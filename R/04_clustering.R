# v1.2.2 changes
#  - run_dimred(): PCA on top HVGs (scran::modelGeneVar, blocked by sample) instead of
#    all genes; UMAP is optional (off by default, nothing downstream uses it).
#  - cluster_spatial(): no dense n x n distance matrices. Spatial and expression
#    information are combined in one augmented feature space and neighbours are found
#    with BiocNeighbors::findKNN (memory O(n*k), feasible for Visium HD bins).
#  - Coordinates are scaled by ONE common factor, so tissue geometry is not distorted
#    (previously x and y were standardised separately).

select_hvgs <- function(spe, n = 2000L, block_col = "sample_id") {
  if (!"logcounts" %in% SummarizedExperiment::assayNames(spe)) stop("Normalize first")
  cd <- SummarizedExperiment::colData(spe)
  block <- if (block_col %in% names(cd) && length(unique(cd[[block_col]])) > 1L)
    factor(cd[[block_col]]) else NULL
  # v1.4.2: computed directly (same idea as scran::modelGeneVar + getTopHVGs, which are
  # deprecated). Per block: variance of logcounts minus a lowess trend of variance on mean
  # (the technical part); the excess ("biological") variance is averaged across blocks.
  bio <- gene_bio_variance(SummarizedExperiment::assay(spe, "logcounts"), block)
  names(bio) <- rownames(spe)
  bio <- bio[is.finite(bio)]
  names(sort(bio, decreasing = TRUE))[seq_len(min(as.integer(n), length(bio)))]
}

# Technical-variance trend: median variance within equal-count bins of mean expression,
# interpolated. Medians cannot be dragged by a minority of truly variable genes, unlike
# lowess, which gets pulled through high-leverage genes at the top of the expression range.
robust_trend <- function(m, v, genes_per_bin = 20L, max_bins = 50L) {
  nb <- max(1L, min(as.integer(max_bins), floor(length(m) / genes_per_bin)))
  if (nb < 2L) return(rep(stats::median(v), length(m)))
  br <- unique(stats::quantile(m, probs = seq(0, 1, length.out = nb + 1L), names = FALSE))
  bin <- cut(m, br, include.lowest = TRUE)
  bx <- tapply(m, bin, stats::median)
  by <- tapply(v, bin, stats::median)
  keep <- is.finite(bx) & is.finite(by)
  if (sum(keep) < 2L) return(rep(stats::median(v), length(m)))
  stats::approx(bx[keep], by[keep], xout = m, rule = 2, ties = mean)$y
}

gene_bio_variance <- function(lc, block = NULL) {
  groups <- if (is.null(block)) list(seq_len(ncol(lc))) else split(seq_len(ncol(lc)), block)
  per_block <- vapply(groups, function(j) {
    x <- lc[, j, drop = FALSE]
    nj <- length(j)
    if (nj < 3L) return(rep(NA_real_, nrow(lc)))
    m <- as.numeric(Matrix::rowMeans(x))
    v <- (as.numeric(Matrix::rowMeans(x^2)) - m^2) * nj / (nj - 1)
    ok <- m > 0 & is.finite(v)
    out <- rep(NA_real_, length(m))
    if (sum(ok) < 3L) { out[ok] <- v[ok]; return(out) }
    out[ok] <- v[ok] - robust_trend(m[ok], v[ok])
    out
  }, numeric(nrow(lc)))
  if (is.null(dim(per_block))) per_block <- matrix(per_block, ncol = 1)
  avg <- rowMeans(per_block, na.rm = TRUE)
  avg[is.nan(avg)] <- NA_real_
  avg
}

run_dimred <- function(spe, n_pcs = 30, n_hvgs = 2000L, run_umap = FALSE) {
  if (!"logcounts" %in% SummarizedExperiment::assayNames(spe)) stop("Normalize first")
  hvgs <- select_hvgs(spe, n = n_hvgs)
  n_pcs <- min(as.integer(n_pcs), length(hvgs) - 1L, ncol(spe) - 1L)
  if (n_pcs < 2L) stop("Need at least 3 informative genes and 3 observations")
  spe <- scater::runPCA(spe, ncomponents = n_pcs, subset_row = hvgs)
  if (run_umap) spe <- scater::runUMAP(spe, dimred = "PCA")
  S4Vectors::metadata(spe)$hvgs <- hvgs
  spe
}

cluster_standard <- function(spe, k = 15, method = "louvain") {
  if (!"PCA" %in% SingleCellExperiment::reducedDimNames(spe)) spe <- run_dimred(spe)
  k <- min(max(2L, as.integer(k)), ncol(spe) - 1L)
  graph <- scran::buildSNNGraph(spe, k = k, use.dimred = "PCA")
  membership <- if (method == "leiden") igraph::cluster_leiden(graph, objective_function = "modularity")$membership
                else igraph::cluster_louvain(graph)$membership
  spe$cluster_louvain <- factor(membership)
  spe
}

# Root-mean-square spread of a matrix's columns: puts blocks on a comparable scale
.block_scale <- function(m) {
  s <- sqrt(sum(apply(m, 2, stats::var)))
  if (!is.finite(s) || s <= 0) 1 else s
}

cluster_spatial <- function(spe, k = 7, method = "spatial_SNN", spatial_weight = 0.5) {
  if (method == "louvain") {
    spe <- cluster_standard(spe, k = k)
    spe$STclust <- spe$cluster_louvain
    return(spe)
  }
  if (method != "spatial_SNN") stop("Only spatial_SNN and louvain are implemented; STclust/BayesSpace are not yet validated")
  if (!is.finite(spatial_weight) || spatial_weight < 0 || spatial_weight > 1) stop("spatial_weight must be in [0, 1]")
  if (!"PCA" %in% SingleCellExperiment::reducedDimNames(spe)) spe <- run_dimred(spe)
  coords <- SpatialExperiment::spatialCoords(spe)
  if (ncol(spe) < 3L || ncol(coords) < 2L || any(!is.finite(coords[, 1:2])))
    stop("Valid coordinates and at least 3 spots are required")
  sid <- if ("sample_id" %in% names(SummarizedExperiment::colData(spe))) as.character(spe$sample_id) else rep("sample", ncol(spe))
  labels <- integer(ncol(spe)); offset <- 0L
  for (s in unique(sid)) {
    idx <- which(sid == s)
    if (length(idx) < 3L) stop("At least 3 spots required in sample ", s)
    xy <- coords[idx, 1:2, drop = FALSE]
    xy <- sweep(xy, 2, colMeans(xy))
    xy <- xy / .block_scale(xy)              # single factor: geometry preserved
    pc <- SingleCellExperiment::reducedDim(spe, "PCA")[idx, , drop = FALSE]
    pc <- pc / .block_scale(pc)
    # Squared distance in this space = w * d_spatial^2 + (1 - w) * d_expression^2
    feat <- cbind(sqrt(spatial_weight) * xy, sqrt(1 - spatial_weight) * pc)
    kk <- min(max(1L, as.integer(k)), length(idx) - 1L)
    nn <- BiocNeighbors::findKNN(feat, k = kk)$index
    edge <- cbind(rep(seq_along(idx), each = kk), as.vector(t(nn)))
    graph <- igraph::simplify(igraph::graph_from_edgelist(edge, directed = FALSE))
    memberships <- igraph::cluster_louvain(graph)$membership
    labels[idx] <- memberships + offset; offset <- max(labels)
  }
  spe$STclust <- factor(labels)
  S4Vectors::metadata(spe)$spatial_clustering <- list(method = method, k = k, spatial_weight = spatial_weight,
                                                      sample_specific = TRUE, neighbours = "BiocNeighbors::findKNN")
  spe
}

# ---------------------------------------------------------------------------
# v1.3.0: joint (cross-sample) spatial domains
#
# Per-sample clustering gives labels that exist in one section only, so
# per-domain pseudobulk DE can never have replicates. cluster_joint():
#   1. PCA on HVGs selected with sample blocking (run_dimred).
#   2. Harmony correction of the PCs by sample_id (if harmony is installed).
#   3. Spatial smoothing in the spirit of BANKSY: each spot's corrected PCs are
#      concatenated with the mean PCs of its spatial neighbours, which are found
#      WITHIN its own section (coordinates are never compared across sections).
#   4. One kNN graph over ALL spots and one Louvain run, so the same label means
#      the same domain in every section.
# spatial_weight = 0 -> expression only; 1 -> neighbourhood average only.
# ---------------------------------------------------------------------------
cluster_joint <- function(spe, k = 15, spatial_weight = 0.5, spatial_k = 6,
                          batch_col = "sample_id", use_harmony = TRUE) {
  if (!is.finite(spatial_weight) || spatial_weight < 0 || spatial_weight > 1) stop("spatial_weight must be in [0, 1]")
  if (!"PCA" %in% SingleCellExperiment::reducedDimNames(spe)) spe <- run_dimred(spe)
  cd <- SummarizedExperiment::colData(spe)
  sid <- if (batch_col %in% names(cd)) as.character(cd[[batch_col]]) else rep("S1", ncol(spe))
  n_samples <- length(unique(sid))

  dimred <- "PCA"; corrected <- FALSE
  if (use_harmony && n_samples > 1L) {
    if (requireNamespace("harmony", quietly = TRUE)) {
      spe <- correct_batch(spe, batch_col = batch_col,
                           condition_col = if ("condition" %in% names(cd)) "condition" else NULL)  # 03_normalisation.R
      dimred <- "PCA_harmony"; corrected <- TRUE
    } else {
      warning("harmony not installed: joint clustering uses uncorrected PCs, so domains may split by section. ",
              "Install with install.packages('harmony')")
    }
  }
  pc <- SingleCellExperiment::reducedDim(spe, dimred)
  pc <- pc / .block_scale(pc)

  coords <- SpatialExperiment::spatialCoords(spe)
  if (any(!is.finite(coords[, 1:2]))) stop("Non-finite spatial coordinates")
  nbr_mean <- matrix(0, nrow(pc), ncol(pc))
  for (s in unique(sid)) {
    idx <- which(sid == s)
    if (length(idx) < 3L) stop("At least 3 spots required in sample ", s)
    sk <- min(max(1L, as.integer(spatial_k)), length(idx) - 1L)
    nn <- BiocNeighbors::findKNN(coords[idx, 1:2, drop = FALSE], k = sk)$index
    p <- pc[idx, , drop = FALSE]
    nbr_mean[idx, ] <- vapply(seq_len(ncol(p)), function(j) rowMeans(matrix(p[nn, j], ncol = sk)), numeric(length(idx)))
  }
  feat <- cbind(sqrt(1 - spatial_weight) * pc, sqrt(spatial_weight) * nbr_mean)
  SingleCellExperiment::reducedDim(spe, "joint_features") <- feat

  kk <- min(max(2L, as.integer(k)), ncol(spe) - 1L)
  nn <- BiocNeighbors::findKNN(feat, k = kk)$index
  edge <- cbind(rep(seq_len(ncol(spe)), each = kk), as.vector(t(nn)))
  graph <- igraph::simplify(igraph::graph_from_edgelist(edge, directed = FALSE))
  spe$STclust <- factor(igraph::cluster_louvain(graph)$membership)

  S4Vectors::metadata(spe)$spatial_clustering <- list(
    method = "joint", k = k, spatial_weight = spatial_weight, spatial_k = spatial_k,
    sample_specific = FALSE, harmony = corrected, dimred = dimred)
  spe
}

# Domain x sample spot counts: shows whether each domain is shared across sections
domain_sample_table <- function(spe, cluster_col = "STclust", sample_col = "sample_id") {
  as.data.frame.matrix(table(domain = spe[[cluster_col]], sample = spe[[sample_col]]))
}

find_markers <- function(spe, cluster_col = "STclust") {
  if (!cluster_col %in% names(SummarizedExperiment::colData(spe))) stop("Cluster annotations absent")
  scran::findMarkers(spe, groups = spe[[cluster_col]])
}
