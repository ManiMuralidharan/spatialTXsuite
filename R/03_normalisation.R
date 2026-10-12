# v1.2.1: normalize_spe() clears previous size factors, logcounts, PCA and
# clusters first. Previously, switching scran -> logNormCounts silently reused
# the scran size factors, and old PCA/clusters survived re-normalisation.
# v1.4.1: log-normalisation computed directly instead of the deprecated scuttle
# librarySizeFactors/normalizeCounts. Same result as scuttle::logNormCounts():
# size factors centred to mean 1, logcounts = log2(counts / size_factor + 1), sparsity kept.
log_norm_counts <- function(spe) {
  m <- SummarizedExperiment::assay(spe, "counts")
  sf <- BiocGenerics::sizeFactors(spe)
  if (is.null(sf)) sf <- Matrix::colSums(m)
  sf <- as.numeric(sf)
  if (any(!is.finite(sf)) || any(sf <= 0)) stop("Size factors must be positive; remove spots with zero counts first")
  sf <- sf / mean(sf)
  BiocGenerics::sizeFactors(spe) <- sf
  if (methods::is(m, "sparseMatrix")) {
    lc <- methods::as(m %*% Matrix::Diagonal(x = 1 / sf), "CsparseMatrix")
    lc@x <- log2(lc@x + 1)
  } else {
    lc <- log2(sweep(as.matrix(m), 2, sf, "/") + 1)
  }
  dimnames(lc) <- dimnames(m)
  SummarizedExperiment::assay(spe, "logcounts") <- lc
  spe
}

normalize_spe <- function(spe, method = "logNormCounts") {
  stopifnot("counts" %in% SummarizedExperiment::assayNames(spe))
  spe <- clear_downstream(spe)  # defined in 02_qc.R
  if (method == "logNormCounts") {
    spe <- log_norm_counts(spe)
  } else if (method == "scran") {
    if (ncol(spe) < 50L) stop("scran normalization requires more observations; use logNormCounts")
    cl <- scran::quickCluster(spe)
    spe <- scran::computeSumFactors(spe, clusters = cl)
    spe <- log_norm_counts(spe)
  } else if (method == "SCT") {
    stop("SCT integration is not validated for this SpatialExperiment workflow; use logNormCounts or scran")
  } else stop("Unknown normalization method: ", method)
  S4Vectors::metadata(spe)$norm_method <- method
  spe
}

scale_spots <- function(mat, method = "zscore") {
  m <- as.matrix(mat)
  if (method == "zscore") {
    out <- t(scale(t(m)))
    out[!is.finite(out)] <- 0
  } else if (method == "0-1") {
    mn <- apply(m, 1, min, na.rm = TRUE)
    mx <- apply(m, 1, max, na.rm = TRUE)
    range <- pmax(mx - mn, .Machine$double.eps)
    out <- sweep(sweep(m, 1, mn, "-"), 1, range, "/")
  } else stop("Unknown method")
  dimnames(out) <- dimnames(m)
  out
}

diag_normalisation <- function(spe, max_spots = 50, max_genes = 2000) {
  if (!"logcounts" %in% SummarizedExperiment::assayNames(spe)) stop("Normalize first")
  i <- seq_len(min(nrow(spe), max_genes)); j <- seq_len(min(ncol(spe), max_spots))
  raw <- log2(as.matrix(SummarizedExperiment::assay(spe, "counts")[i, j, drop = FALSE]) + 1)
  after <- as.matrix(SummarizedExperiment::assay(spe, "logcounts")[i, j, drop = FALSE])
  rle <- function(x) sweep(x, 1, apply(x, 1, median, na.rm = TRUE), "-")
  list(rle_before = rle(raw), rle_after = rle(after))
}

correct_batch <- function(spe, batch_col = "sample_id", condition_col = NULL, method = "Harmony") {
  if (method != "Harmony") stop("ComBat is disabled until confounding-aware design and validation are implemented")
  if (!requireNamespace("harmony", quietly = TRUE)) stop("Install harmony to enable batch correction")
  if (!batch_col %in% names(SummarizedExperiment::colData(spe))) stop("Batch column absent")
  if (!"PCA" %in% SingleCellExperiment::reducedDimNames(spe))
    spe <- scater::runPCA(spe, ncomponents = min(20, ncol(spe) - 1, nrow(spe) - 1))
  if (!is.null(condition_col)) {
    tab <- table(spe[[batch_col]], spe[[condition_col]])
    # Condition is normally constant within a sample, so that is not confounding.
    # Confounding = a condition level observed in only one batch (no replication).
    if (any(colSums(tab > 0) < 2L)) warning("A condition level occurs in only one batch: condition is confounded with batch")
  }
  z <- harmony::RunHarmony(SingleCellExperiment::reducedDim(spe, "PCA"),
                           meta_data = as.data.frame(SummarizedExperiment::colData(spe)),
                           vars_use = batch_col, verbose = FALSE)
  SingleCellExperiment::reducedDim(spe, "PCA_harmony") <- z
  spe
}
