# v1.5.0: logcounts densified 25 genes at a time (large imaging sections).
# v1.2.2 changes
#  - Default gene set: top HVGs from scran::modelGeneVar (sparse-aware) instead of
#    as.matrix() on the whole logcounts assay (~1.3 GB dense for one Visium section).
#    Reuses metadata(spe)$hvgs from run_dimred() when present.
#  - Only the selected genes x one sample are ever densified.
#  - alternative = "greater" by default: SVG detection asks for positive spatial
#    autocorrelation (Moran's I > E[I], Geary's C < 1). Two-sided is still available.

build_spatial_weights <- function(spe, method = "kNN", k = 6, d2 = NULL) {
  xy <- SpatialExperiment::spatialCoords(spe)
  ids <- if ("sample_id" %in% names(SummarizedExperiment::colData(spe))) as.character(spe$sample_id) else rep("S1", ncol(spe))
  out <- list()
  for (s in unique(ids)) {
    ii <- which(ids == s); sub <- xy[ii, 1:2, drop = FALSE]
    if (length(ii) < 4L) stop("At least four observations needed for spatial testing in sample ", s)
    if (any(!is.finite(sub))) stop("Non-finite spatial coordinates in ", s)
    if (method == "kNN") {
      nb <- spdep::knn2nb(spdep::knearneigh(sub, k = min(max(1L, as.integer(k)), length(ii) - 1L)), sym = TRUE)
    } else if (method == "distance") {
      if (is.null(d2) || !is.finite(d2) || d2 <= 0) stop("Distance threshold d2 must be provided in coordinate units")
      nb <- spdep::dnearneigh(sub, 0, d2)
      if (any(spdep::card(nb) == 0)) stop("Isolated observations: increase d2 or switch to kNN")
    } else stop("Unknown neighbor method")
    out[[s]] <- list(index = ii, weights = spdep::nb2listw(nb, style = "W", zero.policy = FALSE))
  }
  out
}

compute_moran_geary <- function(spe, genes = NULL, method = "kNN", k = 6, d2 = NULL,
                                adjust_method = "BH", n_genes = 200L, alternative = "greater") {
  if (!"logcounts" %in% SummarizedExperiment::assayNames(spe)) stop("Normalize first")
  alternative <- match.arg(alternative, c("greater", "two.sided", "less"))
  if (is.null(genes)) {
    hv <- S4Vectors::metadata(spe)$hvgs
    hv <- if (length(hv)) intersect(hv, rownames(spe)) else select_hvgs(spe, n = n_genes)  # select_hvgs: 04_clustering.R
    genes <- head(hv, n_genes)
  }
  genes <- intersect(unique(genes), rownames(spe))
  if (!length(genes)) stop("No matching genes")
  weights <- build_spatial_weights(spe, method, k, d2)
  out <- lapply(names(weights), function(s) {
    info <- weights[[s]]
    lc <- SummarizedExperiment::assay(spe, "logcounts")
    # densify 25 genes at a time: a Xenium section can have 100k+ cells
    chunk_of <- (seq_along(genes) - 1L) %/% 25L
    cache <- new.env()
    rows <- lapply(seq_along(genes), function(i) {
      ch <- as.character(chunk_of[i])
      if (!exists(ch, envir = cache, inherits = FALSE)) {
        rm(list = ls(cache), envir = cache)
        assign(ch, as.matrix(lc[genes[chunk_of == chunk_of[i]], info$index, drop = FALSE]), envir = cache)
      }
      x <- as.numeric(get(ch, envir = cache)[genes[i], ]); m <- g <- NULL
      if (all(is.finite(x)) && stats::sd(x) > 0) {
        m <- tryCatch(spdep::moran.test(x, info$weights, alternative = alternative), error = function(e) NULL)
        g <- tryCatch(spdep::geary.test(x, info$weights, alternative = alternative), error = function(e) NULL)
      }
      data.frame(sample = s, gene = genes[i],
                 moran_I = if (is.null(m)) NA_real_ else unname(m$estimate[1]),
                 moran_p = if (is.null(m)) NA_real_ else m$p.value,
                 geary_C = if (is.null(g)) NA_real_ else unname(g$estimate[1]),
                 geary_p = if (is.null(g)) NA_real_ else g$p.value)
    })
    df <- do.call(rbind, rows)
    df$moran_FDR <- p.adjust(df$moran_p, method = adjust_method)
    df$geary_FDR <- p.adjust(df$geary_p, method = adjust_method)
    df
  })
  df <- do.call(rbind, out)
  sym <- get_gene_symbols(spe); names(sym) <- rownames(spe)   # get_gene_symbols: 02_qc.R
  df$symbol <- unname(sym[df$gene])
  df <- df[, c("sample", "gene", "symbol", setdiff(names(df), c("sample", "gene", "symbol")))]
  attr(df, "alternative") <- alternative
  S4Vectors::metadata(spe)$moran_results <- df
  list(spe = spe, results = df, weights = weights)
}

find_svg <- function(spe, genes = NULL, fdr_cutoff = 0.05, ...) {
  res <- compute_moran_geary(spe, genes = genes, ...)
  res$svg <- subset(res$results, !is.na(moran_FDR) & moran_FDR < fdr_cutoff & moran_I > 0)
  res
}
