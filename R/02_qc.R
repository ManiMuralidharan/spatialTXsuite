# QC and filtering for SpatialExperiment objects.
#
# v1.2.1 changes
#  - Gene symbols are taken from rowData (symbol / gene_name ...) when present.
#    read10xVisium() stores Ensembl IDs in rownames, so grepping rownames for
#    "^MT-" found nothing on real data and silently disabled the mito filter.
#  - compute_qc() removes stale QC columns before recomputing (re-loaded RDS).
#  - clear_downstream() drops results computed on a different spot set;
#    filter_spots() calls it so normalisation/PCA/clusters are never stale.

get_gene_symbols <- function(spe) {
  rd <- SummarizedExperiment::rowData(spe)
  for (col in c("symbol", "Symbol", "SYMBOL", "gene_name", "gene_symbol")) {
    if (col %in% names(rd)) {
      sym <- as.character(rd[[col]])
      if (any(!is.na(sym) & nzchar(sym))) {
        miss <- is.na(sym) | !nzchar(sym)
        sym[miss] <- rownames(spe)[miss]
        return(sym)
      }
    }
  }
  rownames(spe)
}

detect_mito_genes <- function(spe, organism = "auto") {
  genes <- get_gene_symbols(spe)
  if (organism == "human") return(grep("^MT-", genes))
  if (organism == "mouse") return(grep("^mt-", genes))
  grep("^(MT-|mt-)", genes)
}

detect_ribo_genes <- function(spe) {
  grep("^(RP[SL][0-9]|Rp[sl][0-9])", get_gene_symbols(spe))
}

# Remove anything derived from a previous spot set or normalisation.
# keep_logcounts = TRUE is used when only clustering/stat results must go.
clear_downstream <- function(spe, keep_logcounts = FALSE) {
  if (!keep_logcounts) {
    if ("logcounts" %in% SummarizedExperiment::assayNames(spe))
      SummarizedExperiment::assay(spe, "logcounts") <- NULL
    # logNormCounts() silently reuses existing size factors, so they must go too
    BiocGenerics::sizeFactors(spe) <- NULL
  }
  SingleCellExperiment::reducedDims(spe) <- list()
  for (x in intersect(c("STclust", "cluster_louvain"),
                      names(SummarizedExperiment::colData(spe))))
    SummarizedExperiment::colData(spe)[[x]] <- NULL
  md <- S4Vectors::metadata(spe)
  md[intersect(c("spatial_clustering", "moran_results",
                 if (!keep_logcounts) "norm_method"), names(md))] <- NULL
  S4Vectors::metadata(spe) <- md
  spe
}

# v1.4.1: per-spot QC metrics computed directly (sparse-aware) instead of the deprecated
# scater::addPerCellQC. Column names are identical to scater's, so nothing downstream changes.
add_qc_metrics <- function(spe, subsets = list()) {
  m <- SummarizedExperiment::assay(spe, "counts")
  total <- Matrix::colSums(m)
  spe$sum <- as.numeric(total)
  spe$detected <- as.integer(Matrix::colSums(m > 0))
  for (nm in names(subsets)) {
    sub <- m[subsets[[nm]], , drop = FALSE]
    s <- as.numeric(Matrix::colSums(sub))
    spe[[paste0("subsets_", nm, "_sum")]] <- s
    spe[[paste0("subsets_", nm, "_detected")]] <- as.integer(Matrix::colSums(sub > 0))
    spe[[paste0("subsets_", nm, "_percent")]] <- ifelse(total > 0, 100 * s / total, NA_real_)
  }
  spe$total <- as.numeric(total)
  spe
}

compute_qc <- function(spe, organism = "auto") {
  stopifnot(inherits(spe, "SingleCellExperiment"),
            "counts" %in% SummarizedExperiment::assayNames(spe))
  cd <- SummarizedExperiment::colData(spe)
  old <- grep("^(sum|detected|total|subsets_.*)$", names(cd), value = TRUE)
  for (x in old) SummarizedExperiment::colData(spe)[[x]] <- NULL

  mito <- detect_mito_genes(spe, organism)
  ribo <- detect_ribo_genes(spe)
  subsets <- list()
  if (length(mito)) subsets$Mito <- mito
  if (length(ribo)) subsets$Ribo <- ribo
  spe <- add_qc_metrics(spe, subsets)
  # Imaging panels (Xenium/MERFISH/CosMx) normally contain no mitochondrial genes: not a problem
  if (!length(mito) && !(exists("is_imaging", mode = "function") && is_imaging(spe)))   # 01_import.R
    warning("No mitochondrial genes detected (checked rowData symbols and rownames); ",
            "mitochondrial filter will be unavailable")
  S4Vectors::metadata(spe)$qc_params <- list(organism = organism,
                                             mito_genes = length(mito),
                                             ribo_genes = length(ribo),
                                             time = Sys.time())
  spe
}

get_qc_thresholds <- function(platform = "Visium") {
  switch(platform,
         Visium   = list(min_counts = 500, min_features = 250, max_mito = 20),
         Xenium   = list(min_counts = 20,  min_features = 10,  max_mito = NA_real_),
         CosMx    = list(min_counts = 20,  min_features = 10,  max_mito = NA_real_),
         VisiumHD = list(min_counts = 50,  min_features = 20,  max_mito = 20),
         list(min_counts = 500, min_features = 250, max_mito = 20))
}

filter_spots <- function(spe, platform = "Visium", min_counts = NULL, min_features = NULL,
                         max_mito = NULL, return_audit = TRUE) {
  cd <- SummarizedExperiment::colData(spe)
  if (!"sum" %in% names(cd)) stop("Run compute_qc() first")
  defaults <- get_qc_thresholds(platform)
  if (is.null(min_counts))   min_counts   <- defaults$min_counts
  if (is.null(min_features)) min_features <- defaults$min_features
  if (is.null(max_mito))     max_mito     <- defaults$max_mito

  has_mito <- "subsets_Mito_percent" %in% names(cd)
  mito <- if (has_mito) spe$subsets_Mito_percent else rep(NA_real_, ncol(spe))
  imaging <- exists("is_imaging", mode = "function") && is_imaging(spe)   # 01_import.R
  if (!has_mito && !is.na(max_mito) && !imaging)
    warning("Mitochondrial threshold requested but no mito genes were detected; ",
            "mito filter not applied")

  low_counts   <- is.na(spe$sum) | spe$sum < min_counts
  low_features <- is.na(spe$detected) | spe$detected < min_features
  high_mito    <- if (is.na(max_mito)) rep(FALSE, ncol(spe)) else !is.na(mito) & mito > max_mito
  discard <- low_counts | low_features | high_mito

  audit <- data.frame(
    spot = colnames(spe),
    sample_id = if ("sample_id" %in% names(cd)) as.character(spe$sample_id) else "S1",
    total = spe$sum, detected = spe$detected, mito_percent = mito,
    low_counts = low_counts, low_features = low_features, high_mito = high_mito,
    discard = discard)

  if (all(discard)) stop("QC filtering removed all spots; adjust thresholds")
  out <- spe[, !discard, drop = FALSE]
  out <- clear_downstream(out)  # anything computed before filtering is now invalid
  S4Vectors::metadata(out)$filter_audit <- audit
  S4Vectors::metadata(out)$filter_thresholds <- list(platform = platform,
    min_counts = min_counts, min_features = min_features, max_mito = max_mito,
    mito_filter_applied = has_mito && !is.na(max_mito))
  if (return_audit) list(spe = out, audit = audit) else out
}

pseudo_bulk_cor <- function(spe) {
  if (!"sample_id" %in% names(SummarizedExperiment::colData(spe))) stop("sample_id required")
  pb <- scuttle::aggregateAcrossCells(spe, ids = S4Vectors::DataFrame(sample_id = spe$sample_id),
                                      use.assay.type = "counts")
  pb <- scuttle::logNormCounts(pb)
  stats::cor(SummarizedExperiment::assay(pb, "logcounts"))
}
