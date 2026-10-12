# Biological-replicate-level pseudobulk, one model per domain.
#
# v1.3.0: min_spots - a replicate contributes to a domain's pseudobulk only if it has at
#   least min_spots spots/cells in that domain (tiny pseudobulks are very noisy). Domains are
#   meaningful across sections only with joint clustering (cluster_joint, 04_clustering.R).
# v1.5.0: replicate_col - the unit of biological replication. Counts are summed per
#   replicate x domain, so several sections from the same animal/donor are pooled into ONE
#   pseudobulk column instead of being counted as independent samples (pseudoreplication).
#   Default is sample_id (one section = one replicate). The function warns when another
#   column looks like it groups sections from the same animal/donor.
# v1.6.0: covariates - technical factors that span both conditions (e.g. Visium slide,
#   batch) are added to the edgeR design (~ covariates + condition), so their variation is
#   removed from the comparison. The condition coefficient is tested. Per-replicate
#   pseudobulk expression of any gene can be plotted with plot_gene_by_replicate().

# Columns that could be the replication unit: constant within each section, >= 2 values,
# each value within a single condition
replicate_candidates <- function(spe, sample_col = "sample_id", condition_col = "condition") {
  cd <- SummarizedExperiment::colData(spe)
  if (!all(c(sample_col, condition_col) %in% names(cd))) return(sample_col)
  samp <- as.character(cd[[sample_col]]); cond <- as.character(cd[[condition_col]])
  skip <- c(condition_col, "STclust", "cluster_louvain", "platform", "barcode", "in_tissue")
  ok <- vapply(setdiff(names(cd), skip), function(cc) {
    v <- cd[[cc]]
    if (!is.atomic(v) || anyNA(v)) return(FALSE)
    v <- as.character(v)
    map <- unique(data.frame(s = samp, v = v))
    if (anyDuplicated(map$s)) return(FALSE)                       # varies within a section
    if (length(unique(v)) < 2L) return(FALSE)
    all(tapply(cond, v, function(z) length(unique(z))) == 1L)      # nested in condition
  }, logical(1))
  unique(c(sample_col, names(ok)[ok]))
}

# Columns that group several sections into fewer units within a condition: likely
# animal/donor IDs, meaning the sections are not independent replicates.
detect_pseudoreplication <- function(spe, replicate_col, sample_col = "sample_id", condition_col = "condition") {
  cd <- SummarizedExperiment::colData(spe)
  rep_n <- length(unique(cd[[replicate_col]]))
  cond <- as.character(cd[[condition_col]])
  cands <- setdiff(replicate_candidates(spe, sample_col, condition_col), c(sample_col, replicate_col))
  cands[vapply(cands, function(cc) {
    v <- as.character(cd[[cc]])
    n <- length(unique(v))
    # coarser than the chosen unit, and not just a relabelling of condition itself
    n < rep_n && n > length(unique(cond))
  }, logical(1))]
}

# v1.6.1: explicit reference (baseline) condition. Previously the alphabetically first level
#   was the baseline, so "6OHDA" (digit first) became the reference and every logFC was
#   reversed relative to the usual "treated vs control" reading. logFC is now always
#   <other condition> vs <reference>, recorded in the `contrast` column of the results.
CONTROL_LIKE <- "^(control|controls|ctrl|ctl|vehicle|veh|sham|wt|wild[ _-]?type|untreated|baseline|naive|healthy|normal|mock)$"

guess_reference <- function(conds) {
  conds <- sort(unique(as.character(conds)))
  hit <- conds[grepl(CONTROL_LIKE, conds, ignore.case = TRUE)]
  if (!length(hit)) hit <- conds[grepl("control|ctrl|vehicle|sham", conds, ignore.case = TRUE)]
  if (length(hit)) hit[1] else conds[1]
}

# Columns usable as covariates: constant within each section, >= 2 values, and present in
# both conditions (a technical factor such as slide or batch, not a replicate ID)
covariate_candidates <- function(spe, sample_col = "sample_id", condition_col = "condition") {
  cd <- SummarizedExperiment::colData(spe)
  if (!all(c(sample_col, condition_col) %in% names(cd))) return(character(0))
  samp <- as.character(cd[[sample_col]]); cond <- as.character(cd[[condition_col]])
  skip <- c(sample_col, condition_col, "STclust", "cluster_louvain", "platform", "barcode", "in_tissue")
  ok <- vapply(setdiff(names(cd), skip), function(cc) {
    v <- cd[[cc]]
    if (!is.atomic(v) || anyNA(v)) return(FALSE)
    v <- as.character(v)
    if (anyDuplicated(unique(data.frame(s = samp, v = v))$s)) return(FALSE)
    n <- length(unique(v))
    n >= 2L && n < length(unique(samp)) && any(tapply(cond, v, function(z) length(unique(z))) > 1L)
  }, logical(1))
  names(ok)[ok]
}

run_pseudobulk_de <- function(spe, condition_col = "condition", sample_col = "sample_id",
                              cluster_col = NULL, min_spots = 10L, replicate_col = sample_col,
                              covariates = character(0), reference = NULL) {
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  if (!all(c(condition_col, sample_col) %in% names(cd))) stop("condition and sample columns required")
  if (!replicate_col %in% names(cd)) stop("Replicate column '", replicate_col, "' not found")
  if (!is.null(cluster_col) && !cluster_col %in% names(cd)) stop("Cluster column missing")
  if (!is.null(cluster_col) && isTRUE(S4Vectors::metadata(spe)$spatial_clustering$sample_specific))
    warning("Domains were clustered per sample, so labels are not shared across sections; ",
            "use joint clustering for per-domain DE")
  # replicate must be constant within each section
  if (anyDuplicated(unique(data.frame(s = as.character(cd[[sample_col]]), r = as.character(cd[[replicate_col]])))$s))
    stop("Replicate column '", replicate_col, "' varies within a section; it must identify a section, animal or donor")
  pr <- detect_pseudoreplication(spe, replicate_col, sample_col, condition_col)
  if (length(pr))
    warning("Column(s) ", paste(pr, collapse = ", "), " group several sections within a condition. ",
            "If they identify the animal/donor, those sections are not independent replicates: ",
            "choose that column as the replicate column.")

  covariates <- setdiff(as.character(covariates), "")
  bad <- setdiff(covariates, names(cd))
  if (length(bad)) stop("Covariate column(s) not found: ", paste(bad, collapse = ", "))
  if (length(intersect(covariates, c(condition_col, replicate_col))))
    stop("The condition or replicate column cannot also be a covariate")
  df <- data.frame(replicate = as.character(cd[[replicate_col]]), condition = as.character(cd[[condition_col]]),
                   stringsAsFactors = FALSE)
  for (cv in covariates) df[[cv]] <- as.character(cd[[cv]])
  if (anyNA(df) || any(!nzchar(df$replicate)) || any(!nzchar(df$condition))) stop("Missing replicate/condition/covariate labels")
  mapping <- unique(df)
  for (cv in covariates) if (anyDuplicated(unique(mapping[, c("replicate", cv)])$replicate))
    stop("Covariate '", cv, "' varies within a replicate (", replicate_col, "); a covariate must be constant per replicate")
  if (anyDuplicated(mapping$replicate)) stop("Each replicate must map to exactly one condition")
  conds <- sort(unique(mapping$condition))
  if (length(conds) != 2L) stop("Exactly two conditions are currently supported")
  if (is.null(reference) || !nzchar(reference)) reference <- guess_reference(conds)
  if (!reference %in% conds) stop("Reference condition '", reference, "' not found; conditions are: ",
                                  paste(conds, collapse = ", "))
  conds <- c(reference, setdiff(conds, reference))   # baseline first: logFC = conds[2] vs conds[1]
  contrast <- paste(conds[2], "vs", conds[1])

  counts_by_group <- function(sub) {
    m <- SummarizedExperiment::assay(sub, "counts")
    idx <- split(seq_len(ncol(sub)), as.character(sub[[replicate_col]]))
    idx <- idx[lengths(idx) >= min_spots]
    if (!length(idx)) return(NULL)
    z <- vapply(idx, function(j) Matrix::rowSums(m[, j, drop = FALSE]), numeric(nrow(sub)))
    if (is.null(dim(z))) z <- matrix(z, ncol = 1, dimnames = list(rownames(sub), names(idx)))
    z
  }

  levels <- if (is.null(cluster_col)) "all" else sort(unique(as.character(cd[[cluster_col]])))
  result <- list(); bulk <- list(); summary_rows <- list()
  for (cl in levels) {
    sub <- if (is.null(cluster_col)) spe else spe[, !is.na(cd[[cluster_col]]) & as.character(cd[[cluster_col]]) == cl, drop = FALSE]
    mat <- if (ncol(sub)) counts_by_group(sub) else NULL
    meta <- if (is.null(mat)) mapping[0, ] else mapping[match(colnames(mat), mapping$replicate), , drop = FALSE]
    tc <- table(factor(meta$condition, levels = conds))
    tested <- all(tc >= 2L)
    row <- data.frame(domain = cl, spots = ncol(sub))
    row[[paste0("n_", conds[1])]] <- unname(tc[1]); row[[paste0("n_", conds[2])]] <- unname(tc[2])
    row$tested <- tested
    summary_rows[[cl]] <- row
    if (!tested) {
      warning("Skipping domain ", cl, ": need >=2 replicates (", replicate_col, ") per condition with >=",
              min_spots, " spots each")
      next
    }
    dmeta <- data.frame(group = factor(meta$condition, levels = conds))
    for (cv in covariates) dmeta[[cv]] <- droplevels(factor(meta[[cv]]))
    terms <- c(covariates[vapply(covariates, function(cv) nlevels(dmeta[[cv]]) > 1L, logical(1))], "group")
    design <- model.matrix(stats::as.formula(paste("~", paste(terms, collapse = " + "))), data = dmeta)
    coef <- ncol(design)   # group is the last term
    if (qr(design)$rank < ncol(design) || nrow(design) - ncol(design) < 1L) {
      warning("Skipping domain ", cl, ": too few replicates for the design (", paste(terms, collapse = " + "),
              "); remove covariates or lower min_spots")
      summary_rows[[cl]]$tested <- FALSE
      next
    }
    y <- edgeR::DGEList(counts = mat)
    keep <- edgeR::filterByExpr(y, design = design)
    if (sum(keep) < 2L) next
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- if (exists("normLibSizes", envir = asNamespace("edgeR"))) edgeR::normLibSizes(y) else edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design, robust = FALSE)
    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit, coef = coef)
    tt <- edgeR::topTags(qlf, n = Inf)$table
    tt$gene <- rownames(tt); tt$domain <- cl; tt$contrast <- contrast
    result[[cl]] <- tt
    bulk[[cl]] <- list(counts = mat, replicate_metadata = meta)
  }
  domain_summary <- do.call(rbind, summary_rows); rownames(domain_summary) <- NULL
  if (!length(result)) stop("No domains support a valid replicated comparison (>=2 replicates of '", replicate_col,
                            "' per condition, >=", min_spots, " spots each)")
  combined <- do.call(rbind, result)
  sym <- get_gene_symbols(spe); names(sym) <- rownames(spe)   # get_gene_symbols: 02_qc.R
  combined$symbol <- unname(sym[combined$gene])
  combined <- combined[, c("gene", "symbol", "domain", setdiff(names(combined), c("gene", "symbol", "domain")))]
  combined$FDR_within_domain <- combined$FDR
  combined$FDR <- p.adjust(combined$PValue, method = "BH")
  rownames(combined) <- NULL
  list(results = combined, pseudobulk = bulk, domain_summary = domain_summary, replicate_col = replicate_col,
       covariates = covariates, reference = conds[1], contrast = contrast)
}

# Per-replicate pseudobulk expression (log2 CPM) of one gene in one domain, by condition.
# Shows whether a DE result is consistent across replicates or driven by a few of them.
plot_gene_by_replicate <- function(de, gene, domain = names(de$pseudobulk)[1], label = gene) {
  pb <- de$pseudobulk[[domain]]
  if (is.null(pb)) stop("No pseudobulk for domain ", domain)
  if (!gene %in% rownames(pb$counts)) stop("Gene not in pseudobulk: ", gene)
  lib <- colSums(pb$counts)
  d <- data.frame(replicate = colnames(pb$counts),
                  logCPM = log2((pb$counts[gene, ] + 0.5) / (lib + 1) * 1e6),
                  condition = factor(pb$replicate_metadata$condition,
                                     levels = c(de$reference, setdiff(unique(pb$replicate_metadata$condition), de$reference))))
  ggplot2::ggplot(d, ggplot2::aes(condition, logCPM, colour = condition)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey50") +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1), size = 2.5) +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = paste0(label, " (domain: ", domain, ")"),
                  subtitle = "Each point is one replicate's pseudobulk", y = "log2 CPM", x = NULL)
}

plot_volcano <- function(de_res, fdr_cutoff = 0.05, lfc_cutoff = 0.5) {
  de_res$significant <- !is.na(de_res$FDR) & de_res$FDR < fdr_cutoff & abs(de_res$logFC) > lfc_cutoff
  ggplot2::ggplot(de_res, ggplot2::aes(x = logFC, y = -log10(pmax(PValue, 1e-300)), colour = significant)) +
    ggplot2::geom_point(alpha = 0.6) + ggplot2::theme_minimal() +
    ggplot2::labs(title = "Pseudobulk differential expression")
}
