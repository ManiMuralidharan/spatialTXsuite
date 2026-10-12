# Run from project root with Rscript tests/smoke_test.R
library(SpatialExperiment)
for (f in c("01_import.R", "02_qc.R", "03_normalisation.R", "04_clustering.R", "05_spatial_stats.R", "07_diff_analysis.R"))
  source(file.path("R", f))
set.seed(1)
expect_error <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))
capture_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  w
}

# Real Visium layout: Ensembl IDs in rownames, symbols in rowData$symbol
mat <- matrix(rpois(300 * 200, lambda = 6), nrow = 300)  # 300 genes: enough for a variance trend
symbols <- c("MT-CO1", "RPL3", paste0("GENE", seq_len(298)))
rownames(mat) <- sprintf("ENSG%011d", seq_len(nrow(mat)))
colnames(mat) <- paste0("SP", seq_len(ncol(mat)))
xy <- cbind(x = rep(seq_len(20), 10), y = rep(seq_len(10), each = 20))
# Planted signal: gene 3 follows a strong left-to-right gradient in every sample
mat[3, ] <- rpois(ncol(mat), lambda = 2 + 2 * xy[, "x"])
gradient_gene <- rownames(mat)[3]
# Planted domain: genes 4-8 strongly up in the left half (x <= 10) of every section
left <- xy[, "x"] <= 10
mat[4:8, left] <- rpois(5 * sum(left), lambda = 30)
# sample_id must be supplied in colData at construction: SpatialExperiment refuses to
# change the number of unique sample_ids afterwards with spe$sample_id <- ...
spe <- SpatialExperiment(
  assays = list(counts = mat), spatialCoords = xy,
  rowData = S4Vectors::DataFrame(symbol = symbols),
  colData = S4Vectors::DataFrame(
    sample_id = rep(c("A1", "A2", "B1", "B2"), each = 50),
    condition = rep(c("control", "control", "treated", "treated"), each = 50),
    animal = rep(c("m1", "m1", "m2", "m3"), each = 50),   # A1 and A2 are sections of the same animal
    slide = rep(c("s1", "s2", "s1", "s2"), each = 50),     # each slide carries one section per condition
    row.names = colnames(mat)))

# Fix 1: mito/ribo detected via rowData symbols, not rownames
spe <- compute_qc(spe)
stopifnot(length(detect_mito_genes(spe)) == 1L,
          length(detect_ribo_genes(spe)) == 1L,
          "subsets_Mito_percent" %in% names(colData(spe)))
# Re-running QC must not duplicate columns
spe <- compute_qc(spe)
stopifnot(!anyDuplicated(names(colData(spe))))

filtered <- filter_spots(spe, min_counts = 0, min_features = 0, max_mito = 100)$spe
stopifnot(ncol(filtered) == 200)
spe <- normalize_spe(filtered)
stopifnot("logcounts" %in% SummarizedExperiment::assayNames(spe))
# v1.4.1: direct QC metrics and log-normalisation match their definitions
cnt <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
same <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b)))
stopifnot(same(spe$sum, colSums(cnt)), same(spe$detected, colSums(cnt > 0)),
          same(spe$subsets_Mito_percent, 100 * cnt[1, ] / colSums(cnt)))
sf <- colSums(cnt) / mean(colSums(cnt))
stopifnot(isTRUE(all.equal(as.matrix(SummarizedExperiment::assay(spe, "logcounts")),
                           log2(sweep(cnt, 2, sf, "/") + 1), check.attributes = FALSE)))
stopifnot(identical(dim(scale_spots(matrix(1:12, 3, 4), "0-1")), c(3L, 4L)))
spe <- run_dimred(spe, n_pcs = 10)
# HVG-based PCA: HVGs recorded, UMAP off by default
stopifnot(length(S4Vectors::metadata(spe)$hvgs) > 0L,
          !"UMAP" %in% SingleCellExperiment::reducedDimNames(spe))
# v1.4.2: planted gradient and domain genes rank among the most variable genes
stopifnot(all(c(gradient_gene, rownames(spe)[4:8]) %in% head(S4Vectors::metadata(spe)$hvgs, 10)))
spe <- cluster_spatial(spe, k = 5)
stopifnot("STclust" %in% colnames(SummarizedExperiment::colData(spe)))
# kNN-based clustering works at both extremes of spatial weighting
for (w in c(0, 1)) stopifnot(nlevels(cluster_spatial(spe, k = 5, spatial_weight = w)$STclust) >= 1L)

# Default gene selection (HVGs, no full densification); planted gradient must be detected
auto <- compute_moran_geary(spe, k = 4)$results
stopifnot(attr(auto, "alternative") == "greater")
grad <- auto[auto$gene == gradient_gene, ]
stopifnot(nrow(grad) == 4L, all(grad$moran_FDR < 0.05), all(grad$moran_I > 0.3))

# Fix 5: re-filtering after clustering clears stale downstream state
refiltered <- filter_spots(spe, min_counts = 0, min_features = 0, max_mito = 100)$spe
stopifnot(!"logcounts" %in% SummarizedExperiment::assayNames(refiltered),
          is.null(BiocGenerics::sizeFactors(refiltered)),
          !length(SingleCellExperiment::reducedDimNames(refiltered)),
          !"STclust" %in% names(colData(refiltered)))

# Fix 5b: switching normalisation method does not reuse old size factors
sf_scran <- BiocGenerics::sizeFactors(normalize_spe(refiltered, "scran"))
sf_lib   <- BiocGenerics::sizeFactors(normalize_spe(normalize_spe(refiltered, "scran"), "logNormCounts"))
stopifnot(!isTRUE(all.equal(sf_scran, sf_lib)))

# v1.3.0 joint clustering: domains must be shared across sections
joint <- suppressWarnings(cluster_joint(spe, k = 15, spatial_weight = 0.5))
stopifnot(isFALSE(S4Vectors::metadata(joint)$spatial_clustering$sample_specific))
tab <- domain_sample_table(joint)
shared <- rownames(tab)[rowSums(tab >= 10) == 4L]
stopifnot(length(shared) >= 1L)
# Left and right halves of the tissue should mostly land in different domains
maj <- function(v) names(which.max(table(v)))
stopifnot(maj(joint$STclust[left]) != maj(joint$STclust[!left]))
# Per-domain pseudobulk DE now has replicates
de_dom <- suppressWarnings(run_pseudobulk_de(joint, cluster_col = "STclust", min_spots = 10))
stopifnot(nrow(de_dom$results) > 0L, any(de_dom$domain_summary$tested))
# min_spots larger than any sample's domain size -> nothing testable
stopifnot(inherits(try(suppressWarnings(run_pseudobulk_de(joint, cluster_col = "STclust", min_spots = 1000)),
                       silent = TRUE), "try-error"))

spatial <- compute_moran_geary(spe, genes = rownames(spe)[1:2], k = 4)$results
stopifnot(nrow(spatial) == 8L)
stopifnot("symbol" %in% names(spatial), all(spatial$symbol %in% symbols))
# v1.5.0 replicate column and pseudoreplication check
stopifnot(setequal(replicate_candidates(spe), c("sample_id", "animal")))
w <- capture_warnings(de <- run_pseudobulk_de(spe)$results)
stopifnot(nrow(de) > 0L, any(grepl("animal", w)))                 # warns: sections share an animal
expect_error(suppressWarnings(run_pseudobulk_de(spe, replicate_col = "animal")))  # control = 1 animal only
expect_error(run_pseudobulk_de(spe, replicate_col = "STclust"))   # varies within a section

# Windows paths
stopifnot(.is_abs(c("C:\\data", "C:/data", "\\\\server\\share", "/home/x", "~/x")),
          !.is_abs(c("data/raw", "raw")),
          identical(.clean_path(c("C:\\data\\run1\\", "data/raw/", "C:\\", " x ")), c("C:\\data\\run1", "data/raw", "C:\\", "x")))

# v1.6.1 reference level: digits sort first, but "control" must still be the baseline
stopifnot(guess_reference(c("6OHDA", "control")) == "control", guess_reference(c("KO", "WT")) == "WT",
          guess_reference(c("b", "a")) == "a")
ref_c <- suppressWarnings(run_pseudobulk_de(spe))
ref_t <- suppressWarnings(run_pseudobulk_de(spe, reference = "treated"))
stopifnot(ref_c$contrast == "treated vs control", ref_t$contrast == "control vs treated",
          all(ref_c$results$contrast == "treated vs control"))
mm <- match(ref_c$results$gene, ref_t$results$gene)
stopifnot(isTRUE(all.equal(ref_c$results$logFC, -ref_t$results$logFC[mm], tolerance = 1e-4)))  # sign flips
expect_error(suppressWarnings(run_pseudobulk_de(spe, reference = "placebo")))

# v1.6.0 covariates: slide spans both conditions, so it is a covariate (animal is not)
stopifnot(identical(covariate_candidates(spe), "slide"))
res_cov <- suppressWarnings(run_pseudobulk_de(spe, covariates = "slide"))
stopifnot(nrow(res_cov$results) > 0L, identical(res_cov$covariates, "slide"))
expect_error(suppressWarnings(run_pseudobulk_de(spe, covariates = "STclust")))      # varies within a section
expect_error(suppressWarnings(run_pseudobulk_de(spe, covariates = "condition")))    # cannot adjust for the contrast
# Per-replicate plot builds and shows one point per replicate
p <- plot_gene_by_replicate(res_cov, gradient_gene, "all")
stopifnot(inherits(p, "ggplot"), nrow(ggplot2::ggplot_build(p)$data[[2]]) == 4L)
# v1.4.0 sample-sheet import (fake reader: tests merging/validation without Space Ranger files)
fake_reader <- function(genes) function(path, sample_id, ...) {
  m <- matrix(rpois(length(genes) * 20, 5), nrow = length(genes),
              dimnames = list(genes, paste0("AAAC", 1:20, "-1")))  # same barcodes in every slide
  SpatialExperiment(assays = list(counts = m),
                    spatialCoords = cbind(x = 1:20, y = rep(1:4, 5)),
                    rowData = S4Vectors::DataFrame(symbol = genes),
                    colData = S4Vectors::DataFrame(sample_id = rep(sample_id, 20), row.names = colnames(m)))
}
dirs <- file.path(tempdir(), c("A1", "A2", "B1", "B2")); for (d in dirs) dir.create(d, showWarnings = FALSE)
sheet <- data.frame(path = basename(dirs), sample_id = c("A1", "A2", "B1", "B2"),
                    condition = c("control", "control", "treated", "treated"), donor = c("d1", "d2", "d3", "d4"))
genes <- paste0("G", 1:50)
merged <- import_sample_sheet(sheet, base_dir = tempdir(), reader = fake_reader(genes))
stopifnot(ncol(merged) == 80L, !anyDuplicated(colnames(merged)),
          setequal(unique(merged$sample_id), sheet$sample_id),
          all(c("condition", "donor") %in% names(colData(merged))))
rs <- replicate_summary(merged)
stopifnot(all(rs$sections == 2L), all(rs$replicates == 2L), all(rs$de_ready), all(merged$platform == "Visium"))
# Validation failures must stop with an error
expect_error(validate_sample_sheet(transform(sheet, sample_id = c("A1", "A1", "B1", "B2"))))   # duplicate id
expect_error(validate_sample_sheet(sheet[, c("path", "sample_id")]))                            # no condition
expect_error(import_sample_sheet(transform(sheet, path = "no_such_dir"), base_dir = tempdir(),
                                 reader = fake_reader(genes)))                                  # missing folder
mixed_reader <- function(path, sample_id, ...)
  fake_reader(if (sample_id == "A1") genes else paste0("X", 1:50))(path, sample_id)
expect_error(import_sample_sheet(sheet, base_dir = tempdir(), reader = mixed_reader))           # gene mismatch

# v1.5.0 platforms: names normalised, unknown or mixed platforms refused
stopifnot(all(validate_sample_sheet(sheet)$platform == "Visium"),
          all(validate_sample_sheet(transform(sheet, platform = "xenium"))$platform == "Xenium"))
expect_error(validate_sample_sheet(transform(sheet, platform = "SlideSeq")))
expect_error(validate_sample_sheet(transform(sheet, platform = c("Visium", "Xenium", "Visium", "Visium"))))
expect_error(validate_sample_sheet(transform(sheet, platform = c("Xenium", "MERFISH", "Xenium", "Xenium"))))

# Imaging import: control probes removed, no mito warning (panels have no MT genes)
img_reader <- function(path, sample_id, ...) fake_reader(c(genes, "NegControlProbe_1", "BLANK_0001"))(path, sample_id)
img <- import_sample_sheet(transform(sheet, platform = "Xenium"), base_dir = tempdir(), reader = img_reader)
stopifnot(!any(c("NegControlProbe_1", "BLANK_0001") %in% rownames(img)), nrow(img) == 50L, is_imaging(img))
stopifnot(!any(grepl("mitochondrial", capture_warnings(compute_qc(img)))))

# GEO-style Visium files: real 10x files written to disk, then read back
if (requireNamespace("DropletUtils", quietly = TRUE)) {
  geo <- file.path(tempdir(), "geo_visium"); unlink(geo, recursive = TRUE); dir.create(geo)
  bc <- paste0(c("AAAC", "AAAG", "AACC", "AACG", "AAGC"), "TTTTGGGG-1")
  write_sample <- function(prefix, dir_name) {
    cm <- Matrix::Matrix(matrix(rpois(10 * 5, 5), 10, 5), sparse = TRUE)
    DropletUtils::write10xCounts(file.path(geo, dir_name), cm, barcodes = bc,
                                 gene.id = paste0("ENSG", 1:10), gene.symbol = c("MT-ND1", paste0("G", 2:10)),
                                 version = "3")
    pos <- data.frame(bc = c(bc, "GGGGTTTTGGGG-1"), in_tissue = c(1, 1, 1, 1, 0, 1), array_row = 1:6,
                      array_col = 1:6, pxl_row = 101:106, pxl_col = 201:206)
    utils::write.table(pos, file.path(geo, paste0(prefix, "tissue_positions_list.csv")), sep = ",",
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  write_sample("GSM1_A_", "GSM1_A_filtered_feature_bc_matrix")
  write_sample("GSM2_B_", "GSM2_B_filtered_feature_bc_matrix")
  # Two samples in one flat GEO folder: auto-detection must refuse to guess...
  expect_error(read_visium_files(geo, "A"))
  # ...and explicit file names work
  vf <- read_visium_files(geo, "A", counts_file = "GSM1_A_filtered_feature_bc_matrix",
                          positions_file = "GSM1_A_tissue_positions_list.csv")
  stopifnot(ncol(vf) == 4L, nrow(vf) == 10L,                         # out-of-tissue spot dropped
            same(SpatialExperiment::spatialCoords(vf)[1, ], c(201, 101)),
            length(detect_mito_genes(vf)) == 1L)                      # symbols read from features file
  geo_sheet <- data.frame(path = geo, sample_id = c("A", "B"), condition = c("control", "treated"),
                          platform = "VisiumFiles",
                          counts_file = c("GSM1_A_filtered_feature_bc_matrix", "GSM2_B_filtered_feature_bc_matrix"),
                          positions_file = c("GSM1_A_tissue_positions_list.csv", "GSM2_B_tissue_positions_list.csv"))
  geo_spe <- import_sample_sheet(geo_sheet)
  stopifnot(ncol(geo_spe) == 8L, all(geo_spe$platform == "VisiumFiles"),
            !any(c("counts_file", "positions_file") %in% names(colData(geo_spe))))
}

cat("Core synthetic-data smoke tests passed. Real-data validation still required.\n")
