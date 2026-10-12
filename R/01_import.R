# Core import (v1.5.0): several spatial platforms, one sample or many via a sample sheet.
#
# Supported platforms (sample sheet column `platform`, default "Visium"):
#   Visium       Space Ranger output folder (run folder or its outs/)
#   VisiumFiles  GEO-style Visium files: a count matrix (.h5, or matrix.mtx + features + barcodes)
#                plus a tissue_positions(_list).csv[.gz]. Found automatically in the folder, or
#                named in optional sheet columns counts_file / positions_file (for GEO folders
#                that hold many samples with prefixed file names).
#   Xenium       Xenium Onboard Analysis output folder      (needs SpatialExperimentIO)
#   MERFISH      Vizgen MERSCOPE output folder               (needs SpatialExperimentIO)
#   CosMx        NanoString CosMx flat-file export folder    (needs SpatialExperimentIO)
# Anything else can be loaded as a SpatialExperiment .rds.
#
# Sample sheet (CSV), one row per tissue section:
#   path,sample_id,condition[,platform][,counts_file,positions_file][,any extra columns e.g. animal,donor]
# Extra columns are copied onto every spot/cell of that section. If several sections come from
# the same animal/donor, add that column and choose it as the replicate column in the DE tab.

PLATFORMS <- c("Visium", "VisiumFiles", "Xenium", "MERFISH", "CosMx")
PLATFORM_FAMILY <- c(Visium = "sequencing", VisiumFiles = "sequencing",
                     Xenium = "imaging", MERFISH = "imaging", CosMx = "imaging")
SHEET_REQUIRED <- c("path", "sample_id", "condition")
SHEET_RESERVED <- c("path", "sample_id", "platform", "counts_file", "positions_file")

# Negative-control and unassigned probes in imaging panels: not genes, removed on import
CONTROL_PATTERN <- "^(NegControl|NegPrb|Negative|SystemControl|BLANK|Blank-|UnassignedCodeword|DeprecatedCodeword|Intergenic|antisense_)"

is_imaging <- function(spe) {
  "platform" %in% names(SummarizedExperiment::colData(spe)) &&
    all(as.character(spe$platform) %in% names(PLATFORM_FAMILY)[PLATFORM_FAMILY == "imaging"])
}

normalise_platform <- function(x) {
  x <- trimws(as.character(x)); x[is.na(x) | !nzchar(x)] <- "Visium"
  hit <- match(tolower(x), tolower(PLATFORMS))
  if (anyNA(hit)) stop("Unknown platform(s): ", paste(unique(x[is.na(hit)]), collapse = ", "),
                       ". Use one of: ", paste(PLATFORMS, collapse = ", "))
  PLATFORMS[hit]
}

read_sample_sheet <- function(file) {
  if (!file.exists(file)) stop("Sample sheet not found: ", file)
  sheet <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE,
                           strip.white = TRUE, na.strings = c("", "NA"))
  names(sheet) <- tolower(trimws(names(sheet)))
  validate_sample_sheet(sheet)
}

validate_sample_sheet <- function(sheet) {
  miss <- setdiff(SHEET_REQUIRED, names(sheet))
  if (length(miss)) stop("Sample sheet is missing column(s): ", paste(miss, collapse = ", "))
  if (!nrow(sheet)) stop("Sample sheet has no rows")
  for (col in SHEET_REQUIRED) {
    bad <- which(is.na(sheet[[col]]) | !nzchar(trimws(as.character(sheet[[col]]))))
    if (length(bad)) stop("Empty '", col, "' in sample sheet row(s): ", paste(bad, collapse = ", "))
    sheet[[col]] <- trimws(as.character(sheet[[col]]))
  }
  dup <- unique(sheet$sample_id[duplicated(sheet$sample_id)])
  if (length(dup)) stop("Duplicate sample_id(s): ", paste(dup, collapse = ", "),
                        ". Each row must be a distinct tissue section.")
  if (!all(grepl("^[A-Za-z0-9._-]+$", sheet$sample_id)))
    stop("sample_id may contain only letters, digits, '.', '_' and '-'")
  sheet$platform <- normalise_platform(if ("platform" %in% names(sheet)) sheet$platform else rep("Visium", nrow(sheet)))
  fam <- unique(PLATFORM_FAMILY[sheet$platform])
  if (length(fam) > 1L)
    stop("Sequencing-based (Visium) and imaging-based (Xenium/MERFISH/CosMx) samples cannot be merged: ",
         "they measure different gene sets at different resolutions. Analyse them separately.")
  if (length(unique(sheet$platform[PLATFORM_FAMILY[sheet$platform] == "imaging"])) > 1L)
    stop("Different imaging platforms use different gene panels; import one platform at a time")
  sheet
}

# ---- helpers -------------------------------------------------------------------------

to_dgc <- function(m) {
  if (methods::is(m, "dgCMatrix")) return(m)
  if (is.matrix(m)) return(methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix"))
  tryCatch(methods::as(m, "dgCMatrix"),
           error = function(e) methods::as(methods::as(m, "sparseMatrix"), "CsparseMatrix"))
}

drop_control_probes <- function(spe) {
  ctrl <- grepl(CONTROL_PATTERN, get_gene_symbols(spe)) | grepl(CONTROL_PATTERN, rownames(spe))  # get_gene_symbols: 02_qc.R
  if (any(ctrl)) {
    S4Vectors::metadata(spe)$control_probes_removed <- rownames(spe)[ctrl]
    spe <- spe[!ctrl, ]
  }
  spe
}

# absolute: /unix, ~home, C:\ or C:/ drive, \\server\share UNC (Windows network path)
.is_abs <- function(p) grepl("^(/|~|[A-Za-z]:|\\\\)", p, perl = TRUE)
# Windows: dir.exists("C:/data/run1/") is FALSE with a trailing separator, so strip it
.clean_path <- function(p) {
  p <- trimws(p)
  keep <- grepl("^[A-Za-z]:[/\\\\]?$", p, perl = TRUE) | p %in% c("/", "\\")   # leave drive roots alone
  p[!keep] <- sub("[/\\\\]+$", "", p[!keep], perl = TRUE)
  p
}


.find_one <- function(dir, pattern, what, prefer = NULL) {
  f <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(f) > 1L && !is.null(prefer)) {
    p <- f[grepl(prefer, f, ignore.case = TRUE)]
    if (length(p)) f <- p
  }
  if (length(f) > 1L) stop("Several ", what, " files found in ", dir, ": ", paste(basename(f), collapse = ", "),
                           ". Name the files in sample sheet columns counts_file / positions_file.")
  if (length(f)) f else NULL
}

# ---- readers ---------------------------------------------------------------------------

# Space Ranger folder. Accepts the run folder or its outs/ folder.
read_one_visium <- function(path, sample_id) {
  attempt <- function(p) tryCatch(
    SpatialExperiment::read10xVisium(samples = p, sample_id = sample_id, type = "sparse",
                                     data = "filtered", images = NULL, load = TRUE),
    error = function(e) e)
  spe <- attempt(path)
  if (inherits(spe, "error")) {
    alt <- if (basename(normalizePath(path, mustWork = FALSE)) == "outs") dirname(path)
           else if (dir.exists(file.path(path, "outs"))) file.path(path, "outs") else NULL
    first <- conditionMessage(spe)
    if (!is.null(alt)) spe <- attempt(alt)
    if (inherits(spe, "error")) stop("Could not read Visium data for ", sample_id, " at ", path, ": ", first,
                                     ". If this is a GEO download without the Space Ranger folder layout, use platform VisiumFiles.")
  }
  spe
}

read_positions <- function(file) {
  first <- readLines(file, n = 1L)
  pos <- utils::read.csv(file, header = grepl("barcode", first, ignore.case = TRUE), stringsAsFactors = FALSE)
  if (ncol(pos) < 6L) stop("Positions file should have 6 columns (barcode, in_tissue, array_row, array_col, ",
                           "pxl_row_in_fullres, pxl_col_in_fullres): ", basename(file))
  pos <- pos[, 1:6]
  names(pos) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  pos$barcode <- as.character(pos$barcode)
  pos
}

# GEO-style Visium: count matrix + positions table, no Space Ranger folder layout needed
read_visium_files <- function(path, sample_id, counts_file = NULL, positions_file = NULL) {
  resolve <- function(f) if (is.null(f) || is.na(f) || !nzchar(f)) NULL else if (.is_abs(f)) path.expand(f) else file.path(path, f)
  cf <- resolve(counts_file); pf <- resolve(positions_file)
  if (is.null(cf)) {
    cf <- .find_one(path, "\\.h5$", "count matrix (.h5)", prefer = "filtered")
    if (is.null(cf)) {
      mtx <- .find_one(path, "matrix\\.mtx(\\.gz)?$", "matrix.mtx", prefer = "filtered")
      if (!is.null(mtx)) cf <- dirname(mtx)
    }
  }
  if (is.null(cf) || !file.exists(cf)) stop("No count matrix (.h5 or matrix.mtx) found for ", sample_id, " in ", path)
  if (is.null(pf)) pf <- .find_one(path, "tissue_positions(_list)?\\.csv(\\.gz)?$", "tissue positions")
  if (is.null(pf) || !file.exists(pf)) stop("No tissue_positions file found for ", sample_id, " in ", path)

  sce <- DropletUtils::read10xCounts(cf, col.names = TRUE)
  rd <- SummarizedExperiment::rowData(sce)
  if ("Type" %in% names(rd)) sce <- sce[rd$Type == "Gene Expression", ]  # drop antibody/other features
  pos <- read_positions(pf)
  pos <- pos[pos$in_tissue == 1, , drop = FALSE]
  keep <- intersect(colnames(sce), pos$barcode)
  if (!length(keep)) {   # some exports drop or add the "-1" suffix on one side only
    strip <- function(b) sub("-[0-9]+$", "", b)
    m <- match(strip(colnames(sce)), strip(pos$barcode))
    pos$barcode[m[!is.na(m)]] <- colnames(sce)[!is.na(m)]
    keep <- intersect(colnames(sce), pos$barcode)
  }
  if (!length(keep)) stop("No barcodes in the count matrix match in-tissue spots in ", basename(pf))
  pos <- pos[match(keep, pos$barcode), ]
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = to_dgc(SingleCellExperiment::counts(sce)[, keep])),
    rowData = SummarizedExperiment::rowData(sce),
    colData = S4Vectors::DataFrame(barcode = keep, in_tissue = TRUE, array_row = pos$array_row,
                                   array_col = pos$array_col, sample_id = sample_id, row.names = keep),
    spatialCoords = cbind(pxl_col_in_fullres = as.numeric(pos$pxl_col_in_fullres),
                          pxl_row_in_fullres = as.numeric(pos$pxl_row_in_fullres)))
}

# Imaging platforms via SpatialExperimentIO (cell centroids as coordinates)
read_imaging <- function(path, sample_id, platform) {
  if (!requireNamespace("SpatialExperimentIO", quietly = TRUE))
    stop(platform, " import needs SpatialExperimentIO: BiocManager::install('SpatialExperimentIO')")
  spe <- switch(platform,
    Xenium  = SpatialExperimentIO::readXeniumSXE(dirName = path, returnType = "SPE"),
    MERFISH = SpatialExperimentIO::readMerfishSXE(dirName = path, returnType = "SPE"),
    CosMx   = SpatialExperimentIO::readCosmxSXE(dirName = path, returnType = "SPE"))
  if (!inherits(spe, "SpatialExperiment")) stop(platform, " reader did not return a SpatialExperiment for ", sample_id)
  spe
}

# Default reader used by import_sample_sheet(): dispatches on the sheet row's platform
read_spatial_sample <- function(path, sample_id, info = NULL) {
  platform <- if (!is.null(info) && "platform" %in% names(info)) info$platform[1] else "Visium"
  get <- function(col) if (!is.null(info) && col %in% names(info)) info[[col]][1] else NULL
  switch(platform,
    Visium      = read_one_visium(path, sample_id),
    VisiumFiles = read_visium_files(path, sample_id, get("counts_file"), get("positions_file")),
    read_imaging(path, sample_id, platform))
}

# Common clean-up for every reader's output
standardise_sample <- function(s, sample_id, platform) {
  if (!inherits(s, "SpatialExperiment")) stop("Reader did not return a SpatialExperiment for ", sample_id)
  if (!ncol(s)) stop("No spots/cells read for sample ", sample_id)
  if (!"counts" %in% SummarizedExperiment::assayNames(s)) stop("No counts assay for sample ", sample_id)
  if (ncol(SpatialExperiment::spatialCoords(s)) < 2L) stop("No spatial coordinates for sample ", sample_id)
  SummarizedExperiment::assay(s, "counts", withDimnames = FALSE) <- to_dgc(SummarizedExperiment::assay(s, "counts"))
  # keep only the counts assay: other assays differ between readers and would block cbind
  SummarizedExperiment::assays(s) <- SummarizedExperiment::assays(s)["counts"]
  SingleCellExperiment::altExps(s) <- list()
  SingleCellExperiment::reducedDims(s) <- list()
  s$sample_id <- rep(sample_id, ncol(s))
  s$platform <- rep(platform, ncol(s))
  if (PLATFORM_FAMILY[[platform]] == "imaging") s <- drop_control_probes(s)
  s
}

# ---- import -----------------------------------------------------------------------------

# reader is injectable so the merge/validation logic can be tested without real files
import_sample_sheet <- function(sheet, base_dir = NULL, reader = read_spatial_sample,
                                min_gene_overlap = 0.8, progress = NULL) {
  sheet <- validate_sample_sheet(sheet)
  paths <- .clean_path(sheet$path)
  if (!is.null(base_dir) && nzchar(base_dir)) {
    base_dir <- .clean_path(base_dir)
    rel <- !.is_abs(paths)
    paths[rel] <- file.path(base_dir, paths[rel])
  }
  paths <- path.expand(paths)
  missing <- !dir.exists(paths)
  if (any(missing)) stop("Folder(s) not found on the server: ", paste(paths[missing], collapse = "; "))

  extra <- setdiff(names(sheet), SHEET_RESERVED)
  spes <- vector("list", nrow(sheet))
  for (i in seq_len(nrow(sheet))) {
    if (is.function(progress)) progress(i, nrow(sheet), sheet$sample_id[i])
    s <- reader(paths[i], sheet$sample_id[i], sheet[i, , drop = FALSE])
    s <- standardise_sample(s, sheet$sample_id[i], sheet$platform[i])
    for (col in extra) s[[col]] <- rep(sheet[[col]][i], ncol(s))
    colnames(s) <- paste(sheet$sample_id[i], colnames(s), sep = "_")  # barcodes/cell IDs repeat across slides
    spes[[i]] <- s
  }

  if (length(spes) == 1L) {
    S4Vectors::metadata(spes[[1]])$sample_sheet <- sheet
    return(spes[[1]])
  }
  gene_sets <- lapply(spes, rownames)
  common <- Reduce(intersect, gene_sets)
  frac <- length(common) / max(lengths(gene_sets))
  if (frac < min_gene_overlap)
    stop(sprintf("Only %.0f%% of genes are shared across samples (%d genes). Were these run with the same reference or panel?",
                 100 * frac, length(common)))
  if (length(common) < max(lengths(gene_sets)))
    warning(max(lengths(gene_sets)) - length(common), " genes not present in every sample were dropped")
  spes <- lapply(spes, function(s) s[common, ])
  rd_cols <- Reduce(intersect, lapply(spes, function(s) names(SummarizedExperiment::rowData(s))))
  spes <- lapply(spes, function(s) {
    SummarizedExperiment::rowData(s) <- SummarizedExperiment::rowData(spes[[1]])[, rd_cols, drop = FALSE]
    s
  })
  cd_cols <- Reduce(intersect, lapply(spes, function(s) names(SummarizedExperiment::colData(s))))
  spes <- lapply(spes, function(s) {
    SummarizedExperiment::colData(s) <- SummarizedExperiment::colData(s)[, cd_cols, drop = FALSE]
    s
  })
  sc_names <- colnames(SpatialExperiment::spatialCoords(spes[[1]]))
  spes <- lapply(spes, function(s) {
    xy <- SpatialExperiment::spatialCoords(s); colnames(xy) <- sc_names
    SpatialExperiment::spatialCoords(s) <- xy
    s
  })
  merged <- do.call(cbind, spes)
  S4Vectors::metadata(merged) <- list(sample_sheet = sheet, imported = Sys.time())
  merged
}

# Sections (and, if given, true replicates) per condition: tells the user whether DE is possible
replicate_summary <- function(spe, sample_col = "sample_id", condition_col = "condition", replicate_col = sample_col) {
  cd <- SummarizedExperiment::colData(spe)
  if (!condition_col %in% names(cd)) return(NULL)
  cond <- as.character(cd[[condition_col]])
  conds <- sort(unique(cond))
  out <- data.frame(condition = conds,
                    sections = vapply(conds, function(k) length(unique(cd[[sample_col]][cond == k])), integer(1)),
                    replicates = vapply(conds, function(k) length(unique(cd[[replicate_col]][cond == k])), integer(1)),
                    spots_or_cells = vapply(conds, function(k) sum(cond == k), integer(1)),
                    row.names = NULL)
  out$de_ready <- out$replicates >= 2L
  out
}
