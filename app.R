# spatialTXsuite revised scientific-core application (v1.6.1).
# Use Space Ranger output directory on the SERVER (not client-side folder upload).
#
# v1.2.1 changes
#  - Upload limit raised from Shiny's 5 MB default (real SPE RDS files are large).
#  - Re-filtering / re-normalising clears stale logcounts, PCA, clusters and results
#    (logic lives in clear_downstream(), R/02_qc.R).
#  - Progress indicators on long steps; fixed seed for reproducible Louvain clustering.
#  - v1.3.0: joint cross-sample clustering mode (shared domains) and min-spots rule for DE.
#  - v1.4.0: multi-sample import from a sample sheet (R/01_import.R).
#  - v1.5.0: Xenium/MERFISH/CosMx and GEO-style Visium import; replicate column for DE.
#  - v1.6.0: covariates (e.g. slide) in the DE model; per-replicate gene expression plot.
#  - v1.6.1: reference (baseline) condition chosen explicitly; logFC direction labelled.
#  - QC tab reports how many mito genes were found, so a disabled mito filter is visible.
library(shiny)
library(SpatialExperiment)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
for (f in c("01_import.R", "02_qc.R", "03_normalisation.R", "04_clustering.R", "05_spatial_stats.R", "07_diff_analysis.R"))
  source(file.path("R", f), local = FALSE)

options(shiny.maxRequestSize = as.numeric(Sys.getenv("SPATIALTX_MAX_UPLOAD_MB", "4000")) * 1024^2)

ui <- fluidPage(
  titlePanel("spatialTXsuite: spatial transcriptomics analysis"),
  p("Analyzes actual SpatialExperiment data. No simulated analyses. Additional experimental modules are not presented as validated."),
  tabsetPanel(
    tabPanel("Import",
      radioButtons("import_mode", "Import from",
        c("Sample sheet (several samples)" = "sheet",
          "Single sample folder" = "single",
          "SpatialExperiment .rds" = "rds")),
      conditionalPanel("input.import_mode == 'sheet'",
        fileInput("sheet", "Sample sheet CSV", accept = ".csv"),
        textInput("base_dir", "Base folder for relative paths (on the server, optional)", ""),
        helpText("One row per tissue section. Required columns: path, sample_id, condition. Optional: platform ",
                 "(Visium, VisiumFiles, Xenium, MERFISH, CosMx), counts_file and positions_file (VisiumFiles), ",
                 "and any metadata such as animal or donor. Paths are folders on the machine running the app.")),
      conditionalPanel("input.import_mode == 'single'",
        selectInput("platform", "Platform",
          c("Visium (Space Ranger folder)" = "Visium",
            "Visium files (GEO: .h5 or matrix.mtx + tissue_positions)" = "VisiumFiles",
            "Xenium" = "Xenium", "MERFISH (Vizgen)" = "MERFISH", "CosMx" = "CosMx")),
        textInput("path", "Sample output folder on the server"),
        textInput("sample", "Sample ID", "S1"),
        textInput("condition", "Condition label", "control")),
      conditionalPanel("input.import_mode == 'rds'",
        fileInput("rds", "SpatialExperiment .rds", accept = ".rds")),
      actionButton("load", "Import dataset"),
      verbatimTextOutput("status"),
      h5("Spots or cells per section"), tableOutput("samples"),
      h5("Sections and replicates per condition"), tableOutput("replicates"),
      helpText("Replicates are counted with the replicate column chosen on the Differential expression tab.")),
    tabPanel("QC",
      numericInput("mincounts", "Minimum UMIs", 500, min = 0),
      numericInput("minfeatures", "Minimum genes", 250, min = 0),
      numericInput("maxmito", "Maximum mitochondrial percent", 20, min = 0, max = 100),
      actionButton("qc", "Filter spots"),
      helpText("Filtering discards any earlier normalization, clustering and statistics."),
      verbatimTextOutput("qc_status"), plotOutput("qc_plot")),
    tabPanel("Normalization",
      selectInput("method", "Method", c("logNormCounts", "scran")),
      actionButton("norm", "Normalize"), verbatimTextOutput("norm_status")),
    tabPanel("Clusters",
      radioButtons("cluster_mode", "Mode",
        c("Joint across samples (shared domains, needed for per-domain DE)" = "joint",
          "Per sample (labels not comparable between sections)" = "per_sample")),
      numericInput("clusters_k", "Neighbors", 15, min = 1),
      sliderInput("weight", "Spatial weighting", 0, 1, 0.5),
      numericInput("seed", "Random seed", 1, min = 1),
      actionButton("cluster", "Run spatial clustering"),
      verbatimTextOutput("cluster_status"), plotOutput("clusters_plot"),
      h5("Spots per domain and sample"), tableOutput("domain_table")),
    tabPanel("Spatial statistics",
      numericInput("moran_k", "Spatial neighbors", 6, min = 1),
      actionButton("spatial", "Calculate Moran's I and Geary's C"),
      verbatimTextOutput("spatial_status"), tableOutput("spatial_table")),
    tabPanel("Differential expression",
      helpText("Requires at least two independent biological replicates PER condition. Import several samples with a sample sheet (Import tab), then use joint clustering for per-domain DE."),
      selectInput("replicate_col", "Biological replicate (unit of replication)", choices = "sample_id"),
      helpText("If several sections come from the same animal or donor, choose the animal/donor column: ",
               "its sections are pooled into one replicate. Treating them as separate replicates inflates significance."),
      selectInput("domain", "Stratify", c("All spots" = "all", "Spatial domain" = "STclust")),
      selectInput("reference", "Reference (baseline) condition", choices = NULL),
      helpText("logFC is reported as the other condition relative to this one (positive = higher than reference)."),
      selectizeInput("covariates", "Adjust for (technical covariates, optional)", choices = NULL, multiple = TRUE),
      helpText("Factors present in both conditions, such as Visium slide or batch. Each costs one or more ",
               "degrees of freedom, so use them only with enough replicates."),
      numericInput("min_spots", "Minimum spots/cells per replicate in a domain", 10, min = 1),
      actionButton("de", "Run pseudobulk edgeR"),
      verbatimTextOutput("de_status"), h5("Replicates used per domain"), tableOutput("de_summary"),
      tableOutput("de_table"), plotOutput("volcano"),
      h4("Expression per replicate"),
      fluidRow(column(6, selectizeInput("plot_gene", "Gene", choices = NULL)),
               column(6, selectInput("plot_domain", "Domain", choices = NULL))),
      plotOutput("gene_plot")),
    tabPanel("Export",
      h4("Save the processed object"),
      helpText("Recommended for large datasets: writes the .rds directly to a folder on the machine running the app. ",
               "Browser downloads of multi-GB objects can stall, especially in the RStudio viewer pane."),
      textInput("save_path", "File path", file.path(getwd(), "results", "spatialTXsuite_processed.rds"), width = "100%"),
      actionButton("save_disk", "Save to folder"), verbatimTextOutput("save_status"),
      hr(), h4("Downloads"),
      downloadButton("save", "Download processed SpatialExperiment RDS"),
      downloadButton("dl_spatial", "Download spatial statistics CSV"),
      downloadButton("dl_de", "Download DE results CSV"))
  ))

server <- function(input, output, session) {
  rv <- reactiveValues(spe = NULL, spatial = NULL, de = NULL, de_summary = NULL, de_full = NULL)
  status <- reactiveVal("No data loaded")
  output$status <- renderText(status())

  notify_error <- function(e) showNotification(conditionMessage(e), type = "error", duration = NULL)
  # Surface warnings (e.g. "no mito genes") instead of losing them in the server log
  with_warnings <- function(expr) withCallingHandlers(expr, warning = function(w) {
    showNotification(conditionMessage(w), type = "warning", duration = 15)
    invokeRestart("muffleWarning")
  })

  observeEvent(input$load, {
    tryCatch(withProgress(message = "Importing dataset", value = 0.2, {
      spe <- switch(input$import_mode,
        rds = {
          if (is.null(input$rds)) stop("Choose an .rds file")
          x <- readRDS(input$rds$datapath)
          if (!inherits(x, "SpatialExperiment")) stop("RDS must contain SpatialExperiment")
          if ("counts" %in% assayNames(x)) assay(x, "counts", withDimnames = FALSE) <- to_dgc(assay(x, "counts"))
          x
        },
        single = {
          path <- .clean_path(input$path)   # Windows-safe: trailing \ removed
          if (!nzchar(path) || !dir.exists(path)) stop("Enter a valid folder path on the Shiny server")
          sheet <- data.frame(path = path, sample_id = trimws(input$sample), condition = trimws(input$condition),
                              platform = input$platform)
          with_warnings(import_sample_sheet(sheet))
        },
        sheet = {
          if (is.null(input$sheet)) stop("Choose a sample sheet CSV")
          sheet <- read_sample_sheet(input$sheet$datapath)
          with_warnings(import_sample_sheet(sheet, base_dir = trimws(input$base_dir),
            progress = function(i, n, id) setProgress(0.1 + 0.6 * i / n, detail = paste("Reading", id, "-", i, "of", n))))
        })
      if (!"sample_id" %in% names(colData(spe))) stop("Missing sample_id metadata")
      if (!"counts" %in% assayNames(spe)) stop("Missing raw counts assay")
      if (ncol(SpatialExperiment::spatialCoords(spe)) < 2L) stop("Missing spatial coordinates")
      incProgress(0.5, detail = "Computing QC metrics")
      spe <- with_warnings(compute_qc(spe))
      rv$spe <- spe; rv$spatial <- NULL; rv$de <- NULL; rv$de_summary <- NULL; rv$de_full <- NULL
      plat <- if ("platform" %in% names(colData(spe))) unique(as.character(spe$platform)) else "Visium"
      # platform-appropriate QC defaults (imaging panels: low counts per cell, no mito genes)
      th <- get_qc_thresholds(if (length(plat) == 1L && plat != "VisiumFiles") plat else "Visium")
      updateNumericInput(session, "mincounts", value = th$min_counts)
      updateNumericInput(session, "minfeatures", value = th$min_features)
      updateNumericInput(session, "maxmito", value = if (is.na(th$max_mito)) 100 else th$max_mito)
      cands <- if ("condition" %in% names(colData(spe))) replicate_candidates(spe) else "sample_id"
      updateSelectInput(session, "replicate_col", choices = cands, selected = "sample_id")
      if ("condition" %in% names(colData(spe))) {
        conds <- sort(unique(as.character(spe$condition)))
        updateSelectInput(session, "reference", choices = conds, selected = guess_reference(conds))
      }
      updateSelectizeInput(session, "covariates",
        choices = if ("condition" %in% names(colData(spe))) covariate_candidates(spe) else character(0), selected = character(0))
      status(paste0("Loaded ", nrow(spe), " genes by ", ncol(spe), " spots/cells across ",
                    length(unique(spe$sample_id)), " sections (", paste(plat, collapse = ", "), ")"))
    }), error = function(e) status(paste("IMPORT ERROR:", conditionMessage(e))))
  })

  output$samples <- renderTable({ req(rv$spe); as.data.frame(table(sample = rv$spe$sample_id)) })
  output$replicates <- renderTable({
    req(rv$spe)
    rc <- if (!is.null(input$replicate_col) && input$replicate_col %in% names(colData(rv$spe))) input$replicate_col else "sample_id"
    replicate_summary(rv$spe, replicate_col = rc)
  })
  output$qc_plot <- renderPlot({
    req(rv$spe)
    hist(rv$spe$sum, breaks = 40, main = "Raw counts per spot", xlab = "UMIs")
  })
  output$qc_status <- renderText({
    req(rv$spe)
    n_mito <- S4Vectors::metadata(rv$spe)$qc_params$mito_genes
    paste0(ncol(rv$spe), " retained spots\n",
           if (isTRUE(n_mito > 0)) paste(n_mito, "mitochondrial genes detected")
           else "No mitochondrial genes detected: mito filter is NOT applied")
  })

  observeEvent(input$qc, {
    req(rv$spe)
    tryCatch(withProgress(message = "Filtering spots", {
      z <- with_warnings(filter_spots(rv$spe, min_counts = input$mincounts,
                                      min_features = input$minfeatures, max_mito = input$maxmito))
      rv$spe <- z$spe; rv$spatial <- NULL; rv$de <- NULL; rv$de_summary <- NULL; rv$de_full <- NULL
      showNotification(paste(sum(z$audit$discard), "spots excluded; downstream results reset"), type = "message")
    }), error = notify_error)
  })

  output$norm_status <- renderText({ req(rv$spe); paste("Assays:", paste(assayNames(rv$spe), collapse = ", ")) })
  observeEvent(input$norm, {
    req(rv$spe)
    tryCatch(withProgress(message = "Normalizing", {
      rv$spe <- normalize_spe(rv$spe, input$method); rv$spatial <- NULL; rv$de <- NULL; rv$de_summary <- NULL; rv$de_full <- NULL
    }), error = notify_error)
  })

  observeEvent(input$cluster, {
    req(rv$spe)
    tryCatch(withProgress(message = "Running PCA and spatial clustering", {
      set.seed(input$seed)
      rv$spe <- with_warnings(if (input$cluster_mode == "joint")
        cluster_joint(rv$spe, k = input$clusters_k, spatial_weight = input$weight)
        else cluster_spatial(rv$spe, k = input$clusters_k, spatial_weight = input$weight))
      rv$de <- NULL; rv$de_summary <- NULL; rv$de_full <- NULL
    }), error = notify_error)
  })
  output$cluster_status <- renderText({
    req(rv$spe)
    if (!"STclust" %in% names(colData(rv$spe))) "Clustering not yet run"
    else {
      sc <- S4Vectors::metadata(rv$spe)$spatial_clustering
      paste0(length(unique(rv$spe$STclust)), " domains detected (",
             if (isTRUE(sc$sample_specific)) "per sample" else
               paste0("joint", if (isTRUE(sc$harmony)) ", Harmony-corrected" else ", NOT batch-corrected"), ")")
    }
  })
  output$domain_table <- renderTable({
    req(rv$spe); req("STclust" %in% names(colData(rv$spe)))
    domain_sample_table(rv$spe)
  }, rownames = TRUE)
  output$clusters_plot <- renderPlot({
    req(rv$spe); req("STclust" %in% names(colData(rv$spe)))
    xy <- spatialCoords(rv$spe)
    d <- data.frame(x = xy[, 1], y = xy[, 2], cluster = rv$spe$STclust, sample = rv$spe$sample_id)
    # Shift each section to start at 0 so panels share axes; free scales are not
    # allowed with coord_equal() in current ggplot2, and equal aspect keeps true tissue shape
    d$x <- d$x - ave(d$x, d$sample, FUN = min)
    d$y <- d$y - ave(d$y, d$sample, FUN = min)
    pt <- if (nrow(d) > 50000) 0.05 else if (nrow(d) > 10000) 0.3 else 0.9   # imaging data: many cells
    ggplot(d, aes(x, y, color = cluster)) + geom_point(size = pt) +
      facet_wrap(~sample) + coord_equal() + scale_y_reverse() + theme_minimal() +
      guides(color = guide_legend(override.aes = list(size = 3)))
  })

  observeEvent(input$spatial, {
    req(rv$spe)
    tryCatch(withProgress(message = "Computing Moran's I and Geary's C", {
      rv$spatial <- compute_moran_geary(rv$spe, k = input$moran_k)$results
    }), error = notify_error)
  })
  output$spatial_status <- renderText({
    if (is.null(rv$spatial)) "Spatial statistics not yet run" else paste(nrow(rv$spatial), "gene-sample tests")
  })
  output$spatial_table <- renderTable({ req(rv$spatial); head(rv$spatial[order(rv$spatial$moran_FDR), ], 20) }, digits = 4)

  observeEvent(input$de, {
    req(rv$spe)
    tryCatch(withProgress(message = "Running pseudobulk edgeR", {
      res <- with_warnings(run_pseudobulk_de(rv$spe,
                 cluster_col = if (input$domain == "all") NULL else "STclust",
                 min_spots = input$min_spots, replicate_col = input$replicate_col,
                 covariates = if (length(input$covariates)) input$covariates else character(0),
                 reference = input$reference))
      rv$de <- res$results; rv$de_summary <- res$domain_summary; rv$de_full <- res
      # every gene is searchable (not only DE-tested ones, e.g. Xist for sex checks); DE genes listed first
      all_genes <- rownames(res$pseudobulk[[1]]$counts)
      ranked <- unique(c(res$results$gene[order(res$results$PValue)], all_genes))
      sym <- get_gene_symbols(rv$spe); names(sym) <- rownames(rv$spe)
      lab <- ifelse(is.na(sym[ranked]) | sym[ranked] == ranked, ranked, paste0(sym[ranked], " (", ranked, ")"))
      updateSelectizeInput(session, "plot_gene", choices = stats::setNames(ranked, lab), selected = ranked[1], server = TRUE)
      updateSelectInput(session, "plot_domain", choices = names(res$pseudobulk))
    }), error = notify_error)
  })
  output$de_status <- renderText({
    if (is.null(rv$de)) "DE not yet run"
    else paste0(nrow(rv$de), " gene/domain tests. logFC = ", rv$de_full$contrast,
                " (positive = higher in ", sub(" vs .*", "", rv$de_full$contrast), ")")
  })
  output$de_summary <- renderTable({ req(rv$de_summary); rv$de_summary })
  output$gene_plot <- renderPlot({
    req(rv$de_full, input$plot_gene, input$plot_domain)
    req(input$plot_domain %in% names(rv$de_full$pseudobulk))
    sym <- get_gene_symbols(rv$spe)[match(input$plot_gene, rownames(rv$spe))]
    plot_gene_by_replicate(rv$de_full, input$plot_gene, input$plot_domain,
                           label = if (is.na(sym)) input$plot_gene else sym)
  })
  output$de_table <- renderTable({ req(rv$de); head(rv$de[order(rv$de$FDR), ], 20) }, digits = 4)
  output$volcano <- renderPlot({ req(rv$de); plot_volcano(rv$de) })

  save_msg <- reactiveVal("")
  output$save_status <- renderText(save_msg())
  observeEvent(input$save_disk, {
    req(rv$spe)
    path <- path.expand(trimws(input$save_path))
    tryCatch(withProgress(message = "Saving .rds (large objects take a few minutes)", value = 0.3, {
      if (!grepl("\\.rds$", path, ignore.case = TRUE)) path <- paste0(path, ".rds")
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      tmp <- paste0(path, ".partial")
      saveRDS(rv$spe, tmp)                 # write then rename: no half-written file on failure
      file.rename(tmp, path)
      save_msg(sprintf("Saved %s (%.0f MB)", path, file.size(path) / 1e6))
    }), error = function(e) save_msg(paste("SAVE ERROR:", conditionMessage(e))))
  })
  output$save <- downloadHandler(filename = "spatialTXsuite_processed.rds",
    content = function(file) { req(rv$spe); saveRDS(rv$spe, file) })
  output$dl_spatial <- downloadHandler(filename = "spatial_stats.csv",
    content = function(file) { req(rv$spatial); write.csv(rv$spatial, file, row.names = FALSE) })
  output$dl_de <- downloadHandler(filename = "pseudobulk_DE.csv",
    content = function(file) { req(rv$de); write.csv(rv$de, file, row.names = FALSE) })
}
shinyApp(ui, server)
