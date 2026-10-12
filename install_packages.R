# One-time setup: installs everything spatialTXsuite needs.
# Windows and macOS get pre-built binaries, so Rtools / Xcode are not required.
options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("BiocManager", "shiny", "ggplot2", "Matrix", "igraph", "harmony")
bioc <- c("SpatialExperiment", "scater", "scuttle", "scran", "edgeR", "spdep",
          "BiocNeighbors", "DropletUtils", "SpatialExperimentIO")
need <- setdiff(cran, rownames(installed.packages()))
if (length(need)) install.packages(need)
need <- setdiff(bioc, rownames(installed.packages()))
if (length(need)) BiocManager::install(need, ask = FALSE, update = FALSE)
missing <- setdiff(c(cran, bioc), rownames(installed.packages()))
if (length(missing)) {
  message("Could not install: ", paste(missing, collapse = ", "))
  quit(status = 1)
}
message("All packages installed. Run the tests next (run_tests.bat on Windows).")
