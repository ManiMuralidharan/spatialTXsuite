# spatialTXsuite

An R Shiny application for replicate-aware spatial transcriptomics analysis of SpatialExperiment data: import, QC, normalisation, spatial domains, spatial autocorrelation and pseudobulk differential expression.

The synthetic-data test suite passes, and the full workflow has been run end to end on a public 24-section mouse Visium dataset (GSE298680). The Xenium, MERFISH and CosMx readers have so far been tested only on simulated input, not real instrument exports. Check results against known biology before relying on them.

## Supported input

| Platform | What to point at | Notes |
|---|---|---|
| `Visium` | Space Ranger run folder or its `outs/` | Standard 10x layout |
| `VisiumFiles` | A count matrix (`.h5`, or `matrix.mtx` + features + barcodes) and `tissue_positions(_list).csv[.gz]` | For GEO downloads. If one folder holds several samples, name the files in `counts_file` / `positions_file` |
| `Xenium` | Xenium Onboard Analysis output folder | Needs SpatialExperimentIO; control probes removed on import |
| `MERFISH` | Vizgen MERSCOPE output folder | Needs SpatialExperimentIO |
| `CosMx` | CosMx flat-file export folder | Needs SpatialExperimentIO |
| `.rds` | A saved SpatialExperiment | Needs `counts`, spatial coordinates, `sample_id` |

Sequencing-based (Visium) and imaging-based (Xenium, MERFISH, CosMx) data are not merged in one analysis, and different imaging platforms are not mixed.

## Sample sheet

One row per tissue section:

```
path,sample_id,condition,animal,platform
/data/study/sec01/outs,sec01,lesioned,mouse1,Visium
/data/study/sec02/outs,sec02,lesioned,mouse1,Visium
/data/study/sec03/outs,sec03,control,mouse2,Visium
```

Required: `path`, `sample_id`, `condition`. Optional: `platform` (default Visium), `counts_file` and `positions_file` (VisiumFiles), and any metadata such as animal, donor, sex, slide or batch, which is copied onto every spot or cell. Examples: `example_sample_sheet.csv`, `example_sample_sheet_geo.csv` and `GSE298680_sample_sheet.csv` (a worked GEO example; edit the path to your download folder).

## Workflow

1. **Import** several sections at once with a sample sheet. The replicates table shows whether DE is possible (`de_ready`).
2. **QC**: thresholds are set per platform on import. Imaging panels have no mitochondrial genes, so the mito filter is off for them.
3. **Normalise**: library size (default) or scran.
4. **Cluster** (optional): Joint mode gives domains shared across sections (Harmony-corrected PCs on highly variable genes plus within-section neighbourhood smoothing). Required for per-domain DE.
5. **Spatial statistics** (optional): Moran's I and Geary's C per section, one-sided, BH-corrected within section.
6. **Differential expression**: edgeR quasi-likelihood tests on pseudobulk counts per replicate (and per domain).
7. **Export**: "Save to folder" writes the processed object to disk. Load it later in `.rds` mode to continue. CSV downloads are available for spatial statistics and DE results.

## Differential expression design

- **Reference condition.** logFC is always *other condition vs reference*. The reference is detected automatically (control, vehicle, sham, WT and similar) and can be changed. Every result row records the comparison in a `contrast` column.
- **Replicates, not sections.** `sample_id` is a tissue section. If several sections come from the same animal or donor, choose that column as the replicate column: its sections are pooled into one pseudobulk replicate. Treating sections as independent replicates inflates significance. The app warns when a metadata column looks like it groups sections from the same animal or donor.
- **Covariates.** Technical factors present in both conditions, such as Visium slide or batch, can be added to the model (`~ covariates + condition`).
- **Checks.** A per-replicate expression plot shows whether an effect is consistent across replicates or driven by a few of them. A minimum number of spots or cells per replicate and domain prevents noisy tiny pseudobulks.

## Install and run

```r
install.packages(c("shiny", "ggplot2", "Matrix", "igraph", "harmony"))
BiocManager::install(c("SpatialExperiment", "scater", "scuttle", "scran", "edgeR", "spdep",
                       "BiocNeighbors", "DropletUtils", "SpatialExperimentIO"))
shiny::runApp(".")
```

Run the tests from the project root: `Rscript tests/smoke_test.R`

A Dockerfile is included for a containerised deployment.

## Repository layout

```
app.R                  Shiny application
R/01_import.R          Import, sample sheets, platform readers
R/02_qc.R              QC metrics and spot filtering
R/03_normalisation.R   Normalisation and batch correction
R/04_clustering.R      HVGs, PCA, per-section and joint spatial domains
R/05_spatial_stats.R   Moran's I and Geary's C
R/07_diff_analysis.R   Pseudobulk DE, replicates, covariates, plots
R/archive/             Earlier experimental modules (not loaded, not validated)
tests/smoke_test.R     Synthetic-data tests
```

## Citation

If you use spatialTXsuite, please cite it using the metadata in `CITATION.cff` (GitHub shows a "Cite this repository" button).

Mani M. *spatialTXsuite: an R Shiny application for replicate-aware spatial transcriptomics analysis.* Version 1.6.1, 2026. https://github.com/ManiMuralidharan/spatialTXsuite

Please also cite the methods the app relies on, in particular SpatialExperiment, edgeR, Harmony and spdep.

## License

MIT. See `LICENSE`.
