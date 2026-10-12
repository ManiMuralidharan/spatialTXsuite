# spatialTXsuite news

## 1.6.1 (2026-10-12)
- Windows: launchers (`install_packages.bat`, `run_tests.bat`, `run_app.bat`) that locate R without PATH changes; Windows drive, backslash and UNC paths handled in sample sheets (trailing separators no longer break folder checks).
- Differential expression uses an explicit reference (baseline) condition, detected automatically from names such as control, vehicle, sham or WT. Previously the alphabetically first level was the baseline, which reversed logFC when a treatment name began with a digit (e.g. "6OHDA"). Results carry a `contrast` column.
- Export tab: "Save to folder" writes the processed object directly to disk (reliable for multi-GB objects).
- Per-replicate expression plot searches all genes.

## 1.6.0
- Technical covariates (e.g. Visium slide, batch) in the edgeR design.
- Per-replicate pseudobulk expression plot.

## 1.5.0
- Import of GEO-style Visium files (`.h5` or `matrix.mtx` plus `tissue_positions`), Xenium, MERFISH and CosMx; control probes removed for imaging panels.
- Replicate column for DE: sections from the same animal or donor are pooled; likely pseudoreplication is flagged.
- Moran's I computed in gene chunks for large imaging sections.

## 1.4.x
- Multi-sample import from a sample sheet.
- QC metrics, log-normalisation and highly variable gene selection computed directly (no deprecated scater/scuttle/scran calls); robust binned variance trend.
- Gene symbols in spatial statistics and DE results.

## 1.3.0
- Joint spatial domains shared across sections (Harmony plus neighbourhood smoothing), required for per-domain DE.
- Minimum spots per replicate in pseudobulk.

## 1.2.x
- Mitochondrial and ribosomal genes detected from gene symbols (real Visium data has Ensembl row names).
- Stale results cleared when QC or normalisation is re-run.
- HVG-based PCA; kNN neighbour search instead of dense distance matrices; one-sided Moran's I.
