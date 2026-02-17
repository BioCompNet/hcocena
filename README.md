# hCoCena
Horizontal integration and analysis of transcriptomics datasets  
Paper: https://doi.org/10.1093/bioinformatics/btac589

hCoCena is an R package for network-based transcriptomics analysis.  
It supports both:
- multi-layer integration (e.g. RNA-seq + array), and
- single-layer analyses (seq-only workflows).

The repository is maintained at: https://github.com/BioCompNet/hcocena

## Current Version
This repository currently tracks **hcocena v1.28**.

![hCoCena overview v1.28](assets/hCocena_updated_1.27.png)

## What hCoCena Can Do (v1.28)
- S4-first workflow (`HCoCenaExperiment`) with `hc_*` API wrappers.
- Backward compatibility with legacy workflows.
- Correlation cutoff tuning with tiered selection (`hc_tune_cutoff` / `hc_auto_tune_cutoff_tiered`).
- Automatic cutoff application via `hc_set_cutoff(auto = TRUE, fallback_cutoff = ...)`.
- Resolution tuning for clustering (`hc_tune_resolution` / `hc_auto_tune`).
- Leiden default partition set to `RBConfigurationVertexPartition`.
- Module splitting and undo (`hc_split_modules`, `hc_unsplit_modules`) including resolution testing.
- Advanced functional enrichment:
  - per DB plots,
  - mixed/combined DB plots,
  - selected terms and all-significant exports,
  - merged Excel outputs across DBs.
- Upstream inference (DoRothEA TF + PROGENy pathway via decoupleR):
  - `activity_input = "gfc"` or `"fc"` (with user-defined comparisons),
  - per-comparison plotting,
  - consistent-term mode for comparability.
- Module knowledge network plot combining enrichment + upstream results.
- Automatic module cell-type annotation with Enrichr DBs (`hc_celltype_annotation`), plus DB discovery/preview helpers:
  - `hc_list_celltype_databases`
  - `hc_preview_celltype_database`

## Repository Layout
- `hCoCena-r-package/`: R package source.
- `hcocena_main.Rmd`: main multi-layer workflow.
- `hcocena_main_seq_only.Rmd`: main single-layer workflow.
- `hcocena_satellite.Rmd`: optional analyses and extended utilities.
- `legacy_rmd/`: archived legacy script variants.
- `STAR_protocol/`: STAR protocol-oriented materials.
- `reference_files/`: reference GMT/TF files for enrichment workflows.
- `scripts/`: maintenance/helpers (e.g. reinstall helpers).

## Installation

### Option A: Reproducible install using provided scripts
Run, in order:
1. `install_versioned_dependencies.R`
2. `install_hcocena.R`

### Option B: Direct install from GitHub
```r
install.packages("remotes")
remotes::install_github("BioCompNet/hcocena", subdir = "hCoCena-r-package", dependencies = TRUE)
```

### Option C: Install from local checkout
```r
install.packages("remotes")
remotes::install_local("hCoCena-r-package", dependencies = TRUE, upgrade = "never")
```

## Minimal S4 Workflow
```r
library(hcocena)

hc <- hc_init()
# ... set paths/layers/settings/import ...
hc <- hc_run_expression_analysis_1(hc, corr_method = "pearson")
hc <- hc_set_cutoff(hc, auto = TRUE, fallback_cutoff = 0.982)
hc <- hc_run_expression_analysis_2(hc)
hc <- hc_build_integrated_network(hc, mode = "u", multi_edges = "min")
hc <- hc_cluster_calculation(hc, cluster_algo = "cluster_leiden")
hc <- hc_plot_cluster_heatmap(hc)
hc <- hc_functional_enrichment(hc, gene_sets = c("Hallmark", "Kegg", "Go"))
```

For end-to-end examples, use:
- `hcocena_main.Rmd`
- `hcocena_main_seq_only.Rmd`
- `hcocena_satellite.Rmd`

## Legacy Workflow
Legacy behavior is still available.  
If needed, use conversion helpers:
- `as_hcocena()` (legacy -> S4)
- `as_hcobject()` (S4 -> legacy)

## Citation
Marie Oestreich, Lisa Holsten, Shobhit Agrawal, Kilian Dahm, Philipp Koch, Han Jin, Matthias Becker, Thomas Ulas,  
**hCoCena: horizontal integration and analysis of transcriptomics datasets**, Bioinformatics, 38(20), 2022, 4727-4734.  
https://doi.org/10.1093/bioinformatics/btac589

STAR Protocol reference:  
https://star-protocols.cell.com/protocols/3341
