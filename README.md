# hcocena

[![bioc-check](https://github.com/BioCompNet/hcocena/actions/workflows/bioc-check.yaml/badge.svg)](https://github.com/BioCompNet/hcocena/actions/workflows/bioc-check.yaml)

`hcocena` is an R package for horizontal integration and downstream analysis of
transcriptomics datasets. It combines a modern S4 workflow built around
`HCoCenaExperiment` with compatibility for the historical `hcobject` workflow,
so new analyses and older projects can live in the same package.

![hcocena overview](.github/assets/hcocena-overview.png)

The package supports both multi-layer integration, such as RNA-seq plus array
data, and single-layer analyses using the same API. The focus is a
module-centric workflow: from data import and correlation-based network
construction to clustering, heatmaps, functional enrichment, upstream
inference, cell-type annotation, longitudinal analysis, and optional
LLM-assisted module interpretation.

## What hcocena provides

- S4-first workflow with `HCoCenaExperiment`, `MultiAssayExperiment`, and
  `SummarizedExperiment`
- Backward-compatible support for the legacy `hcobject` workflow
- Correlation cutoff tuning and automatic cutoff selection helpers
- Clustering, integrated network construction, module splitting, and hCoCena
  heatmaps
- Functional enrichment across multiple databases with export helpers
- Upstream inference with DoRothEA and PROGENy via `decoupleR`
- Cell-type annotation helpers and reference-data preview utilities
- Longitudinal module and endotype analyses
- A Docker workflow with bundled `reference_files` for a ready-to-run setup

## Repository structure

- Package source is at the repository root and follows a Bioconductor-style
  layout
- Docker support lives in [`docker/`](docker), including bundled
  `reference_files`
- CI for package checks is defined in [`.github/workflows/bioc-check.yaml`](.github/workflows/bioc-check.yaml)

## Installation

After Bioconductor acceptance:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("hcocena")
```

For local development or pre-submission testing from a checkout:

```r
install.packages("remotes")
remotes::install_local(".", dependencies = TRUE, upgrade = "never")
```

## Docker

To build a ready-to-use RStudio image from this repository:

```bash
docker build -f docker/Dockerfile -t hcocena .
docker run --rm -p 8787:8787 -e PASSWORD=hcocena hcocena
```

The container prepares a workspace at `/home/rstudio/hcocena` and includes:

- the local `hcocena` installation
- bundled `reference_files/`
- preinstalled optional packages for common workflows, including
  `CALIBERrfimpute`, `RCy3`, `SpatialExperiment`, and `GSVA`
- empty `count_data`, `annotation_data`, and `output` directories

See [`docker/README.md`](docker/README.md) for the Docker-specific notes.

## Minimal S4 workflow

```r
library(hcocena)

hc <- hc_init()
hc <- hc_set_paths(
  hc,
  dir_count_data = "/path/to/counts/",
  dir_annotation = "/path/to/annotation/",
  dir_reference_files = "/path/to/reference/",
  dir_output = tempdir()
)
hc <- hc_define_layers(
  hc,
  data_sets = list(
    Layer1 = c("counts.tsv", "anno.tsv")
  )
)
hc <- hc_read_data(
  hc,
  gene_symbol_col = "SYMBOL",
  sample_col = "SampleID",
  count_has_rn = FALSE,
  anno_has_rn = FALSE
)
hc <- hc_run_expression_analysis_1(hc, export = FALSE)
hc <- hc_plot_cutoffs(hc, interactive = FALSE)
```

For a reproducible package-based example, install `hcocena` and run:

```r
browseVignettes("hcocena")
```

The package ships toy data and prepared example objects in `inst/extdata` to
support documentation, testing, and manual smoke tests.

## Legacy workflow

Legacy APIs are still available. For mixed projects or migration work, the
conversion helpers below bridge between the old and new object models:

- `as_hcocena()` converts legacy objects to the S4 workflow
- `as_hcobject()` converts S4 objects back to the legacy format

One practical note: when using longitudinal imputation with
`impute_method = "rfcont"`, attach `CALIBERrfimpute` in the session first:

```r
library(CALIBERrfimpute)
```

## Documentation and references

- Bioinformatics paper:
  https://doi.org/10.1093/bioinformatics/btac589
- STAR Protocol:
  https://star-protocols.cell.com/protocols/3341

## Citation

Marie Oestreich, Lisa Holsten, Shobhit Agrawal, Kilian Dahm, Philipp Koch,
Han Jin, Matthias Becker, Thomas Ulas (2022). "hCoCena: horizontal integration
and analysis of transcriptomics datasets." *Bioinformatics* 38(20):4727-4734.
doi:10.1093/bioinformatics/btac589
