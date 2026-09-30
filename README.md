# hcocena

[![R-CMD-check](https://github.com/BioCompNet/hcocena/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BioCompNet/hcocena/actions/workflows/R-CMD-check.yaml)

`hcocena` is an R package for horizontal integration and downstream analysis of
transcriptomics datasets. It provides a modern S4 workflow built around
`HCoCenaExperiment` for reproducible network-centric transcriptomics analyses.

![hcocena overview](.github/assets/hcocena-overview.jpg)

The package supports both multi-layer integration, such as RNA-seq plus array
data, and single-layer analyses using the same API. The focus is a
module-centric workflow: from data import and correlation-based network
construction to clustering, heatmaps, functional enrichment, upstream
inference, cell-type annotation, longitudinal analysis, and optional
LLM-assisted module interpretation.

## What hcocena provides

- S4-first workflow with `HCoCenaExperiment`, `MultiAssayExperiment`, and
  `SummarizedExperiment`
- Correlation cutoff tuning and automatic cutoff selection helpers
- Clustering, integrated network construction, module splitting, and hCoCena
  heatmaps
- Functional enrichment across multiple databases with export helpers
- Upstream inference with DoRothEA and PROGENy via `decoupleR`
- Cell-type annotation helpers and reference-data preview utilities
- Longitudinal module and endotype analyses

## Repository structure

- Package source is at the repository root and follows the Bioconductor
  layout
- Analysis templates are kept in [`inst/scripts/workflows/`](inst/scripts/workflows/)
- CI for R CMD check and BiocCheck (release and Bioconductor devel) is defined
  in [`.github/workflows/R-CMD-check.yaml`](.github/workflows/R-CMD-check.yaml)
- The Docker image (RStudio with hcocena and bundled reference files) and the
  development history before the Bioconductor submission are maintained in
  [BioCompNet/hcocena-dev](https://github.com/BioCompNet/hcocena-dev)

## Contributors

- [Waqar Hanif](https://github.com/waqarhanif-biocode)

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

## Analysis templates

The longer walkthroughs the package is normally driven from are installed with
it, so they are available from an installed copy and not only from a clone:

- `inst/scripts/workflows/hcocena_main.Rmd` -- full analysis, import to enrichment
- `inst/scripts/workflows/hcocena_satellite.Rmd` -- optional downstream analyses

```r
dir(system.file("scripts", "workflows", package = "hcocena"))

file.copy(
  system.file("scripts", "workflows", "hcocena_main.Rmd", package = "hcocena"),
  "hcocena_main.Rmd"
)
```

They are templates, not reproducible documents: they point at your own count
and annotation files, and several steps contact external services (ChEA3,
Enrichr, DoRothEA, Cytoscape, LLM providers). The built vignette
(`vignette("hcocena-s4-workflow")`) is the runnable short version.

One practical note: when using longitudinal imputation with
`impute_method = "rfcont"`, attach `CALIBERrfimpute` in the session first:

```r
library(CALIBERrfimpute)
```

## Documentation and references

- Method paper: Oestreich et al. (2022), *Bioinformatics*,
  https://doi.org/10.1093/bioinformatics/btac589
- Step-by-step protocol: Holsten et al. (2024), *STAR Protocols*,
  https://doi.org/10.1016/j.xpro.2024.102922

### Coming from the papers

Both papers use the earlier hCoCena interface, in which the functions worked
on a global `hcobject` and had no `hc_` prefix (repository
`MarieOestreich/hCoCena`, Docker image `mo126/hcocena`). In this package every
step takes and returns an `HCoCenaExperiment` (`hc <- hc_step(hc, ...)`), and
the current Docker image is `therealtomek/hcocena`. The functions map as
follows:

| In the papers | In this package |
|---|---|
| `init_wd()`, `check_dirs()`, `init_save_folder()` | `hc_init()`, `hc_set_paths()`, `hc_check_dirs()`, `hc_init_save_folder()` |
| `define_layers()`, `read_data()` | `hc_define_layers()`, `hc_read_data()` |
| `set_supp_files()`, `read_supplementary()` | `hc_set_supp_files()`, `hc_read_supplementary()` |
| `set_global_settings()`, `set_layer_settings()` | `hc_set_global_settings()`, `hc_set_layer_settings()` |
| `suggest_topvar()` | `hc_suggest_topvar()` |
| `run_expression_analysis_1()`, `plot_cutoffs()`, `set_cutoff()`, `plot_deg_dist()` | `hc_run_expression_analysis_1()`, `hc_plot_cutoffs()`, `hc_set_cutoff()`, `hc_plot_deg_dist()` |
| `run_expression_analysis_2()` | `hc_run_expression_analysis_2()` |
| `build_integrated_network()`, `plot_integrated_network()` | `hc_build_integrated_network()`, `hc_plot_integrated_network()` |
| `cluster_calculation()`, `plot_cluster_heatmap()` | `hc_cluster_calculation()`, `hc_plot_cluster_heatmap()` |
| `functional_enrichment()` | `hc_functional_enrichment()` |
| `TF_overrep_module()`, `TF_overrep_network()`, `check_tf()` | `hc_tf_overrep_module()`, `hc_tf_overrep_network()`, `hc_check_tf()` |
| `find_hubs()` | `hc_find_hubs()` |
| `change_grouping_parameter()`, `cut_hclust()` | `hc_change_grouping_parameter()`, `hc_cut_hclust()` |
| `write_session_info()` | `hc_write_session_info()` |

Results that the papers read from `hcobject` are available through
accessors such as `hc_clusters()`, `hc_gene_to_cluster()`, `hc_graph()` and
`hc_satellite()`.

## Citation

Please cite the method paper, and the protocol where you followed it:

Oestreich M, Holsten L, Agrawal S, Dahm K, Koch P, Jin H, Becker M, Ulas T
(2022). "hCoCena: horizontal integration and analysis of transcriptomics
datasets." *Bioinformatics* 38(20):4727-4734.
doi:10.1093/bioinformatics/btac589

Holsten L, Dahm K, Oestreich M, Becker M, Ulas T (2024). "hCoCena: A toolbox
for network-based co-expression analysis and horizontal integration of
transcriptomic datasets." *STAR Protocols* 5(1):102922.
doi:10.1016/j.xpro.2024.102922

In R: `citation("hcocena")`.
