# hcocena

`hcocena` is a Bioconductor-oriented software package for horizontal
integration and downstream analysis of transcriptomics datasets. It provides an
S4 workflow centered on `HCoCenaExperiment`, with support for
`MultiAssayExperiment` / `SummarizedExperiment`, while keeping compatibility
with the historical `hcobject` workflow.

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
```

The container workspace is prepared at `/home/rstudio/hcocena` and includes
`reference_files/` so users can start without downloading those files
separately.

## Minimal workflow

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
```

See the package vignettes for a reproducible example based on the bundled toy
data.

## Citation

Marie Oestreich, Lisa Holsten, Shobhit Agrawal, Kilian Dahm, Philipp Koch,
Han Jin, Matthias Becker, Thomas Ulas (2022). "hCoCena: horizontal integration
and analysis of transcriptomics datasets." *Bioinformatics* 38(20):4727-4734.
doi:10.1093/bioinformatics/btac589
