## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----install, eval = FALSE----------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install("hcocena")

## ----libraries----------------------------------------------------------------
library(hcocena)
library(MultiAssayExperiment)
library(SummarizedExperiment)

## ----toy-workflow-------------------------------------------------------------
extdir <- normalizePath(system.file("extdata", package = "hcocena"), winslash = "/")
extdir <- paste0(extdir, "/")
outdir <- file.path(tempdir(), "hcocena-vignette-output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

hc <- hc_init()

hc <- hc_set_paths(
  hc,
  dir_count_data = extdir,
  dir_annotation = extdir,
  dir_reference_files = extdir,
  dir_output = outdir
)

hc <- hc_define_layers(
  hc,
  data_sets = list(
    Layer1 = c("toy_layer1_counts.tsv", "toy_layer1_anno.tsv"),
    Layer2 = c("toy_layer2_counts.tsv", "toy_layer2_anno.tsv")
  )
)

hc <- hc_read_data(
  hc,
  sep_counts = "\t",
  sep_anno = "\t",
  gene_symbol_col = "SYMBOL",
  sample_col = "SampleID",
  count_has_rn = FALSE,
  anno_has_rn = FALSE
)

hc <- hc_set_global_settings(
  hc,
  organism = "human",
  control_keyword = "control",
  variable_of_interest = "group",
  data_in_log = TRUE
)

## ----container-inspection-----------------------------------------------------
names(experiments(hc_mae(hc)))
lapply(experiments(hc_mae(hc)), dim)
as.data.frame(SummarizedExperiment::colData(experiments(hc_mae(hc))[[1L]]))

## ----config-inspection--------------------------------------------------------
as.data.frame(methods::slot(hc_config(hc), "global"))
as.data.frame(methods::slot(hc_config(hc), "layer"))[, c("layer_id", "count_source", "annotation_source")]

## ----session-info-------------------------------------------------------------
sessionInfo()

