## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----libraries----------------------------------------------------------------
library(hcocena)

## ----roundtrip-conversion-----------------------------------------------------
hc <- hc_init()
legacy <- as_hcobject(hc)
hc_roundtrip <- as_hcocena(legacy)

class(legacy)
sort(names(legacy))[1:8]
methods::is(hc_roundtrip, "HCoCenaExperiment")
slotNames(hc_roundtrip)

## ----workspace-example, eval=FALSE--------------------------------------------
# # load a workspace that contains a legacy `hcobject`
# load("legacy_workspace_with_hcobject.RData")
# 
# hc <- as_hcocena(hcobject)
# 
# # use object-based API
# hc <- hc_set_paths(
#   hc,
#   dir_count_data = FALSE,
#   dir_annotation = FALSE,
#   dir_reference_files = tempdir(),
#   dir_output = tempdir()
# )
# 
# # if needed: convert back for legacy custom scripts
# hcobject <- as_hcobject(hc)

## ----session-info-------------------------------------------------------------
sessionInfo()

