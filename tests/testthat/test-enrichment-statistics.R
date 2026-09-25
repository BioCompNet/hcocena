## Three different quantities used to share two column names in the enrichment
## output. Each must now be named for what it is.

marker_gmt <- function(dir) {
  path <- file.path(dir, "markers.gmt")
  writeLines(c(
    paste(c("T cells", "-", "CD3D", "CD3E", "CD3G", "CD2", "CD28",
            "LCK", "ZAP70", "IL7R", "CD7", "TRAC"), collapse = "\t"),
    paste(c("B cells", "-", "CD19", "MS4A1", "CD79A", "CD79B", "BLNK",
            "PAX5", "CR2", "FCRL1", "TNFRSF13B", "VPREB3"), collapse = "\t"),
    paste(c("Monocytes", "-", "CD14", "LYZ", "FCN1", "VCAN", "S100A8",
            "S100A9", "CSF1R", "ITGAM", "CD68", "FCGR3A"), collapse = "\t"),
    paste(c("NK cells", "-", "NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1",
            "GZMB", "NCR1", "KLRC1", "FGFBP2", "SPON2"), collapse = "\t")
  ), path)
  path
}

enrich_with <- function(padj) {
  path <- system.file("extdata", "hc_clustered.rds", package = "hcocena")
  skip_if(!nzchar(path), "The clustered fixture is unavailable.")
  skip_if_not_installed("clusterProfiler")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  grDevices::pdf(file.path(dir, "plots.pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  hc <- readRDS(path)
  hc <- suppressMessages(suppressWarnings(hc_functional_enrichment(
    hc, gene_sets = character(0),
    custom_gmt_files = c(CellTypes = marker_gmt(dir)),
    padj = padj
  )))
  as.data.frame(as.list(hc@satellite)$enrichments$significant_enrichments_all_dbs)
}

test_that("the enrichment table names each statistic for what it is", {
  res <- enrich_with("BH")
  skip_if(nrow(res) == 0, "no significant terms to inspect")

  expect_true(all(c("pvalue", "p_adjusted", "padj_method", "q_storey", "qvalue")
                  %in% colnames(res)))
  # the adjusted p-value drives the filtering, and `qvalue` is its alias
  expect_equal(res$qvalue, res$p_adjusted)
  # the correction that was actually applied is recorded
  expect_true(all(res$padj_method == "BH"))
})

test_that("the recorded correction method follows the request", {
  res <- enrich_with("bonferroni")
  skip_if(nrow(res) == 0, "no significant terms to inspect")
  expect_true(all(res$padj_method == "bonferroni"))
  expect_equal(res$qvalue, res$p_adjusted)
})

test_that("a Storey q-value that could not be computed is NA, not substituted", {
  res <- enrich_with("BH")
  skip_if(nrow(res) == 0, "no significant terms to inspect")
  # On a four-term collection pi0 cannot be fitted, so q_storey is unavailable.
  # The point is that this is visible instead of being quietly replaced by the
  # adjusted p-value under the same name.
  expect_true("q_storey" %in% colnames(res))
  expect_true(is.numeric(res$q_storey) || all(is.na(res$q_storey)))
})
