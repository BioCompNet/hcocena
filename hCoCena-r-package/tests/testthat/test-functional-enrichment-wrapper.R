test_that("functional-enrichment panel keys keep single DBs before combined views", {
  enrich <- list(
    top_all_dbs_mixed = list(),
    top_Go = list(),
    panel_order = c("top_Hallmark", "top_Go", "top_all_dbs", "top_all_dbs_mixed"),
    top_all_dbs = list(),
    top_Hallmark = list()
  )

  keys <- hcocena:::.hc_functional_enrichment_panel_keys(enrich)
  expect_identical(
    keys,
    c("top_Hallmark", "top_Go", "top_all_dbs", "top_all_dbs_mixed")
  )
})

test_that("functional-enrichment panel keys are deduplicated", {
  enrich <- list(
    top_Hallmark = list(),
    top_all_dbs = list(),
    top_all_dbs_mixed = list(),
    panel_order = c("top_Hallmark", "top_Hallmark", "top_all_dbs", "top_all_dbs_mixed")
  )

  keys <- hcocena:::.hc_functional_enrichment_panel_keys(enrich)
  expect_identical(
    keys,
    c("top_Hallmark", "top_all_dbs", "top_all_dbs_mixed")
  )
})
