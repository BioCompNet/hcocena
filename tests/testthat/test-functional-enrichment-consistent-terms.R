test_that("default term selection carries union-top terms into other modules", {
  select_rows <- get(".hc_select_enrichment_rows", asNamespace("hcocena"))

  significant_rows <- list(
    blue = data.frame(
      Description = c("TermA", "TermB", "TermC"),
      qvalue = c(0.001, 0.010, 0.020),
      rank = 1:3,
      stringsAsFactors = FALSE
    ),
    red = data.frame(
      Description = c("TermD", "TermA", "TermE", "TermB"),
      qvalue = c(0.002, 0.003, 0.004, 0.020),
      rank = 1:4,
      stringsAsFactors = FALSE
    )
  )

  out <- select_rows(
    significant_rows = significant_rows,
    top = 2L,
    cluster_levels = c("blue", "red")
  )

  expect_identical(out$top_enr$blue, c("TermA", "TermB"))
  expect_identical(out$top_enr$red, c("TermD", "TermA", "TermB"))
  expect_identical(out$selected_rows$red$rank, c(1L, 2L, 4L))
})

test_that("legacy term selection keeps per-module top terms only", {
  select_rows <- get(".hc_select_enrichment_rows", asNamespace("hcocena"))

  significant_rows <- list(
    blue = data.frame(
      Description = c("TermA", "TermB", "TermC"),
      qvalue = c(0.001, 0.010, 0.020),
      rank = 1:3,
      stringsAsFactors = FALSE
    ),
    red = data.frame(
      Description = c("TermD", "TermA", "TermE", "TermB"),
      qvalue = c(0.002, 0.003, 0.004, 0.020),
      rank = 1:4,
      stringsAsFactors = FALSE
    )
  )

  out <- select_rows(
    significant_rows = significant_rows,
    top = 2L,
    cluster_levels = c("blue", "red"),
    consistent_terms = FALSE
  )

  expect_identical(out$top_enr$blue, c("TermA", "TermB"))
  expect_identical(out$top_enr$red, c("TermD", "TermA"))
  expect_identical(out$selected_rows$red$rank, c(1L, 2L))
})
