test_that("regression: utils no longer references working_director", {
  src <- paste(deparse(get("run_expression_analysis_1_body", asNamespace("hcocena"))), collapse = "\n")
  expect_false(grepl("working_director", src, fixed = TRUE))
})


test_that("regression: rho path now guards optional propr dependency", {
  src <- paste(deparse(get("run_expression_analysis_1_body", asNamespace("hcocena"))), collapse = "\n")
  expect_true(grepl("requireNamespace\\(\"propr\"", src))
})
