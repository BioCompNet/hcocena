test_that("rfcont imputation works with CALIBERrfimpute installed but not attached", {
  skip_if_not_installed("mice")
  skip_if_not_installed("CALIBERrfimpute")

  rf_pkg_search <- "package:CALIBERrfimpute"
  if (rf_pkg_search %in% search()) {
    detach(rf_pkg_search, unload = FALSE, character.only = TRUE)
  }

  expect_false(rf_pkg_search %in% search())
  set.seed(42)
  out <- suppressWarnings(
    hcocena:::.hc_longitudinal_impute_time_data(
      time_data = data.frame(
        donor = c("d1", "d2", "d3", "d4"),
        `1` = c(1, 2, 3, 4),
        `2` = c(2, 3, NA, 5),
        `3` = c(3, 4, 5, 6),
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      donor_col = "donor",
      method = "rfcont",
      ntree = 10,
      m = 2,
      maxit = 1,
      seed = 42
    )
  )
  expect_s3_class(out, "data.frame")
  expect_equal(out$donor, c("d1", "d2", "d3", "d4"))
  expect_false(rf_pkg_search %in% search())
})
