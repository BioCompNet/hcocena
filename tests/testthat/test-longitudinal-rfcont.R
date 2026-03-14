test_that("rfcont imputation requires CALIBERrfimpute to be attached", {
  skip_if_not_installed("mice")
  skip_if_not_installed("CALIBERrfimpute")

  rf_pkg_search <- "package:CALIBERrfimpute"
  if (rf_pkg_search %in% search()) {
    detach(rf_pkg_search, unload = FALSE, character.only = TRUE)
  }

  expect_false(rf_pkg_search %in% search())
  expect_error(
    hcocena:::.hc_legacy_impute_time_data(
      time_data = data.frame(
        donor = c("d1", "d2"),
        `1` = c(1, NA),
        `2` = c(2, 3),
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      donor_col = "donor",
      method = "rfcont",
      ntree = 10,
      m = 2,
      maxit = 1,
      seed = 42
    ),
    "library\\(CALIBERrfimpute\\)"
  )
})
