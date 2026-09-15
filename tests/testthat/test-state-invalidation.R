## Results computed from inputs or parameters that have since changed must not
## stay in the object. Each of these reproduces a state the package used to
## accept silently.

fixture <- function() {
  path <- system.file("extdata", "hc_clustered.rds", package = "hcocena")
  skip_if(!nzchar(path), "The clustered fixture is unavailable.")
  readRDS(path)
}

has_modules <- function(hc) {
  cl <- tryCatch(as.list(hc@integration@cluster)[["cluster_information"]],
                 error = function(e) NULL)
  !is.null(cl) && nrow(as.data.frame(cl)) > 0
}
has_network <- function(hc) !is.null(hc@integration@graph)
n_part <- function(hc, part) {
  length(as.list(methods::slot(hc@layer_results[[1L]], part)))
}

publish_layers <- function(hc, env = parent.frame(), mutate = identity) {
  legacy <- hcocena:::as_hcobject(hc)
  for (nm in c("set1_counts", "set1_anno", "set2_counts", "set2_anno")) {
    value <- legacy$data[[nm]]
    if (grepl("_counts$", nm)) value <- mutate(value)
    assign(nm, value, envir = env)
  }
  invisible(NULL)
}

test_that("re-reading changed counts discards network and modules", {
  hc <- fixture()
  env <- new.env(parent = globalenv())
  set.seed(4)
  publish_layers(hc, env = env, mutate = function(m) {
    m[] <- matrix(stats::runif(length(m), 5, 50), nrow = nrow(m))
    m
  })

  changed <- suppressMessages(
    hcocena:::.hc_read_data_impl(hc, count_has_rn = TRUE, anno_has_rn = TRUE,
                                 auto_setup_output = FALSE, source_env = env)
  )
  expect_false(has_network(changed))
  expect_false(has_modules(changed))
  expect_identical(n_part(changed, "part1"), 0L)
})

test_that("re-reading identical counts keeps the analysis", {
  hc <- fixture()
  env <- new.env(parent = globalenv())
  publish_layers(hc, env = env)

  same <- suppressMessages(
    hcocena:::.hc_read_data_impl(hc, count_has_rn = TRUE, anno_has_rn = TRUE,
                                 auto_setup_output = FALSE, source_env = env)
  )
  expect_true(has_network(same))
  expect_true(has_modules(same))
})

test_that("changing layer settings discards everything computed from them", {
  hc <- fixture()
  changed <- suppressMessages(hc_set_layer_settings(
    hc, top_var = c(10, 10), min_corr = c(0.9, 0.9),
    range_cutoff_length = c(20, 20), print_distribution_plots = c(FALSE, FALSE)
  ))
  expect_identical(n_part(changed, "part1"), 0L)
  expect_false(has_network(changed))
  expect_false(has_modules(changed))
})

test_that("re-applying the same layer settings keeps the analysis", {
  hc <- fixture()
  cfg <- hc@config@layer
  same <- suppressMessages(hc_set_layer_settings(
    hc,
    top_var = as.numeric(cfg$top_var),
    min_corr = as.numeric(cfg$min_corr),
    range_cutoff_length = as.numeric(cfg$range_cutoff_length),
    print_distribution_plots = as.logical(cfg$print_distribution_plots)
  ))
  expect_true(has_network(same))
  expect_true(has_modules(same))
})

test_that("a new cutoff discards the network but keeps the correlations", {
  hc <- fixture()
  before_part1 <- n_part(hc, "part1")
  changed <- suppressMessages(hc_set_cutoff(hc, cutoff_vector = c(0.95, 0.95)))

  # correlations do not depend on the cutoff and are expensive, so they stay
  expect_identical(n_part(changed, "part1"), before_part1)
  expect_identical(n_part(changed, "part2"), 0L)
  expect_false(has_network(changed))
  expect_false(has_modules(changed))
  # and the object must not report a cutoff its network was not built with
  expect_equal(as.numeric(changed@config@layer$cutoff), c(0.95, 0.95))
})

test_that("re-applying the same cutoff keeps the analysis", {
  hc <- fixture()
  same <- suppressMessages(
    hc_set_cutoff(hc, cutoff_vector = as.numeric(hc@config@layer$cutoff))
  )
  expect_true(has_network(same))
  expect_true(has_modules(same))
})

test_that("downstream analysis on a discarded state fails with an explanation", {
  hc <- fixture()
  changed <- suppressMessages(hc_set_cutoff(hc, cutoff_vector = c(0.95, 0.95)))
  expect_error(hc_get_module_scores(changed), "No module assignment available")
})
