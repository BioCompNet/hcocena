## A named cutoff must reach the layer it names, and nothing else. Every case
## here is one the package used to get wrong or accept without comment.

fixture <- function() {
  path <- system.file("extdata", "hc_after_part1.rds", package = "hcocena")
  skip_if(!nzchar(path), "The after-part1 fixture is unavailable.")
  readRDS(path)
}
cutoffs <- function(hc) as.numeric(hc@config@layer$cutoff)
set_cut <- function(hc, v) suppressMessages(hc_set_cutoff(hc, cutoff_vector = v))

test_that("a layer can be addressed by the name the user gave it", {
  hc <- fixture()
  before <- cutoffs(hc)
  display <- as.character(hc@config@layer$layer_name)

  only_second <- set_cut(hc, stats::setNames(0.8, display[[2L]]))
  expect_equal(cutoffs(only_second), c(before[[1L]], 0.8))

  only_first <- set_cut(hc, stats::setNames(0.8, display[[1L]]))
  expect_equal(cutoffs(only_first), c(0.8, before[[2L]]))
})

test_that("names decide the assignment, not the order they are written in", {
  hc <- fixture()
  display <- as.character(hc@config@layer$layer_name)
  reversed <- set_cut(hc, stats::setNames(c(0.7, 0.9), c(display[[2L]], display[[1L]])))
  expect_equal(cutoffs(reversed), c(0.9, 0.7))
})

test_that("internal layer ids address the same layers", {
  hc <- fixture()
  ids <- as.character(hc@config@layer$layer_id)
  display <- as.character(hc@config@layer$layer_name)
  by_id <- set_cut(hc, stats::setNames(c(0.7, 0.9), c(ids[[2L]], ids[[1L]])))
  by_name <- set_cut(hc, stats::setNames(c(0.7, 0.9), c(display[[2L]], display[[1L]])))
  expect_equal(cutoffs(by_id), cutoffs(by_name))
})

test_that("a named value does not spill onto the layers it did not name", {
  hc <- fixture()
  before <- cutoffs(hc)
  skip_if(isTRUE(all.equal(before[[1L]], before[[2L]])),
          "needs two different starting cutoffs to be meaningful")
  one <- set_cut(hc, stats::setNames(0.8, as.character(hc@config@layer$layer_name)[[2L]]))
  expect_equal(cutoffs(one)[[1L]], before[[1L]])
})

test_that("an unknown layer name is refused", {
  hc <- fixture()
  expect_error(set_cut(hc, c(no_such_layer = 0.8)), "does not exist")
})

test_that("naming the same layer twice is refused", {
  hc <- fixture()
  nm <- as.character(hc@config@layer$layer_name)[[1L]]
  expect_error(set_cut(hc, stats::setNames(c(0.7, 0.8), c(nm, nm))), "more than once")
})

test_that("cutoffs outside the correlation range are refused", {
  hc <- fixture()
  expect_error(set_cut(hc, c(2.5, 2.5)), "\\[-1, 1\\]")
  expect_error(set_cut(hc, c(-5, 0.5)), "\\[-1, 1\\]")
})

test_that("a vector of the wrong length is refused", {
  hc <- fixture()
  expect_error(set_cut(hc, c(0.6, 0.7, 0.8)), "3 values for 2 layers")
})

test_that("a single unnamed value applies to every layer", {
  hc <- fixture()
  all_layers <- set_cut(hc, 0.75)
  expect_equal(cutoffs(all_layers), rep(0.75, nrow(hc@config@layer)))
})

test_that("NA means 'not set for this layer' rather than shifting the rest", {
  hc <- fixture()
  before <- cutoffs(hc)
  partial <- set_cut(hc, c(NA, 0.5))
  expect_equal(cutoffs(partial), c(before[[1L]], 0.5))
})
