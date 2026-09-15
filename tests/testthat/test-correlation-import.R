## Imported correlation and p-value matrices are written independently and may
## be sorted differently. Aligning them by gene name is the difference between
## correct p-values and p-values attached to the wrong gene pairs.

write_matrix <- function(m, path) {
  utils::write.table(m, path, quote = FALSE, row.names = FALSE,
                     col.names = TRUE, sep = "\t")
  path
}

toy_pair <- function(dir) {
  genes <- c("GA", "GB", "GC", "GD")
  set.seed(11)
  r <- matrix(stats::runif(16, -1, 1), 4, dimnames = list(genes, genes))
  r[lower.tri(r)] <- t(r)[lower.tri(r)]
  diag(r) <- 1
  p <- matrix(stats::runif(16, 0, 1), 4, dimnames = list(genes, genes))
  p[lower.tri(p)] <- t(p)[lower.tri(p)]
  diag(p) <- 0
  list(genes = genes, r = r, p = p,
       r_path = write_matrix(r, file.path(dir, "r.txt")),
       p_path = write_matrix(p, file.path(dir, "p.txt")))
}

test_that("a differently sorted p-value matrix is aligned, not misread", {
  dir <- withr::local_tempdir()
  tp <- toy_pair(dir)

  shuffled <- tp$p[rev(tp$genes), rev(tp$genes), drop = FALSE]
  shuffled_path <- write_matrix(shuffled, file.path(dir, "p_rev.txt"))

  straight <- hcocena:::.hc_read_correlation_import(c(tp$r_path, tp$p_path), layer = 1L)
  reversed <- hcocena:::.hc_read_correlation_import(c(tp$r_path, shuffled_path), layer = 1L)

  expect_equal(straight$r, reversed$r)
  expect_equal(straight$P, reversed$P)
  # and the p-values really do follow the correlation matrix's gene order
  expect_identical(rownames(reversed$P), rownames(reversed$r))
  expect_equal(reversed$P["GA", "GC"], tp$p["GA", "GC"])
})

test_that("matrices describing different genes are refused", {
  dir <- withr::local_tempdir()
  tp <- toy_pair(dir)
  wrong <- tp$p
  colnames(wrong)[[2]] <- "GX"
  rownames(wrong)[[2]] <- "GX"
  wrong_path <- write_matrix(wrong, file.path(dir, "p_wrong.txt"))

  expect_error(
    hcocena:::.hc_read_correlation_import(c(tp$r_path, wrong_path), layer = 1L),
    "different genes"
  )
})

test_that("a single path, a missing file and impossible values are refused", {
  dir <- withr::local_tempdir()
  tp <- toy_pair(dir)

  expect_error(hcocena:::.hc_read_correlation_import(tp$r_path, layer = 1L),
               "exactly two file paths")
  expect_error(
    hcocena:::.hc_read_correlation_import(c(tp$r_path, file.path(dir, "nope.txt")), layer = 1L),
    "file not found"
  )

  too_big <- tp$r
  too_big[1, 2] <- 4
  too_big[2, 1] <- 4
  big_path <- write_matrix(too_big, file.path(dir, "r_big.txt"))
  expect_error(
    hcocena:::.hc_read_correlation_import(c(big_path, tp$p_path), layer = 1L),
    "\\[-1, 1\\]"
  )
})

test_that("duplicate gene names are refused", {
  dir <- withr::local_tempdir()
  tp <- toy_pair(dir)
  dup <- tp$r
  colnames(dup)[[2]] <- "GA"
  dup_path <- write_matrix(dup, file.path(dir, "r_dup.txt"))
  expect_error(
    hcocena:::.hc_read_correlation_import(c(dup_path, tp$p_path), layer = 1L),
    "repeats gene name"
  )
})
