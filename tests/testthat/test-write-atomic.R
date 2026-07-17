# Robust export helpers: atomic write + post-write verification. These guard the
# module tables and heatmap PDF/PNG exports against truncated files on synced
# output folders (Sciebo / OneDrive).

write_atomic <- get(".hc_write_atomic", asNamespace("hcocena"))
verify_output <- get(".hc_verify_output_file", asNamespace("hcocena"))

test_that("atomic write produces the final file and leaves no temp behind", {
  dir <- withr::local_tempdir()
  final <- file.path(dir, "sub", "out.txt")
  res <- write_atomic(final, function(tmp) writeLines("hello", tmp))
  expect_identical(res, final)
  expect_true(file.exists(final))
  expect_identical(readLines(final), "hello")
  # no leftover *.part-* temp files
  leftovers <- list.files(dirname(final), pattern = "\\.part-", full.names = TRUE)
  expect_length(leftovers, 0)
})

test_that("atomic write overwrites an existing file", {
  dir <- withr::local_tempdir()
  final <- file.path(dir, "out.txt")
  writeLines("old", final)
  write_atomic(final, function(tmp) writeLines("new", tmp))
  expect_identical(readLines(final), "new")
})

test_that("atomic write errors and keeps the old file when the producer writes nothing", {
  dir <- withr::local_tempdir()
  final <- file.path(dir, "out.txt")
  writeLines("keep-me", final)
  # producer does not create tmp -> should error, not silently succeed
  expect_error(write_atomic(final, function(tmp) invisible(NULL)))
  # the pre-existing file is only removed once a valid temp exists, so it survives
  expect_true(file.exists(final))
  expect_identical(readLines(final), "keep-me")
  # no leftover temp files
  expect_length(list.files(dir, pattern = "\\.part-"), 0)
})

test_that("atomic write works for a real .xlsx payload", {
  skip_if_not_installed("openxlsx")
  dir <- withr::local_tempdir()
  final <- file.path(dir, "tbl.xlsx")
  tbl <- data.frame(genes = c("A", "B"), module = c("M1", "M1"), stringsAsFactors = FALSE)
  write_atomic(final, function(tmp) {
    openxlsx::write.xlsx(list(module_gene_list = tbl), file = tmp, overwrite = TRUE)
  })
  back <- openxlsx::read.xlsx(final)
  expect_equal(nrow(back), 2)
  expect_identical(colnames(back), c("genes", "module"))
})

test_that("verify_output_file warns on a missing or empty file", {
  dir <- withr::local_tempdir()
  missing <- file.path(dir, "nope.pdf")
  expect_warning(res <- verify_output(missing, label = "PDF"), "missing or empty")
  expect_false(res)

  empty <- file.path(dir, "empty.pdf")
  file.create(empty)
  expect_warning(verify_output(empty), "missing or empty")

  ok <- file.path(dir, "ok.txt")
  writeLines("content", ok)
  expect_silent(res_ok <- verify_output(ok))
  expect_true(res_ok)
})
