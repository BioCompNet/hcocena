# Regression: the fast correlation backend (.hc_fast_rcorr) must reproduce
# Hmisc::rcorr exactly for the values the pipeline consumes, i.e. the
# off-diagonal (upper.tri) correlation and p-value entries.

fast_rcorr <- get(".hc_fast_rcorr", asNamespace("hcocena"))

# Compare only what run_expression_analysis_1_body() actually reads back:
# the strict upper triangle of both r and P.
expect_matches_rcorr <- function(x, type) {
  ref <- Hmisc::rcorr(as.matrix(x), type = type)
  new <- fast_rcorr(x, type = type, backend = "auto")

  ut <- upper.tri(ref$r, diag = FALSE)
  # values agree to floating-point precision
  expect_equal(new$r[ut], ref$r[ut], tolerance = 1e-10)
  expect_equal(new$P[ut], ref$P[ut], tolerance = 1e-10)
  # missing-value patterns agree
  expect_identical(is.na(new$r[ut]), is.na(ref$r[ut]))
  expect_identical(is.na(new$P[ut]), is.na(ref$P[ut]))
  # gene names carried through
  expect_identical(colnames(new$r), colnames(ref$r))
  # pairwise observation counts retain the rcorr matrix contract
  expect_identical(new$n, ref$n)
}

set.seed(42)
n <- 60L
p <- 120L
x <- matrix(stats::rnorm(n * p), n, p)
colnames(x) <- paste0("g", seq_len(p))

test_that("fast backend matches Hmisc::rcorr for Pearson", {
  expect_matches_rcorr(x, "pearson")
})

test_that("fast backend matches Hmisc::rcorr for Spearman", {
  expect_matches_rcorr(x, "spearman")
})

test_that("fast backend matches rcorr with a constant (zero-variance) column", {
  xz <- x
  xz[, 5] <- 3.0
  expect_matches_rcorr(xz, "pearson")
  expect_matches_rcorr(xz, "spearman")
})

test_that("fast backend falls back to Hmisc::rcorr when NAs are present", {
  xna <- x
  xna[3, 7] <- NA
  xna[10, 20] <- NA
  ref <- Hmisc::rcorr(as.matrix(xna), type = "pearson")
  new <- suppressMessages(fast_rcorr(xna, type = "pearson", backend = "auto"))
  # exact fallback: identical to rcorr, including pairwise-complete p-values
  expect_equal(new$r, ref$r)
  expect_equal(new$P, ref$P)
})

test_that("backend = 'rcorr' forces the original path", {
  ref <- Hmisc::rcorr(as.matrix(x), type = "pearson")
  new <- fast_rcorr(x, type = "pearson", backend = "rcorr")
  expect_equal(new$r, ref$r)
  expect_equal(new$P, ref$P)
})

test_that("fast backend rejects invalid arguments", {
  expect_error(fast_rcorr(x, type = "kendall"), "pearson")
  expect_error(fast_rcorr(x, backend = "nope"))
})

test_that("pwcorr tolerates constant genes when finite pairs remain", {
  pwcorr_fun <- get("pwcorr", asNamespace("hcocena"))
  x_constant <- cbind(
    g1 = seq_len(12),
    g2 = seq_len(12) + c(rep(0, 11), 0.1),
    constant = rep(3, 12)
  )

  out <- pwcorr_fun(
    dd2 = x_constant,
    layer_set = list(min_corr = 0.2, range_cutoff_length = 5L),
    bayes = FALSE,
    prior = 2,
    alpha = 0.5,
    padj = "none",
    export = FALSE,
    layer = 1L,
    import = NULL,
    corr_method = "pearson"
  )

  expect_true(all(is.finite(out$range_cutoff)))
  expect_false(anyNA(out$correlation_df_filt))
  expect_true(all(out$correlation_df_filt$V1 != "constant"))
  expect_true(all(out$correlation_df_filt$V2 != "constant"))
})

test_that("pwcorr reports an informative error if every pair is undefined", {
  pwcorr_fun <- get("pwcorr", asNamespace("hcocena"))
  x_constant <- cbind(a = rep(1, 8), b = rep(2, 8), c = rep(3, 8))

  expect_error(
    pwcorr_fun(
      dd2 = x_constant,
      layer_set = list(min_corr = 0.2, range_cutoff_length = 5L),
      bayes = FALSE,
      prior = 2,
      alpha = 0.5,
      padj = "none",
      export = FALSE,
      layer = 1L,
      import = NULL,
      corr_method = "pearson"
    ),
    "No finite pairwise correlations"
  )
})
