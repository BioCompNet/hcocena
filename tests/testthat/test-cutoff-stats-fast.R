# The fast cutoff-statistics loop (.hc_cutoff_stats_fast) must reproduce the
# original per-cutoff cutoff_prep() loop exactly. Component statistics are
# order- and vertex-label-invariant, which is what makes the integer/prefix
# optimisation valid, so this test pins that equivalence.

fast_loop <- get(".hc_cutoff_stats_fast", asNamespace("hcocena"))
cutoff_prep <- get("cutoff_prep", asNamespace("hcocena"))

old_loop <- function(df, range_cutoff, min_nodes) {
  do.call(rbind, lapply(
    range_cutoff,
    cutoff_prep,
    corrdf_r = df,
    print.all.plots = FALSE,
    x = 1L,
    min_nodes = min_nodes
  ))
}

make_corr_df <- function(n_genes, n_edges, seed = 1L) {
  set.seed(seed)
  genes <- sprintf("G%04d", seq_len(n_genes))
  i <- sample.int(n_genes, n_edges, replace = TRUE)
  j <- sample.int(n_genes, n_edges, replace = TRUE)
  keep <- i != j
  i <- i[keep]; j <- j[keep]
  df <- data.frame(
    V1 = genes[i],
    V2 = genes[j],
    rval = round(runif(length(i), 0.30, 0.99), 4),
    pval = runif(length(i), 0, 0.049),
    stringsAsFactors = FALSE
  )
  # drop duplicate undirected edges the way an upper-triangle correlation df has none
  df[!duplicated(paste(pmin(i, j), pmax(i, j))), ]
}

test_that("fast cutoff loop matches cutoff_prep loop (small graph)", {
  df <- make_corr_df(n_genes = 150, n_edges = 1200, seed = 42)
  range_cutoff <- round(seq(0.30, 0.95, length.out = 10), 3)
  for (mn in c(1L, 3L, 5L)) {
    expect_equal(
      fast_loop(df, range_cutoff, min_nodes = mn),
      old_loop(df, range_cutoff, min_nodes = mn),
      info = paste("min_nodes =", mn)
    )
  }
})

test_that("fast cutoff loop matches for a denser graph and wider cutoff range", {
  df <- make_corr_df(n_genes = 400, n_edges = 8000, seed = 7)
  range_cutoff <- round(seq(0.30, 0.98, length.out = 20), 3)
  expect_equal(
    fast_loop(df, range_cutoff, min_nodes = 4L),
    old_loop(df, range_cutoff, min_nodes = 4L)
  )
})

test_that("fast cutoff loop handles cutoffs above the maximum correlation (zero-edge rows)", {
  df <- make_corr_df(n_genes = 80, n_edges = 300, seed = 3)
  # include a cutoff higher than every rval -> zero passing edges
  range_cutoff <- c(0.5, 0.8, 0.999)
  fast <- fast_loop(df, range_cutoff, min_nodes = 3L)
  old <- old_loop(df, range_cutoff, min_nodes = 3L)
  expect_equal(fast, old)
  # the 0.999 cutoff row is the all-zero placeholder
  zrow <- fast[fast$cutoff == 0.999, ]
  expect_equal(zrow$no_edges, 0)
  expect_equal(zrow$no_nodes, 0)
})

union_find <- get(".hc_union_find_components", asNamespace("hcocena"))

test_that("union-find components match igraph::components", {
  skip_if_not_installed("igraph")
  canon <- function(membership) ave(seq_along(membership), membership, FUN = min)

  for (params in list(
    c(n = 200L, e = 500L, s = 1L),   # fragmented
    c(n = 200L, e = 4000L, s = 2L),  # near-fully-connected
    c(n = 500L, e = 300L, s = 3L)    # very sparse, many singletons
  )) {
    set.seed(params[["s"]])
    from <- sample.int(params[["n"]], params[["e"]], replace = TRUE)
    to <- sample.int(params[["n"]], params[["e"]], replace = TRUE)
    keep <- from != to
    from <- from[keep]; to <- to[keep]

    # Compact vertices to those present in edges (exactly what the cutoff loop
    # does via `unique(c(from_k, to_k))`), so there are no isolated vertices.
    # igraph::graph_from_edgelist otherwise drops trailing isolated vertices,
    # which would make the counts diverge on inputs the real code never produces.
    present <- sort(unique(c(from, to)))
    from <- match(from, present)
    to <- match(to, present)
    nv <- length(present)

    uf <- union_find(from, to, nv)
    g <- igraph::graph_from_edgelist(matrix(c(from, to), ncol = 2), directed = FALSE)
    ig <- igraph::components(g)

    # same partition (labels may differ) and same multiset of component sizes
    expect_identical(canon(uf$membership), canon(ig$membership))
    expect_identical(sort(uf$csize), sort(as.integer(ig$csize)))
  }
})

test_that("union-find handles isolated vertices (no edges)", {
  uf <- union_find(integer(0), integer(0), 5L)
  expect_equal(sort(uf$csize), rep(1L, 5))       # five singletons
  expect_equal(length(unique(uf$membership)), 5)
})

test_that("fast cutoff loop returns zero rows for an empty correlation table", {
  df <- data.frame(V1 = character(0), V2 = character(0), rval = numeric(0), pval = numeric(0))
  range_cutoff <- c(0.5, 0.7)
  fast <- fast_loop(df, range_cutoff, min_nodes = 3L)
  expect_equal(nrow(fast), 2)
  expect_true(all(fast$no_edges == 0))
  expect_equal(fast$cutoff, range_cutoff)
})
