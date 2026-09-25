test_that("hc_auto_tune applies layer-wise optimal cutoffs", {
  hc <- hc_init()
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = c("set1", "set2"),
    layer_name = c("A", "B"),
    count_source = c("set1_counts", "set2_counts"),
    annotation_source = c("set1_anno", "set2_anno")
  )

  hc@layer_results <- S4Vectors::SimpleList(
    set1 = new(
      "HCoCenaLayerResult",
      part1 = S4Vectors::SimpleList(cutoff_calc_out = list(optimal_cutoff = 0.81)),
      part2 = S4Vectors::SimpleList()
    ),
    set2 = new(
      "HCoCenaLayerResult",
      part1 = S4Vectors::SimpleList(cutoff_calc_out = list(optimal_cutoff = 0.93)),
      part2 = S4Vectors::SimpleList()
    )
  )

  tuned <- hc_auto_tune(
    hc,
    tune_cutoff = TRUE,
    tune_clustering = FALSE,
    apply = TRUE,
    verbose = FALSE
  )

  expect_true("cutoff" %in% colnames(tuned@config@layer))
  expect_equal(round(as.numeric(tuned@config@layer$cutoff), 2), c(0.81, 0.93))
  expect_true("auto_tune" %in% names(tuned@satellite))
})

test_that("hc_auto_tune evaluates and applies clustering recommendation", {
  edge_df <- data.frame(
    from = c("g1", "g1", "g2", "g4", "g4", "g5", "g3"),
    to = c("g2", "g3", "g3", "g5", "g6", "g6", "g4"),
    weight = c(0.91, 0.87, 0.88, 0.93, 0.95, 0.92, 0.25),
    stringsAsFactors = FALSE
  )
  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
  igraph::E(g)$weight <- edge_df$weight

  hc <- hc_init()
  hc@config@global <- S4Vectors::DataFrame(
    min_nodes_number_for_cluster = 1L,
    min_nodes_number_for_network = 1L,
    range_GFC = 2,
    voi = "group"
  )
  counts <- matrix(
    c(10, 12, 9, 8, 7, 11),
    nrow = 6,
    ncol = 1,
    dimnames = list(igraph::V(g)$name, "s1")
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(group = "A", row.names = "s1")
  )
  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(set1 = se)
  )
  hc@config@layer <- S4Vectors::DataFrame(
    layer_id = "set1",
    layer_name = "Layer 1",
    count_source = "set1_counts",
    annotation_source = "set1_anno"
  )
  hc@integration@graph <- g
  hc@integration@gfc <- S4Vectors::DataFrame(
    Gene = igraph::V(g)$name,
    baseline = c(-1.2, -1.1, -0.9, 1.0, 1.1, 1.2),
    stim = c(1.2, 1.1, 0.9, -1.0, -1.1, -1.2)
  )

  tuned <- hc_auto_tune(
    hc,
    tune_cutoff = FALSE,
    tune_clustering = TRUE,
    apply = TRUE,
    cluster_algorithms = c("cluster_louvain"),
    resolution_grid = 0.1,
    module_count_bounds = c(1, 50),
    no_of_iterations = 1,
    verbose = FALSE
  )

  expect_true("cluster_information" %in% names(tuned@integration@cluster))
  ci <- tuned@integration@cluster[["cluster_information"]]
  expect_true(nrow(as.data.frame(ci)) > 0)
  expect_true("auto_tune" %in% names(tuned@satellite))
})
