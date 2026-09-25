.make_regroup_hc <- function(two_layers = FALSE) {
  expression <- matrix(
    c(
      12, 11, 2, 2,
      10, 12, 2, 3,
      2, 2, 11, 12,
      3, 2, 12, 10
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("g", 1:4),
      paste0("s", 1:4)
    )
  )
  anno <- S4Vectors::DataFrame(
    group = c("control", "control", "case", "case"),
    row.names = colnames(expression)
  )
  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = expression),
    colData = anno
  )
  experiments <- S4Vectors::SimpleList(set1 = se1)

  layer_cfg <- S4Vectors::DataFrame(
    layer_id = "set1",
    layer_name = "Layer 1",
    count_source = "set1_counts",
    annotation_source = "set1_anno"
  )
  if (isTRUE(two_layers)) {
    expression2 <- expression
    colnames(expression2) <- paste0("t", 1:4)
    anno2 <- anno
    rownames(anno2) <- colnames(expression2)
    se2 <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = expression2),
      colData = anno2
    )
    experiments[["set2"]] <- se2
    layer_cfg <- S4Vectors::DataFrame(
      layer_id = c("set1", "set2"),
      layer_name = c("Layer 1", "Layer 2"),
      count_source = c("set1_counts", "set2_counts"),
      annotation_source = c("set1_anno", "set2_anno")
    )
  }

  hc <- hc_init()
  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = experiments)
  hc@config@layer <- layer_cfg
  hc@config@global <- S4Vectors::DataFrame(
    voi = "group",
    control = "control",
    range_GFC = 2,
    data_in_log = FALSE
  )
  hc@integration@cluster <- S4Vectors::SimpleList(
    cluster_information = data.frame(
      color = c("red", "blue"),
      gene_n = c("g1,g2", "g3,g4"),
      gene_no = c(2L, 2L),
      cluster_included = c("yes", "yes"),
      stringsAsFactors = FALSE
    ),
    module_label_map = c(red = "M1", blue = "M2"),
    heatmap_row_order = c("red", "blue"),
    gfc_colors = c("navy", "white", "firebrick"),
    gfc_scale_limits = c(-2, 2)
  )
  methods::validObject(hc)
  hc
}

test_that("hc_regroup_samples previews sample clusters without mutating hc", {
  hc <- .make_regroup_hc()
  before <- serialize(hc, NULL)

  preview <- hc_regroup_samples(
    hc,
    layer = "set1",
    k = 2,
    silent = TRUE
  )

  expect_s3_class(preview, "hc_sample_regrouping")
  expect_false(preview$applied)
  expect_identical(serialize(hc, NULL), before)
  expect_identical(serialize(preview$hc, NULL), before)
  expect_equal(dim(preview$results$set1$sample_gfc), c(2, 4))
  expect_identical(rownames(preview$results$set1$sample_gfc), c("M1", "M2"))
  expect_s3_class(preview$results$set1$tree, "hclust")
  expect_s3_class(preview$results$set1$plot, "pheatmap")

  groups <- stats::setNames(
    preview$results$set1$clusters$cluster,
    preview$results$set1$clusters$sample
  )
  expect_identical(groups[["s1"]], groups[["s2"]])
  expect_identical(groups[["s3"]], groups[["s4"]])
  expect_false(groups[["s1"]] == groups[["s3"]])

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  expect_invisible(plot(preview))
  grDevices::dev.off()
  expect_gt(file.info(plot_file)$size, 0)
  expect_error(plot(preview, layer = "missing"), "Unknown regrouping plot layer")
})

test_that("hc_regroup_samples applies a new annotation column without replacing VOI", {
  hc <- .make_regroup_hc()
  original_group <- as.character(
    SummarizedExperiment::colData(
      MultiAssayExperiment::experiments(hc@mae)[["set1"]]
    )$group
  )

  applied <- hc_regroup_samples(
    hc,
    layer = 1,
    k = 2,
    group_col = "module_profile_group",
    apply = TRUE,
    silent = TRUE
  )
  out_anno <- SummarizedExperiment::colData(
    MultiAssayExperiment::experiments(applied$hc@mae)[["set1"]]
  )
  input_anno <- SummarizedExperiment::colData(
    MultiAssayExperiment::experiments(hc@mae)[["set1"]]
  )

  expect_true(applied$applied)
  expect_true("module_profile_group" %in% colnames(out_anno))
  expect_false("module_profile_group" %in% colnames(input_anno))
  expect_identical(as.character(out_anno$group), original_group)
  expect_identical(as.character(applied$hc@config@global$voi[[1]]), "group")
  expect_true("sample_regrouping" %in% names(applied$hc@satellite))
  expect_identical(
    applied$hc@satellite$sample_regrouping$parameters$group_col,
    "module_profile_group"
  )
  expect_true(methods::validObject(applied$hc))
})

test_that("hc_regroup_samples resolves split labels and selected modules", {
  hc <- .make_regroup_hc()
  cluster_calc <- as.list(hc@integration@cluster)
  cluster_calc$module_label_map <- c(red = "M1.1", blue = "M2")
  hc@integration@cluster <- S4Vectors::SimpleList(cluster_calc)

  preview <- hc_regroup_samples(
    hc,
    modules = "M1.1",
    k = 2,
    silent = TRUE
  )

  expect_identical(rownames(preview$results$set1$sample_gfc), "M1.1")
  expect_identical(preview$results$set1$module_table$color, "red")
  expect_identical(preview$parameters$modules, "M1.1")
})

test_that("hc_regroup_samples supports layer-specific k values", {
  hc <- .make_regroup_hc(two_layers = TRUE)
  preview <- hc_regroup_samples(
    hc,
    layer = "all",
    k = c("Layer 1" = 2, set2 = 3),
    silent = TRUE
  )

  expect_identical(names(preview$results), c("set1", "set2"))
  expect_equal(length(unique(preview$results$set1$clusters$cluster)), 2)
  expect_equal(length(unique(preview$results$set2$clusters$cluster)), 3)
})

test_that("hc_regroup_samples validates destructive and reference choices", {
  hc <- .make_regroup_hc()
  se <- MultiAssayExperiment::experiments(hc@mae)[["set1"]]
  anno <- SummarizedExperiment::colData(se)
  anno$hc_regroup <- "existing"
  SummarizedExperiment::colData(se) <- anno
  hc@mae[["set1"]] <- se

  expect_error(
    hc_regroup_samples(hc, apply = TRUE, silent = TRUE),
    "already exists"
  )
  expect_error(
    hc_regroup_samples(hc, k = 5, silent = TRUE),
    "cannot exceed"
  )

  control_preview <- hc_regroup_samples(
    hc,
    reference = "control",
    k = 2,
    silent = TRUE
  )
  expect_equal(dim(control_preview$results$set1$sample_gfc), c(2, 4))
})
