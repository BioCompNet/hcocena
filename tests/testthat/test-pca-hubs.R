.make_pca_hub_legacy <- function(single_annotation = FALSE,
                                 include_tf_reference = FALSE) {
  anno1 <- data.frame(group = c("A", "B"), row.names = c("s1", "s2"))
  anno2 <- data.frame(group = c("A", "B"), row.names = c("s3", "s4"))
  if (!isTRUE(single_annotation)) {
    anno1$batch <- c("x", "y")
    anno2$batch <- c("x", "y")
  }

  supplementary_data <- list()
  if (isTRUE(include_tf_reference)) {
    supplementary_data$TF <- data.frame(
      human = c("g1", "g2"),
      category = c("Innate", "Adaptive"),
      stringsAsFactors = FALSE
    )
  }

  list(
    working_directory = list(dir_output = tempdir()),
    global_settings = list(
      voi = "group",
      save_folder = "",
      organism = "human"
    ),
    layers = list(
      set1 = c("set1_counts", "set1_anno"),
      set2 = c("set2_counts", "set2_anno")
    ),
    layers_names = c("Layer1", "Layer2"),
    data = list(
      set1_counts = matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("g1", "g2"), c("s1", "s2"))
      ),
      set1_anno = anno1,
      set2_counts = matrix(
        c(5, 6),
        nrow = 1,
        dimnames = list("g2", c("s3", "s4"))
      ),
      set2_anno = anno2
    ),
    integrated_output = list(
      cluster_calc = list(
        cluster_information = data.frame(
          gene_n = "g1,g2",
          color = "red",
          stringsAsFactors = FALSE
        )
      )
    ),
    satellite_outputs = list(),
    supplementary_data = supplementary_data,
    layer_settings = list(),
    supplement = list(),
    cutoff_vec = NULL
  )
}

.with_pca_hub_bridge <- function(hcobject, code) {
  state <- hcocena:::.hc_bind_bridge_hcobject(
    hcobject,
    envo = asNamespace("hcocena")
  )
  on.exit(hcocena:::.hc_restore_bridge_hcobject(state), add = TRUE)
  force(code)
}

test_that("PCA accepts sparse topvar expression matrices", {
  set.seed(1)
  dense <- matrix(
    stats::rnorm(30),
    nrow = 5,
    dimnames = list(paste0("g", 1:5), paste0("s", 1:6))
  )
  sparse <- Matrix::Matrix(dense, sparse = TRUE)

  prepared <- hcocena:::.hc_pca_prepare_expression(
    sparse,
    layer_label = "sparse",
    scale_features = TRUE
  )

  expect_true(is.matrix(prepared))
  expect_type(prepared, "double")
  expect_equal(dim(prepared), c(6, 5))
  expect_s3_class(stats::prcomp(prepared, scale. = TRUE), "prcomp")
})

test_that("PCA removes constant features before scaling", {
  expression <- rbind(
    variable = c(1, 2, 4, 8),
    constant = c(3, 3, 3, 3)
  )

  prepared <- NULL
  expect_warning(
    prepared <- hcocena:::.hc_pca_prepare_expression(
      expression,
      layer_label = "constant",
      scale_features = TRUE
    ),
    "removed 1 constant feature"
  )

  expect_equal(dim(prepared), c(4, 1))
  expect_s3_class(stats::prcomp(prepared, scale. = TRUE), "prcomp")
})

test_that("PCA plotting supports a single principal component", {
  expression <- matrix(
    c(1, 2, 4, 8),
    ncol = 1,
    dimnames = list(paste0("s", 1:4), "g1")
  )
  pca <- stats::prcomp(expression, scale. = FALSE)

  plot <- hcocena:::.hc_pca_individual_plot(
    res.pca = pca,
    groups = c("A", "A", "B", "B"),
    palette = c("#0072B2", "#D55E00"),
    ellipses = FALSE,
    title = "One-component PCA"
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(unique(plot$data$PC2), 0)
})

test_that("public topvar PCA works with a sparse stored matrix", {
  fixture <- system.file(
    "extdata",
    "hc_after_part1.rds",
    package = "hcocena"
  )
  skip_if(!nzchar(fixture), "The after-part1 fixture is unavailable.")

  hc <- readRDS(fixture)
  legacy <- hcocena:::as_hcobject(hc)
  legacy$layer_specific_outputs$set1$part1$topvar <- Matrix::Matrix(
    legacy$layer_specific_outputs$set1$part1$topvar,
    sparse = TRUE
  )
  output_dir <- tempfile("hc-pca-sparse-")
  dir.create(output_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  legacy$working_directory$dir_output <- paste0(
    output_dir,
    .Platform$file.sep
  )
  legacy$global_settings$save_folder <- ""
  hc <- hcocena:::as_hcocena(legacy)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- hc_pca(
    hc,
    which = "topvar",
    color_by = "none",
    ellipses = FALSE
  )

  expect_s4_class(result, "HCoCenaExperiment")
  expect_true("pca" %in% names(as.list(result@satellite)))
  expect_true(file.exists(file.path(output_dir, "PCA_topvar_Layer1.pdf")))
  expect_true(file.exists(file.path(output_dir, "PCA_topvar_Layer1.png")))
})

test_that("hub expression skips genes absent from individual layers", {
  legacy <- .make_pca_hub_legacy()

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_message(
    expect_no_error(
      .with_pca_hub_bridge(
        legacy,
        hcocena:::.hc_visualize_gene_expression_driver(
          genes = "g1",
          name = "hub_probe",
          save = FALSE
        )
      )
    ),
    "absent from layer `Layer2`"
  )
})

test_that("hub expression preserves one-column annotations", {
  legacy <- .make_pca_hub_legacy(single_annotation = TRUE)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(
    .with_pca_hub_bridge(
      legacy,
      hcocena:::.hc_visualize_gene_expression_driver(
        genes = "g1",
        name = "hub_probe",
        save = FALSE
      )
    )
  )
})

test_that("TF_only validates filters and references", {
  expect_identical(
    formals(hcocena:::.hc_find_hubs_driver)$TF_only,
    FALSE
  )
  expect_identical(
    formals(hcocena:::get_hub_nodes)$TF_only,
    FALSE
  )

  no_reference <- .make_pca_hub_legacy()
  expect_null(
    .with_pca_hub_bridge(
      no_reference,
      hcocena:::.hc_hub_tf_filter_genes(FALSE)
    )
  )
  expect_error(
    .with_pca_hub_bridge(
      no_reference,
      hcocena:::.hc_hub_tf_filter_genes("all")
    ),
    "requires a non-empty transcription-factor reference"
  )
  expect_error(
    .with_pca_hub_bridge(
      no_reference,
      hcocena:::.hc_hub_tf_filter_genes(TRUE)
    ),
    "must be FALSE"
  )

  with_reference <- .make_pca_hub_legacy(include_tf_reference = TRUE)
  expect_equal(
    .with_pca_hub_bridge(
      with_reference,
      hcocena:::.hc_hub_tf_filter_genes("all")
    ),
    c("g1", "g2")
  )
  expect_equal(
    .with_pca_hub_bridge(
      with_reference,
      hcocena:::.hc_hub_tf_filter_genes("Innate")
    ),
    "g1"
  )
  expect_error(
    .with_pca_hub_bridge(
      with_reference,
      hcocena:::.hc_hub_tf_filter_genes("Unknown")
    ),
    "Unknown `TF_only` category"
  )
})
