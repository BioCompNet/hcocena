test_that("S4 API config helpers update object", {
  hc <- hc_init()
  hc <- hc_set_paths(
    hc,
    dir_count_data = "counts/",
    dir_annotation = "anno/",
    dir_reference_files = "ref/",
    dir_output = "out/"
  )

  hc <- hc_define_layers(
    hc,
    data_sets = list(
      RNA = c("rna_counts.tsv", "rna_anno.tsv"),
      ARR = c("arr_counts.tsv", "arr_anno.tsv")
    )
  )

  expect_equal(nrow(hc@config@paths), 1)
  expect_equal(nrow(hc@config@layer), 2)
  expect_equal(as.character(hc@config@layer$layer_id), c("set1", "set2"))
})

test_that("extended S4 wrappers are exported", {
  expected <- c(
    "hc_check_dirs",
    "hc_init_save_folder",
    "hc_plot_cutoffs",
    "hc_plot_cluster_heatmap",
    "hc_plot_integrated_network",
    "hc_tf_overrep_module",
    "hc_write_session_info",
    "hc_suggest_topvar",
    "hc_plot_sample_distributions",
    "hc_pca",
    "hc_meta_plot",
    "hc_export_clusters",
    "hc_get_module_scores",
    "hc_algo_alluvial",
    "hc_pca_algo_compare",
    "hc_update_clustering_algorithm",
    "hc_export_to_local_folder",
    "hc_import_layout_from_local_folder",
    "hc_export_to_cytoscape",
    "hc_import_layout_from_cytoscape",
    "hc_plot_gfc_network",
    "hc_tf_overrep_network",
    "hc_check_tf",
    "hc_find_hubs",
    "hc_visualize_gene_expression",
    "hc_highlight_geneset",
    "hc_colour_single_cluster",
    "hc_col_anno_categorical",
    "hc_meta_correlation_cat"
  )

  exported <- getNamespaceExports("hcocena")
  expect_true(all(expected %in% exported))
})

test_that("output directory setup can create missing root and use it directly", {
  out_dir <- file.path(tempdir(), paste0("hcocena_out_", as.integer(Sys.time())))
  if (dir.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }

  hc <- hc_init()
  hc <- hc_set_paths(
    hc,
    dir_count_data = FALSE,
    dir_annotation = FALSE,
    dir_reference_files = tempdir(),
    dir_output = out_dir
  )

  hc <- hc_check_dirs(hc, create_output_dir = TRUE)
  expect_true(dir.exists(out_dir))

  hc <- hc_init_save_folder(hc, name = "", use_output_dir = TRUE)
  legacy <- as_hcobject(hc)
  expect_equal(legacy$global_settings$save_folder, "")
})

test_that("hc_read_data auto-initializes output setup", {
  work_dir <- file.path(tempdir(), paste0("hcocena_read_", as.integer(Sys.time())))
  if (dir.exists(work_dir)) {
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  count_file <- file.path(work_dir, "counts.tsv")
  anno_file <- file.path(work_dir, "anno.tsv")
  out_dir <- file.path(work_dir, "output_root")

  counts <- data.frame(
    gene = c("G1", "G2"),
    sample1 = c(10, 20),
    sample2 = c(15, 25),
    stringsAsFactors = FALSE
  )
  anno <- data.frame(
    sample = c("sample1", "sample2"),
    merged = c("A", "B"),
    stringsAsFactors = FALSE
  )

  utils::write.table(counts, file = count_file, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(anno, file = anno_file, sep = "\t", row.names = FALSE, quote = FALSE)

  hc <- hc_init()
  hc <- hc_set_paths(
    hc,
    dir_count_data = paste0(work_dir, "/"),
    dir_annotation = paste0(work_dir, "/"),
    dir_reference_files = paste0(work_dir, "/"),
    dir_output = out_dir
  )
  hc <- hc_define_layers(hc, data_sets = list(Layer = c("counts.tsv", "anno.tsv")))

  hc <- hc_read_data(
    hc,
    gene_symbol_col = "gene",
    sample_col = "sample",
    count_has_rn = FALSE,
    anno_has_rn = FALSE
  )

  legacy <- as_hcobject(hc)
  expect_true(dir.exists(out_dir))
  expect_identical(legacy$global_settings$save_folder, "")
  expect_true(all(c("set1_counts", "set1_anno") %in% names(legacy$data)))

  unlink(work_dir, recursive = TRUE, force = TRUE)
})
