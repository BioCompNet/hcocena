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
    "hc_plot_enrichment_panels",
    "hc_plot_integrated_network",
    "hc_resolve_paths",
    "hc_auto_set_paths",
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
    "hc_functional_enrichment",
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

test_that("hc_read_data supports object sources with FALSE path settings", {
  count_name <- paste0("hcocena_counts_obj_", as.integer(Sys.time()))
  anno_name <- paste0("hcocena_anno_obj_", as.integer(Sys.time()))

  assign(
    count_name,
    data.frame(
      SYMBOL = c("G1", "G2"),
      sample1 = c(10, 20),
      sample2 = c(15, 25),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    envir = .GlobalEnv
  )
  assign(
    anno_name,
    data.frame(
      SampleID = c("sample1", "sample2"),
      merged = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    envir = .GlobalEnv
  )
  on.exit({
    if (exists(count_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = count_name, envir = .GlobalEnv)
    }
    if (exists(anno_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = anno_name, envir = .GlobalEnv)
    }
  }, add = TRUE)

  hc <- hc_init()
  hc <- hc_set_paths(
    hc,
    dir_count_data = FALSE,
    dir_annotation = FALSE,
    dir_reference_files = tempdir(),
    dir_output = tempdir()
  )
  hc <- hc_define_layers(hc, data_sets = list(Layer = c(count_name, anno_name)))

  expect_error(
    hc <- hc_read_data(
      hc,
      gene_symbol_col = "SYMBOL",
      sample_col = "SampleID",
      count_has_rn = FALSE,
      anno_has_rn = FALSE
    ),
    NA
  )

  legacy <- as_hcobject(hc)
  expect_true(all(c("set1_counts", "set1_anno") %in% names(legacy$data)))
  expect_identical(rownames(legacy$data$set1_counts), c("G1", "G2"))
  expect_identical(colnames(legacy$data$set1_counts), c("sample1", "sample2"))
})

test_that("hc_read_data reports missing object sources clearly with FALSE paths", {
  hc <- hc_init()
  hc <- hc_set_paths(
    hc,
    dir_count_data = FALSE,
    dir_annotation = FALSE,
    dir_reference_files = tempdir(),
    dir_output = tempdir()
  )
  hc <- hc_define_layers(
    hc,
    data_sets = list(Layer = c("definitely_missing_counts", "definitely_missing_anno"))
  )

  expect_error(
    hc_read_data(
      hc,
      gene_symbol_col = "SYMBOL",
      sample_col = "SampleID",
      count_has_rn = FALSE,
      anno_has_rn = FALSE
    ),
    "was not found as an object"
  )
})

test_that("hc_set_paths enforces scalar path-like inputs", {
  hc <- hc_init()
  fancy_path <- structure(tempdir(), class = c("fs_path", "character"))

  hc <- hc_set_paths(
    hc,
    dir_count_data = fancy_path,
    dir_annotation = fancy_path,
    dir_reference_files = fancy_path,
    dir_output = fancy_path
  )

  legacy <- as_hcobject(hc)
  expect_type(legacy$working_directory$dir_count_data, "character")

  expect_error(
    hc_set_paths(
      hc_init(),
      dir_count_data = c("a", "b"),
      dir_annotation = "a",
      dir_reference_files = "a",
      dir_output = "a"
    ),
    "scalar path string or FALSE"
  )
})

test_that("hc_check_dirs accepts FALSE path entries", {
  out_dir <- file.path(tempdir(), paste0("hcocena_out_false_", as.integer(Sys.time())))
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

  expect_error(hc <- hc_check_dirs(hc, create_output_dir = TRUE), NA)
  expect_true(dir.exists(out_dir))
})

test_that("hc_resolve_paths prefers local project folders", {
  root <- file.path(tempdir(), paste0("hcocena_paths_", as.integer(Sys.time())))
  dir.create(file.path(root, "counts_local"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "anno_local"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "ref_local"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "out_local"), recursive = TRUE, showWarnings = FALSE)

  paths <- hc_resolve_paths(
    project_root_hint = root,
    local_count_subdir = "counts_local",
    local_annotation_subdir = "anno_local",
    local_reference_subdir = "ref_local",
    local_output_subdir = "out_local",
    docker_reference_dir = file.path(root, "no_such_docker_reference")
  )

  expect_false(paths$docker_mode)
  expect_true(grepl("/$", paths$dir_count_data))
  expect_true(grepl("/counts_local/$", paths$dir_count_data))
  expect_true(grepl("/anno_local/$", paths$dir_annotation))
  expect_true(grepl("/ref_local/$", paths$dir_reference_files))
  expect_true(grepl("/out_local/$", paths$dir_output))
})

test_that("hc_auto_set_paths applies overrides", {
  root <- file.path(tempdir(), paste0("hcocena_paths_override_", as.integer(Sys.time())))
  dir.create(file.path(root, "reference_files"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "output"), recursive = TRUE, showWarnings = FALSE)

  hc <- hc_init()
  hc <- hc_auto_set_paths(
    hc,
    project_root_hint = root,
    docker_reference_dir = file.path(root, "no_such_docker_reference"),
    dir_count_data = FALSE,
    dir_annotation = FALSE
  )

  legacy <- as_hcobject(hc)
  expect_identical(legacy$working_directory$dir_count_data, FALSE)
  expect_identical(legacy$working_directory$dir_annotation, FALSE)
  expect_true(grepl("/reference_files/$", legacy$working_directory$dir_reference_files))
  expect_true(grepl("/output/$", legacy$working_directory$dir_output))
})

test_that("S4 wrappers do not leak hcobject into .GlobalEnv", {
  had_existing <- exists("hcobject", envir = .GlobalEnv, inherits = FALSE)
  old_hcobject <- NULL
  if (had_existing) {
    old_hcobject <- get("hcobject", envir = .GlobalEnv, inherits = FALSE)
    rm(list = "hcobject", envir = .GlobalEnv)
  }
  on.exit({
    if (had_existing) {
      assign("hcobject", old_hcobject, envir = .GlobalEnv)
    } else if (exists("hcobject", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "hcobject", envir = .GlobalEnv)
    }
  }, add = TRUE)

  hc <- hc_init()
  hc <- hc_set_global_settings(hc, variable_of_interest = "merged")

  expect_false(exists("hcobject", envir = .GlobalEnv, inherits = FALSE))
})

test_that("S4 wrappers do not overwrite an existing global hcobject", {
  had_existing <- exists("hcobject", envir = .GlobalEnv, inherits = FALSE)
  old_hcobject <- NULL
  if (had_existing) {
    old_hcobject <- get("hcobject", envir = .GlobalEnv, inherits = FALSE)
  }
  sentinel <- list(sentinel = TRUE)
  assign("hcobject", sentinel, envir = .GlobalEnv)
  on.exit({
    if (had_existing) {
      assign("hcobject", old_hcobject, envir = .GlobalEnv)
    } else if (exists("hcobject", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "hcobject", envir = .GlobalEnv)
    }
  }, add = TRUE)

  hc <- hc_init()
  hc <- hc_set_global_settings(hc, variable_of_interest = "merged")

  expect_identical(get("hcobject", envir = .GlobalEnv, inherits = FALSE), sentinel)
})

test_that("hc_col_anno_categorical can remove stored annotations", {
  hc <- hc_init()
  hc@satellite <- S4Vectors::SimpleList()
  hc@satellite[["column_annos_categorical"]] <- list(
    cohort = matrix(
      c(1, 0),
      nrow = 1,
      dimnames = list("grpA", c("A", "B"))
    ),
    site = matrix(
      c(0, 1),
      nrow = 1,
      dimnames = list("grpA", c("X", "Y"))
    )
  )

  hc_one <- hc_col_anno_categorical(hc, variables = NULL, variable_label = "cohort")
  expect_true("column_annos_categorical" %in% names(hc_one@satellite))
  expect_false("cohort" %in% names(hc_one@satellite$column_annos_categorical))
  expect_true("site" %in% names(hc_one@satellite$column_annos_categorical))

  hc_all <- hc_col_anno_categorical(hc_one, variables = character(0))
  expect_false("column_annos_categorical" %in% names(hc_all@satellite))
})
