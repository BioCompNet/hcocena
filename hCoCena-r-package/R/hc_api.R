#' Initialize an empty `HCoCenaExperiment`
#'
#' @return A valid, empty `HCoCenaExperiment`.
#' @export
hc_init <- function() {
  new("HCoCenaExperiment")
}

#' Set input/output paths in an `HCoCenaExperiment`
#'
#' @param hc A `HCoCenaExperiment`.
#' @param dir_count_data Path to count data.
#' @param dir_annotation Path to annotation data.
#' @param dir_reference_files Path to reference files.
#' @param dir_output Output directory.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_set_paths <- function(hc, dir_count_data, dir_annotation, dir_reference_files, dir_output) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  hc@config@paths <- S4Vectors::DataFrame(
    dir_count_data = dir_count_data,
    dir_annotation = dir_annotation,
    dir_reference_files = dir_reference_files,
    dir_output = dir_output
  )
  methods::validObject(hc)
  hc
}

#' Define layers in an `HCoCenaExperiment`
#'
#' @param hc A `HCoCenaExperiment`.
#' @param data_sets Named list of layers, each value a pair (count, annotation).
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_define_layers <- function(hc, data_sets = list()) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (length(data_sets) == 0) {
    hc@config@layer <- S4Vectors::DataFrame()
    methods::validObject(hc)
    return(hc)
  }

  layer_names <- names(data_sets)
  if (is.null(layer_names)) {
    stop("`data_sets` must be a named list.")
  }

  rows <- vector("list", length(data_sets))
  for (i in seq_along(data_sets)) {
    pair <- data_sets[[i]]
    rows[[i]] <- list(
      layer_id = paste0("set", i),
      layer_name = as.character(layer_names[[i]]),
      count_source = if (length(pair) >= 1) as.character(pair[[1]]) else NA_character_,
      annotation_source = if (length(pair) >= 2) as.character(pair[[2]]) else NA_character_
    )
  }
  hc@config@layer <- .hc_rows_to_data_frame(rows)
  methods::validObject(hc)
  hc
}

#' Read expression and annotation data (S4 API)
#'
#' @inheritParams read_data
#' @param hc A `HCoCenaExperiment`.
#' @param auto_setup_output Boolean. If `TRUE`, run `hc_check_dirs()` and
#'   initialize the save folder before reading data (if `dir_output` is set).
#' @param create_output_dir Boolean passed to `hc_check_dirs()` when
#'   `auto_setup_output = TRUE`.
#' @param project_folder Optional save subfolder name. Use `""` to write
#'   directly into `dir_output`. If `NULL`, keep an already configured
#'   `save_folder` or default to `""`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_read_data <- function(hc,
                         sep_counts = "\t",
                         sep_anno = "\t",
                         gene_symbol_col = NULL,
                         sample_col = NULL,
                         count_has_rn = TRUE,
                         anno_has_rn = TRUE,
                         auto_setup_output = TRUE,
                         create_output_dir = TRUE,
                         project_folder = NULL) {
  if (isTRUE(auto_setup_output)) {
    has_output_dir <- FALSE
    paths_cfg <- hc@config@paths
    if (base::nrow(paths_cfg) > 0 && "dir_output" %in% base::colnames(paths_cfg)) {
      out_dir <- as.character(paths_cfg$dir_output[[1]])
      has_output_dir <- (
        base::length(out_dir) == 1 &&
          !is.na(out_dir) &&
          base::nzchar(out_dir) &&
          !identical(out_dir, "FALSE")
      )
    }

    if (isTRUE(has_output_dir)) {
      hc <- hc_check_dirs(hc, create_output_dir = create_output_dir)

      save_folder <- project_folder
      if (is.null(save_folder)) {
        global_cfg <- hc@config@global
        if (base::nrow(global_cfg) > 0 && "save_folder" %in% base::colnames(global_cfg)) {
          current <- as.character(global_cfg$save_folder[[1]])
          if (base::length(current) == 1 && !is.na(current)) {
            save_folder <- current
          }
        }
      }

      if (is.null(save_folder) || base::length(save_folder) == 0 || is.na(save_folder)) {
        save_folder <- ""
      }
      save_folder <- as.character(save_folder[[1]])

      hc <- hc_init_save_folder(
        hc,
        name = save_folder,
        use_output_dir = identical(save_folder, "")
      )
    }
  }

  .hc_run_legacy(
    hc = hc,
    fun = "read_data",
    sep_counts = sep_counts,
    sep_anno = sep_anno,
    gene_symbol_col = gene_symbol_col,
    sample_col = sample_col,
    count_has_rn = count_has_rn,
    anno_has_rn = anno_has_rn
  )
}

#' Set global settings (S4 API)
#'
#' @inheritParams set_global_settings
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_set_global_settings <- function(hc,
                                   organism,
                                   control_keyword,
                                   variable_of_interest,
                                   min_nodes_number_for_network = 15,
                                   min_nodes_number_for_cluster = 15,
                                   range_GFC = 2.0,
                                   layout_algorithm = "layout_with_fr",
                                   data_in_log) {
  .hc_run_legacy(
    hc = hc,
    fun = "set_global_settings",
    organism = organism,
    control_keyword = control_keyword,
    variable_of_interest = variable_of_interest,
    min_nodes_number_for_network = min_nodes_number_for_network,
    min_nodes_number_for_cluster = min_nodes_number_for_cluster,
    range_GFC = range_GFC,
    layout_algorithm = layout_algorithm,
    data_in_log = data_in_log
  )
}

#' Set layer settings (S4 API)
#'
#' @inheritParams set_layer_settings
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_set_layer_settings <- function(hc,
                                  top_var,
                                  min_corr = 0.7,
                                  range_cutoff_length,
                                  print_distribution_plots = FALSE) {
  .hc_run_legacy(
    hc = hc,
    fun = "set_layer_settings",
    top_var = top_var,
    min_corr = min_corr,
    range_cutoff_length = range_cutoff_length,
    print_distribution_plots = print_distribution_plots
  )
}

#' Register supplementary references (S4 API)
#'
#' @inheritParams set_supp_files
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_set_supp_files <- function(hc, Tf = NULL, Hallmark = NULL, Go = NULL, Kegg = NULL, Reactome = NULL, ...) {
  .hc_run_legacy(
    hc = hc,
    fun = "set_supp_files",
    Tf = Tf,
    Hallmark = Hallmark,
    Go = Go,
    Kegg = Kegg,
    Reactome = Reactome,
    ...
  )
}

#' Load supplementary references (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_read_supplementary <- function(hc) {
  .hc_run_legacy(
    hc = hc,
    fun = "read_supplementary"
  )
}

#' Run expression analysis part I (S4 API)
#'
#' @inheritParams run_expression_analysis_1
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_run_expression_analysis_1 <- function(hc,
                                         padj = "none",
                                         export = FALSE,
                                         import = NULL,
                                         bayes = FALSE,
                                         prior = 2,
                                         alpha = 0.5,
                                         corr_method = "pearson") {
  .hc_run_legacy(
    hc = hc,
    fun = "run_expression_analysis_1",
    padj = padj,
    export = export,
    import = import,
    bayes = bayes,
    prior = prior,
    alpha = alpha,
    corr_method = corr_method
  )
}

#' Set network cutoffs (S4 API)
#'
#' @inheritParams set_cutoff
#' @param hc A `HCoCenaExperiment`.
#' @param auto Logical; if `TRUE`, automatically select cutoffs from available
#'   tuning outputs in this priority:
#'   1) `hc@satellite$cutoff_tuning$applied_cutoff_vector`,
#'   2) `hc@satellite$cutoff_tuning$recommended_cutoff_vector`,
#'   3) `hc@satellite$auto_tune$cutoff$recommended_cutoff_vector`,
#'   4) internal simple auto cutoff extraction from layer results,
#'   5) existing `hc@config@layer$cutoff`,
#'   6) user-provided `cutoff_vector`,
#'   7) `fallback_cutoff`.
#' @param fallback_cutoff Numeric fallback cutoff used for still-missing layers
#'   (for both `auto = TRUE` and `auto = FALSE`).
#' @param verbose Logical; print additional source details.
#'   The final applied cutoff vector is always printed.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_set_cutoff <- function(hc,
                          cutoff_vector = base::c(),
                          auto = FALSE,
                          fallback_cutoff = 0.982,
                          verbose = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!is.logical(auto) || length(auto) != 1 || is.na(auto)) {
    stop("`auto` must be TRUE or FALSE.")
  }
  if (!is.numeric(fallback_cutoff) || length(fallback_cutoff) != 1 || !is.finite(fallback_cutoff)) {
    stop("`fallback_cutoff` must be a finite numeric scalar.")
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  layer_ids <- character(0)
  if (base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer)) {
    layer_ids <- as.character(hc@config@layer$layer_id)
  } else if (base::length(hc@layer_results) > 0) {
    layer_ids <- base::names(hc@layer_results)
  }

  n_layers <- base::length(layer_ids)
  if (n_layers == 0 && base::nrow(hc@config@layer) > 0) {
    n_layers <- base::nrow(hc@config@layer)
    layer_ids <- paste0("set", seq_len(n_layers))
  }
  if (n_layers == 0 && base::length(hc@layer_results) > 0) {
    n_layers <- base::length(hc@layer_results)
    layer_ids <- if (!is.null(base::names(hc@layer_results))) {
      base::names(hc@layer_results)
    } else {
      paste0("set", seq_len(n_layers))
    }
  }
  if (n_layers == 0 && base::length(cutoff_vector) > 0) {
    n_layers <- base::length(cutoff_vector)
    layer_ids <- paste0("set", seq_len(n_layers))
  }
  if (n_layers == 0) {
    stop("No layers found to set cutoffs for.")
  }

  align_vec <- function(x) {
    xv <- suppressWarnings(as.numeric(x))
    out <- rep(NA_real_, n_layers)
    if (length(xv) == 0) {
      return(out)
    }

    nms <- names(x)
    if (!is.null(nms) && length(nms) == length(x) && any(nzchar(nms))) {
      idx <- match(layer_ids, nms)
      ok <- is.finite(idx)
      if (any(ok)) {
        out[ok] <- suppressWarnings(as.numeric(x[idx[ok]]))
      }
    }

    miss <- which(!is.finite(out))
    seq_vals <- xv[is.finite(xv)]
    if (length(miss) > 0 && length(seq_vals) > 0) {
      n_copy <- min(length(miss), length(seq_vals))
      out[miss[seq_len(n_copy)]] <- seq_vals[seq_len(n_copy)]
    }
    out
  }

  chosen <- rep(NA_real_, n_layers)
  source_per_layer <- rep(NA_character_, n_layers)

  fill_missing <- function(vec, source_name) {
    aligned <- align_vec(vec)
    idx <- which(!is.finite(chosen) & is.finite(aligned))
    if (length(idx) > 0) {
      chosen[idx] <<- aligned[idx]
      source_per_layer[idx] <<- source_name
    }
  }

  if (isTRUE(auto)) {
    sat <- as.list(hc@satellite)
    if ("cutoff_tuning" %in% names(sat) && is.list(sat[["cutoff_tuning"]])) {
      ct <- sat[["cutoff_tuning"]]
      if ("applied_cutoff_vector" %in% names(ct)) {
        fill_missing(ct[["applied_cutoff_vector"]], "cutoff_tuning.applied")
      }
      if ("recommended_cutoff_vector" %in% names(ct)) {
        fill_missing(ct[["recommended_cutoff_vector"]], "cutoff_tuning.recommended")
      }
    }
    if ("auto_tune" %in% names(sat) &&
        is.list(sat[["auto_tune"]]) &&
        "cutoff" %in% names(sat[["auto_tune"]]) &&
        is.list(sat[["auto_tune"]][["cutoff"]]) &&
        "recommended_cutoff_vector" %in% names(sat[["auto_tune"]][["cutoff"]])) {
      fill_missing(sat[["auto_tune"]][["cutoff"]][["recommended_cutoff_vector"]], "auto_tune.recommended")
    }

    if (base::exists(".hc_auto_collect_cutoffs", mode = "function", inherits = TRUE)) {
      simple_auto <- tryCatch(.hc_auto_collect_cutoffs(hc), error = function(e) numeric(0))
      fill_missing(simple_auto, "simple_auto")
    }

    if (base::nrow(hc@config@layer) > 0 && "cutoff" %in% base::colnames(hc@config@layer)) {
      fill_missing(hc@config@layer$cutoff, "config.layer.cutoff")
    }

    fill_missing(cutoff_vector, "user.cutoff_vector")
  } else {
    fill_missing(cutoff_vector, "user.cutoff_vector")
    if (base::nrow(hc@config@layer) > 0 && "cutoff" %in% base::colnames(hc@config@layer)) {
      fill_missing(hc@config@layer$cutoff, "config.layer.cutoff")
    }
  }

  miss <- which(!is.finite(chosen))
  if (length(miss) > 0) {
    chosen[miss] <- as.numeric(fallback_cutoff)
    source_per_layer[miss] <- "fallback"
  }

  message(
    "hc_set_cutoff(auto=", if (isTRUE(auto)) "TRUE" else "FALSE", "): applying cutoffs = ",
    paste0(format(chosen, digits = 4), collapse = ", ")
  )
  if (isTRUE(verbose)) {
    src_df <- data.frame(
      layer_id = layer_ids,
      source = as.character(source_per_layer),
      stringsAsFactors = FALSE
    )
    message(
      "hc_set_cutoff(): source per layer -> ",
      paste0(src_df$layer_id, ":", src_df$source, collapse = ", ")
    )
  }

  hc <- .hc_run_legacy(
    hc = hc,
    fun = "set_cutoff",
    cutoff_vector = chosen
  )

  sat2 <- as.list(hc@satellite)
  sat2[["cutoff_selection"]] <- list(
    created_at = as.character(base::Sys.time()),
    mode = if (isTRUE(auto)) "auto" else "manual",
    fallback_cutoff = as.numeric(fallback_cutoff),
    applied_cutoff_vector = chosen,
    layer_source = data.frame(
      layer_id = layer_ids,
      cutoff = as.numeric(chosen),
      source = as.character(source_per_layer),
      stringsAsFactors = FALSE
    )
  )
  hc@satellite <- S4Vectors::SimpleList(sat2)
  methods::validObject(hc)
  hc
}

#' Run expression analysis part II (S4 API)
#'
#' @inheritParams run_expression_analysis_2
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_run_expression_analysis_2 <- function(hc,
                                         grouping_v = NULL,
                                         plot_HM = TRUE,
                                         method = "complete",
                                         additional_anno = NULL,
                                         cols = NULL) {
  .hc_run_legacy(
    hc = hc,
    fun = "run_expression_analysis_2",
    grouping_v = grouping_v,
    plot_HM = plot_HM,
    method = method,
    additional_anno = additional_anno,
    cols = cols
  )
}

#' Build integrated network (S4 API)
#'
#' @inheritParams build_integrated_network
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_build_integrated_network <- function(hc,
                                        mode = "u",
                                        with = NULL,
                                        multi_edges = "min",
                                        GFC_when_missing = NULL) {
  if (is.null(GFC_when_missing) &&
      base::nrow(hc@config@global) > 0 &&
      "range_GFC" %in% base::colnames(hc@config@global)) {
    GFC_when_missing <- -hc@config@global$range_GFC[[1]]
  }
  if (is.null(GFC_when_missing)) {
    GFC_when_missing <- -2.0
  }
  .hc_run_legacy(
    hc = hc,
    fun = "build_integrated_network",
    mode = mode,
    with = with,
    multi_edges = multi_edges,
    GFC_when_missing = GFC_when_missing
  )
}

#' Cluster integrated network (S4 API)
#'
#' @inheritParams cluster_calculation
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_cluster_calculation <- function(hc,
                                   cluster_algo = "cluster_leiden",
                                   no_of_iterations = 2,
                                   resolution = 0.1,
                                   partition_type = "RBConfigurationVertexPartition",
                                   max_cluster_count_per_gene = 1,
                                   return_result = FALSE) {
  .hc_run_legacy(
    hc = hc,
    fun = "cluster_calculation",
    cluster_algo = cluster_algo,
    no_of_iterations = no_of_iterations,
    resolution = resolution,
    partition_type = partition_type,
    max_cluster_count_per_gene = max_cluster_count_per_gene,
    return_result = return_result
  )
}

#' Merge modules based on module-heatmap similarity (S4 API)
#'
#' @inheritParams merge_clusters
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_merge_clusters <- function(hc,
                              k = "auto",
                              save = TRUE,
                              method = "complete",
                              k_min = 2,
                              k_max = NULL,
                              auto_parsimony_penalty = 1e-04,
                              verbose = TRUE) {
  .hc_run_legacy(
    hc = hc,
    fun = "merge_clusters",
    k = k,
    save = save,
    method = method,
    k_min = k_min,
    k_max = k_max,
    auto_parsimony_penalty = auto_parsimony_penalty,
    verbose = verbose
  )
}

#' Split one or multiple modules into submodules (S4 API)
#'
#' Re-clusters genes within selected modules and replaces those modules by
#' submodules (with color shades of the parent module).
#'
#' @param hc A `HCoCenaExperiment`.
#' @param modules Character/numeric vector of modules to split. Accepts module
#'   labels (e.g. `"M3"`), module colors, or module indices.
#' @inheritParams split_modules
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_split_modules <- function(hc,
                             modules,
                             cluster_algo = "cluster_leiden",
                             no_of_iterations = 2,
                             resolution = 0.1,
                             resolution_grid = NULL,
                             resolution_test_only = FALSE,
                             partition_type = "RBConfigurationVertexPartition",
                             seed = 168575,
                             drop_small_submodules = TRUE,
                             min_submodule_size = NULL,
                             min_module_size = NULL,
                             verbose = TRUE) {
  .hc_run_legacy(
    hc = hc,
    fun = "split_modules",
    modules = modules,
    cluster_algo = cluster_algo,
    no_of_iterations = no_of_iterations,
    resolution = resolution,
    resolution_grid = resolution_grid,
    resolution_test_only = resolution_test_only,
    partition_type = partition_type,
    seed = seed,
    drop_small_submodules = drop_small_submodules,
    min_submodule_size = min_submodule_size,
    min_module_size = min_module_size,
    verbose = verbose
  )
}

#' Undo module splitting (S4 API)
#'
#' Restores module structure from split history.
#'
#' @param hc A `HCoCenaExperiment`.
#' @inheritParams unsplit_modules
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_unsplit_modules <- function(hc,
                               which = c("last", "all"),
                               verbose = TRUE) {
  which <- base::match.arg(which)
  .hc_run_legacy(
    hc = hc,
    fun = "unsplit_modules",
    which = which,
    verbose = verbose
  )
}

#' Functional enrichment (S4 API)
#'
#' @inheritParams functional_enrichment
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_functional_enrichment <- function(hc,
                                     gene_sets = c("Go", "Kegg", "Hallmark", "Reactome"),
                                     top = 5,
                                     clusters = c("all"),
                                     padj = "BH",
                                     qval = 0.05,
                                     heatmap_side = "left",
                                     heatmap_cluster_rows = FALSE,
                                     heatmap_cluster_columns = FALSE,
                                     heatmap_show_row_dend = FALSE,
                                     heatmap_show_column_dend = FALSE,
                                     heatmap_order = NULL,
                                     heatmap_module_label_mode = "same",
                                     heatmap_show_gene_counts = FALSE,
                                     ...) {
  .hc_run_legacy(
    hc = hc,
    fun = "functional_enrichment",
    gene_sets = gene_sets,
    top = top,
    clusters = clusters,
    padj = padj,
    qval = qval,
    heatmap_side = heatmap_side,
    heatmap_cluster_rows = heatmap_cluster_rows,
    heatmap_cluster_columns = heatmap_cluster_columns,
    heatmap_show_row_dend = heatmap_show_row_dend,
    heatmap_show_column_dend = heatmap_show_column_dend,
    heatmap_order = heatmap_order,
    heatmap_module_label_mode = heatmap_module_label_mode,
    heatmap_show_gene_counts = heatmap_show_gene_counts,
    ...
  )
}

#' Upstream regulator/pathway inference (S4 API)
#'
#' @inheritParams upstream_inference
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_upstream_inference <- function(hc,
                                  resources = c("TF", "Pathway"),
                                  top = 5,
                                  clusters = c("all"),
                                  padj = "BH",
                                  qval = 0.05,
                                  tf_confidence = c("A", "B", "C"),
                                  minsize = 5,
                                  method = "ulm",
                                  activity_input = "gfc",
                                  fc_comparisons = NULL,
                                  heatmap_side = "left",
                                  plot = TRUE,
                                  save_pdf = TRUE,
                                  plot_per_comparison = TRUE,
                                  consistent_terms = FALSE,
                                  overall_plot_scale = 1) {
  .hc_run_legacy(
    hc = hc,
    fun = "upstream_inference",
    resources = resources,
    top = top,
    clusters = clusters,
    padj = padj,
    qval = qval,
    tf_confidence = tf_confidence,
    minsize = minsize,
    method = method,
    activity_input = activity_input,
    fc_comparisons = fc_comparisons,
    heatmap_side = heatmap_side,
    plot = plot,
    save_pdf = save_pdf,
    plot_per_comparison = plot_per_comparison,
    consistent_terms = consistent_terms,
    overall_plot_scale = overall_plot_scale
  )
}

#' Module cell-type annotation from Enrichr (S4 API)
#'
#' @inheritParams celltype_annotation
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_celltype_annotation <- function(hc,
                                   databases = c("Descartes_Cell_Types_and_Tissue_2021", "Human_Gene_Atlas"),
                                   clusters = c("all"),
                                   mode = c("coarse", "fine"),
                                   top = 3,
                                   qval = 0.1,
                                   padj = "BH",
                                   min_term_genes = 5,
                                   min_gs_size = 10,
                                   max_gs_size = 5000,
                                   annotation_slot = c("auto", "enriched_per_cluster", "enriched_per_cluster2"),
                                   coarse_map = NULL,
                                   coarse_include_other = TRUE,
                                   refresh_db = FALSE,
                                   export_excel = TRUE,
                                   excel_file = "Module_Celltype_Annotation.xlsx",
                                   plot_heatmap = FALSE,
                                   heatmap_file_name = "Heatmap_modules_celltype_annotation.pdf",
                                   ...) {
  .hc_run_legacy(
    hc = hc,
    fun = "celltype_annotation",
    databases = databases,
    clusters = clusters,
    mode = mode,
    top = top,
    qval = qval,
    padj = padj,
    min_term_genes = min_term_genes,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    annotation_slot = annotation_slot,
    coarse_map = coarse_map,
    coarse_include_other = coarse_include_other,
    refresh_db = refresh_db,
    export_excel = export_excel,
    excel_file = excel_file,
    plot_heatmap = plot_heatmap,
    heatmap_file_name = heatmap_file_name,
    ...
  )
}

#' Plot module knowledge network from enrichment + upstream inference (S4 API)
#'
#' @inheritParams plot_enrichment_upstream_network
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_enrichment_upstream_network <- function(hc,
                                                enrichment_mode = "selected",
                                                upstream_mode = "selected",
                                                clusters = c("all"),
                                                max_enrichment_per_module = NULL,
                                                max_upstream_per_module = NULL,
                                                label_mode = "both",
                                                show_plot = TRUE,
                                                save_pdf = TRUE,
                                                pdf_name = "Module_Knowledge_Network.pdf",
                                                overall_plot_scale = 1) {
  .hc_run_legacy(
    hc = hc,
    fun = "plot_enrichment_upstream_network",
    enrichment_mode = enrichment_mode,
    upstream_mode = upstream_mode,
    clusters = clusters,
    max_enrichment_per_module = max_enrichment_per_module,
    max_upstream_per_module = max_upstream_per_module,
    label_mode = label_mode,
    show_plot = show_plot,
    save_pdf = save_pdf,
    pdf_name = pdf_name,
    overall_plot_scale = overall_plot_scale
  )
}

#' Run `check_dirs()` with S4 state synchronization
#'
#' @param hc A `HCoCenaExperiment`.
#' @param create_output_dir Boolean. If TRUE and `dir_output` is missing, create it.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_check_dirs <- function(hc, create_output_dir = TRUE) {
  .hc_run_legacy(hc = hc, fun = "check_dirs", create_output_dir = create_output_dir)
}

#' Run `init_save_folder()` with S4 state synchronization
#'
#' @param hc A `HCoCenaExperiment`.
#' @param name Folder name. Use `""` to write directly into `dir_output`.
#' @param use_output_dir Boolean. If TRUE, ignore `name` and use `dir_output` directly.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_init_save_folder <- function(hc, name, use_output_dir = FALSE) {
  if (isTRUE(use_output_dir)) {
    name <- ""
  }
  .hc_run_legacy(hc = hc, fun = "init_save_folder", name = name)
}

#' Plot cut-off diagnostics (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `plot_cutoffs()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_cutoffs <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "plot_cutoffs", ...)
}

#' Plot degree distributions (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_deg_dist <- function(hc) {
  .hc_run_legacy(hc = hc, fun = "plot_deg_dist")
}

#' Plot module heatmap (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `plot_cluster_heatmap()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_cluster_heatmap <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "plot_cluster_heatmap", ...)
}

#' Plot integrated network (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `plot_integrated_network()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_integrated_network <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "plot_integrated_network", ...)
}

#' Plot network colored by GFC (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_gfc_network <- function(hc) {
  .hc_run_legacy(hc = hc, fun = "plot_GFC_network")
}

#' TF enrichment per module (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `TF_overrep_module()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_tf_overrep_module <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "TF_overrep_module", ...)
}

#' TF enrichment network-wide (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `TF_overrep_network()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_tf_overrep_network <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "TF_overrep_network", ...)
}

#' Check TF targets (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param TF Transcription factor symbol.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_check_tf <- function(hc, TF) {
  .hc_run_legacy(hc = hc, fun = "check_tf", TF = TF)
}

#' Write session info (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_write_session_info <- function(hc) {
  .hc_run_legacy(hc = hc, fun = "write_session_info")
}

#' Suggest top variable genes (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_suggest_topvar <- function(hc) {
  .hc_run_legacy(hc = hc, fun = "suggest_topvar")
}

#' Plot sample distributions (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `plot_sample_distributions()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_sample_distributions <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "plot_sample_distributions", ...)
}

#' PCA plotting (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `PCA()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_pca <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "PCA", ...)
}

#' Meta-data plotting (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `meta_plot()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_meta_plot <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "meta_plot", ...)
}

#' Export clusters (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_export_clusters <- function(hc) {
  .hc_run_legacy(hc = hc, fun = "export_clusters")
}

#' Module scores (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `get_module_scores()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_get_module_scores <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "get_module_scores", ...)
}

#' Alluvial comparison plots (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_algo_alluvial <- function(hc) {
  .hc_run_legacy(hc = hc, fun = "algo_alluvial")
}

#' PCA algorithm comparison (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `PCA_algo_compare()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_pca_algo_compare <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "PCA_algo_compare", ...)
}

#' Update clustering algorithm (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `update_clustering_algorithm()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_update_clustering_algorithm <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "update_clustering_algorithm", ...)
}

#' Export network to local folder (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `export_to_local_folder()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_export_to_local_folder <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "export_to_local_folder", ...)
}

#' Import Cytoscape layout from local folder (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `import_layout_from_local_folder()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_import_layout_from_local_folder <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "import_layout_from_local_folder", ...)
}

#' Export network to Cytoscape (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `export_to_cytoscape()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_export_to_cytoscape <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "export_to_cytoscape", ...)
}

#' Import layout from Cytoscape (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_import_layout_from_cytoscape <- function(hc) {
  .hc_run_legacy(hc = hc, fun = "import_layout_from_cytoscape")
}

#' Hub detection (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `find_hubs()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_find_hubs <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "find_hubs", ...)
}

#' Visualize gene expression (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `visualize_gene_expression()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_visualize_gene_expression <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "visualize_gene_expression", ...)
}

#' Highlight gene set in network (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `highlight_geneset()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_highlight_geneset <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "highlight_geneset", ...)
}

#' Highlight single cluster in network (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `colour_single_cluster()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_colour_single_cluster <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "colour_single_cluster", ...)
}

#' Add categorical module heatmap annotations (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `col_anno_categorical()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_col_anno_categorical <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "col_anno_categorical", ...)
}

#' Correlate categorical metadata with modules (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `meta_correlation_cat()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_meta_correlation_cat <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "meta_correlation_cat", ...)
}
