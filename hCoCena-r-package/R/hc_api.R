#' Initialize an empty `HCoCenaExperiment`
#'
#' @return A valid, empty `HCoCenaExperiment`.
#' @export
hc_init <- function() {
  methods::new("HCoCenaExperiment")
}

.hc_normalize_path_value <- function(x, arg_name) {
  if (base::isFALSE(x) || base::identical(x, FALSE)) {
    return(FALSE)
  }

  if (base::length(x) != 1) {
    stop("`", arg_name, "` must be a scalar path string or FALSE.")
  }

  x_chr <- base::as.character(x[[1]])
  if (base::is.na(x_chr) || !base::nzchar(x_chr)) {
    stop("`", arg_name, "` must be a non-empty path string or FALSE.")
  }

  x_chr
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

  dir_count_data <- .hc_normalize_path_value(dir_count_data, "dir_count_data")
  dir_annotation <- .hc_normalize_path_value(dir_annotation, "dir_annotation")
  dir_reference_files <- .hc_normalize_path_value(dir_reference_files, "dir_reference_files")
  dir_output <- .hc_normalize_path_value(dir_output, "dir_output")

  hc@config@paths <- S4Vectors::DataFrame(
    dir_count_data = dir_count_data,
    dir_annotation = dir_annotation,
    dir_reference_files = dir_reference_files,
    dir_output = dir_output
  )
  methods::validObject(hc)
  hc
}

.hc_path_with_trailing_slash <- function(x) {
  x <- base::as.character(x[[1]])
  x <- base::gsub("\\\\", "/", x)
  if (base::grepl("/$", x)) x else base::paste0(x, "/")
}

.hc_scalar_or_null <- function(x, arg_name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (base::length(x) != 1) {
    stop("`", arg_name, "` must be NULL or a scalar string.")
  }
  out <- base::as.character(x[[1]])
  if (base::is.na(out) || !base::nzchar(out)) {
    return(NULL)
  }
  out
}

.hc_find_local_dir <- function(project_root, subdir) {
  if (is.null(subdir) || !base::nzchar(subdir)) {
    return(NULL)
  }
  candidate <- base::file.path(project_root, subdir)
  if (!base::dir.exists(candidate)) {
    return(NULL)
  }
  .hc_path_with_trailing_slash(base::normalizePath(candidate, winslash = "/", mustWork = TRUE))
}

#' Resolve hCoCena input/output directories for Docker or local projects
#'
#' @param project_root_hint Optional local project root to prioritize if it exists.
#'   If missing or not found, the current working directory is used.
#' @param local_count_subdir Local count-data subfolder below `project_root_hint`.
#' @param local_annotation_subdir Local annotation-data subfolder below `project_root_hint`.
#' @param local_reference_subdir Local reference-files subfolder below `project_root_hint`.
#' @param local_output_subdir Local output subfolder below `project_root_hint`.
#' @param docker_reference_dir Docker reference-files directory used to detect Docker mode.
#' @param docker_count_dir Docker count-data directory.
#' @param docker_annotation_dir Docker annotation-data directory.
#' @param docker_output_dir Docker output directory.
#' @param create_docker_dirs Boolean. If TRUE and Docker mode is detected, create the
#'   Docker directories when missing.
#' @param fallback_count_data Fallback count-data directory string used outside Docker
#'   when no local folder was detected.
#' @param fallback_annotation Fallback annotation directory string used outside Docker
#'   when no local folder was detected.
#' @param fallback_reference_files Fallback reference-files directory string used outside
#'   Docker when no local folder was detected.
#' @param fallback_output Fallback output directory string used outside Docker when no
#'   local folder was detected.
#' @return Named list with `docker_mode`, `project_root`, `dir_count_data`,
#'   `dir_annotation`, `dir_reference_files`, and `dir_output`.
#' @export
hc_resolve_paths <- function(
    project_root_hint = NULL,
    local_count_subdir = "count_data",
    local_annotation_subdir = "annotation_data",
    local_reference_subdir = "reference_files",
    local_output_subdir = "output",
    docker_reference_dir = "/home/rstudio/reference_files",
    docker_count_dir = "/home/rstudio/count_data",
    docker_annotation_dir = "/home/rstudio/annotation_data",
    docker_output_dir = "/home/rstudio/output",
    create_docker_dirs = TRUE,
    fallback_count_data = "PATH_TO_FOLDER_THAT_HOLDS_YOUR_GENE_EXPRESSION_TABLES/",
    fallback_annotation = "PATH_TO_FOLDER_THAT_HOLDS_YOUR_ANNOTATION_TABLES/",
    fallback_reference_files = "PATH_TO_REFERENCE_FILES_FOLDER/",
    fallback_output = "PATH_TO_FOLDER_WHERE_HCOCENA_SHOULD_SAVE_ALL_ANALYSIS_OUTPUTS/"
) {
  project_root_hint <- .hc_scalar_or_null(project_root_hint, "project_root_hint")
  local_count_subdir <- .hc_scalar_or_null(local_count_subdir, "local_count_subdir")
  local_annotation_subdir <- .hc_scalar_or_null(local_annotation_subdir, "local_annotation_subdir")
  local_reference_subdir <- .hc_scalar_or_null(local_reference_subdir, "local_reference_subdir")
  local_output_subdir <- .hc_scalar_or_null(local_output_subdir, "local_output_subdir")
  docker_reference_dir <- .hc_scalar_or_null(docker_reference_dir, "docker_reference_dir")
  docker_count_dir <- .hc_scalar_or_null(docker_count_dir, "docker_count_dir")
  docker_annotation_dir <- .hc_scalar_or_null(docker_annotation_dir, "docker_annotation_dir")
  docker_output_dir <- .hc_scalar_or_null(docker_output_dir, "docker_output_dir")
  fallback_count_data <- .hc_scalar_or_null(fallback_count_data, "fallback_count_data")
  fallback_annotation <- .hc_scalar_or_null(fallback_annotation, "fallback_annotation")
  fallback_reference_files <- .hc_scalar_or_null(fallback_reference_files, "fallback_reference_files")
  fallback_output <- .hc_scalar_or_null(fallback_output, "fallback_output")

  docker_mode <- !is.null(docker_reference_dir) && base::dir.exists(docker_reference_dir)
  if (isTRUE(docker_mode) && isTRUE(create_docker_dirs)) {
    for (x in base::c(docker_count_dir, docker_annotation_dir, docker_reference_dir, docker_output_dir)) {
      if (!is.null(x) && base::nzchar(x)) {
        base::dir.create(x, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }

  project_root <- if (!is.null(project_root_hint) && base::dir.exists(project_root_hint)) {
    base::normalizePath(project_root_hint, winslash = "/", mustWork = TRUE)
  } else {
    base::normalizePath(".", winslash = "/", mustWork = TRUE)
  }

  local_count_dir <- .hc_find_local_dir(project_root, local_count_subdir)
  local_annotation_dir <- .hc_find_local_dir(project_root, local_annotation_subdir)
  local_reference_dir <- .hc_find_local_dir(project_root, local_reference_subdir)
  local_output_dir <- .hc_find_local_dir(project_root, local_output_subdir)

  dir_count_data <- if (!is.null(local_count_dir)) {
    local_count_dir
  } else if (isTRUE(docker_mode) && !is.null(docker_count_dir)) {
    .hc_path_with_trailing_slash(docker_count_dir)
  } else {
    .hc_path_with_trailing_slash(fallback_count_data)
  }

  dir_annotation <- if (!is.null(local_annotation_dir)) {
    local_annotation_dir
  } else if (isTRUE(docker_mode) && !is.null(docker_annotation_dir)) {
    .hc_path_with_trailing_slash(docker_annotation_dir)
  } else {
    .hc_path_with_trailing_slash(fallback_annotation)
  }

  dir_reference_files <- if (!is.null(local_reference_dir)) {
    local_reference_dir
  } else if (isTRUE(docker_mode) && !is.null(docker_reference_dir)) {
    .hc_path_with_trailing_slash(docker_reference_dir)
  } else {
    .hc_path_with_trailing_slash(fallback_reference_files)
  }

  dir_output <- if (!is.null(local_output_dir)) {
    local_output_dir
  } else if (isTRUE(docker_mode) && !is.null(docker_output_dir)) {
    .hc_path_with_trailing_slash(docker_output_dir)
  } else {
    .hc_path_with_trailing_slash(fallback_output)
  }

  list(
    docker_mode = isTRUE(docker_mode),
    project_root = project_root,
    dir_count_data = dir_count_data,
    dir_annotation = dir_annotation,
    dir_reference_files = dir_reference_files,
    dir_output = dir_output
  )
}

#' Resolve and apply standard paths in one call
#'
#' @param hc A `HCoCenaExperiment`.
#' @param dir_count_data Optional explicit override for count-data directory.
#' @param dir_annotation Optional explicit override for annotation directory.
#' @param dir_reference_files Optional explicit override for reference-files directory.
#' @param dir_output Optional explicit override for output directory.
#' @inheritParams hc_resolve_paths
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_auto_set_paths <- function(
    hc,
    project_root_hint = NULL,
    local_count_subdir = "count_data",
    local_annotation_subdir = "annotation_data",
    local_reference_subdir = "reference_files",
    local_output_subdir = "output",
    docker_reference_dir = "/home/rstudio/reference_files",
    docker_count_dir = "/home/rstudio/count_data",
    docker_annotation_dir = "/home/rstudio/annotation_data",
    docker_output_dir = "/home/rstudio/output",
    create_docker_dirs = TRUE,
    fallback_count_data = "PATH_TO_FOLDER_THAT_HOLDS_YOUR_GENE_EXPRESSION_TABLES/",
    fallback_annotation = "PATH_TO_FOLDER_THAT_HOLDS_YOUR_ANNOTATION_TABLES/",
    fallback_reference_files = "PATH_TO_REFERENCE_FILES_FOLDER/",
    fallback_output = "PATH_TO_FOLDER_WHERE_HCOCENA_SHOULD_SAVE_ALL_ANALYSIS_OUTPUTS/",
    dir_count_data = NULL,
    dir_annotation = NULL,
    dir_reference_files = NULL,
    dir_output = NULL
) {
  resolved <- hc_resolve_paths(
    project_root_hint = project_root_hint,
    local_count_subdir = local_count_subdir,
    local_annotation_subdir = local_annotation_subdir,
    local_reference_subdir = local_reference_subdir,
    local_output_subdir = local_output_subdir,
    docker_reference_dir = docker_reference_dir,
    docker_count_dir = docker_count_dir,
    docker_annotation_dir = docker_annotation_dir,
    docker_output_dir = docker_output_dir,
    create_docker_dirs = create_docker_dirs,
    fallback_count_data = fallback_count_data,
    fallback_annotation = fallback_annotation,
    fallback_reference_files = fallback_reference_files,
    fallback_output = fallback_output
  )

  final_count <- if (is.null(dir_count_data)) {
    resolved$dir_count_data
  } else {
    .hc_normalize_path_value(dir_count_data, "dir_count_data")
  }
  final_anno <- if (is.null(dir_annotation)) {
    resolved$dir_annotation
  } else {
    .hc_normalize_path_value(dir_annotation, "dir_annotation")
  }
  final_ref <- if (is.null(dir_reference_files)) {
    resolved$dir_reference_files
  } else {
    .hc_normalize_path_value(dir_reference_files, "dir_reference_files")
  }
  final_out <- if (is.null(dir_output)) {
    resolved$dir_output
  } else {
    .hc_normalize_path_value(dir_output, "dir_output")
  }

  hc_set_paths(
    hc,
    dir_count_data = final_count,
    dir_annotation = final_anno,
    dir_reference_files = final_ref,
    dir_output = final_out
  )
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
    if (length(pair) < 2) {
      stop(
        "Layer `", layer_names[[i]], "` must provide exactly two entries: ",
        "count source and annotation source."
      )
    }

    count_src <- pair[[1]]
    anno_src <- pair[[2]]

    if (!is.character(count_src) || length(count_src) != 1 || is.na(count_src) || !nzchar(count_src)) {
      stop(
        "Layer `", layer_names[[i]], "`: count source must be a single string.\n",
        "If you use objects from the environment, pass quoted object names, e.g. ",
        "`c(\"counts_df\", \"anno_df\")` (not `c(counts_df, anno_df)`)."
      )
    }
    if (!is.character(anno_src) || length(anno_src) != 1 || is.na(anno_src) || !nzchar(anno_src)) {
      stop(
        "Layer `", layer_names[[i]], "`: annotation source must be a single string.\n",
        "If you use objects from the environment, pass quoted object names, e.g. ",
        "`c(\"counts_df\", \"anno_df\")` (not `c(counts_df, anno_df)`)."
      )
    }

    rows[[i]] <- list(
      layer_id = paste0("set", i),
      layer_name = as.character(layer_names[[i]]),
      count_source = as.character(count_src),
      annotation_source = as.character(anno_src)
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
                                   organism = "human",
                                   control_keyword = "none",
                                   variable_of_interest = "merged",
                                   min_nodes_number_for_network = 50,
                                   min_nodes_number_for_cluster = 50,
                                   range_GFC = 2.0,
                                   layout_algorithm = "layout_with_stress",
                                   data_in_log = TRUE) {
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
#' @param ... Additional supplementary files passed through to `set_supp_files()`.
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

# Internal helpers for S4 functional-enrichment panel redraw.
# @noRd
.hc_functional_enrichment_entries <- function(hc) {
  sat <- as.list(hc@satellite)
  enrich <- sat[["enrichments"]]
  if (is.null(enrich)) {
    sat_keys <- names(sat)
    if (!is.null(sat_keys) && any(grepl("^top_", sat_keys))) {
      enrich <- sat
    }
  }
  if (is.null(enrich)) {
    return(list())
  }
  enrich_names <- names(enrich)
  enrich <- as.list(enrich)
  if (is.null(names(enrich)) && !is.null(enrich_names) && length(enrich_names) == length(enrich)) {
    names(enrich) <- enrich_names
  }
  if (length(enrich) == 0) {
    return(list())
  }
  enrich
}

.hc_functional_enrichment_panel_keys <- function(enrich) {
  if (is.null(enrich) || length(enrich) == 0) {
    return(character(0))
  }

  panel_keys <- names(enrich)
  panel_keys <- panel_keys[grepl("^top_", panel_keys)]
  if (length(panel_keys) == 0) {
    return(character(0))
  }
  panel_keys <- unique(panel_keys)

  requested_order <- enrich[["panel_order"]]
  if (!is.null(requested_order)) {
    requested_order <- base::as.character(requested_order)
    requested_order <- requested_order[!base::is.na(requested_order) & base::nzchar(requested_order)]
    requested_order <- requested_order[requested_order %in% panel_keys]
    panel_keys <- unique(base::c(requested_order, panel_keys))
  }

  single_db_keys <- panel_keys[!panel_keys %in% c("top_all_dbs", "top_all_dbs_mixed")]
  combined_keys <- intersect(c("top_all_dbs", "top_all_dbs_mixed"), panel_keys)
  panel_keys <- unique(c(single_db_keys, combined_keys))

  panel_keys
}

.hc_functional_enrichment_has_draw_object <- function(entry) {
  entry <- tryCatch(as.list(entry), error = function(e) NULL)
  if (is.null(entry)) {
    return(FALSE)
  }
  for (nm in c("p", "hc_heatmap", "enrichment_plot")) {
    if (nm %in% names(entry) && !is.null(entry[[nm]])) {
      return(TRUE)
    }
  }
  FALSE
}

.hc_functional_enrichment_draw_object <- function(entry, heatmap_side = "left") {
  if (is.null(entry)) {
    return(NULL)
  }
  entry <- tryCatch(as.list(entry), error = function(e) NULL)
  if (is.null(entry)) {
    return(NULL)
  }

  clone_obj <- function(x) {
    .hc_safe_deep_clone(x, context = "stored enrichment panel object")
  }

  draw_obj <- NULL
  has_components <- (
    "hc_heatmap" %in% names(entry) &&
      "enrichment_plot" %in% names(entry)
  )
  if (isTRUE(has_components)) {
    hc_ht <- clone_obj(entry[["hc_heatmap"]])
    enr_ht <- clone_obj(entry[["enrichment_plot"]])
    if (inherits(hc_ht, "Heatmap") && inherits(enr_ht, "Heatmap")) {
      draw_obj <- if (identical(heatmap_side, "right")) {
        enr_ht + hc_ht
      } else {
        hc_ht + enr_ht
      }
    }
  }
  if (is.null(draw_obj)) {
    p <- clone_obj(entry[["p"]])
    if (inherits(p, "HeatmapList")) {
      draw_obj <- p
    }
  }
  draw_obj
}

.hc_functional_enrichment_panel_title <- function(panel_key) {
  if (is.null(panel_key) || length(panel_key) != 1 || is.na(panel_key)) {
    return(NULL)
  }
  panel_key <- as.character(panel_key[[1]])
  if (identical(panel_key, "top_all_dbs")) {
    return("Combined enrichment (all DBs)")
  }
  if (identical(panel_key, "top_all_dbs_mixed")) {
    return("Combined enrichment (all DBs, mixed)")
  }
  db_name <- sub("^top_", "", panel_key)
  if (!nzchar(db_name) || identical(db_name, panel_key)) {
    return(NULL)
  }
  db_name <- gsub("_", " ", db_name, fixed = TRUE)
  paste0(db_name, " enrichment")
}

.hc_functional_enrichment_panel_target <- function(panel_key) {
  if (is.null(panel_key) || length(panel_key) != 1 || is.na(panel_key)) {
    return(NULL)
  }
  panel_key <- as.character(panel_key[[1]])
  if (identical(panel_key, "top_all_dbs")) {
    return("enrichment_all")
  }
  if (identical(panel_key, "top_all_dbs_mixed")) {
    return("enrichment_all_mixed")
  }
  "enrichment"
}

.hc_draw_functional_enrichment_panel_title <- function(panel_key, fontsize = 16) {
  panel_title <- .hc_functional_enrichment_panel_title(panel_key)
  panel_target <- .hc_functional_enrichment_panel_target(panel_key)
  if (is.null(panel_title) || !nzchar(panel_title) || is.null(panel_target)) {
    return(invisible(NULL))
  }

  panel_key_chr <- as.character(panel_key[[1]])
  title_drawn <- FALSE
  slice_candidates <- if (identical(panel_key_chr, "top_all_dbs")) c(2L, 1L) else 1L
  title_offset_mm <- max(4, 0.35 * as.numeric(fontsize))
  components <- try(ComplexHeatmap::list_components(), silent = TRUE)
  if (inherits(components, "try-error")) {
    components <- character(0)
  }

  if (identical(panel_key_chr, "top_all_dbs") && length(components) > 0) {
    body_pattern <- paste0("^", panel_target, "_heatmap_body_[0-9]+_[0-9]+$")
    body_components <- components[grepl(body_pattern, components)]
    if (length(body_components) >= 1) {
      draw_res <- try(
        {
          left_edges <- numeric(0)
          right_edges <- numeric(0)
          top_edges <- numeric(0)

          for (vp_name in body_components) {
            grid::seekViewport(vp_name)
            loc_left <- grid::deviceLoc(
              x = grid::unit(0, "npc"),
              y = grid::unit(1, "npc")
            )
            loc_right <- grid::deviceLoc(
              x = grid::unit(1, "npc"),
              y = grid::unit(1, "npc")
            )
            grid::upViewport(0)
            left_edges <- c(left_edges, grid::convertX(loc_left[["x"]], "inches", valueOnly = TRUE))
            right_edges <- c(right_edges, grid::convertX(loc_right[["x"]], "inches", valueOnly = TRUE))
            top_edges <- c(top_edges, grid::convertY(loc_left[["y"]], "inches", valueOnly = TRUE))
          }

          loc_x <- grid::unit((min(left_edges) + max(right_edges)) / 2, "inches")
          loc_y <- grid::unit(max(top_edges), "inches")
          grid::grid.text(
            label = panel_title,
            x = loc_x,
            y = loc_y + grid::unit(title_offset_mm, "mm"),
            just = c("center", "bottom"),
            gp = grid::gpar(fontsize = fontsize, fontface = "bold")
          )
        },
        silent = TRUE
      )
      try(grid::upViewport(0), silent = TRUE)
      if (!inherits(draw_res, "try-error")) {
        return(invisible(NULL))
      }
    }
  }

  for (slice_idx in slice_candidates) {
    body_vp_name <- paste0(panel_target, "_heatmap_body_", slice_idx, "_1")
    if (length(components) > 0 && !(body_vp_name %in% components)) {
      next
    }
    draw_res <- try(
      {
        grid::seekViewport(body_vp_name)
        loc <- grid::deviceLoc(
          x = grid::unit(0.5, "npc"),
          y = grid::unit(1, "npc")
        )
        grid::upViewport(0)
        loc_x <- loc[["x"]]
        loc_y <- loc[["y"]]
        grid::grid.text(
          label = panel_title,
          x = loc_x,
          y = loc_y + grid::unit(title_offset_mm, "mm"),
          just = c("center", "bottom"),
          gp = grid::gpar(fontsize = fontsize, fontface = "bold")
        )
      },
      silent = TRUE
    )
    try(grid::upViewport(0), silent = TRUE)
    if (!inherits(draw_res, "try-error")) {
      title_drawn <- TRUE
      break
    }
  }
  if (!title_drawn) {
    return(invisible(NULL))
  }
  invisible(NULL)
}

.hc_draw_functional_enrichment_panels <- function(hc, heatmap_side = "left", record_history = FALSE) {
  enrich <- .hc_functional_enrichment_entries(hc)
  panel_keys <- .hc_functional_enrichment_panel_keys(enrich)
  if (length(panel_keys) == 0) {
    return(invisible(NULL))
  }

  if (isTRUE(record_history) && grDevices::dev.cur() > 1) {
    try(grDevices::dev.control(displaylist = "enable"), silent = TRUE)
  }

  draw_errors <- character(0)
  for (k in panel_keys) {
    draw_obj <- .hc_functional_enrichment_draw_object(enrich[[k]], heatmap_side = heatmap_side)
    if (is.null(draw_obj)) {
      next
    }

    draw_res <- try(
      ComplexHeatmap::draw(
        draw_obj,
        newpage = TRUE,
        merge_legends = TRUE,
        show_annotation_legend = TRUE,
        show_heatmap_legend = TRUE
      ),
      silent = TRUE
    )
    if (inherits(draw_res, "try-error")) {
      draw_errors <- c(draw_errors, k)
      next
    }
    try(.hc_draw_functional_enrichment_panel_title(k), silent = TRUE)
    if (isTRUE(record_history) && grDevices::dev.cur() > 1) {
      try(grDevices::recordPlot(), silent = TRUE)
    }
  }

  if (length(draw_errors) > 0) {
    warning(
      "Could not redraw enrichment panel(s): ",
      paste(unique(draw_errors), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.hc_emit_functional_enrichment_knitr_panels <- function(hc,
                                                        heatmap_side = "left",
                                                        panel_keys = NULL,
                                                        res = 200) {
  if (!isTRUE(getOption("knitr.in.progress")) &&
      !isTRUE(getOption("rstudio.notebook.executing"))) {
    return(invisible(NULL))
  }
  if (requireNamespace("knitr", quietly = TRUE)) {
    old_fig_keep <- try(knitr::opts_current$get("fig.keep"), silent = TRUE)
    if (!inherits(old_fig_keep, "try-error")) {
      on.exit(
        try(knitr::opts_current$set(fig.keep = old_fig_keep), silent = TRUE),
        add = TRUE
      )
    }
    try(knitr::opts_current$set(fig.keep = "all"), silent = TRUE)
  }

  enrich <- .hc_functional_enrichment_entries(hc)
  available_panel_keys <- .hc_functional_enrichment_panel_keys(enrich)
  combined_panel_keys <- intersect(c("top_all_dbs", "top_all_dbs_mixed"), available_panel_keys)
  if (is.null(panel_keys)) {
    panel_keys <- available_panel_keys
  } else {
    panel_keys <- base::as.character(panel_keys)
    panel_keys <- panel_keys[!base::is.na(panel_keys) & base::nzchar(panel_keys)]
    panel_keys <- panel_keys[panel_keys %in% available_panel_keys]
  }
  panel_keys <- unique(panel_keys)
  panel_keys <- c(
    panel_keys[!panel_keys %in% combined_panel_keys],
    intersect(combined_panel_keys, panel_keys),
    setdiff(panel_keys[panel_keys %in% combined_panel_keys], combined_panel_keys)
  )
  if (length(panel_keys) == 0) {
    return(invisible(NULL))
  }
  panel_keys <- c(
    panel_keys[!panel_keys %in% c("top_all_dbs", "top_all_dbs_mixed")],
    intersect(c("top_all_dbs", "top_all_dbs_mixed"), panel_keys)
  )

  tmp_plot_dir <- tempfile(pattern = "hc_functional_enrichment_panels_")
  dir.create(tmp_plot_dir, recursive = TRUE, showWarnings = FALSE)
  png_paths <- character(0)
  failed_panels <- character(0)

  for (idx in base::seq_along(panel_keys)) {
    k <- panel_keys[[idx]]
    draw_obj <- .hc_functional_enrichment_draw_object(enrich[[k]], heatmap_side = heatmap_side)
    if (is.null(draw_obj)) {
      next
    }

    is_combined_panel <- k %in% c("top_all_dbs", "top_all_dbs_mixed")
    panel_width_in <- if (is_combined_panel) 18 else 14
    panel_height_in <- if (is_combined_panel) 10 else 9
    safe_key <- gsub("[^A-Za-z0-9._-]+", "_", k)
    png_file <- base::file.path(tmp_plot_dir, sprintf("%03d_%s.png", idx, safe_key))
    opened <- FALSE
    open_res <- try(
      {
        grDevices::png(
          filename = png_file,
          width = panel_width_in,
          height = panel_height_in,
          units = "in",
          res = res
        )
        opened <- TRUE
      },
      silent = TRUE
    )
    if (!isTRUE(opened) || inherits(open_res, "try-error")) {
      failed_panels <- c(failed_panels, k)
      next
    }

    draw_ok <- FALSE
    draw_res <- try(
      {
        ComplexHeatmap::draw(
          draw_obj,
          newpage = TRUE,
          merge_legends = TRUE,
          show_annotation_legend = TRUE,
          show_heatmap_legend = TRUE
        )
        .hc_draw_functional_enrichment_panel_title(k)
        draw_ok <- TRUE
      },
      silent = TRUE
    )
    try(grDevices::dev.off(), silent = TRUE)

    if (inherits(draw_res, "try-error") || !isTRUE(draw_ok) || !file.exists(png_file)) {
      failed_panels <- c(failed_panels, k)
      if (file.exists(png_file)) {
        try(unlink(png_file, force = TRUE), silent = TRUE)
      }
      next
    }
    png_paths <- c(png_paths, png_file)
  }

  if (length(png_paths) > 0) {
    if (requireNamespace("knitr", quietly = TRUE)) {
      print(knitr::include_graphics(png_paths))
    } else {
      warning(
        "Package `knitr` is required for notebook panel rendering.",
        call. = FALSE
      )
    }
  }

  if (length(failed_panels) > 0) {
    warning(
      "Could not render enrichment panel(s) for knitr output: ",
      paste(unique(failed_panels), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' Draw saved functional-enrichment panels from the current object
#'
#' This is a convenience plotting helper for notebook/interactive usage. It
#' draws the enrichment panels already stored in `hc@satellite$enrichments`
#' after `hc_functional_enrichment()`, without recomputing enrichment.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param panels Optional character vector of panel names to draw.
#'   Default (`NULL`) draws all available panels in standard order.
#'   Accepted values include `"top_Hallmark"`, `"top_Kegg"`,
#'   `"top_all_dbs"`, `"top_all_dbs_mixed"` (or without `"top_"` prefix).
#' @param heatmap_side One of `"left"` (default) or `"right"`.
#' @return Invisibly returns `hc`.
#' @export
hc_plot_enrichment_panels <- function(hc,
                                      panels = NULL,
                                      heatmap_side = "left") {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  heatmap_side <- base::match.arg(heatmap_side, choices = c("left", "right"))

  enrich <- .hc_functional_enrichment_entries(hc)
  if (length(enrich) == 0) {
    stop(
      "No enrichment panels found in `hc@satellite$enrichments`. ",
      "Run `hc_functional_enrichment()` first."
    )
  }

  available_panel_keys <- .hc_functional_enrichment_panel_keys(enrich)
  if (length(available_panel_keys) == 0) {
    stop(
      "No drawable enrichment panels found. ",
      "Expected entries like `top_Hallmark`, `top_all_dbs`, ... in ",
      "`hc@satellite$enrichments`."
    )
  }

  combined_panel_keys <- intersect(c("top_all_dbs", "top_all_dbs_mixed"), available_panel_keys)
  single_panel_keys <- available_panel_keys[!available_panel_keys %in% combined_panel_keys]
  panel_keys <- c(single_panel_keys, combined_panel_keys)
  if (!is.null(panels)) {
    panels <- base::as.character(panels)
    panels <- panels[!base::is.na(panels) & base::nzchar(panels)]
    panels <- ifelse(grepl("^top_", panels), panels, paste0("top_", panels))
    missing_panels <- panels[!panels %in% available_panel_keys]
    if (length(missing_panels) > 0) {
      warning(
        "Skipping unknown enrichment panel(s): ",
        paste(unique(missing_panels), collapse = ", "),
        call. = FALSE
      )
    }
    requested <- panels[panels %in% available_panel_keys]
    if (length(requested) == 0) {
      stop(
        "None of the requested panels are available. Available panels: ",
        paste(available_panel_keys, collapse = ", ")
      )
    }

    panel_keys <- unique(requested)
  }

  panel_keys <- unique(panel_keys)
  drawable_mask <- base::vapply(
    panel_keys,
    function(k) .hc_functional_enrichment_has_draw_object(enrich[[k]]),
    FUN.VALUE = base::logical(1)
  )
  if (!base::any(drawable_mask)) {
    stop(
      "Selected enrichment panels were not stored as drawable objects in `hc`. ",
      "This usually happens in memory-saving mode for multi-database enrichment runs. ",
      "Rerun `hc_functional_enrichment(..., store_panel_objects = 'always')` if you need redraw from the object."
    )
  }
  message(
    "Drawing enrichment panels in order: ",
    paste(panel_keys, collapse = ", ")
  )

  in_notebook <- isTRUE(getOption("knitr.in.progress")) ||
    isTRUE(getOption("rstudio.notebook.executing"))
  if (isTRUE(in_notebook)) {
    for (k in panel_keys) {
      .hc_emit_functional_enrichment_knitr_panels(
        hc = hc,
        heatmap_side = heatmap_side,
        panel_keys = k
      )
    }
    return(invisible(hc))
  }

  draw_errors <- character(0)
  for (k in panel_keys) {
    draw_obj <- .hc_functional_enrichment_draw_object(enrich[[k]], heatmap_side = heatmap_side)
    if (is.null(draw_obj)) {
      draw_errors <- c(draw_errors, k)
      next
    }

    draw_res <- try({
      grid::grid.newpage()
      ComplexHeatmap::draw(
        draw_obj,
        newpage = FALSE,
        merge_legends = TRUE,
        show_annotation_legend = TRUE,
        show_heatmap_legend = TRUE
      )
    }, silent = TRUE)
    if (inherits(draw_res, "try-error")) {
      draw_errors <- c(draw_errors, k)
      next
    }
    try(.hc_draw_functional_enrichment_panel_title(k), silent = TRUE)
  }

  if (length(draw_errors) > 0) {
    warning(
      "Could not draw enrichment panel(s): ",
      paste(unique(draw_errors), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(hc)
}

#' Functional enrichment (S4 API)
#'
#' @inheritParams functional_enrichment
#' @param hc A `HCoCenaExperiment`.
#' @param ... Additional plotting arguments forwarded to the legacy
#'   `functional_enrichment()` implementation.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_functional_enrichment <- function(hc,
                                     gene_sets = c("Go", "Kegg", "Hallmark", "Reactome"),
                                     custom_gmt_files = NULL,
                                     top = 5,
                                     clusters = c("all"),
                                     padj = "BH",
                                     qval = 0.05,
                                     heatmap_side = "left",
                                     heatmap_cluster_rows = FALSE,
                                     heatmap_cluster_columns = FALSE,
                                     heatmap_show_row_dend = FALSE,
                                     heatmap_show_column_dend = FALSE,
                                     heatmap_col_order = NULL,
                                     heatmap_order = NULL,
                                     heatmap_module_label_mode = "same",
                                     heatmap_show_gene_counts = FALSE,
                                     heatmap_column_label_fontsize = NULL,
                                     heatmap_module_label_fontsize = NULL,
                                     legend_fontsize = NULL,
                                     enrichment_label_fontsize = NULL,
                                     enrichment_db_header_fontsize = NULL,
                                     enrichment_label_wrap = FALSE,
                                     enrichment_label_wrap_width = 30,
                                     gfc_scale_limits = NULL,
                                     store_panel_objects = c("auto", "always", "never"),
                                     pdf_width = NULL,
                                     pdf_height = NULL,
                                     pdf_pointsize = 11,
                                     ...) {
  .hc_run_legacy(
    hc = hc,
    fun = "functional_enrichment",
    gene_sets = gene_sets,
    custom_gmt_files = custom_gmt_files,
    top = top,
    clusters = clusters,
    padj = padj,
    qval = qval,
    heatmap_side = heatmap_side,
    heatmap_cluster_rows = heatmap_cluster_rows,
    heatmap_cluster_columns = heatmap_cluster_columns,
    heatmap_show_row_dend = heatmap_show_row_dend,
    heatmap_show_column_dend = heatmap_show_column_dend,
    heatmap_col_order = heatmap_col_order,
    heatmap_order = heatmap_order,
    heatmap_module_label_mode = heatmap_module_label_mode,
    heatmap_show_gene_counts = heatmap_show_gene_counts,
    heatmap_column_label_fontsize = heatmap_column_label_fontsize,
    heatmap_module_label_fontsize = heatmap_module_label_fontsize,
    legend_fontsize = legend_fontsize,
    enrichment_label_fontsize = enrichment_label_fontsize,
    enrichment_db_header_fontsize = enrichment_db_header_fontsize,
    enrichment_label_wrap = enrichment_label_wrap,
    enrichment_label_wrap_width = enrichment_label_wrap_width,
    gfc_scale_limits = gfc_scale_limits,
    store_panel_objects = store_panel_objects,
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    pdf_pointsize = pdf_pointsize,
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
                                  custom_pathway_gmt = NULL,
                                  heatmap_side = "left",
                                  heatmap_cluster_columns = FALSE,
                                  heatmap_col_order = NULL,
                                  gfc_scale_limits = NULL,
                                  plot = TRUE,
                                  save_pdf = TRUE,
                                  pdf_width = NULL,
                                  pdf_height = NULL,
                                  pdf_pointsize = 11,
                                  plot_per_comparison = TRUE,
                                  consistent_terms = TRUE,
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
    custom_pathway_gmt = custom_pathway_gmt,
    heatmap_side = heatmap_side,
    heatmap_cluster_columns = heatmap_cluster_columns,
    heatmap_col_order = heatmap_col_order,
    gfc_scale_limits = gfc_scale_limits,
    plot = plot,
    save_pdf = save_pdf,
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    pdf_pointsize = pdf_pointsize,
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
                                   custom_gmt_files = NULL,
                                   clusters = c("all"),
                                   mode = c("coarse", "fine"),
                                   top = 3,
                                   qval = 0.1,
                                   padj = "BH",
                                   min_term_genes = 5,
                                   min_gs_size = 10,
                                   max_gs_size = 5000,
                                   annotation_slot = c("auto", "enriched_per_cluster", "enriched_per_cluster2"),
                                   slot_suffix = NULL,
                                   clear_previous_slots = TRUE,
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
    custom_gmt_files = custom_gmt_files,
    clusters = clusters,
    mode = mode,
    top = top,
    qval = qval,
    padj = padj,
    min_term_genes = min_term_genes,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    annotation_slot = annotation_slot,
    slot_suffix = slot_suffix,
    clear_previous_slots = clear_previous_slots,
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

#' Module cell-type activity from Enrichr markers via decoupleR (S4 API)
#'
#' @inheritParams celltype_activity_decoupler
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_celltype_activity_decoupler <- function(hc,
                                           databases = c("Descartes_Cell_Types_and_Tissue_2021", "Human_Gene_Atlas"),
                                           custom_gmt_files = NULL,
                                           clusters = c("all"),
                                           mode = c("coarse", "fine"),
                                           top = 3,
                                           qval = 0.1,
                                           padj = "BH",
                                           activity_input = c("gfc", "fc", "expression"),
                                           fc_comparisons = NULL,
                                           method = "ulm",
                                           minsize = 5,
                                           min_term_genes = 5,
                                           annotation_slot = c("auto", "enriched_per_cluster", "enriched_per_cluster2"),
                                           slot_suffix = "decoupler",
                                           clear_previous_slots = FALSE,
                                           coarse_map = NULL,
                                           coarse_include_other = TRUE,
                                           refresh_db = FALSE,
                                           export_excel = TRUE,
                                           excel_file = "Module_Celltype_Activity_decoupler.xlsx",
                                           plot_heatmap = FALSE,
                                           heatmap_file_name = "Heatmap_modules_celltype_activity_decoupler.pdf",
                                           ...) {
  .hc_run_legacy(
    hc = hc,
    fun = "celltype_activity_decoupler",
    databases = databases,
    custom_gmt_files = custom_gmt_files,
    clusters = clusters,
    mode = mode,
    top = top,
    qval = qval,
    padj = padj,
    activity_input = activity_input,
    fc_comparisons = fc_comparisons,
    method = method,
    minsize = minsize,
    min_term_genes = min_term_genes,
    annotation_slot = annotation_slot,
    slot_suffix = slot_suffix,
    clear_previous_slots = clear_previous_slots,
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
                                                gfc_scale_limits = NULL,
                                                heatmap_col_order = NULL,
                                                heatmap_cluster_columns = FALSE,
                                                pdf_width = NULL,
                                                pdf_height = NULL,
                                                pdf_pointsize = 11,
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
    gfc_scale_limits = gfc_scale_limits,
    heatmap_col_order = heatmap_col_order,
    heatmap_cluster_columns = heatmap_cluster_columns,
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    pdf_pointsize = pdf_pointsize,
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
#' @param file_name Optional file name passed to `plot_cluster_heatmap()`.
#' @param ... Passed to `plot_cluster_heatmap()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_cluster_heatmap <- function(hc, file_name = NULL, ...) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  legacy_envo <- .hc_legacy_state_env()
  legacy_state <- .hc_bind_legacy_hcobject(
    .hc_as_hcobject_for_cluster_plot(hc),
    envo = legacy_envo
  )
  on.exit(.hc_restore_legacy_hcobject(legacy_state), add = TRUE)

  old_opt <- getOption("hcocena.suppress_legacy_warning", FALSE)
  options(hcocena.suppress_legacy_warning = TRUE)
  on.exit(options(hcocena.suppress_legacy_warning = old_opt), add = TRUE)

  base::do.call(plot_cluster_heatmap, c(list(file_name = file_name), list(...)))
  .hc_update_hc_from_cluster_plot(
    hc = hc,
    hcobject = base::get("hcobject", envir = legacy_envo, inherits = FALSE)
  )
}

#' Recalculate grouped GFCs and replot module heatmap (S4 API)
#'
#' Uses the existing module definitions from the current object and only
#' recalculates grouped GFCs for a different grouping variable. No integrated
#' network rebuild and no module re-clustering is performed.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param group_by Grouping column to use (must be present in all annotation tables).
#' @param col_order Optional heatmap column order.
#' @param cluster_columns Whether to cluster heatmap columns. Default is
#'   `FALSE`, so the stored main hCoCena column order is reused when available.
#' @param row_order Optional heatmap row order.
#' @param cluster_rows Whether to cluster heatmap rows.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_change_grouping_parameter <- function(hc,
                                         group_by,
                                         col_order = NULL,
                                         cluster_columns = FALSE,
                                         row_order = NULL,
                                         cluster_rows = TRUE) {
  .hc_run_legacy(
    hc = hc,
    fun = "change_grouping_parameter",
    group_by = group_by,
    col_order = col_order,
    cluster_columns = cluster_columns,
    row_order = row_order,
    cluster_rows = cluster_rows
  )
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

#' Test module differences between conditions (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `module_condition_significance()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_module_condition_significance <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = "module_condition_significance", ...)
}
