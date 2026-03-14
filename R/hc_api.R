#' Initialize an empty `HCoCenaExperiment`
#'
#' @examples
#' hc <- hc_init()
#' methods::is(hc, "HCoCenaExperiment")
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
#' @examples
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = FALSE,
#'   dir_annotation = FALSE,
#'   dir_reference_files = tempdir(),
#'   dir_output = tempdir()
#' )
#' as.data.frame(methods::slot(hc_config(hc), "paths"))
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
#' @examples
#' root <- file.path(tempdir(), "hcocena-paths")
#' dir.create(file.path(root, "reference_files"), recursive = TRUE, showWarnings = FALSE)
#' dir.create(file.path(root, "output"), recursive = TRUE, showWarnings = FALSE)
#' paths <- hc_resolve_paths(
#'   project_root_hint = root,
#'   local_count_subdir = NULL,
#'   local_annotation_subdir = NULL
#' )
#' names(paths)
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
#' @examples
#' root <- file.path(tempdir(), "hcocena-auto-paths")
#' dir.create(file.path(root, "reference_files"), recursive = TRUE, showWarnings = FALSE)
#' dir.create(file.path(root, "output"), recursive = TRUE, showWarnings = FALSE)
#' hc <- hc_init()
#' hc <- hc_auto_set_paths(
#'   hc,
#'   project_root_hint = root,
#'   local_count_subdir = NULL,
#'   local_annotation_subdir = NULL
#' )
#' as.data.frame(methods::slot(hc_config(hc), "paths"))
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
#' @rdname define_layers
#' @param hc A `HCoCenaExperiment`.
#' @param data_sets Named list of layers, each value a pair (count, annotation).
#' @examples
#' hc <- hc_init()
#' hc <- hc_define_layers(
#'   hc,
#'   data_sets = list(
#'     Layer1 = c("counts.tsv", "anno.tsv")
#'   )
#' )
#' as.data.frame(methods::slot(hc_config(hc), "layer"))
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
# Internal helpers shared by modern and bridged legacy import code.
.hc_has_output_dir <- function(hc) {
  paths_cfg <- hc@config@paths
  if (!(base::nrow(paths_cfg) > 0 && "dir_output" %in% base::colnames(paths_cfg))) {
    return(FALSE)
  }

  out_dir <- as.character(paths_cfg$dir_output[[1]])
  base::length(out_dir) == 1 &&
    !base::is.na(out_dir) &&
    base::nzchar(out_dir) &&
    !base::identical(out_dir, "FALSE")
}

.hc_auto_prepare_output <- function(hc, create_output_dir = TRUE, project_folder = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!isTRUE(.hc_has_output_dir(hc))) {
    return(hc)
  }

  hc <- hc_check_dirs(hc, create_output_dir = create_output_dir)

  save_folder <- project_folder
  if (is.null(save_folder)) {
    global_cfg <- hc@config@global
    if (base::nrow(global_cfg) > 0 && "save_folder" %in% base::colnames(global_cfg)) {
      current <- as.character(global_cfg$save_folder[[1]])
      if (base::length(current) == 1 && !base::is.na(current)) {
        save_folder <- current
      }
    }
  }

  if (is.null(save_folder) || base::length(save_folder) == 0 || base::is.na(save_folder)) {
    save_folder <- ""
  }
  save_folder <- as.character(save_folder[[1]])

  hc_init_save_folder(
    hc,
    name = save_folder,
    use_output_dir = base::identical(save_folder, "")
  )
}

.hc_normalize_layer_source <- function(src, kind, layer_idx) {
  if (!base::is.character(src) || base::length(src) != 1 || base::is.na(src) || !base::nzchar(src)) {
    stop(
      "Layer ", layer_idx, ": ", kind, " source must be a single non-empty string.\n",
      "If you provide objects, pass quoted names in `hc_define_layers()`, e.g. ",
      "`c(\"counts_df\", \"anno_df\")`."
    )
  }
  base::as.character(src)
}

.hc_legacy_dir_is_false <- function(x) {
  base::isFALSE(x) || (base::is.character(x) && base::identical(x, "FALSE"))
}

.hc_resolve_source_path <- function(dir_path, source_name) {
  if (base::file.exists(source_name)) {
    return(source_name)
  }
  if (.hc_legacy_dir_is_false(dir_path)) {
    return(source_name)
  }
  base::file.path(base::as.character(dir_path), source_name)
}

.hc_normalize_count_object <- function(obj, source_name, gene_symbol_col = NULL) {
  if (base::is.matrix(obj)) {
    obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!base::is.data.frame(obj)) {
    stop("Count source `", source_name, "` exists but is not a data.frame or matrix.")
  }
  obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)

  if (!base::is.null(gene_symbol_col) && gene_symbol_col %in% base::colnames(obj)) {
    obj <- make_rownames_unique(counts = obj, gene_symbol_col = gene_symbol_col)
  } else {
    if (base::is.null(base::rownames(obj)) || base::any(!base::nzchar(base::rownames(obj)))) {
      stop(
        "Count object `", source_name, "` has no usable rownames.\n",
        "Provide `gene_symbol_col` (matching a column in the count object) ",
        "or set rownames to gene symbols before calling `hc_read_data()`."
      )
    }
    numeric_cols <- base::vapply(obj, base::is.numeric, logical(1))
    if (!base::all(numeric_cols)) {
      obj <- obj[, numeric_cols, drop = FALSE]
    }
    if (base::ncol(obj) == 0) {
      stop("Count object `", source_name, "` has no numeric sample columns.")
    }
  }

  obj
}

.hc_normalize_anno_object <- function(obj, source_name, sample_col = NULL) {
  if (base::is.matrix(obj)) {
    obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!base::is.data.frame(obj)) {
    stop("Annotation source `", source_name, "` exists but is not a data.frame or matrix.")
  }
  obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)

  if (!base::is.null(sample_col) && sample_col %in% base::colnames(obj)) {
    base::rownames(obj) <- base::as.character(obj[[sample_col]])
  }
  if (base::is.null(base::rownames(obj)) || base::any(!base::nzchar(base::rownames(obj)))) {
    stop(
      "Annotation object `", source_name, "` has no usable rownames.\n",
      "Provide `sample_col` (matching a column in the annotation object) ",
      "or set rownames to sample IDs before calling `hc_read_data()`."
    )
  }

  obj[] <- base::lapply(obj, base::factor)
  obj
}

.hc_load_count_source <- function(source_name,
                                  dir_count_data,
                                  gene_symbol_col,
                                  count_has_rn,
                                  sep_counts,
                                  source_env) {
  if (base::exists(source_name, envir = source_env, inherits = TRUE)) {
    return(.hc_normalize_count_object(
      base::get(source_name, envir = source_env, inherits = TRUE),
      source_name = source_name,
      gene_symbol_col = gene_symbol_col
    ))
  }

  if (base::is.null(gene_symbol_col)) {
    stop("You must provide the 'gene_symbol_col' parameter.")
  }

  count_file <- .hc_resolve_source_path(dir_count_data, source_name)
  if (!base::file.exists(count_file)) {
    stop(
      "Count source `", source_name, "` was not found as an object, and file `",
      count_file, "` does not exist.\n",
      "If you use object input, pass quoted object names in `hc_define_layers()`."
    )
  }

  read_expression_data(
    file = count_file,
    rown = count_has_rn,
    sep = sep_counts,
    gene_symbol_col = gene_symbol_col
  )
}

.hc_load_anno_source <- function(source_name,
                                 dir_annotation,
                                 sample_col,
                                 anno_has_rn,
                                 sep_anno,
                                 source_env) {
  if (base::exists(source_name, envir = source_env, inherits = TRUE)) {
    return(.hc_normalize_anno_object(
      base::get(source_name, envir = source_env, inherits = TRUE),
      source_name = source_name,
      sample_col = sample_col
    ))
  }

  if (base::is.null(sample_col)) {
    stop("You must provide the 'sample_col' parameter.")
  }

  anno_file <- .hc_resolve_source_path(dir_annotation, source_name)
  if (!base::file.exists(anno_file)) {
    stop(
      "Annotation source `", source_name, "` was not found as an object, and file `",
      anno_file, "` does not exist.\n",
      "If you use object input, pass quoted object names in `hc_define_layers()`."
    )
  }

  read_anno(
    file = anno_file,
    rown = anno_has_rn,
    sep = sep_anno,
    sample_col = sample_col
  )
}

.hc_read_data_impl <- function(hc,
                               sep_counts = "\t",
                               sep_anno = "\t",
                               gene_symbol_col = NULL,
                               sample_col = NULL,
                               count_has_rn = TRUE,
                               anno_has_rn = TRUE,
                               auto_setup_output = TRUE,
                               create_output_dir = TRUE,
                               project_folder = NULL,
                               source_env = parent.frame()) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (isTRUE(auto_setup_output)) {
    hc <- .hc_auto_prepare_output(
      hc,
      create_output_dir = create_output_dir,
      project_folder = project_folder
    )
  }

  layer_cfg <- hc@config@layer
  required_cols <- c("layer_id", "count_source", "annotation_source")
  if (!(base::nrow(layer_cfg) > 0 && base::all(required_cols %in% base::colnames(layer_cfg)))) {
    stop("No valid layers found. Run `hc_define_layers()` before `hc_read_data()`.")
  }

  paths <- .hc_row_to_list(hc@config@paths)
  data <- list()
  for (i in base::seq_len(base::nrow(layer_cfg))) {
    lid <- as.character(layer_cfg$layer_id[[i]])
    count_source <- .hc_normalize_layer_source(layer_cfg$count_source[[i]], "count", i)
    anno_source <- .hc_normalize_layer_source(layer_cfg$annotation_source[[i]], "annotation", i)

    counts <- .hc_load_count_source(
      source_name = count_source,
      dir_count_data = paths[["dir_count_data"]],
      gene_symbol_col = gene_symbol_col,
      count_has_rn = count_has_rn,
      sep_counts = sep_counts,
      source_env = source_env
    )

    var.df <- rank_variance(counts)
    if (base::all(base::is.na(var.df$variance))) {
      stop(
        "Count data for layer ", i, " contains no valid numeric expression values after preprocessing."
      )
    }
    zero_var <- var.df$gene[!base::is.na(var.df$variance) & var.df$variance == 0]
    if (base::length(zero_var) > 0) {
      message(base::paste0("Detected genes with 0 variance in dataset ", i, "."))
      counts <- counts[!(base::rownames(counts) %in% zero_var), , drop = FALSE]
      message(base::paste0(base::length(zero_var), " gene(s) were removed from dataset ", i, "."))
    }

    anno <- .hc_load_anno_source(
      source_name = anno_source,
      dir_annotation = paths[["dir_annotation"]],
      sample_col = sample_col,
      anno_has_rn = anno_has_rn,
      sep_anno = sep_anno,
      source_env = source_env
    )

    if (!base::ncol(counts) == base::nrow(anno)) {
      stop(
        base::paste0(
          "The count table has ", base::ncol(counts), " columns but the annotation has ",
          base::nrow(anno), " rows. These values are required to be the same since they\n",
          "                 should correspond to the number of samples. THE LOADING OF THE DATA WILL BE TERMINATED."
        )
      )
    }
    if (!base::all(base::as.character(base::colnames(counts)) %in% base::as.character(base::rownames(anno)))) {
      stop("The column names of the count file do not all match the rownames of the annotation. Please make sure they contain the same samples.")
    }
    anno <- anno[base::colnames(counts), , drop = FALSE]

    data[[base::paste0(lid, "_counts")]] <- counts
    data[[base::paste0(lid, "_anno")]] <- anno
  }

  hc@mae <- .hc_build_mae(list(data = data), layer_cfg)
  methods::validObject(hc)
  hc
}
#
#' @rdname read_data
#' @inheritParams read_data
#' @param hc A `HCoCenaExperiment`.
#' @param auto_setup_output Boolean. If `TRUE`, run `hc_check_dirs()` and
#'   initialize the save folder before reading data (if `dir_output` is set).
#' @param create_output_dir Boolean passed to `hc_check_dirs()` when
#'   `auto_setup_output = TRUE`.
#' @param project_folder Optional save subfolder name. Use `""` to write
#'   directly into `dir_output`. If `NULL`, keep an already configured
#'   `save_folder` or default to `""`.
#' @examples
#' extdir <- paste0(
#'   normalizePath(system.file("extdata", package = "hcocena"), winslash = "/"),
#'   "/"
#' )
#' outdir <- file.path(tempdir(), "hcocena-read-data")
#' dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = extdir,
#'   dir_annotation = extdir,
#'   dir_reference_files = extdir,
#'   dir_output = outdir
#' )
#' hc <- hc_define_layers(
#'   hc,
#'   data_sets = list(
#'     Layer1 = c("toy_layer1_counts.tsv", "toy_layer1_anno.tsv"),
#'     Layer2 = c("toy_layer2_counts.tsv", "toy_layer2_anno.tsv")
#'   )
#' )
#' hc <- hc_read_data(
#'   hc,
#'   sep_counts = "\t",
#'   sep_anno = "\t",
#'   gene_symbol_col = "SYMBOL",
#'   sample_col = "SampleID",
#'   count_has_rn = FALSE,
#'   anno_has_rn = FALSE
#' )
#' names(MultiAssayExperiment::experiments(hc_mae(hc)))
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
  .hc_read_data_impl(
    hc = hc,
    sep_counts = sep_counts,
    sep_anno = sep_anno,
    gene_symbol_col = gene_symbol_col,
    sample_col = sample_col,
    count_has_rn = count_has_rn,
    anno_has_rn = anno_has_rn,
    auto_setup_output = auto_setup_output,
    create_output_dir = create_output_dir,
    project_folder = project_folder,
    source_env = parent.frame()
  )
}

#' Set global settings (S4 API)
#'
# Internal implementation shared by S4 and legacy entry points.
.hc_set_global_settings_impl <- function(hc,
                                         organism = "human",
                                         control_keyword = "none",
                                         variable_of_interest = "merged",
                                         min_nodes_number_for_network = 50,
                                         min_nodes_number_for_cluster = 50,
                                         range_GFC = 2.0,
                                         layout_algorithm = "layout_with_stress",
                                         data_in_log = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  layout_algorithm <- .hc_normalize_layout_algorithm(layout_algorithm)
  if (layout_algorithm %in% c("layout_with_stress", "layout_with_sparse_stress") &&
      !requireNamespace("graphlayouts", quietly = TRUE)) {
    warning(
      "Selected `layout_algorithm = '", layout_algorithm, "'`, but package `graphlayouts` is not installed.\n",
      "Network plotting will fall back to `layout_with_fr` until `graphlayouts` is available.",
      call. = FALSE
    )
  }

  global_cfg <- .hc_row_to_list(hc@config@global)
  global_cfg[["organism"]] <- organism
  global_cfg[["control"]] <- control_keyword
  global_cfg[["voi"]] <- variable_of_interest
  global_cfg[["min_nodes_number_for_network"]] <- min_nodes_number_for_network
  global_cfg[["min_nodes_number_for_cluster"]] <- min_nodes_number_for_cluster
  global_cfg[["range_GFC"]] <- range_GFC
  global_cfg[["layout_algorithm"]] <- layout_algorithm
  global_cfg[["data_in_log"]] <- data_in_log

  hc@config@global <- .hc_to_data_frame(global_cfg)
  methods::validObject(hc)
  hc
}
#
#' @rdname set_global_settings
#' @inheritParams set_global_settings
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' extdir <- paste0(
#'   normalizePath(system.file("extdata", package = "hcocena"), winslash = "/"),
#'   "/"
#' )
#' outdir <- file.path(tempdir(), "hcocena-global-settings")
#' dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = extdir,
#'   dir_annotation = extdir,
#'   dir_reference_files = extdir,
#'   dir_output = outdir
#' )
#' hc <- hc_define_layers(
#'   hc,
#'   data_sets = list(
#'     Layer1 = c("toy_layer1_counts.tsv", "toy_layer1_anno.tsv"),
#'     Layer2 = c("toy_layer2_counts.tsv", "toy_layer2_anno.tsv")
#'   )
#' )
#' hc <- hc_read_data(
#'   hc,
#'   sep_counts = "\t",
#'   sep_anno = "\t",
#'   gene_symbol_col = "SYMBOL",
#'   sample_col = "SampleID",
#'   count_has_rn = FALSE,
#'   anno_has_rn = FALSE
#' )
#' hc <- hc_set_global_settings(
#'   hc,
#'   organism = "human",
#'   control_keyword = "control",
#'   variable_of_interest = "group",
#'   data_in_log = TRUE
#' )
#' as.data.frame(methods::slot(hc_config(hc), "global"))
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
  .hc_set_global_settings_impl(
    hc = hc,
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
# Internal implementation shared by S4 and legacy entry points.
.hc_recycle_layer_arg <- function(x, n_layers, arg_name) {
  if (base::length(x) == 1 && n_layers > 1) {
    return(base::rep(x, n_layers))
  }
  if (base::length(x) != n_layers) {
    stop("`", arg_name, "` must have length 1 or match the number of layers (", n_layers, ").")
  }
  x
}

.hc_set_layer_settings_impl <- function(hc,
                                        top_var,
                                        min_corr = 0.7,
                                        range_cutoff_length,
                                        print_distribution_plots = FALSE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (base::nrow(hc@config@layer) == 0 || !("layer_id" %in% base::colnames(hc@config@layer))) {
    stop("No layers found. Run `hc_define_layers()` before `hc_set_layer_settings()`.")
  }

  n_layers <- base::nrow(hc@config@layer)
  top_var <- .hc_recycle_layer_arg(top_var, n_layers, "top_var")
  min_corr <- .hc_recycle_layer_arg(min_corr, n_layers, "min_corr")
  range_cutoff_length <- .hc_recycle_layer_arg(range_cutoff_length, n_layers, "range_cutoff_length")
  print_distribution_plots <- .hc_recycle_layer_arg(
    print_distribution_plots,
    n_layers,
    "print_distribution_plots"
  )

  rows <- vector("list", n_layers)
  for (i in base::seq_len(n_layers)) {
    row <- list()
    for (nm in base::colnames(hc@config@layer)) {
      row[[nm]] <- hc@config@layer[[nm]][[i]]
    }
    row[["top_var"]] <- top_var[[i]]
    row[["min_corr"]] <- min_corr[[i]]
    row[["range_cutoff_length"]] <- range_cutoff_length[[i]]
    row[["print_distribution_plots"]] <- print_distribution_plots[[i]]
    rows[[i]] <- row
  }

  hc@config@layer <- .hc_rows_to_data_frame(rows)
  methods::validObject(hc)
  hc
}
#
#' @rdname set_layer_settings
#' @inheritParams set_layer_settings
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' extdir <- paste0(
#'   normalizePath(system.file("extdata", package = "hcocena"), winslash = "/"),
#'   "/"
#' )
#' outdir <- file.path(tempdir(), "hcocena-layer-settings")
#' dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = extdir,
#'   dir_annotation = extdir,
#'   dir_reference_files = extdir,
#'   dir_output = outdir
#' )
#' hc <- hc_define_layers(
#'   hc,
#'   data_sets = list(
#'     Layer1 = c("toy_layer1_counts.tsv", "toy_layer1_anno.tsv"),
#'     Layer2 = c("toy_layer2_counts.tsv", "toy_layer2_anno.tsv")
#'   )
#' )
#' hc <- hc_read_data(
#'   hc,
#'   sep_counts = "\t",
#'   sep_anno = "\t",
#'   gene_symbol_col = "SYMBOL",
#'   sample_col = "SampleID",
#'   count_has_rn = FALSE,
#'   anno_has_rn = FALSE
#' )
#' hc <- hc_set_layer_settings(
#'   hc,
#'   top_var = c(2, 2),
#'   min_corr = 0.1,
#'   range_cutoff_length = c(2, 2),
#'   print_distribution_plots = FALSE
#' )
#' as.data.frame(methods::slot(hc_config(hc), "layer"))
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_set_layer_settings <- function(hc,
                                  top_var,
                                  min_corr = 0.7,
                                  range_cutoff_length,
                                  print_distribution_plots = FALSE) {
  .hc_set_layer_settings_impl(
    hc = hc,
    top_var = top_var,
    min_corr = min_corr,
    range_cutoff_length = range_cutoff_length,
    print_distribution_plots = print_distribution_plots
  )
}

#' Register supplementary references (S4 API)
#'
# Internal implementation shared by S4 and legacy entry points.
.hc_set_supp_files_impl <- function(hc,
                                    Tf = NULL,
                                    Hallmark = NULL,
                                    Go = NULL,
                                    Kegg = NULL,
                                    Reactome = NULL,
                                    ...) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  hc@references@registry <- .hc_reference_registry(
    base::list(Tf = Tf, Hallmark = Hallmark, Go = Go, Kegg = Kegg, Reactome = Reactome, ...)
  )
  methods::validObject(hc)
  hc
}
#
#' @rdname set_supp_files
#' @inheritParams set_supp_files
#' @param hc A `HCoCenaExperiment`.
#' @param ... Additional supplementary files passed through to `set_supp_files()`.
#' @examples
#' hc <- hc_init()
#' hc <- hc_set_supp_files(
#'   hc,
#'   Hallmark = "hallmark.gmt",
#'   Go = "go.gmt"
#' )
#' sort(names(as_hcobject(hc)$supplement))
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_set_supp_files <- function(hc, Tf = NULL, Hallmark = NULL, Go = NULL, Kegg = NULL, Reactome = NULL, ...) {
  .hc_set_supp_files_impl(
    hc = hc,
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
# Internal implementation shared by S4 and legacy entry points.
.hc_read_supplementary_target_name <- function(name) {
  if (base::identical(name, "Tf")) {
    return("TF")
  }
  if (name %in% c("Hallmark", "Go", "Kegg", "Reactome")) {
    return(name)
  }
  stringr::str_to_title(name)
}

.hc_read_supplementary_one <- function(name, source, dir_reference_files) {
  source_path <- .hc_resolve_source_path(dir_reference_files, source)
  if (!base::file.exists(source_path)) {
    stop("Supplementary source `", source, "` for `", name, "` does not exist.")
  }

  if (base::identical(name, "Tf")) {
    return(utils::read.delim(source_path, header = TRUE, check.names = FALSE))
  }
  if (name %in% c("Hallmark", "Go", "Kegg", "Reactome")) {
    return(clusterProfiler::read.gmt(source_path))
  }
  if (base::grepl("\\.csv$", source, ignore.case = TRUE)) {
    return(utils::read.csv(source_path, header = TRUE, stringsAsFactors = FALSE, quote = ""))
  }
  if (base::grepl("\\.gmt$", source, ignore.case = TRUE)) {
    return(clusterProfiler::read.gmt(source_path))
  }

  print(base::paste0("invalid input format of database: ", name, ". Valid inputs are: .csv and .gmt files!"))
  NULL
}

.hc_read_supplementary_impl <- function(hc) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  registry <- hc@references@registry
  if (!(base::nrow(registry) > 0 && base::all(c("name", "source") %in% base::colnames(registry)))) {
    return(hc)
  }

  paths <- .hc_row_to_list(hc@config@paths)
  dir_reference_files <- paths[["dir_reference_files"]]
  if (base::is.null(dir_reference_files) || .hc_legacy_dir_is_false(dir_reference_files)) {
    stop("`dir_reference_files` is not set. Please run `init_wd()`/`hc_set_paths()` first.")
  }

  loaded <- as.list(hc@references@data)
  for (i in base::seq_len(base::nrow(registry))) {
    name <- as.character(registry$name[[i]])
    source <- as.character(registry$source[[i]])
    target_name <- .hc_read_supplementary_target_name(name)
    value <- .hc_read_supplementary_one(name, source, dir_reference_files)
    if (!base::is.null(value)) {
      loaded[[target_name]] <- value
    }
  }

  hc@references@data <- S4Vectors::SimpleList(loaded)
  methods::validObject(hc)
  hc
}
#
#' @rdname read_supplementary
#' @param hc A `HCoCenaExperiment`.
#' @examples
#' refdir <- file.path(tempdir(), "hcocena-reference-files")
#' dir.create(refdir, recursive = TRUE, showWarnings = FALSE)
#' writeLines("ToySet\\tdescription\\tG1\\tG2", file.path(refdir, "hallmark.gmt"))
#' hc <- hc_init()
#' hc <- hc_set_paths(
#'   hc,
#'   dir_count_data = FALSE,
#'   dir_annotation = FALSE,
#'   dir_reference_files = paste0(normalizePath(refdir, winslash = "/"), "/"),
#'   dir_output = tempdir()
#' )
#' hc <- hc_set_supp_files(hc, Hallmark = "hallmark.gmt")
#' hc <- hc_read_supplementary(hc)
#' names(as_hcobject(hc)$supplementary_data)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_read_supplementary <- function(hc) {
  .hc_read_supplementary_impl(hc)
}

#' Run expression analysis part I (S4 API)
#'
# Internal implementation shared by S4 and legacy entry points.
.hc_run_expression_analysis_1_legacy_driver <- function(padj,
                                                        export,
                                                        import,
                                                        bayes,
                                                        prior,
                                                        alpha,
                                                        corr_method) {
  for (x in base::seq_len(base::length(hcobject[["layers"]]))) {
    hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]] <<- run_expression_analysis_1_body(
      x = x,
      bayes = bayes,
      prior = prior,
      alpha = alpha,
      padj = padj,
      export = export,
      import = import,
      corr_method = corr_method
    )
  }
}

.hc_run_expression_analysis_1_impl <- function(hc,
                                               padj = "none",
                                               export = FALSE,
                                               import = NULL,
                                               bayes = FALSE,
                                               prior = 2,
                                               alpha = 0.5,
                                               corr_method = "pearson") {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!corr_method %in% c("pearson", "spearman", "rho")) {
    stop("Parameter 'corr_method' must be either 'pearson', 'spearman' or 'rho'.")
  }
  if (!(base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer))) {
    stop("No layers found. Run `hc_define_layers()` before `hc_run_expression_analysis_1()`.")
  }

  .hc_run_legacy(
    hc = hc,
    fun = .hc_run_expression_analysis_1_legacy_driver,
    padj = padj,
    export = export,
    import = import,
    bayes = bayes,
    prior = prior,
    alpha = alpha,
    corr_method = corr_method
  )
}
#
#' @rdname run_expression_analysis_1
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
  .hc_run_expression_analysis_1_impl(
    hc = hc,
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
#' @rdname set_cutoff
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
.hc_set_cutoff_legacy_driver <- function(cutoff_vector) {
  hcobject[["cutoff_vec"]] <<- cutoff_vector
}

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
    fun = .hc_set_cutoff_legacy_driver,
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
# Internal implementation shared by S4 and legacy entry points.
.hc_run_expression_analysis_2_legacy_driver <- function(grouping_v,
                                                        plot_HM,
                                                        method,
                                                        additional_anno,
                                                        cols) {
  for (x in base::seq_len(base::length(hcobject[["layers"]]))) {
    extra_anno <- if (base::is.null(additional_anno) || base::length(additional_anno) < x) {
      NULL
    } else {
      additional_anno[[x]]
    }

    run_expression_analysis_2_body(
      x = x,
      grouping_v = grouping_v,
      plot_HM = plot_HM,
      method = method,
      additional_anno = extra_anno,
      title = hcobject[["layers_names"]][x],
      cols = cols
    )
  }
}

.hc_run_expression_analysis_2_impl <- function(hc,
                                               grouping_v = NULL,
                                               plot_HM = TRUE,
                                               method = "complete",
                                               additional_anno = NULL,
                                               cols = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!(base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer))) {
    stop("No layers found. Run `hc_define_layers()` before `hc_run_expression_analysis_2()`.")
  }

  .hc_run_legacy(
    hc = hc,
    fun = .hc_run_expression_analysis_2_legacy_driver,
    grouping_v = grouping_v,
    plot_HM = plot_HM,
    method = method,
    additional_anno = additional_anno,
    cols = cols
  )
}
#
#' @rdname run_expression_analysis_2
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
  .hc_run_expression_analysis_2_impl(
    hc = hc,
    grouping_v = grouping_v,
    plot_HM = plot_HM,
    method = method,
    additional_anno = additional_anno,
    cols = cols
  )
}

#' Build integrated network (S4 API)
#'
# Internal implementation shared by S4 and legacy entry points.
.hc_build_integrated_network_legacy_driver <- function(mode,
                                                       with,
                                                       multi_edges,
                                                       GFC_when_missing) {
  if (mode == "u") {
    message("Intergrating network based on union.")
    get_union()
  } else if (mode == "i") {
    message("Intergrating network based on intersection.")
    if (is.null(with)) {
      stop("The 'with' parameter must be specified.")
    }
    get_intersection(with = with)
  } else {
    stop("No valid choice of 'mode'. Must be either 'u' for integration by union or 'i' for integration by intersection.")
  }

  merged_net <- igraph::graph_from_data_frame(
    hcobject[["integrated_output"]][["combined_edgelist"]],
    directed = FALSE
  )
  merged_net <- igraph::simplify(
    merged_net,
    edge.attr.comb = list(weight = multi_edges, "ignore")
  )
  new_edgelist <- base::cbind(
    igraph::get.edgelist(merged_net),
    base::round(igraph::E(merged_net)$weight, 7)
  ) %>%
    base::as.data.frame()
  base::colnames(new_edgelist) <- base::colnames(
    hcobject[["integrated_output"]][["combined_edgelist"]]
  )
  hcobject[["integrated_output"]][["combined_edgelist"]] <<- new_edgelist
  hcobject[["integrated_output"]][["merged_net"]] <<- merged_net
  hcobject[["integrated_output"]][["GFC_all_layers"]] <<- merge_GFCs(
    GFC_when_missing = GFC_when_missing
  )

  # Rebuilding the integrated graph invalidates prior cluster-dependent outputs.
  hcobject[["integrated_output"]][["cluster_calc"]] <<- list()
  hcobject[["integrated_output"]][["enrichments"]] <<- NULL
  hcobject[["integrated_output"]][["upstream_inference"]] <<- NULL
  hcobject[["integrated_output"]][["knowledge_network"]] <<- NULL
  hcobject[["satellite_outputs"]][["enrichments"]] <<- NULL
  hcobject[["satellite_outputs"]][["upstream_inference"]] <<- NULL
  hcobject[["satellite_outputs"]][["knowledge_network"]] <<- NULL
  hcobject[["satellite_outputs"]][["labelled_network"]] <<- NULL
  hcobject[["satellite_outputs"]][["network_col_by_module"]] <<- NULL
}
#
#' @rdname build_integrated_network
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
    fun = .hc_build_integrated_network_legacy_driver,
    mode = mode,
    with = with,
    multi_edges = multi_edges,
    GFC_when_missing = GFC_when_missing
  )
}

#' Cluster integrated network (S4 API)
#'
# Internal implementation shared by S4 and legacy entry points.
.hc_cluster_calculation_legacy_driver <- function(cluster_algo,
                                                  no_of_iterations,
                                                  resolution,
                                                  partition_type,
                                                  max_cluster_count_per_gene,
                                                  return_result) {
  hcobject[["global_settings"]][["chosen_clustering_algo"]] <<- cluster_algo

  if (cluster_algo == "cluster_leiden") {
    hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- leiden_clustering(
      g = hcobject[["integrated_output"]][["merged_net"]],
      num_it = no_of_iterations,
      resolution = resolution,
      partition_type = partition_type
    )
    hcobject[["integrated_output"]][["cluster_calc"]][["labelled_network"]] <<- NULL
    hcobject[["integrated_output"]][["cluster_calc"]][["network_col_by_module"]] <<- NULL
    hcobject[["satellite_outputs"]][["labelled_network"]] <<- NULL
    hcobject[["satellite_outputs"]][["network_col_by_module"]] <<- NULL
    return(invisible(NULL))
  }

  cluster_algo_list <- c(
    "cluster_label_prop",
    "cluster_fast_greedy",
    "cluster_louvain",
    "cluster_infomap",
    "cluster_walktrap",
    "cluster_leiden"
  )

  if (cluster_algo == "auto") {
    algos_to_use <- cluster_algo_list
  } else {
    algos_to_use <- cluster_algo
  }

  cluster_run <- .hc_with_seed(168575L, {
    if (cluster_algo == "auto") {
      print(algos_to_use)
      df_modularity_score <- base::do.call(
        "rbind",
        base::lapply(algos_to_use, function(algo_now) {
          cluster_calculation_internal(
            graph_obj = hcobject[["integrated_output"]][["merged_net"]],
            algo = algo_now,
            case = "test",
            resolution = resolution,
            partition_type = partition_type
          )
        })
      )

      cluster_algo_used <- df_modularity_score %>%
        dplyr::filter(modularity_score == base::max(modularity_score)) %>%
        dplyr::select(cluster_algorithm) %>%
        base::as.character()

      print(base::paste(cluster_algo_used, " will be used based on the highest modularity score."))
    } else {
      cluster_algo_used <- cluster_algo
      print(base::paste(cluster_algo_used, " will be used based on your input."))
    }

    gene_which_cluster <- base::do.call(
      "cbind",
      base::lapply(base::seq_len(no_of_iterations), function(iteration_now) {
        cluster_calculation_internal(
          graph_obj = hcobject[["integrated_output"]][["merged_net"]],
          algo = cluster_algo_used,
          case = "best",
          resolution = resolution,
          partition_type = partition_type,
          it = iteration_now
        )
      })
    )

    list(
      cluster_algo_used = cluster_algo_used,
      gene_which_cluster = gene_which_cluster
    )
  })
  cluster_algo_used <- cluster_run$cluster_algo_used
  gene_which_cluster <- cluster_run$gene_which_cluster

  if (base::ncol(gene_which_cluster) > 1) {
    gene_cluster_ident <- base::apply(gene_which_cluster, 1, function(x) {
      if (base::length(base::unique(x)) > max_cluster_count_per_gene) {
        0
      } else {
        base::names(base::which(base::table(x) == base::max(base::table(x))))[1]
      }
    })
  } else {
    gene_cluster_ident <- gene_which_cluster[, 1]
  }

  white_genes_clustercounts <- base::as.integer(
    base::length(base::grep(gene_cluster_ident, pattern = "\\b0\\b"))
  )

  print(
    base::paste(
      white_genes_clustercounts,
      "genes were assigned to more than",
      max_cluster_count_per_gene,
      "cluster(s). These genes are assigned to Cluster 0 (white) and will be left out of the network and further analyses."
    )
  )

  cluster_data <- base::data.frame(
    genes = igraph::vertex_attr(hcobject[["integrated_output"]][["merged_net"]], "name"),
    clusters = base::paste0("Cluster ", gene_cluster_ident),
    stringsAsFactors = FALSE
  )

  dfk <- cluster_data %>%
    dplyr::count(clusters, genes) %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(
      gene_no = base::sum(n),
      gene_n = base::paste0(genes, collapse = ",")
    ) %>%
    dplyr::mutate(
      cluster_included = base::ifelse(
        gene_no >= hcobject[["global_settings"]][["min_nodes_number_for_cluster"]],
        "yes",
        "no"
      ),
      color = "white"
    )

  color.cluster <- get_cluster_colours()
  plot_clusters <- ggplot2::ggplot(
    data = dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", ],
    ggplot2::aes(x = clusters)
  ) +
    ggplot2::geom_bar(ggplot2::aes(fill = clusters)) +
    ggplot2::scale_fill_manual(values = color.cluster)

  plot_clust <- ggplot2::ggplot_build(plot_clusters)
  dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", "color"] <- plot_clust$data[[1]][["fill"]]

  white_genes_clustersize <- base::as.integer(
    dfk %>%
      dplyr::filter(cluster_included == "no") %>%
      dplyr::summarise(n = base::sum(gene_no)) %>%
      purrr::map(1)
  )

  print(
    base::paste0(
      white_genes_clustersize,
      " genes were assigned to clusters with a smaller size than the defined minimal cluster size of ",
      hcobject[["global_settings"]][["min_nodes_number_for_cluster"]],
      " genes per cluster. These genes will also be assigned to Cluster 0 (white) and left out of the network and further analyses."
    )
  )

  dfk_allinfo <- base::do.call(
    "rbind",
    base::lapply(
      1:base::nrow(dfk),
      gfc_mean_clustergene,
      cluster_df = dfk,
      gfc_dat = hcobject[["integrated_output"]][["GFC_all_layers"]]
    )
  )
  dfk_allinfo$vertexsize <- base::ifelse(dfk_allinfo$cluster_included == "yes", 3, 1)

  if (isTRUE(return_result)) {
    return(dfk_allinfo)
  }

  hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- dfk_allinfo
  hcobject[["integrated_output"]][["cluster_calc"]][["labelled_network"]] <<- NULL
  hcobject[["integrated_output"]][["cluster_calc"]][["network_col_by_module"]] <<- NULL
  hcobject[["satellite_outputs"]][["labelled_network"]] <<- NULL
  hcobject[["satellite_outputs"]][["network_col_by_module"]] <<- NULL
  invisible(NULL)
}

.hc_cluster_calculation_impl <- function(hc,
                                         cluster_algo = "cluster_leiden",
                                         no_of_iterations = 2,
                                         resolution = 0.1,
                                         partition_type = "RBConfigurationVertexPartition",
                                         max_cluster_count_per_gene = 1,
                                         return_result = FALSE) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = .hc_cluster_calculation_legacy_driver,
    cluster_algo = cluster_algo,
    no_of_iterations = no_of_iterations,
    resolution = resolution,
    partition_type = partition_type,
    max_cluster_count_per_gene = max_cluster_count_per_gene,
    return_result = return_result
  )

  if (isTRUE(return_result) && !base::is.null(out$result)) {
    return(list(hc = out$hc, result = out$result))
  }

  out$hc
}
#
#' @rdname cluster_calculation
#' @inheritParams cluster_calculation
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`. If `return_result = TRUE` and the
#'   clustering backend returns a cluster table, that table is returned
#'   instead.
#' @export
hc_cluster_calculation <- function(hc,
                                   cluster_algo = "cluster_leiden",
                                   no_of_iterations = 2,
                                   resolution = 0.1,
                                   partition_type = "RBConfigurationVertexPartition",
                                   max_cluster_count_per_gene = 1,
                                   return_result = FALSE) {
  out <- .hc_cluster_calculation_impl(
    hc = hc,
    cluster_algo = cluster_algo,
    no_of_iterations = no_of_iterations,
    resolution = resolution,
    partition_type = partition_type,
    max_cluster_count_per_gene = max_cluster_count_per_gene,
    return_result = return_result
  )

  if (base::is.list(out) &&
      "hc" %in% base::names(out) &&
      inherits(out[["hc"]], "HCoCenaExperiment")) {
    return(out[["result"]])
  }

  out
}

#' Merge modules based on module-heatmap similarity (S4 API)
#'
.hc_merge_clusters_impl <- function(hc,
                                    k = "auto",
                                    save = TRUE,
                                    method = "complete",
                                    k_min = 2,
                                    k_max = NULL,
                                    auto_parsimony_penalty = 1e-04,
                                    verbose = TRUE) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = .hc_merge_clusters_legacy_driver,
    k = k,
    save = save,
    method = method,
    k_min = k_min,
    k_max = k_max,
    auto_parsimony_penalty = auto_parsimony_penalty,
    verbose = verbose
  )

  if (!base::is.null(out$result)) {
    return(list(hc = out$hc, result = out$result))
  }

  out$hc
}
#
#' @rdname merge_clusters
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
  out <- .hc_merge_clusters_impl(
    hc = hc,
    k = k,
    save = save,
    method = method,
    k_min = k_min,
    k_max = k_max,
    auto_parsimony_penalty = auto_parsimony_penalty,
    verbose = verbose
  )

  if (base::is.list(out) &&
      "hc" %in% base::names(out) &&
      inherits(out[["hc"]], "HCoCenaExperiment")) {
    return(out[["hc"]])
  }

  out
}

#' Split one or multiple modules into submodules (S4 API)
#'
#' @rdname split_modules
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
    fun = .hc_split_modules_legacy_driver,
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
#' @rdname unsplit_modules
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
    fun = .hc_unsplit_modules_legacy_driver,
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
# Internal implementation shared by S4 and legacy entry points.
.hc_functional_enrichment_impl <- function(hc, ...) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = ".hc_functional_enrichment_legacy_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
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
  .hc_functional_enrichment_impl(
    hc = hc,
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
  )[["hc"]]
}

#' Upstream regulator/pathway inference (S4 API)
#'
#' @inheritParams upstream_inference
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
# Internal implementation shared by S4 and legacy entry points.
.hc_upstream_inference_impl <- function(hc, ...) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = ".hc_upstream_inference_legacy_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
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
  .hc_upstream_inference_impl(
    hc = hc,
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
  )[["hc"]]
}

#' Module cell-type annotation from Enrichr (S4 API)
#'
#' @inheritParams celltype_annotation
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
# Internal implementation shared by S4 and legacy entry points.
.hc_celltype_annotation_impl <- function(hc, ...) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = ".hc_celltype_annotation_legacy_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
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
  .hc_celltype_annotation_impl(
    hc = hc,
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
  )[["hc"]]
}

#' Module cell-type activity from Enrichr markers via decoupleR (S4 API)
#'
#' @inheritParams celltype_activity_decoupler
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
# Internal implementation shared by S4 and legacy entry points.
.hc_celltype_activity_decoupler_impl <- function(hc, ...) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = ".hc_celltype_activity_decoupler_legacy_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
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
  .hc_celltype_activity_decoupler_impl(
    hc = hc,
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
  )[["hc"]]
}

#' Plot module knowledge network from enrichment + upstream inference (S4 API)
#'
#' @inheritParams plot_enrichment_upstream_network
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
# Internal implementation shared by S4 and legacy entry points.
.hc_plot_enrichment_upstream_network_impl <- function(hc, ...) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = ".hc_plot_enrichment_upstream_network_legacy_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
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
  .hc_plot_enrichment_upstream_network_impl(
    hc = hc,
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
  )[["hc"]]
}

#' Run `check_dirs()` with S4 state synchronization
#'
# Internal implementation shared by S4 and legacy entry points.
.hc_check_dirs_impl <- function(hc, create_output_dir = TRUE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }

  paths <- .hc_row_to_list(hc@config@paths)
  if (base::length(paths) == 0) {
    return(hc)
  }

  for (nm in base::names(paths)) {
    current_dir <- paths[[nm]]
    if (nm == "dir_output" &&
        !base::identical(current_dir, FALSE) &&
        isTRUE(create_output_dir) &&
        !base::dir.exists(current_dir)) {
      base::dir.create(current_dir, recursive = TRUE, showWarnings = FALSE)
      message("Created missing output directory: ", current_dir)
    }
    paths[[nm]] <- fix_dir(current_dir)
  }

  hc@config@paths <- .hc_to_data_frame(paths)
  methods::validObject(hc)
  hc
}
#
#' @rdname check_dirs
#' @param hc A `HCoCenaExperiment`.
#' @param create_output_dir Boolean. If TRUE and `dir_output` is missing, create it.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_check_dirs <- function(hc, create_output_dir = TRUE) {
  .hc_check_dirs_impl(hc = hc, create_output_dir = create_output_dir)
}

#' Run `init_save_folder()` with S4 state synchronization
#'
# Internal implementation shared by S4 and legacy entry points.
.hc_init_save_folder_impl <- function(hc, name, use_output_dir = FALSE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (isTRUE(use_output_dir)) {
    name <- ""
  }
  if (base::is.null(name) || base::length(name) != 1 || base::is.na(name)) {
    stop("`name` must be a non-NA character scalar. Use \"\" to skip a subfolder.")
  }

  paths <- .hc_row_to_list(hc@config@paths)
  out_dir <- paths[["dir_output"]]
  if (base::is.null(out_dir) || base::identical(out_dir, FALSE) || !base::nzchar(base::as.character(out_dir))) {
    stop("`dir_output` is not set. Please run `init_wd()`/`hc_set_paths()` first.")
  }

  if (!base::dir.exists(out_dir)) {
    base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("Created output directory: ", out_dir)
  }

  save_folder <- base::as.character(name)
  if (base::identical(save_folder, "")) {
    message("Using output directory directly (no additional save subfolder): ", out_dir)
  } else {
    target_dir <- base::file.path(out_dir, save_folder)
    if (!base::dir.exists(target_dir)) {
      base::dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
      message("Created save folder: ", target_dir)
    } else {
      message("Using existing save folder: ", target_dir)
    }
  }

  global_cfg <- .hc_row_to_list(hc@config@global)
  global_cfg[["save_folder"]] <- name
  hc@config@global <- .hc_to_data_frame(global_cfg)
  methods::validObject(hc)
  hc
}
#
#' @rdname init_save_folder
#' @param hc A `HCoCenaExperiment`.
#' @param name Folder name. Use `""` to write directly into `dir_output`.
#' @param use_output_dir Boolean. If TRUE, ignore `name` and use `dir_output` directly.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_init_save_folder <- function(hc, name, use_output_dir = FALSE) {
  .hc_init_save_folder_impl(hc = hc, name = name, use_output_dir = use_output_dir)
}

#' Plot cut-off diagnostics (S4 API)
#'
#' @rdname plot_cutoffs
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `plot_cutoffs()`.
#' @template example-hc-after-part1
#' hc <- hc_plot_cutoffs(hc)
#' @return Updated `HCoCenaExperiment`.
# Internal implementation shared by S4 and legacy entry points.
.hc_plot_cutoffs_impl <- function(hc, ...) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = ".hc_plot_cutoffs_legacy_driver",
    ...
  )
  list(hc = out$hc, result = out$result)
}
#
#' @export
hc_plot_cutoffs <- function(hc, ...) {
  .hc_plot_cutoffs_impl(hc = hc, ...)[["hc"]]
}

#' Plot degree distributions (S4 API)
#'
#' @rdname plot_deg_dist
#' @param hc A `HCoCenaExperiment`.
#' @template example-hc-after-part1
#' hc <- hc_plot_deg_dist(hc)
#' @return Updated `HCoCenaExperiment`.
# Internal implementation shared by S4 and legacy entry points.
.hc_plot_deg_dist_impl <- function(hc) {
  out <- .hc_run_legacy_capture(
    hc = hc,
    fun = ".hc_plot_deg_dist_legacy_driver"
  )
  list(hc = out$hc, result = out$result)
}
#
#' @export
hc_plot_deg_dist <- function(hc) {
  .hc_plot_deg_dist_impl(hc = hc)[["hc"]]
}

#' Plot module heatmap (S4 API)
#'
#' @rdname plot_cluster_heatmap
#' @param hc A `HCoCenaExperiment`.
#' @param file_name Optional file name passed to `plot_cluster_heatmap()`.
#' @param ... Passed to `plot_cluster_heatmap()`.
#' @template example-hc-clustered
#' hc <- hc_plot_cluster_heatmap(hc, file_name = FALSE)
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

  base::do.call(.hc_plot_cluster_heatmap_legacy_driver, c(list(file_name = file_name), list(...)))
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
#' @rdname change_grouping_parameter
#' @param hc A `HCoCenaExperiment`.
#' @param group_by Grouping column to use (must be present in all annotation tables).
#' @param col_order Optional heatmap column order.
#' @param cluster_columns Whether to cluster heatmap columns. Default is
#'   `FALSE`, so the stored main hCoCena column order is reused when available.
#' @param row_order Optional heatmap row order.
#' @param cluster_rows Whether to cluster heatmap rows.
#' @template example-hc-clustered
#' hc <- hc_change_grouping_parameter(hc, group_by = "batch")
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
    fun = .hc_change_grouping_parameter_legacy_driver,
    group_by = group_by,
    col_order = col_order,
    cluster_columns = cluster_columns,
    row_order = row_order,
    cluster_rows = cluster_rows
  )
}

#' Plot integrated network (S4 API)
#'
#' @rdname plot_integrated_network
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `plot_integrated_network()`.
#' @template example-hc-clustered
#' hc <- hc_plot_integrated_network(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_integrated_network <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_plot_integrated_network_legacy_driver, ...)
}

#' Plot network colored by GFC (S4 API)
#'
#' @rdname plot_GFC_network
#' @param hc A `HCoCenaExperiment`.
#' @template example-hc-clustered
#' hc <- hc_plot_gfc_network(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_gfc_network <- function(hc) {
  .hc_run_legacy(hc = hc, fun = .hc_plot_GFC_network_legacy_driver)
}

#' TF enrichment per module (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `TF_overrep_module()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_tf_overrep_module <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_TF_overrep_module_legacy_driver, ...)
}

#' TF enrichment network-wide (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `TF_overrep_network()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_tf_overrep_network <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_TF_overrep_network_legacy_driver, ...)
}

#' Check TF targets (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param TF Transcription factor symbol.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_check_tf <- function(hc, TF) {
  .hc_run_legacy(hc = hc, fun = .hc_check_tf_legacy_driver, TF = TF)
}

#' Write session info (S4 API)
#'
#' @rdname write_session_info
#' @param hc A `HCoCenaExperiment`.
#' @template example-hc-prepared
#' hc <- hc_write_session_info(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_write_session_info <- function(hc) {
  .hc_run_legacy(hc = hc, fun = .hc_write_session_info_legacy_driver)
}

#' Suggest top variable genes (S4 API)
#'
#' @rdname suggest_topvar
#' @param hc A `HCoCenaExperiment`.
#' @template example-hc-prepared
#' hc <- hc_suggest_topvar(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_suggest_topvar <- function(hc) {
  .hc_run_legacy(hc = hc, fun = .hc_suggest_topvar_legacy_driver)
}

#' Plot sample distributions (S4 API)
#'
#' @rdname plot_sample_distributions
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `plot_sample_distributions()`.
#' @template example-hc-prepared
#' hc <- hc_plot_sample_distributions(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_plot_sample_distributions <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_plot_sample_distributions_legacy_driver, ...)
}

#' PCA plotting (S4 API)
#'
#' @rdname PCA
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `PCA()`.
#' @template example-hc-prepared
#' hc <- hc_pca(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_pca <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_PCA_legacy_driver, ...)
}

#' Meta-data plotting (S4 API)
#'
#' @rdname meta_plot
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `meta_plot()`.
#' @template example-hc-prepared
#' hc <- hc_meta_plot(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_meta_plot <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_meta_plot_legacy_driver, ...)
}

#' Export clusters (S4 API)
#'
#' @rdname export_clusters
#' @param hc A `HCoCenaExperiment`.
#' @template example-hc-clustered
#' hc <- hc_export_clusters(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_export_clusters <- function(hc) {
  .hc_run_legacy(hc = hc, fun = .hc_export_clusters_legacy_driver)
}

#' Module scores (S4 API)
#'
#' @rdname get_module_scores
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `get_module_scores()`.
#' @template example-hc-clustered
#' hc <- hc_get_module_scores(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_get_module_scores <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_get_module_scores_legacy_driver, ...)
}

#' Alluvial comparison plots (S4 API)
#'
#' @rdname algo_alluvial
#' @param hc A `HCoCenaExperiment`.
#' @template example-hc-clustered
#' hc <- hc_algo_alluvial(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_algo_alluvial <- function(hc) {
  .hc_run_legacy(hc = hc, fun = .hc_algo_alluvial_legacy_driver)
}

#' PCA algorithm comparison (S4 API)
#'
#' @rdname PCA_algo_compare
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `PCA_algo_compare()`.
#' @template example-hc-prepared
#' hc <- hc_pca_algo_compare(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_pca_algo_compare <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_PCA_algo_compare_legacy_driver, ...)
}

#' Update clustering algorithm (S4 API)
#'
#' @rdname update_clustering_algorithm
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `update_clustering_algorithm()`.
#' @template example-hc-clustered
#' hc <- hc_update_clustering_algorithm(hc, cluster_algo = "cluster_louvain")
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_update_clustering_algorithm <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_update_clustering_algorithm_legacy_driver, ...)
}

#' Export network to local folder (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `export_to_local_folder()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_export_to_local_folder <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_export_to_local_folder_legacy_driver, ...)
}

#' Import Cytoscape layout from local folder (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `import_layout_from_local_folder()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_import_layout_from_local_folder <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_import_layout_from_local_folder_legacy_driver, ...)
}

#' Export network to Cytoscape (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `export_to_cytoscape()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_export_to_cytoscape <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_export_to_cytoscape_legacy_driver, ...)
}

#' Import layout from Cytoscape (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_import_layout_from_cytoscape <- function(hc) {
  .hc_run_legacy(hc = hc, fun = .hc_import_layout_from_cytoscape_legacy_driver)
}

#' Hub detection (S4 API)
#'
#' @rdname find_hubs
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `find_hubs()`.
#' @template example-hc-clustered
#' hc <- hc_find_hubs(hc)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_find_hubs <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_find_hubs_legacy_driver, ...)
}

#' Visualize gene expression (S4 API)
#'
#' @rdname visualize_gene_expression
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `visualize_gene_expression()`.
#' @template example-hc-prepared
#' hc <- hc_visualize_gene_expression(hc, gene = "G1")
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_visualize_gene_expression <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_visualize_gene_expression_legacy_driver, ...)
}

#' Highlight gene set in network (S4 API)
#'
#' @rdname highlight_geneset
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `highlight_geneset()`.
#' @template example-hc-clustered
#' hc <- hc_highlight_geneset(hc, geneset = c("G1", "G2"))
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_highlight_geneset <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_highlight_geneset_legacy_driver, ...)
}

#' Highlight single cluster in network (S4 API)
#'
#' @rdname colour_single_cluster
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `colour_single_cluster()`.
#' @template example-hc-clustered
#' hc <- hc_colour_single_cluster(hc, selected_cluster = 1)
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_colour_single_cluster <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_colour_single_cluster_legacy_driver, ...)
}

#' Add categorical module heatmap annotations (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `col_anno_categorical()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_col_anno_categorical <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_col_anno_categorical_legacy_driver, ...)
}

#' Correlate categorical metadata with modules (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `meta_correlation_cat()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_meta_correlation_cat <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_meta_correlation_cat_legacy_driver, ...)
}

#' Test module differences between conditions (S4 API)
#'
#' @param hc A `HCoCenaExperiment`.
#' @param ... Passed to `module_condition_significance()`.
#' @return Updated `HCoCenaExperiment`.
#' @export
hc_module_condition_significance <- function(hc, ...) {
  .hc_run_legacy(hc = hc, fun = .hc_module_condition_significance_legacy_driver, ...)
}
