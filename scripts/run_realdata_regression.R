#!/usr/bin/env Rscript

# Real-data regression runner for local hCoCena development.
#
# The runner is intentionally opt-in: real data and generated outputs are kept
# outside package tests unless HCOCENA_RUN_REALDATA is enabled.

.hcr_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) {
    return(default)
  }
  tolower(trimws(as.character(x[[1]]))) %in% c("1", "true", "yes", "y", "on")
}

.hcr_scalar <- function(x, default = NULL) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(as.character(x[[1]]))) {
    return(default)
  }
  as.character(x[[1]])
}

.hcr_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }
    key <- sub("^--", "", arg)
    value <- "true"
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      key <- parts[[1]]
      value <- paste(parts[-1], collapse = "=")
    } else if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
      i <- i + 1L
      value <- args[[i]]
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- value
    i <- i + 1L
  }
  out
}

.hcr_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    desc <- file.path(path, "DESCRIPTION")
    if (file.exists(desc)) {
      dcf <- tryCatch(read.dcf(desc), error = function(e) NULL)
      if (!is.null(dcf) && identical(unname(dcf[1, "Package"]), "hcocena")) {
        return(path)
      }
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the hcocena repository root from: ", start, call. = FALSE)
    }
    path <- parent
  }
}

.hcr_slash <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (grepl("/$", path)) path else paste0(path, "/")
}

.hcr_default_data_dir <- function(repo_root) {
  env <- Sys.getenv("HCOCENA_REALDATA_DIR", unset = "")
  if (nzchar(env)) {
    return(env)
  }
  file.path(dirname(repo_root), "data")
}

.hcr_default_reference_dir <- function(repo_root, mode) {
  env <- Sys.getenv("HCOCENA_REALDATA_REFERENCE_DIR", unset = "")
  if (nzchar(env)) {
    return(file.path(env, mode))
  }
  file.path(repo_root, "realdata-reference", mode)
}

.hcr_load_package <- function(repo_root, load_source = TRUE) {
  if (isTRUE(load_source) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo_root, quiet = TRUE)
    return(invisible("source"))
  }
  suppressPackageStartupMessages(library(hcocena))
  invisible("installed")
}

.hcr_mode_config <- function(mode) {
  mode <- match.arg(mode, c("quick", "full"))
  if (identical(mode, "quick")) {
    return(list(
      mode = mode,
      top_var = c(2000, 2000),
      min_corr = c(0.90, 0.90),
      range_cutoff_length = c(5, 5),
      cutoffs = c(0.90, 0.90),
      min_nodes_number_for_network = 5,
      min_nodes_number_for_cluster = 5,
      cluster_algo = "cluster_fast_greedy",
      no_of_iterations = 1,
      resolution = 0.1,
      read_supplementary = TRUE,
      run_enrichment = TRUE,
      enrichment_gene_sets = c("Hallmark", "Kegg"),
      enrichment_top = 3,
      run_module_split = TRUE,
      split_cluster_algo = "cluster_fast_greedy",
      split_resolution = 0.1,
      split_resolution_grid = c(0.05, 0.1, 0.2),
      split_min_submodule_size = 5,
      plot_layer_heatmaps = FALSE,
      module_label_fontsize = 11,
      module_label_pt_size = 0.55,
      module_box_width_cm = 0.80
    ))
  }

  list(
    mode = mode,
    top_var = c("all", "all"),
    min_corr = c(0.90, 0.90),
    range_cutoff_length = c(100, 100),
    cutoffs = c(0.90, 0.90),
    min_nodes_number_for_network = 15,
    min_nodes_number_for_cluster = 15,
    cluster_algo = "cluster_leiden",
    no_of_iterations = 2,
    resolution = 0.1,
    read_supplementary = TRUE,
    run_enrichment = .hcr_bool(Sys.getenv("HCOCENA_REALDATA_FULL_ENRICHMENT", unset = "true"), TRUE),
    enrichment_gene_sets = c("Hallmark", "Kegg", "Go"),
    enrichment_top = 5,
    run_module_split = TRUE,
    split_cluster_algo = "cluster_fast_greedy",
    split_resolution = 0.1,
    split_resolution_grid = c(0.05, 0.1, 0.2),
    split_min_submodule_size = 15,
    plot_layer_heatmaps = FALSE,
    module_label_fontsize = 11,
    module_label_pt_size = 0.55,
    module_box_width_cm = 0.80
  )
}

.hcr_check_inputs <- function(data_dir, reference_dir) {
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  required <- c(
    "data_seq_processed.txt",
    "annotation_seq.txt",
    "data_array_processed.txt",
    "annotation_array.txt"
  )
  missing <- required[!file.exists(file.path(data_dir, required))]
  if (length(missing) > 0) {
    stop(
      "Real-data directory is missing required files: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(list(data_dir = data_dir, reference_dir = reference_dir))
}

.hcr_as_df <- function(x) {
  if (is.null(x)) {
    return(data.frame())
  }
  if (inherits(x, "DataFrame")) {
    return(as.data.frame(x))
  }
  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

.hcr_write_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = path, row.names = row.names, quote = TRUE, na = "")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.hcr_write_json <- function(x, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package `jsonlite` is required to write the regression manifest.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path = path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.hcr_sort_df <- function(x, preferred = character()) {
  x <- .hcr_as_df(x)
  if (nrow(x) == 0) {
    return(x)
  }
  cols <- unique(c(preferred[preferred %in% names(x)], names(x)))
  ord_data <- x[cols]
  ord <- do.call(order, lapply(ord_data, function(col) {
    if (is.numeric(col)) col else tolower(as.character(col))
  }))
  x[ord, , drop = FALSE]
}

.hcr_table_metrics <- function(hc) {
  layer_cfg <- .hcr_as_df(hc@config@layer)
  experiments <- MultiAssayExperiment::experiments(hc@mae)
  layer_rows <- lapply(seq_along(experiments), function(i) {
    se <- experiments[[i]]
    data.frame(
      metric = paste0("layer_", names(experiments)[[i]], "_dimensions"),
      value = paste0(nrow(se), " genes x ", ncol(se), " samples"),
      stringsAsFactors = FALSE
    )
  })
  cfg_rows <- if (nrow(layer_cfg) > 0) {
    data.frame(
      metric = paste0("configured_layer_", seq_len(nrow(layer_cfg))),
      value = apply(layer_cfg, 1, function(row) paste(names(row), row, sep = "=", collapse = "; ")),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }
  do.call(rbind, c(layer_rows, list(cfg_rows)))
}

.hcr_parse_module_genes <- function(gene_n) {
  parser <- tryCatch(
    get(".hc_parse_genes_from_gene_n", envir = asNamespace("hcocena"), inherits = FALSE),
    error = function(e) NULL
  )
  if (is.function(parser)) {
    out <- tryCatch(parser(gene_n), error = function(e) character())
    return(unique(as.character(out[!is.na(out) & nzchar(out)])))
  }

  gene_n <- as.character(gene_n)
  gene_n <- gene_n[!is.na(gene_n)]
  if (length(gene_n) == 0) {
    return(character())
  }
  out <- unlist(strsplit(gene_n, "[#,;|[:space:]]+"), use.names = FALSE)
  unique(out[!is.na(out) & nzchar(out)])
}

.hcr_choose_split_module <- function(hc) {
  cluster <- as.list(hc@integration@cluster)
  cluster_info <- .hcr_as_df(cluster[["cluster_information"]])
  label_map <- cluster[["module_label_map"]]
  if (nrow(cluster_info) == 0 || !all(c("color", "gene_n") %in% names(cluster_info))) {
    stop("Cannot choose a split module: cluster information is missing `color` or `gene_n`.", call. = FALSE)
  }

  included <- if ("cluster_included" %in% names(cluster_info)) {
    as.character(cluster_info$cluster_included) == "yes"
  } else {
    rep(TRUE, nrow(cluster_info))
  }
  candidates <- cluster_info[included, , drop = FALSE]
  if (nrow(candidates) == 0) {
    stop("Cannot choose a split module: no included modules found.", call. = FALSE)
  }

  gene_counts <- vapply(candidates$gene_n, function(x) length(.hcr_parse_module_genes(x)), integer(1))
  target_idx <- which.max(gene_counts)
  target_color <- as.character(candidates$color[[target_idx]])
  target_label <- if (length(label_map) > 0 && target_color %in% names(label_map)) {
    as.character(label_map[[target_color]])
  } else {
    target_color
  }

  list(
    module_color = target_color,
    module_label = target_label,
    gene_count = as.integer(gene_counts[[target_idx]])
  )
}

.hcr_empty_df <- function(cols) {
  stats::setNames(as.data.frame(rep(list(character()), length(cols))), cols)
}

.hcr_enrichment_outputs <- function(satellite) {
  enrich <- satellite[["enrichments"]]
  if (is.null(enrich) || !is.list(enrich)) {
    enrich <- list()
  }
  all <- .hcr_sort_df(enrich[["all_enrichments_all_dbs"]], c("database", "module", "module_label", "ID", "Description", "p.adjust"))
  selected <- .hcr_sort_df(enrich[["selected_enrichments_all_dbs"]], c("database", "module", "module_label", "ID", "Description", "p.adjust"))
  significant <- .hcr_sort_df(enrich[["significant_enrichments_all_dbs"]], c("database", "module", "module_label", "ID", "Description", "p.adjust"))
  if (ncol(all) == 0) all <- .hcr_empty_df(c("database", "module", "module_label", "ID", "Description", "p.adjust"))
  if (ncol(selected) == 0) selected <- .hcr_empty_df(c("database", "module", "module_label", "ID", "Description", "p.adjust"))
  if (ncol(significant) == 0) significant <- .hcr_empty_df(c("database", "module", "module_label", "ID", "Description", "p.adjust"))
  list(all = all, selected = selected, significant = significant)
}

.hcr_split_outputs <- function(satellite) {
  resolution <- satellite[["module_split_resolution_test_last"]]
  last <- satellite[["module_split_last"]]

  by_resolution <- if (!is.null(resolution) && is.list(resolution)) {
    .hcr_sort_df(resolution[["by_resolution"]], c("module_label", "resolution"))
  } else {
    data.frame()
  }
  by_module <- if (!is.null(resolution) && is.list(resolution)) {
    .hcr_sort_df(resolution[["by_module"]], c("module_label", "resolution"))
  } else {
    data.frame()
  }
  resolved <- if (!is.null(last) && is.list(last)) {
    .hcr_sort_df(last[["resolved_modules"]], c("resolved_index", "resolved_label", "resolved_color"))
  } else {
    data.frame()
  }
  summary <- if (!is.null(last) && is.list(last)) {
    .hcr_sort_df(last[["split_summary"]], c("parent_label", "status"))
  } else {
    data.frame()
  }

  if (ncol(by_resolution) == 0) by_resolution <- .hcr_empty_df(c("module_label", "resolution", "n_submodules"))
  if (ncol(by_module) == 0) by_module <- .hcr_empty_df(c("module_label", "resolution", "n_submodules"))
  if (ncol(resolved) == 0) resolved <- .hcr_empty_df(c("input", "resolved_index", "resolved_color", "resolved_label", "status"))
  if (ncol(summary) == 0) summary <- .hcr_empty_df(c("parent_color", "parent_label", "parent_genes", "n_submodules", "status"))

  list(
    by_resolution = by_resolution,
    by_module = by_module,
    resolved = resolved,
    summary = summary
  )
}

.hcr_visual_text_page <- function(title, lines = character(), footer = NULL) {
  grid::grid.newpage()
  grid::grid.text(
    title,
    x = grid::unit(0.05, "npc"),
    y = grid::unit(0.94, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )
  lines <- unlist(lapply(as.character(lines), strwrap, width = 120), use.names = FALSE)
  if (length(lines) > 0) {
    y <- 0.86
    for (line in lines[seq_len(min(length(lines), 36L))]) {
      grid::grid.text(
        line,
        x = grid::unit(0.05, "npc"),
        y = grid::unit(y, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 9)
      )
      y <- y - 0.023
    }
  }
  if (!is.null(footer)) {
    grid::grid.text(
      as.character(footer),
      x = grid::unit(0.05, "npc"),
      y = grid::unit(0.04, "npc"),
      just = c("left", "bottom"),
      gp = grid::gpar(fontsize = 8, col = "#555555")
    )
  }
  invisible(NULL)
}

.hcr_visual_table_page <- function(title, df, n = 18L, footer = NULL) {
  df <- .hcr_as_df(df)
  if (nrow(df) == 0) {
    return(.hcr_visual_text_page(title, "No rows available.", footer = footer))
  }
  lines <- utils::capture.output(print(utils::head(df, n), row.names = FALSE, right = FALSE))
  .hcr_visual_text_page(title, lines, footer = footer)
}

.hcr_visual_parameters_df <- function(parameters) {
  if (is.null(parameters) || length(parameters) == 0) {
    return(data.frame())
  }
  values <- vapply(parameters, .hcr_value_label, FUN.VALUE = character(1), USE.NAMES = FALSE)
  data.frame(
    parameter = names(parameters),
    value = values,
    stringsAsFactors = FALSE
  )
}

.hcr_value_label <- function(value) {
  if (is.null(value)) {
    return("NULL")
  }
  if (is.atomic(value) && length(value) <= 12) {
    return(paste(as.character(value), collapse = ", "))
  }
  paste(utils::capture.output(str(value, give.attr = FALSE)), collapse = " ")
}

.hcr_values_equal <- function(value, default) {
  if (is.null(value) && is.null(default)) {
    return(TRUE)
  }
  if (is.null(value) || is.null(default)) {
    return(FALSE)
  }
  if (is.numeric(value) && is.numeric(default) && length(value) == length(default)) {
    return(isTRUE(all.equal(as.numeric(value), as.numeric(default), tolerance = 1e-12)))
  }
  identical(as.character(value), as.character(default))
}

.hcr_parameter_diffs_df <- function(parameters, defaults = list(), scope = "") {
  if (is.null(parameters) || length(parameters) == 0) {
    return(data.frame())
  }
  rows <- lapply(names(parameters), function(parameter) {
    value <- parameters[[parameter]]
    default <- defaults[[parameter]]
    if (.hcr_values_equal(value, default)) {
      return(NULL)
    }
    data.frame(
      scope = scope,
      parameter = parameter,
      value = .hcr_value_label(value),
      default = .hcr_value_label(default),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

.hcr_heatmap_plot_defaults <- function() {
  list(
    module_label_preset = "balanced",
    module_label_fontsize = NULL,
    module_label_pt_size = NULL,
    module_box_width_cm = NULL,
    gene_count_mode = "text"
  )
}

.hcr_enrichment_plot_defaults <- function() {
  list(
    gene_sets = "Hallmark",
    top = 5,
    consistent_terms = TRUE,
    cluster_columns = FALSE,
    store_panel_objects = "auto",
    heatmap_module_label_fontsize = NULL
  )
}

.hcr_plot_catalog_df <- function(plot_variants) {
  if (is.null(plot_variants) || length(plot_variants) == 0) {
    return(data.frame())
  }
  rows <- lapply(plot_variants, function(variant) {
    data.frame(
      plot = .hcr_scalar(variant$plot, ""),
      variant = .hcr_scalar(variant$variant, ""),
      pdf_file = .hcr_scalar(variant$pdf_file, ""),
      png_file = .hcr_scalar(variant$png_file, ""),
      stage = .hcr_scalar(variant$stage, ""),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.hcr_plot_parameter_diffs_df <- function(plot_variants, font_probe = data.frame()) {
  if (is.null(plot_variants) || length(plot_variants) == 0) {
    return(data.frame())
  }
  rows <- list()
  for (variant in plot_variants) {
    variant_name <- paste(
      c(.hcr_scalar(variant$plot, ""), .hcr_scalar(variant$variant, "")),
      collapse = " | "
    )
    params <- variant$parameters
    defaults <- variant$defaults
    diff <- .hcr_parameter_diffs_df(params, defaults, scope = variant_name)
    if (nrow(diff) > 0) {
      diff$source <- "non-default parameter"
      rows[[length(rows) + 1L]] <- diff
    }
  }

  font_probe_scope <- "Cluster heatmap | module-label font-size probe"
  font_probe_idx <- which(vapply(plot_variants, function(variant) {
    identical(.hcr_scalar(variant$plot, ""), "Cluster heatmap") &&
      grepl("font-size probe", .hcr_scalar(variant$variant, ""), ignore.case = TRUE)
  }, FUN.VALUE = logical(1)))
  if (length(font_probe_idx) > 0) {
    font_probe_variant <- plot_variants[[font_probe_idx[[1]]]]
    font_probe_scope <- paste(
      c(.hcr_scalar(font_probe_variant$plot, ""), .hcr_scalar(font_probe_variant$variant, "")),
      collapse = " | "
    )
  }

  if (nrow(font_probe) > 0 &&
    "module_label_auto_fit_shrunk" %in% names(font_probe) &&
    isTRUE(tolower(as.character(font_probe$module_label_auto_fit_shrunk[[1]])) %in% c("true", "1", "yes"))) {
    reason <- paste(
      c(
        if ("module_label_auto_fit_width_limited" %in% names(font_probe) &&
          tolower(as.character(font_probe$module_label_auto_fit_width_limited[[1]])) %in% c("true", "1", "yes")) "width",
        if ("module_label_auto_fit_height_limited" %in% names(font_probe) &&
          tolower(as.character(font_probe$module_label_auto_fit_height_limited[[1]])) %in% c("true", "1", "yes")) "height"
      ),
      collapse = "+"
    )
    if (!nzchar(reason)) {
      reason <- "fit guard"
    }
    rows[[length(rows) + 1L]] <- data.frame(
      scope = font_probe_scope,
      parameter = "effective_module_label_font_size_pt",
      value = paste0(
        .hcr_scalar(font_probe$effective_module_label_pt_size_min_pt, ""),
        " to ",
        .hcr_scalar(font_probe$effective_module_label_pt_size_max_pt, "")
      ),
      default = paste0("requested max ", .hcr_scalar(font_probe$requested_module_label_pt_size_max_pt, "")),
      source = paste0("runtime auto-fit (", reason, ")"),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

.hcr_run_parameter_defaults <- function(mode) {
  c(
    list(
      mode = mode,
      seed = 168575,
      compare_reference = TRUE,
      update_reference = FALSE,
      tolerance = 1e-8
    ),
    .hcr_mode_config(mode)
  )
}

.hcr_plot_file <- function(path, base_dir) {
  rel <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (startsWith(rel, prefix)) {
    rel <- substring(rel, nchar(prefix) + 1L)
  }
  rel
}

.hcr_plot_variant <- function(plot,
                              variant,
                              pdf_file,
                              base_dir,
                              stage,
                              parameters,
                              defaults) {
  png_file <- sub("\\.pdf$", ".png", pdf_file, ignore.case = TRUE)
  list(
    plot = plot,
    variant = variant,
    pdf_file = .hcr_plot_file(pdf_file, base_dir),
    png_file = if (file.exists(png_file)) .hcr_plot_file(png_file, base_dir) else "",
    stage = stage,
    parameters = parameters,
    defaults = defaults
  )
}

.hcr_run_parameter_diffs_df <- function(parameters, mode) {
  keep <- setdiff(names(parameters), c("data_dir", "output_dir", "reference_dir", "package_source", "package_version"))
  .hcr_parameter_diffs_df(parameters[keep], .hcr_run_parameter_defaults(mode), scope = "runner")
}

.hcr_visual_parameter_table <- function(title, df, empty_message = "No non-default parameters.") {
  df <- .hcr_as_df(df)
  if (nrow(df) == 0) {
    return(.hcr_visual_text_page(title, empty_message))
  }
  .hcr_visual_table_page(title, df, n = 48L)
}

.hcr_is_abs_path <- function(path) {
  path <- as.character(path)
  grepl("^[A-Za-z]:[/\\\\]", path) || grepl("^[/\\\\]{2}", path) || startsWith(path, "/")
}

.hcr_abs_file <- function(path, base_dir) {
  path <- .hcr_scalar(path, "")
  if (!nzchar(path)) {
    return("")
  }
  if (.hcr_is_abs_path(path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, path), winslash = "/", mustWork = FALSE)
}

.hcr_plot_png_file <- function(variant, output_dir) {
  png_file <- .hcr_scalar(variant$png_file, "")
  if (!nzchar(png_file)) {
    pdf_file <- .hcr_scalar(variant$pdf_file, "")
    if (nzchar(pdf_file)) {
      png_file <- sub("\\.pdf$", ".png", pdf_file, ignore.case = TRUE)
    }
  }
  .hcr_abs_file(png_file, output_dir)
}

.hcr_plot_detail_lines <- function(variant, output_dir, font_probe = data.frame()) {
  pdf_file <- .hcr_scalar(variant$pdf_file, "")
  png_file <- .hcr_scalar(variant$png_file, "")
  if (!nzchar(png_file) && nzchar(pdf_file)) {
    png_file <- sub("\\.pdf$", ".png", pdf_file, ignore.case = TRUE)
  }
  diff <- .hcr_parameter_diffs_df(variant$parameters, variant$defaults)
  parameter_line <- if (nrow(diff) == 0) {
    "Non-default parameters: none"
  } else {
    paste0(
      "Non-default parameters: ",
      paste0(diff$parameter, "=", diff$value, " (default ", diff$default, ")", collapse = "; ")
    )
  }
  lines <- c(
    paste0("Stage: ", .hcr_scalar(variant$stage, "")),
    paste0("PDF: ", pdf_file),
    paste0("PNG: ", png_file),
    parameter_line
  )

  is_font_probe <- grepl("font-size probe", .hcr_scalar(variant$variant, ""), ignore.case = TRUE)
  if (is_font_probe && nrow(font_probe) > 0) {
    lines <- c(
      lines,
      paste0(
        "Runtime module label size: ",
        .hcr_scalar(font_probe$effective_module_label_pt_size_min_pt, ""),
        " to ",
        .hcr_scalar(font_probe$effective_module_label_pt_size_max_pt, ""),
        " pt (requested max ",
        .hcr_scalar(font_probe$requested_module_label_pt_size_max_pt, ""),
        " pt; auto-fit shrunk=",
        .hcr_scalar(font_probe$module_label_auto_fit_shrunk, ""),
        ")"
      )
    )
  }
  lines
}

.hcr_visual_plot_page <- function(title, png_file, lines = character(), footer = NULL) {
  if (!nzchar(.hcr_scalar(png_file, "")) || !file.exists(png_file)) {
    return(.hcr_visual_text_page(
      title,
      c("Plot image is unavailable.", paste0("PNG: ", png_file), lines),
      footer = footer
    ))
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    return(.hcr_visual_text_page(
      title,
      c("Package `png` is not installed, so the plot image cannot be embedded.", paste0("PNG: ", png_file), lines),
      footer = footer
    ))
  }

  image <- tryCatch(
    png::readPNG(png_file),
    error = function(e) e
  )
  if (inherits(image, "error")) {
    return(.hcr_visual_text_page(
      title,
      c("Plot image could not be read.", paste0("PNG: ", png_file), paste0("Error: ", conditionMessage(image)), lines),
      footer = footer
    ))
  }

  grid::grid.newpage()
  grid::grid.text(
    title,
    x = grid::unit(0.05, "npc"),
    y = grid::unit(0.985, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 11, fontface = "bold")
  )

  detail_lines <- unlist(lapply(as.character(lines), strwrap, width = 150), use.names = FALSE)
  detail_lines <- detail_lines[seq_len(min(length(detail_lines), 6L))]

  detail_bottom <- if (is.null(footer)) 0.035 else 0.055
  detail_line_height <- 0.017
  detail_height <- if (length(detail_lines) > 0) {
    (length(detail_lines) * detail_line_height) + 0.014
  } else {
    0
  }
  image_top <- 0.955
  image_bottom <- detail_bottom + detail_height
  max_w <- 0.96
  max_h <- max(0.25, image_top - image_bottom)
  aspect <- dim(image)[[2]] / dim(image)[[1]]
  if (aspect >= max_w / max_h) {
    width <- max_w
    height <- max_w / aspect
  } else {
    height <- max_h
    width <- max_h * aspect
  }
  grid::grid.raster(
    image,
    x = grid::unit(0.5, "npc"),
    y = grid::unit(image_bottom + (max_h / 2), "npc"),
    width = grid::unit(width, "npc"),
    height = grid::unit(height, "npc"),
    interpolate = TRUE
  )

  if (length(detail_lines) > 0) {
    y <- detail_bottom + (length(detail_lines) * detail_line_height)
    for (line in detail_lines) {
      grid::grid.text(
        line,
        x = grid::unit(0.05, "npc"),
        y = grid::unit(y, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 6.4, col = "#333333")
      )
      y <- y - detail_line_height
    }
  }

  if (!is.null(footer)) {
    grid::grid.text(
      as.character(footer),
      x = grid::unit(0.05, "npc"),
      y = grid::unit(0.028, "npc"),
      just = c("left", "bottom"),
      gp = grid::gpar(fontsize = 6.8, col = "#555555")
    )
  }
  invisible(NULL)
}

.hcr_visual_plot_pages <- function(plot_variants, output_dir, font_probe = data.frame()) {
  if (is.null(plot_variants) || length(plot_variants) == 0) {
    return(invisible(NULL))
  }
  footer <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  for (variant in plot_variants) {
    title <- paste(
      c("Plot image", .hcr_scalar(variant$plot, ""), .hcr_scalar(variant$variant, "")),
      collapse = " - "
    )
    .hcr_visual_plot_page(
      title = title,
      png_file = .hcr_plot_png_file(variant, output_dir),
      lines = .hcr_plot_detail_lines(variant, output_dir, font_probe = font_probe),
      footer = footer
    )
  }
  invisible(NULL)
}

.hcr_open_visual_report_pdf <- function(report_file, width = 13, height = 9.5) {
  candidates <- report_file
  fallback <- file.path(
    dirname(report_file),
    paste0(
      tools::file_path_sans_ext(basename(report_file)),
      "_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".pdf"
    )
  )
  candidates <- unique(c(candidates, fallback))
  errors <- character()
  for (candidate in candidates) {
    ok <- tryCatch({
      grDevices::pdf(candidate, width = width, height = height, onefile = TRUE)
      TRUE
    }, error = function(e) {
      errors <<- c(errors, paste0(candidate, ": ", conditionMessage(e)))
      FALSE
    })
    if (isTRUE(ok)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  stop(
    "Could not open a visual-check PDF device. Tried:\n",
    paste(errors, collapse = "\n"),
    call. = FALSE
  )
}

.hcr_visual_parameters_df_legacy <- function(parameters) {
  if (is.null(parameters) || length(parameters) == 0) {
    return(data.frame())
  }
  values <- vapply(parameters, function(value) {
    if (is.null(value)) {
      return("NULL")
    }
    if (is.atomic(value) && length(value) <= 12) {
      return(paste(as.character(value), collapse = ", "))
    }
    paste(utils::capture.output(str(value, give.attr = FALSE)), collapse = " ")
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
  data.frame(
    parameter = names(parameters),
    value = values,
    stringsAsFactors = FALSE
  )
}

.hcr_write_visual_report <- function(hc,
                                     output_dir,
                                     mode,
                                     font_probe_file,
                                     split_target = NULL,
                                     parameters = NULL,
                                     plot_variants = NULL) {
  report_file <- file.path(output_dir, paste0("visual_check_report_", mode, ".pdf"))
  dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)

  read_export <- function(name) {
    path <- file.path(output_dir, name)
    if (!file.exists(path)) {
      return(data.frame())
    }
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }

  qc <- read_export("qc_summary.csv")
  split_summary <- read_export("split_summary.csv")
  split_resolution <- read_export("split_resolution_summary.csv")
  enrichment_selected <- read_export("enrichment_selected.csv")
  font_probe <- read_export("font_probe.csv")
  plot_catalog <- .hcr_plot_catalog_df(plot_variants)
  plot_parameter_diffs <- .hcr_plot_parameter_diffs_df(plot_variants, font_probe = font_probe)
  run_parameter_diffs <- .hcr_run_parameter_diffs_df(parameters, mode = mode)

  report_file <- .hcr_open_visual_report_pdf(report_file, width = 13, height = 9.5)
  on.exit(grDevices::dev.off(), add = TRUE)

  qc_lines <- c(
    paste0("Mode: ", mode),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    if (!is.null(split_target)) {
      paste0(
        "Split target: ", split_target$module_label,
        " (", split_target$gene_count, " genes)"
      )
    } else {
      "Split target: none"
    },
    "",
    utils::capture.output(print(qc, row.names = FALSE, right = FALSE))
  )
  .hcr_visual_text_page(
    title = "hCoCena real-data visual check",
    lines = qc_lines,
    footer = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  )

  .hcr_visual_table_page("Plot variants", plot_catalog, n = 32L)

  .hcr_visual_parameter_table(
    "Plot parameter deviations",
    plot_parameter_diffs,
    empty_message = "No plot parameters differ from their defaults."
  )

  .hcr_visual_parameter_table(
    "Run parameter deviations",
    run_parameter_diffs,
    empty_message = "No runner parameters differ from the mode defaults."
  )

  .hcr_visual_plot_pages(plot_variants, output_dir, font_probe = font_probe)

  .hcr_visual_table_page("Module split summary", split_summary)

  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(split_resolution) > 0 &&
    all(c("resolution", "total_created_submodules", "resulting_included_modules") %in% names(split_resolution))) {
    p <- ggplot2::ggplot(split_resolution, ggplot2::aes(x = resolution)) +
      ggplot2::geom_line(ggplot2::aes(y = total_created_submodules, color = "created submodules"), linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(y = total_created_submodules, color = "created submodules"), size = 2) +
      ggplot2::geom_line(ggplot2::aes(y = resulting_included_modules, color = "included modules"), linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(y = resulting_included_modules, color = "included modules"), size = 2) +
      ggplot2::labs(
        title = "Split resolution probe",
        x = "Resolution",
        y = "Count",
        color = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12)
    print(p)
  } else {
    .hcr_visual_table_page("Split resolution probe", split_resolution)
  }

  if (requireNamespace("ggplot2", quietly = TRUE) && nrow(enrichment_selected) > 0 &&
    all(c("database", "module_label", "term", "p.adjust") %in% names(enrichment_selected))) {
    enrich_plot <- enrichment_selected
    enrich_plot$p_adjust_num <- suppressWarnings(as.numeric(enrich_plot$p.adjust))
    enrich_plot <- enrich_plot[is.finite(enrich_plot$p_adjust_num), , drop = FALSE]
    enrich_plot <- enrich_plot[order(enrich_plot$p_adjust_num), , drop = FALSE]
    enrich_plot <- utils::head(enrich_plot, 18L)
    enrich_plot$label <- paste0(enrich_plot$module_label, " | ", enrich_plot$database, " | ", enrich_plot$term)
    enrich_plot$label <- factor(enrich_plot$label, levels = rev(enrich_plot$label))
    p <- ggplot2::ggplot(enrich_plot, ggplot2::aes(x = -log10(p_adjust_num), y = label)) +
      ggplot2::geom_col(fill = "#4C78A8") +
      ggplot2::labs(
        title = "Selected enrichment terms",
        x = "-log10(adjusted p)",
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
    print(p)
  } else {
    .hcr_visual_table_page("Selected enrichment terms", enrichment_selected)
  }

  .hcr_visual_table_page("Module label font-size probe", font_probe, n = 8L)

  cluster <- as.list(hc@integration@cluster)
  heatmap_obj <- if (!is.null(cluster[["heatmap_cluster"]])) {
    cluster[["heatmap_cluster"]]
  } else {
    cluster[["heatmap_cluster_raw"]]
  }
  if (!is.null(heatmap_obj) && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    ok <- tryCatch({
      ComplexHeatmap::draw(heatmap_obj, newpage = TRUE, merge_legends = TRUE)
      TRUE
    }, error = function(e) {
      .hcr_visual_text_page(
        "Cluster heatmap",
        c(
          "The heatmap object could not be redrawn in the visual report.",
          paste0("Standalone heatmap file: ", file.path("hcocena", basename(font_probe_file))),
          paste0("Error: ", conditionMessage(e))
        )
      )
      FALSE
    })
    invisible(ok)
  } else {
    .hcr_visual_text_page(
      "Cluster heatmap",
      c(
        "No drawable heatmap object was stored.",
        paste0("Standalone heatmap file: ", file.path("hcocena", basename(font_probe_file)))
      )
    )
  }

  normalizePath(report_file, winslash = "/", mustWork = TRUE)
}

.hcr_export_outputs <- function(hc, output_dir, mode, font_probe_file, plot_variants = NULL) {
  cluster <- as.list(hc@integration@cluster)
  satellite <- as.list(hc@satellite)
  enrichment <- .hcr_enrichment_outputs(satellite)
  split <- .hcr_split_outputs(satellite)
  cluster_info <- .hcr_sort_df(cluster[["cluster_information"]], c("cluster", "color", "gene_n"))
  edge_list <- .hcr_sort_df(hc@integration@combined_edgelist, c("V1", "V2", "weight"))
  gfc <- .hcr_sort_df(hc@integration@gfc, c("Gene"))
  module_gene_src <- if (!is.null(satellite[["module_gene_list"]])) {
    satellite[["module_gene_list"]]
  } else {
    cluster[["module_gene_list"]]
  }
  module_gene_list <- .hcr_sort_df(module_gene_src, c("module", "module_color", "gene"))
  if (ncol(module_gene_list) == 0) {
    module_gene_list <- data.frame(module = character(), module_color = character(), gene = character())
  }

  label_map <- cluster[["module_label_map"]]
  label_map_df <- if (length(label_map) > 0) {
    data.frame(
      module_color = names(label_map),
      module_label = as.character(label_map),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(module_color = character(), module_label = character())
  }
  label_map_df <- .hcr_sort_df(label_map_df, c("module_label", "module_color"))

  heatmap_matrix <- cluster[["heatmap_matrix"]]
  heatmap_df <- if (!is.null(heatmap_matrix) && length(heatmap_matrix) > 0) {
    data.frame(module = rownames(heatmap_matrix), as.data.frame(heatmap_matrix, check.names = FALSE), check.names = FALSE)
  } else {
    data.frame()
  }

  graph <- hc@integration@graph
  graph_metrics <- data.frame(
    metric = c("nodes", "edges"),
    value = c(
      if (inherits(graph, "igraph")) igraph::vcount(graph) else NA_integer_,
      if (inherits(graph, "igraph")) igraph::ecount(graph) else NA_integer_
    ),
    stringsAsFactors = FALSE
  )
  module_metrics <- data.frame(
    metric = c(
      "included_modules",
      "cluster_rows",
      "module_label_fontsize",
      "module_label_pt_size",
      "module_label_pt_size_effective_pt_min",
      "module_label_pt_size_effective_pt_max",
      "module_label_pt_size_requested_pt_max",
      "module_label_auto_fit_shrunk",
      "module_label_auto_fit_width_limited",
      "module_label_auto_fit_height_limited",
      "module_box_width_cm",
      "split_summary_rows",
      "enrichment_all_rows",
      "enrichment_selected_rows",
      "enrichment_significant_rows",
      "heatmap_rows",
      "heatmap_columns"
    ),
    value = c(
      if (nrow(cluster_info) > 0 && "cluster_included" %in% names(cluster_info)) {
        length(unique(cluster_info$color[cluster_info$cluster_included == "yes"]))
      } else {
        NA_integer_
      },
      nrow(cluster_info),
      .hcr_scalar(cluster[["module_label_fontsize"]], NA_character_),
      .hcr_scalar(cluster[["module_label_pt_size"]], NA_character_),
      .hcr_scalar(cluster[["module_label_pt_size_effective_pt_min"]], NA_character_),
      .hcr_scalar(cluster[["module_label_pt_size_effective_pt_max"]], NA_character_),
      .hcr_scalar(cluster[["module_label_pt_size_requested_pt_max"]], NA_character_),
      .hcr_scalar(cluster[["module_label_auto_fit_shrunk"]], NA_character_),
      .hcr_scalar(cluster[["module_label_auto_fit_width_limited"]], NA_character_),
      .hcr_scalar(cluster[["module_label_auto_fit_height_limited"]], NA_character_),
      .hcr_scalar(cluster[["module_box_width_cm"]], NA_character_),
      nrow(split$summary),
      nrow(enrichment$all),
      nrow(enrichment$selected),
      nrow(enrichment$significant),
      nrow(heatmap_df),
      max(0L, ncol(heatmap_df) - 1L)
    ),
    stringsAsFactors = FALSE
  )
  font_probe <- data.frame(
    requested_module_label_fontsize = .hcr_scalar(cluster[["module_label_fontsize"]], NA_character_),
    requested_module_label_pt_size = .hcr_scalar(cluster[["module_label_pt_size"]], NA_character_),
    requested_module_box_width_cm = .hcr_scalar(cluster[["module_box_width_cm"]], NA_character_),
    requested_module_label_pt_size_max_pt = .hcr_scalar(cluster[["module_label_pt_size_requested_pt_max"]], NA_character_),
    effective_module_label_pt_size_min_pt = .hcr_scalar(cluster[["module_label_pt_size_effective_pt_min"]], NA_character_),
    effective_module_label_pt_size_max_pt = .hcr_scalar(cluster[["module_label_pt_size_effective_pt_max"]], NA_character_),
    module_label_auto_fit_shrunk = .hcr_scalar(cluster[["module_label_auto_fit_shrunk"]], NA_character_),
    module_label_auto_fit_width_limited = .hcr_scalar(cluster[["module_label_auto_fit_width_limited"]], NA_character_),
    module_label_auto_fit_height_limited = .hcr_scalar(cluster[["module_label_auto_fit_height_limited"]], NA_character_),
    effective_size_control = "module_label_pt_size_with_box_fit_guard",
    heatmap_file = file.path("hcocena", basename(font_probe_file)),
    heatmap_file_exists = file.exists(font_probe_file),
    stringsAsFactors = FALSE
  )
  plot_catalog <- .hcr_plot_catalog_df(plot_variants)
  plot_parameter_diffs <- .hcr_plot_parameter_diffs_df(plot_variants, font_probe = font_probe)
  qc <- rbind(.hcr_table_metrics(hc), graph_metrics, module_metrics)

  files <- list(
    qc_summary = .hcr_write_csv(qc, file.path(output_dir, "qc_summary.csv")),
    combined_edgelist = .hcr_write_csv(edge_list, file.path(output_dir, "combined_edgelist.csv")),
    gfc_matrix = .hcr_write_csv(gfc, file.path(output_dir, "gfc_matrix.csv")),
    cluster_information = .hcr_write_csv(cluster_info, file.path(output_dir, "cluster_information.csv")),
    module_label_map = .hcr_write_csv(label_map_df, file.path(output_dir, "module_label_map.csv")),
    module_gene_list = .hcr_write_csv(module_gene_list, file.path(output_dir, "module_gene_list.csv")),
    heatmap_matrix = .hcr_write_csv(heatmap_df, file.path(output_dir, "heatmap_matrix.csv")),
    font_probe = .hcr_write_csv(font_probe, file.path(output_dir, "font_probe.csv")),
    plot_variants = .hcr_write_csv(plot_catalog, file.path(output_dir, "plot_variants.csv")),
    plot_variant_parameters = .hcr_write_csv(plot_parameter_diffs, file.path(output_dir, "plot_variant_parameters.csv")),
    split_resolution_summary = .hcr_write_csv(split$by_resolution, file.path(output_dir, "split_resolution_summary.csv")),
    split_resolution_by_module = .hcr_write_csv(split$by_module, file.path(output_dir, "split_resolution_by_module.csv")),
    split_resolved_modules = .hcr_write_csv(split$resolved, file.path(output_dir, "split_resolved_modules.csv")),
    split_summary = .hcr_write_csv(split$summary, file.path(output_dir, "split_summary.csv")),
    enrichment_all = .hcr_write_csv(enrichment$all, file.path(output_dir, "enrichment_all.csv")),
    enrichment_selected = .hcr_write_csv(enrichment$selected, file.path(output_dir, "enrichment_selected.csv")),
    enrichment_significant = .hcr_write_csv(enrichment$significant, file.path(output_dir, "enrichment_significant.csv"))
  )

  list(
    files = files,
    metrics = qc,
    font_probe = font_probe,
    plot_variants = plot_catalog,
    plot_variant_parameters = plot_parameter_diffs
  )
}

.hcr_compare_table <- function(current_path, reference_path, tolerance = 1e-8) {
  if (!file.exists(reference_path)) {
    return(list(status = "missing_reference", details = "reference file not found"))
  }
  cur <- utils::read.csv(current_path, check.names = FALSE, stringsAsFactors = FALSE)
  ref <- utils::read.csv(reference_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(dim(cur), dim(ref))) {
    return(list(status = "different", details = paste0("dimension ", paste(dim(cur), collapse = "x"), " != ", paste(dim(ref), collapse = "x"))))
  }
  if (!identical(names(cur), names(ref))) {
    return(list(status = "different", details = "column names differ"))
  }
  for (nm in names(cur)) {
    cx <- cur[[nm]]
    rx <- ref[[nm]]
    c_num <- suppressWarnings(as.numeric(cx))
    r_num <- suppressWarnings(as.numeric(rx))
    numeric_like <- !all(is.na(c_num)) || !all(is.na(r_num))
    if (numeric_like) {
      bad <- xor(is.na(c_num), is.na(r_num)) |
        (!is.na(c_num) & !is.na(r_num) & abs(c_num - r_num) > tolerance)
      if (any(bad)) {
        return(list(status = "different", details = paste0("numeric column differs: ", nm)))
      }
    } else if (!identical(as.character(cx), as.character(rx))) {
      return(list(status = "different", details = paste0("character column differs: ", nm)))
    }
  }
  list(status = "ok", details = "")
}

.hcr_compare_outputs <- function(output_dir, reference_dir, tolerance = 1e-8) {
  files <- c(
    "qc_summary.csv",
    "combined_edgelist.csv",
    "gfc_matrix.csv",
    "cluster_information.csv",
    "module_label_map.csv",
    "module_gene_list.csv",
    "heatmap_matrix.csv",
    "font_probe.csv",
    "plot_variants.csv",
    "plot_variant_parameters.csv",
    "split_resolution_summary.csv",
    "split_resolution_by_module.csv",
    "split_resolved_modules.csv",
    "split_summary.csv",
    "enrichment_all.csv",
    "enrichment_selected.csv",
    "enrichment_significant.csv"
  )
  rows <- lapply(files, function(file) {
    res <- .hcr_compare_table(
      current_path = file.path(output_dir, file),
      reference_path = file.path(reference_dir, file),
      tolerance = tolerance
    )
    data.frame(file = file, status = res$status, details = res$details, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.hcr_update_reference <- function(output_dir, reference_dir) {
  dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(output_dir, pattern = "\\.(csv|json)$", full.names = TRUE)
  ok <- file.copy(files, reference_dir, overwrite = TRUE)
  if (!all(ok)) {
    stop("Could not update all reference files in: ", reference_dir, call. = FALSE)
  }
  invisible(normalizePath(reference_dir, winslash = "/", mustWork = TRUE))
}

run_hcocena_realdata_regression <- function(mode = c("quick", "full"),
                                            data_dir = NULL,
                                            output_dir = NULL,
                                            reference_dir = NULL,
                                            update_reference = FALSE,
                                            compare_reference = TRUE,
                                            load_source = TRUE,
                                            seed = 168575,
                                            tolerance = 1e-8,
                                            repo_root = NULL) {
  mode <- match.arg(mode)
  repo_root <- if (is.null(repo_root)) .hcr_repo_root() else normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  cfg <- .hcr_mode_config(mode)
  data_dir <- .hcr_scalar(data_dir, .hcr_default_data_dir(repo_root))
  output_dir <- .hcr_scalar(output_dir, file.path(repo_root, "realdata-output", mode))
  reference_dir <- .hcr_scalar(reference_dir, .hcr_default_reference_dir(repo_root, mode))

  .hcr_check_inputs(data_dir, reference_dir)
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  reference_dir <- normalizePath(reference_dir, winslash = "/", mustWork = FALSE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  hc_output_dir <- file.path(output_dir, "hcocena")
  dir.create(hc_output_dir, recursive = TRUE, showWarnings = FALSE)

  package_source <- .hcr_load_package(repo_root, load_source = load_source)
  set.seed(as.integer(seed))

  reference_files <- file.path(repo_root, "docker", "reference_files")
  if (!dir.exists(reference_files)) {
    reference_files <- file.path(repo_root, "inst", "extdata")
  }

  message("Running hCoCena real-data regression: mode=", mode)
  message("Data dir: ", data_dir)
  message("Output dir: ", output_dir)

  hc <- hcocena::hc_init()
  hc <- hcocena::hc_set_paths(
    hc,
    dir_count_data = .hcr_slash(data_dir),
    dir_annotation = .hcr_slash(data_dir),
    dir_reference_files = .hcr_slash(reference_files),
    dir_output = .hcr_slash(hc_output_dir)
  )
  hc <- hcocena::hc_define_layers(
    hc,
    data_sets = list(
      RNA_Seq = c("data_seq_processed.txt", "annotation_seq.txt"),
      Array = c("data_array_processed.txt", "annotation_array.txt")
    )
  )
  hc <- hcocena::hc_read_data(
    hc,
    sep_counts = "\t",
    sep_anno = "\t",
    gene_symbol_col = "SYMBOL",
    sample_col = "SampleID",
    count_has_rn = FALSE,
    anno_has_rn = FALSE,
    project_folder = ""
  )

  if (isTRUE(cfg$read_supplementary)) {
    hc <- hcocena::hc_set_supp_files(
      hc,
      Tf = "TFcat.txt",
      Hallmark = "h.all.v2023.1.Hs.symbols.gmt",
      Go = "c5.go.bp.v2023.1.Hs.symbols.gmt",
      Kegg = "c2.cp.kegg.v2023.1.Hs.symbols.gmt",
      Reactome = "c2.cp.reactome.v2023.1.Hs.symbols.gmt"
    )
    hc <- hcocena::hc_read_supplementary(hc)
  }

  hc <- hcocena::hc_set_global_settings(
    hc,
    organism = "human",
    control_keyword = "baseline",
    variable_of_interest = "merged",
    min_nodes_number_for_network = cfg$min_nodes_number_for_network,
    min_nodes_number_for_cluster = cfg$min_nodes_number_for_cluster,
    range_GFC = 2.0,
    layout_algorithm = "layout_with_fr",
    data_in_log = TRUE
  )
  hc <- hcocena::hc_set_layer_settings(
    hc,
    top_var = cfg$top_var,
    min_corr = cfg$min_corr,
    range_cutoff_length = cfg$range_cutoff_length,
    print_distribution_plots = c(FALSE, FALSE)
  )
  hc <- hcocena::hc_run_expression_analysis_1(hc, corr_method = "pearson")
  hc <- hcocena::hc_set_cutoff(hc, cutoff_vector = cfg$cutoffs, verbose = FALSE)
  hc <- hcocena::hc_run_expression_analysis_2(hc, plot_HM = isTRUE(cfg$plot_layer_heatmaps))
  hc <- hcocena::hc_build_integrated_network(hc, mode = "u", multi_edges = "min")
  hc <- hcocena::hc_cluster_calculation(
    hc,
    cluster_algo = cfg$cluster_algo,
    no_of_iterations = cfg$no_of_iterations,
    resolution = cfg$resolution
  )

  plot_variants <- list()
  heatmap_plot_defaults <- .hcr_heatmap_plot_defaults()
  default_heatmap <- file.path(hc_output_dir, paste0("cluster_heatmap_", mode, ".pdf"))
  default_heatmap_parameters <- list(
    module_label_preset = "balanced",
    module_label_fontsize = NULL,
    module_label_pt_size = NULL,
    module_box_width_cm = NULL,
    gene_count_mode = "text"
  )
  hc <- hcocena::hc_plot_cluster_heatmap(
    hc,
    file_name = basename(default_heatmap),
    module_label_preset = default_heatmap_parameters$module_label_preset,
    module_label_fontsize = default_heatmap_parameters$module_label_fontsize,
    module_label_pt_size = default_heatmap_parameters$module_label_pt_size,
    module_box_width_cm = default_heatmap_parameters$module_box_width_cm,
    gene_count_mode = default_heatmap_parameters$gene_count_mode,
    return_HM = FALSE
  )
  plot_variants[[length(plot_variants) + 1L]] <- .hcr_plot_variant(
    plot = "Cluster heatmap",
    variant = "baseline before module split",
    pdf_file = default_heatmap,
    base_dir = output_dir,
    stage = "post clustering",
    parameters = default_heatmap_parameters,
    defaults = heatmap_plot_defaults
  )

  split_target <- NULL
  if (isTRUE(cfg$run_module_split)) {
    split_target <- .hcr_choose_split_module(hc)
    message(
      "Testing module split on ",
      split_target$module_label,
      " (", split_target$gene_count, " genes)."
    )
    hc <- hcocena::hc_split_modules(
      hc,
      modules = split_target$module_label,
      cluster_algo = cfg$split_cluster_algo,
      no_of_iterations = 1,
      resolution = cfg$split_resolution,
      resolution_grid = cfg$split_resolution_grid,
      resolution_test_only = TRUE,
      seed = seed,
      min_submodule_size = cfg$split_min_submodule_size,
      verbose = FALSE
    )
    hc <- hcocena::hc_split_modules(
      hc,
      modules = split_target$module_label,
      cluster_algo = cfg$split_cluster_algo,
      no_of_iterations = 1,
      resolution = cfg$split_resolution,
      resolution_grid = cfg$split_resolution_grid,
      resolution_test_only = FALSE,
      seed = seed,
      min_submodule_size = cfg$split_min_submodule_size,
      verbose = FALSE
    )
  }

  font_probe_suffix <- if (isTRUE(cfg$run_module_split)) {
    "_post_split_module_label_size_probe.pdf"
  } else {
    "_module_label_size_probe.pdf"
  }
  font_probe_file <- file.path(hc_output_dir, paste0("cluster_heatmap_", mode, font_probe_suffix))
  font_probe_heatmap_parameters <- list(
    module_label_preset = "balanced",
    module_label_fontsize = cfg$module_label_fontsize,
    module_label_pt_size = cfg$module_label_pt_size,
    module_box_width_cm = cfg$module_box_width_cm,
    gene_count_mode = "text"
  )
  hc <- hcocena::hc_plot_cluster_heatmap(
    hc,
    file_name = basename(font_probe_file),
    module_label_preset = font_probe_heatmap_parameters$module_label_preset,
    module_label_fontsize = font_probe_heatmap_parameters$module_label_fontsize,
    module_label_pt_size = font_probe_heatmap_parameters$module_label_pt_size,
    module_box_width_cm = font_probe_heatmap_parameters$module_box_width_cm,
    gene_count_mode = font_probe_heatmap_parameters$gene_count_mode,
    return_HM = TRUE
  )
  plot_variants[[length(plot_variants) + 1L]] <- .hcr_plot_variant(
    plot = "Cluster heatmap",
    variant = if (isTRUE(cfg$run_module_split)) {
      "module-label font-size probe after split"
    } else {
      "module-label font-size probe"
    },
    pdf_file = font_probe_file,
    base_dir = output_dir,
    stage = if (isTRUE(cfg$run_module_split)) "post module split" else "post clustering",
    parameters = font_probe_heatmap_parameters,
    defaults = heatmap_plot_defaults
  )

  if (isTRUE(cfg$run_enrichment)) {
    enrichment_parameters <- list(
      gene_sets = cfg$enrichment_gene_sets,
      top = cfg$enrichment_top,
      consistent_terms = TRUE,
      cluster_columns = FALSE,
      store_panel_objects = "never",
      heatmap_module_label_fontsize = cfg$module_label_fontsize
    )
    hc <- hcocena::hc_functional_enrichment(
      hc,
      gene_sets = enrichment_parameters$gene_sets,
      top = enrichment_parameters$top,
      consistent_terms = enrichment_parameters$consistent_terms,
      cluster_columns = enrichment_parameters$cluster_columns,
      store_panel_objects = enrichment_parameters$store_panel_objects,
      heatmap_module_label_fontsize = enrichment_parameters$heatmap_module_label_fontsize
    )
    plot_variants[[length(plot_variants) + 1L]] <- .hcr_plot_variant(
      plot = "Functional enrichment",
      variant = "combined all DBs",
      pdf_file = file.path(hc_output_dir, paste0("Enrichment_All_DBs_top_", cfg$enrichment_top, ".pdf")),
      base_dir = output_dir,
      stage = if (isTRUE(cfg$run_module_split)) "post module split" else "post clustering",
      parameters = enrichment_parameters,
      defaults = .hcr_enrichment_plot_defaults()
    )
  }

  export <- .hcr_export_outputs(
    hc,
    output_dir,
    mode = mode,
    font_probe_file = font_probe_file,
    plot_variants = plot_variants
  )
  visual_report <- .hcr_write_visual_report(
    hc = hc,
    output_dir = output_dir,
    mode = mode,
    font_probe_file = font_probe_file,
    split_target = split_target,
    parameters = c(
      list(
        mode = mode,
        seed = seed,
        data_dir = data_dir,
        output_dir = output_dir,
        reference_dir = reference_dir,
        package_source = package_source,
        package_version = as.character(utils::packageVersion("hcocena")),
        compare_reference = compare_reference,
        update_reference = update_reference,
        tolerance = tolerance
      ),
      cfg
    ),
    plot_variants = plot_variants
  )
  export$files[["visual_report"]] <- visual_report
  compare <- NULL
  if (isTRUE(update_reference)) {
    .hcr_update_reference(output_dir, reference_dir)
  }
  if (isTRUE(compare_reference) && dir.exists(reference_dir) && !isTRUE(update_reference)) {
    compare <- .hcr_compare_outputs(output_dir, reference_dir, tolerance = tolerance)
    .hcr_write_csv(compare, file.path(output_dir, "reference_comparison.csv"))
  }

  manifest <- list(
    mode = mode,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    repo_root = repo_root,
    data_dir = data_dir,
    output_dir = output_dir,
    reference_dir = reference_dir,
    package_source = package_source,
    package_version = as.character(utils::packageVersion("hcocena")),
    seed = as.integer(seed),
    config = cfg,
    split_target = split_target,
    files = export$files,
    reference_updated = isTRUE(update_reference),
    reference_compared = !is.null(compare),
    reference_ok = if (is.null(compare)) NA else all(compare$status == "ok")
  )
  manifest_path <- .hcr_write_json(manifest, file.path(output_dir, "manifest.json"))
  md5_files <- unlist(c(export$files, list(manifest = manifest_path)), use.names = TRUE)
  md5 <- data.frame(
    file = names(md5_files),
    path = as.character(md5_files),
    md5 = as.character(tools::md5sum(md5_files)),
    stringsAsFactors = FALSE
  )
  .hcr_write_csv(md5, file.path(output_dir, "checksums.csv"))

  result <- list(
    hc = hc,
    manifest = manifest,
    output_dir = output_dir,
    reference_dir = reference_dir,
    comparison = compare,
    files = export$files
  )
  class(result) <- c("hcocena_realdata_regression", class(result))
  result
}

.hcr_main <- function() {
  args <- .hcr_parse_args()
  if (.hcr_bool(args$help, FALSE)) {
    cat(
      "Usage: Rscript scripts/run_realdata_regression.R [--mode quick|full]\n",
      "       [--data-dir PATH] [--output-dir PATH] [--reference-dir PATH]\n",
      "       [--update-reference] [--no-compare-reference] [--load-source true|false]\n",
      "\n",
      "Environment:\n",
      "  HCOCENA_REALDATA_DIR             Real data directory override.\n",
      "  HCOCENA_REALDATA_REFERENCE_DIR   Parent directory for references.\n",
      "  HCOCENA_REALDATA_FULL_ENRICHMENT true/false for full-mode enrichment.\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  mode <- .hcr_scalar(args$mode, Sys.getenv("HCOCENA_REALDATA_MODE", unset = "quick"))
  update_reference <- .hcr_bool(args$update_reference, FALSE)
  compare_reference <- !.hcr_bool(args$no_compare_reference, FALSE)
  load_source <- .hcr_bool(args$load_source, TRUE)
  seed <- as.integer(.hcr_scalar(args$seed, Sys.getenv("HCOCENA_REALDATA_SEED", unset = "168575")))
  tolerance <- as.numeric(.hcr_scalar(args$tolerance, "1e-8"))

  res <- run_hcocena_realdata_regression(
    mode = mode,
    data_dir = args$data_dir,
    output_dir = args$output_dir,
    reference_dir = args$reference_dir,
    update_reference = update_reference,
    compare_reference = compare_reference,
    load_source = load_source,
    seed = seed,
    tolerance = tolerance
  )

  cat("Real-data regression completed.\n")
  cat("Output: ", res$output_dir, "\n", sep = "")
  if (!is.null(res$comparison)) {
    cat("Reference comparison:\n")
    print(res$comparison, row.names = FALSE)
    if (!all(res$comparison$status == "ok")) {
      quit(status = 1)
    }
  }
  invisible(res)
}

if (identical(sys.nframe(), 0L)) {
  .hcr_main()
}
