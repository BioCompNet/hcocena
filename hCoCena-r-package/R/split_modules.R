#' Split one or multiple modules into submodules
#'
#' Re-clusters genes inside selected modules and replaces each selected module
#' with module-specific submodules (color shades of the parent module).
#'
#' This is useful for very large modules that likely contain multiple
#' biological programs.
#'
#' Accepted values in `modules`:
#' - module labels from `module_label_map` (e.g. `"M3"`),
#' - module colors (e.g. `"#FFD700"`),
#' - numeric module indices (based on current included-module order).
#'
#' The function stores an undo snapshot in
#' `hcobject$satellite_outputs$module_split_history`.
#'
#' @param modules Character/numeric vector of modules to split.
#' @param cluster_algo Clustering algorithm for within-module splitting.
#'  One of `"cluster_leiden"` (default), `"cluster_louvain"`,
#'  `"cluster_fast_greedy"`, `"cluster_infomap"`, `"cluster_walktrap"`,
#'  `"cluster_label_prop"` or `"auto"`.
#' @param no_of_iterations Number of Leiden iterations (used only for Leiden).
#' @param resolution Leiden resolution (used only for Leiden).
#' @param resolution_grid Optional numeric vector of candidate resolutions to
#'  test before splitting. For each candidate, hCoCena reports how many
#'  submodules would be retained after size filtering.
#' @param resolution_test_only Logical; if `TRUE`, only run the resolution test
#'  (when `resolution_grid` is set) and do not apply any split.
#' @param partition_type Leiden partition type (used only for Leiden).
#' @param seed Random seed used for deterministic clustering.
#' @param drop_small_submodules Logical; if `TRUE` (default), submodules with
#'  fewer than `min_submodule_size` genes are dropped from the module
#'  annotation (their genes become unassigned/white in network visualizations).
#' @param min_submodule_size Optional minimum size for retained split
#'  submodules. If `NULL`, uses `global_settings$min_nodes_number_for_cluster`.
#' @param min_module_size Optional alias for `min_submodule_size`. If set, it
#'  takes precedence.
#' @param verbose Logical; print progress messages.
#' @export
split_modules <- function(modules,
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
  .hc_legacy_warning("split_modules")

  if (missing(modules) || is.null(modules) || length(modules) == 0) {
    stop("`modules` must contain at least one module identifier.")
  }
  if (!is.character(cluster_algo) || length(cluster_algo) != 1 || is.na(cluster_algo)) {
    stop("`cluster_algo` must be a single character string.")
  }
  if (!is.numeric(no_of_iterations) || length(no_of_iterations) != 1 || !is.finite(no_of_iterations) || no_of_iterations < 1) {
    stop("`no_of_iterations` must be a positive numeric scalar.")
  }
  if (!is.numeric(resolution) || length(resolution) != 1 || !is.finite(resolution) || resolution <= 0) {
    stop("`resolution` must be a positive numeric scalar.")
  }
  if (!is.null(resolution_grid) &&
      (!is.numeric(resolution_grid) || length(resolution_grid) == 0 || any(!is.finite(resolution_grid)) || any(resolution_grid <= 0))) {
    stop("`resolution_grid` must be NULL or a numeric vector with positive finite values.")
  }
  if (!is.logical(resolution_test_only) || length(resolution_test_only) != 1 || is.na(resolution_test_only)) {
    stop("`resolution_test_only` must be TRUE or FALSE.")
  }
  if (!is.character(partition_type) || length(partition_type) != 1 || is.na(partition_type)) {
    stop("`partition_type` must be a single character string.")
  }
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
    stop("`seed` must be a finite numeric scalar.")
  }
  if (!is.logical(drop_small_submodules) || length(drop_small_submodules) != 1 || is.na(drop_small_submodules)) {
    stop("`drop_small_submodules` must be TRUE or FALSE.")
  }
  if (!is.null(min_submodule_size) &&
      (!is.numeric(min_submodule_size) || length(min_submodule_size) != 1 || !is.finite(min_submodule_size) || min_submodule_size < 1)) {
    stop("`min_submodule_size` must be NULL or a positive numeric scalar.")
  }
  if (!is.null(min_module_size) &&
      (!is.numeric(min_module_size) || length(min_module_size) != 1 || !is.finite(min_module_size) || min_module_size < 1)) {
    stop("`min_module_size` must be NULL or a positive numeric scalar.")
  }
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  merged_net <- hcobject[["integrated_output"]][["merged_net"]]
  if (is.null(merged_net) || !inherits(merged_net, "igraph")) {
    stop("Integrated graph missing. Run `build_integrated_network()` first.")
  }

  cluster_calc <- hcobject[["integrated_output"]][["cluster_calc"]]
  cluster_info <- cluster_calc[["cluster_information"]]
  if (is.null(cluster_info) || !is.data.frame(cluster_info) || nrow(cluster_info) == 0) {
    stop("No cluster information found. Run `cluster_calculation()` first.")
  }
  required_cols <- c("color", "gene_n")
  if (!all(required_cols %in% colnames(cluster_info))) {
    stop("`cluster_information` is missing required columns: ", paste(required_cols, collapse = ", "))
  }

  cluster_info <- as.data.frame(cluster_info, stringsAsFactors = FALSE)
  cluster_included <- if ("cluster_included" %in% colnames(cluster_info)) {
    as.character(cluster_info$cluster_included) == "yes"
  } else {
    rep(TRUE, nrow(cluster_info))
  }
  incl_info <- cluster_info[cluster_included, , drop = FALSE]
  if (nrow(incl_info) == 0) {
    stop("No included modules found in `cluster_information`.")
  }
  if (is.null(min_submodule_size)) {
    min_submodule_size <- suppressWarnings(
      as.numeric(hcobject[["global_settings"]][["min_nodes_number_for_cluster"]])
    )
  }
  if (!is.null(min_module_size)) {
    min_submodule_size <- min_module_size
  }
  if (!is.finite(min_submodule_size) || min_submodule_size < 1) {
    min_submodule_size <- 15
  }
  min_submodule_size <- as.integer(round(min_submodule_size))

  module_prefix <- if (!is.null(cluster_calc[["module_prefix"]]) && length(cluster_calc[["module_prefix"]]) == 1) {
    as.character(cluster_calc[["module_prefix"]])
  } else {
    "M"
  }

  module_label_map <- cluster_calc[["module_label_map"]]
  if (is.null(module_label_map) || length(module_label_map) == 0) {
    module_label_map <- stats::setNames(
      paste0(module_prefix, seq_len(nrow(incl_info))),
      as.character(incl_info$color)
    )
  } else {
    module_label_map <- as.character(module_label_map)
    names(module_label_map) <- as.character(names(cluster_calc[["module_label_map"]]))
  }

  resolved <- .hc_resolve_modules_for_split(
    modules = modules,
    available_colors = as.character(incl_info$color),
    module_label_map = module_label_map
  )
  unresolved_tbl <- resolved$resolution_table[
    as.character(resolved$resolution_table$status) != "ok",
    ,
    drop = FALSE
  ]
  if (nrow(unresolved_tbl) > 0) {
    unresolved_inputs <- unique(as.character(unresolved_tbl$input))
    available_labels <- unique(as.character(module_label_map))
    if (length(available_labels) > 0) {
      available_preview <- paste(utils::head(available_labels, 12L), collapse = ", ")
      if (length(available_labels) > 12L) {
        available_preview <- paste0(available_preview, ", ...")
      }
      stop(
        "Could not resolve the following module identifier(s): ",
        paste(unresolved_inputs, collapse = ", "),
        "\nUse exact current module labels (prefix-sensitive), module colors, or numeric indices.",
        "\nAvailable module labels: ", available_preview
      )
    }
    stop(
      "Could not resolve the following module identifier(s): ",
      paste(unresolved_inputs, collapse = ", "),
      "\nUse exact current module labels (prefix-sensitive), module colors, or numeric indices."
    )
  }
  target_colors <- resolved$target_colors
  if (length(target_colors) == 0) {
    stop(
      "Could not match any requested module in `modules`.\n",
      "Use module labels (e.g. M3), module colors, or module indices."
    )
  }

  if (isTRUE(verbose)) {
    message(
      "split_modules(): splitting ",
      paste(resolved$resolved_labels, collapse = ", "),
      " (", length(target_colors), " module(s))."
    )
  }

  all_graph_nodes <- as.character(igraph::V(merged_net)$name)

  if (!is.null(resolution_grid)) {
    resolution_grid <- unique(as.numeric(resolution_grid))
    resolution_grid <- resolution_grid[is.finite(resolution_grid) & resolution_grid > 0]

    preview <- .hc_preview_split_resolutions(
      merged_net = merged_net,
      cluster_info = cluster_info,
      target_colors = target_colors,
      module_label_map = module_label_map,
      cluster_algo = cluster_algo,
      no_of_iterations = as.integer(round(no_of_iterations)),
      partition_type = partition_type,
      seed = as.integer(round(seed)),
      resolution_grid = resolution_grid,
      drop_small_submodules = drop_small_submodules,
      min_submodule_size = min_submodule_size,
      all_graph_nodes = all_graph_nodes
    )

    sat <- hcobject[["satellite_outputs"]]
    if (is.null(sat) || !is.list(sat)) {
      sat <- list()
    }
    sat[["module_split_resolution_test_last"]] <- list(
      timestamp = as.character(Sys.time()),
      parameters = list(
        modules = as.character(modules),
        cluster_algo = cluster_algo,
        no_of_iterations = as.integer(round(no_of_iterations)),
        resolution_grid = resolution_grid,
        partition_type = partition_type,
        seed = as.integer(round(seed)),
        drop_small_submodules = drop_small_submodules,
        min_submodule_size = min_submodule_size
      ),
      resolved_modules = resolved$resolution_table,
      by_module = preview$by_module,
      by_resolution = preview$by_resolution
    )
    hcobject[["satellite_outputs"]] <<- sat

    if (isTRUE(verbose)) {
      message(
        "split_modules(): tested ", length(resolution_grid),
        " resolution value(s); summary stored in ",
        "`hcobject$satellite_outputs$module_split_resolution_test_last`."
      )
      print(preview$by_resolution)
    }

    if (isTRUE(resolution_test_only)) {
      if (isTRUE(verbose)) {
        message("split_modules(): `resolution_test_only = TRUE`; no split applied.")
      }
      return(invisible(NULL))
    }
  }

  gfc_all <- hcobject[["integrated_output"]][["GFC_all_layers"]]
  if (is.null(gfc_all) || !is.data.frame(gfc_all) || !("Gene" %in% colnames(gfc_all))) {
    stop("`integrated_output$GFC_all_layers` is missing or malformed.")
  }

  before_cluster_info <- cluster_info
  before_module_label_map <- module_label_map

  new_rows <- vector("list", nrow(cluster_info))
  row_ptr <- 0L
  split_summary_rows <- list()
  used_colors <- as.character(cluster_info$color)
  total_created_submodules <- 0L
  total_removed_small_submodules <- 0L
  total_removed_small_genes <- 0L
  total_removed_small_parent_modules <- 0L

  for (i in seq_len(nrow(cluster_info))) {
    row_i <- cluster_info[i, , drop = FALSE]
    this_color <- as.character(row_i$color[[1]])

    if (!(this_color %in% target_colors)) {
      row_ptr <- row_ptr + 1L
      new_rows[[row_ptr]] <- row_i
      next
    }

    parent_label <- if (this_color %in% names(module_label_map)) {
      as.character(module_label_map[[this_color]])
    } else {
      as.character(this_color)
    }

    genes <- .hc_parse_genes_from_gene_n(row_i$gene_n[[1]])
    genes <- genes[genes %in% all_graph_nodes]
    genes <- unique(genes)

    if (length(genes) < 3) {
      row_ptr <- row_ptr + 1L
      new_rows[[row_ptr]] <- row_i
      split_summary_rows[[length(split_summary_rows) + 1L]] <- data.frame(
        parent_color = this_color,
        parent_label = parent_label,
        parent_genes = length(genes),
        n_submodules_raw = 1L,
        n_submodules = 1L,
        removed_small_submodules = 0L,
        removed_small_genes = 0L,
        status = "skipped_too_few_genes",
        stringsAsFactors = FALSE
      )
      next
    }

    subgraph <- igraph::induced_subgraph(merged_net, vids = genes)
    membership <- .hc_split_membership(
      graph_obj = subgraph,
      cluster_algo = cluster_algo,
      no_of_iterations = as.integer(round(no_of_iterations)),
      resolution = resolution,
      partition_type = partition_type,
      seed = as.integer(round(seed))
    )

    if (length(membership) == 0) {
      row_ptr <- row_ptr + 1L
      new_rows[[row_ptr]] <- row_i
      split_summary_rows[[length(split_summary_rows) + 1L]] <- data.frame(
        parent_color = this_color,
        parent_label = parent_label,
        parent_genes = length(genes),
        n_submodules_raw = 1L,
        n_submodules = 1L,
        removed_small_submodules = 0L,
        removed_small_genes = 0L,
        status = "skipped_failed_membership",
        stringsAsFactors = FALSE
      )
      next
    }

    member_tbl <- sort(table(membership), decreasing = TRUE)
    member_levels <- names(member_tbl)
    member_sizes <- as.integer(member_tbl)
    n_sub_raw <- length(member_levels)
    kept_idx <- seq_along(member_levels)
    removed_small_submodules <- 0L
    removed_small_genes <- 0L
    if (isTRUE(drop_small_submodules)) {
      keep_mask <- member_sizes >= min_submodule_size
      kept_idx <- which(keep_mask)
      removed_small_submodules <- sum(!keep_mask)
      if (removed_small_submodules > 0) {
        removed_small_genes <- sum(member_sizes[!keep_mask])
      }
    }
    if (length(kept_idx) < 2) {
      if (isTRUE(drop_small_submodules) && length(genes) < min_submodule_size) {
        total_removed_small_parent_modules <- total_removed_small_parent_modules + 1L
        total_removed_small_genes <- total_removed_small_genes + length(genes)
        split_summary_rows[[length(split_summary_rows) + 1L]] <- data.frame(
          parent_color = this_color,
          parent_label = parent_label,
          parent_genes = length(genes),
          n_submodules_raw = n_sub_raw,
          n_submodules = 0L,
          removed_small_submodules = removed_small_submodules,
          removed_small_genes = length(genes),
          status = "dropped_parent_below_min_size",
          stringsAsFactors = FALSE
        )
        next
      }
      row_ptr <- row_ptr + 1L
      new_rows[[row_ptr]] <- row_i
      split_summary_rows[[length(split_summary_rows) + 1L]] <- data.frame(
        parent_color = this_color,
        parent_label = parent_label,
        parent_genes = length(genes),
        n_submodules_raw = n_sub_raw,
        n_submodules = 1L,
        removed_small_submodules = removed_small_submodules,
        removed_small_genes = removed_small_genes,
        status = if (isTRUE(drop_small_submodules)) "skipped_after_small_module_filter" else "skipped_single_submodule",
        stringsAsFactors = FALSE
      )
      next
    }
    member_levels <- member_levels[kept_idx]
    n_sub <- length(member_levels)
    total_removed_small_submodules <- total_removed_small_submodules + removed_small_submodules
    total_removed_small_genes <- total_removed_small_genes + removed_small_genes

    child_colors <- .hc_generate_module_shades(
      base_color = this_color,
      n = n_sub,
      avoid = used_colors
    )
    used_colors <- c(used_colors, child_colors)
    child_labels <- paste0(parent_label, ".", seq_len(n_sub))
    names(child_labels) <- child_colors

    child_rows <- .hc_build_child_cluster_rows(
      template_row = row_i,
      membership = membership,
      member_levels = member_levels,
      child_colors = child_colors,
      child_labels = child_labels,
      gfc_all = gfc_all
    )
    for (j in seq_len(nrow(child_rows))) {
      row_ptr <- row_ptr + 1L
      new_rows[[row_ptr]] <- child_rows[j, , drop = FALSE]
    }
    total_created_submodules <- total_created_submodules + n_sub

    module_label_map <- module_label_map[setdiff(names(module_label_map), this_color)]
    module_label_map <- c(module_label_map, child_labels)

    split_summary_rows[[length(split_summary_rows) + 1L]] <- data.frame(
      parent_color = this_color,
      parent_label = parent_label,
      parent_genes = length(genes),
      n_submodules_raw = n_sub_raw,
      n_submodules = n_sub,
      removed_small_submodules = removed_small_submodules,
      removed_small_genes = removed_small_genes,
      status = "split_ok",
      stringsAsFactors = FALSE
    )
  }

  if (row_ptr == 0L) {
    stop("No module rows were retained after split; this should never happen.")
  }

  new_cluster_info <- do.call(rbind, new_rows[seq_len(row_ptr)])
  rownames(new_cluster_info) <- NULL

  included_new <- if ("cluster_included" %in% colnames(new_cluster_info)) {
    as.character(new_cluster_info$color[as.character(new_cluster_info$cluster_included) == "yes"])
  } else {
    as.character(new_cluster_info$color)
  }
  included_new <- unique(included_new)
  if (length(included_new) == 0) {
    stop("Split produced no included modules.")
  }

  missing_map <- setdiff(included_new, names(module_label_map))
  if (length(missing_map) > 0) {
    map_start <- length(module_label_map) + 1L
    module_label_map <- c(
      module_label_map,
      stats::setNames(
        paste0(module_prefix, seq.int(map_start, length.out = length(missing_map))),
        missing_map
      )
    )
  }
  module_label_map <- module_label_map[included_new]

  hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- new_cluster_info
  hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]] <<- module_label_map
  hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]] <<- NULL

  .hc_clear_outputs_after_module_change()

  split_summary <- if (length(split_summary_rows) > 0) {
    do.call(rbind, split_summary_rows)
  } else {
    data.frame()
  }

  sat <- hcobject[["satellite_outputs"]]
  if (is.null(sat) || !is.list(sat)) {
    sat <- list()
  }
  hist <- sat[["module_split_history"]]
  if (is.null(hist) || !is.list(hist)) {
    hist <- list()
  }

  history_entry <- list(
    timestamp = as.character(Sys.time()),
    parameters = list(
      modules = as.character(modules),
      cluster_algo = cluster_algo,
      no_of_iterations = as.integer(round(no_of_iterations)),
      resolution = resolution,
      partition_type = partition_type,
      seed = as.integer(round(seed))
      ,
      drop_small_submodules = drop_small_submodules,
      min_submodule_size = min_submodule_size,
      min_module_size = min_module_size,
      resolution_grid = resolution_grid,
      resolution_test_only = resolution_test_only
    ),
    resolved_modules = resolved$resolution_table,
    split_summary = split_summary,
    before = list(
      cluster_information = before_cluster_info,
      module_label_map = before_module_label_map
    ),
    after = list(
      cluster_information = new_cluster_info,
      module_label_map = module_label_map
    )
  )
  hist[[length(hist) + 1L]] <- history_entry
  sat[["module_split_history"]] <- hist
  sat[["module_split_last"]] <- history_entry
  hcobject[["satellite_outputs"]] <<- sat

  if (isTRUE(verbose)) {
    split_ok <- if (nrow(split_summary) > 0) {
      sum(split_summary$status == "split_ok")
    } else {
      0L
    }
    total_modules_after <- if ("cluster_included" %in% colnames(new_cluster_info)) {
      sum(as.character(new_cluster_info$cluster_included) == "yes")
    } else {
      nrow(new_cluster_info)
    }
    message(
      "split_modules(): completed. ",
      split_ok, " parent module(s) split; ",
      total_created_submodules, " submodule(s) created; ",
      total_removed_small_submodules, " submodule(s) removed (< ", min_submodule_size, " genes); ",
      total_removed_small_parent_modules, " parent module(s) removed (< ", min_submodule_size, " genes); ",
      total_removed_small_genes, " gene(s) dropped by size filter; ",
      "resulting included modules = ", total_modules_after, "; ",
      "history depth = ", length(hist), "."
    )
  }
}

.hc_preview_split_resolutions <- function(merged_net,
                                          cluster_info,
                                          target_colors,
                                          module_label_map,
                                          cluster_algo,
                                          no_of_iterations,
                                          partition_type,
                                          seed,
                                          resolution_grid,
                                          drop_small_submodules,
                                          min_submodule_size,
                                          all_graph_nodes) {
  if (length(target_colors) == 0 || length(resolution_grid) == 0) {
    return(list(
      by_module = data.frame(),
      by_resolution = data.frame()
    ))
  }

  cluster_info <- as.data.frame(cluster_info, stringsAsFactors = FALSE)
  cluster_included <- if ("cluster_included" %in% colnames(cluster_info)) {
    as.character(cluster_info$cluster_included) == "yes"
  } else {
    rep(TRUE, nrow(cluster_info))
  }
  base_included_modules <- sum(cluster_included)

  rows <- list()
  for (color in target_colors) {
    idx <- which(as.character(cluster_info$color) == as.character(color))
    if (length(idx) == 0) {
      next
    }
    row_i <- cluster_info[idx[[1]], , drop = FALSE]
    parent_label <- if (color %in% names(module_label_map)) {
      as.character(module_label_map[[color]])
    } else {
      as.character(color)
    }
    genes <- .hc_parse_genes_from_gene_n(row_i$gene_n[[1]])
    genes <- genes[genes %in% all_graph_nodes]
    genes <- unique(genes)

    for (res in resolution_grid) {
      n_sub_raw <- 1L
      n_sub_kept <- 1L
      removed_small_submodules <- 0L
      removed_small_genes <- 0L
      would_split <- FALSE
      status <- "skipped_too_few_genes"

      if (length(genes) >= 3) {
        subgraph <- igraph::induced_subgraph(merged_net, vids = genes)
        mem <- .hc_split_membership(
          graph_obj = subgraph,
          cluster_algo = cluster_algo,
          no_of_iterations = no_of_iterations,
          resolution = as.numeric(res),
          partition_type = partition_type,
          seed = seed
        )

        if (!is.null(mem) && length(mem) > 0) {
          member_sizes <- as.integer(sort(table(mem), decreasing = TRUE))
          n_sub_raw <- length(member_sizes)
          n_sub_kept <- n_sub_raw
          if (isTRUE(drop_small_submodules)) {
            keep_mask <- member_sizes >= min_submodule_size
            n_sub_kept <- sum(keep_mask)
            removed_small_submodules <- sum(!keep_mask)
            if (removed_small_submodules > 0) {
              removed_small_genes <- sum(member_sizes[!keep_mask])
            }
          }
          would_split <- n_sub_kept >= 2
          status <- if (would_split) "would_split" else "would_not_split"
        } else {
          status <- "skipped_failed_membership"
        }
      }

      rows[[length(rows) + 1L]] <- data.frame(
        parent_color = as.character(color),
        parent_label = parent_label,
        parent_genes = length(genes),
        resolution = as.numeric(res),
        n_submodules_raw = as.integer(n_sub_raw),
        n_submodules_kept = as.integer(n_sub_kept),
        removed_small_submodules = as.integer(removed_small_submodules),
        removed_small_genes = as.integer(removed_small_genes),
        would_split = isTRUE(would_split),
        status = status,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    return(list(
      by_module = data.frame(),
      by_resolution = data.frame()
    ))
  }
  by_module <- do.call(rbind, rows)
  rownames(by_module) <- NULL

  res_levels <- unique(by_module$resolution)
  by_res_rows <- vector("list", length(res_levels))
  for (i in seq_along(res_levels)) {
    res <- res_levels[[i]]
    sub <- by_module[by_module$resolution == res, , drop = FALSE]
    n_modules_tested <- nrow(sub)
    n_modules_would_split <- sum(sub$would_split)
    total_created_submodules <- sum(ifelse(sub$would_split, sub$n_submodules_kept, 0L))
    total_removed_small_submodules <- sum(sub$removed_small_submodules)
    total_removed_small_genes <- sum(sub$removed_small_genes)
    resulting_included_modules <- base_included_modules - n_modules_would_split + total_created_submodules

    by_res_rows[[i]] <- data.frame(
      resolution = as.numeric(res),
      n_modules_tested = as.integer(n_modules_tested),
      n_modules_would_split = as.integer(n_modules_would_split),
      total_created_submodules = as.integer(total_created_submodules),
      total_removed_small_submodules = as.integer(total_removed_small_submodules),
      total_removed_small_genes = as.integer(total_removed_small_genes),
      resulting_included_modules = as.integer(resulting_included_modules),
      stringsAsFactors = FALSE
    )
  }
  by_resolution <- do.call(rbind, by_res_rows)
  rownames(by_resolution) <- NULL

  list(
    by_module = by_module,
    by_resolution = by_resolution
  )
}

#' Undo module splitting
#'
#' Restores cluster assignments from the split history created by
#' [split_modules()].
#'
#' @param which Either `"last"` (undo one split step) or `"all"` (restore the
#'  original pre-split cluster state).
#' @param verbose Logical; print progress messages.
#' @export
unsplit_modules <- function(which = c("last", "all"), verbose = TRUE) {
  .hc_legacy_warning("unsplit_modules")
  which <- base::match.arg(which)

  sat <- hcobject[["satellite_outputs"]]
  if (is.null(sat) || !is.list(sat) || !("module_split_history" %in% names(sat))) {
    stop("No module split history found.")
  }
  hist <- sat[["module_split_history"]]
  if (is.null(hist) || !is.list(hist) || length(hist) == 0) {
    stop("No module split history found.")
  }

  if (which == "last") {
    entry <- hist[[length(hist)]]
    if (length(hist) == 1) {
      hist <- list()
    } else {
      hist <- hist[-length(hist)]
    }
  } else {
    entry <- hist[[1]]
    hist <- list()
  }

  before <- entry[["before"]]
  if (is.null(before) || !is.list(before) || is.null(before[["cluster_information"]])) {
    stop("Split history entry is malformed: missing `before$cluster_information`.")
  }

  hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- before[["cluster_information"]]
  if (!is.null(before[["module_label_map"]])) {
    hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]] <<- before[["module_label_map"]]
  } else {
    hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]] <<- NULL
  }
  hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]] <<- NULL

  .hc_clear_outputs_after_module_change()

  sat[["module_split_history"]] <- hist
  sat[["module_split_last_restore"]] <- list(
    timestamp = as.character(Sys.time()),
    restored_which = which,
    remaining_history_depth = length(hist)
  )
  if (length(hist) == 0) {
    sat[["module_split_last"]] <- NULL
  } else {
    sat[["module_split_last"]] <- hist[[length(hist)]]
  }
  hcobject[["satellite_outputs"]] <<- sat

  if (isTRUE(verbose)) {
    message(
      "unsplit_modules(): restored `", which, "` split state; ",
      "history depth = ", length(hist), "."
    )
  }
}

.hc_resolve_modules_for_split <- function(modules, available_colors, module_label_map) {
  available_colors <- unique(as.character(available_colors))
  label_to_color <- stats::setNames(names(module_label_map), as.character(module_label_map))

  resolved_colors <- character(0)
  resolution_rows <- list()
  input_vec <- modules
  if (is.numeric(input_vec)) {
    input_vec <- as.list(input_vec)
  } else {
    input_vec <- as.list(as.character(input_vec))
  }

  for (i in seq_along(input_vec)) {
    x <- input_vec[[i]]
    x_chr <- as.character(x)
    status <- "ok"
    resolved <- NA_character_

    if (is.numeric(x) && is.finite(x)) {
      idx <- as.integer(round(x))
      if (idx >= 1 && idx <= length(available_colors)) {
        resolved <- available_colors[[idx]]
      } else {
        status <- "not_found"
      }
    } else {
      x_chr <- trimws(x_chr)
      if (x_chr %in% available_colors) {
        resolved <- x_chr
      } else if (x_chr %in% names(label_to_color)) {
        resolved <- as.character(label_to_color[[x_chr]])
      } else {
        status <- "not_found"
      }
    }

    if (!is.na(resolved) && !(resolved %in% available_colors)) {
      status <- "not_found"
      resolved <- NA_character_
    }

    if (!is.na(resolved)) {
      resolved_colors <- c(resolved_colors, resolved)
    }

    resolved_label <- if (!is.na(resolved) && resolved %in% names(module_label_map)) {
      as.character(module_label_map[[resolved]])
    } else {
      NA_character_
    }

    resolution_rows[[length(resolution_rows) + 1L]] <- data.frame(
      input = x_chr,
      resolved_color = resolved,
      resolved_label = resolved_label,
      status = status,
      stringsAsFactors = FALSE
    )
  }

  resolution_table <- do.call(rbind, resolution_rows)
  resolved_colors <- unique(resolved_colors[!is.na(resolved_colors)])
  resolved_labels <- if (length(resolved_colors) == 0) {
    character(0)
  } else {
    out <- module_label_map[resolved_colors]
    out[is.na(out)] <- resolved_colors[is.na(out)]
    as.character(out)
  }

  list(
    target_colors = resolved_colors,
    resolved_labels = resolved_labels,
    resolution_table = resolution_table
  )
}

.hc_split_membership <- function(graph_obj,
                                 cluster_algo = "cluster_leiden",
                                 no_of_iterations = 2L,
                                 resolution = 0.1,
                                 partition_type = "RBConfigurationVertexPartition",
                                 seed = 168575L) {
  if (igraph::vcount(graph_obj) <= 1) {
    out <- rep(1L, igraph::vcount(graph_obj))
    names(out) <- as.character(igraph::V(graph_obj)$name)
    return(out)
  }

  get_membership_single <- function(algo) {
    if (algo == "cluster_leiden") {
      set.seed(seed)
      part <- leidenbase::leiden_find_partition(
        igraph = graph_obj,
        partition_type = partition_type,
        edge_weights = igraph::E(graph_obj)$weight,
        resolution = resolution,
        num_iter = no_of_iterations,
        seed = seed
      )
      mem <- suppressWarnings(as.integer(part$membership)) + 1L
      names(mem) <- as.character(igraph::V(graph_obj)$name)
      return(mem)
    }

    if (!exists(algo, mode = "function", where = asNamespace("igraph"), inherits = FALSE)) {
      stop("Unknown clustering algorithm: ", algo)
    }
    fun <- getExportedValue("igraph", algo)
    fit <- tryCatch(fun(graph_obj), error = function(e) NULL)
    if (is.null(fit)) {
      return(NULL)
    }
    mem <- tryCatch(igraph::membership(fit), error = function(e) NULL)
    if (is.null(mem) && !is.null(fit$membership)) {
      mem <- fit$membership
    }
    if (is.null(mem)) {
      return(NULL)
    }

    if (is.list(mem) && !is.numeric(mem)) {
      # Older/newer igraph variants can expose a list-of-members representation.
      vec <- rep(NA_integer_, igraph::vcount(graph_obj))
      for (ii in seq_along(mem)) {
        idx <- suppressWarnings(as.integer(mem[[ii]]))
        idx <- idx[is.finite(idx) & idx >= 1 & idx <= length(vec)]
        vec[idx] <- ii
      }
      mem <- vec
    }

    mem <- suppressWarnings(as.integer(mem))
    if (all(is.na(mem))) {
      return(NULL)
    }
    if (is.null(names(mem)) || all(!nzchar(names(mem)))) {
      names(mem) <- as.character(igraph::V(graph_obj)$name)
    }
    mem
  }

  if (cluster_algo == "auto") {
    algos <- c(
      "cluster_leiden",
      "cluster_louvain",
      "cluster_fast_greedy",
      "cluster_infomap",
      "cluster_walktrap",
      "cluster_label_prop"
    )
    best_mem <- NULL
    best_mod <- -Inf
    for (algo in algos) {
      mem <- tryCatch(get_membership_single(algo), error = function(e) NULL)
      if (is.null(mem) || length(mem) != igraph::vcount(graph_obj)) {
        next
      }
      mod <- tryCatch(igraph::modularity(graph_obj, mem), error = function(e) NA_real_)
      if (is.finite(mod) && mod > best_mod) {
        best_mod <- mod
        best_mem <- mem
      }
    }
    if (is.null(best_mem)) {
      out <- rep(1L, igraph::vcount(graph_obj))
      names(out) <- as.character(igraph::V(graph_obj)$name)
      return(out)
    }
    out <- as.integer(as.factor(best_mem))
    names(out) <- names(best_mem)
    return(out)
  }

  mem <- get_membership_single(cluster_algo)
  if (is.null(mem)) {
    out <- rep(1L, igraph::vcount(graph_obj))
    names(out) <- as.character(igraph::V(graph_obj)$name)
    return(out)
  }
  out <- as.integer(as.factor(mem))
  names(out) <- names(mem)
  out
}

.hc_build_child_cluster_rows <- function(template_row,
                                         membership,
                                         member_levels,
                                         child_colors,
                                         child_labels,
                                         gfc_all) {
  out <- template_row[rep(1, length(member_levels)), , drop = FALSE]
  gfc_cols <- setdiff(colnames(gfc_all), "Gene")

  for (k in seq_along(member_levels)) {
    lvl <- member_levels[[k]]
    genes_k <- names(membership)[membership == as.integer(lvl)]
    genes_k <- genes_k[!is.na(genes_k) & genes_k != ""]
    genes_k <- unique(genes_k)

    if ("clusters" %in% colnames(out)) {
      out$clusters[[k]] <- child_labels[[k]]
    }
    if ("gene_no" %in% colnames(out)) {
      out$gene_no[[k]] <- length(genes_k)
    }
    if ("gene_n" %in% colnames(out)) {
      out$gene_n[[k]] <- paste0(genes_k, collapse = ",")
    }
    if ("cluster_included" %in% colnames(out)) {
      out$cluster_included[[k]] <- "yes"
    }
    if ("color" %in% colnames(out)) {
      out$color[[k]] <- child_colors[[k]]
    }
    if ("vertexsize" %in% colnames(out)) {
      out$vertexsize[[k]] <- 3
    }

    if (length(gfc_cols) > 0) {
      if ("conditions" %in% colnames(out)) {
        out$conditions[[k]] <- paste0(gfc_cols, collapse = "#")
      }
      if ("grp_means" %in% colnames(out)) {
        sub <- gfc_all[gfc_all$Gene %in% genes_k, gfc_cols, drop = FALSE]
        means <- if (nrow(sub) > 0) {
          colMeans(sub, na.rm = TRUE)
        } else {
          rep(NA_real_, length(gfc_cols))
        }
        out$grp_means[[k]] <- paste0(round(means, 3), collapse = ",")
      }
    }
  }

  out
}

.hc_parse_genes_from_gene_n <- function(gene_n_value) {
  if (is.null(gene_n_value) || length(gene_n_value) == 0 || is.na(gene_n_value)) {
    return(character(0))
  }
  genes <- unlist(strsplit(as.character(gene_n_value), ",", fixed = TRUE), use.names = FALSE)
  genes <- trimws(as.character(genes))
  genes[genes != ""]
}

.hc_generate_module_shades <- function(base_color, n, avoid = character(0)) {
  if (n <= 1) {
    return(base_color)
  }

  mix_color <- function(col_a, col_b, alpha) {
    a <- grDevices::col2rgb(col_a) / 255
    b <- grDevices::col2rgb(col_b) / 255
    m <- (1 - alpha) * a + alpha * b
    grDevices::rgb(m[1], m[2], m[3])
  }

  dark <- mix_color(base_color, "#000000", 0.25)
  light <- mix_color(base_color, "#FFFFFF", 0.35)
  cols <- grDevices::colorRampPalette(c(dark, base_color, light))(n)

  avoid <- unique(as.character(avoid))
  for (i in seq_along(cols)) {
    if (!(cols[[i]] %in% avoid) && !duplicated(cols)[[i]]) {
      next
    }

    base_hsv <- grDevices::rgb2hsv(grDevices::col2rgb(cols[[i]]) / 255)
    h <- as.numeric(base_hsv[1, 1])
    s <- as.numeric(base_hsv[2, 1])
    v <- as.numeric(base_hsv[3, 1])

    trial <- 0L
    repeat {
      trial <- trial + 1L
      h_new <- (h + (trial * 0.03)) %% 1
      s_new <- min(1, max(0.3, s + ((trial %% 2) * 0.08)))
      v_new <- min(1, max(0.25, v - ((trial %% 3) * 0.06)))
      candidate <- grDevices::hsv(h = h_new, s = s_new, v = v_new)
      if (!(candidate %in% c(avoid, cols[seq_len(max(1L, i - 1L))]))) {
        cols[[i]] <- candidate
        break
      }
      if (trial > 25L) {
        cols[[i]] <- candidate
        break
      }
    }
  }
  cols
}

.hc_clear_outputs_after_module_change <- function() {
  hcobject[["integrated_output"]][["enrichments"]] <<- NULL
  hcobject[["integrated_output"]][["upstream_inference"]] <<- NULL
  hcobject[["integrated_output"]][["knowledge_network"]] <<- NULL

  sat <- hcobject[["satellite_outputs"]]
  if (is.null(sat) || !is.list(sat)) {
    sat <- list()
  }

  sat[["enrichments"]] <- NULL
  sat[["upstream_inference"]] <- NULL
  sat[["knowledge_network"]] <- NULL

  sat_names <- names(sat)
  if (!is.null(sat_names) && length(sat_names) > 0) {
    drop_nm <- sat_names[grepl("^celltype_annotation", sat_names)]
    if (length(drop_nm) > 0) {
      for (nm in drop_nm) {
        sat[[nm]] <- NULL
      }
    }
  }

  hcobject[["satellite_outputs"]] <<- sat
}
