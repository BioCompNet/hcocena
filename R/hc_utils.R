#' Function To Provide Cluster Colours
#' @noRd

get_cluster_colours <- function() {
  # NOTE: colours must stay unique. Modules are coloured by index here, but
  # `merge_clusters()` rebuilds its module table keyed by colour, so a repeated
  # colour would silently collapse two distinct groups into one. The second
  # "slategray" (position 39) was replaced with "purple" to keep the palette
  # unique without shifting the indices of earlier colours. `unique()` guards
  # against accidental duplicates being reintroduced later.
  col_vec <- c(
    "coral", "gold", "steelblue", "lightgreen", "turquoise", "plum", "maroon", "seagreen", "wheat", "slategray", "lightblue",
    "orchid", "darkgreen", "darkorange", "darkgrey", "indianred", "pink", "sandybrown", "khaki", "darkblue", "cadetblue",
    "greenyellow", "cyan", "thistle", "darkmagenta", "red", "blue", "green", "yellow", "brown", "black", "darkgoldenrod",
    "cornsilk", "firebrick", "deeppink", "dodgerblue", "lightpink", "midnightblue", "purple", "aquamarine", "chocolate",
    "darkred", "navy", "olivedrab", "peachpuff", "tomato", "snow"
  )
  return(base::unique(col_vec))
}

.hc_match_axis_indices_with_duplicates <- function(axis_ids, requested_order) {
  axis_ids <- base::as.character(axis_ids)
  requested_order <- base::as.character(requested_order)
  if (base::length(axis_ids) == 0) {
    return(base::integer())
  }
  if (base::length(requested_order) == 0) {
    return(base::seq_along(axis_ids))
  }

  source_idx <- base::split(base::seq_along(axis_ids), axis_ids)
  used <- stats::setNames(base::integer(base::length(source_idx)), base::names(source_idx))
  out_idx <- base::integer()

  for (val in requested_order) {
    idxs <- source_idx[[val]]
    if (base::is.null(idxs) || base::length(idxs) == 0) {
      next
    }
    next_pos <- used[[val]] + 1L
    if (next_pos <= base::length(idxs)) {
      out_idx <- base::c(out_idx, idxs[[next_pos]])
      used[[val]] <- next_pos
    }
  }

  remaining_idx <- base::seq_along(axis_ids)[!(base::seq_along(axis_ids) %in% out_idx)]
  base::c(out_idx, remaining_idx)
}

.hc_match_axis_order_with_duplicates <- function(axis_ids, requested_order) {
  axis_ids <- base::as.character(axis_ids)
  axis_ids[.hc_match_axis_indices_with_duplicates(axis_ids, requested_order)]
}

.hc_subset_matrix_cols_with_duplicates <- function(mat, requested_order = NULL) {
  if (base::is.null(mat)) {
    return(mat)
  }
  if (base::is.null(base::colnames(mat)) || base::ncol(mat) == 0) {
    return(mat)
  }
  idx <- .hc_match_axis_indices_with_duplicates(base::colnames(mat), requested_order)
  mat[, idx, drop = FALSE]
}

.hc_resolve_heatmap_col_order <- function(mat_cols,
                                          requested_order = NULL,
                                          context = "heatmap",
                                          warn_on_missing = TRUE) {
  mat_cols <- base::as.character(mat_cols)
  if (base::is.null(requested_order) || base::length(requested_order) == 0) {
    return(mat_cols)
  }

  requested_order <- base::as.character(requested_order)
  requested_order <- requested_order[!base::is.na(requested_order) & base::nzchar(base::trimws(requested_order))]
  if (base::length(requested_order) == 0) {
    return(mat_cols)
  }

  keep <- requested_order[requested_order %in% mat_cols]
  dropped <- base::setdiff(base::unique(requested_order), base::unique(mat_cols))
  if (isTRUE(warn_on_missing) && base::length(dropped) > 0) {
    preview <- base::paste(utils::head(dropped, 8L), collapse = ", ")
    if (base::length(dropped) > 8L) {
      preview <- base::paste0(preview, ", ...")
    }
    warning(
      "Ignoring ", base::length(dropped),
      " `col_order` entries not present in the current ", context, ": ",
      preview,
      call. = FALSE
    )
  }

  if (base::length(keep) == 0) {
    return(mat_cols)
  }

  .hc_match_axis_order_with_duplicates(mat_cols, keep)
}

.hc_resolve_col_order_alias <- function(col_order = NULL,
                                        heatmap_col_order = NULL,
                                        col_order_missing = FALSE,
                                        heatmap_col_order_missing = TRUE,
                                        context = "heatmap plot") {
  if (!base::is.null(col_order)) {
    col_order <- base::as.character(col_order)
  }
  alias_provided <- !isTRUE(heatmap_col_order_missing) && !base::is.null(heatmap_col_order)
  direct_provided <- !isTRUE(col_order_missing) && !base::is.null(col_order)

  if (!alias_provided) {
    return(col_order)
  }

  heatmap_col_order <- base::as.character(heatmap_col_order)
  if (direct_provided && !base::identical(col_order, heatmap_col_order)) {
    stop(
      "Use either `col_order` or legacy `heatmap_col_order` in ", context,
      ", not both with different values.",
      call. = FALSE
    )
  }

  if (direct_provided) {
    return(col_order)
  }
  heatmap_col_order
}

.hc_resolve_cluster_columns_alias <- function(cluster_columns = FALSE,
                                              heatmap_cluster_columns = NULL,
                                              cluster_columns_missing = FALSE,
                                              heatmap_cluster_columns_missing = TRUE,
                                              context = "heatmap plot") {
  alias_provided <- !isTRUE(heatmap_cluster_columns_missing) &&
    !base::is.null(heatmap_cluster_columns)
  direct_provided <- !isTRUE(cluster_columns_missing)

  if (!alias_provided) {
    return(cluster_columns)
  }

  if (direct_provided && !base::identical(cluster_columns, heatmap_cluster_columns)) {
    stop(
      "Use either `cluster_columns` or legacy `heatmap_cluster_columns` in ", context,
      ", not both with different values.",
      call. = FALSE
    )
  }

  if (direct_provided) {
    return(cluster_columns)
  }
  heatmap_cluster_columns
}


#' Leiden Clustering
#'
#' Applies the Leiden community detection to the network.
#' @param g The network, an igraph object.
#' @param num_it The number of iteration the algorithm is supposed to run.
#' @param resolution The resolution of the leiden clustering, higher values result in more clusters and vice versa (Default: 0.1).
#' @param partition_type Name of the partition type. Select from 'CPMVertexPartition', 'ModularityVertexPartition', 'RBConfigurationVertexPartition' and 'RBERVertexPartition' (Default: 'RBConfigurationVertexPartition').
#' @noRd

#' Community detection algorithms that return the same partition on every run
#'
#' For these, repeating the clustering `no_of_iterations` times only costs time:
#' every replicate is byte-identical, so the stability vote is a no-op.
#' @noRd

.hc_deterministic_cluster_algos <- function() {
  base::c("cluster_fast_greedy", "cluster_walktrap")
}


#' Align the community labels of one partition to a reference partition
#'
#' Community detection returns arbitrary label *names*: run twice and the same
#' group of genes may be called "3" once and "7" the next time. Comparing the
#' raw labels across iterations - as the stability vote used to do - therefore
#' measures label permutation rather than clustering instability, and discards
#' most genes even when the partitions agree almost perfectly.
#'
#' This relabels `x` onto `reference` by greedily matching the pair of
#' communities with the largest overlap, then the next largest among what is
#' left, and so on. Source communities with no counterpart keep a distinct
#' `unmatched_*` label so genuine disagreement still registers.
#'
#' @param x Vector of community labels to relabel.
#' @param reference Vector of community labels to align to (same length).
#' @return A character vector of `x` expressed in `reference`'s labels.
#' @noRd

.hc_match_partition_labels <- function(x, reference) {
  x <- base::as.character(x)
  reference <- base::as.character(reference)
  if (base::length(x) != base::length(reference) || base::length(x) == 0) {
    return(x)
  }

  overlap <- base::table(x, reference)
  if (base::length(overlap) == 0) {
    return(x)
  }

  map <- stats::setNames(
    base::rep(NA_character_, base::nrow(overlap)),
    base::rownames(overlap)
  )
  remaining <- overlap
  while (base::any(remaining > 0)) {
    idx <- base::which(remaining == base::max(remaining), arr.ind = TRUE)
    r <- idx[1L, 1L]
    cc <- idx[1L, 2L]
    map[[base::rownames(remaining)[[r]]]] <- base::colnames(remaining)[[cc]]
    remaining[r, ] <- 0L
    remaining[, cc] <- 0L
  }

  unmatched <- base::is.na(map)
  if (base::any(unmatched)) {
    map[unmatched] <- base::paste0("unmatched_", base::names(map)[unmatched])
  }

  base::unname(map[x])
}


#' Leiden partition as a membership vector
#'
#' Runs the Leiden algorithm and returns a minimal `communities`-like list with
#' a 1-based `membership` vector named by vertex, so it can be used the same way
#' as the return value of the `igraph::cluster_*` functions.
#' @noRd

.hc_leiden_membership <- function(g, num_it, resolution, partition_type,
                                  seed = 168575L) {
  tmp.partition <- leidenbase::leiden_find_partition(
    igraph = g,
    partition_type = partition_type,
    edge_weights = igraph::E(g)$weight,
    resolution = resolution,
    num_iter = num_it,
    seed = seed
  )

  # leidenbase already enumerates communities from 1 (the comment in the old
  # code claiming 0-based enumeration was wrong). Re-index through a factor so
  # the ids stay a gapless 1..K run even if the backend skips a number.
  membership <- base::as.integer(base::as.factor(tmp.partition$membership))
  base::names(membership) <- igraph::V(g)$name

  base::list(
    membership = membership,
    algorithm = "leiden",
    resolution = resolution,
    n.iter = num_it,
    names = igraph::V(g)$name
  )
}


leiden_clustering <- function(g, num_it, resolution, partition_type) {
  color.cluster <- get_cluster_colours()

  # run Leiden algorithm ion network:
  partition <- .hc_leiden_membership(
    g = g, num_it = num_it, resolution = resolution,
    partition_type = partition_type
  )

  # extract found clusters and the genes belonging to them.
  # NB: `.hc_leiden_membership()` returns gapless 1-based ids. The previous code
  # added 1 on top of the already 1-based factor codes, so module ids started at
  # 2 and the first palette colour (coral) was never used.
  clusters_df <- base::data.frame(
    cluster = base::as.integer(partition$membership),
    gene = partition$names,
    stringsAsFactors = FALSE
  )

  # get gene counts per cluster:
  cluster_frequencies <- base::table(clusters_df$cluster) %>% base::as.data.frame()

  # detect clusters large enough to be kept:
  clusters_to_keep <- dplyr::filter(cluster_frequencies, Freq >= hcobject[["global_settings"]][["min_nodes_number_for_cluster"]]) %>%
    dplyr::pull(., "Var1")

  # define white clusters (those that are too small to be kept):
  clusters_df_white <- dplyr::filter(clusters_df, !cluster %in% clusters_to_keep)

  # remove white clusters:
  clusters_df <- dplyr::filter(clusters_df, cluster %in% clusters_to_keep)

  # inform how many clusters and accordingly how many genes were lost due to insufficient cluster size:
  message(
    base::length(base::unique(clusters_df_white$cluster)),
    " cluster/s was/were smaller than the set minimum cluster size and therefore discarded.",
    " This removes ",
    base::nrow(clusters_df_white),
    " genes."
  )


  out <- base::lapply(base::unique(clusters_df$cluster), function(x) {
    # extract genes present in current cluster:
    genes <- dplyr::filter(clusters_df, cluster == x) %>%
      dplyr::pull(., "gene")

    # cluster name:
    col_clusters <- base::paste0("cluster ", x)
    # number of genes in cluster:
    col_gene_no <- dplyr::filter(clusters_df, cluster == x) %>%
      base::nrow()
    # comma separated list of all genes in the cluster:
    col_gene_n <- genes %>% base::paste0(., collapse = ",")
    # is the cluster included in the network (i.e., is it large enough):
    col_cluster_included <- "yes"
    # cluster color:
    col_color <- color.cluster[x]
    # order of conditions for following GFC values:
    col_conditions <- .hc_gfc_condition_names(hcobject[["integrated_output"]][["GFC_all_layers"]]) %>%
      base::paste0(., collapse = "#")
    # mean GFCs of the cluster genes per sample group:
    col_grp_means <- .hc_gfc_colmeans_for_genes(
      hcobject[["integrated_output"]][["GFC_all_layers"]],
      genes = genes
    ) %>%
      base::round(., digits = 3) %>%
      base::paste0(., collapse = ",")
    # collect all information:
    out <- base::data.frame(
      clusters = col_clusters,
      gene_no = col_gene_no,
      gene_n = col_gene_n,
      cluster_included = col_cluster_included,
      color = col_color,
      conditions = col_conditions,
      grp_means = col_grp_means,
      vertexsize = 3
    )
    return(out)
  }) %>% rlist::list.rbind()
  return(out)
}


#' Internal Function Used In cluster_calculation()
#' @noRd

cluster_calculation_internal <- function(graph_obj,
                                         algo,
                                         case,
                                         resolution,
                                         partition_type = "RBConfigurationVertexPartition",
                                         it = 2L,
                                         seed_offset = 0L) {
  # `it` used to default to the free variable `no_of_iterations`, which does not
  # exist in the package namespace. The default is only forced in the Leiden
  # branch, so `cluster_algo = "auto"` ran five algorithms and then died with
  # "object 'no_of_iterations' not found" on cluster_leiden.
  if (algo == "cluster_leiden") {
    # leiden_clustering() returns the finished module *table*, which has no
    # $membership -- modularity() and case = "best" both need a membership
    # vector, so derive one here instead.
    #
    # `seed_offset` makes replicate runs differ. The caller used to vary the
    # *number of Leiden iterations* per replicate instead, which produced runs
    # under systematically different algorithm settings rather than repeats of
    # the same one. Varying the seed keeps replicates comparable and the whole
    # thing reproducible.
    cfg <- .hc_leiden_membership(
      g = graph_obj, num_it = it, resolution = resolution,
      partition_type = partition_type,
      seed = 168575L + base::as.integer(seed_offset)
    )
  } else {
    cfg <- base::getExportedValue("igraph", algo)(graph_obj)
  }

  mod_score <- igraph::modularity(graph_obj, base::as.numeric(cfg$membership))

  mod_df <- base::data.frame(
    modularity_score = mod_score,
    cluster_algorithm = algo,
    stringsAsFactors = FALSE
  )

  # making switch so that in the end when only the best algorithm is to be used then the same function can be used
  output <- base::switch(case,
    best = cfg$membership,
    test = mod_df,
    final = cfg
  )

  message(algo, " algorithm tested")
  return(output)
}

#' Resolve GFC value-column indices while preserving duplicate condition names
#' @noRd
.hc_gfc_value_col_idx <- function(gfc_df) {
  if (base::is.null(gfc_df) || !base::is.data.frame(gfc_df) || base::ncol(gfc_df) == 0) {
    return(base::integer())
  }
  gene_idx <- base::which(base::colnames(gfc_df) %in% "Gene")
  if (base::length(gene_idx) == 0) {
    stop("`GFC_all_layers` must contain a `Gene` column.")
  }
  base::setdiff(base::seq_len(base::ncol(gfc_df)), gene_idx)
}

#' Return GFC condition names without triggering name repair on duplicates
#' @noRd
.hc_gfc_condition_names <- function(gfc_df) {
  idx <- .hc_gfc_value_col_idx(gfc_df)
  if (base::length(idx) == 0) {
    return(base::character())
  }
  base::colnames(gfc_df)[idx]
}

#' Extract the numeric GFC value frame without the `Gene` column
#' @noRd
.hc_gfc_value_frame <- function(gfc_df) {
  idx <- .hc_gfc_value_col_idx(gfc_df)
  if (base::length(idx) == 0) {
    return(base::data.frame())
  }
  out <- gfc_df[, idx, drop = FALSE]
  if (inherits(out, "DataFrame")) {
    return(.hc_to_base_data_frame_preserve_names(out))
  }
  base::data.frame(base::lapply(out, base::identity), check.names = FALSE)
}

#' Compute mean GFC values for a gene set while preserving duplicate condition names
#' @noRd
.hc_gfc_colmeans_for_genes <- function(gfc_df, genes) {
  if (base::is.null(gfc_df) || !base::is.data.frame(gfc_df) || !"Gene" %in% base::colnames(gfc_df)) {
    stop("`GFC_all_layers` must be a data.frame with a `Gene` column.")
  }
  cond_names <- .hc_gfc_condition_names(gfc_df)
  sub_df <- gfc_df[gfc_df[["Gene"]] %in% genes, , drop = FALSE]
  val_df <- .hc_gfc_value_frame(sub_df)
  if (base::ncol(val_df) == 0) {
    return(stats::setNames(base::numeric(), base::character()))
  }
  if (base::nrow(val_df) == 0) {
    return(stats::setNames(base::rep(NA_real_, base::ncol(val_df)), cond_names))
  }
  val_mat <- base::as.matrix(val_df)
  storage.mode(val_mat) <- "numeric"
  stats::setNames(base::colMeans(val_mat, na.rm = TRUE), cond_names)
}

.hc_group_values_from_annotation_for_gfc <- function(anno_df, voi = NULL) {
  if (base::is.null(anno_df) || !base::is.data.frame(anno_df) || base::nrow(anno_df) == 0) {
    return(base::character())
  }

  candidate_cols <- base::intersect(base::as.character(voi), base::colnames(anno_df))
  grp <- if (base::length(candidate_cols) > 1) {
    do.call(base::paste, base::c(anno_df[, candidate_cols, drop = FALSE], sep = "-"))
  } else if (base::length(candidate_cols) == 1) {
    anno_df[[candidate_cols[[1]]]]
  } else {
    anno_df[[1]]
  }

  grp <- base::trimws(base::as.character(grp))
  grp[grp %in% c("", "NA", "<NA>", "[NA]", "[<NA>]")] <- NA_character_
  grp[!base::is.na(grp) & base::nzchar(grp)]
}

.hc_layer_name_map <- function(hcobject) {
  layer_ids <- base::names(hcobject[["layers"]])
  if (base::is.null(layer_ids) || base::length(layer_ids) == 0) {
    anno_keys <- base::grep("_anno$", base::names(hcobject[["data"]]), value = TRUE)
    layer_ids <- base::sub("_anno$", "", anno_keys)
  }
  if (base::is.null(layer_ids) || base::length(layer_ids) == 0) {
    return(stats::setNames(base::character(), base::character()))
  }

  layer_names <- hcobject[["layers_names"]]
  if (base::is.null(layer_names) || base::length(layer_names) != base::length(layer_ids)) {
    layer_names <- layer_ids
  }
  stats::setNames(base::as.character(layer_names), base::as.character(layer_ids))
}

.hc_gfc_layer_source_rows <- function(hcobject) {
  layer_map <- .hc_layer_name_map(hcobject)
  layer_ids <- base::names(layer_map)
  if (base::length(layer_ids) == 0) {
    return(base::data.frame())
  }

  voi <- tryCatch(hcobject[["global_settings"]][["voi"]], error = function(e) NULL)
  layer_specific <- hcobject[["layer_specific_outputs"]]
  out_rows <- base::vector("list", base::length(layer_ids))

  for (i in base::seq_along(layer_ids)) {
    lid <- layer_ids[[i]]
    layer_label <- layer_map[[lid]]
    anno_df <- hcobject[["data"]][[base::paste0(lid, "_anno")]]

    gfc_layer <- NULL
    if (!base::is.null(layer_specific) && base::length(layer_specific) > 0) {
      gfc_layer <- tryCatch(layer_specific[[lid]][["part2"]][["GFC_all_genes"]], error = function(e) NULL)
      if (base::is.null(gfc_layer) && base::length(layer_specific) >= i) {
        gfc_layer <- tryCatch(layer_specific[[i]][["part2"]][["GFC_all_genes"]], error = function(e) NULL)
      }
    }

    cond_cols <- if (!base::is.null(gfc_layer) && base::is.data.frame(gfc_layer)) {
      .hc_gfc_condition_names(gfc_layer)
    } else {
      base::sort(base::unique(.hc_group_values_from_annotation_for_gfc(anno_df, voi = voi)))
    }
    cond_cols <- base::as.character(cond_cols)
    cond_cols <- cond_cols[!base::is.na(cond_cols) & base::nzchar(cond_cols)]

    grp_vals <- .hc_group_values_from_annotation_for_gfc(anno_df, voi = voi)
    grp_tbl <- if (base::length(grp_vals) > 0) base::table(grp_vals) else base::integer()
    sample_count <- .hc_as_integer_safely(grp_tbl[cond_cols])
    sample_count[base::is.na(sample_count)] <- 0L

    out_rows[[i]] <- base::data.frame(
      raw_condition = cond_cols,
      layer_id = lid,
      layer_name = layer_label,
      sample_count = sample_count,
      stringsAsFactors = FALSE
    )
  }

  out_rows <- out_rows[base::vapply(out_rows, function(x) base::is.data.frame(x) && base::nrow(x) > 0, FUN.VALUE = base::logical(1))]
  if (base::length(out_rows) == 0) {
    return(base::data.frame())
  }

  out <- base::do.call(base::rbind, out_rows)
  out$occurrence <- stats::ave(base::seq_len(base::nrow(out)), out$raw_condition, FUN = base::seq_along)
  base::rownames(out) <- NULL
  out
}

.hc_gfc_column_display_metadata <- function(hcobject, cols) {
  cols <- base::as.character(cols)
  out <- base::data.frame(
    raw_condition = cols,
    display_label = cols,
    count_label = cols,
    layer_id = NA_character_,
    layer_name = NA_character_,
    sample_count = NA_integer_,
    stringsAsFactors = FALSE
  )
  if (base::length(cols) == 0) {
    return(out)
  }

  dup_raw <- base::duplicated(cols) | base::duplicated(cols, fromLast = TRUE)
  source <- .hc_gfc_layer_source_rows(hcobject)
  if (base::nrow(source) == 0) {
    return(out)
  }

  if (base::nrow(source) == base::length(cols) &&
    base::all(base::as.character(source$raw_condition) == cols)) {
    mapped <- source
  } else {
    occ <- stats::ave(base::seq_along(cols), cols, FUN = base::seq_along)
    mapped_idx <- base::vapply(
      base::seq_along(cols),
      function(i) {
        hit <- base::which(source$raw_condition == cols[[i]] & source$occurrence == occ[[i]])
        if (base::length(hit) > 0) hit[[1]] else NA_integer_
      },
      FUN.VALUE = base::integer(1)
    )
    mapped <- source[mapped_idx, , drop = FALSE]
  }

  has_map <- !base::is.na(mapped$layer_name) & base::nzchar(mapped$layer_name)
  display <- cols
  display[dup_raw & has_map] <- base::paste0(mapped$layer_name[dup_raw & has_map], ": ", cols[dup_raw & has_map])

  count_label <- display
  has_count <- !base::is.na(mapped$sample_count) & mapped$sample_count > 0L
  count_label[has_count] <- base::paste0(display[has_count], "  [", mapped$sample_count[has_count], "]")

  out$display_label <- display
  out$count_label <- count_label
  out$layer_id <- mapped$layer_id
  out$layer_name <- mapped$layer_name
  out$sample_count <- mapped$sample_count
  out
}

.hc_gfc_display_col_labels <- function(hcobject, cols) {
  .hc_gfc_column_display_metadata(hcobject, cols)$display_label
}

.hc_gfc_display_count_labels <- function(hcobject, cols) {
  .hc_gfc_column_display_metadata(hcobject, cols)$count_label
}

.hc_gfc_duplicate_condition_width_scale <- function(hcobject, cols) {
  meta <- .hc_gfc_column_display_metadata(hcobject, cols)
  if (!base::is.data.frame(meta) || base::nrow(meta) <= 1) {
    return(1)
  }

  raw_labels <- base::as.character(meta$raw_condition)
  display_labels <- base::as.character(meta$display_label)
  keep <- !base::is.na(raw_labels) & base::nzchar(raw_labels)
  raw_labels <- raw_labels[keep]
  display_labels <- display_labels[keep]
  if (base::length(raw_labels) <= 1) {
    return(1)
  }

  has_layer_prefixed_duplicates <- base::any(
    !base::is.na(display_labels) &
      base::nzchar(display_labels) &
      display_labels != raw_labels
  )
  if (!isTRUE(has_layer_prefixed_duplicates)) {
    return(1)
  }

  unique_n <- base::length(base::unique(raw_labels))
  total_n <- base::length(raw_labels)
  if (unique_n <= 0 || unique_n >= total_n) {
    return(1)
  }

  dup_factor <- total_n / unique_n
  width_scale <- 1 + (0.12 * (dup_factor - 1))
  width_scale <- base::max(1, base::min(1.25, width_scale))
  as.numeric(width_scale[[1]])
}

.hc_parse_heatmap_col_layer_suffix <- function(hcobject, cols) {
  cols <- base::as.character(cols)
  if (base::length(cols) == 0) {
    return(NULL)
  }

  layer_names <- base::unique(base::unname(.hc_layer_name_map(hcobject)))
  layer_names <- base::as.character(layer_names)
  layer_names <- layer_names[!base::is.na(layer_names) & base::nzchar(layer_names)]
  if (base::length(layer_names) == 0) {
    return(NULL)
  }

  layer_names <- layer_names[base::order(base::nchar(layer_names), decreasing = TRUE)]
  parsed_layer <- base::rep(NA_character_, base::length(cols))
  parsed_prefix <- base::rep(NA_character_, base::length(cols))

  for (i in base::seq_along(cols)) {
    current_col <- cols[[i]]
    suffix_hits <- layer_names[base::endsWith(current_col, base::paste0("_", layer_names))]
    if (base::length(suffix_hits) == 0) {
      next
    }
    matched_layer <- suffix_hits[[1]]
    suffix_txt <- base::paste0("_", matched_layer)
    prefix_txt <- base::substr(
      current_col,
      1,
      base::nchar(current_col) - base::nchar(suffix_txt)
    )
    parsed_layer[[i]] <- matched_layer
    parsed_prefix[[i]] <- if (base::nzchar(prefix_txt)) prefix_txt else NA_character_
  }

  if (!base::any(!base::is.na(parsed_layer) & base::nzchar(parsed_layer))) {
    return(NULL)
  }

  list(
    layer_name = parsed_layer,
    prefix = parsed_prefix
  )
}

.hc_heatmap_column_metadata_values <- function(hcobject,
                                               cols,
                                               metadata_column) {
  if (base::is.null(metadata_column) || base::length(metadata_column) == 0) {
    stop("`column_gap_by` must be NULL or a non-empty metadata column name.")
  }
  metadata_column <- base::trimws(base::as.character(metadata_column[[1]]))
  if (base::is.na(metadata_column) || !base::nzchar(metadata_column)) {
    stop("`column_gap_by` must be NULL or a non-empty metadata column name.")
  }

  cols <- base::as.character(cols)
  if (base::length(cols) == 0) {
    return(base::character())
  }

  clean_chr <- function(x) {
    x <- base::trimws(base::as.character(x))
    x[x %in% c("", "NA", "<NA>", "[NA]", "[<NA>]")] <- NA_character_
    x
  }
  group_vector <- function(anno_df, voi) {
    candidate_cols <- base::intersect(base::as.character(voi), base::colnames(anno_df))
    grp <- if (base::length(candidate_cols) > 1) {
      do.call(base::paste, base::c(anno_df[, candidate_cols, drop = FALSE], sep = "-"))
    } else if (base::length(candidate_cols) == 1) {
      anno_df[[candidate_cols[[1]]]]
    } else {
      anno_df[[1]]
    }
    clean_chr(grp)
  }

  layer_map <- .hc_layer_name_map(hcobject)
  layer_ids <- base::names(layer_map)
  if (base::length(layer_ids) == 0) {
    stop("`column_gap_by` cannot be used because no annotation layers were found.")
  }

  meta <- .hc_gfc_column_display_metadata(hcobject, cols)
  parsed_layer_suffix <- .hc_parse_heatmap_col_layer_suffix(hcobject, cols)
  voi <- tryCatch(hcobject[["global_settings"]][["voi"]], error = function(e) NULL)
  data_list <- hcobject[["data"]]
  out <- base::rep(NA_character_, base::length(cols))
  metadata_found_anywhere <- FALSE

  for (i in base::seq_along(cols)) {
    candidate_layer_ids <- base::character()
    mapped_layer_id <- base::as.character(meta$layer_id[[i]])
    if (!base::is.na(mapped_layer_id) && base::nzchar(mapped_layer_id)) {
      candidate_layer_ids <- mapped_layer_id
    }
    if (base::length(candidate_layer_ids) == 0 &&
      !base::is.null(parsed_layer_suffix)) {
      parsed_layer <- base::as.character(parsed_layer_suffix$layer_name[[i]])
      if (!base::is.na(parsed_layer) && base::nzchar(parsed_layer)) {
        candidate_layer_ids <- base::names(layer_map)[base::unname(layer_map) == parsed_layer]
      }
    }
    if (base::length(candidate_layer_ids) == 0) {
      candidate_layer_ids <- layer_ids
    }
    candidate_layer_ids <- base::unique(candidate_layer_ids)

    condition_candidates <- base::unique(clean_chr(base::c(
      cols[[i]],
      meta$raw_condition[[i]],
      if (!base::is.null(parsed_layer_suffix)) parsed_layer_suffix$prefix[[i]] else NA_character_
    )))
    condition_candidates <- condition_candidates[!base::is.na(condition_candidates)]

    values_i <- base::character()
    for (lid in candidate_layer_ids) {
      anno_df <- data_list[[base::paste0(lid, "_anno")]]
      if (base::is.null(anno_df)) {
        layer_pos <- base::match(lid, layer_ids)
        if (!base::is.na(layer_pos)) {
          anno_df <- data_list[[base::paste0("set", layer_pos, "_anno")]]
        }
      }
      if (base::is.null(anno_df) || !base::is.data.frame(anno_df)) {
        next
      }
      if (!(metadata_column %in% base::colnames(anno_df))) {
        next
      }
      metadata_found_anywhere <- TRUE
      grp <- group_vector(anno_df, voi)
      hit <- !base::is.na(grp) & grp %in% condition_candidates
      if (!base::any(hit)) {
        next
      }
      values_i <- base::c(values_i, clean_chr(anno_df[[metadata_column]][hit]))
    }

    values_i <- base::unique(values_i[!base::is.na(values_i)])
    if (base::length(values_i) == 1) {
      out[[i]] <- values_i[[1]]
    } else if (base::length(values_i) > 1) {
      stop(
        "`column_gap_by = \"", metadata_column, "\"` maps heatmap column `",
        cols[[i]], "` to multiple metadata values: ",
        base::paste(values_i, collapse = ", "),
        ". Use a metadata column that is constant within each heatmap column.",
        call. = FALSE
      )
    }
  }

  if (!isTRUE(metadata_found_anywhere)) {
    stop(
      "`column_gap_by = \"", metadata_column,
      "\"` was not found in any annotation table.",
      call. = FALSE
    )
  }
  if (base::any(base::is.na(out))) {
    missing_cols <- cols[base::is.na(out)]
    stop(
      "`column_gap_by = \"", metadata_column,
      "\"` could not be resolved for heatmap column(s): ",
      base::paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  out
}

.hc_heatmap_column_gap_spec <- function(hcobject,
                                        cols,
                                        cluster_columns = FALSE,
                                        gap_mm = 0.6,
                                        enabled = FALSE,
                                        metadata_column = NULL) {
  cols <- base::as.character(cols)
  metadata_column <- if (base::is.null(metadata_column)) {
    NULL
  } else if (base::length(metadata_column) == 0) {
    NULL
  } else {
    base::trimws(base::as.character(metadata_column[[1]]))
  }
  if (!base::is.null(metadata_column) &&
    (base::is.na(metadata_column) || !base::nzchar(metadata_column))) {
    metadata_column <- NULL
  }
  metadata_split_requested <- !base::is.null(metadata_column)
  empty_out <- list(
    column_split = NULL,
    column_gap = NULL,
    total_gap_mm = 0,
    source = NULL,
    slice_count = 1L,
    slice_titles = NULL
  )
  if (!isTRUE(enabled) && !isTRUE(metadata_split_requested)) {
    return(empty_out)
  }
  if (base::length(cols) <= 1 || isTRUE(cluster_columns)) {
    return(empty_out)
  }
  if (base::length(cols) <= 3 && !isTRUE(metadata_split_requested)) {
    return(empty_out)
  }

  add_candidate <- function(store,
                            keys,
                            source,
                            priority,
                            show_titles = FALSE,
                            allow_singleton_runs = FALSE) {
    keys <- base::as.character(keys)
    if (base::length(keys) != base::length(cols)) {
      return(store)
    }
    keys[base::is.na(keys) | !base::nzchar(keys)] <- NA_character_
    if (base::any(base::is.na(keys))) {
      return(store)
    }
    runs <- base::rle(keys)
    if (base::length(runs$lengths) <= 1) {
      return(store)
    }
    if (base::all(runs$lengths == 1) && !isTRUE(allow_singleton_runs)) {
      return(store)
    }
    if (base::length(base::unique(keys)) <= 1) {
      return(store)
    }

    store[[base::length(store) + 1]] <- list(
      keys = keys,
      source = source,
      priority = as.integer(priority[[1]]),
      max_run = base::max(runs$lengths),
      mean_run = base::mean(runs$lengths),
      covered_cols = base::sum(runs$lengths[runs$lengths > 1]),
      n_runs = base::length(runs$lengths),
      show_titles = isTRUE(show_titles)
    )
    store
  }

  candidates <- list()
  if (isTRUE(metadata_split_requested)) {
    metadata_keys <- .hc_heatmap_column_metadata_values(
      hcobject = hcobject,
      cols = cols,
      metadata_column = metadata_column
    )
    candidates <- add_candidate(
      candidates,
      metadata_keys,
      base::paste0("metadata:", metadata_column),
      0L,
      show_titles = TRUE,
      allow_singleton_runs = TRUE
    )
  } else {
    meta <- .hc_gfc_column_display_metadata(hcobject, cols)
    parsed_layer_suffix <- .hc_parse_heatmap_col_layer_suffix(hcobject, cols)

    if (!base::is.null(parsed_layer_suffix)) {
      candidates <- add_candidate(candidates, parsed_layer_suffix$prefix, "prefix_before_layer", 1L)
      candidates <- add_candidate(candidates, parsed_layer_suffix$layer_name, "layer_suffix", 4L, show_titles = TRUE)
    }
    if (base::is.data.frame(meta) && base::nrow(meta) == base::length(cols)) {
      candidates <- add_candidate(candidates, meta$raw_condition, "raw_condition", 2L)
      candidates <- add_candidate(candidates, meta$layer_name, "layer_name", 3L, show_titles = TRUE)
    }
    generic_prefix <- ifelse(base::grepl("_", cols), base::sub("_[^_]+$", "", cols), NA_character_)
    candidates <- add_candidate(candidates, generic_prefix, "prefix_before_last_underscore", 5L)
  }

  if (base::length(candidates) == 0) {
    return(empty_out)
  }

  ordering <- base::order(
    -base::vapply(candidates, `[[`, numeric(1), "max_run"),
    -base::vapply(candidates, `[[`, numeric(1), "mean_run"),
    -base::vapply(candidates, `[[`, numeric(1), "covered_cols"),
    base::vapply(candidates, `[[`, integer(1), "n_runs"),
    base::vapply(candidates, `[[`, integer(1), "priority")
  )
  best <- candidates[[ordering[[1]]]]
  runs <- base::rle(best$keys)
  slice_count <- base::length(runs$lengths)
  if (slice_count <= 1) {
    return(empty_out)
  }

  slice_titles <- base::as.character(runs$values)
  slice_titles[base::is.na(slice_titles) | !base::nzchar(slice_titles)] <- base::paste0("Part ", base::seq_len(slice_count))
  if (isTRUE(best$show_titles)) {
    split_values <- base::make.unique(slice_titles, sep = " ")
    split_ids <- base::inverse.rle(list(
      values = split_values,
      lengths = runs$lengths
    ))
    split_ids <- base::factor(split_ids, levels = split_values)
  } else {
    slice_titles[] <- ""
    split_ids <- base::inverse.rle(list(
      values = base::seq_len(slice_count),
      lengths = runs$lengths
    ))
    split_ids <- base::factor(split_ids, levels = base::seq_len(slice_count))
  }
  gap_mm_use <- as.numeric(gap_mm[[1]])
  if (!base::is.finite(gap_mm_use) || gap_mm_use <= 0) {
    gap_mm_use <- 0.6
  }
  gap_mm_use <- base::min(1.2, base::max(0.35, gap_mm_use))

  list(
    column_split = split_ids,
    column_gap = grid::unit(base::rep(gap_mm_use, slice_count - 1L), "mm"),
    total_gap_mm = (slice_count - 1L) * gap_mm_use,
    source = best$source,
    slice_count = as.integer(slice_count),
    slice_titles = slice_titles
  )
}

.hc_heatmap_add_column_gap_args <- function(hm_args,
                                            column_gap_spec,
                                            title_gp = NULL) {
  if (base::is.null(hm_args) || !base::is.list(hm_args)) {
    return(hm_args)
  }
  if (base::is.null(column_gap_spec) ||
    base::is.null(column_gap_spec$column_split) ||
    base::is.null(column_gap_spec$column_gap) ||
    base::length(column_gap_spec$column_split) == 0) {
    return(hm_args)
  }

  hm_args$column_split <- column_gap_spec$column_split
  hm_args$column_gap <- column_gap_spec$column_gap
  hm_args$cluster_column_slices <- FALSE
  hm_args$column_title <- column_gap_spec$slice_titles
  if (!base::is.null(title_gp) &&
    !base::is.null(column_gap_spec$slice_titles) &&
    base::any(base::nzchar(base::as.character(column_gap_spec$slice_titles)))) {
    hm_args$column_title_gp <- title_gp
  }
  hm_args
}

.hc_heatmap_ggplot_column_layout <- function(cols,
                                             column_gap_spec = NULL,
                                             default_cell_mm = 5) {
  cols <- base::as.character(cols)
  n_cols <- base::length(cols)
  if (n_cols == 0) {
    return(list(
      x = base::numeric(0),
      limits = c(0.5, 0.5),
      slice_df = base::data.frame(
        title = base::character(0),
        x = base::numeric(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  has_gap <- !base::is.null(column_gap_spec) &&
    !base::is.null(column_gap_spec$column_split) &&
    base::length(column_gap_spec$column_split) == n_cols
  if (!isTRUE(has_gap)) {
    x <- base::seq_len(n_cols)
    return(list(
      x = x,
      limits = c(0.5, n_cols + 0.5),
      slice_df = base::data.frame(
        title = base::character(0),
        x = base::numeric(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  split_chr <- base::as.character(column_gap_spec$column_split)
  runs <- base::rle(split_chr)
  if (base::length(runs$lengths) <= 1) {
    x <- base::seq_len(n_cols)
    return(list(
      x = x,
      limits = c(0.5, n_cols + 0.5),
      slice_df = base::data.frame(
        title = base::character(0),
        x = base::numeric(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  default_cell_mm <- .hc_first_numeric_value(default_cell_mm[[1]])
  if (!base::is.finite(default_cell_mm) || default_cell_mm <= 0) {
    default_cell_mm <- 5
  }
  total_gap_mm <- .hc_first_numeric_value(column_gap_spec$total_gap_mm[[1]])
  gap_mm_each <- if (base::is.finite(total_gap_mm) && (base::length(runs$lengths) > 1)) {
    total_gap_mm / (base::length(runs$lengths) - 1L)
  } else {
    0.6
  }
  gap_units <- base::max(0.08, base::min(0.24, gap_mm_each / default_cell_mm))

  x <- base::numeric(n_cols)
  slice_centers <- base::numeric(base::length(runs$lengths))
  idx_start <- 1L
  cursor <- 1
  for (i in base::seq_along(runs$lengths)) {
    len_i <- runs$lengths[[i]]
    idx_end <- idx_start + len_i - 1L
    pos_i <- cursor + base::seq.int(0, len_i - 1L)
    x[idx_start:idx_end] <- pos_i
    slice_centers[[i]] <- base::mean(pos_i)
    cursor <- base::max(pos_i) + 1 + if (i < base::length(runs$lengths)) gap_units else 0
    idx_start <- idx_end + 1L
  }

  slice_titles <- column_gap_spec$slice_titles
  if (base::is.null(slice_titles) || base::length(slice_titles) != base::length(runs$lengths)) {
    slice_titles <- base::as.character(runs$values)
  } else {
    slice_titles <- base::as.character(slice_titles)
  }
  keep_titles <- !base::is.na(slice_titles) & base::nzchar(slice_titles)

  list(
    x = x,
    limits = c(base::min(x) - 0.5, base::max(x) + 0.5),
    slice_df = base::data.frame(
      title = slice_titles[keep_titles],
      x = slice_centers[keep_titles],
      stringsAsFactors = FALSE
    )
  )
}

#' Internal Function Used In cluster_calculation()
#' @noRd

gfc_mean_clustergene <- function(rownum, cluster_df, gfc_dat) {
  d1 <- cluster_df[rownum, , drop = FALSE]
  gene_names <- d1["gene_n"] %>%
    stringi::stri_split_regex(pattern = ",") %>%
    base::unlist()

  gfc_means <- .hc_gfc_colmeans_for_genes(gfc_dat, genes = gene_names)

  d1$conditions <- base::paste0(base::names(gfc_means), collapse = "#")
  d1$grp_means <- base::paste0(base::round(gfc_means, 3), collapse = ",")
  return(d1)
}


#' Function That Reads Gene Expression Matrices
#'
#' Reads in the gene expression data from a matrix format. Files can be provided in most common file formats (.txt, .csv).
#'  Rows must correspond to genes, columns must correspond to samples.
#' @param file A string defining the file path.
#' @param rown A Boolean. Whether or not the file has rownames. Default is TRUE.
#' @param sep The separator of the file. Default is "\t" for tab separated files.
#' @param gene_symbol_col A String. Name of the column that contains the gene symbols.
#' @noRd

# function to read expression files:
read_expression_data <- function(file, rown = TRUE, sep = "\t", gene_symbol_col) {
  if (rown) {
    expression_data <- utils::read.table(
      file = file, row.names = 1,
      stringsAsFactors = FALSE, sep = sep, check.names = FALSE, header = TRUE, quote = ""
    )
  } else {
    expression_data <- utils::read.table(
      file = file,
      stringsAsFactors = FALSE, sep = sep, check.names = FALSE, header = TRUE, quote = ""
    )
  }
  expression_data <- make_rownames_unique(counts = expression_data, gene_symbol_col = gene_symbol_col)

  return(expression_data)
}


#' Remove Rows With Duplicate Rownames
#'
#' The function removes all duplicate genes and only keeps the first occurance. It also drops all non-numeric columns.
#' @param counts The count matrix.
#' @param gene_symbol_col A String. Name of the column that contains the gene symbols.
#' @noRd

make_rownames_unique <- function(counts, gene_symbol_col) {
  counts <- counts[!base::duplicated(counts[gene_symbol_col]), ] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., gene_symbol_col)

  # remove all non-numeric columns (description, gene id, etc.):
  for (x in base::colnames(counts)) {
    if (!base::is.numeric(counts[[x]])) {
      counts[[x]] <- NULL
    }
  }

  if (base::ncol(counts) == 0) {
    warning(
      "All columns were deleted when removing non-numeric columns. Please check the data type of your expression values.",
      call. = FALSE
    )
  }

  return(counts)
}


#' Function That Reads The Annotation File
#'
#' Reads in the annotation data from a matrix format. Files can be provided in most common file formats (.txt, .csv).
#'  Rows must correspond to samples, columns must correspond to meta information categories.
#'  Transforms all columns to factors.
#' @param file A string defining the file path.
#' @param rown A Boolean. Whether or not the file has rownames. Default is TRUE.
#' @param sep The separator of the file. Default is "\t" for tab separated files.
#' @param sample_col A String. Name of the column that contains the sample IDs.
#' @noRd

read_anno <- function(file, rown = TRUE, sep = "\t", sample_col) {
  if (rown) {
    anno <- utils::read.table(
      file = file, row.names = 1,
      stringsAsFactors = FALSE, sep = sep, check.names = TRUE, header = TRUE
    )
  } else {
    anno <- utils::read.table(
      file = file,
      stringsAsFactors = FALSE, sep = sep, check.names = TRUE, header = TRUE
    )
  }
  base::rownames(anno) <- dplyr::pull(anno, sample_col)

  anno[] <- base::lapply(anno, base::factor)

  return(anno)
}

#' Internal Implementation Of run_expression_analysis_1()
#'
#' @param x An Integer. Gives the dataset that is currently processed.
#' @param padj A String. Defines the method to be used for p-value adjustment. Valid values are "none" (default) or values for "method" in stats::p.adjust.
#' @param export A Boolean. If TRUE, correlation values and p-values will be exported.
#'  This can save time if you plan on re-running the analysis since computing pari-wise correlations is a bottleneck of the analysis. Default is FALSE.
#' @param import A list. Each slot represents a dataset and contains a vector of two file names, the first is the name of the correlation file of that dataset exported in previous runs.
#'  The second is the name of the p-value file exported in a previous run. For details, see the information file found in the repository. Default is NULL.
#' @param bayes Sanchez-Taltavull et al. (2016) suggest superiority of Bayesian correlation analysis to Pearson correlation in some cases.
#'  Therefore, the Pearson correlation values can be weighted with Bayesian correlation values. To do so, set the "bayes"-parameter to TRUE. Default is FALSE, using only Pearson correlations.
#' @param alpha A numeric value in `[0,1]`. Allows to adjust the strength of the Bayes weighting: For alpha = 0 the Pearson correlation values remain unaltered, for alpha = 1 the Pearson correlation value and the Bayesian correlation value contribute equally to the final correlation.
#' @param prior An integer, either 2 or 3, using prior 2 or 3 for the Bayes weighting as described in "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016).
#' @param corr_method Method for the correlation calculation if bayes is FALSE.
#'   Supported values are 'pearson' and 'spearman'.
#' @noRd


run_expression_analysis_1_body <- function(
  x,
  bayes,
  prior,
  alpha,
  padj,
  export,
  import,
  corr_method,
  corr_backend = "auto"
) {
  message("Currently processed dataset: ", hcobject[["layers_names"]][x])
  output <- list()

  # retrieve gene count matrix for current layer:
  count_table <- hcobject[["data"]][[base::paste0("set", x, "_counts")]]

  # retrieve meta information for current layer:
  anno_table <- hcobject[["data"]][[base::paste0("set", x, "_anno")]]

  # filter for top most variant genes if required:

  top_var <- hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]]
  if (!top_var == "all") {
    message("...extracting ", top_var, " top variant genes...")
  }

  # sort the count data based on the genes' decreasing variance:
  ds <- count_table[base::order(base::apply(count_table, 1, stats::var), decreasing = TRUE), ]

  output[["ds"]] <- ds

  # filter:
  if (!hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]] == "all") {
    dd2 <- utils::head(ds, hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]])
  } else {
    dd2 <- ds
  }

  output[["topvar"]] <- dd2
  dd2 <- t(dd2)


  # calculate pair-wise correlations:
  corr_calc_out <- pwcorr(
    dd2 = dd2,
    layer_set = hcobject[["layer_settings"]][[base::paste0("set", x)]],
    bayes = bayes,
    prior = prior,
    alpha = alpha,
    padj = padj,
    export = export,
    layer = x,
    import = import,
    corr_method = corr_method,
    corr_backend = corr_backend
  )

  output[["corr_calc_out"]] <- corr_calc_out

  # calculate cut-off statistics:
  message("...calculating cutoff statistics...")
  cutoff_stats <- .hc_cutoff_stats_fast(
    correlation_df_filt = corr_calc_out[["correlation_df_filt"]],
    range_cutoff = corr_calc_out[["range_cutoff"]],
    print.all.plots = hcobject[["layer_settings"]][[base::paste0("set", x)]][["print_distribution_plots"]],
    x = x,
    min_nodes = hcobject[["global_settings"]][["min_nodes_number_for_network"]]
  )


  output[["cutoff_stats"]] <- cutoff_stats

  # reshape cutoff stats:
  cutoff_calc_out <- reshape_cutoff_stats(cutoff_stats = cutoff_stats)

  output[["cutoff_calc_out"]] <- cutoff_calc_out


  return(output)
}


#' Calculate p-value
#'
#' Calculates the p-value ofd a given value based on a reference distribution.
#' @param x Value for which to calculate the p-value
#' @param mu The mean of the reference population
#' @param sigma The standard deviation of the reference population
#' @param n The size of the reference population
#' @noRd

calc_pval <- function(x, mu, sigma, n) {
  z <- (x - mu) / (sigma / base::sqrt(n))
  p <- stats::pnorm(-base::abs(z))
  return(p)
}


#' Fast drop-in replacement for `Hmisc::rcorr()`
#'
#' Computes a dense gene-by-gene correlation matrix together with the matching
#' two-sided p-value matrix, returning the same `list(r, P, n)` structure as
#' [Hmisc::rcorr()]. Correlations are obtained from the cross-product of the
#' column-standardised matrix (a BLAS matmul, which is ~13-25x faster than
#' `rcorr`'s single-threaded Fortran), and p-values from the analytic
#' t-distribution formula that `rcorr` itself uses, so results are numerically
#' identical to `rcorr` (to floating-point precision).
#'
#' The fast path assumes complete observations (constant `n` across all pairs).
#' If the input contains `NA`s, `rcorr`'s pairwise-complete semantics change the
#' effective `n` per pair, so `backend = "auto"` transparently falls back to
#' [Hmisc::rcorr()] to keep results exact. Use `backend = "rcorr"` to force the
#' original behaviour.
#'
#' @param x A samples-by-genes numeric matrix (correlations across columns), as
#'   passed to [Hmisc::rcorr()].
#' @param type Either "pearson" (default) or "spearman". Spearman ranks each
#'   column first, then applies the Pearson path, matching `rcorr`.
#' @param backend "auto" (fast cross-product path with automatic `rcorr`
#'   fallback on `NA`s) or "rcorr" (always use [Hmisc::rcorr()]).
#' @return A list with elements `r` (correlation matrix), `P` (p-value matrix,
#'   `NA` on the diagonal) and `n` (pairwise observation-count matrix).
#' @noRd

.hc_fast_rcorr <- function(x, type = "pearson", backend = "auto") {
  x <- base::as.matrix(x)
  backend <- base::match.arg(backend, c("auto", "rcorr"))
  if (!type %in% c("pearson", "spearman")) {
    stop("Parameter 'type' must be either 'pearson' or 'spearman'.", call. = FALSE)
  }

  has_na <- base::anyNA(x)
  if (backend == "rcorr" || has_na) {
    if (has_na && backend == "auto") {
      message("...NAs detected in expression matrix; using Hmisc::rcorr for exact pairwise-complete p-values...")
    }
    return(Hmisc::rcorr(x, type = type))
  }

  n <- base::nrow(x)
  if (n < 3L) {
    stop("Need at least 3 observations (samples) to compute correlation p-values.", call. = FALSE)
  }

  # Spearman: rank each column, then treat as Pearson (matches Hmisc::rcorr).
  if (type == "spearman") {
    x <- base::apply(x, 2, base::rank)
  }

  # Correlation via cross-product of column-standardised data (BLAS matmul).
  # `scale` already yields the per-column SDs, so reuse them to flag constant
  # columns instead of a second pass.
  xs <- base::scale(x)
  col_sd <- base::attr(xs, "scaled:scale")
  zero_var <- !base::is.finite(col_sd) | col_sd == 0
  r <- base::crossprod(xs) / (n - 1)
  r[r > 1] <- 1
  r[r < -1] <- -1

  # Analytic two-sided t p-values -- the exact formula Hmisc::rcorr uses.
  # (Constant columns give NaN here; they are overwritten with NA below.)
  df <- n - 2
  tstat <- base::suppressWarnings(r * base::sqrt(df / (1 - r^2)))
  P <- 2 * stats::pt(-base::abs(tstat), df)

  base::diag(r) <- 1
  base::diag(P) <- NA_real_

  # Constant columns have undefined correlation; match rcorr by returning NA.
  if (base::any(zero_var)) {
    r[zero_var, ] <- NA_real_
    r[, zero_var] <- NA_real_
    P[zero_var, ] <- NA_real_
    P[, zero_var] <- NA_real_
  }

  gene_names <- base::colnames(x)
  base::rownames(r) <- base::colnames(r) <- gene_names
  base::rownames(P) <- base::colnames(P) <- gene_names
  n_matrix <- base::matrix(
    base::as.integer(n),
    nrow = base::ncol(x),
    ncol = base::ncol(x),
    dimnames = base::list(gene_names, gene_names)
  )

  base::list(r = r, P = P, n = n_matrix)
}


#' Calculate Pair-Wise Correlations
#'
#' The function calculates the pair-wise correlations for all genes in each dataset and performs multiple-testing correction.
#'  Alternatively, previously calculated correlations and p-values can be imported. Pearson correlation coefficient may be weighted with Bayesian correlations.
#'  Correlations and p-values may be exported.
#' @param dd2 Transposed gene expression matrix.
#' @param padj A String. Defines the method to be used for p-value adjustment. Valid values are "none" (default) or values for "method" in stats::p.adjust.
#' @param export A Boolean. If TRUE, correlation values and p-values will be exported.
#'  This can save time if you plan on re-running the analysis since computing pari-wise correlations is a bottleneck of the analysis. Default is FALSE.
#' @param import A list. Each slot represents a dataset and contains a vector of two file names, the first is the name of the correlation file of that dataset exported in previous runs.
#'  The second is the name of the p-value file exported in a previous run. For details, see the information file found in the repository. Default is NULL.
#' @param bayes Sanchez-Taltavull et al. (2016) suggest superiority of Bayesian correlation analysis to Pearson correlation in some cases.
#'  Therefore, the Pearson correlation values can be weighted with Bayesian correlation values. To do so, set the "bayes"-parameter to TRUE. Default is FALSE, using only Pearson correlations.
#' @param alpha A numeric value in `[0,1]`. Allows to adjust the strength of the Bayes weighting: For alpha = 0 the Pearson correlation values remain unaltered, for alpha = 1 the Pearson correlation value and the Bayesian correlation value contribute equally to the final correlation.
#' @param prior An integer, either 2 or 3, using prior 2 or 3 for the Bayes weighting as described in "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016).
#' @param layer_set The layer specific settings for this layer.
#' @param layer An Integer indicating the currently processed dataset.
#' @param corr_backend Correlation backend: "auto" (fast cross-product path via
#'   [.hc_fast_rcorr()], with automatic `Hmisc::rcorr` fallback on `NA`s) or
#'   "rcorr" (always use `Hmisc::rcorr`). Default "auto".
#' @noRd

pwcorr <- function(
  dd2,
  layer_set,
  bayes,
  prior,
  alpha,
  padj,
  export,
  layer,
  import,
  corr_method,
  corr_backend = "auto"
) {
  message("...calculating pairwise correlations...")

  if (!base::is.character(corr_method) || base::length(corr_method) != 1L ||
    base::is.na(corr_method) || !corr_method %in% c("pearson", "spearman")) {
    stop("Parameter 'corr_method' must be either 'pearson' or 'spearman'.", call. = FALSE)
  }

  output <- list()

  # import of pre-calculated correlation values and their p-values:
  if (base::length(import) > 1) {
    if (!base::is.na(import[layer])) {
      # import matrix
      message("...importing correlation matrix from file...")
      correlation_matrix <- list()
      correlation_matrix[["r"]] <- utils::read.table(import[[layer]][1], header = TRUE, check.names = FALSE) %>% base::as.matrix()
      correlation_matrix[["P"]] <- utils::read.table(import[[layer]][2], header = TRUE, check.names = FALSE) %>% base::as.matrix()
      base::rownames(correlation_matrix[["r"]]) <- base::colnames(correlation_matrix[["r"]])
      base::rownames(correlation_matrix[["P"]]) <- base::colnames(correlation_matrix[["P"]])
    } else {
      correlation_matrix <- .hc_fast_rcorr(base::as.matrix(dd2), type = corr_method, backend = corr_backend)
    }
  } else if (base::length(import) == 1 & !base::is.null(import)) {
    # import matrix
    message("...importing correlation matrix from file...")
    correlation_matrix <- list()
    correlation_matrix[["r"]] <- utils::read.table(import[[layer]][1], header = TRUE, check.names = FALSE) %>% base::as.matrix()
    correlation_matrix[["P"]] <- utils::read.table(import[[layer]][2], header = TRUE, check.names = FALSE) %>% base::as.matrix()
    base::rownames(correlation_matrix[["r"]]) <- base::colnames(correlation_matrix[["r"]])
    base::rownames(correlation_matrix[["P"]]) <- base::colnames(correlation_matrix[["P"]])
  } else {
    correlation_matrix <- .hc_fast_rcorr(base::as.matrix(dd2), type = corr_method, backend = corr_backend)
  }

  # export correlations and p-values for future re-runs to avoid the bottleneck:
  if (export) {
    utils::write.table(
      x = correlation_matrix[["r"]],
      file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/correlation_matrix_", hcobject[["layers_names"]][layer], "_correlations.txt"),
      quote = FALSE, row.names = FALSE, col.names = TRUE, dec = "."
    )
    utils::write.table(
      x = correlation_matrix[["P"]],
      file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/correlation_matrix_", hcobject[["layers_names"]][layer], "_pvalues.txt"),
      quote = FALSE, row.names = FALSE, col.names = TRUE, dec = "."
    )
  }

  # Reshape the upper triangle without materialising a second full
  # `upper.tri()` matrix for the p-values.
  ind <- base::which(base::upper.tri(correlation_matrix[["r"]], diag = FALSE), arr.ind = TRUE)
  gene_names <- base::colnames(correlation_matrix[["r"]])
  correlation_df <- base::data.frame(
    V1 = gene_names[ind[, 1L]],
    V2 = gene_names[ind[, 2L]],
    rval = base::as.numeric(correlation_matrix[["r"]][ind]),
    pval = base::as.numeric(correlation_matrix[["P"]][ind]),
    stringsAsFactors = FALSE
  )

  # Bayes weighting:
  if (bayes) {
    correlation_df[["rval"]] <- bayes_weighting(dd2 = dd2, alpha = alpha, prior = prior, pearson_rval = correlation_df[["rval"]])
  }

  # multiple testing correction:
  if (!padj == "none") {
    message("...conducting multiple testing correction using method ", padj, "...")
    correlation_df[["pval"]] <- stats::p.adjust(correlation_df[["pval"]], method = padj)
  }


  # Retain finite, significant positive correlations. Undefined pairs can occur
  # for constant genes and must not turn into all-NA rows during subsetting.
  finite_pairs <- base::is.finite(correlation_df[["pval"]]) &
    base::is.finite(correlation_df[["rval"]])
  correlation_df_filt <- correlation_df[
    finite_pairs & correlation_df[["pval"]] < 0.05 & correlation_df[["rval"]] > 0,
    ,
    drop = FALSE
  ]


  # range of cutoff min to max (correlation)
  finite_correlations <- correlation_df[["rval"]][base::is.finite(correlation_df[["rval"]])]
  if (base::length(finite_correlations) == 0L) {
    stop(
      "No finite pairwise correlations could be calculated. ",
      "Check whether at least two genes have non-constant expression profiles.",
      call. = FALSE
    )
  }
  range_cutoff <- base::seq(
    from = layer_set[["min_corr"]], to = base::max(finite_correlations),
    length.out = layer_set[["range_cutoff_length"]]
  )
  range_cutoff <- base::round(range_cutoff, 3)
  if (base::length(range_cutoff) > layer_set[["range_cutoff_length"]]) {
    range_cutoff <- range_cutoff[base::seq_len(layer_set[["range_cutoff_length"]])]
  } else {
    range_cutoff <- range_cutoff
  }


  output[["corr_mat"]] <- correlation_matrix
  output[["correlation_df"]] <- correlation_df
  output[["correlation_df_filt"]] <- correlation_df_filt
  output[["range_cutoff"]] <- range_cutoff

  return(output)
}


#' Function That Weights Given Pearson Correlation Coefficients With Bayes Correlation Values
#'
#' @param dd2 Transposed gene expression matrix.
#' @param alpha A numeric value in `[0,1]`. Allows to adjust the strength of the Bayes weighting: For alpha = 0 the Pearson correlation values remain unaltered, for alpha = 1 the Pearson correlation value and the Bayesian correlation value contribute equally to the final correlation.
#' @param prior An integer, either 2 or 3, using prior 2 or 3 for the Bayes weighting as described in "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016).
#' @param pearson_rval A vector of Pearson Correlation Coefficients that are to be weighted.
#' @noRd

bayes_weighting <- function(dd2,
                            alpha,
                            prior,
                            pearson_rval) {
  # use prior 2 as described in "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016).
  if (prior == 2) {
    message("---weighting pearson with bayes using prior: 2---")
    # bayes:
    corr_b <- Bayes_Corr_Prior2(X = base::t(dd2))
    # use prior 3 as described in "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016).
  } else {
    message("---weighting pearson with bayes using prior: 3---")
    # bayes:
    corr_b <- Bayes_Corr_Prior3(X = base::t(dd2))
  }

  # reshape matrix:
  ind <- base::which(base::upper.tri(corr_b, diag = FALSE), arr.ind = TRUE)
  correlation_df2 <- base::cbind(ind, corr_b[ind]) %>% base::as.data.frame()
  correlation_df2[, 1] <- base::colnames(dd2)[correlation_df2[, 1]]
  correlation_df2[, 2] <- base::colnames(dd2)[correlation_df2[, 2]]
  base::colnames(correlation_df2) <- c("V1", "V2", "rval")
  weighted_pearson <- (pearson_rval + (alpha * correlation_df2[["rval"]])) / (1 + alpha)

  return(weighted_pearson)
}

#' Computing The Bayesian Correlations Assuming Second (Dirichlet-marginalized) Prior
#'
#' This function has been taken from  "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016), https://doi.org/10.1371/journal.pone.0163595.s001
#' @noRd

Bayes_Corr_Prior2 <- function(X) {
  d <- base::dim(X)
  alpha0 <- base::rep(1 / d[1], d[2])
  beta0 <- base::rep(1 - 1 / d[1], d[2])
  Bcorrvals <- Bayes_Corr(alpha0, beta0, X)
  return(Bcorrvals)
}


#' Computing The Bayesian Correlations Assuming Third (zero count-motivated) Prior
#'
#' This function has been taken from  "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016), https://doi.org/10.1371/journal.pone.0163595.s001
#' @noRd

Bayes_Corr_Prior3 <- function(X) {
  d <- base::dim(X)
  cs <- base::colSums(X)
  alpha0 <- (cs + 1) / (base::max(cs) + 1)
  beta0 <- base::rep(1, d[2])
  Bcorrvals <- Bayes_Corr(alpha0, beta0, X)
  return(Bcorrvals)
}

#' Function To Compute Bayesian Correlation Coefficients
#'
#' This function has been taken from  "Bayesian correlation analysis for sequence count data" by Sanchez-Taltavull et al. (2016), https://doi.org/10.1371/journal.pone.0163595.s001
#' @noRd

Bayes_Corr <- function(alpha0, beta0, X) {
  nrowsX <- base::nrow(X)
  k <- base::ncol(X)
  cs <- base::colSums(X)
  alphas <- base::matrix(base::rep(alpha0, nrowsX), nrow = nrowsX, byrow = TRUE) + X
  betas <- base::matrix(base::rep(beta0, nrowsX), nrow = nrowsX, byrow = TRUE) + base::matrix(base::rep(cs, nrowsX), nrow = nrowsX, byrow = TRUE) - X
  alphasPLUSbetas <- alphas + betas

  # First BIG product term for covariance formula
  Psi <- alphas / alphasPLUSbetas - base::matrix(base::rep(base::rowSums(alphas / alphasPLUSbetas) / k, k), ncol = k, byrow = FALSE)

  # Covariance matrix
  cov_mtrx <- Psi %*% base::t(Psi) / k

  # Variances (this is a column vector of length = nrowsX)
  var_vec <- base::as.matrix((base::rowSums((alphas * betas) / ((alphasPLUSbetas^2) * (alphasPLUSbetas + 1))) + base::rowSums(Psi^2)) / k)

  Bcorrvals <- cov_mtrx / base::sqrt(var_vec %*% base::t(var_vec))
  base::diag(Bcorrvals) <- 1
  return(Bcorrvals)
}


#' Calculate Cutoff Statistics
#'
#' @param cutoff The minimum correlation fro which to filter the data.
#' @param corrdf_r Data frame holding the correlation coefficients.
#' @param print.all.plots Boolean. Whether or not to print the degree distribution plots for all cutoffs.
#' @param x An Integer. Gives the dataset that is currently processed.
#' @noRd

.hc_union_find_components <- function(from_idx, to_idx, n_vertices) {
  # Weighted union-find with in-place mutation. The `parent`/`rank` vectors live
  # in this frame and are updated via `<<-`, which R performs in place because
  # each vector is single-referenced. `find` only *reads* `parent` (walking to
  # the root), so no per-operation vector copies occur -- unlike the previous
  # implementation, which returned the whole `parent` vector from every find and
  # rebound it, making each edge O(n_vertices) instead of O(alpha). Union by rank
  # keeps trees shallow (O(log n) finds). This is the large-graph fallback used
  # when igraph construction would be too memory-hungry.
  parent <- seq_len(n_vertices)
  rank <- integer(n_vertices)

  # `find` is inlined as a read-only while-loop (rather than a closure) so that
  # `parent` keeps a single reference and `parent[x] <- y` mutates it in place.
  from_idx <- as.integer(from_idx)
  to_idx <- as.integer(to_idx)
  for (k in seq_along(from_idx)) {
    root_from <- from_idx[k]
    while (parent[root_from] != root_from) {
      root_from <- parent[root_from]
    }
    root_to <- to_idx[k]
    while (parent[root_to] != root_to) {
      root_to <- parent[root_to]
    }
    if (root_from == root_to) {
      next
    }
    if (rank[root_from] < rank[root_to]) {
      parent[root_from] <- root_to
    } else if (rank[root_from] > rank[root_to]) {
      parent[root_to] <- root_from
    } else {
      parent[root_to] <- root_from
      rank[root_from] <- rank[root_from] + 1L
    }
  }

  roots <- integer(n_vertices)
  for (i in seq_len(n_vertices)) {
    r <- i
    while (parent[r] != r) {
      r <- parent[r]
    }
    roots[i] <- r
  }
  comp_ids <- match(roots, unique(roots))
  list(
    membership = comp_ids,
    csize = tabulate(comp_ids, nbins = max(comp_ids))
  )
}

.hc_safe_deep_clone <- function(x,
                                context = "object",
                                max_bytes = NULL) {
  if (base::is.null(x)) {
    return(NULL)
  }

  max_bytes_num <- if (base::is.null(max_bytes) || base::length(max_bytes) == 0) {
    NA_real_
  } else {
    .hc_first_numeric_value(max_bytes[[1]])
  }
  if (base::length(max_bytes_num) == 1 && base::is.finite(max_bytes_num) && max_bytes_num > 0) {
    est_size <- tryCatch(
      as.numeric(utils::object.size(x)),
      error = function(e) NA_real_
    )
    if (base::is.finite(est_size) && est_size > max_bytes_num) {
      warning(
        "Skipping deep clone of ", context, "; estimated size exceeds the safe clone threshold. ",
        "Reusing the original object instead.",
        call. = FALSE
      )
      return(x)
    }
  }

  out <- tryCatch(
    base::unserialize(base::serialize(x, NULL)),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    warning(
      "Could not deep-clone ", context, "; reusing the original object instead. ",
      "Reason: ", conditionMessage(out),
      call. = FALSE
    )
    return(x)
  }

  out
}

.hc_resolve_panel_storage_mode <- function(store_panel_objects = c("auto", "always", "never"),
                                           n_databases = NA_integer_) {
  mode <- base::match.arg(store_panel_objects)
  if (!identical(mode, "auto")) {
    return(mode)
  }

  n_databases <- .hc_first_numeric_value(.hc_as_integer_safely(n_databases[[1]]))
  if (!base::is.finite(n_databases) || n_databases <= 0L) {
    return("never")
  }
  if (n_databases <= 1L) {
    return("always")
  }
  "never"
}

.hc_cutoff_component_summary <- function(graph_df, min_nodes) {
  v1 <- base::as.character(graph_df[["V1"]])
  v2 <- base::as.character(graph_df[["V2"]])
  vertices <- base::unique(base::c(v1, v2))
  .hc_cutoff_component_summary_idx(
    from_idx = match(v1, vertices),
    to_idx = match(v2, vertices),
    n_vertices = length(vertices),
    min_nodes = min_nodes
  )
}

# Integer-index core of the cutoff component summary. `from_idx`/`to_idx` must be
# 1-based vertex indices into a compact vertex set of size `n_vertices` (i.e.
# every vertex 1..n_vertices appears in at least one edge), so the `tabulate`
# degree count stays consistent. Splitting this out lets the fast cutoff loop map
# gene names to integers once and reuse them for every cutoff instead of
# re-deriving them from character vectors each time.
#
# Components come from the in-place weighted union-find (`.hc_union_find_components`),
# which is faster than building an igraph object per cutoff (the cutoff loop
# calls this once per cutoff, so the repeated graph construction dominated) and
# uses less memory. Its component partition is verified identical to
# `igraph::components()` in the tests.
.hc_cutoff_component_summary_idx <- function(from_idx, to_idx, n_vertices, min_nodes) {
  components <- .hc_union_find_components(
    from_idx = from_idx,
    to_idx = to_idx,
    n_vertices = n_vertices
  )

  keep_components <- components[["csize"]] >= min_nodes
  keep_vertices <- keep_components[components[["membership"]]]
  degrees_all <- tabulate(base::c(from_idx, to_idx), nbins = n_vertices)
  keep_edges <- keep_vertices[from_idx] & keep_vertices[to_idx]

  list(
    keep_vertices = keep_vertices,
    degrees = degrees_all[keep_vertices],
    num_nodes = base::sum(keep_vertices),
    num_edges = base::sum(keep_edges),
    num_networks = base::sum(keep_components)
  )
}

cutoff_prep <- function(cutoff, corrdf_r, print.all.plots, x, min_nodes = hcobject[["global_settings"]][["min_nodes_number_for_network"]]) {
  ### filter correlations above the cutoff
  filteredmatrix <- corrdf_r[corrdf_r[["rval"]] >= cutoff, ]

  ## create expected df with initial values
  output <- base::data.frame(
    R.squared = 0,
    degree = 0,
    Probs = 0,
    cutoff = cutoff,
    no_edges = 0,
    no_nodes = 0,
    no_of_networks = 0
  )


  rownums <- base::ifelse(base::nrow(filteredmatrix) == 0, "zero_rows", "gtzero")

  base::switch(rownums,
    zero_rows = {
      output
    },
    gtzero = {
      rsquaredfun(graph_df = filteredmatrix, cutoff = cutoff, print.all.plots = print.all.plots, min_nodes = min_nodes, x = x)
    }
  )
}


#' Vectorised cutoff-statistics loop
#'
#' Result-identical, faster replacement for
#' `do.call(rbind, lapply(range_cutoff, cutoff_prep, corrdf_r = ...))`.
#'
#' The original loop re-derives the vertex set for every cutoff by converting
#' millions of gene-name strings with `as.character()` + `unique()` + `match()`
#' and by re-subsetting a multi-million-row data.frame. This version maps gene
#' names to integer vertex ids **once**, sorts the edges by correlation once, and
#' for each cutoff processes the integer prefix of passing edges. Component
#' statistics are order- and label-invariant, so the output is identical to the
#' per-cutoff `cutoff_prep()` (verified in tests).
#' @noRd
.hc_cutoff_stats_fast <- function(correlation_df_filt,
                                  range_cutoff,
                                  print.all.plots = FALSE,
                                  x = NULL,
                                  min_nodes = hcobject[["global_settings"]][["min_nodes_number_for_network"]]) {
  zero_row <- function(cutoff) {
    base::data.frame(
      R.squared = 0,
      degree = 0,
      Probs = 0,
      cutoff = cutoff,
      no_edges = 0,
      no_nodes = 0,
      no_of_networks = 0
    )
  }

  if (base::is.null(correlation_df_filt) || base::nrow(correlation_df_filt) == 0) {
    return(base::do.call(base::rbind, base::lapply(range_cutoff, zero_row)))
  }

  valid_edges <- !base::is.na(correlation_df_filt[["V1"]]) &
    !base::is.na(correlation_df_filt[["V2"]]) &
    base::nzchar(base::as.character(correlation_df_filt[["V1"]])) &
    base::nzchar(base::as.character(correlation_df_filt[["V2"]])) &
    base::is.finite(base::as.numeric(correlation_df_filt[["rval"]]))
  correlation_df_filt <- correlation_df_filt[valid_edges, , drop = FALSE]
  if (base::nrow(correlation_df_filt) == 0L) {
    return(base::do.call(base::rbind, base::lapply(range_cutoff, zero_row)))
  }

  # Map gene names to integer vertex ids once, then sort edges by correlation so
  # the edges passing any cutoff are a prefix of the sorted arrays.
  v1 <- base::as.character(correlation_df_filt[["V1"]])
  v2 <- base::as.character(correlation_df_filt[["V2"]])
  vertices <- base::unique(base::c(v1, v2))
  from_all <- match(v1, vertices)
  to_all <- match(v2, vertices)
  rval <- base::as.numeric(correlation_df_filt[["rval"]])

  ord <- base::order(rval, decreasing = TRUE)
  from_all <- from_all[ord]
  to_all <- to_all[ord]
  rval_sorted <- rval[ord]

  # `-rval_sorted` is ascending, so findInterval gives all edge-prefix sizes
  # without rescanning every edge for every cutoff.
  cutoff_values <- base::as.numeric(range_cutoff)
  k_by_cutoff <- base::integer(base::length(cutoff_values))
  valid_cutoffs <- !base::is.na(cutoff_values)
  k_by_cutoff[valid_cutoffs] <- base::findInterval(
    -cutoff_values[valid_cutoffs],
    -rval_sorted
  )
  summary_cache <- base::new.env(parent = base::emptyenv(), hash = TRUE)

  base::do.call(base::rbind, base::lapply(base::seq_along(cutoff_values), function(i) {
    cutoff <- cutoff_values[[i]]
    k <- k_by_cutoff[[i]]
    if (k == 0) {
      return(zero_row(cutoff))
    }
    cache_key <- base::as.character(k)
    if (base::exists(cache_key, envir = summary_cache, inherits = FALSE)) {
      stats_summary <- base::get(cache_key, envir = summary_cache, inherits = FALSE)
    } else {
      idx <- base::seq_len(k)
      from_k <- from_all[idx]
      to_k <- to_all[idx]
      # Compact the vertex ids to those present at this cutoff so the degree
      # tabulation remains consistent with cutoff_prep().
      present <- base::unique(base::c(from_k, to_k))
      stats_summary <- .hc_cutoff_component_summary_idx(
        from_idx = match(from_k, present),
        to_idx = match(to_k, present),
        n_vertices = base::length(present),
        min_nodes = min_nodes
      )
      base::assign(cache_key, stats_summary, envir = summary_cache)
    }
    .hc_rsquared_from_summary(
      stats_summary = stats_summary,
      cutoff = cutoff,
      print.all.plots = print.all.plots,
      x = x
    )
  }))
}


#' Function For Calculating Network Statistics
#'
#' @param graph_df Description of the network in matrix format with two columns giving edges between nodes and a third column giving edge weights.
#' @param cutoff The minimum correlation for which the edges are filtered.
#' @noRd


rsquaredfun <- function(graph_df, cutoff, print.all.plots, min_nodes = hcobject[["global_settings"]][["min_nodes_number_for_network"]], x = NULL) {
  stats_summary <- .hc_cutoff_component_summary(graph_df = graph_df, min_nodes = min_nodes)
  .hc_rsquared_from_summary(
    stats_summary = stats_summary,
    cutoff = cutoff,
    print.all.plots = print.all.plots,
    x = x
  )
}

# Scale-free-fit statistics for one cutoff from a precomputed component summary
# (see `.hc_cutoff_component_summary`). Identical to the tail of the original
# `rsquaredfun`; factored out so the fast cutoff loop can reuse it.
.hc_rsquared_from_summary <- function(stats_summary, cutoff, print.all.plots, x = NULL) {
  num_networks <- stats_summary[["num_networks"]]
  num_nodes <- stats_summary[["num_nodes"]]
  num_edges <- stats_summary[["num_edges"]]
  kept_degrees <- stats_summary[["degrees"]]


  ## calculate stats
  if (num_nodes == 0) {
    R.square <- NA
    degree <- NA
    probability <- NA
  } else {
    degree_counts <- tabulate(kept_degrees, nbins = base::max(kept_degrees))
    degree <- seq_len(base::length(degree_counts))
    probability <- degree_counts / base::sum(degree_counts)
    nonzero.position <- base::which(probability != 0)
    probability <- probability[nonzero.position]
    degree <- degree[nonzero.position]

    if (base::length(probability) == 0) {
      R.square <- 0
    } else {
      forplot <- base::data.frame(probability = probability, degree = degree)
      reg <- stats::lm(base::log(probability) ~ base::log(degree))
      cozf <- stats::coef(reg)
      power.law.fit <- function(x) base::exp(cozf[[1]] + cozf[[2]] * base::log(x))
      alpha <- -cozf[[2]]
      R.square <- base::summary(reg)$r.squared


      if (print.all.plots) {
        degree_distribution_wd <- base::paste0("dir_DegreeDistribution_", hcobject[["layer_settings"]][[base::paste0("set", x)]][["top_var"]])

        if (!degree_distribution_wd %in% base::list.dirs(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]]))) {
          base::dir.create(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], degree_distribution_wd))
        }

        dd_plot <- ggplot2::ggplot(forplot, ggplot2::aes(x = base::log(degree), y = base::log(probability))) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "lm") +
          ggplot2::theme_bw()

        .hc_export_ggplot_file(
          file = .hc_output_file(
            base::paste0("Degree_distribution_plot_", cutoff, "_set_", x, ".pdf"),
            degree_distribution_wd
          ),
          plot = dd_plot,
          width = 7,
          height = 7
        )
      }
    }
  }

  output <- base::data.frame(
    R.squared = R.square,
    degree = degree,
    Probs = probability,
    cutoff = cutoff,
    no_edges = num_edges,
    no_nodes = num_nodes,
    no_of_networks = num_networks
  )

  return(output)
}

####### ADD FUNCTIONS FROM optimal_cutoff_MultiOmics.R #######
#' Normalize criteria by column maxima
#'
#' Percentage-of-max normalization for criteria tables.
#' @param perf_tbl A numeric matrix/data.frame with criteria in columns.
#' @noRd
.normalize_by_max <- function(perf_tbl) {
  mat <- base::as.matrix(perf_tbl)
  mode(mat) <- "numeric"

  out <- base::apply(mat, 2, function(x) {
    if (base::all(base::is.na(x))) {
      return(base::rep(NA_real_, base::length(x)))
    }
    max_val <- base::max(x, na.rm = TRUE)
    if (!base::is.finite(max_val) || max_val == 0) {
      return(base::rep(0, base::length(x)))
    }
    x / max_val
  })

  if (base::is.null(base::dim(out))) {
    out <- base::matrix(out, ncol = 1)
    base::colnames(out) <- base::colnames(mat)
  }
  base::rownames(out) <- base::rownames(mat)
  base::as.data.frame(out, check.names = FALSE)
}

# Default GFC palette used across heatmaps.
# Matches the saved heatmap snapshot appearance.
.hc_default_gfc_colors <- function() {
  base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
}

#' Weighted sum over normalized criteria
#'
#' Weighted sum over normalized criteria.
#' @param perf_tbl A numeric matrix/data.frame with criteria in columns.
#' @param weights Named numeric vector with one weight per criterion.
#' @noRd
.weighted_sum <- function(perf_tbl, weights) {
  if (base::is.null(base::names(weights))) {
    stop("`weights` must be a named numeric vector.")
  }
  if (!base::all(base::names(weights) %in% base::colnames(perf_tbl))) {
    stop("`weights` names must be present in `perf_tbl` columns.")
  }

  mat <- base::as.matrix(perf_tbl[, base::names(weights), drop = FALSE])
  mode(mat) <- "numeric"
  ws <- mat %*% base::as.numeric(weights)
  ws <- base::as.vector(ws)
  base::names(ws) <- base::rownames(perf_tbl)
  ws
}

#' Reshape Cutoff Statistics
#'
#' The function reshapes the cutoff statistics for downstream usage.
#' @param cutoff_stats The native representation of the cutoff statistics.
#' @noRd

reshape_cutoff_stats <- function(cutoff_stats) {
  output <- list()

  cutoff_stats_concise <- cutoff_stats %>%
    dplyr::select(R.squared, cutoff, no_edges, no_nodes, no_of_networks) %>%
    dplyr::distinct()

  cutoff_stats_concise <- cutoff_stats_concise %>% dplyr::filter(no_of_networks != 0)
  base::rownames(cutoff_stats_concise) <- cutoff_stats_concise[["cutoff"]]
  cutoff_stats_concise <- cutoff_stats_concise[, -2]

  if (base::nrow(cutoff_stats_concise) == 0) {
    warning("cutoff_stats_concise is empty", call. = FALSE)
    return(base::data.frame())
  }
  nPT <- .normalize_by_max(cutoff_stats_concise[, c("R.squared", "no_edges", "no_nodes", "no_of_networks")])
  w <- base::c(0.5, 0.1, 0.5, -1)
  base::names(w) <- base::colnames(nPT)
  ws <- .weighted_sum(nPT, w)
  ranked_ws <- base::rank(-ws) %>% base::sort()


  calculated_optimal_cutoff <- base::as.numeric(base::names(ranked_ws[1]))


  stats_calculated_optimal_cutoff <- cutoff_stats[cutoff_stats[["cutoff"]] == calculated_optimal_cutoff, c("degree", "Probs")]

  dd_plot_calculated_optimal <- ggplot2::ggplot(stats_calculated_optimal_cutoff, ggplot2::aes(x = base::log(degree), y = base::log(Probs))) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(base::paste0("Calculated optimal correlation cut-off [", calculated_optimal_cutoff, "]"))


  output[["cutoff_stats_concise"]] <- cutoff_stats_concise
  output[["dd_plot_calculated_optimal"]] <- dd_plot_calculated_optimal
  output[["optimal_cutoff"]] <- calculated_optimal_cutoff
  return(output)
}

#' Internal Realization Of run_expression_analysis_2()
#'
#' Iterates over all datasets to perform the actions described in run_expression_analysis_2().
#' @param x An integer giving the number of the dataset currently processed.
#' @param grouping_v A string giving a column name present in all annotation files, if this variable shall be used for grouping the samles isntead of the variable of interest. Default is NULL.
#' @param plot_HM A Boolean. Whether or not to plot the heatmap (for networks with many genes this may be very demanding for your computer if you are running the analysis locally). Default is TRUE.
#' @param method The method used for clustering the heatmap in the pheatmap function. Default is "complete".
#' @param additional_anno A list, with one slot per data set. A slot contains a vector of column names from that data set's annotation file that you wish to annotate with.
#'  If for some of the data sets you don't wish any further annotation, you can set the corresponding list slot to NULL. Default is NULL.
#' @param cols A named list of color vectors. The list names need to match the chosen annotation column names. Default is NULL which uses implemented colors.
#' @noRd

run_expression_analysis_2_body <- function(x, grouping_v, plot_HM, method, additional_anno, title, cols) {
  message("...Currently processed dataset: ", hcobject[["layers_names"]][x], "...")

  .hc_set_bridge_hcobject_slot(
    c("layer_specific_outputs", base::paste0("set", x), "part2", "heatmap_out"),
    heatmap_network_genes(
      x = x,
      plot_HM = plot_HM,
      method = method,
      additional_anno = additional_anno,
      title = title,
      cols = cols
    )
  )


  .hc_set_bridge_hcobject_slot(
    c("layer_specific_outputs", base::paste0("set", x), "part2", "GFC_all_genes"),
    GFC_calculation(
      info_dataset = hcobject[["data"]][[base::paste0("set", x, "_anno")]],
      grouping_v = grouping_v,
      x = x
    )
  )
}


#' Plot Heatmap Of Network Genes
#'
#' Plots a heatmap of the genes present in the network, thus those genes left after filtering for the minimum correlation and removing too smal graph components.
#' @param x An integer giving the number of the dataset currently processed.
#' @param plot_HM A Boolean. Whether or not to plot the heatmap.
#' @param method The method used for clustering in the pheatmap function.
#' @param additional_anno A vector of strings giving other columns from the annotation to be annotated in the heatmap in addition to the variable of interest.
#' @param cols A named list of color vectors. The list names need to match the chosen annotation column names. Default is NULL which uses implemented colors.
#' @noRd

heatmap_network_genes <- function(x, plot_HM, method, additional_anno, title, cols) {
  # extract annotation data for current data layer:
  info_dataset <- hcobject[["data"]][[base::paste0("set", x, "_anno")]]

  message("...creating heatmap of network genes...")

  # collect function output:
  output <- list()

  # conditions to annotate in the heatmap:
  all_conditions <- base::unique(base::c(additional_anno, hcobject[["global_settings"]][["voi"]]))
  all_conditions <- all_conditions[!base::is.null(all_conditions)]

  # extract data exceeding set cutoff:
  filt_cutoff_data <- hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["corr_calc_out"]][["correlation_df_filt"]] %>%
    dplyr::filter(., rval >= hcobject[["cutoff_vec"]][x])

  # build network based on filtered data:
  filt_cutoff_graph <- igraph::graph_from_data_frame(filt_cutoff_data, directed = FALSE)

  # remove too small components:
  graph_components <- igraph::components(filt_cutoff_graph)

  # remove nodes from too small components from the network:
  gene_to_comp <- base::data.frame(gene = base::names(graph_components[["membership"]]), component = graph_components[["membership"]])
  comps_to_keep <- base::which(graph_components[["csize"]] >= hcobject[["global_settings"]][["min_nodes_number_for_network"]])
  nodes_to_remove <- dplyr::filter(gene_to_comp, !component %in% comps_to_keep) %>% dplyr::pull(., "gene")
  filt_cutoff_graph <- igraph::delete.vertices(filt_cutoff_graph, nodes_to_remove)

  # remove edges from removed nodes from the cutoff data:
  filt_cutoff_data <- dplyr::filter(filt_cutoff_data, V1 %in% igraph::V(filt_cutoff_graph)$name & V2 %in% igraph::V(filt_cutoff_graph)$name)

  filt_cutoff_counts <- hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["ds"]][base::row.names(hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["ds"]]) %in%
    base::names(igraph::V(filt_cutoff_graph)), ]
  corresp_info <- info_dataset[base::rownames(base::t(hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["topvar"]])) %in% base::rownames(info_dataset), ]

  output[["filt_cutoff_graph"]] <- filt_cutoff_graph
  output[["filt_cutoff_data"]] <- filt_cutoff_data


  message(
    "...After using the optimal cutoff of ", hcobject[["cutoff_vec"]][x], " the number of edges = ",
    base::nrow(filt_cutoff_data), " and the number of nodes = ", base::nrow(filt_cutoff_counts), "..."
  )

  col_list <- list()
  if (base::is.null(cols)) {
    for (i in all_conditions) {
      tmp_col <- ggsci::pal_nejm(alpha = 1)(8)
      # if there are more groups than colours, expand palette:
      if (base::length(base::unique(hcobject[["data"]][[paste0("set", x, "_anno")]][, i])) > base::length(tmp_col)) {
        tmp_col <- grDevices::colorRampPalette(tmp_col)(base::length(base::unique(hcobject[["data"]][[base::paste0("set", x, "_anno")]][, i])))
      } else {
        tmp_col <- tmp_col[base::seq_len(base::length(base::unique(hcobject[["data"]][[base::paste0("set", x, "_anno")]][, i])))]
      }
      base::names(tmp_col) <- base::unique(hcobject[["data"]][[base::paste0("set", x, "_anno")]][, i])
      col_list[[i]] <- tmp_col
    }
  } else {
    col_list <- cols
  }

  # Large matrices are rasterized explicitly to avoid the default
  # magick-based temp-file roundtrip, which can fail on some Windows setups.
  hm_use_raster <- base::nrow(filt_cutoff_counts) > 2000
  hm_raster_device <- if (isTRUE(base::capabilities("cairo"))) "CairoPNG" else "png"

  # filter annotation for pheatmap:
  anno_df <- dplyr::select(hcobject[["data"]][[base::paste0("set", x, "_anno")]], tidyselect::all_of(all_conditions))

  heatmap_filtered_counts <- ComplexHeatmap::pheatmap(
    mat = base::as.matrix(filt_cutoff_counts),
    color = .hc_default_gfc_colors(),
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_colors = col_list,
    annotation_col = anno_df,
    fontsize = 8,
    show_rownames = FALSE,
    show_colnames = TRUE,
    main = title,
    annotation_names_col = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_method = method,
    use_raster = hm_use_raster,
    raster_by_magick = FALSE,
    raster_device = hm_raster_device,
    heatmap_legend_param = list(title = "scaled expr.")
  )

  if (plot_HM) {
    ComplexHeatmap::plot.Heatmap(heatmap_filtered_counts)
  }

  output[["heatmap"]] <- heatmap_filtered_counts

  heatmap_export_file <- base::paste0(
    hcobject[["working_directory"]][["dir_output"]],
    hcobject[["global_settings"]][["save_folder"]],
    "/",
    "Heatmap_topvar_genes_",
    title,
    ".pdf"
  )
  output[["files"]] <- .hc_export_single_page_plot(
    file = heatmap_export_file,
    width = 7,
    height = 10,
    draw_fun = function() {
      ComplexHeatmap::plot.Heatmap(heatmap_filtered_counts)
    }
  )

  return(output)
}


#' Truncated Fold Changes
#'
#' Function calculates fold changes between given columns and a reference and truncates them if the FCs exceed the set maximum GFC range.
#' @param grp A vector of column names for which the fold changes shall be calculated.
#' @param group_means Name of the column containing the reference level from which the FC shall be calculated.
#' @param trans_norm The dataframe containing the data.
#' @noRd

gfc_calc <- function(grp, trans_norm, group_means) {
  df1 <- trans_norm[, grp]
  df2 <- gtools::foldchange(df1, group_means) %>%
    base::ifelse(. > hcobject[["global_settings"]][["range_GFC"]], hcobject[["global_settings"]][["range_GFC"]], .) %>%
    base::ifelse(. < (-hcobject[["global_settings"]][["range_GFC"]]), -hcobject[["global_settings"]][["range_GFC"]], .) %>%
    base::as.data.frame()
  base::colnames(df2) <- base::paste0("", grp)
  return(df2)
}


#' Calculate Group Fold Changes
#'
#' Function that computes fold changes for the define sample groups either with reference to a control group or to the mean of all groups.
#' @param info_dataset The annotation dataframe for the dataset that contains the group labels.
#' @param grouping_v A string giving a column name present in all annotation files, if this variable shall be used for grouping the samples instead of the variable of interest.
#' @param x An integer giving the number of the currently processed dataset.
#' @noRd


GFC_calculation <- function(info_dataset, grouping_v, x) {
  message("...calculate Group-Fold-Changes...")

  if (!base::is.null(grouping_v)) {
    info_dataset[["grpvar"]] <- info_dataset[, base::c(grouping_v)]
    message("User-defined variable ", grouping_v, " will be used for grouping the data.")
  } else if (base::intersect(hcobject[["global_settings"]][["voi"]], base::colnames(info_dataset)) %>% base::length() > 0) {
    message(
      "...Variable: '", hcobject[["global_settings"]][["voi"]],
      "' will be used as grouping variables..."
    )

    info_dataset[["grpvar"]] <- purrr::pmap(info_dataset[base::intersect(hcobject[["global_settings"]][["voi"]], base::colnames(info_dataset))],
      paste,
      sep = "-"
    ) %>% base::unlist()
  } else {
    message(
      "...The first column in the metadata will be used as the grouping variable",
      "since the voi_id is not present in the metadata..."
    )

    info_dataset[["grpvar"]] <- info_dataset[, 1]
  }

  # Normalize grouping labels and remove missing/NA-like entries to avoid
  # creating artificial "[NA]" groups in downstream heatmaps.
  info_dataset[["grpvar"]] <- base::trimws(base::as.character(info_dataset[["grpvar"]]))
  info_dataset[["grpvar"]][info_dataset[["grpvar"]] %in% c("", "NA", "<NA>", "[NA]", "[<NA>]")] <- NA_character_

  missing_grp <- base::is.na(info_dataset[["grpvar"]]) | !base::nzchar(info_dataset[["grpvar"]])
  if (base::any(missing_grp)) {
    warning(
      "Dropping ", base::sum(missing_grp),
      " sample(s) with missing grouping labels before GFC calculation.",
      call. = FALSE
    )
    info_dataset <- info_dataset[!missing_grp, , drop = FALSE]
  }
  if (base::nrow(info_dataset) == 0) {
    stop("No samples with valid grouping labels available for GFC calculation.")
  }


  if (hcobject[["global_settings"]][["control"]] == "none") {
    message("...GFC calculation with foldchange from mean...")

    # GFC calculation with fold-change from mean
    norm_data_anno <- base::merge(base::t(hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["topvar"]]), base::subset(info_dataset, select = base::c("grpvar")), by = "row.names", all.x = TRUE)

    norm_data_anno <- norm_data_anno[, -1]

    norm_data_anno <- norm_data_anno[, base::c(base::ncol(norm_data_anno), base::seq_len(base::ncol(norm_data_anno) - 1L))]

    trans_norm <- stats::setNames(base::data.frame(base::t(norm_data_anno[, -1])), norm_data_anno[, 1])

    if (hcobject[["global_settings"]][["data_in_log"]] == TRUE) {
      trans_norm <- .hc_antilog_impl(trans_norm, 2)
    }


    trans_norm <- base::t(base::apply(trans_norm, 1, function(i) base::tapply(i, base::colnames(trans_norm), base::mean)))
    trans_norm <- base::cbind(trans_norm, base::rowMeans(trans_norm))

    base::colnames(trans_norm)[base::ncol(trans_norm)] <- "group_mean"
    grplist <- base::colnames(trans_norm)[-(base::ncol(trans_norm))]


    GFC_all_genes <- base::do.call("cbind", base::lapply(grplist, gfc_calc, trans_norm = trans_norm, group_means = trans_norm[, "group_mean"]))
    GFC_all_genes <- base::round(GFC_all_genes, 3)
    GFC_all_genes$Gene <- base::rownames(GFC_all_genes)


    return(GFC_all_genes)
  } else {
    message("...GFC calculation with foldchange from control...")
    # GFC calculation with foldchange from control
    norm_data_anno <- base::merge(base::t(hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["topvar"]]), base::subset(info_dataset, select = base::c("grpvar")), by = "row.names", all.x = TRUE)


    # contains a column called Row.names (1st col)
    norm_data_anno <- norm_data_anno[, -1]

    norm_data_anno <- norm_data_anno[, base::c(base::ncol(norm_data_anno), base::seq_len(base::ncol(norm_data_anno) - 1L))]


    trans_norm <- stats::setNames(base::data.frame(base::t(norm_data_anno[, -1])), norm_data_anno[, 1])

    if (hcobject[["global_settings"]][["data_in_log"]] == TRUE) {
      trans_norm[trans_norm == 0] <- 1
      trans_norm <- .hc_antilog_impl(trans_norm, 2)
    }

    trans_norm <- base::t(base::apply(trans_norm, 1, function(i) base::tapply(i, base::colnames(trans_norm), base::mean)))


    # Identify the control group. `control` is documented as a substring of the
    # control group's label, so keep substring matching -- but use fixed = TRUE
    # (group labels are data, not regular expressions) and fail loudly when the
    # keyword does not resolve to exactly one group. Previously a typo matched
    # nothing, `cbind()` appended no column, and the rename below silently
    # turned the last (alphabetically sorted) condition into the reference, so
    # every GFC was computed against an arbitrary group with no warning.
    control_keyword <- base::as.character(hcobject[["global_settings"]][["control"]])
    is_ctrl <- base::grepl(control_keyword, base::colnames(trans_norm),
      ignore.case = TRUE, fixed = FALSE
    )
    if (base::sum(is_ctrl) == 0) {
      stop(
        "The control keyword '", control_keyword, "' does not match any sample group. ",
        "Available groups: ", base::paste(base::colnames(trans_norm), collapse = ", "),
        ". Set `control_keyword` to a string contained in the control group's label, ",
        "or use \"none\" to compute group fold changes against the mean of all groups.",
        call. = FALSE
      )
    }
    if (base::sum(is_ctrl) > 1) {
      stop(
        "The control keyword '", control_keyword, "' matches ", base::sum(is_ctrl),
        " sample groups (", base::paste(base::colnames(trans_norm)[is_ctrl], collapse = ", "),
        "). It must identify exactly one control group; choose a more specific keyword.",
        call. = FALSE
      )
    }

    ctrl_name <- base::colnames(trans_norm)[is_ctrl]
    trans_norm_no_ctrl <- trans_norm[, !is_ctrl, drop = FALSE]
    if (base::ncol(trans_norm_no_ctrl) == 0) {
      stop(
        "The only sample group ('", ctrl_name, "') was identified as the control group, ",
        "so no fold changes can be computed.",
        call. = FALSE
      )
    }

    trans_norm <- base::cbind(trans_norm_no_ctrl, trans_norm[, is_ctrl, drop = FALSE])

    base::rownames(trans_norm) <- base::rownames(trans_norm_no_ctrl)
    base::colnames(trans_norm)[base::ncol(trans_norm)] <- "group_mean"
    grplist <- base::colnames(trans_norm)[-(base::ncol(trans_norm))]

    GFC_all_genes <- base::do.call("cbind", base::lapply(grplist, gfc_calc, trans_norm = trans_norm, group_means = trans_norm[, "group_mean"]))
    GFC_all_genes <- base::round(GFC_all_genes, 3)
    base::rownames(GFC_all_genes) <- base::rownames(trans_norm)
    GFC_all_genes$Gene <- base::rownames(GFC_all_genes)
    # The control column was already excluded from `grplist`; drop it by exact
    # name rather than by substring, which used to also delete any non-control
    # condition (and potentially the "Gene" column) containing the keyword.
    GFC_all_genes <- GFC_all_genes[, base::colnames(GFC_all_genes) != ctrl_name, drop = FALSE]

    return(GFC_all_genes)
  }
}


#' Get Union Of Graphs
#'
#' The function builds a multigraph from the union of all layer-specific networks.
#'  Thus, also network parts that are unique to some layers will be present in the resulting integrated network, creating an integrated network that provides a wholistic, cross-dataset view on gene co-expression.
#' @noRd

get_union <- function() {
  combined_edgelist <- NULL

  for (x in base::seq_along(hcobject[["layer_specific_outputs"]])) {
    combined_edgelist <- base::rbind(combined_edgelist, hcobject[["layer_specific_outputs"]][[x]][["part2"]][["heatmap_out"]][["filt_cutoff_data"]])
  }
  combined_edgelist$weight <- combined_edgelist$rval
  combined_edgelist$rval <- NULL
  combined_edgelist$pval <- NULL


  .hc_set_bridge_hcobject_slot(c("integrated_output", "combined_edgelist"), combined_edgelist)
}


#' Get Intersection With Reference Graph
#'
#' It generates a multigraph of the reference network.
#'  The vertices of the resulting network will be identical to the reference network, but the edges connecting them will be greatly impacted by the other datasets.
#' @param with Either an integer giving the number of the dataset to be used as reference (e.g., 1) or the name given to the layer.
#' @noRd

get_intersection <- function(with) {
  # change to numeric representation of the dataset in case it was given as the name of a dataset:
  if (with %in% hcobject[["layers_names"]]) {
    with <- match(with, hcobject[["layers_names"]])
  }

  # change to numeric representation of the dataset in case it was given as a string:
  if (startsWith(with, "set")) {
    with <- base::strsplit(base::as.character(with), split = "set")[[1]][2] %>% base::as.numeric()
  }

  # get edgelist of the reference network:
  combined_edgelist <- hcobject[["layer_specific_outputs"]][[with]][["part2"]][["heatmap_out"]][["filt_cutoff_data"]]

  # Edge keys need a separator that cannot occur in a gene symbol: pasting the
  # names directly made ("MT","CO1") and ("M","TCO1") collide, which produced
  # spurious intersection edges.
  .edge_key <- function(a, b) {
    base::paste0(base::as.character(a), "\r", base::as.character(b))
  }

  combined_edgelist$merged <- .edge_key(combined_edgelist$V1, combined_edgelist$V2)
  combined_edgelist$revmerged <- .edge_key(combined_edgelist$V2, combined_edgelist$V1)
  # iterate over datasets:
  for (x in base::seq_along(hcobject[["layer_specific_outputs"]])) {
    if (!x == with) {
      # get edgelist of current dataset:
      tmp <- hcobject[["layer_specific_outputs"]][[x]][["part2"]][["heatmap_out"]][["filt_cutoff_data"]]

      tmp$merged <- .edge_key(tmp$V1, tmp$V2)

      # check which edges overlap with reference network:
      tmp <- dplyr::filter(tmp, merged %in% combined_edgelist$merged | merged %in% combined_edgelist$revmerged)
      combined_edgelist <- base::rbind(combined_edgelist, tmp)
    }
  }
  combined_edgelist$weight <- combined_edgelist$rval
  combined_edgelist$rval <- NULL
  combined_edgelist$pval <- NULL
  combined_edgelist$merged <- NULL
  combined_edgelist$revmerged <- NULL

  .hc_set_bridge_hcobject_slot(c("integrated_output", "combined_edgelist"), combined_edgelist)
}

#' Merge GFCs From Different Datasets
#'
#' Combines the GFC values per gene from the different datasets and substitutes missing values.
#' @param GFC_when_missing The GFC value to enter when a gene has not been measured in a dataset, but in others. Default is the lower bound of the set GFC range.
#' @noRd

merge_GFCs <- function(GFC_when_missing = -hcobject[["global_settings"]][["range_GFC"]]) {
  col_names_new_GFC <- NULL
  for (x in base::seq_along(hcobject[["layer_specific_outputs"]])) {
    new_col <- base::colnames(hcobject[["layer_specific_outputs"]][[x]][["part2"]][["GFC_all_genes"]])[-base::ncol(hcobject[["layer_specific_outputs"]][[x]][["part2"]][["GFC_all_genes"]])]

    col_names_new_GFC <- base::c(col_names_new_GFC, new_col)
  }
  new_GFC <- NULL
  for (y in igraph::get.vertex.attribute(hcobject[["integrated_output"]][["merged_net"]])$name) {
    line <- NULL


    for (z in base::seq_along(hcobject[["layer_specific_outputs"]])) {
      if (y %in% hcobject[["layer_specific_outputs"]][[z]][["part2"]][["GFC_all_genes"]][["Gene"]]) {
        GFC_tmp <- hcobject[["layer_specific_outputs"]][[z]][["part2"]][["GFC_all_genes"]]

        line <- base::c(line, GFC_tmp[GFC_tmp$Gene == y, base::colnames(GFC_tmp)[base::seq_len(base::ncol(GFC_tmp) - 1L)]]) %>%
          base::unlist(.)
      } else {
        line <- base::c(line, base::rep(GFC_when_missing, (base::length(base::colnames(hcobject[["layer_specific_outputs"]][[z]][["part2"]][["GFC_all_genes"]])) - 1)))
      }
    }

    line <- base::as.data.frame(line) %>% base::t()
    new_GFC <- base::rbind(new_GFC, line)
  }
  new_GFC <- base::as.data.frame(new_GFC)

  base::colnames(new_GFC) <- base::c(col_names_new_GFC)
  new_GFC$Gene <- igraph::get.vertex.attribute(hcobject[["integrated_output"]][["merged_net"]])$name
  rownames(new_GFC) <- new_GFC$Gene
  return(new_GFC)
}


#' Fixes Missing Slashes
#'
#' Adds slash to the end of provided directory if it is missing and gives and error if the provided directory does not exist.
#' @param directory A string describing a directory path.
#' @noRd

fix_dir <- function(directory) {
  if (base::isFALSE(directory) || base::identical(directory, FALSE)) {
    return(FALSE)
  }
  if (base::length(directory) != 1) {
    stop("`directory` must be a scalar path string or FALSE.")
  }
  directory <- base::as.character(directory[[1]])
  if (base::is.na(directory) || !base::nzchar(directory)) {
    stop("`directory` must be a non-empty path string or FALSE.")
  }
  # check existance of path:
  if (!base::dir.exists(directory)) {
    stop("The directory '", directory, "' does not exist.")
  }
  # add slash if missing at the end:
  split_dir <- base::strsplit(directory, split = "")[[1]]
  if (split_dir[base::length(split_dir)] == "/") {
    return(directory)
  } else {
    return(base::paste0(directory, "/"))
  }
}


#' Non-Interactive Plot of Cutoff Statistics
#'
#' Plots the R-squared vlaue, number of edges, number of genes and number of networks for different cut-offs as a ggplot.
#' @param cutoff_stats A dataframe of cutoff statistics generated in previous steps.
#' @param hline A list with four slots ("R.squared", "no_edges", "no_nodes", "no_networks") each of which can be set either to NULL (default) or a number to introduce a horizontal line for orientation at that value in the respective plot.
#' @param x An integer giving the number of the currently processed layer.
#' @param tuning_marker Optional list with selected/strict/relaxed/simple cutoff marker metadata.
#' @noRd

plot_cutoffs_internal_static <- function(cutoff_stats,
                                         hline = list("R.squared" = NULL, "no_edges" = NULL, "no_nodes" = NULL, "no_networks" = NULL),
                                         x,
                                         tuning_marker = NULL) {
  as_num1 <- function(z) {
    out <- .hc_as_numeric_safely(z)
    if (base::length(out) < 1) {
      return(NA_real_)
    }
    out[[1]]
  }
  as_chr1 <- function(z) {
    out <- base::as.character(z)
    if (base::length(out) < 1) {
      return(NA_character_)
    }
    out[[1]]
  }
  normalize_policy <- function(z) {
    pol <- base::tolower(as_chr1(z))
    if (!base::is.character(pol) || base::is.na(pol) || !base::nzchar(pol)) {
      return(NA_character_)
    }
    base::switch(pol,
      tier1 = "strict",
      tier2 = "relaxed",
      fallback = "best_available",
      pol
    )
  }
  same_cutoff <- function(a, b, tol = 1e-8) {
    base::is.finite(a) && base::is.finite(b) && base::abs(a - b) <= tol
  }
  format_marker_label <- function(value, metric) {
    if (!base::is.finite(value)) {
      return(NA_character_)
    }
    if (metric == "R.squared") {
      return(base::formatC(value, format = "f", digits = 3))
    }
    if (metric == "no_of_networks") {
      return(base::formatC(base::round(value), format = "f", digits = 0))
    }
    scales::comma(base::round(value))
  }
  extract_marker_stats <- function(cutoff_value) {
    if (!base::is.finite(cutoff_value)) {
      return(NULL)
    }
    tmp <- cutoff_stats[, c("corr", "R.squared", "no_edges", "no_nodes", "no_of_networks"), drop = FALSE]
    tmp[["corr"]] <- .hc_as_numeric_safely(tmp[["corr"]])
    tmp <- tmp[base::order(tmp[["corr"]]), , drop = FALSE]
    metric_at_cutoff <- function(metric) {
      y <- .hc_as_numeric_safely(tmp[[metric]])
      ok <- base::is.finite(tmp[["corr"]]) & base::is.finite(y)
      if (base::sum(ok) < 1) {
        return(NA_real_)
      }
      x_ok <- tmp[["corr"]][ok]
      y_ok <- y[ok]
      if (base::length(base::unique(x_ok)) >= 2) {
        return(as.numeric(stats::approx(x = x_ok, y = y_ok, xout = cutoff_value, ties = base::mean, rule = 2)$y))
      }
      y_ok[[1]]
    }
    base::data.frame(
      corr = cutoff_value,
      R.squared = metric_at_cutoff("R.squared"),
      no_edges = metric_at_cutoff("no_edges"),
      no_nodes = metric_at_cutoff("no_nodes"),
      no_of_networks = metric_at_cutoff("no_of_networks")
    )
  }

  cutoff_stats[["corr"]] <- base::rownames(cutoff_stats) %>% base::as.numeric()
  layer_title <- base::as.character(hcobject[["layers_names"]][x])

  selected_cutoff <- as_num1(tuning_marker[["selected_cutoff"]])
  if (!base::is.finite(selected_cutoff)) {
    selected_cutoff <- as_num1(tuning_marker[["tiered_cutoff"]])
  }
  selected_policy <- normalize_policy(tuning_marker[["selected_policy"]])
  if (!base::is.character(selected_policy) || base::is.na(selected_policy) || !base::nzchar(selected_policy)) {
    selected_policy <- normalize_policy(tuning_marker[["tier"]])
  }
  strict_cutoff <- as_num1(tuning_marker[["strict_cutoff"]])
  relaxed_cutoff <- as_num1(tuning_marker[["relaxed_cutoff"]])
  best_available_cutoff <- as_num1(tuning_marker[["best_available_cutoff"]])
  simple_cutoff <- as_num1(tuning_marker[["simple_cutoff"]])

  if (!base::is.finite(strict_cutoff) && base::identical(selected_policy, "strict") && base::is.finite(selected_cutoff)) {
    strict_cutoff <- selected_cutoff
  }
  if (!base::is.finite(relaxed_cutoff) && base::identical(selected_policy, "relaxed") && base::is.finite(selected_cutoff)) {
    relaxed_cutoff <- selected_cutoff
  }
  if (!base::is.finite(best_available_cutoff) && base::identical(selected_policy, "best_available") && base::is.finite(selected_cutoff)) {
    best_available_cutoff <- selected_cutoff
  }

  selected_dash <- base::switch(selected_policy,
    strict = "solid",
    relaxed = "dashed",
    best_available = "dotdash",
    "dotted"
  )
  selected_color <- base::switch(selected_policy,
    strict = "#1B9E77",
    relaxed = "#E69F00",
    best_available = "#D73027",
    "#7B2CBF"
  )
  x_vals <- cutoff_stats[["corr"]]
  nudge_base <- 0
  if (base::length(x_vals) > 1) {
    x_range <- base::max(x_vals, na.rm = TRUE) - base::min(x_vals, na.rm = TRUE)
    if (base::is.finite(x_range) && x_range > 0) {
      nudge_base <- 0.02 * x_range
    }
  }

  add_cutoff_markers <- function(p, show_legend = FALSE) {
    if (base::is.finite(relaxed_cutoff) && !isTRUE(show_legend)) {
      p <- p + ggplot2::geom_vline(
        xintercept = relaxed_cutoff,
        color = "#E69F00",
        linetype = "dashed",
        linewidth = if (base::identical(selected_policy, "relaxed") && same_cutoff(selected_cutoff, relaxed_cutoff)) 1.4 else 0.9
      )
    }
    if (base::is.finite(strict_cutoff) && !isTRUE(show_legend)) {
      p <- p + ggplot2::geom_vline(
        xintercept = strict_cutoff,
        color = "#1B9E77",
        linetype = "solid",
        linewidth = if (base::identical(selected_policy, "strict") && same_cutoff(selected_cutoff, strict_cutoff)) 1.4 else 0.9
      )
    }
    if (base::is.finite(best_available_cutoff)) {
      p <- p + ggplot2::geom_vline(
        xintercept = best_available_cutoff,
        color = "#D73027",
        linetype = "dotdash",
        linewidth = if (base::identical(selected_policy, "best_available") && same_cutoff(selected_cutoff, best_available_cutoff)) 1.4 else 0.9
      )
    }
    if (base::is.finite(simple_cutoff)) {
      p <- p + ggplot2::geom_vline(
        xintercept = simple_cutoff,
        color = "#2D2D2D",
        linetype = "longdash",
        linewidth = 0.8
      )
    }

    has_selected_reference <- same_cutoff(selected_cutoff, strict_cutoff) ||
      same_cutoff(selected_cutoff, relaxed_cutoff) ||
      same_cutoff(selected_cutoff, best_available_cutoff) ||
      same_cutoff(selected_cutoff, simple_cutoff)
    if (base::is.finite(selected_cutoff) && !has_selected_reference) {
      p <- p + ggplot2::geom_vline(
        xintercept = selected_cutoff,
        color = selected_color,
        linetype = selected_dash,
        linewidth = 1.3
      )
    }
    if (isTRUE(show_legend)) {
      legend_df <- base::data.frame(ref = base::character(0), cutoff = base::numeric(0), stringsAsFactors = FALSE)
      if (base::is.finite(strict_cutoff)) {
        legend_df <- base::rbind(
          legend_df,
          base::data.frame(ref = "strict", cutoff = strict_cutoff, stringsAsFactors = FALSE)
        )
      }
      if (base::is.finite(relaxed_cutoff)) {
        legend_df <- base::rbind(
          legend_df,
          base::data.frame(ref = "relaxed", cutoff = relaxed_cutoff, stringsAsFactors = FALSE)
        )
      }
      if (base::nrow(legend_df) > 0) {
        legend_breaks <- base::intersect(c("strict", "relaxed"), legend_df$ref)
        p <- p +
          ggplot2::geom_vline(
            data = legend_df,
            mapping = ggplot2::aes(xintercept = cutoff, color = ref, linetype = ref),
            inherit.aes = FALSE,
            linewidth = 1.1,
            show.legend = TRUE
          ) +
          ggplot2::scale_color_manual(
            name = "Cutoff lines",
            values = c(strict = "#1B9E77", relaxed = "#E69F00"),
            breaks = legend_breaks
          ) +
          ggplot2::scale_linetype_manual(
            name = "Cutoff lines",
            values = c(strict = "solid", relaxed = "dashed"),
            breaks = legend_breaks
          ) +
          ggplot2::guides(
            color = ggplot2::guide_legend(order = 1),
            linetype = ggplot2::guide_legend(order = 1)
          ) +
          ggplot2::theme(legend.position = "right")
      }
    } else {
      p <- p + ggplot2::guides(color = "none", linetype = "none")
    }
    p
  }

  add_policy_marker_values <- function(p, metric_col) {
    policies <- list(
      list(id = "strict", cutoff = strict_cutoff, color = "#1B9E77", prefix = "S ", nudge_x = nudge_base, vjust = -0.55),
      list(id = "relaxed", cutoff = relaxed_cutoff, color = "#E69F00", prefix = "R ", nudge_x = -nudge_base, vjust = 1.35)
    )
    for (pol in policies) {
      if (!base::is.finite(pol$cutoff)) {
        next
      }
      st <- extract_marker_stats(pol$cutoff)
      if (base::is.null(st)) {
        next
      }
      x_val <- .hc_first_numeric_value(st[["corr"]][[1]])
      y_val <- .hc_first_numeric_value(st[[metric_col]][[1]])
      if (!base::is.finite(x_val) || !base::is.finite(y_val)) {
        next
      }
      label_txt <- format_marker_label(y_val, metric = metric_col)
      if (!base::is.character(label_txt) || base::is.na(label_txt) || !base::nzchar(label_txt)) {
        next
      }
      label_txt <- base::paste0(pol$prefix, label_txt)
      pt_df <- base::data.frame(corr = x_val, yval = y_val, label = label_txt)
      p <- p +
        ggplot2::geom_point(
          data = pt_df,
          mapping = ggplot2::aes(x = corr, y = yval),
          inherit.aes = FALSE,
          color = pol$color,
          size = 2.5
        ) +
        ggplot2::geom_label(
          data = pt_df,
          mapping = ggplot2::aes(x = corr, y = yval, label = label),
          inherit.aes = FALSE,
          nudge_x = pol$nudge_x,
          vjust = pol$vjust,
          size = 2.8,
          color = pol$color,
          fill = "white",
          linewidth = 0.15,
          label.padding = grid::unit(0.08, "lines")
        )
    }
    p
  }

  # plot R-squared values:
  p1 <- ggplot2::ggplot(cutoff_stats, ggplot2::aes(x = corr)) +
    ggplot2::geom_vline(xintercept = base::c(base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.01)), color = "darkgrey") +
    ggplot2::geom_hline(yintercept = hline[[1]], color = "darkgrey", linewidth = 1.7) +
    ggplot2::geom_line(ggplot2::aes(y = R.squared), color = "#374E55FF") +
    ggplot2::geom_point(ggplot2::aes(y = R.squared), color = "#374E55FF", size = 2) +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.005)) +
    ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
    ggplot2::ggtitle(layer_title)
  p1 <- add_cutoff_markers(p1, show_legend = TRUE)
  p1 <- add_policy_marker_values(p1, metric_col = "R.squared")

  # plot number of edges:
  p2 <- ggplot2::ggplot(cutoff_stats, ggplot2::aes(x = corr)) +
    ggplot2::geom_vline(xintercept = base::c(base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.01)), color = "darkgrey") +
    ggplot2::geom_hline(yintercept = hline[[2]], color = "darkgrey", linewidth = 1.7) +
    ggplot2::geom_line(ggplot2::aes(y = no_edges), color = "#DF8F44FF") +
    ggplot2::geom_point(ggplot2::aes(y = no_edges), color = "#DF8F44FF", size = 2) +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.005)) +
    ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0.05, 0.18)))
  p2 <- add_cutoff_markers(p2, show_legend = FALSE)
  p2 <- add_policy_marker_values(p2, metric_col = "no_edges")

  # plot number of nodes:
  p3 <- ggplot2::ggplot(cutoff_stats, ggplot2::aes(x = corr)) +
    ggplot2::geom_vline(xintercept = base::c(base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.01)), color = "darkgrey") +
    ggplot2::geom_hline(yintercept = hline[[3]], color = "darkgrey", linewidth = 1.7) +
    ggplot2::geom_line(ggplot2::aes(y = no_nodes), color = "#00A1D5FF") +
    ggplot2::geom_point(ggplot2::aes(y = no_nodes), color = "#00A1D5FF", size = 2) +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.005)) +
    ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0.05, 0.18)))
  p3 <- add_cutoff_markers(p3, show_legend = FALSE)
  p3 <- add_policy_marker_values(p3, metric_col = "no_nodes")

  # plot number of networks:
  p4 <- ggplot2::ggplot(cutoff_stats, ggplot2::aes(x = corr)) +
    ggplot2::geom_vline(xintercept = base::c(base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.01)), color = "darkgrey") +
    ggplot2::geom_hline(yintercept = hline[[4]], color = "darkgrey", linewidth = 1.7) +
    ggplot2::geom_point(ggplot2::aes(y = no_of_networks), color = "#B24745FF", size = 2) +
    ggplot2::scale_x_continuous(breaks = base::seq(from = base::min(cutoff_stats[["corr"]]), to = 1, by = 0.005)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18)))
  p4 <- add_cutoff_markers(p4, show_legend = FALSE)
  p4 <- add_policy_marker_values(p4, metric_col = "no_of_networks")


  p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1, align = "v")
  return(p)
}

#' Interactive Plot of Cutoff Statistics
#'
#' Plots the R-squared vlaue, number of edges, number of genes and number of networks for different cut-offs as an interactive widget using plotly.
#' @param cutoff_stats A dataframe of cutoff statistics generated in previous steps.
#' @param x An integer giving the number of the currently processed layer.
#' @param tuning_marker Optional list with selected/strict/relaxed/simple cutoff marker metadata.
#' @noRd

plot_cutoffs_internal_interactive <- function(cutoff_stats,
                                              x,
                                              tuning_marker = NULL) {
  as_num1 <- function(z) {
    out <- .hc_as_numeric_safely(z)
    if (base::length(out) < 1) {
      return(NA_real_)
    }
    out[[1]]
  }
  as_chr1 <- function(z) {
    out <- base::as.character(z)
    if (base::length(out) < 1) {
      return(NA_character_)
    }
    out[[1]]
  }
  normalize_policy <- function(z) {
    pol <- base::tolower(as_chr1(z))
    if (!base::is.character(pol) || base::is.na(pol) || !base::nzchar(pol)) {
      return(NA_character_)
    }
    base::switch(pol,
      tier1 = "strict",
      tier2 = "relaxed",
      fallback = "best_available",
      pol
    )
  }
  same_cutoff <- function(a, b, tol = 1e-8) {
    base::is.finite(a) && base::is.finite(b) && base::abs(a - b) <= tol
  }
  format_marker_label <- function(value, metric) {
    if (!base::is.finite(value)) {
      return(NA_character_)
    }
    if (metric == "R.squared") {
      return(base::formatC(value, format = "f", digits = 3))
    }
    if (metric == "no_of_networks") {
      return(base::formatC(base::round(value), format = "f", digits = 0))
    }
    scales::comma(base::round(value))
  }
  extract_marker_stats <- function(cutoff_value) {
    if (!base::is.finite(cutoff_value)) {
      return(NULL)
    }
    tmp <- cutoff_stats[, c("corr", "R.squared", "no_edges", "no_nodes", "no_of_networks"), drop = FALSE]
    tmp[["corr"]] <- .hc_as_numeric_safely(tmp[["corr"]])
    tmp <- tmp[base::order(tmp[["corr"]]), , drop = FALSE]
    metric_at_cutoff <- function(metric) {
      y <- .hc_as_numeric_safely(tmp[[metric]])
      ok <- base::is.finite(tmp[["corr"]]) & base::is.finite(y)
      if (base::sum(ok) < 1) {
        return(NA_real_)
      }
      x_ok <- tmp[["corr"]][ok]
      y_ok <- y[ok]
      if (base::length(base::unique(x_ok)) >= 2) {
        return(as.numeric(stats::approx(x = x_ok, y = y_ok, xout = cutoff_value, ties = base::mean, rule = 2)$y))
      }
      y_ok[[1]]
    }
    base::data.frame(
      corr = cutoff_value,
      R.squared = metric_at_cutoff("R.squared"),
      no_edges = metric_at_cutoff("no_edges"),
      no_nodes = metric_at_cutoff("no_nodes"),
      no_of_networks = metric_at_cutoff("no_of_networks")
    )
  }

  cutoff_stats[["corr"]] <- base::rownames(cutoff_stats) %>% base::as.numeric()
  selected_cutoff <- as_num1(tuning_marker[["selected_cutoff"]])
  if (!base::is.finite(selected_cutoff)) {
    selected_cutoff <- as_num1(tuning_marker[["tiered_cutoff"]])
  }
  selected_policy <- normalize_policy(tuning_marker[["selected_policy"]])
  if (!base::is.character(selected_policy) || base::is.na(selected_policy) || !base::nzchar(selected_policy)) {
    selected_policy <- normalize_policy(tuning_marker[["tier"]])
  }
  strict_cutoff <- as_num1(tuning_marker[["strict_cutoff"]])
  relaxed_cutoff <- as_num1(tuning_marker[["relaxed_cutoff"]])
  best_available_cutoff <- as_num1(tuning_marker[["best_available_cutoff"]])
  simple_cutoff <- as_num1(tuning_marker[["simple_cutoff"]])
  if (!base::is.finite(strict_cutoff) && base::identical(selected_policy, "strict") && base::is.finite(selected_cutoff)) {
    strict_cutoff <- selected_cutoff
  }
  if (!base::is.finite(relaxed_cutoff) && base::identical(selected_policy, "relaxed") && base::is.finite(selected_cutoff)) {
    relaxed_cutoff <- selected_cutoff
  }
  if (!base::is.finite(best_available_cutoff) && base::identical(selected_policy, "best_available") && base::is.finite(selected_cutoff)) {
    best_available_cutoff <- selected_cutoff
  }
  selected_dash <- base::switch(selected_policy,
    strict = "solid",
    relaxed = "dash",
    best_available = "dashdot",
    "dot"
  )
  selected_color <- base::switch(selected_policy,
    strict = "#1B9E77",
    relaxed = "#E69F00",
    best_available = "#D73027",
    "#7B2CBF"
  )
  title_txt <- base::paste0("Cut-off selection guide: ", hcobject[["layers_names"]][x])

  # plot R-squared value:
  p1 <- plotly::plot_ly(cutoff_stats,
    x = ~corr, y = ~R.squared, type = "scatter",
    mode = "lines+markers", name = "R^2", line = list(color = "#374E55FF"), marker = list(color = "#374E55FF")
  )
  # plot number of edges:
  p2 <- plotly::plot_ly(cutoff_stats,
    x = ~corr, y = ~no_edges, type = "scatter",
    mode = "lines+markers", name = "no. edges", line = list(color = "#DF8F44FF"), marker = list(color = "#DF8F44FF")
  )
  # plot number of nodes:
  p3 <- plotly::plot_ly(cutoff_stats,
    x = ~corr, y = ~no_nodes, type = "scatter",
    mode = "lines+markers", name = "no. nodes", line = list(color = "#00A1D5FF"), marker = list(color = "#00A1D5FF")
  )
  # plot number of networks:
  p4 <- plotly::plot_ly(cutoff_stats,
    x = ~corr, y = ~no_of_networks, type = "scatter",
    mode = "markers", name = "no. networks", marker = list(color = "#B24745FF")
  )
  # combine plots:
  p <- plotly::subplot(p1, p2, p3, p4, nrows = 4, shareX = TRUE)
  add_reference_legend_trace <- function(p_obj, cutoff_value, nm, color, dash) {
    if (!base::is.finite(cutoff_value)) {
      return(p_obj)
    }
    x_vals <- .hc_as_numeric_safely(cutoff_stats[["corr"]])
    x_vals <- x_vals[base::is.finite(x_vals)]
    if (base::length(x_vals) >= 2) {
      x_ref <- base::range(x_vals)
    } else if (base::length(x_vals) == 1) {
      x_ref <- base::rep(x_vals[[1]], 2)
    } else {
      x_ref <- c(0, 1)
    }
    y_ref <- .hc_first_numeric_value(cutoff_stats[["R.squared"]][[1]])
    if (!base::is.finite(y_ref)) {
      y_ref <- 0
    }
    legend_df <- base::data.frame(corr = x_ref, y = base::rep(y_ref, 2))
    p_obj %>% plotly::add_trace(
      data = legend_df,
      x = ~corr,
      y = ~y,
      type = "scatter",
      mode = "lines",
      line = list(color = color, dash = dash, width = 1.4),
      name = nm,
      showlegend = TRUE,
      visible = "legendonly",
      hoverinfo = "skip",
      inherit = FALSE
    )
  }
  p <- add_reference_legend_trace(p, strict_cutoff, "strict cutoff", "#1B9E77", "solid")
  p <- add_reference_legend_trace(p, relaxed_cutoff, "relaxed cutoff", "#E69F00", "dash")

  # adjust layout:
  steps <- list()
  for (i in base::seq_along(cutoff_stats[["corr"]])) {
    step <- list(
      args = list("marker.color", list(
        base::rep("#374E55FF", base::length(cutoff_stats[["corr"]])),
        base::rep("#DF8F44FF", base::length(cutoff_stats[["corr"]])),
        base::rep("#00A1D5FF", base::length(cutoff_stats[["corr"]])),
        base::rep("#B24745FF", base::length(cutoff_stats[["corr"]]))
      )),
      label = base::paste0(
        base::as.character(cutoff_stats[["corr"]][i]), ", R^2: ",
        base::round(cutoff_stats[["R.squared"]][i], 3),
        "; no. edges: ", cutoff_stats[["no_edges"]][i], "; no. nodes: ",
        cutoff_stats[["no_nodes"]][i], "; no. networks: ", cutoff_stats[["no_of_networks"]][i]
      ),
      method = "restyle"
    )

    for (j in base::seq_len(4L)) {
      step[["args"]][[2]][[j]][i] <- "red"
    }

    steps[[i]] <- step
  }
  add_shape <- function(shape_list, x0, color, dash, width) {
    if (!base::is.finite(x0)) {
      return(shape_list)
    }
    shape_list[[base::length(shape_list) + 1]] <- list(
      type = "line",
      x0 = x0,
      x1 = x0,
      y0 = 0,
      y1 = 1,
      xref = "x",
      yref = "paper",
      line = list(color = color, dash = dash, width = width)
    )
    shape_list
  }
  shape_list <- list()
  shape_list <- add_shape(
    shape_list = shape_list,
    x0 = relaxed_cutoff,
    color = "#E69F00",
    dash = "dash",
    width = if (base::identical(selected_policy, "relaxed") && same_cutoff(selected_cutoff, relaxed_cutoff)) 1.5 else 1.0
  )
  shape_list <- add_shape(
    shape_list = shape_list,
    x0 = strict_cutoff,
    color = "#1B9E77",
    dash = "solid",
    width = if (base::identical(selected_policy, "strict") && same_cutoff(selected_cutoff, strict_cutoff)) 1.5 else 1.0
  )
  shape_list <- add_shape(
    shape_list = shape_list,
    x0 = best_available_cutoff,
    color = "#D73027",
    dash = "dashdot",
    width = if (base::identical(selected_policy, "best_available") && same_cutoff(selected_cutoff, best_available_cutoff)) 1.5 else 1.0
  )
  shape_list <- add_shape(
    shape_list = shape_list,
    x0 = simple_cutoff,
    color = "#2D2D2D",
    dash = "longdash",
    width = 1.2
  )
  has_selected_reference <- same_cutoff(selected_cutoff, strict_cutoff) ||
    same_cutoff(selected_cutoff, relaxed_cutoff) ||
    same_cutoff(selected_cutoff, best_available_cutoff) ||
    same_cutoff(selected_cutoff, simple_cutoff)
  if (base::is.finite(selected_cutoff) && !has_selected_reference) {
    shape_list <- add_shape(
      shape_list = shape_list,
      x0 = selected_cutoff,
      color = selected_color,
      dash = selected_dash,
      width = 1.4
    )
  }
  annotation_list <- list()
  policy_ann <- list(
    list(id = "strict", cutoff = strict_cutoff, color = "#1B9E77", prefix = "S ", ax = -30, ay = -18),
    list(id = "relaxed", cutoff = relaxed_cutoff, color = "#E69F00", prefix = "R ", ax = 30, ay = 18)
  )
  marker_specs <- list(
    list(metric = "R.squared", yref = "y"),
    list(metric = "no_edges", yref = "y2"),
    list(metric = "no_nodes", yref = "y3"),
    list(metric = "no_of_networks", yref = "y4")
  )
  for (pol in policy_ann) {
    if (!base::is.finite(pol$cutoff)) {
      next
    }
    pol_stats <- extract_marker_stats(pol$cutoff)
    if (base::is.null(pol_stats)) {
      next
    }
    x_val <- .hc_first_numeric_value(pol_stats[["corr"]][[1]])
    if (!base::is.finite(x_val)) {
      next
    }
    for (sp in marker_specs) {
      y_val <- .hc_first_numeric_value(pol_stats[[sp$metric]][[1]])
      if (!base::is.finite(y_val)) {
        next
      }
      label_txt <- format_marker_label(y_val, metric = sp$metric)
      if (!base::is.character(label_txt) || base::is.na(label_txt) || !base::nzchar(label_txt)) {
        next
      }
      annotation_list[[base::length(annotation_list) + 1]] <- list(
        x = x_val,
        y = y_val,
        xref = "x",
        yref = sp$yref,
        text = base::paste0(pol$prefix, label_txt),
        showarrow = TRUE,
        arrowhead = 1,
        arrowsize = 0.8,
        arrowwidth = 1,
        arrowcolor = pol$color,
        ax = pol$ax,
        ay = pol$ay,
        font = list(color = pol$color, size = 10),
        bgcolor = "rgba(255,255,255,0.85)",
        bordercolor = pol$color,
        borderwidth = 1,
        borderpad = 2
      )
    }
  }

  p <- p %>%
    plotly::layout(hovermode = "x unified") %>%
    plotly::layout(
      title = title_txt,
      sliders = list(
        list(
          pad = list(t = 60),
          active = 2,
          currentvalue = list(prefix = "Cut-off: ", font = list(color = "black", size = 14)),
          steps = steps,
          font = list(color = "white", size = 0)
        )
      ),
      shapes = shape_list,
      annotations = annotation_list
    )
  return(p)
}

.hc_heatmap_extract_matrix <- function(heatmap_obj) {
  if (inherits(heatmap_obj, "Heatmap")) {
    return(tryCatch(heatmap_obj@matrix, error = function(e) NULL))
  }
  if (inherits(heatmap_obj, "HeatmapList")) {
    return(tryCatch(heatmap_obj@ht_list[[1]]@matrix, error = function(e) NULL))
  }
  NULL
}

.hc_normalize_heatmap_axis_order <- function(order_value, axis_ids) {
  if (is.null(axis_ids) || length(axis_ids) == 0) {
    return(NULL)
  }
  if (is.list(order_value) && length(order_value) > 0) {
    order_value <- order_value[[1]]
  }

  ord <- NULL
  if (is.numeric(order_value) && length(order_value) == length(axis_ids)) {
    idx <- as.integer(order_value)
    idx <- idx[!is.na(idx) & idx >= 1L & idx <= length(axis_ids)]
    ord <- axis_ids[idx]
  } else if (is.character(order_value) && length(order_value) > 0) {
    ord <- as.character(order_value)
  }

  if (is.null(ord)) {
    return(as.character(axis_ids))
  }

  ord <- ord[!is.na(ord) & nzchar(ord) & ord %in% axis_ids]
  .hc_match_axis_order_with_duplicates(as.character(axis_ids), ord)
}

.hc_heatmap_cache_info <- function(cluster_calc) {
  stored_raw <- tryCatch(cluster_calc[["heatmap_cluster_raw"]], error = function(e) NULL)
  stored_hm <- tryCatch(cluster_calc[["heatmap_cluster"]], error = function(e) NULL)
  mat <- tryCatch(cluster_calc[["heatmap_matrix"]], error = function(e) NULL)

  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  } else if (!is.null(mat) && !is.matrix(mat)) {
    mat <- tryCatch(as.matrix(mat), error = function(e) NULL)
  }

  if (is.null(mat)) {
    source_obj <- if (!is.null(stored_raw)) stored_raw else stored_hm
    mat <- .hc_heatmap_extract_matrix(source_obj)
  }

  row_ids <- if (!is.null(mat)) rownames(mat) else NULL
  col_ids <- if (!is.null(mat)) colnames(mat) else NULL

  row_order <- tryCatch(cluster_calc[["heatmap_row_order"]], error = function(e) NULL)
  col_order <- tryCatch(cluster_calc[["heatmap_column_order"]], error = function(e) NULL)

  if (is.null(row_order) && !is.null(stored_hm)) {
    row_order <- tryCatch(ComplexHeatmap::row_order(stored_hm), error = function(e) NULL)
  }
  if (is.null(col_order) && !is.null(stored_hm)) {
    col_order <- tryCatch(ComplexHeatmap::column_order(stored_hm), error = function(e) NULL)
  }

  list(
    matrix = mat,
    row_order = .hc_normalize_heatmap_axis_order(row_order, row_ids),
    col_order = .hc_normalize_heatmap_axis_order(col_order, col_ids),
    heatmap_obj = stored_hm,
    raw_heatmap_obj = stored_raw
  )
}

.hc_select_heatmap_col_order <- function(available_cols,
                                         plot_order = NULL,
                                         main_order = NULL,
                                         fallback_order = NULL,
                                         context = "heatmap") {
  available_cols <- as.character(available_cols)
  if (length(available_cols) == 0) {
    return(available_cols)
  }

  resolve_candidate <- function(candidate, warn = FALSE) {
    if (is.null(candidate) || length(candidate) == 0) {
      return(NULL)
    }
    candidate <- as.character(candidate)
    candidate <- candidate[!is.na(candidate) & nzchar(candidate)]
    if (length(candidate) == 0 || !any(candidate %in% available_cols)) {
      return(NULL)
    }
    if (isTRUE(warn)) {
      return(.hc_resolve_heatmap_col_order(
        mat_cols = available_cols,
        requested_order = candidate,
        context = context
      ))
    }
    .hc_resolve_heatmap_col_order(
      mat_cols = available_cols,
      requested_order = candidate,
      context = context,
      warn_on_missing = FALSE
    )
  }

  resolved <- resolve_candidate(plot_order, warn = TRUE)
  if (!is.null(resolved)) {
    return(resolved)
  }

  resolved <- resolve_candidate(main_order, warn = FALSE)
  if (!is.null(resolved)) {
    return(resolved)
  }

  resolved <- resolve_candidate(fallback_order, warn = FALSE)
  if (!is.null(resolved)) {
    return(resolved)
  }

  available_cols
}

.hc_prepare_plot_heatmap_columns <- function(mat,
                                             cluster_columns = FALSE,
                                             plot_order = NULL,
                                             main_order = NULL,
                                             fallback_order = NULL,
                                             context = "heatmap") {
  mat <- base::as.matrix(mat)
  if (base::is.null(base::colnames(mat)) || base::ncol(mat) == 0) {
    return(list(
      mat = mat,
      col_order = base::character(0),
      col_dend = NULL,
      clustered = FALSE
    ))
  }

  if (!isTRUE(cluster_columns) || base::ncol(mat) <= 1) {
    selected_col_order <- .hc_select_heatmap_col_order(
      available_cols = base::colnames(mat),
      plot_order = plot_order,
      main_order = main_order,
      fallback_order = fallback_order,
      context = context
    )
    if (base::length(selected_col_order) > 0) {
      mat <- .hc_subset_matrix_cols_with_duplicates(mat, selected_col_order)
    }
    return(list(
      mat = mat,
      col_order = base::colnames(mat),
      col_dend = NULL,
      clustered = FALSE
    ))
  }

  mat_num <- mat
  mode(mat_num) <- "numeric"
  mat_num[!base::is.finite(mat_num)] <- 0

  clustered <- tryCatch(
    {
      hc <- stats::hclust(stats::dist(base::t(mat_num)), method = "complete")
      ord_idx <- hc$order
      ord <- base::colnames(mat_num)[ord_idx]
      list(
        mat = mat[, ord_idx, drop = FALSE],
        col_order = ord,
        col_dend = stats::as.dendrogram(hc),
        clustered = TRUE
      )
    },
    error = function(e) NULL
  )

  if (!base::is.null(clustered)) {
    return(clustered)
  }

  selected_col_order <- .hc_select_heatmap_col_order(
    available_cols = base::colnames(mat),
    plot_order = plot_order,
    main_order = main_order,
    fallback_order = fallback_order,
    context = context
  )
  if (base::length(selected_col_order) > 0) {
    mat <- .hc_subset_matrix_cols_with_duplicates(mat, selected_col_order)
  }
  warning(
    "Could not cluster heatmap columns for ", context, ". Falling back to explicit/stored order.",
    call. = FALSE
  )
  list(
    mat = mat,
    col_order = base::colnames(mat),
    col_dend = NULL,
    clustered = FALSE
  )
}

#' Helper Function To Plot TF Enrichment
#' @noRd

plot_TF <- function(hubs_df) {
  # get column order from module heatmap:
  heatmap_info <- .hc_heatmap_cache_info(hcobject[["integrated_output"]][["cluster_calc"]])
  co <- heatmap_info$col_order
  if (is.null(co) || length(co) == 0) {
    co <- base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]][, -base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]]), drop = FALSE])
  }

  colname_state <- new.env(parent = emptyenv())
  colname_state$final_colnames <- NULL
  hub_exp <- base::apply(hubs_df, 2, function(x) {
    genes <- x[!x == " "]
    if (base::length(genes) == 0) {
      return(NULL)
    }

    out <- base::lapply(base::seq_along(hcobject[["layers"]]), function(l) {
      exp <- hcobject[["data"]][[base::paste0("set", l, "_counts")]][genes, ]
      mean_exp <- base::lapply(base::unique(hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]), function(c) {
        samples <- base::rownames(hcobject[["data"]][[base::paste0("set", l, "_anno")]][hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]] == c, ])
        tmp <- dplyr::select(exp, samples) %>% base::rowMeans()
      }) %>%
        rlist::list.cbind() %>%
        base::as.data.frame()
      base::colnames(mean_exp) <- base::unique(hcobject[["data"]][[base::paste0("set", l, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]) %>% base::as.character()
      l_co <- co[co %in% base::colnames(mean_exp)]
      colname_state$final_colnames <- c(colname_state$final_colnames, l_co)

      mean_exp <- mean_exp[, l_co]

      return(mean_exp)
    }) %>% rlist::list.cbind()

    return(out)
  })
  title_plot <- base::names(hub_exp)[1]
  hub_exp <- hub_exp[[1]][stats::complete.cases(hub_exp[[1]]), ]
  base::colnames(hub_exp) <- colname_state$final_colnames


  p <- ComplexHeatmap::Heatmap(base::as.matrix(hub_exp),
    show_column_names = TRUE,
    show_row_names = TRUE,
    border = "black",
    col = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 7, name = "BrBG")))(50),
    cluster_columns = FALSE, cluster_rows = FALSE, rect_gp = grid::gpar(col = "black"),
    heatmap_legend_param = list(title = ""), width = grid::unit(5, "cm"), height = grid::unit(7, "cm"),
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8), column_title = title_plot
  )
  return(p)
}

#' Change Colour Transparency
#'
#' @param color String. Name of the colour.
#' @param alpha Integer between 1 and 100, giving the alpha value (saturation) of the target colour. Default is 100.
#' @noRd

makeTransparent <- function(color, alpha = 100) {
  newColor <- grDevices::col2rgb(color)
  transp <- base::apply(newColor, 2, function(rgbdata) {
    grDevices::rgb(
      red = rgbdata[1], green = rgbdata[2],
      blue = rgbdata[3], alpha = alpha, maxColorValue = 255
    )
  })
  return(transp)
}


#' Plot Boxplots From A List Of Dataframes As Input
#' @noRd

boxplot_from_list <- function(data, log_2, bool_plot) {
  plts <- list()

  for (x in base::seq_along(data)) {
    if (methods::is(data[[x]], "DataFrame")) {
      data[[x]] <- base::as.data.frame(data[[x]])
    } else if (base::is.matrix(data[[x]])) {
      data[[x]] <- base::as.data.frame(data[[x]])
    }

    if (!base::is.data.frame(data[[x]])) {
      stop("Your list needs to contain data frames only!")
    }

    plts[[x]] <- boxplot_from_df(data[[x]], it = hcobject[["layers_names"]][x], log_2 = log_2, bool_plot = bool_plot)
  }
  return(plts)
}


#' Plot Boxplots From A Dataframe As Input
#' @noRd

boxplot_from_df <- function(data, it = NULL, log_2, bool_plot) {
  # in case of large data set inform the user about the plots being split:
  if (base::ncol(data) > 50) {
    message("To avoid over-crowded plots, several plots will be created with up to 50 samples each.")
  }

  if (log_2 == TRUE) {
    data[data < 1] <- NA

    data <- base::log(data, base = 2)
  }

  # if there are less than 50 samples, the data can be plotted in one plot:
  if (base::ncol(data) <= 50) {
    bp_df <- base::data.frame(values = base::unlist(data), sample = base::rep(base::colnames(data), each = base::nrow(data)))

    p <- ggplot2::ggplot(bp_df, ggplot2::aes(x = sample, y = values)) +
      ggplot2::geom_boxplot(fill = "gray90", color = "deepskyblue4") +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(it) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = ggplot2::element_text(size = 10)
      )

    .hc_export_ggplot_file(
      file = .hc_output_file(base::paste0("Sample_distribution_bp_", it, ".pdf")),
      plot = p,
      width = 17,
      height = 7.8
    )

    return(p)

    # if there are more samples, the plots need to be split up to improve readability
  } else {
    # find the best suitable number of boxplots per plot to make the figures as balanced as possible
    # (e.g. when there are 51 columns, we would not want a figure with 50 boxplots and another figure with just one boxplot
    # but rather one with 26 boxplots and one with 25)
    n <- find_best_mod(nc = base::ncol(data), ref = 50)

    plts <- list()
    # plot each plot:
    for (x in base::seq_len((base::ncol(data) %/% n) + 1)) {
      start <- ((x - 1) * n) + 1

      end <- x * n

      if (start > base::ncol(data)) {
        break()
      }

      if (end > base::ncol(data)) {
        end <- base::ncol(data)
      }

      bp_df <- base::data.frame(values = base::unlist(data[, start:end]), sample = base::rep(base::colnames(data[, start:end]), each = base::nrow(data)))

      p <- ggplot2::ggplot(bp_df, ggplot2::aes(x = sample, y = values)) +
        ggplot2::geom_boxplot(fill = "gray90", color = "deepskyblue4") +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(it) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
          text = ggplot2::element_text(size = 10)
        )

      .hc_export_ggplot_file(
        file = .hc_output_file(base::paste0("Sample_distribution_bp_", it, "_", x, ".pdf")),
        plot = p,
        width = 10,
        height = 8
      )

      plts[[x]] <- p
    }
    return(plts)
  }
}


#' Plot Frequency Distributions From A List Of Dataframes As Input
#' @noRd

freqdist_plot_from_list <- function(data, log_2, bool_plot) {
  plts <- list()

  for (x in base::seq_along(data)) {
    if (methods::is(data[[x]], "DataFrame")) {
      data[[x]] <- base::as.data.frame(data[[x]])
    } else if (base::is.matrix(data[[x]])) {
      data[[x]] <- base::as.data.frame(data[[x]])
    }

    if (!base::is.data.frame(data[[x]])) {
      stop("Your list needs to contain data frames only!")
    }

    plts[[x]] <- freqdist_plot_from_df(data = data[[x]], log_2 = log_2, it = hcobject[["layers_names"]][x], bool_plot = bool_plot)
  }

  return(plts)
}


#' Plot Frequency Distributions From A Dataframe As Input
#' @noRd

freqdist_plot_from_df <- function(data, log_2, bool_plot, it = NULL) {
  if (base::ncol(data) > 42) {
    message("A high number of samples was detected. Frequency distributions are saved per sample.")
  }

  if (log_2 == TRUE) {
    data[data < 1] <- NA

    data <- base::log(data, base = 2)
  }

  plts <- list()
  out_path <- .hc_output_dir()

  for (x in base::seq_len(base::ncol(data))) {
    sample_name <- base::colnames(data)[x]
    vals <- data[[x]]
    xmax <- .hc_max_finite(vals)
    if (!is.finite(xmax) || xmax <= 0) {
      xmax <- 1
    }

    p <- ggplot2::ggplot(data[, x, drop = FALSE], ggplot2::aes(x = data[[x]])) +
      ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(count)), fill = "lightgray") +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = base::mean(data[[x]][stats::complete.cases(data[[x]])])),
        linetype = "dashed",
        linewidth = 0.6,
        color = "#FC4E07"
      ) +
      ggplot2::xlab(sample_name) +
      ggplot2::xlim(base::c(0, xmax)) +
      ggplot2::theme_bw() +
      ggplot2::theme(text = ggplot2::element_text(size = 10)) +
      ggplot2::ggtitle(it)

    safe_sample <- gsub("[^A-Za-z0-9._-]+", "_", sample_name)
    .hc_export_ggplot_file(
      file = base::file.path(out_path, base::paste0("Sample_distribution_freq_", it, "_", safe_sample, ".pdf")),
      plot = p,
      width = 7,
      height = 5
    )

    plts[[sample_name]] <- p
  }

  return(plts)
}


#' Plot List Of Plots
#' @noRd

plot_list_of_plots <- function(plts) {
  for (j in base::seq_along(plts)) {
    tmp <- plts[[j]]

    if (inherits(tmp, "ggplot")) {
      graphics::plot(tmp)
    } else {
      if (base::is.list(tmp)) {
        for (k in base::seq_along(tmp)) {
          if (inherits(tmp[[k]], "ggplot")) {
            graphics::plot(tmp[[k]])
          } else {
            .hc_display_object(patchwork::wrap_plots(tmp[[k]]))
          }
        }
      } else {
        .hc_display_object(patchwork::wrap_plots(tmp))
      }
    }
  }
}


#' Split A Number In Well-Balanced Way
#' @noRd

find_best_mod <- function(nc, ref) {
  mods <- NULL

  for (x in base::seq_len(ref)) {
    mods <- dplyr::bind_rows(mods, base::data.frame(x = x, mod = nc %% x))
  }

  mods <- mods[!base::duplicated(mods$mod), ]

  mods <- mods[mods$mod == base::max(mods$mod), "x"][1]

  return(mods)
}


#' Filter Integrated Network For Non-White Genes
#' @noRd

network_filt <- function() {
  gtc <- .hc_gene_to_cluster_impl()
  network <- hcobject[["integrated_output"]][["merged_net"]]
  del_v <- igraph::V(network)$name[!igraph::V(network)$name %in% gtc$gene]
  network <- igraph::delete.vertices(network, del_v)
  return(network)
}


#' For each sample calculates the mean expression per cluster
#' @param set An integer. Number of the dataset, for which samples the means should be calculated.
#' @noRd

sample_wise_cluster_expression <- function(set) {
  gtc <- .hc_gene_to_cluster_impl()
  counts <- hcobject[["data"]][[base::paste0("set", set, "_counts")]]
  modules <- base::unique(gtc$color[!gtc$color == "white"])

  # Modules are built on the *integrated* network, so they routinely contain
  # genes that were not measured in this particular layer. `counts` is a matrix
  # in the S4 pipeline, where indexing by an unknown row name is a hard
  # "subscript out of bounds" error (a data.frame would have yielded NA rows,
  # which is what the complete.cases() guard below was written for). Restrict to
  # the genes actually present, and keep the matrix 2-dimensional so modules
  # with a single remaining gene do not collapse to a vector.
  available <- base::rownames(counts)
  n_samples <- base::ncol(counts)

  cluster_means <- base::lapply(modules, function(c) {
    genes <- dplyr::filter(gtc, color == c) %>% dplyr::pull(., "gene")
    genes <- base::intersect(genes, available)
    if (base::length(genes) == 0) {
      return(base::rep(NA_real_, n_samples))
    }
    tmp <- counts[genes, , drop = FALSE]
    tmp <- tmp[stats::complete.cases(tmp), , drop = FALSE]
    if (base::nrow(tmp) == 0) {
      return(base::rep(NA_real_, n_samples))
    }
    base::colMeans(tmp)
  }) %>%
    rlist::list.rbind() %>%
    base::as.data.frame()

  base::rownames(cluster_means) <- modules
  base::colnames(cluster_means) <- base::colnames(counts)

  return(cluster_means)
}


#' Hub Node Detection
#'
#' Subroutine to .hc_find_hubs_driver().
#' @noRd


hub_node_detection <- function(cluster, top, save, tree_layout, TF_only, plot,
                               label = NULL) {
  # `cluster` selects the subnetwork by colour; `label` is what the user sees
  # (module label such as "M2.1"). Falls back to the colour when not supplied.
  if (base::is.null(label) || base::length(label) == 0 ||
    base::is.na(label[[1]]) || !base::nzchar(base::as.character(label[[1]]))) {
    label <- cluster
  }
  label <- base::as.character(label[[1]])
  # extract chosen cluster as an isolated network:
  g <- cluster_to_network(cluster = cluster)
  # determine hub nodes:
  hub_out <- get_hub_nodes(network = g, top = top, TF_only = TF_only)

  # rank_df <- hub_out$rank_df
  if (base::nrow(hub_out$rank_df) == 0) {
    return(hub_out)
  }

  # vertex size: non-hubs small & uniform, hubs scaled by their combined
  # centrality rank so the strongest hub reads as the largest node.
  rank_lookup <- stats::setNames(hub_out$rank_df$sum, hub_out$rank_df$node)
  hub_sums <- rank_lookup[hub_out$hub_nodes]
  hub_sizes <- if (base::length(hub_sums) > 0 &&
    base::diff(base::range(hub_sums, na.rm = TRUE)) > 0) {
    scales::rescale(hub_sums, to = base::c(6, 14))
  } else {
    base::rep(9, base::length(hub_sums))
  }
  base::names(hub_sizes) <- hub_out$hub_nodes
  vertex_size <- base::vapply(igraph::V(g)$name, function(node) {
    if (node %in% hub_out$hub_nodes) hub_sizes[[node]] else 2.5
  }, FUN.VALUE = base::numeric(1))

  # layout:
  if (cluster == "all") {
    l <- hcobject[["integrated_output"]][["cluster_calc"]][["layout"]]
  } else {
    l <- .hc_hub_network_layout(g, tree_layout = tree_layout)
    base::rownames(l) <- igraph::V(g)$name
  }


  # set some new node attributes:
  igraph::V(g)$size <- vertex_size
  igraph::V(g)$label <- NA
  igraph::V(g)$color <- hub_out$colour_df$colour
  # fade the non-hub background so the (fully opaque) hub nodes stand out.
  is_hub <- igraph::V(g)$name %in% hub_out$hub_nodes
  if (base::any(!is_hub)) {
    igraph::V(g)$color[!is_hub] <- grDevices::adjustcolor(
      igraph::V(g)$color[!is_hub],
      alpha.f = 0.5
    )
  }
  # plot network: modern ggplot2 style by default, with a graceful fallback to
  # the classic base-graphics renderer when ggrepel is unavailable or the user
  # opts out via options(hcocena.hub_network_style = "classic").
  use_modern <- !base::identical(
    base::getOption("hcocena.hub_network_style", "modern"), "classic"
  ) && .hc_has_modern_hub_plot_deps()
  if (use_modern) {
    .hc_hub_network_modern_plot(
      network = g,
      hub_nodes = hub_out$hub_nodes,
      gene_ranks = base::seq_along(hub_out$hub_nodes),
      layout = l,
      centrality = stats::setNames(hub_out$rank_df$sum, hub_out$rank_df$node),
      title = base::c(label, top),
      save = save,
      plot = plot
    )
  } else {
    network_with_labels(
      network = g,
      gene_labels = hub_out$hub_nodes,
      gene_ranks = base::seq_along(hub_out$hub_nodes),
      l = l,
      label_offset = 10,
      title = base::c(label, top),
      save = save,
      plot = plot
    )
  }


  return(hub_out)
}


#' Cluster To Network
#'
#' Extracts a cluster as a standalone network.
#' @noRd

cluster_to_network <- function(cluster) {
  gtc <- .hc_gene_to_cluster_impl() %>% dplyr::filter(., color == cluster)
  g <- hcobject[["integrated_output"]][["merged_net"]]
  g <- igraph::delete_vertices(g, igraph::V(g)$name[!igraph::V(g)$name %in% gtc$gene])

  return(g)
}


#' Get Hub Nodes
#'
#' Subroutine to .hc_find_hubs_driver().
#' @noRd

.hc_hub_tf_filter_genes <- function(TF_only = FALSE) {
  if (base::identical(TF_only, FALSE)) {
    return(NULL)
  }
  if (!base::is.character(TF_only) ||
    base::length(TF_only) != 1 ||
    base::is.na(TF_only) ||
    !base::nzchar(base::trimws(TF_only))) {
    stop(
      "`TF_only` must be FALSE, \"all\", or one transcription-factor category.",
      call. = FALSE
    )
  }
  TF_only <- base::trimws(TF_only)

  tf_table <- hcobject[["supplementary_data"]][["TF"]]
  if (base::is.null(tf_table) ||
    base::length(tf_table) == 0 ||
    base::is.null(base::dim(tf_table)) ||
    base::nrow(tf_table) == 0 ||
    base::ncol(tf_table) == 0) {
    stop(
      "`TF_only = \"",
      TF_only,
      "\"` requires a non-empty transcription-factor reference in ",
      "`hc@supplementary$TF`.",
      call. = FALSE
    )
  }
  tf_table <- base::as.data.frame(tf_table, stringsAsFactors = FALSE)

  organism <- hcobject[["global_settings"]][["organism"]]
  if (base::is.null(organism) ||
    base::length(organism) != 1 ||
    base::is.na(organism[[1]]) ||
    !base::nzchar(base::trimws(base::as.character(organism[[1]])))) {
    stop(
      "A single non-empty organism setting is required for `TF_only` filtering.",
      call. = FALSE
    )
  }
  organism <- base::trimws(base::as.character(organism[[1]]))
  organism_columns <- base::grep(
    organism,
    base::colnames(tf_table),
    ignore.case = TRUE,
    value = TRUE
  )
  exact_columns <- base::colnames(tf_table)[
    base::tolower(base::colnames(tf_table)) == base::tolower(organism)
  ]
  if (base::length(exact_columns) == 1) {
    organism_column <- exact_columns[[1]]
  } else if (base::length(organism_columns) == 1) {
    organism_column <- organism_columns[[1]]
  } else if (base::length(organism_columns) == 0) {
    stop(
      "No transcription-factor reference column matches organism `",
      organism,
      "`.",
      call. = FALSE
    )
  } else {
    stop(
      "Multiple transcription-factor reference columns match organism `",
      organism,
      "`: ",
      base::paste(organism_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  tf_rows <- tf_table
  if (!base::identical(TF_only, "all")) {
    if (base::ncol(tf_table) < 2) {
      stop(
        "A categorized `TF_only` filter requires a category column in the ",
        "transcription-factor reference.",
        call. = FALSE
      )
    }
    category_column <- base::colnames(tf_table)[base::ncol(tf_table)]
    available_categories <- base::unique(base::as.character(tf_table[[category_column]]))
    available_categories <- available_categories[
      !base::is.na(available_categories) & base::nzchar(available_categories)
    ]
    if (!(TF_only %in% available_categories)) {
      stop(
        "Unknown `TF_only` category `",
        TF_only,
        "`. Available categories are: ",
        base::paste(base::sort(available_categories), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    tf_rows <- tf_table[
      base::as.character(tf_table[[category_column]]) == TF_only,
      ,
      drop = FALSE
    ]
  }

  genes <- base::trimws(base::as.character(tf_rows[[organism_column]]))
  genes <- base::unique(genes[!base::is.na(genes) & base::nzchar(genes)])
  if (base::length(genes) == 0) {
    stop(
      "The requested transcription-factor filter contains no genes for organism `",
      organism,
      "`.",
      call. = FALSE
    )
  }
  genes
}


get_hub_nodes <- function(network = hcobject[["integrated_output"]][["merged_net"]],
                          top = 10,
                          TF_only = FALSE) {
  tf_genes <- .hc_hub_tf_filter_genes(TF_only)
  rank_df <- combined_centrality(network = network)
  if (!"node" %in% base::colnames(rank_df) ||
    base::all(base::is.na(rank_df$node) | !base::nzchar(base::as.character(rank_df$node)))) {
    rank_df$node <- base::rownames(rank_df)
  }

  if (!base::is.null(tf_genes)) {
    rank_df <- dplyr::filter(rank_df, node %in% tf_genes)
  }

  if (top > base::nrow(rank_df)) {
    hub_nodes <- rank_df %>% dplyr::pull(., "node")
  } else {
    hub_nodes <- rank_df[base::seq_len(top), ] %>% dplyr::pull(., "node")
  }

  if (base::nrow(rank_df) == 0) {
    message("No hubs found for this cluster.")
    return(list(rank_df = rank_df, colour_df = NULL, hub_nodes = hub_nodes))
  }
  colour_df <- centrality_colours(rank_df = rank_df, network = network)

  return(list(rank_df = rank_df, colour_df = colour_df, hub_nodes = hub_nodes))
}


#' Combined Centrality Measure
#'
#' Subroutine to .hc_find_hubs_driver().
#' @noRd

combined_centrality <- function(network) {
  dc <- weighted_DC(network)
  cc <- weighted_CC(network)
  bc <- weighted_BC(network)
  node_names <- igraph::V(network)$name
  if (base::is.null(node_names) || base::length(node_names) != igraph::vcount(network)) {
    node_names <- igraph::V(network) %>% base::as.character()
  }

  rank_df <- base::data.frame(
    dc = base::rank(base::unname(dc), ties.method = "average"),
    cc = base::rank(base::unname(cc), ties.method = "average"),
    bc = base::rank(base::unname(bc), ties.method = "average"),
    node = base::as.character(node_names),
    stringsAsFactors = FALSE
  )

  rank_df$sum <- base::apply(rank_df[, c("dc", "cc", "bc"), drop = FALSE], 1, base::sum)
  rank_df$id <- igraph::V(network) %>% base::as.character()

  # order based on highest rank (strongest hub candidates):
  rank_df <- rank_df[base::order(rank_df$sum, decreasing = TRUE), ]
  base::rownames(rank_df) <- rank_df$node
  return(rank_df)
}

.hc_graph_edge_weights <- function(network, default = 1) {
  edge_n <- igraph::ecount(network)
  if (edge_n == 0) {
    return(base::numeric())
  }
  w <- igraph::edge_attr(network, "weight", index = igraph::E(network))
  if (base::is.null(w)) {
    return(base::rep(default, edge_n))
  }
  w <- .hc_as_numeric_safely(w)
  if (base::length(w) != edge_n) {
    w <- base::rep(default, edge_n)
  }
  w[!base::is.finite(w)] <- default
  w
}

.hc_graph_distance_weights <- function(network,
                                       weights = .hc_graph_edge_weights(network),
                                       eps = 1e-8) {
  if (base::length(weights) == 0) {
    return(NULL)
  }

  if (base::all(weights >= 0, na.rm = TRUE) && base::all(weights <= 1, na.rm = TRUE)) {
    dist_w <- 1 - weights
  } else {
    w_min <- .hc_min_finite(weights)
    w_max <- .hc_max_finite(weights)
    if (!base::is.finite(w_min) || !base::is.finite(w_max) || w_min == w_max) {
      dist_w <- base::rep(1, base::length(weights))
    } else {
      scaled_w <- (weights - w_min) / (w_max - w_min)
      dist_w <- 1 - scaled_w
    }
  }

  dist_w[!base::is.finite(dist_w)] <- 1
  base::pmax(dist_w, eps)
}


#' Hub Network Layout
#'
#' Layout for the per-cluster hub networks. Prefers graphlayouts' stress layout
#' (deterministic, cleanly separates components, avoids the hairball/cutoff
#' issues of `layout.lgl`), falling back to a weighted Fruchterman-Reingold and
#' finally `layout.lgl` so behaviour degrades gracefully when graphlayouts is
#' unavailable or a layout cannot be computed.
#' @noRd

.hc_hub_network_layout <- function(g, tree_layout = FALSE) {
  if (base::isTRUE(tree_layout)) {
    return(igraph::layout_as_tree(g))
  }
  weights <- .hc_graph_edge_weights(g)
  l <- .hc_with_seed(1L, tryCatch(
    {
      if (base::requireNamespace("graphlayouts", quietly = TRUE)) {
        graphlayouts::layout_with_stress(g)
      } else {
        igraph::layout_with_fr(g, weights = weights)
      }
    },
    error = function(e) NULL
  ))
  if (base::is.null(l) || !base::is.matrix(l) || base::nrow(l) != igraph::vcount(g)) {
    l <- igraph::layout.lgl(g)
  }
  l
}


#' Rescale A Network Layout To A Stable Coordinate Range
#'
#' Maps layout coordinates into a `[0, target]` box while preserving the aspect
#' ratio (the larger of the two spans is scaled to `target`). This decouples
#' downstream label placement -- which uses an absolute offset -- from the
#' native coordinate scale of the layout backend, so stress/FR layouts are no
#' longer squished into a thin strip the way a fixed offset against their small
#' coordinate range would cause.
#' @noRd

.hc_rescale_layout <- function(l, target = 100) {
  l <- base::as.matrix(l)
  if (base::nrow(l) == 0) {
    return(l)
  }
  rng_x <- base::range(l[, 1], na.rm = TRUE)
  rng_y <- base::range(l[, 2], na.rm = TRUE)
  span <- base::max(base::diff(rng_x), base::diff(rng_y))
  if (!base::is.finite(span) || span <= 0) {
    return(l)
  }
  scale <- target / span
  l[, 1] <- (l[, 1] - rng_x[1]) * scale
  l[, 2] <- (l[, 2] - rng_y[1]) * scale
  l
}


#' Weighted Degree Centrality
#'
#' Caclulates the weighted degree centrality of a node. Modified to be weighted from https://doi.org/10.1155/2019/9728742.
#' @noRd

weighted_DC <- function(network) {
  message("Calculating weighted degree centrality.")
  if (igraph::vcount(network) == 0) {
    return(base::numeric())
  }

  weights <- .hc_graph_edge_weights(network)
  total_weight <- base::sum(weights)
  if (!base::is.finite(total_weight) || total_weight <= 0) {
    total_weight <- base::max(1, igraph::ecount(network))
    weights <- if (igraph::ecount(network) > 0) base::rep(1, igraph::ecount(network)) else base::numeric()
  }

  strengths <- igraph::strength(
    graph = network,
    vids = igraph::V(network),
    mode = "all",
    loops = FALSE,
    weights = if (base::length(weights) == 0) NULL else weights
  )
  out <- strengths / total_weight
  stats::setNames(base::as.numeric(out), base::as.character(igraph::V(network)))
}


#' Weighted Closeness Centrality
#'
#' Caclulates the weighted closeness centrality of a node. Modified to be weighted from https://doi.org/10.1155/2019/9728742.
#' @noRd

weighted_CC <- function(network) {
  message("Calculating weighted closeness centrality.")
  if (igraph::vcount(network) == 0) {
    return(base::numeric())
  }
  if (igraph::vcount(network) == 1) {
    return(stats::setNames(0, base::as.character(igraph::V(network))))
  }

  if (igraph::count_components(network) > 1) {
    message("Module consists of disconnected components. Closeness is computed on reachable vertices only.")
  }

  dist_weights <- .hc_graph_distance_weights(network)
  out <- igraph::closeness(
    graph = network,
    vids = igraph::V(network),
    mode = "all",
    weights = dist_weights,
    normalized = FALSE
  )
  out[!base::is.finite(out)] <- 0
  stats::setNames(base::as.numeric(out), base::as.character(igraph::V(network)))
}


#' Weighted Betweenness Centrality
#'
#' Caclulates the weighted betweenness centrality of a node. Modified to be weighted from https://doi.org/10.1155/2019/9728742.
#' @noRd

weighted_BC <- function(network) {
  message("Calculating weighted betweenness centrality.")
  if (igraph::vcount(network) == 0) {
    return(base::numeric())
  }
  if (igraph::ecount(network) == 0 || igraph::vcount(network) <= 2) {
    out <- base::rep(0, igraph::vcount(network))
    return(stats::setNames(out, base::as.character(igraph::V(network))))
  }

  dist_weights <- .hc_graph_distance_weights(network)
  out <- igraph::betweenness(
    graph = network,
    v = igraph::V(network),
    directed = FALSE,
    weights = dist_weights,
    normalized = FALSE
  )
  out[!base::is.finite(out)] <- 0
  stats::setNames(base::as.numeric(out), base::as.character(igraph::V(network)))
}


#' Colours Based On Centrality
#'
#' Subroutine to .hc_find_hubs_driver().
#' @noRd

centrality_colours <- function(rank_df, network) {
  mypalette <- grDevices::colorRampPalette(base::c("#FEE8C8", "#FC8D59", "#B30000"))
  colour_df <- base::data.frame(
    sum = rank_df$sum,
    id = rank_df$id,
    colour = (mypalette(base::max(rank_df$sum))[rank_df$sum])
  )
  colour_df <- colour_df[base::match(igraph::V(network), colour_df$id), ]
  return(colour_df)
}


#' Modern Hub-Network Renderer Dependencies
#'
#' TRUE when the optional packages for the ggplot2 hub-network style are present.
#' @noRd

.hc_has_modern_hub_plot_deps <- function() {
  base::requireNamespace("ggplot2", quietly = TRUE) &&
    base::requireNamespace("ggrepel", quietly = TRUE)
}


#' Plot A Hub Network (Modern ggplot2 Style)
#'
#' Renders the per-cluster hub network as a clean ggplot2 figure: the full
#' network sits faintly in the background, the hub genes are drawn as bright
#' nodes sized and coloured by their combined centrality, and only the hubs are
#' labelled (directly at the node, with a white halo and short repelled
#' connectors via ggrepel). `coord_equal()` keeps the network undistorted.
#' Falls back to [network_with_labels()] is handled by the caller.
#' @noRd

.hc_hub_network_modern_plot <- function(network,
                                        hub_nodes,
                                        gene_ranks,
                                        layout,
                                        centrality,
                                        title,
                                        save,
                                        plot,
                                        label_fontsize = 3.6,
                                        node_size_range = base::c(3, 8)) {
  plot_title <- base::paste0("cluster '", title[1], "' hub genes [top ", title[2], "]")

  node_names <- igraph::V(network)$name
  coords <- base::as.matrix(layout)
  if (!base::is.null(base::rownames(coords)) &&
    base::all(node_names %in% base::rownames(coords))) {
    coords <- coords[node_names, , drop = FALSE]
  }

  is_hub <- node_names %in% hub_nodes
  nodes_df <- base::data.frame(
    x = coords[, 1],
    y = coords[, 2],
    cent = base::as.numeric(centrality[node_names]),
    is_hub = is_hub,
    label = node_names,
    stringsAsFactors = FALSE
  )
  hubs_df <- nodes_df[nodes_df$is_hub, , drop = FALSE]
  bg_df <- nodes_df[!nodes_df$is_hub, , drop = FALSE]

  el <- igraph::as_edgelist(network, names = FALSE)
  edges_df <- if (base::nrow(el) > 0) {
    base::data.frame(
      x = coords[el[, 1], 1], y = coords[el[, 1], 2],
      xend = coords[el[, 2], 1], yend = coords[el[, 2], 2]
    )
  } else {
    base::data.frame(x = base::numeric(0), y = base::numeric(0),
      xend = base::numeric(0), yend = base::numeric(0))
  }

  p <- ggplot2::ggplot()
  if (base::nrow(edges_df) > 0) {
    p <- p + ggplot2::geom_segment(
      data = edges_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = "grey80", alpha = 0.18, linewidth = 0.15
    )
  }
  if (base::nrow(bg_df) > 0) {
    p <- p + ggplot2::geom_point(
      data = bg_df, ggplot2::aes(x = x, y = y),
      colour = "grey75", alpha = 0.45, size = 0.9
    )
  }
  p <- p +
    ggplot2::geom_point(
      data = hubs_df,
      ggplot2::aes(x = x, y = y, fill = cent, size = cent),
      shape = 21, colour = "white", stroke = 0.9
    ) +
    ggplot2::scale_fill_viridis_c(option = "rocket", direction = -1, end = 0.92, name = "centrality") +
    ggplot2::scale_size(range = node_size_range, guide = "none") +
    ggrepel::geom_text_repel(
      data = hubs_df,
      ggplot2::aes(x = x, y = y, label = label),
      fontface = "bold", size = label_fontsize, colour = "grey15",
      bg.color = "white", bg.r = 0.18,
      box.padding = 0.6, point.padding = 0.3, min.segment.length = 0,
      segment.colour = "grey55", segment.size = 0.3,
      max.overlaps = Inf, seed = 7
    ) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = plot_title) +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
      plot.margin = ggplot2::margin(12, 18, 12, 18),
      legend.position = "right"
    )

  if (isTRUE(save)) {
    .hc_ggsave_pdf_png(
      filename = .hc_output_file(base::paste0("Hub_genes_", title[1], "_module_network.pdf")),
      plot = p, width = 10, height = 9, bg = "white"
    )
  }
  if (isTRUE(plot)) {
    base::print(p)
  }
  base::invisible(p)
}


#' Plots A Network With Labels
#'
#' Subroutine to .hc_find_hubs_driver().
#' @noRd

network_with_labels <- function(network, gene_labels, gene_ranks, l, label_offset, title, save, plot) {
  plot_title <- base::paste0("cluster '", title[1], "' hub genes [top ", title[2], "]")

  # Normalise the layout to a stable coordinate scale and derive the side-label
  # offset from the node spread. A fixed absolute offset squished small-scale
  # layouts (stress/FR) into a thin vertical strip because it dwarfed the actual
  # network extent; a proportional offset keeps a sensible aspect ratio for any
  # layout backend.
  l <- .hc_rescale_layout(l)
  label_offset <- 0.55 * base::diff(base::range(l[, 1]))

  new_genes_df <- base::data.frame(name = base::paste0("label_", base::seq_along(gene_labels)), label = gene_labels %>% base::as.character())
  new_indeces <- base::match(new_genes_df$label, igraph::get.vertex.attribute(network)$name)
  new_genes_df$coords_x <- l[new_indeces, 1]
  new_genes_df$coords_y <- l[new_indeces, 2]
  # Split the labels into the left/right (and up/down) columns around the median
  # of the *hub* positions instead of all nodes, so the two columns stay roughly
  # balanced even when the hubs cluster on one side of the layout.
  mean_coord_x <- stats::median(new_genes_df$coords_x, na.rm = TRUE)
  mean_coord_y <- stats::median(new_genes_df$coords_y, na.rm = TRUE)
  left_up <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y >= mean_coord_y) %>% dplyr::pull(., label)
  left_down <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y < mean_coord_y) %>% dplyr::pull(., label)
  right_up <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y >= mean_coord_y) %>% dplyr::pull(., label)
  right_down <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y < mean_coord_y) %>% dplyr::pull(., label)

  new_position_l <- base::matrix(base::cbind(
    base::rep(base::ceiling(base::min(l[, 1])) - label_offset, base::length(left_up) + base::length(left_down)),
    base::seq(
      from = base::ceiling(base::max(l[, 2])), to = base::ceiling(base::min(l[, 2])),
      length.out = base::length(left_up) + base::length(left_down)
    )
  ), ncol = 2)
  new_position_r <- base::matrix(base::cbind(
    base::rep(base::ceiling(base::max(l[, 1])) + label_offset, base::length(right_up) + base::length(right_down)),
    base::seq(
      from = base::ceiling(base::max(l[, 2])), to = base::ceiling(base::min(l[, 2])),
      length.out = base::length(right_up) + base::length(right_down)
    )
  ), ncol = 2)
  new_position <- base::rbind(new_position_l, new_position_r)

  base::colnames(new_position) <- base::colnames(l)
  l2 <- base::rbind(l, new_position) %>% base::as.matrix()
  new_genes_df_l <- new_genes_df[new_genes_df$label %in% base::c(left_up, left_down), ]
  new_genes_df_l <- new_genes_df_l[base::order(new_genes_df_l$coords_y, decreasing = TRUE), ]
  new_genes_df_r <- new_genes_df[new_genes_df$label %in% base::c(right_up, right_down), ]
  new_genes_df_r <- new_genes_df_r[base::order(new_genes_df_r$coords_y, decreasing = TRUE), ]
  new_genes_df <- base::rbind(new_genes_df_l, new_genes_df_r)

  network2 <- igraph::add.vertices(network,
    nv = base::length(new_genes_df$name),
    attr = list(name = new_genes_df$name)
  )
  new_labels <- base::lapply(igraph::get.vertex.attribute(network2)$name, function(x) {
    if (x %in% igraph::get.vertex.attribute(network)$name) {
      NA
    } else {
      tmp_label <- dplyr::filter(new_genes_df, name == x) %>%
        dplyr::pull(., "label")
      tmp_rank <- gene_ranks[base::match(tmp_label, gene_labels)]
      base::paste0(tmp_label, " [", tmp_rank, "]")
    }
  }) %>% base::unlist()

  igraph::V(network2)$label <- new_labels


  name_to_id <- base::data.frame(name = igraph::V(network2)$name, id = base::seq_along(igraph::V(network2)$name))


  nti1 <- name_to_id[name_to_id$name %in% new_genes_df[, 1], ]
  base::rownames(nti1) <- nti1$name
  nti1 <- nti1[new_genes_df[, 1], ]

  nti2 <- name_to_id[name_to_id$name %in% new_genes_df[, 2], ]
  base::rownames(nti2) <- nti2$name
  nti2 <- nti2[new_genes_df[, 2], ]

  new_edges <- base::matrix(base::c(nti1$id, nti2$id), ncol = 2, byrow = FALSE)

  new_edges <- base::as.vector(base::t(new_edges))

  network2 <- igraph::add.edges(network2, new_edges)

  new_edge_color <- base::apply(igraph::get.edgelist(network2), 1, function(x) {
    if (x[2] %in% base::as.character(new_genes_df$name)) {
      "black"
    } else {
      "lightgrey"
    }
  })

  # edge aesthetics: thin straight leader lines to the labels, co-expression
  # edges scaled by weight and gently curved to reduce overplotting.
  el <- igraph::get.edgelist(network2)
  is_leader_edge <- el[, 2] %in% base::as.character(new_genes_df$name)
  edge_weights <- .hc_graph_edge_weights(network2, default = NA_real_)
  net_w <- edge_weights[!is_leader_edge]
  new_edge_width <- base::rep(0.6, base::nrow(el))
  if (base::any(!is_leader_edge) &&
    base::diff(base::range(net_w, na.rm = TRUE)) > 0) {
    new_edge_width[!is_leader_edge] <- scales::rescale(net_w, to = base::c(0.3, 2.2))
  } else if (base::any(!is_leader_edge)) {
    new_edge_width[!is_leader_edge] <- 0.8
  }
  new_edge_curved <- base::ifelse(is_leader_edge, 0, 0.12)

  vertex_shape <- base::lapply(igraph::V(network2)$name, function(x) {
    if (x %in% igraph::V(network)$name) {
      "circle"
    } else {
      "none"
    }
  }) %>% base::unlist()

  new_label_color <- base::lapply(igraph::get.vertex.attribute(network2)$name, function(x) {
    if (x %in% new_genes_df$name) {
      tmp_gene_n <- dplyr::filter(new_genes_df, name == x) %>%
        dplyr::pull(., label)
      tfs <- hcobject[["supplementary_data"]][["TF"]][, base::grep(base::colnames(hcobject[["supplementary_data"]][["TF"]]), pattern = hcobject[["global_settings"]][["organism"]], ignore.case = TRUE)]
      if (tmp_gene_n %in% tfs) {
        "black"
      } else {
        igraph::get.vertex.attribute(network2, name = "color", index = tmp_gene_n)
      }
    } else {
      NA
    }
  }) %>% base::unlist()

  name_to_size <- base::data.frame(name = igraph::V(network)$name, size = igraph::V(network)$size)

  vertex_size <- base::lapply(igraph::V(network2)$name, function(x) {
    if (x %in% igraph::V(network)$name) {
      dplyr::filter(name_to_size, name == x) %>% dplyr::pull(., "size")
    } else {
      0
    }
  }) %>% base::unlist()

  # outline the hub nodes (the labelled genes) so they stand out from the faded
  # background; leave the rest borderless to avoid speckling dense networks.
  is_hub_node <- igraph::V(network2)$name %in% base::as.character(gene_labels)
  vertex_frame_color <- base::ifelse(is_hub_node, "black", NA)

  # Normalise l2 into a centred square and plot with rescale = FALSE / asp = 1.
  # igraph's default per-axis rescale stretches each axis to [-1, 1]
  # independently; because the side labels widen only the x-range, that made
  # the network render roughly twice as tall as wide. Scaling both axes by the
  # same factor preserves the network's true aspect ratio.
  # Centre the frame on the network-node centroid (not the bounding-box midpoint)
  # so the dense core sits in the middle and a few peripheral nodes no longer pull
  # the composition off to one side. Scale both axes by the same factor to keep
  # the network undistorted (plotted with rescale = FALSE / asp = 1 below).
  net_rows <- base::seq_len(igraph::vcount(network))
  l2_cx <- stats::median(l2[net_rows, 1], na.rm = TRUE)
  l2_cy <- stats::median(l2[net_rows, 2], na.rm = TRUE)
  l2_scale <- base::max(base::diff(base::range(l2[, 1])), base::diff(base::range(l2[, 2]))) / 2
  if (base::is.finite(l2_scale) && l2_scale > 0) {
    l2[, 1] <- (l2[, 1] - l2_cx) / l2_scale
    l2[, 2] <- (l2[, 2] - l2_cy) / l2_scale
  }
  # window reaches the farthest label on each axis (extra horizontal room for the
  # label text) so nothing is clipped while the core stays centred.
  x_reach <- base::max(base::abs(l2[, 1]), na.rm = TRUE)
  y_reach <- base::max(base::abs(l2[, 2]), na.rm = TRUE)
  plot_xlim <- base::c(-x_reach, x_reach) * 1.18
  plot_ylim <- base::c(-y_reach, y_reach) * 1.08
  if (!base::all(base::is.finite(plot_xlim)) || base::diff(plot_xlim) <= 0) {
    plot_xlim <- base::c(-1.2, 1.2)
  }
  if (!base::all(base::is.finite(plot_ylim)) || base::diff(plot_ylim) <= 0) {
    plot_ylim <- base::c(-1.2, 1.2)
  }

  if (save == TRUE) {
    .hc_export_single_page_plot(
      file = .hc_output_file(base::paste0("Hub_genes_", title[1], "_module_network.pdf")),
      width = 20,
      height = 15,
      draw_fun = function() {
        igraph::plot.igraph(
          network2,
          vertex.size = vertex_size,
          vertex.label = new_labels,
          vertex.label.cex = 1.5,
          layout = l2,
          rescale = FALSE,
          xlim = plot_xlim,
          ylim = plot_ylim,
          asp = 1,
          vertex.label.dist = 1,
          vertex.shape = vertex_shape,
          vertex.frame.color = vertex_frame_color,
          edge.color = new_edge_color,
          edge.width = new_edge_width,
          edge.curved = new_edge_curved,
          vertex.label.color = new_label_color
        )
        graphics::title(plot_title, cex.main = 3)
      }
    )
  }
  if (plot) {
    igraph::plot.igraph(network2,
      vertex.size = vertex_size, vertex.label = new_labels, vertex.label.cex = 0.75,
      layout = l2, main = plot_title, vertex.label.dist = 1, vertex.shape = vertex_shape,
      rescale = FALSE, xlim = plot_xlim, ylim = plot_ylim, asp = 1,
      vertex.frame.color = vertex_frame_color, edge.color = new_edge_color,
      edge.width = new_edge_width, edge.curved = new_edge_curved,
      vertex.label.color = new_label_color
    )
  }
}


#' Run All Clustering Algorithms
#'
#' Runs all clustering algorithms on the network and returns the clustering.
#' @noRd

run_all_cluster_algos <- function() {
  output <- list()

  network <- hcobject[["integrated_output"]][["merged_net"]]

  cluster_algo_list <- base::c(
    "cluster_louvain",
    "cluster_fast_greedy",
    "cluster_infomap",
    "cluster_walktrap",
    "cluster_label_prop",
    "cluster_leiden"
  )


  color.cluster <- get_cluster_colours()


  current_algo <- NULL

  for (c in base::unique(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]$color)) {
    genes <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], color == c) %>%
      dplyr::pull(., "gene_n") %>%
      base::strsplit(., split = ",") %>%
      base::unlist(.)
    current_algo <- base::rbind(current_algo, base::data.frame(gene = genes, cluster = base::rep(c, base::length(genes))))
  }

  current_algo$cluster <- base::as.factor(current_algo$cluster)

  output[[hcobject[["global_settings"]][["chosen_clustering_algo"]]]] <- current_algo

  for (a in cluster_algo_list[!cluster_algo_list == hcobject[["global_settings"]][["chosen_clustering_algo"]]]) {
    partition_df <- .hc_with_seed(1L, tryCatch(
      {
        if (a == "cluster_leiden") {
          partition <- leidenAlg::leiden.community(graph = network)
          partition_df <- base::data.frame(gene = partition$names, cluster = base::as.numeric(partition$membership))
          partition_df$cluster <- partition_df$cluster + 1
        } else {
          cfg <- base::getExportedValue("igraph", a)(network)
          partition_df <- base::data.frame(gene = igraph::get.vertex.attribute(network)$name, cluster = cfg$membership)
        }
        partition_df
      },
      error = function(e) {
        warning(
          "Skipping `", a, "` in `run_all_cluster_algos()`: ",
          base::conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    ))

    if (base::is.null(partition_df) || base::nrow(partition_df) == 0) {
      next
    }

    partition_df$cluster <- base::lapply(partition_df$cluster, function(x) {
      color.cluster[x]
    }) %>% base::unlist()


    partition_df_table <- base::table(partition_df$cluster) %>%
      base::as.data.frame()

    partition_df <- partition_df[partition_df$cluster %in% (dplyr::filter(partition_df_table, Freq >= hcobject[["global_settings"]][["min_nodes_number_for_cluster"]]) %>%
      dplyr::pull(., "Var1")), ]

    partition_df$cluster <- base::as.factor(partition_df$cluster)

    output[[a]] <- partition_df
  }
  if (base::length(output) <= 1) {
    warning("No alternative clustering algorithm completed in `run_all_cluster_algos()`.", call. = FALSE)
  }
  return(output)
}


#' Plot PCA of top most variant
#'
#' Subroutine to .hc_PCA_algo_compare_driver().
#' @noRd

plot_PCA_topvar <- function(PCA_save_folder, cols = cols) {
  plotlist <- list()
  pca_list <- list()
  for (i in base::seq_along(hcobject[["layers"]])) {
    pca_input <- .hc_pca_prepare_expression(
      hcobject[["layer_specific_outputs"]][[base::paste0("set", i)]][["part1"]][["topvar"]],
      layer_label = hcobject[["layers_names"]][i],
      scale_features = TRUE
    )
    pca <- stats::prcomp(pca_input, scale. = TRUE)
    pca.var <- pca$sdev^2
    pca.var.per <- base::data.frame(pc = base::seq_along(pca.var), val = base::round(pca.var / base::sum(pca.var) * 100, 1))
    pc2 <- if (base::ncol(pca$x) >= 2) pca$x[, 2] else base::rep(0, base::nrow(pca$x))
    pc2_var <- if (base::nrow(pca.var.per) >= 2) pca.var.per[2, 2] else 0
    pca.data <- base::data.frame(
      Sample = base::rownames(pca$x),
      X = pca$x[, 1],
      Y = pc2,
      Group = hcobject[["data"]][[base::paste0("set", i, "_anno")]][[hcobject[["global_settings"]][["voi"]]]]
    )

    num_colours <- base::length(base::unique(dplyr::pull(hcobject[["data"]][[base::paste0("set", i, "_anno")]], hcobject[["global_settings"]][["voi"]])))
    my_palette <- ggsci::pal_d3("category20")(20)

    if (base::is.null(cols)) {
      if (base::length(num_colours) > base::length(my_palette)) {
        my_colours <- grDevices::colorRampPalette(my_palette)(num_colours)
      } else {
        my_colours <- my_palette[base::seq_len(num_colours)]
      }
    } else {
      my_colours <- cols
    }


    p <- ggplot2::ggplot(pca.data, ggplot2::aes(x = X, y = Y, col = Group, label = Sample)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::ylab(base::paste0("PC 2", " (", pc2_var, "%)")) +
      ggplot2::xlab(base::paste0("PC 1", " (", pca.var.per[1, 2], "%)")) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(base::paste0(hcobject[["layers_names"]][i], " by topvar")) +
      ggplot2::scale_color_manual(values = my_colours)


    folder <- PCA_save_folder
    .hc_output_dir(folder)

    .hc_export_ggplot_file(
      file = .hc_output_file(
        base::paste0("PCA_topvar_", hcobject[["layers_names"]][i], ".pdf"),
        folder
      ),
      plot = p,
      width = 10,
      height = 7
    )
    plotlist[[i]] <- p
    pca_list[[i]] <- pca
  }
  if (base::length(hcobject[["layers"]]) == 1) {
    cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 1)
  } else {
    cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 2, align = "hv")
  }
  graphics::plot(cp)
  return(pca_list)
}


#' Plot PCA based on cluster expressions
#'
#' Subroutine to .hc_PCA_algo_compare_driver().
#' @noRd

plot_PCA_cluster <- function(gtc = NULL, algo = NULL, PCA_save_folder, cols = cols) {
  plotlist <- list()
  pca_list <- list()
  for (l in base::seq_along(hcobject[["layers"]])) {
    if (base::is.null(gtc)) {
      FC <- intra_sample_FC(l)
    } else {
      FC <- intra_sample_FC(l, gtc = gtc)
    }

    pca <- stats::prcomp(base::t(FC), scale = TRUE)
    pca.var <- pca$sdev^2
    pca.var.per <- base::data.frame(pc = base::seq_along(pca.var), val = base::round(pca.var / base::sum(pca.var) * 100, 1))
    anno <- hcobject[["data"]][[base::paste0("set", l, "_anno")]]
    if (base::ncol(pca$x) == 1) {
      pca.data <- base::data.frame(
        Sample = base::rownames(pca$x),
        X = pca$x[, 1],
        Y = 0,
        Group = dplyr::pull(anno, hcobject[["global_settings"]][["voi"]])
      )
    } else {
      pca.data <- base::data.frame(
        Sample = base::rownames(pca$x),
        X = pca$x[, 1],
        Y = pca$x[, 2],
        Group = dplyr::pull(anno, hcobject[["global_settings"]][["voi"]])
      )
    }
    num_colours <- base::length(base::unique(dplyr::pull(hcobject[["data"]][[base::paste0("set", l, "_anno")]], hcobject[["global_settings"]][["voi"]])))
    my_palette <- ggsci::pal_d3("category20")(20)

    if (base::is.null(cols)) {
      if (base::length(num_colours) > base::length(my_palette)) {
        my_colours <- grDevices::colorRampPalette(my_palette)(num_colours)
      } else {
        my_colours <- my_palette[base::seq_len(num_colours)]
      }
    } else {
      my_colours <- cols
    }


    p <- ggplot2::ggplot(pca.data, ggplot2::aes(x = X, y = Y, label = Sample, color = Group)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::scale_color_manual(values = my_colours) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(base::paste0(hcobject[["layers_names"]][l], " by module - ", algo)) +
      ggplot2::ylab(base::paste0("PC 2", " (", pca.var.per[2, 2], "%)")) +
      ggplot2::xlab(base::paste0("PC 1", " (", pca.var.per[1, 2], "%)"))
    plotlist[[l]] <- p
    pca_list[[l]] <- pca

    folder <- PCA_save_folder
    .hc_output_dir(folder)

    .hc_export_ggplot_file(
      file = .hc_output_file(
        base::paste0("PCA_module_", hcobject[["layers_names"]][l], "_", algo, ".pdf"),
        folder
      ),
      plot = p,
      width = 10,
      height = 7
    )
  }
  if (base::length(hcobject[["layers"]]) == 1) {
    cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 1)
  } else {
    cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 2, align = "hv")
  }
  graphics::plot(cp)
  return(pca_list)
}


#' Calculate Intra-Sample Fold Changes
#'
#' Mean cluster exression from mean sample expression.
#' Subroutine to .hc_PCA_algo_compare_driver().
#' @noRd

intra_sample_FC <- function(l, gtc = NULL) {
  if (base::is.null(gtc)) {
    gene_to_cluster <- .hc_gene_to_cluster_impl()
  } else {
    base::colnames(gtc) <- base::c("gene", "cluster")
    gene_to_cluster <- gtc
    gene_to_cluster$gene <- base::as.character(gene_to_cluster$gene)
    gene_to_cluster$cluster <- base::as.character(gene_to_cluster$cluster)
    base::colnames(gene_to_cluster) <- base::c("gene", "color")
  }


  counts <- hcobject[["data"]][[base::paste0("set", l, "_counts")]]

  mean_expression_per_sample <- base::apply(counts, 2, base::mean)

  mean_expression_per_cluster <- NULL

  for (c in base::unique(gene_to_cluster$color)) {
    if (!c == "white") {
      genes <- gene_to_cluster[gene_to_cluster$color == c, ] %>%
        dplyr::pull(., "gene")

      # drop = FALSE: a module with a single gene present in this layer would
      # otherwise collapse to a vector and break apply()
      filt_counts <- counts[base::rownames(counts) %in% genes, , drop = FALSE]
      if (base::nrow(filt_counts) == 0) {
        next
      }
      tmp <- base::colMeans(filt_counts) %>%
        base::as.data.frame() %>%
        base::t() %>%
        base::as.data.frame()
      base::rownames(tmp) <- c
      base::colnames(tmp) <- base::colnames(counts)
      mean_expression_per_cluster <- base::rbind(mean_expression_per_cluster, tmp)
    }
  }
  mean_expression_per_cluster <- base::rbind(mean_expression_per_cluster, mean_expression_per_sample)
  FC_from_mean_per_cluster <- base::apply(mean_expression_per_cluster, 2, function(x) {
    x / x[base::length(x)]
  }) %>% base::as.data.frame()
  base::colnames(FC_from_mean_per_cluster) <- base::colnames(mean_expression_per_cluster)
  base::rownames(FC_from_mean_per_cluster) <- base::rownames(mean_expression_per_cluster)
  FC_from_mean_per_cluster <- FC_from_mean_per_cluster[base::seq_len(base::nrow(FC_from_mean_per_cluster) - 1L), ]

  return(FC_from_mean_per_cluster)
}


#' Find New Control
#'
#' Detects sample subset with highest control content. Subroutine to .hc_cut_hclust_impl().
#' @noRd

find_new_ctrl <- function(anno, l) {
  ctrl_samples <- anno[base::grepl(hcobject[["global_settings"]][["control"]], anno[[base::paste0(hcobject[["global_settings"]][["voi"]], "_old")]], ignore.case = TRUE), ] %>% base::rownames()
  new_cons <- base::unique(dplyr::pull(anno, hcobject[["global_settings"]][["voi"]]))

  maxvec <- base::lapply(new_cons, function(x) {
    tmp <- anno[anno[[hcobject[["global_settings"]][["voi"]]]] == x, ]
    l <- base::intersect(ctrl_samples, base::rownames(tmp)) %>% base::length()
    return(l)
  }) %>% base::unlist()
  new_ctrl <- new_cons[base::which(maxvec == base::max(maxvec))[1]]
  message("New control: ", new_ctrl)
  anno[hcobject[["global_settings"]][["voi"]]][anno[hcobject[["global_settings"]][["voi"]]] == new_ctrl] <- base::paste0(hcobject[["layers_names"]][l], "_", hcobject[["global_settings"]][["control"]])
  return(anno)
}


#' Get Annotation Matrix
#'
#' Subroutine to .hc_col_anno_categorical_driver().
#' @noRd

get_anno_matrix <- function(variables) {
  all_mats <- list()
  for (v in base::seq_along(variables)) {
    if (base::is.na(variables[v])) {
      names <- dplyr::pull(hcobject[["data"]][[base::paste0("set", v, "_anno")]], hcobject[["global_settings"]][["voi"]]) %>% base::unique()
      layer_mat <- base::matrix(base::rep(0, base::length(names[names %in% base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])])), ncol = 1)
      base::rownames(layer_mat) <- names[names %in% base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])]
    } else {
      anno <- hcobject[["data"]][[base::paste0("set", v, "_anno")]]
      if (!variables[v] %in% colnames(anno)) {
        stop(variables[v], " is not a column name found in the annotation of dataset ", v, ". Please check the spelling.")
      }
      layer_mat <- NULL
      for (i in base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])[!base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]]) == "Gene"]) {
        if (i %in% dplyr::pull(anno, hcobject[["global_settings"]][["voi"]])) {
          tmp <- dplyr::filter(anno, base::as.vector(anno[hcobject[["global_settings"]][["voi"]]] == i)) %>%
            dplyr::pull(., variables[v])

          tmp <- base::table(tmp) %>%
            base::t() %>%
            base::as.data.frame()

          tmp$Var1 <- NULL
          base::colnames(tmp) <- c("var", "freq")
          mat <- base::matrix(tmp$freq, nrow = 1, byrow = TRUE)
          base::colnames(mat) <- tmp$var
          base::rownames(mat) <- i
          layer_mat <- dplyr::bind_rows(base::as.data.frame(layer_mat), base::as.data.frame(mat))
        }
      }
    }
    layer_mat[base::is.na(layer_mat)] <- 0
    layer_mat <- base::as.matrix(layer_mat)
    all_mats[[v]] <- layer_mat
  }
  return(all_mats)
}


#' Unify Matrices
#'
#' Subroutine to .hc_col_anno_categorical_driver().
#' @noRd

unify_mats <- function(mat_list) {
  merged_mat <- NULL
  new_rn <- NULL
  for (x in mat_list) {
    if (!base::is.null(merged_mat)) {
      if (base::all(base::colnames(merged_mat) %in% base::colnames(x))) {
        new_rn <- base::c(new_rn, base::rownames(x))
        x <- x[, base::colnames(merged_mat)]
        merged_mat <- base::rbind(merged_mat, x)
      } else {
        new_rn <- base::c(new_rn, base::rownames(x))
        merged_mat <- base::merge(merged_mat, x, by = "row.names", all = TRUE)
        merged_mat$Row.names <- NULL
      }
    } else {
      new_rn <- base::c(new_rn, base::rownames(x))
      merged_mat <- base::merge(merged_mat, x, by = "row.names", all = TRUE)
      merged_mat$Row.names <- NULL
    }
  }
  merged_mat[base::is.na(merged_mat)] <- 0
  merged_mat <- merged_mat[, base::colSums(merged_mat) > 0]
  base::rownames(merged_mat) <- new_rn
  return(merged_mat)
}


#' Adapted Version Of plot_cluster_heatmap()
#'
#' Subroutine to regroup_cluster_heatmap().
#' @noRd

replot_cluster_heatmap <- function(col_order = NULL,
                                   row_order = NULL,
                                   cluster_columns = FALSE,
                                   cluster_rows = TRUE,
                                   k = 0,
                                   return_HM = TRUE,
                                   cat_as_bp = NULL,
                                   file_name = "module_heatmap.pdf",
                                   GFCs,
                                   group,
                                   data) {
  # set user specific enrichments if they exist:
  if ("enriched_per_cluster" %in% base::names(hcobject[["satellite_outputs"]])) {
    if (!base::is.null(hcobject[["satellite_outputs"]][["enriched_per_cluster"]])) {
      user_enrichment_1 <- hcobject[["satellite_outputs"]][["enriched_per_cluster"]][["categories_per_cluster"]]
    }
  } else {
    user_enrichment_1 <- NULL
  }

  if ("enriched_per_cluster2" %in% base::names(hcobject[["satellite_outputs"]])) {
    if (!base::is.null(hcobject[["satellite_outputs"]][["enriched_per_cluster2"]])) {
      user_enrichment_2 <- hcobject[["satellite_outputs"]][["enriched_per_cluster2"]][["categories_per_cluster"]]
    }
  } else {
    user_enrichment_2 <- NULL
  }

  column_anno_categorical <- NULL

  column_anno_numerical <- NULL

  if (base::is.null(cat_as_bp)) {
    if (!base::is.null(column_anno_categorical)) {
      cat_as_bp <- base::rep(FALSE, base::length(column_anno_categorical))
    }
  }

  base::gc()

  # filter for included clusters (non-white)
  c_df <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], cluster_included == "yes")
  mat_heatmap <- NULL


  if (!base::is.null(row_order)) {
    for (c in row_order) {
      # get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)

      # GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(GFCs, Gene %in% genes)
      c_GFC_means <- base::apply(c_GFCs[, base::seq_len(base::ncol(c_GFCs) - 1L)], 2, base::mean)

      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)
    }
    base::rownames(mat_heatmap) <- row_order
  } else {
    for (c in base::unique(c_df$color)) {
      # get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)


      # GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(GFCs, Gene %in% genes)

      if (base::is.vector(c_GFCs)) {
        c_GFC_means <- cGFCs
      } else {
        c_GFC_means <- base::apply(c_GFCs[, base::seq_len(base::ncol(c_GFCs) - 1L)] %>%
          base::as.data.frame(), 2, base::mean)
      }


      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)
    }
    base::rownames(mat_heatmap) <- c_df$color
  }


  base::colnames(mat_heatmap) <- base::colnames(GFCs)[base::seq_len(base::ncol(GFCs) - 1L)]

  existing_heatmap_col_order <- .hc_heatmap_cache_info(
    hcobject[["integrated_output"]][["cluster_calc"]]
  )$col_order
  if (!isTRUE(cluster_columns)) {
    selected_col_order <- .hc_select_heatmap_col_order(
      available_cols = base::colnames(mat_heatmap),
      plot_order = col_order,
      main_order = existing_heatmap_col_order,
      fallback_order = base::colnames(mat_heatmap),
      context = "regrouped cluster heatmap"
    )
    if (base::length(selected_col_order) > 0) {
      mat_heatmap <- .hc_subset_matrix_cols_with_duplicates(
        mat_heatmap,
        selected_col_order
      ) %>% base::as.matrix()
    }
  }
  column_labels_display <- .hc_gfc_display_col_labels(hcobject, base::colnames(mat_heatmap))
  column_gap_spec <- .hc_heatmap_column_gap_spec(
    hcobject = hcobject,
    cols = base::colnames(mat_heatmap),
    cluster_columns = cluster_columns,
    gap_mm = 0.6
  )

  enrich_mat1 <- list()
  enrich_count1 <- list()
  enrich_mat2 <- list()
  enrich_count2 <- list()

  if (!base::is.null(row_order)) {
    if (!base::is.null(user_enrichment_1)) {
      for (x in row_order) {
        enrich_mat1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
          dplyr::pull(., count)
        enrich_count1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
          dplyr::pull(., hits) %>%
          dplyr::first(.)
      }
    }
    if (!base::is.null(user_enrichment_2)) {
      for (x in row_order) {
        enrich_mat2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
          dplyr::pull(., count)
        enrich_count2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
          dplyr::pull(., hits) %>%
          dplyr::first(.)
      }
    }
  } else {
    for (x in base::unique(user_enrichment_1$cluster)) {
      enrich_mat1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
        dplyr::pull(., count)
      enrich_count1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
        dplyr::pull(., hits) %>%
        dplyr::first(.)
    }
    for (x in base::unique(user_enrichment_2$cluster)) {
      enrich_mat2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
        dplyr::pull(., count)
      enrich_count2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
        dplyr::pull(., hits) %>%
        dplyr::first(.)
    }
  }

  if (!base::length(enrich_mat1) == 0) {
    enrich_mat1 <- base::matrix(base::unlist(enrich_mat1), nrow = base::length(enrich_mat1), byrow = TRUE)
    enrich_count1 <- base::unlist(enrich_count1)
  }
  if (base::all(enrich_mat1 == 0) == TRUE) {
    enrich_mat1 <- list()
  }

  if (!base::length(enrich_mat2) == 0) {
    enrich_mat2 <- base::matrix(base::unlist(enrich_mat2), nrow = base::length(enrich_mat2), byrow = TRUE)
    enrich_count2 <- base::unlist(enrich_count2)
  }
  if (base::all(enrich_mat2 == 0) == TRUE) {
    enrich_mat2 <- list()
  }

  if (!base::is.null(row_order)) {
    cluster_colors <- base::factor(row_order)
    base::names(cluster_colors) <- row_order
    c_df <- c_df[base::match(row_order, c_df$color), ]
  } else {
    cluster_colors <- base::factor(c_df$color)
    base::names(cluster_colors) <- c_df$color
    row_order <- base::unique(c_df$color)
  }


  if (base::length(enrich_mat1) == 0 & base::length(enrich_mat2) == 0) {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      modules = ComplexHeatmap::anno_simple(row_order,
        col = cluster_colors,
        simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")
      ),
      genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(2.5, "cm")),
      gene_nums = ComplexHeatmap::anno_text(c_df$gene_no, width = grid::unit(1.5, "cm"), gp = grid::gpar(fontsize = 10)),
      which = "row",
      width = grid::unit(4.5, "cm"),
      annotation_name_side = "top",
      gap = grid::unit(2, "mm"),
      annotation_name_rot = 0,
      annotation_name_gp = grid::gpar(fontsize = 8)
    )

    lgd_list <- list()
  } else if (base::length(enrich_mat1) > 0 & base::length(enrich_mat2) == 0) {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      modules = ComplexHeatmap::anno_simple(row_order,
        col = cluster_colors,
        simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")
      ),
      genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(2.5, "cm")),
      enriched_count = ComplexHeatmap::anno_text(base::paste0(enrich_count1, "/", c_df$gene_no), width = grid::unit(1.5, "cm")),
      enriched = ComplexHeatmap::anno_barplot(enrich_mat1,
        width = grid::unit(3, "cm"),
        gp = grid::gpar(
          fill = RColorBrewer::brewer.pal(n = 12, name = "Paired"),
          col = RColorBrewer::brewer.pal(n = 12, name = "Paired")
        )
      ),
      which = "row",
      width = grid::unit(9, "cm"),
      annotation_name_side = "top",
      gap = grid::unit(2, "mm"),
      annotation_name_rot = 0,
      annotation_name_gp = grid::gpar(fontsize = 8)
    )

    lgd_list <- list(
      ComplexHeatmap::Legend(
        labels = base::unique(user_enrichment_1$cell_type), title = "enriched",
        legend_gp = grid::gpar(col = RColorBrewer::brewer.pal(n = 12, name = "Paired")),
        type = "points", pch = 15
      )
    )
  } else if (base::length(enrich_mat1) == 0 & base::length(enrich_mat2) > 0) {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      modules = ComplexHeatmap::anno_simple(row_order,
        col = cluster_colors,
        simple_anno_size = grid::unit(0.5, "cm"), gp = grid::gpar(col = "black")
      ),
      genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(1.5, "cm")),
      enriched_count = ComplexHeatmap::anno_text(base::paste0(enrich_count2, "/", c_df$gene_no),
        width = grid::unit(1.5, "cm"),
        gp = grid::gpar(fontsize = 8)
      ),
      enriched = ComplexHeatmap::anno_barplot(enrich_mat2,
        width = grid::unit(5, "cm"),
        gp = grid::gpar(
          fill = ggsci::pal_d3(palette = "category20")(base::ncol(enrich_mat2)),
          col = ggsci::pal_d3(palette = "category20")(base::ncol(enrich_mat2))
        )
      ),
      which = "row",
      width = grid::unit(9, "cm"),
      annotation_name_side = "top",
      gap = grid::unit(2, "mm"),
      annotation_name_rot = 0,
      annotation_name_gp = grid::gpar(fontsize = 8)
    )

    lgd_list <- list(
      ComplexHeatmap::Legend(
        labels = base::unique(user_enrichment_2$cell_type), title = "enriched",
        legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)),
        type = "points", pch = 15
      )
    )
  } else {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      modules = ComplexHeatmap::anno_simple(row_order,
        col = cluster_colors,
        simple_anno_size = grid::unit(0.25, "cm"), gp = grid::gpar(col = "black")
      ),
      genes = ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(0.75, "cm")),
      enriched_count_1 = ComplexHeatmap::anno_text(base::paste0(enrich_count1, "/", c_df$gene_no),
        width = grid::unit(0.75, "cm"),
        gp = grid::gpar(fontsize = 8)
      ),
      enriched_1 = ComplexHeatmap::anno_barplot(enrich_mat1,
        width = grid::unit(2, "cm"),
        gp = grid::gpar(
          fill = RColorBrewer::brewer.pal(n = 12, name = "Paired")[base::seq_len(base::ncol(enrich_mat1))],
          col = RColorBrewer::brewer.pal(n = 12, name = "Paired")[base::seq_len(base::ncol(enrich_mat1))]
        ),
        baseline = 0
      ),
      enriched_count_2 = ComplexHeatmap::anno_text(base::paste0(enrich_count2, "/", c_df$gene_no),
        width = grid::unit(0.75, "cm"),
        gp = grid::gpar(fontsize = 8)
      ),
      enriched_2 = ComplexHeatmap::anno_barplot(enrich_mat2,
        width = grid::unit(2, "cm"),
        gp = grid::gpar(
          fill = ggsci::pal_d3(palette = "category20")(20)[base::seq_len(base::ncol(enrich_mat2))],
          col = ggsci::pal_d3(palette = "category20")(20)[base::seq_len(base::ncol(enrich_mat2))]
        ),
        baseline = 0
      ),
      which = "row",
      width = grid::unit(12, "cm"),
      annotation_name_side = "top",
      gap = grid::unit(2, "mm"),
      annotation_name_rot = 0,
      annotation_name_gp = grid::gpar(fontsize = 8)
    )

    lgd_list <- list(
      ComplexHeatmap::Legend(
        labels = base::unique(user_enrichment_1$cell_type), title = "enriched_1",
        legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)[base::seq_len(base::ncol(enrich_mat1))]),
        type = "points", pch = 15
      ),
      ComplexHeatmap::Legend(
        labels = base::unique(user_enrichment_2$cell_type), title = "enriched_2",
        legend_gp = grid::gpar(col = ggsci::pal_d3(palette = "category20")(20)[base::seq_len(base::ncol(enrich_mat2))]),
        type = "points", pch = 15
      )
    )
  }


  anno_list <- NULL


  if (!base::length(column_anno_categorical) == 0) {
    for (a in base::seq_along(column_anno_categorical)) {
      tmp_colour <- grDevices::colorRampPalette(ggsci::pal_d3("category20")(20))(base::ncol(column_anno_categorical[[a]]))
      if (cat_as_bp[a] == TRUE) {
        column_anno_categorical[[a]][base::is.na(column_anno_categorical[[a]])] <- 0
        if (base::is.null(anno_list)) {
          anno_list <- ComplexHeatmap::HeatmapAnnotation(
            col_anno = ComplexHeatmap::anno_barplot(column_anno_categorical[[a]] %>% base::as.matrix(),
              width = grid::unit(2, "cm"),
              gp = grid::gpar(
                fill = tmp_colour,
                col = tmp_colour
              )
            ),
            which = "column",
            height = grid::unit(1, "cm"),
            annotation_name_side = "right",
            gap = grid::unit(2, "mm"),
            annotation_name_rot = 0,
            annotation_name_gp = grid::gpar(fontsize = 8),
            annotation_label = base::names(column_anno_categorical)[a]
          )
        } else {
          anno_list <- ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(
            col_anno = ComplexHeatmap::anno_barplot(column_anno_categorical[[a]] %>% base::as.matrix(),
              width = grid::unit(2, "cm"),
              gp = grid::gpar(
                fill = tmp_colour,
                col = tmp_colour
              )
            ),
            which = "column",
            height = grid::unit(1, "cm"),
            annotation_name_side = "right",
            gap = grid::unit(2, "mm"),
            annotation_name_rot = 0,
            annotation_name_gp = grid::gpar(fontsize = 8),
            annotation_label = base::names(column_anno_categorical)[a]
          ), direction = "vertical")
        }
      } else {
        if (base::is.null(anno_list)) {
          anno_list <- ComplexHeatmap::HeatmapAnnotation(
            col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),
              width = grid::unit(2, "cm"),
              gp = grid::gpar(col = tmp_colour),
              add_points = TRUE,
              pt_gp = grid::gpar(col = tmp_colour), pch = 16
            ),
            which = "column",
            height = grid::unit(1, "cm"),
            annotation_name_side = "right",
            gap = grid::unit(2, "mm"),
            annotation_name_rot = 0,
            annotation_name_gp = grid::gpar(fontsize = 8),
            annotation_label = base::names(column_anno_categorical)[a]
          )
        } else {
          anno_list <- ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(
            col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),
              width = grid::unit(2, "cm"),
              gp = grid::gpar(col = tmp_colour),
              add_points = TRUE,
              pt_gp = grid::gpar(col = tmp_colour), pch = 16
            ),
            which = "column",
            height = grid::unit(1, "cm"),
            annotation_name_side = "right",
            gap = grid::unit(2, "mm"),
            annotation_name_rot = 0,
            annotation_name_gp = grid::gpar(fontsize = 8),
            annotation_label = base::names(column_anno_categorical)[a]
          ), direction = "vertical")
        }
      }


      lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(
        labels = colnames(column_anno_categorical[[a]] %>% as.matrix()), title = names(column_anno_categorical)[a],
        legend_gp = grid::gpar(col = tmp_colour),
        type = "points", pch = 15
      ))
    }
  }


  if (!base::length(column_anno_numerical) == 0) {
    for (a in base::seq_along(column_anno_numerical)) {
      tmp_col_anno_2 <- column_anno_numerical[[a]]
      tmp_col_anno_2 <- tmp_col_anno_2[base::colnames(mat_heatmap)]
      if (base::is.null(anno_list)) {
        anno_list <- ComplexHeatmap::HeatmapAnnotation(
          cont_anno = ComplexHeatmap::anno_boxplot(tmp_col_anno_2, height = grid::unit(1, "cm")),
          which = "column",
          annotation_name_side = "right",
          gap = grid::unit(2, "mm"),
          annotation_name_rot = 0,
          annotation_name_gp = grid::gpar(fontsize = 8),
          annotation_label = base::names(column_anno_numerical)[a], show_legend = FALSE
        )
      } else {
        anno_list <- ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::HeatmapAnnotation(
          cont_anno = ComplexHeatmap::anno_boxplot(tmp_col_anno_2, height = grid::unit(1, "cm")),
          which = "column",
          annotation_name_side = "right",
          gap = grid::unit(2, "mm"),
          annotation_name_rot = 0,
          annotation_name_gp = grid::gpar(fontsize = 8),
          annotation_label = base::names(column_anno_numerical)[a], show_legend = FALSE
        ), direction = "vertical")
      }
    }
  }


  all_conditions <- .hc_gfc_display_count_labels(hcobject, base::colnames(mat_heatmap))

  if (base::is.null(anno_list)) {
    anno_list <- ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions))
  } else {
    anno_list <- ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions)), direction = "vertical")
  }


  # }


  gfc_scale_limits <- .hc_as_numeric_safely(hcobject[["integrated_output"]][["cluster_calc"]][["gfc_scale_limits"]])
  if (base::length(gfc_scale_limits) == 1 && base::is.finite(gfc_scale_limits) && gfc_scale_limits > 0) {
    gfc_scale_limits <- c(-base::abs(gfc_scale_limits), base::abs(gfc_scale_limits))
  }
  if (base::length(gfc_scale_limits) != 2 || any(!base::is.finite(gfc_scale_limits))) {
    fallback_lim <- .hc_first_numeric_value(hcobject[["global_settings"]][["range_GFC"]])
    if (!base::is.finite(fallback_lim) || fallback_lim <= 0) {
      fallback_lim <- 2
    }
    gfc_scale_limits <- c(-base::abs(fallback_lim), base::abs(fallback_lim))
  } else {
    gfc_scale_limits <- base::sort(gfc_scale_limits)
    if (gfc_scale_limits[1] == gfc_scale_limits[2]) {
      lim_abs <- base::abs(gfc_scale_limits[1])
      if (!base::is.finite(lim_abs) || lim_abs <= 0) {
        lim_abs <- 2
      }
      gfc_scale_limits <- c(-lim_abs, lim_abs)
    }
  }
  gfc_scale_breaks <- pretty(gfc_scale_limits, n = 5)
  gfc_scale_breaks <- gfc_scale_breaks[
    gfc_scale_breaks >= gfc_scale_limits[1] - .Machine$double.eps^0.5 &
      gfc_scale_breaks <= gfc_scale_limits[2] + .Machine$double.eps^0.5
  ]
  if (!any(base::abs(gfc_scale_breaks) < .Machine$double.eps^0.5)) {
    gfc_scale_breaks <- base::sort(base::unique(base::c(gfc_scale_breaks, 0)))
  }
  if (base::length(gfc_scale_breaks) < 3) {
    gfc_scale_breaks <- base::seq(gfc_scale_limits[1], gfc_scale_limits[2], length.out = 5)
  }
  gfc_scale_labels <- base::formatC(gfc_scale_breaks, format = "fg", digits = 3)
  gfc_scale_labels <- base::trimws(gfc_scale_labels)
  gfc_label_width <- base::max(base::nchar(gfc_scale_labels), na.rm = TRUE)
  gfc_scale_labels <- base::format(gfc_scale_labels, width = gfc_label_width, justify = "right")
  gfc_palette <- grDevices::colorRampPalette(.hc_default_gfc_colors())(51)
  gfc_col_fun <- circlize::colorRamp2(
    seq(gfc_scale_limits[1], gfc_scale_limits[2], length.out = base::length(gfc_palette)),
    gfc_palette
  )

  hm_args <- list(
    matrix = mat_heatmap,
    right_annotation = ha,
    col = gfc_col_fun,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    column_names_rot = 90,
    column_labels = column_labels_display,
    column_names_centered = FALSE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    rect_gp = grid::gpar(col = "black"),
    heatmap_legend_param = list(
      title = "",
      at = gfc_scale_breaks,
      labels = gfc_scale_labels,
      title_gp = grid::gpar(fontsize = 7.6, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 6.6),
      legend_height = grid::unit(3, "cm")
    ),
    column_km = k
  )
  if (!base::is.null(column_gap_spec$column_split) &&
    !base::is.null(column_gap_spec$column_gap) &&
    (!base::is.numeric(k) || base::length(k) == 0 || base::all(k <= 0))) {
    hm_args$column_split <- column_gap_spec$column_split
    hm_args$column_gap <- column_gap_spec$column_gap
    hm_args$cluster_column_slices <- FALSE
    hm_args$column_title <- column_gap_spec$slice_titles
  }
  hm <- do.call(ComplexHeatmap::Heatmap, hm_args)

  if (base::is.null(anno_list)) {
    anno_list <- hm
  } else {
    anno_list <- ComplexHeatmap::add_heatmap(hm, anno_list, direction = c("vertical"))
  }

  export_draw_fun <- function() {
    ComplexHeatmap::draw(
      object = anno_list,
      annotation_legend_list = lgd_list,
      merge_legends = TRUE,
      padding = grid::unit(c(2, 2, 2, 30), "mm")
    )
  }
  .hc_export_single_page_plot(
    file = paste0(
      hcobject[["working_directory"]][["dir_output"]],
      hcobject[["global_settings"]][["save_folder"]],
      "/",
      file_name
    ),
    width = 50,
    height = 30,
    pointsize = 11,
    res = 300,
    pdf_dpi = 300,
    draw_fun = export_draw_fun
  )

  hm_w_lgd <- export_draw_fun()

  .hc_display_object(hm_w_lgd)
  if (return_HM) {
    return(hm_w_lgd)
  }
}
