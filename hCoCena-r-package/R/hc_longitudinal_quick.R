.hc_normalize_longitudinal_k <- function(k, arg_name) {
  if (is.null(k)) {
    return(NULL)
  }
  k <- base::as.integer(k)
  k <- k[base::is.finite(k) & !base::is.na(k)]
  k <- base::sort(base::unique(k))
  if (base::length(k) == 0) {
    stop("`", arg_name, "` must contain at least one finite integer.")
  }
  if (base::any(k < 2)) {
    stop("`", arg_name, "` must contain integers >= 2.")
  }
  k
}

.hc_normalize_time_levels <- function(time_levels) {
  if (missing(time_levels)) {
    stop("`time_levels` is required and must be an explicit ordered timepoint vector.")
  }
  if (!base::is.atomic(time_levels) || base::length(time_levels) == 0) {
    stop("`time_levels` must be a non-empty vector.")
  }
  time_levels <- base::as.character(time_levels)
  time_levels <- time_levels[!base::is.na(time_levels) & base::nzchar(time_levels)]
  if (base::length(time_levels) == 0) {
    stop("`time_levels` must contain at least one non-empty value.")
  }
  base::unique(time_levels)
}

.hc_longitudinal_quick_preset_defaults <- function(preset = c("standard", "consensus")) {
  preset <- base::match.arg(preset)
  if (identical(preset, "standard")) {
    return(list(
      module_k = 2:4,
      clustering_nstart = 100,
      cap_runs = 50,
      cap_nstart = 30,
      meta_k = 2:8,
      meta_consensus_runs = 250,
      meta_consensus_sample_fraction = 0.8,
      meta_consensus_feature_fraction = 0.85,
      include_consensus = FALSE,
      cap_donor_order = "none",
      meta_save_tables = TRUE
    ))
  }

  list(
    module_k = 2:4,
    clustering_nstart = 100,
    cap_runs = 50,
    cap_nstart = 30,
    meta_k = 3:8,
    meta_consensus_runs = 500,
    meta_consensus_sample_fraction = 0.9,
    meta_consensus_feature_fraction = 0.9,
    include_consensus = TRUE,
    cap_donor_order = "module_cluster_heatmap",
    meta_save_tables = FALSE
  )
}

#' Longitudinal step 1: module/donor clustering
#'
#' Computes longitudinal module means, exact `kml` module/donor clustering, and
#' CAP, then returns the two wave plots and two heatmaps for this step.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param donor_col Annotation column containing donor IDs.
#' @param time_col Annotation column containing timepoint labels.
#' @param layer Optional layer index, layer id, or layer name. `NULL` uses the
#'   first available layer.
#' @param time_levels Required explicit timepoint order vector.
#' @param k Candidate donor-trajectory cluster numbers per module. This is the
#'   main parameter controlling how many donor clusters each module may form.
#' @param rerolls Number of repeated `kml` redrawings per module. Increase this
#'   for more stable results.
#' @param impute Logical. If `TRUE`, impute missing donor-time values before
#'   `kml` clustering.
#' @param impute_method Missing-value method passed to `mice`. Use `"rfcont"`
#'   to match the previous workflow most closely.
#' @param ntree Number of trees for `rfcont` imputation.
#' @param min_cluster_fraction Minimum allowed fraction of donors in the
#'   smallest cluster when scoring candidate `k`.
#' @param score_method Rule used to choose the final per-module clustering from
#'   the `kml` score tables.
#' @param seed Random seed.
#' @param means_slot Satellite slot for module means.
#' @param output_slot Satellite slot for endotype outputs.
#'
#' @return A list with updated `hc`, step plots, and diagnostics.
#' @export
hc_longitudinal_step1_module_donor <- function(hc,
                                               donor_col = "Subject",
                                               time_col = "Time_token",
                                               layer = NULL,
                                               time_levels,
                                               k = 2:6,
                                               rerolls = 1000,
                                               impute = TRUE,
                                               impute_method = "rfcont",
                                               ntree = 10,
                                               min_cluster_fraction = 0.1,
                                               score_method = "median",
                                               seed = 42,
                                               means_slot = "longitudinal_module_means",
                                               output_slot = "longitudinal_endotypes") {
  time_levels <- .hc_normalize_time_levels(time_levels)
  k <- .hc_normalize_longitudinal_k(k, "k")
  if (!base::is.numeric(rerolls) || base::length(rerolls) != 1 || !base::is.finite(rerolls) || rerolls < 1) {
    stop("`rerolls` must be a single integer >= 1.")
  }
  rerolls <- as.integer(rerolls)
  if (!base::is.logical(impute) || base::length(impute) != 1 || base::is.na(impute)) {
    stop("`impute` must be TRUE or FALSE.")
  }
  if (!base::is.character(impute_method) || base::length(impute_method) != 1 || !base::nzchar(impute_method)) {
    stop("`impute_method` must be a non-empty character scalar.")
  }
  if (!base::is.numeric(ntree) || base::length(ntree) != 1 || !base::is.finite(ntree) || ntree < 1) {
    stop("`ntree` must be a single integer >= 1.")
  }
  if (!base::is.numeric(min_cluster_fraction) || base::length(min_cluster_fraction) != 1 ||
      !base::is.finite(min_cluster_fraction) || min_cluster_fraction < 0 || min_cluster_fraction > 1) {
    stop("`min_cluster_fraction` must be between 0 and 1.")
  }
  if (!score_method %in% c("max", "mean", "median", "z_max", "z_mean", "z_median")) {
    stop("`score_method` must be one of `max`, `mean`, `median`, `z_max`, `z_mean`, `z_median`.")
  }

  lmm_args <- list(
    hc = hc,
    donor_col = donor_col,
    time_col = time_col,
    layer = layer,
    group_col = NULL,
    use_module_labels = TRUE,
    time_levels = time_levels,
    impute_missing = "none",
    slot_name = means_slot
  )
  if ("value_label" %in% base::names(base::formals(hc_longitudinal_module_means))) {
    lmm_args$value_label <- NULL
  }
  hc <- do.call(hc_longitudinal_module_means, lmm_args)

  hc <- .hc_run_legacy_step1_exact(
    hc = hc,
    means_slot = means_slot,
    output_slot = output_slot,
    k = k,
    rerolls = rerolls,
    impute = impute,
    impute_method = impute_method,
    ntree = ntree,
    min_cluster_fraction = min_cluster_fraction,
    score_method = score_method,
    seed = seed
  )

  module_means <- hc_plot_longitudinal_module_means(
    hc,
    slot_name = means_slot,
    save_pdf = TRUE,
    file_prefix = "Longitudinal_ModuleMeans",
    square_panels = TRUE,
    save_width = 10,
    save_height = 10
  )

  module_clusters <- hc_plot_longitudinal_module_clusters(
    hc,
    slot_name = output_slot,
    save_pdf = TRUE,
    file_prefix = "Longitudinal_ModuleClusters",
    show_heatmap_numbers = TRUE,
    square_panels = TRUE,
    save_waves_width = 10,
    save_waves_height = 10
  )

  donor_order <- NULL
  heatmap_obj <- module_clusters$module_cluster_heatmap
  if (!base::is.null(heatmap_obj) &&
      !base::is.null(heatmap_obj$data) &&
      "donor" %in% base::colnames(heatmap_obj$data)) {
    donor_vec <- heatmap_obj$data$donor
    donor_levels <- base::levels(donor_vec)
    donor_order <- if (!base::is.null(donor_levels) && base::length(donor_levels) > 0) {
      donor_levels
    } else {
      base::unique(base::as.character(donor_vec))
    }
  }

  cap_args <- list(
    hc = hc,
    slot_name = output_slot,
    save_pdf = TRUE,
    file_prefix = "Longitudinal_CAP",
    show_values = FALSE
  )
  if (!base::is.null(donor_order)) {
    cap_args$donor_order <- donor_order
  }
  cap <- do.call(hc_plot_longitudinal_cap, cap_args)

  list(
    hc = hc,
    plots = list(
      module_means_waves = module_means$module_means,
      module_cluster_waves = module_clusters$module_cluster_waves,
      module_cluster_heatmap = module_clusters$module_cluster_heatmap,
      cap_heatmap = cap$cap_heatmap
    ),
    diagnostics = list(
      k = k,
      rerolls = as.integer(rerolls),
      impute = isTRUE(impute),
      impute_method = impute_method,
      ntree = as.integer(ntree),
      min_cluster_fraction = as.numeric(min_cluster_fraction),
      score_method = score_method,
      module_cluster_score = hc@satellite[[output_slot]]$module_cluster_score,
      module_cluster_best_k = hc@satellite[[output_slot]]$module_cluster_best_k
    )
  )
}

#' Longitudinal step 2: meta-clustering
#'
#' Runs the exact CAP -> PCA -> graph -> Leiden donor meta-clustering and
#' returns the PCA/UMAP embeddings.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing step 1 outputs.
#' @param dimensions Number of PCA dimensions used to build the donor graph.
#' @param graph_method Graph type used before Leiden clustering.
#' @param knn_method Nearest-neighbor backend used for graph construction.
#' @param graph_k Number of neighbors used for the donor graph.
#' @param resolution Leiden resolution parameter.
#' @param leiden_method Leiden partition method.
#' @param cluster_prefix Prefix for meta-cluster labels.
#' @param compute_umap Logical. If `TRUE`, compute an additional UMAP from the
#'   PCA subspace used for graph clustering.
#' @param umap_neighbors UMAP neighborhood size.
#' @param umap_min_dist UMAP minimum distance.
#' @param seed Random seed.
#'
#' @return A list with updated `hc`, step plots, and diagnostics.
#' @export
hc_longitudinal_step2_meta_clustering <- function(hc,
                                                  slot_name = "longitudinal_endotypes",
                                                  dimensions = 4,
                                                  graph_method = c("knn", "snn"),
                                                  knn_method = c("annoy", "rann"),
                                                  graph_k = 7,
                                                  resolution = 0.4,
                                                  leiden_method = "RBConfigurationVertexPartition",
                                                  cluster_prefix = "MC",
                                                  compute_umap = TRUE,
                                                  umap_neighbors = 15,
                                                  umap_min_dist = 0.3,
                                                  seed = 42) {
  graph_method <- base::match.arg(graph_method)
  knn_method <- base::match.arg(knn_method)
  if (!base::is.numeric(dimensions) || base::length(dimensions) != 1 || !base::is.finite(dimensions) || dimensions < 1) {
    stop("`dimensions` must be a single integer >= 1.")
  }
  if (!base::is.numeric(graph_k) || base::length(graph_k) != 1 || !base::is.finite(graph_k) || graph_k < 1) {
    stop("`graph_k` must be a single integer >= 1.")
  }
  if (!base::is.numeric(resolution) || base::length(resolution) != 1 || !base::is.finite(resolution) || resolution <= 0) {
    stop("`resolution` must be a single positive number.")
  }

  hc <- .hc_run_legacy_step2_exact(
    hc = hc,
    slot_name = slot_name,
    dimensions = as.integer(dimensions),
    graph_method = graph_method,
    knn_method = knn_method,
    graph_k = as.integer(graph_k),
    resolution = as.numeric(resolution),
    leiden_method = leiden_method,
    cluster_prefix = cluster_prefix,
    seed = seed,
    compute_umap = compute_umap,
    umap_neighbors = umap_neighbors,
    umap_min_dist = umap_min_dist
  )

  meta_embeddings <- hc_plot_longitudinal_meta_embeddings(
    hc,
    slot_name = slot_name,
    save_pdf = TRUE,
    file_prefix = "Longitudinal_Meta",
    show_endotype_crosstab = FALSE,
    show_cluster_labels = TRUE,
    save_tables = TRUE,
    table_format = "xlsx",
    table_detail = "simple"
  )
  meta_obj <- hc@satellite[[slot_name]]

  list(
    hc = hc,
    plots = list(
      pca = meta_embeddings$pca,
      umap = meta_embeddings$umap
    ),
    diagnostics = list(
      k_criterion = NULL,
      meta_consensus = NULL,
      cross_tab = meta_embeddings$cross_tab,
      tables = meta_embeddings$tables,
      graph_method = graph_method,
      knn_method = knn_method,
      graph_k = as.integer(graph_k),
      resolution = as.numeric(resolution),
      dimensions = as.integer(dimensions),
      leiden_method = leiden_method,
      meta_cluster = meta_obj$meta_cluster,
      meta_method_comparison = meta_obj$meta_method_comparison,
      meta_score_table = meta_obj$meta_score_table
    )
  )
}

#' Longitudinal step 3: module trajectories by meta-cluster
#'
#' Plots module trajectories after meta-clustering using donor and meta-cluster
#' mean waves.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing meta-clustering outputs.
#' @param facet_ncol Number of facet columns.
#' @param free_y Logical. If `TRUE`, use free y-scales per module.
#' @param square_panels Logical. If `TRUE`, draw square panels.
#' @param value_mode One of `"scaled_mean_vst"` or `"expression"`.
#' @param value_range Optional numeric vector of length 2 for the y-axis range
#'   when `value_mode = "scaled_mean_vst"`. Default is `c(-2, 2)`. Use `NULL`
#'   for automatic scaling.
#'
#' @return A list with unchanged `hc` and the step 3 plot.
#' @export
hc_longitudinal_step3_meta_module_trajectories <- function(hc,
                                                           slot_name = "longitudinal_endotypes",
                                                           facet_ncol = 4,
                                                           free_y = TRUE,
                                                           square_panels = TRUE,
                                                           value_mode = c("scaled_mean_vst", "expression"),
                                                           value_range = c(-2, 2)) {
  value_mode <- base::match.arg(value_mode)
  meta_module_waves <- hc_plot_longitudinal_meta_module_waves(
    hc,
    slot_name = slot_name,
    save_pdf = TRUE,
    file_prefix = "Longitudinal_Meta_ModuleWaves",
    facet_ncol = facet_ncol,
    free_y = free_y,
    square_panels = square_panels,
    value_mode = value_mode,
    value_range = value_range,
    save_width = 10,
    save_height = 10
  )

  list(
    hc = hc,
    plots = list(
      meta_module_waves = meta_module_waves$meta_module_waves
    )
  )
}

#' Print longitudinal endotype outputs
#'
#' Convenience printer for outputs from the stepwise longitudinal functions.
#'
#' @param x Output list from one of the longitudinal step functions, or a list
#'   containing nested `step1_module_donor`, `step2_meta_clustering`, and
#'   `step3_meta_module_trajectories` plot blocks.
#' @param show_tables Logical. If `TRUE`, print optional meta-clustering tables
#'   when available.
#'
#' @return Invisibly returns `x`.
#' @export
hc_print_longitudinal_endotypes <- function(x, show_tables = TRUE) {
  if (base::is.null(x$plots)) {
    return(base::invisible(x))
  }
  p <- x$plots

  step1 <- p$step1_module_donor
  step2 <- p$step2_meta_clustering
  step3 <- p$step3_meta_module_trajectories

  if (base::is.null(step1) &&
      base::any(base::c("module_means_waves", "module_cluster_waves", "module_cluster_heatmap", "cap_heatmap") %in% base::names(p))) {
    step1 <- p
  }
  if (base::is.null(step2) &&
      base::any(base::c("pca", "umap") %in% base::names(p))) {
    step2 <- p
  }
  if (base::is.null(step3) &&
      "meta_module_waves" %in% base::names(p)) {
    step3 <- p
  }

  if (!base::is.null(step1$module_means_waves)) print(step1$module_means_waves)
  if (!base::is.null(step1$module_cluster_waves)) print(step1$module_cluster_waves)
  if (!base::is.null(step1$module_cluster_heatmap)) print(step1$module_cluster_heatmap)
  if (!base::is.null(step1$cap_heatmap)) print(step1$cap_heatmap)
  if (!base::is.null(step2$pca)) print(step2$pca)
  if (!base::is.null(step2$umap)) print(step2$umap)
  if (!base::is.null(step3$meta_module_waves)) print(step3$meta_module_waves)

  meta_embeddings <- p$meta_embeddings
  if (base::is.null(meta_embeddings) && !base::is.null(x$diagnostics)) {
    meta_embeddings <- list(
      cross_tab = x$diagnostics$cross_tab,
      tables = x$diagnostics$tables
    )
  }
  k_criterion <- if (!base::is.null(p$k_criterion)) p$k_criterion else x$diagnostics$k_criterion
  meta_consensus <- if (!base::is.null(p$meta_consensus)) p$meta_consensus else x$diagnostics$meta_consensus

  if (!base::is.null(meta_embeddings$cross_tab)) print(meta_embeddings$cross_tab)
  if (isTRUE(show_tables) && !base::is.null(k_criterion$global_k)) print(k_criterion$global_k)
  if (isTRUE(show_tables) && !base::is.null(k_criterion$module_k)) print(k_criterion$module_k)
  if (isTRUE(show_tables) && !base::is.null(meta_consensus$consensus_k)) print(meta_consensus$consensus_k)
  if (isTRUE(show_tables) && !base::is.null(meta_consensus$stability)) print(meta_consensus$stability)
  if (isTRUE(show_tables) && !base::is.null(meta_consensus$consensus_heatmap)) print(meta_consensus$consensus_heatmap)
  if (isTRUE(show_tables) && !base::is.null(meta_embeddings$tables$Method_Clusters)) {
    print(meta_embeddings$tables$Method_Clusters)
  }
  if (isTRUE(show_tables) && !base::is.null(meta_embeddings$tables$Consensus_Decision)) {
    print(meta_embeddings$tables$Consensus_Decision)
  }

  base::invisible(x)
}

#' Add meta-cluster x time grouping column for heatmap regrouping
#'
#' Creates an annotation column that combines donor meta-cluster assignment and
#' timepoint labels (e.g. `MC1__24h`). This column can be used as `grouping_v`
#' in `hc_run_expression_analysis_2()` to build a cluster heatmap by
#' meta-cluster x timepoint.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param donor_col Annotation column containing donor IDs.
#' @param time_col Annotation column containing timepoint labels.
#' @param grouping_col Name of the new annotation column to create.
#' @param slot_name Satellite slot holding longitudinal outputs. Must contain
#'   `meta_cluster` table.
#' @param layer Layer index, layer id/name, or `"all"` (default) to add the
#'   grouping column to all layers.
#' @param time_levels Optional explicit timepoint order vector.
#' @param meta_cluster_levels Optional explicit meta-cluster order vector.
#' @param separator String between meta-cluster and timepoint labels.
#' @param append_layer_suffix Logical. If `TRUE`, append layer id to each group
#'   label (recommended when multiple layers share group names).
#'
#' @return Updated `HCoCenaExperiment`. The computed order is stored in
#'   `hc@satellite[[slot_name]]$meta_time_grouping`.
#' @export
hc_add_meta_time_grouping <- function(hc,
                                      donor_col = "Subject",
                                      time_col = "Time_token",
                                      grouping_col = "meta_cluster_time",
                                      slot_name = "longitudinal_endotypes",
                                      layer = "all",
                                      time_levels = NULL,
                                      meta_cluster_levels = NULL,
                                      separator = "__",
                                      append_layer_suffix = FALSE) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  if (!base::is.character(grouping_col) || base::length(grouping_col) != 1 || !base::nzchar(grouping_col)) {
    stop("`grouping_col` must be a non-empty string.")
  }
  if (!base::is.character(slot_name) || base::length(slot_name) != 1 || !base::nzchar(slot_name)) {
    stop("`slot_name` must be a non-empty string.")
  }
  if (!base::is.character(separator) || base::length(separator) != 1) {
    stop("`separator` must be a single string.")
  }

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (base::is.null(obj) || !base::is.list(obj) || base::is.null(obj$meta_cluster)) {
    stop("Slot `", slot_name, "` must contain `meta_cluster`. Run `hc_longitudinal_step2_meta_clustering()` first.")
  }

  meta_tbl <- base::as.data.frame(obj$meta_cluster, stringsAsFactors = FALSE)
  if (!base::all(c("donor", "meta_cluster") %in% base::colnames(meta_tbl))) {
    stop("`meta_cluster` table must contain columns `donor` and `meta_cluster`.")
  }
  meta_tbl$donor <- base::trimws(base::as.character(meta_tbl$donor))
  meta_tbl$meta_cluster <- base::trimws(base::as.character(meta_tbl$meta_cluster))
  keep_meta <- !base::is.na(meta_tbl$donor) & base::nzchar(meta_tbl$donor) &
    !base::is.na(meta_tbl$meta_cluster) & base::nzchar(meta_tbl$meta_cluster)
  meta_tbl <- meta_tbl[keep_meta, , drop = FALSE]
  meta_tbl <- meta_tbl[!base::duplicated(meta_tbl$donor), , drop = FALSE]
  if (base::nrow(meta_tbl) == 0) {
    stop("`meta_cluster` table does not contain valid donor assignments.")
  }

  if (base::is.null(meta_cluster_levels)) {
    meta_cluster_levels <- base::sort(base::unique(meta_tbl$meta_cluster))
  } else {
    meta_cluster_levels <- base::trimws(base::as.character(meta_cluster_levels))
    meta_cluster_levels <- meta_cluster_levels[!base::is.na(meta_cluster_levels) & base::nzchar(meta_cluster_levels)]
    obs_meta <- base::sort(base::unique(meta_tbl$meta_cluster))
    miss_meta <- obs_meta[!obs_meta %in% meta_cluster_levels]
    if (base::length(miss_meta) > 0) {
      warning("`meta_cluster_levels` missed observed clusters: ", base::paste(miss_meta, collapse = ", "), ". Appending at the end.")
      meta_cluster_levels <- base::c(meta_cluster_levels, miss_meta)
    }
  }
  meta_cluster_levels <- base::unique(meta_cluster_levels)

  exps <- MultiAssayExperiment::experiments(hc@mae)
  exp_names <- base::names(exps)
  if (base::length(exp_names) == 0) {
    stop("No layers found in `hc@mae`.")
  }

  target_layers <- NULL
  if (base::is.character(layer) && base::length(layer) == 1 && identical(layer, "all")) {
    target_layers <- exp_names
  } else if (base::is.null(layer)) {
    target_layers <- exp_names
  } else {
    target_layers <- .hc_resolve_layer_id(hc = hc, layer = layer)
  }

  col_order_by_layer <- list()
  time_levels_by_layer <- list()

  for (lid in target_layers) {
    se <- exps[[lid]]
    anno <- base::as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
    if (!base::all(c(donor_col, time_col) %in% base::colnames(anno))) {
      stop("Layer `", lid, "` is missing required annotation columns: `", donor_col, "` and/or `", time_col, "`.")
    }

    donor_vec <- base::trimws(base::as.character(anno[[donor_col]]))
    time_vec <- base::trimws(base::as.character(anno[[time_col]]))
    donor_vec[!base::nzchar(donor_vec)] <- NA_character_
    time_vec[!base::nzchar(time_vec)] <- NA_character_

    if (base::is.null(time_levels)) {
      raw_time <- SummarizedExperiment::colData(se)[[time_col]]
      if (base::is.factor(raw_time)) {
        tl <- base::trimws(base::levels(raw_time))
        tl <- tl[!base::is.na(tl) & base::nzchar(tl)]
        obs_t <- base::unique(time_vec[!base::is.na(time_vec) & base::nzchar(time_vec)])
        time_levels_use <- tl[tl %in% obs_t]
        if (base::length(time_levels_use) == 0) {
          time_levels_use <- base::sort(obs_t)
        }
      } else {
        time_levels_use <- base::sort(base::unique(time_vec[!base::is.na(time_vec) & base::nzchar(time_vec)]))
      }
    } else {
      time_levels_use <- base::trimws(base::as.character(time_levels))
      time_levels_use <- time_levels_use[!base::is.na(time_levels_use) & base::nzchar(time_levels_use)]
      obs_t <- base::unique(time_vec[!base::is.na(time_vec) & base::nzchar(time_vec)])
      miss_t <- obs_t[!obs_t %in% time_levels_use]
      if (base::length(miss_t) > 0) {
        warning("Layer `", lid, "` has time values missing in `time_levels`: ",
                base::paste(miss_t, collapse = ", "), ". Appending at the end.")
        time_levels_use <- base::c(time_levels_use, miss_t)
      }
    }
    time_levels_use <- base::unique(time_levels_use)
    if (base::length(time_levels_use) == 0) {
      stop("No valid time levels found for layer `", lid, "`.")
    }

    meta_vec <- meta_tbl$meta_cluster[base::match(donor_vec, meta_tbl$donor)]
    time_idx <- base::match(time_vec, time_levels_use)
    time_use <- base::ifelse(base::is.na(time_idx), NA_character_, time_levels_use[time_idx])

    valid_pair <- !base::is.na(meta_vec) & !base::is.na(time_use) & base::nzchar(time_use)
    combined <- base::rep(NA_character_, base::length(time_use))
    combined[valid_pair] <- base::paste0(meta_vec[valid_pair], separator, time_use[valid_pair])

    combo_levels <- base::as.vector(base::outer(
      meta_cluster_levels,
      time_levels_use,
      FUN = function(mc, tm) base::paste0(mc, separator, tm)
    ))

    if (isTRUE(append_layer_suffix)) {
      combined[!base::is.na(combined)] <- base::paste0(combined[!base::is.na(combined)], "_", lid)
      combo_levels <- base::paste0(combo_levels, "_", lid)
    }

    anno[[grouping_col]] <- base::factor(combined, levels = combo_levels, ordered = TRUE)
    n_assigned <- base::sum(!base::is.na(anno[[grouping_col]]))
    if (n_assigned == 0) {
      donor_examples <- base::unique(donor_vec[!base::is.na(donor_vec)])
      donor_examples <- utils::head(donor_examples, 8)
      stop(
        "No sample in layer `", lid, "` could be assigned to `", grouping_col, "`.\n",
        "Likely mismatch between annotation donors/times and longitudinal meta-cluster/time levels.\n",
        "Check donor values in `", donor_col, "` and time values in `", time_col, "`.\n",
        "Example donors in this layer: ", base::paste(donor_examples, collapse = ", ")
      )
    }
    n_unassigned <- base::length(anno[[grouping_col]]) - n_assigned
    if (n_unassigned > 0) {
      warning(
        "Layer `", lid, "`: ", n_unassigned, " sample(s) could not be assigned to `",
        grouping_col, "` and remain NA."
      )
    }

    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(anno, check.names = FALSE)
    exps[[lid]] <- se

    present_levels <- combo_levels[combo_levels %in% base::unique(base::as.character(stats::na.omit(anno[[grouping_col]])))]
    col_order_by_layer[[lid]] <- present_levels
    time_levels_by_layer[[lid]] <- time_levels_use
  }

  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = exps)

  if (base::is.null(sat[[slot_name]]) || !base::is.list(sat[[slot_name]])) {
    sat[[slot_name]] <- list()
  }
  sat[[slot_name]][["meta_time_grouping"]] <- list(
    grouping_col = grouping_col,
    donor_col = donor_col,
    time_col = time_col,
    separator = separator,
    append_layer_suffix = isTRUE(append_layer_suffix),
    layer_ids = target_layers,
    meta_cluster_levels = meta_cluster_levels,
    time_levels_by_layer = time_levels_by_layer,
    col_order_by_layer = col_order_by_layer,
    col_order = base::unlist(col_order_by_layer, use.names = FALSE)
  )
  hc@satellite <- S4Vectors::SimpleList(sat)

  methods::validObject(hc)
  hc
}

#' Get column order for meta-cluster x timepoint heatmaps
#'
#' Returns the column-order vector stored by `hc_add_meta_time_grouping()`.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param slot_name Satellite slot containing `meta_time_grouping`.
#' @param layer Layer index, layer id/name, or `NULL` for combined order.
#'
#' @return Character vector with heatmap column order.
#' @export
hc_get_meta_time_col_order <- function(hc,
                                       slot_name = "longitudinal_endotypes",
                                       layer = NULL) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (base::is.null(obj) || !base::is.list(obj) || base::is.null(obj$meta_time_grouping)) {
    stop("No `meta_time_grouping` found in slot `", slot_name, "`. Run `hc_add_meta_time_grouping()` first.")
  }
  info <- obj$meta_time_grouping

  if (base::is.null(layer)) {
    out <- base::as.character(info$col_order)
    out <- out[!base::is.na(out) & base::nzchar(out)]
    return(base::unique(out))
  }

  lid <- .hc_resolve_layer_id(hc = hc, layer = layer)
  if (base::is.null(info$col_order_by_layer[[lid]])) {
    stop("No stored col-order for layer `", lid, "`.")
  }
  out <- base::as.character(info$col_order_by_layer[[lid]])
  out <- out[!base::is.na(out) & base::nzchar(out)]
  base::unique(out)
}
