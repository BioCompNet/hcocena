.hc_require_namespace <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package `", pkg, "` is required")
    if (!is.null(reason) && nzchar(reason)) {
      msg <- paste0(msg, " for ", reason)
    }
    stop(msg, ". Please install it first.", call. = FALSE)
  }
}

.hc_legacy_impute_time_data <- function(time_data,
                                        donor_col = "donor",
                                        method = "rfcont",
                                        ntree = 10,
                                        m = 50,
                                        maxit = 50,
                                        seed = 42) {
  .hc_require_namespace("mice", "legacy longitudinal imputation")
  if (identical(method, "rfcont")) {
    .hc_require_namespace("CALIBERrfimpute", "legacy `rfcont` imputation")
    if (!"package:CALIBERrfimpute" %in% base::search()) {
      suppressPackageStartupMessages(
        base::library("CALIBERrfimpute", character.only = TRUE)
      )
    }
    CALIBERrfimpute::setRFoptions(ntree_cont = ntree)
  }

  donor_ids <- time_data[[donor_col]]
  data_use <- time_data[, base::setdiff(base::colnames(time_data), donor_col), drop = FALSE]
  base::colnames(data_use) <- base::paste0("t", base::colnames(data_use))

  set.seed(seed)
  mids <- mice::mice(
    data = as.data.frame(data_use, stringsAsFactors = FALSE),
    m = m,
    method = method,
    seed = seed,
    maxit = maxit,
    printFlag = FALSE
  )

  data_comp_list <- lapply(seq_len(m), function(i) {
    comp <- mice::complete(mids, i)
    comp$..row_id <- seq_len(base::nrow(comp))
    comp
  })
  data_comp <- do.call(base::rbind, data_comp_list)

  agg <- stats::aggregate(
    data_comp[, base::setdiff(base::colnames(data_comp), "..row_id"), drop = FALSE],
    by = list(..row_id = data_comp$..row_id),
    FUN = stats::median,
    na.rm = TRUE
  )
  agg <- agg[base::order(agg$..row_id), , drop = FALSE]
  agg$..row_id <- NULL
  base::colnames(agg) <- gsub("^t", "", base::colnames(agg))
  agg[[donor_col]] <- donor_ids
  agg <- agg[, c(donor_col, base::setdiff(base::colnames(agg), donor_col)), drop = FALSE]
  agg
}

.hc_legacy_kml_cluster_one <- function(time_data,
                                       k = 2:6,
                                       rerolls = 1000,
                                       maxNA = 1) {
  .hc_require_namespace("kml", "legacy longitudinal module clustering")
  longData_df <- kml::cld(
    time_data,
    timeInData = 2:base::ncol(time_data),
    maxNA = maxNA
  )
  longData_df["initializationMethod"] <- "kmeans++"
  par_algo <- kml::parALGO(
    startingCond = "kmeans++",
    imputationMethod = "copyMean",
    scale = TRUE,
    saveFreq = Inf
  )
  kml::kml(
    longData_df,
    nbRedrawing = rerolls,
    nbClusters = k,
    parAlgo = par_algo
  )
  longData_df
}

.hc_legacy_kml_cluster_scoring <- function(kml_object,
                                           min_cluster_fraction = 0.1) {
  kml_criterion_df_list <- lapply(kml_object, function(x) {
    criterion <- x["criterionActif"]
    allCrit <- x["criterionValues", criterion]
    nbCriterion <- 1000
    lengthList <- max(vapply(allCrit, length, integer(1)))
    allCrit <- sapply(allCrit, function(y) c(y, rep(NA, lengthList - length(y) + 1)))
    lengthList <- min(lengthList, nbCriterion)
    allCrit_df <- apply(allCrit, 2, unlist) |> as.data.frame(stringsAsFactors = FALSE)
    base::colnames(allCrit_df) <- gsub("c", "", base::colnames(allCrit_df))
    tmp_output <- allCrit_df[-base::nrow(allCrit_df), , drop = FALSE]

    if (!is.null(min_cluster_fraction)) {
      tmp <- vapply(seq_len(base::ncol(tmp_output)), function(y) {
        slot_name <- paste0("c", base::colnames(tmp_output)[[y]])
        part_list <- methods::slot(x, slot_name)
        base::min(base::table(part_list[[1]]@clusters))
      }, numeric(1))
      for (y in seq_len(base::ncol(tmp_output))) {
        if (tmp[[y]] < length(x@idFewNA) * min_cluster_fraction) {
          tmp_output[, y] <- NA
        }
      }
    }
    tmp_output
  })

  output <- list()
  output[["max"]] <- lapply(kml_criterion_df_list, function(x) {
    apply(x, 2, max, na.rm = TRUE) |> sapply(function(y) y * -1) |> rank()
  }) |> dplyr::bind_rows() |> as.data.frame(stringsAsFactors = FALSE) |> magrittr::set_rownames(base::names(kml_object))
  output[["mean"]] <- lapply(kml_criterion_df_list, function(x) {
    apply(x, 2, mean, na.rm = TRUE) |> sapply(function(y) y * -1) |> rank()
  }) |> dplyr::bind_rows() |> as.data.frame(stringsAsFactors = FALSE) |> magrittr::set_rownames(base::names(kml_object))
  output[["median"]] <- lapply(kml_criterion_df_list, function(x) {
    apply(x, 2, stats::median, na.rm = TRUE) |> sapply(function(y) y * -1) |> rank()
  }) |> dplyr::bind_rows() |> as.data.frame(stringsAsFactors = FALSE) |> magrittr::set_rownames(base::names(kml_object))

  kml_criterion_df_list <- lapply(kml_criterion_df_list, function(x) {
    base::t(base::scale(base::t(x)))
  })

  output[["z_max"]] <- lapply(kml_criterion_df_list, function(x) {
    apply(x, 2, max, na.rm = TRUE) |> sapply(function(y) y * -1) |> rank()
  }) |> dplyr::bind_rows() |> as.data.frame(stringsAsFactors = FALSE) |> magrittr::set_rownames(base::names(kml_object))
  output[["z_mean"]] <- lapply(kml_criterion_df_list, function(x) {
    apply(x, 2, mean, na.rm = TRUE) |> sapply(function(y) y * -1) |> rank()
  }) |> dplyr::bind_rows() |> as.data.frame(stringsAsFactors = FALSE) |> magrittr::set_rownames(base::names(kml_object))
  output[["z_median"]] <- lapply(kml_criterion_df_list, function(x) {
    apply(x, 2, stats::median, na.rm = TRUE) |> sapply(function(y) y * -1) |> rank()
  }) |> dplyr::bind_rows() |> as.data.frame(stringsAsFactors = FALSE) |> magrittr::set_rownames(base::names(kml_object))
  output
}

.hc_legacy_get_clustering <- function(kml_object,
                                      score_method = "median",
                                      score) {
  if (!score_method %in% base::names(score)) {
    stop("Unknown `score_method`: ", score_method, call. = FALSE)
  }
  score_tbl <- as.data.frame(score[[score_method]], stringsAsFactors = FALSE)
  k_selected <- apply(score_tbl, 1, function(x) {
    base::names(x)[base::match(base::min(x), x)]
  })

  kml_criterion_list <- lapply(seq_along(kml_object), function(i) {
    criterion <- kml_object[[i]]["criterionActif"]
    allCrit <- kml_object[[i]]["criterionValues", criterion]
    nbCriterion <- 1000
    lengthList <- max(vapply(allCrit, length, integer(1)))
    allCrit <- sapply(allCrit, function(x) c(x, rep(NA, lengthList - length(x) + 1)))
    lengthList <- min(lengthList, nbCriterion)
    allCrit_long <- apply(allCrit, 2, unlist) |> as.data.frame(stringsAsFactors = FALSE)
    base::colnames(allCrit_long) <- gsub("c", "", base::colnames(allCrit_long))
    allCrit_long[-base::nrow(allCrit_long), k_selected[[i]], drop = TRUE]
  })

  position_max_k <- lapply(kml_criterion_list, function(x) {
    x_num <- suppressWarnings(base::as.numeric(x))
    if (base::length(x_num) == 0) {
      return(1L)
    }
    ok <- base::is.finite(x_num)
    if (!base::any(ok)) {
      return(1L)
    }
    base::which.max(replace(x_num, !ok, -Inf))
  }) |> unlist()
  base::names(position_max_k) <- base::names(kml_object)

  kml_clusters <- lapply(seq_along(kml_object), function(i) {
    slot_name <- paste0("c", k_selected[[i]])
    part_list <- methods::slot(kml_object[[i]], slot_name)
    if (base::length(part_list) == 0) {
      candidate_slots <- paste0("c", base::colnames(score_tbl))
      candidate_slots <- candidate_slots[candidate_slots %in% methods::slotNames(kml_object[[i]])]
      candidate_lengths <- vapply(candidate_slots, function(nm) {
        base::length(methods::slot(kml_object[[i]], nm))
      }, integer(1))
      candidate_slots <- candidate_slots[candidate_lengths > 0]
      if (base::length(candidate_slots) == 0) {
        stop(
          "No non-empty kml partition slots found for module `",
          base::names(kml_object)[[i]],
          "`.",
          call. = FALSE
        )
      }
      slot_name <- candidate_slots[[1]]
      part_list <- methods::slot(kml_object[[i]], slot_name)
    }
    pos <- suppressWarnings(as.integer(position_max_k[[i]]))
    if (!base::is.finite(pos) || pos < 1L) {
      pos <- 1L
    }
    pos <- min(pos, base::length(part_list))
    suppressWarnings(as.numeric(part_list[[pos]]@clusters))
  }) |>
    dplyr::bind_cols() |>
    as.data.frame(stringsAsFactors = FALSE) |>
    magrittr::set_rownames(kml_object[[1]]@idFewNA) |>
    magrittr::set_colnames(base::names(kml_object))

  list(clusters = kml_clusters, selected_k = k_selected, selected_run = position_max_k)
}

.hc_legacy_cap_matrix <- function(kml_object, clustering_df) {
  output <- lapply(base::names(kml_object), function(mod) {
    kk <- base::max(clustering_df[[mod]], na.rm = TRUE)
    slot_name <- paste0("c", kk)
    part_list <- methods::slot(kml_object[[mod]], slot_name)
    if (base::length(part_list) == 0) {
      candidate_slots <- grep("^c[0-9]+$", methods::slotNames(kml_object[[mod]]), value = TRUE)
      candidate_lengths <- vapply(candidate_slots, function(nm) {
        base::length(methods::slot(kml_object[[mod]], nm))
      }, integer(1))
      candidate_slots <- candidate_slots[candidate_lengths > 0]
      if (base::length(candidate_slots) == 0) {
        stop("No non-empty kml partition slots found for module `", mod, "`.", call. = FALSE)
      }
      slot_name <- candidate_slots[[1]]
      part_list <- methods::slot(kml_object[[mod]], slot_name)
      kk <- suppressWarnings(as.integer(sub("^c", "", slot_name)))
    }
    runs <- base::length(part_list)
    mat_x <- do.call(base::rbind, lapply(part_list, function(y) {
      y@clusters
    }))
    likelihood <- apply(mat_x, MARGIN = 2, FUN = function(t) {
      base::table(base::factor(t, levels = base::seq_len(kk)))
    }) / runs
    likelihood <- base::t(likelihood)
    likelihood <- as.data.frame(likelihood, stringsAsFactors = FALSE)
    base::colnames(likelihood) <- base::paste0(mod, "__", base::seq_len(kk))
    likelihood
  })
  out <- dplyr::bind_cols(output)
  base::rownames(out) <- kml_object[[1]]@idFewNA
  out
}

.hc_legacy_group_singletons <- function(ids,
                                        snn,
                                        group.singletons = TRUE,
                                        verbose = TRUE) {
  singletons <- base::names(base::which(base::table(ids) == 1))
  singletons <- base::intersect(base::unique(ids), singletons)
  if (!group.singletons) {
    ids[base::which(ids %in% singletons)] <- "singleton"
    return(ids)
  }

  cluster_names <- base::as.character(base::unique(ids))
  cluster_names <- base::setdiff(cluster_names, singletons)
  connectivity <- base::numeric(base::length(cluster_names))
  base::names(connectivity) <- cluster_names

  for (i in singletons) {
    i.cells <- base::names(base::which(ids == i))
    for (j in cluster_names) {
      j.cells <- base::names(base::which(ids == j))
      subSNN <- snn[i.cells, j.cells, drop = FALSE]
      set.seed(1)
      if (methods::is(subSNN, "Matrix")) {
        connectivity[[j]] <- base::sum(subSNN) / (base::nrow(subSNN) * base::ncol(subSNN))
      } else {
        connectivity[[j]] <- base::mean(subSNN)
      }
    }
    m <- base::max(connectivity, na.rm = TRUE)
    mi <- base::which(connectivity == m)
    closest_cluster <- sample(base::names(connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (base::length(singletons) > 0 && isTRUE(verbose)) {
    message(
      base::length(singletons), " singletons identified. ",
      base::length(base::unique(ids)), " final clusters."
    )
  }
  ids
}

.hc_legacy_build_knn_graph <- function(knn_data,
                                       dimensions = 4,
                                       graph_method = c("knn", "snn"),
                                       knn_method = c("annoy", "rann"),
                                       graph_k = 7) {
  graph_method <- base::match.arg(graph_method)
  knn_method <- base::match.arg(knn_method)

  dims_use <- seq_len(min(as.integer(dimensions), base::ncol(knn_data)))
  x <- knn_data[, dims_use, drop = FALSE]
  if (is.null(graph_k) || !is.finite(graph_k) || graph_k <= 0) {
    graph_k <- ifelse(base::sqrt(base::nrow(x)) < 10, 10, base::floor(base::sqrt(base::nrow(x))))
  }
  graph_k <- as.integer(graph_k)

  if (identical(knn_method, "rann")) {
    .hc_require_namespace("RANN", "legacy KNN graph construction")
    nn_result <- RANN::nn2(x, k = graph_k)
    neighbor_indices <- nn_result$nn.idx
  } else {
    .hc_require_namespace("RcppAnnoy", "legacy Annoy KNN graph construction")
    ann_index <- RcppAnnoy::AnnoyAngular$new(base::ncol(x))
    for (i in seq_len(base::nrow(x))) {
      ann_index$addItem(i - 1L, x[i, ])
    }
    ann_index$build(50)
    idx <- base::matrix(nrow = base::nrow(x), ncol = graph_k)
    search_k <- 100L * graph_k
    for (i in seq_len(base::nrow(x))) {
      annoy_res <- ann_index$getNNsByVectorList(x[i, ], graph_k, search_k, TRUE)
      if (base::length(annoy_res$item) != graph_k) {
        stop("Annoy search failed to find the requested number of neighbors.", call. = FALSE)
      }
      idx[i, ] <- annoy_res$item
    }
    neighbor_indices <- idx + 1L
  }

  if (identical(graph_method, "knn")) {
    graph <- base::matrix(0, nrow = base::nrow(x), ncol = base::nrow(x))
    for (i in seq_len(base::nrow(x))) {
      graph[i, neighbor_indices[i, ]] <- 1
    }
    base::rownames(graph) <- base::rownames(knn_data)
    base::colnames(graph) <- base::rownames(knn_data)
    return(graph)
  }

  graph <- .hc_legacy_compute_snn(
    nn_ranked = neighbor_indices[seq_len(base::nrow(x)), , drop = FALSE],
    prune = as.double(1 / 15)
  )
  base::rownames(graph) <- base::rownames(knn_data)
  base::colnames(graph) <- base::rownames(knn_data)
  graph
}

# Seurat computes SNN weights as shared-neighbor overlap scaled by
# overlap / (k + (k - overlap)) and prunes low-weight edges.
.hc_legacy_compute_snn <- function(nn_ranked, prune = 1 / 15) {
  k <- base::ncol(nn_ranked)
  j <- base::as.numeric(base::t(nn_ranked))
  i <- ((seq_along(j) - 1L) %/% k) + 1L
  nn_matrix <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(base::nrow(nn_ranked), base::nrow(nn_ranked))
  )
  snn_matrix <- nn_matrix %*% Matrix::t(nn_matrix)
  snn_matrix@x <- snn_matrix@x / (k + (k - snn_matrix@x))
  snn_matrix@x[snn_matrix@x < prune] <- 0
  base::as.matrix(Matrix::drop0(snn_matrix))
}

.hc_legacy_compute_umap <- function(x,
                                    donor_ids,
                                    labels,
                                    seed = 42,
                                    umap_neighbors = 15,
                                    umap_min_dist = 0.3) {
  meta_umap <- NULL
  if (requireNamespace("uwot", quietly = TRUE)) {
    set.seed(seed)
    um <- uwot::umap(
      x,
      n_neighbors = as.integer(umap_neighbors),
      min_dist = as.numeric(umap_min_dist),
      metric = "euclidean",
      verbose = FALSE
    )
    meta_umap <- base::data.frame(
      donor = donor_ids,
      UMAP1 = um[, 1],
      UMAP2 = um[, 2],
      meta_cluster = labels,
      stringsAsFactors = FALSE
    )
  } else if (requireNamespace("umap", quietly = TRUE)) {
    set.seed(seed)
    um <- umap::umap(x)
    meta_umap <- base::data.frame(
      donor = donor_ids,
      UMAP1 = um$layout[, 1],
      UMAP2 = um$layout[, 2],
      meta_cluster = labels,
      stringsAsFactors = FALSE
    )
  } else {
    warning("UMAP not computed: neither `uwot` nor `umap` is installed.")
  }
  meta_umap
}

.hc_run_legacy_step1_exact <- function(hc,
                                       means_slot = "longitudinal_module_means",
                                       output_slot = "longitudinal_endotypes",
                                       k = 2:6,
                                       rerolls = 1000,
                                       impute = TRUE,
                                       impute_method = "rfcont",
                                       ntree = 10,
                                       min_cluster_fraction = 0.1,
                                       score_method = "median",
                                       seed = 42) {
  .hc_require_namespace("kml", "legacy longitudinal step 1")
  sat <- as.list(hc@satellite)
  inp <- sat[[means_slot]]
  if (is.null(inp) || is.null(inp$module_lookup) || is.null(inp$layer_id)) {
    stop("Step 1 requires outputs from `hc_longitudinal_module_means()`.", call. = FALSE)
  }

  module_lookup <- base::as.data.frame(inp$module_lookup, stringsAsFactors = FALSE)
  layer_id <- as.character(inp$layer_id[[1]])
  donor_col <- as.character(inp$donor_col[[1]])
  time_col <- as.character(inp$time_col[[1]])
  time_levels <- base::as.character(inp$time_levels)
  module_ids <- base::as.character(module_lookup$module)
  module_colors <- base::as.character(module_lookup$module_color)

  se <- MultiAssayExperiment::experiments(hc@mae)[[layer_id]]
  counts <- SummarizedExperiment::assay(se, "counts")
  counts <- base::as.matrix(counts)
  anno <- base::as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  if (is.null(base::rownames(anno))) {
    base::rownames(anno) <- base::colnames(counts)
  }
  sample_ids <- base::colnames(counts)
  anno <- anno[base::match(sample_ids, base::rownames(anno)), , drop = FALSE]
  sample_table <- anno
  sample_table$ID <- sample_ids

  cluster_obj <- as.list(hc@integration@cluster)
  cluster_info <- cluster_obj[["cluster_information"]]
  if (is.null(cluster_info)) {
    stop("No cluster information found. Run `hc_cluster_calculation()` first.", call. = FALSE)
  }
  cluster_info <- base::as.data.frame(cluster_info, stringsAsFactors = FALSE)
  if ("cluster_included" %in% base::colnames(cluster_info)) {
    cluster_info <- dplyr::filter(cluster_info, cluster_included == "yes")
  }
  cluster_map <- .hc_module_gene_map(cluster_info = cluster_info)

  gene_groups <- base::do.call(
    base::rbind,
    base::lapply(base::seq_len(base::nrow(module_lookup)), function(i) {
      mcol <- module_colors[[i]]
      genes <- cluster_map[[mcol]]
      if (base::length(genes) == 0) {
        return(NULL)
      }
      base::data.frame(
        color = mcol,
        gene = genes,
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(gene_groups) || base::nrow(gene_groups) == 0) {
    stop("No module genes available for legacy longitudinal step 1.", call. = FALSE)
  }

  FCs <- NULL
  for (mcol in module_colors) {
    genes <- gene_groups[gene_groups$color == mcol, "gene", drop = TRUE]
    tmp_counts <- counts[base::rownames(counts) %in% genes, , drop = FALSE]
    tmp_fc <- apply(tmp_counts, 2, base::mean) |> t() |> as.data.frame(stringsAsFactors = FALSE)
    base::rownames(tmp_fc) <- mcol
    FCs <- base::rbind(FCs, tmp_fc)
  }
  base::colnames(FCs) <- base::colnames(counts)
  FCs$module <- base::rownames(FCs)
  FC_long <- as.data.frame(base::t(FCs), stringsAsFactors = FALSE) |>
    tibble::rownames_to_column(var = "ID") |>
    dplyr::select("ID", dplyr::everything()) |>
    tidyr::pivot_longer(cols = !"ID", names_to = "module", values_to = "GFC")
  sample_FC <- base::merge(FC_long, sample_table, by = "ID") |>
    dplyr::mutate(GFC = as.numeric(GFC))

  if (!is.numeric(sample_FC[[time_col]])) {
    if (is.factor(sample_FC[[time_col]])) {
      sample_FC[[time_col]] <- as.character(as.numeric(sample_FC[[time_col]]))
    } else {
      stop("`", time_col, "` is neither numeric nor factor in the selected layer annotation.", call. = FALSE)
    }
  }

  donor_time_module <- stats::aggregate(
    stats::as.formula(base::paste0("GFC ~ ", donor_col, " + ", time_col, " + module")),
    data = sample_FC,
    FUN = base::mean
  )
  base::colnames(donor_time_module) <- c("donor", "time", "module_color", "value")
  donor_time_module$module <- module_lookup$module[
    base::match(base::as.character(donor_time_module$module_color), module_lookup$module_color)
  ]
  donor_time_module <- donor_time_module[!base::is.na(donor_time_module$module), , drop = FALSE]
  donor_time_module <- donor_time_module[, c("donor", "time", "module", "value"), drop = FALSE]

  kml_object <- list()
  for (i in base::seq_len(base::nrow(module_lookup))) {
    mod <- module_ids[[i]]
    mcol <- module_colors[[i]]
    tmp_sample_fc <- sample_FC[sample_FC$module == mcol, , drop = FALSE]
    time_data <- reshape2::dcast(
      tmp_sample_fc,
      stats::as.formula(base::paste0(donor_col, "~", time_col)),
      value.var = "GFC"
    )
    base::colnames(time_data) <- base::as.character(base::colnames(time_data))
    keep_cols <- c(donor_col, time_levels)
    if (!base::all(keep_cols %in% base::colnames(time_data))) {
      missing_tp <- base::setdiff(keep_cols, base::colnames(time_data))
      for (tp in missing_tp) {
        time_data[[tp]] <- NA_real_
      }
      time_data <- time_data[, keep_cols, drop = FALSE]
    } else {
      time_data <- time_data[, keep_cols, drop = FALSE]
    }
    if (isTRUE(impute)) {
      time_data <- .hc_legacy_impute_time_data(
        time_data = time_data,
        donor_col = donor_col,
        method = impute_method,
        ntree = ntree,
        seed = seed
      )
    }
    base::colnames(time_data)[base::colnames(time_data) == donor_col] <- "donor"
    base::colnames(time_data)[-1] <- as.character(seq_len(base::length(time_levels)))
    kml_object[[mod]] <- .hc_legacy_kml_cluster_one(
      time_data = time_data,
      k = k,
      rerolls = rerolls
    )
  }

  master_donors <- kml_object[[1]]@idFewNA

  score <- .hc_legacy_kml_cluster_scoring(
    kml_object = kml_object,
    min_cluster_fraction = min_cluster_fraction
  )
  clustering <- .hc_legacy_get_clustering(
    kml_object = kml_object,
    score_method = score_method,
    score = score
  )

  module_cluster_matrix <- base::matrix(
    NA_integer_,
    nrow = base::length(master_donors),
    ncol = base::length(module_ids),
    dimnames = list(master_donors, module_ids)
  )
  common_rows <- base::intersect(master_donors, base::rownames(clustering$clusters))
  common_cols <- base::intersect(module_ids, base::colnames(clustering$clusters))
  module_cluster_matrix[common_rows, common_cols] <- base::as.matrix(
    clustering$clusters[common_rows, common_cols, drop = FALSE]
  )
  module_cluster_matrix <- base::as.data.frame(module_cluster_matrix, stringsAsFactors = FALSE)

  module_best_k <- base::data.frame(
    module = base::names(clustering$selected_k),
    best_k = as.integer(clustering$selected_k),
    selected_run = as.integer(clustering$selected_run[base::names(clustering$selected_k)]),
    stringsAsFactors = FALSE
  )
  module_best_k <- module_best_k[base::match(module_ids, module_best_k$module), , drop = FALSE]

  module_cluster_palette <- .hc_module_cluster_palette(
    module_lookup = module_lookup,
    module_best_k = module_best_k[, c("module", "best_k"), drop = FALSE]
  )

  module_cluster_assignments <- base::as.data.frame(
    base::as.table(base::as.matrix(module_cluster_matrix)),
    stringsAsFactors = FALSE
  )
  base::colnames(module_cluster_assignments) <- c("donor", "module", "cluster")
  module_cluster_assignments$cluster <- suppressWarnings(as.integer(module_cluster_assignments$cluster))
  module_cluster_assignments <- module_cluster_assignments[!base::is.na(module_cluster_assignments$cluster), , drop = FALSE]
  module_cluster_assignments$module_cluster_key <- base::paste0(
    module_cluster_assignments$module, "__", module_cluster_assignments$cluster
  )

  donor_time_module_with_clusters <- base::merge(
    donor_time_module,
    module_cluster_assignments[, c("donor", "module", "cluster", "module_cluster_key"), drop = FALSE],
    by = c("donor", "module"),
    all.x = FALSE,
    sort = FALSE
  )

  module_cluster_trajectory_mean <- stats::aggregate(
    value ~ module + time + module_cluster_key,
    data = donor_time_module_with_clusters,
    FUN = base::mean
  )

  cap_matrix <- .hc_legacy_cap_matrix(
    kml_object = kml_object,
    clustering_df = module_cluster_matrix
  )
  cap_matrix <- base::as.matrix(cap_matrix)
  cap_matrix <- cap_matrix[master_donors, , drop = FALSE]

  cap_long <- base::as.data.frame(base::as.table(cap_matrix), stringsAsFactors = FALSE)
  base::colnames(cap_long) <- c("donor", "module_cluster_key", "proportion")
  cap_long$module <- sub("__.*$", "", base::as.character(cap_long$module_cluster_key))
  cap_long$cluster <- suppressWarnings(as.integer(sub("^.*__", "", base::as.character(cap_long$module_cluster_key))))
  cap_long$proportion <- as.numeric(cap_long$proportion)

  cap_consensus_assignments <- stats::aggregate(
    proportion ~ donor + module,
    data = cap_long,
    FUN = max
  )
  base::colnames(cap_consensus_assignments)[base::colnames(cap_consensus_assignments) == "proportion"] <- "max_prop"

  sat[[output_slot]] <- list(
    created_at = as.character(base::Sys.time()),
    source_slot = means_slot,
    donor_col = inp$donor_col,
    time_col = inp$time_col,
    group_col = inp$group_col,
    time_levels = time_levels,
    module_lookup = module_lookup,
    donor_time_module = donor_time_module,
    value_label = inp$value_label,
    module_cluster_matrix = base::as.matrix(module_cluster_matrix),
    module_cluster_assignments = module_cluster_assignments,
    donor_time_module_with_module_clusters = donor_time_module_with_clusters,
    module_cluster_trajectory_mean = module_cluster_trajectory_mean,
    module_cluster_palette = module_cluster_palette,
    module_cluster_best_k = module_best_k,
    module_cluster_score = score,
    module_cluster_score_method = score_method,
    module_cluster_min_fraction = min_cluster_fraction,
    legacy_kml_object = kml_object,
    legacy_master_donors = master_donors,
    legacy_kml_rerolls = as.integer(rerolls),
    legacy_kml_k = as.integer(k),
    legacy_kml_impute = isTRUE(impute),
    legacy_kml_impute_method = impute_method,
    legacy_kml_ntree = as.integer(ntree),
    cap_matrix = cap_matrix,
    cap_long = cap_long,
    cap_consensus_assignments = cap_consensus_assignments,
    cap_method = "legacy_exact"
  )
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}

.hc_run_legacy_step2_exact <- function(hc,
                                       slot_name = "longitudinal_endotypes",
                                       dimensions = 4,
                                       graph_method = c("knn", "snn"),
                                       knn_method = c("annoy", "rann"),
                                       graph_k = 7,
                                       resolution = 0.4,
                                       leiden_method = "RBConfigurationVertexPartition",
                                       cluster_prefix = "MC",
                                       seed = 42,
                                       compute_umap = TRUE,
                                       umap_neighbors = 15,
                                       umap_min_dist = 0.3) {
  .hc_require_namespace("igraph", "legacy longitudinal step 2")
  .hc_require_namespace("leidenbase", "legacy graph Leiden clustering")

  graph_method <- base::match.arg(graph_method)
  knn_method <- base::match.arg(knn_method)

  sat <- as.list(hc@satellite)
  obj <- sat[[slot_name]]
  if (is.null(obj) || is.null(obj$cap_matrix)) {
    stop("Step 2 requires `cap_matrix` from step 1.", call. = FALSE)
  }

  cap_matrix <- base::as.matrix(obj$cap_matrix)
  donor_ids <- base::rownames(cap_matrix)
  if (base::is.null(donor_ids) || base::length(donor_ids) < 3) {
    stop("Need at least 3 donors in the CAP matrix for meta-clustering.", call. = FALSE)
  }

  pca_fit <- stats::prcomp(cap_matrix)
  meta_pca_var <- 100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))
  knn_data <- pca_fit$x
  graph <- .hc_legacy_build_knn_graph(
    knn_data = knn_data,
    dimensions = dimensions,
    graph_method = graph_method,
    knn_method = knn_method,
    graph_k = graph_k
  )

  igraph_obj <- igraph::graph.adjacency(graph, mode = "undirected", weighted = TRUE, diag = FALSE)
  leiden_clusters <- leidenbase::leiden_find_partition(
    igraph_obj,
    leiden_method,
    resolution_parameter = resolution,
    seed = seed,
    num_iter = 10
  )$membership
  base::names(leiden_clusters) <- base::colnames(graph)
  leiden_clusters <- .hc_legacy_group_singletons(
    ids = leiden_clusters,
    snn = graph,
    group.singletons = TRUE,
    verbose = TRUE
  )
  leiden_clusters <- as.integer(leiden_clusters[donor_ids])
  meta_labels <- base::paste0(cluster_prefix, leiden_clusters)

  meta_cluster <- base::data.frame(
    donor = donor_ids,
    cluster = leiden_clusters,
    meta_cluster = meta_labels,
    stringsAsFactors = FALSE
  )

  meta_pca <- base::data.frame(
    donor = donor_ids,
    PC1 = pca_fit$x[, 1],
    PC2 = if (base::ncol(pca_fit$x) >= 2) pca_fit$x[, 2] else base::rep(0, base::nrow(pca_fit$x)),
    meta_cluster = meta_labels,
    stringsAsFactors = FALSE
  )

  dims_use <- seq_len(min(as.integer(dimensions), base::ncol(knn_data)))
  meta_umap <- NULL
  if (isTRUE(compute_umap)) {
    meta_umap <- .hc_legacy_compute_umap(
      x = knn_data[, dims_use, drop = FALSE],
      donor_ids = donor_ids,
      labels = meta_labels,
      seed = seed,
      umap_neighbors = umap_neighbors,
      umap_min_dist = umap_min_dist
    )
  }

  n_meta <- base::length(base::unique(meta_labels))
  obj$meta_feature_source <- "cap_matrix"
  obj$meta_feature_matrix_used <- cap_matrix
  obj$meta_method <- "legacy_graph_leiden"
  obj$meta_candidate_k <- as.integer(n_meta)
  obj$meta_nstart <- NA_integer_
  obj$meta_seed <- as.integer(seed)
  obj$meta_score_table <- base::data.frame(
    k = as.integer(n_meta),
    ch_index = NA_real_,
    graph_method = graph_method,
    knn_method = knn_method,
    graph_k = as.integer(graph_k),
    resolution = as.numeric(resolution),
    dimensions = as.integer(dimensions),
    stringsAsFactors = FALSE
  )
  obj$meta_best_k <- as.integer(n_meta)
  obj$meta_cluster <- meta_cluster
  obj$meta_donor_clusterings <- base::data.frame(
    donor = donor_ids,
    method = "graph_leiden",
    best_k_ch = as.integer(n_meta),
    cluster_id = leiden_clusters,
    cluster_label = meta_labels,
    stringsAsFactors = FALSE
  )
  obj$meta_pca <- meta_pca
  obj$meta_pca_variance <- meta_pca_var
  obj$meta_umap <- meta_umap
  obj$meta_method_comparison <- base::data.frame(
    method = "graph_leiden",
    best_k_ch = as.integer(n_meta),
    best_ch_index = NA_real_,
    n_clusters_at_best_k = as.integer(n_meta),
    status = "ok",
    selected_method = TRUE,
    selected_k_pipeline = as.integer(n_meta),
    n_meta_clusters_pipeline = as.integer(n_meta),
    selection_mode = "graph_leiden",
    graph_method = graph_method,
    knn_method = knn_method,
    graph_k = as.integer(graph_k),
    resolution = as.numeric(resolution),
    dimensions = as.integer(dimensions),
    stringsAsFactors = FALSE
  )
  obj$meta_consensus <- FALSE
  obj$meta_consensus_matrix <- NULL
  obj$meta_consensus_runs <- NA_integer_
  obj$meta_consensus_sample_fraction <- NA_real_
  obj$meta_consensus_feature_fraction <- NA_real_
  obj$meta_consensus_linkage <- NA_character_
  obj$meta_donor_stability <- NULL
  obj$meta_graph_method <- graph_method
  obj$meta_knn_method <- knn_method
  obj$meta_graph_k <- as.integer(graph_k)
  obj$meta_resolution <- as.numeric(resolution)
  obj$meta_dimensions <- as.integer(dimensions)
  obj$meta_leiden_method <- leiden_method

  sat[[slot_name]] <- obj
  hc@satellite <- S4Vectors::SimpleList(sat)
  methods::validObject(hc)
  hc
}
