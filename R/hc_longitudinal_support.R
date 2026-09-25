.hc_require_namespace <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package `", pkg, "` is required")
    if (!is.null(reason) && nzchar(reason)) {
      msg <- paste0(msg, " for ", reason)
    }
    stop(msg, ". Please install it first.", call. = FALSE)
  }
}

.hc_longitudinal_impute_time_data <- function(time_data,
                                              donor_col = "donor",
                                              method = "rfcont",
                                              ntree = 10,
                                              m = 50,
                                              maxit = 50,
                                              seed = 42) {
  .hc_require_namespace("mice", "longitudinal imputation")
  if (identical(method, "rfcont")) {
    .hc_require_namespace("CALIBERrfimpute", "`rfcont` longitudinal imputation")
    rf_pkg_search <- "package:CALIBERrfimpute"
    attached_here <- FALSE
    if (!(rf_pkg_search %in% search())) {
      base::attachNamespace(base::asNamespace("CALIBERrfimpute"))
      attached_here <- TRUE
    }
    if (isTRUE(attached_here)) {
      on.exit(
        detach(rf_pkg_search, unload = FALSE, character.only = TRUE),
        add = TRUE
      )
    }
    CALIBERrfimpute::setRFoptions(ntree_cont = ntree)
  }

  donor_ids <- time_data[[donor_col]]
  data_use <- time_data[, base::setdiff(base::colnames(time_data), donor_col), drop = FALSE]
  base::colnames(data_use) <- base::paste0("t", base::colnames(data_use))

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

.hc_longitudinal_group_singletons <- function(ids,
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

  # If every cluster is a singleton there is no multi-member cluster to merge
  # into. Bailing out here avoids `max(numeric(0))` (-Inf + warning) and the
  # subsequent `sample(character(0), 1)` error.
  if (base::length(cluster_names) == 0) {
    if (base::length(singletons) > 0 && isTRUE(verbose)) {
      message(
        base::length(singletons),
        " singletons identified but no multi-member cluster to join; left unchanged."
      )
    }
    return(ids)
  }

  connectivity <- base::numeric(base::length(cluster_names))
  base::names(connectivity) <- cluster_names

  for (i in singletons) {
    i.cells <- base::names(base::which(ids == i))
    for (j in cluster_names) {
      j.cells <- base::names(base::which(ids == j))
      subSNN <- snn[i.cells, j.cells, drop = FALSE]
      if (methods::is(subSNN, "Matrix")) {
        connectivity[[j]] <- base::sum(subSNN) / (base::nrow(subSNN) * base::ncol(subSNN))
      } else {
        connectivity[[j]] <- base::mean(subSNN)
      }
    }
    if (!base::any(base::is.finite(connectivity))) {
      # No usable connectivity to any cluster: leave this singleton as-is.
      next
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

.hc_longitudinal_build_knn_graph <- function(knn_data,
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
    .hc_require_namespace("RANN", "longitudinal KNN graph construction")
    nn_result <- RANN::nn2(x, k = graph_k)
    neighbor_indices <- nn_result$nn.idx
  } else {
    .hc_require_namespace("RcppAnnoy", "longitudinal Annoy KNN graph construction")
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

  graph <- .hc_longitudinal_compute_snn(
    nn_ranked = neighbor_indices[seq_len(base::nrow(x)), , drop = FALSE],
    prune = as.double(1 / 15)
  )
  base::rownames(graph) <- base::rownames(knn_data)
  base::colnames(graph) <- base::rownames(knn_data)
  graph
}

.hc_longitudinal_compute_snn <- function(nn_ranked, prune = 1 / 15) {
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

.hc_longitudinal_compute_umap <- function(x,
                                          donor_ids,
                                          labels,
                                          seed = 42,
                                          umap_neighbors = 15,
                                          umap_min_dist = 0.3) {
  meta_umap <- NULL
  if (requireNamespace("uwot", quietly = TRUE)) {
    um <- .hc_with_seed(seed, uwot::umap(
      x,
      n_neighbors = as.integer(umap_neighbors),
      min_dist = as.numeric(umap_min_dist),
      metric = "euclidean",
      verbose = FALSE
    ))
    meta_umap <- base::data.frame(
      donor = donor_ids,
      UMAP1 = um[, 1],
      UMAP2 = um[, 2],
      meta_cluster = labels,
      stringsAsFactors = FALSE
    )
  } else if (requireNamespace("umap", quietly = TRUE)) {
    um <- .hc_with_seed(seed, umap::umap(x))
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

.hc_run_longitudinal_step2_graph <- function(hc,
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
  .hc_require_namespace("igraph", "longitudinal step 2 graph clustering")
  .hc_require_namespace("leidenbase", "graph Leiden clustering")

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
  graph <- .hc_longitudinal_build_knn_graph(
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
  leiden_clusters <- .hc_longitudinal_group_singletons(
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
    meta_umap <- .hc_longitudinal_compute_umap(
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
  obj$meta_method <- "graph_leiden"
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
