#' Merge Clusters
#' 
#' Merge similar clusters in module heatmap.
#' @param k Either an integer (number of merged modules) or `"auto"` to choose
#'   a data-driven value based on the average silhouette score across candidate
#'   cuts of the module dendrogram.
#' @param save Boolean. If FALSE (default) cutting is showcased but not saved.
#'   Set to TRUE once you have found the right clustering and overwrite the old
#'   cluster assignment (cannot be reversed).
#' @param method Method used for clustering. Default is "complete".
#' @param k_min Minimum candidate `k` for automatic selection (`k = "auto"`).
#' @param k_max Maximum candidate `k` for automatic selection (`k = "auto"`).
#'   Defaults to `min(10, n_modules - 1)`.
#' @param auto_parsimony_penalty Small penalty used in auto mode to prefer
#'   simpler solutions when silhouette scores are very similar.
#' @param verbose Boolean. Print selected `k` and selection diagnostics.
#' @export

merge_clusters <- function(k = 1,
                           save = FALSE,
                           method = "complete",
                           k_min = 2,
                           k_max = NULL,
                           auto_parsimony_penalty = 1e-04,
                           verbose = TRUE) {
  m <- hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]]@ht_list[[1]]@matrix
  m <- as.matrix(m)
  if (nrow(m) < 2) {
    stop("Need at least 2 modules to merge.")
  }

  auto_mode <- is.null(k) || (is.character(k) && length(k) == 1 && tolower(k) == "auto")
  if (is.numeric(k) && length(k) == 1 && is.na(k)) {
    auto_mode <- TRUE
  }

  selected_k <- NA_integer_
  auto_report <- NULL

  if (isTRUE(auto_mode)) {
    auto_report <- .hc_select_merge_k_auto(
      m = m,
      method = method,
      k_min = k_min,
      k_max = k_max,
      parsimony_penalty = auto_parsimony_penalty
    )
    selected_k <- auto_report$selected_k
    if (isTRUE(verbose)) {
      message(
        "merge_clusters(): auto-selected k = ", selected_k,
        " (criterion: max adjusted mean silhouette)."
      )
    }
  } else {
    if (!is.numeric(k) || length(k) != 1 || !is.finite(k)) {
      stop("`k` must be an integer >= 1 or \"auto\".")
    }
    selected_k <- as.integer(round(k))
    if (selected_k < 1 || selected_k > nrow(m)) {
      stop("`k` must be between 1 and the number of modules (", nrow(m), ").")
    }
  }

  p <- pheatmap::pheatmap(
    mat = m,
    color = base::rev(RColorBrewer::brewer.pal(11, "RdBu")),
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize = 8,
    show_rownames = FALSE,
    show_colnames = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_method = method,
    cutree_rows = selected_k
  )

  if (isTRUE(save)) {
    cut <- stats::cutree(p$tree_row, k = selected_k)
    new_group <- base::data.frame(old_cluster = base::names(cut), group = cut)
    new_group$new_cluster <- base::lapply(new_group$group, function(x){
      get_cluster_colours()[x]
    }) %>% base::unlist()
    
    gtc <- GeneToCluster()
    base::colnames(gtc) <- base::c("gene", "old_color")
    gtc$color <- base::lapply(gtc$old_color, function(x){
      nc <- dplyr::filter(new_group, as.character(old_cluster) == x) %>% dplyr::pull(., "new_cluster")
      
      return(nc)
    }) %>% base::unlist()
    
    new_cluster_info <- NULL
    
    for(c in base::unique(gtc$color)){
      tmp <- gtc[gtc$color == c, ]
      gene_no <- base::nrow(tmp)
      gene_n <- base::paste0(tmp$gene, collapse = ",")
      if(c == "white"){
        cluster_included <- "no"
        vertexsize <- 1
      }else{
        cluster_included <- "yes"
        vertexsize <- 3
      }
      
      color <- c
      conditions <- base::paste0(base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])[1:(base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]])-1)], collapse = "#")
      gfc_means = hcobject[["integrated_output"]][["GFC_all_layers"]][hcobject[["integrated_output"]][["GFC_all_layers"]][["Gene"]] %in% tmp$gene,] %>%
        dplyr::select(-Gene) %>%
        base::colMeans()
      
      grp_means = base::paste0(base::round(gfc_means,3) , collapse = ",")
      
      new_cluster_info <- rbind(new_cluster_info,
                                data.frame(clusters = c,
                                           gene_no = gene_no,
                                           gene_n = gene_n,
                                           cluster_included = cluster_included,
                                           color = c,
                                           conditions = conditions,
                                           grp_means = grp_means,
                                           vertexsize = vertexsize,
                                           stringsAsFactors = F))
      
    }
    hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]] <<- new_cluster_info
    plot_cluster_heatmap()

    return(invisible(list(
      selected_k = selected_k,
      auto_summary = if (isTRUE(auto_mode)) auto_report$summary else NULL
    )))
  } else {
    message("This is only a preview. To save these new clusters, set `save = TRUE`.")
    return(invisible(list(
      selected_k = selected_k,
      auto_summary = if (isTRUE(auto_mode)) auto_report$summary else NULL,
      plot = p
    )))
  }
}

.hc_select_merge_k_auto <- function(m,
                                    method = "complete",
                                    k_min = 2,
                                    k_max = NULL,
                                    parsimony_penalty = 1e-04) {
  n <- nrow(m)
  if (n < 3) {
    return(list(
      selected_k = 1L,
      summary = data.frame(
        k = 1L,
        mean_silhouette = NA_real_,
        singleton_fraction = NA_real_,
        score = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  if (!is.numeric(k_min) || length(k_min) != 1 || !is.finite(k_min)) {
    stop("`k_min` must be a finite numeric scalar.")
  }
  if (!is.null(k_max) && (!is.numeric(k_max) || length(k_max) != 1 || !is.finite(k_max))) {
    stop("`k_max` must be NULL or a finite numeric scalar.")
  }
  if (!is.numeric(parsimony_penalty) || length(parsimony_penalty) != 1 || !is.finite(parsimony_penalty) || parsimony_penalty < 0) {
    stop("`auto_parsimony_penalty` must be a finite numeric scalar >= 0.")
  }

  k_min <- max(2L, as.integer(round(k_min)))
  if (is.null(k_max)) {
    k_max <- min(10L, n - 1L)
  } else {
    k_max <- as.integer(round(k_max))
  }
  k_max <- min(k_max, n - 1L)
  if (k_min > k_max) {
    k_min <- max(2L, min(k_max, n - 1L))
  }
  if (k_min > k_max) {
    return(list(
      selected_k = 2L,
      summary = data.frame(
        k = 2L,
        mean_silhouette = NA_real_,
        singleton_fraction = NA_real_,
        score = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  d <- stats::dist(m, method = "euclidean")
  dmat <- as.matrix(d)
  tree <- stats::hclust(d, method = method)
  k_grid <- seq.int(k_min, k_max)

  rows <- lapply(k_grid, function(k_now) {
    cl <- stats::cutree(tree, k = k_now)
    sil <- .hc_mean_silhouette_from_dist(dmat = dmat, cluster_id = cl)
    singleton_fraction <- mean(table(cl) <= 1)
    score <- sil - parsimony_penalty * (k_now - k_min)
    data.frame(
      k = k_now,
      mean_silhouette = sil,
      singleton_fraction = singleton_fraction,
      score = score,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$k), , drop = FALSE]

  valid <- is.finite(tab$score)
  if (!any(valid)) {
    return(list(selected_k = k_min, summary = tab))
  }

  best_score <- max(tab$score[valid], na.rm = TRUE)
  best_rows <- which(valid & abs(tab$score - best_score) < 1e-12)
  selected_k <- tab$k[min(best_rows)]

  list(selected_k = as.integer(selected_k), summary = tab)
}

.hc_mean_silhouette_from_dist <- function(dmat, cluster_id) {
  n <- length(cluster_id)
  if (n <= 1 || nrow(dmat) != n || ncol(dmat) != n) {
    return(NA_real_)
  }
  groups <- sort(unique(cluster_id))
  if (length(groups) < 2) {
    return(NA_real_)
  }

  s <- numeric(n)
  idx_all <- seq_len(n)
  for (i in idx_all) {
    own <- cluster_id[[i]]
    same <- idx_all[cluster_id == own & idx_all != i]

    if (length(same) == 0) {
      a_i <- 0
    } else {
      a_i <- mean(dmat[i, same])
    }

    b_i <- Inf
    for (g in groups[groups != own]) {
      other <- idx_all[cluster_id == g]
      if (length(other) > 0) {
        b_i <- min(b_i, mean(dmat[i, other]))
      }
    }
    if (!is.finite(b_i)) {
      b_i <- 0
    }

    den <- max(a_i, b_i)
    s[[i]] <- if (den > 0) (b_i - a_i) / den else 0
  }
  mean(s)
}
