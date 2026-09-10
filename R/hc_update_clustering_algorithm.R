#' Change Or Update The Clustering Algorithm
#'
#' The clustering can be updated to another algortihm (one of "cluster_louvain", "cluster_fast_greedy", "cluster_infomap", "cluster_walktrap", "cluster_label_prop", "cluster_leiden") or to a user defined clsutering using the 'gtc' parameter.
#' @param new_algo A string (name of the new algorithm), or NULL (if gtc is provided)
#' @param gtc A user defined clsutering can be provided as long as it has the same format as the output of .hc_gene_to_cluster_impl(): A data frame with two columns, the first containing gene names as strings, the second containing cluster colours as strings.
#' @noRd

.hc_update_clustering_algorithm_driver <- function(new_algo = NULL, gtc = NULL) {
  # recycle clustering of all algos from alluvial function, if it exists, otherwise conduct clustering now:
  if (!"alluvials" %in% names(hcobject[["integrated_output"]])) {
    all_clusterings <- run_all_cluster_algos()
  } else {
    all_clusterings <- hcobject[["integrated_output"]][["alluvials"]]
  }

  if (base::is.null(gtc)) {
    gtc <- all_clusterings[[new_algo]]
  } else {
    gtc <- gtc
    base::colnames(gtc) <- base::c("gene", "cluster")
    new_algo <- "my_clusters"
  }

  gtc[] <- base::lapply(gtc, base::as.character)
  new_cluster_info <- NULL

  for (c in base::unique(gtc$cluster)) {
    tmp <- gtc[gtc$cluster == c, ]
    gene_no <- base::nrow(tmp)
    gene_n <- base::paste0(tmp$gene, collapse = ",")
    if (c == "white") {
      cluster_included <- "no"
      vertexsize <- 1
    } else {
      cluster_included <- "yes"
      vertexsize <- 3
    }

    color <- c
    conditions <- base::paste0(
      .hc_gfc_condition_names(hcobject[["integrated_output"]][["GFC_all_layers"]]),
      collapse = "#"
    )
    gfc_means <- .hc_gfc_colmeans_for_genes(
      hcobject[["integrated_output"]][["GFC_all_layers"]],
      genes = tmp$gene
    )

    grp_means <- base::paste0(base::round(gfc_means, 3), collapse = ",")

    new_cluster_info <- rbind(
      new_cluster_info,
      data.frame(
        clusters = c,
        gene_no = gene_no,
        gene_n = gene_n,
        cluster_included = cluster_included,
        color = c,
        conditions = conditions,
        grp_means = grp_means,
        vertexsize = vertexsize,
        stringsAsFactors = FALSE
      )
    )
  }

  .hc_set_bridge_hcobject_slot(c("global_settings", "chosen_clustering_algo"), new_algo)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "cluster_information"), new_cluster_info)
  .hc_plot_cluster_heatmap_driver(return_HM = TRUE)
  .hc_plot_integrated_network_driver(layout = hcobject[["integrated_output"]][["cluster_calc"]][["layout"]])
}


