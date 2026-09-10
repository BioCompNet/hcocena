#' Get Cluster Scores
#'
#' For every gene, the ratio of its edges to genes in the same cluster to its total number of edges is determined.
#' 	The corresponding values are returned as a data frame and a box plot is generated showing the scores for each of the clusters.
#' @param save A Boolean. Whether or not to save the plot to PDF (default is TRUE).

.hc_get_module_scores_driver <- function(save = TRUE) {
  gtc <- .hc_gene_to_cluster_impl()
  # remove white cluster since not of interest:
  gtc <- dplyr::filter(gtc, !color == "white")

  e <- hcobject[["integrated_output"]][["combined_edgelist"]]

  # The score is the fraction of a gene's neighbours that sit in its own module.
  # Two things used to go wrong here:
  #  * the numerator counted *unique neighbour genes* while the denominator
  #    counted *edges*. The integrated network is a multigraph (one edge per
  #    layer), so a gene whose neighbours are all in its own module scored
  #    1/n_layers instead of 1. Both sides now count unique neighbours.
  #  * the neighbour lookup rescanned the full edge list once per gene, making
  #    this O(n_genes * n_edges). It is now a single vectorised join.
  v1 <- base::as.character(e[["V1"]])
  v2 <- base::as.character(e[["V2"]])

  # undirected adjacency, both directions, de-duplicated across layers
  adj <- base::unique(base::data.frame(
    gene = base::c(v1, v2),
    neighbour = base::c(v2, v1),
    stringsAsFactors = FALSE
  ))
  adj <- adj[adj$gene != adj$neighbour, , drop = FALSE]

  colour_of <- stats::setNames(base::as.character(gtc$color), base::as.character(gtc$gene))
  adj$gene_colour <- colour_of[adj$gene]
  adj$neighbour_colour <- colour_of[adj$neighbour]

  n_neighbours <- base::table(adj$gene)
  n_same <- base::table(adj$gene[
    !base::is.na(adj$gene_colour) &
      !base::is.na(adj$neighbour_colour) &
      adj$gene_colour == adj$neighbour_colour
  ])

  total_vec <- base::as.numeric(n_neighbours[base::as.character(gtc$gene)])
  same_vec <- base::as.numeric(n_same[base::as.character(gtc$gene)])
  same_vec[base::is.na(same_vec)] <- 0
  gtc$score <- base::ifelse(base::is.na(total_vec) | total_vec == 0,
    NA_real_, same_vec / total_vec
  )

  module_sizes <- base::table(gtc$color)
  gtc$label <- base::paste0(
    gtc$color, " [", base::as.integer(module_sizes[base::as.character(gtc$color)]), "]"
  )


  p <- ggplot2::ggplot(gtc, ggplot2::aes(x = label, y = score, color = color, fill = color)) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_color_manual(values = base::sort(base::unique(gtc$color))) +
    ggplot2::scale_fill_manual(values = grDevices::adjustcolor(base::sort(base::unique(gtc$color)), alpha.f = 0.5)) +
    ggplot2::xlab("modules") +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::ggtitle("Module scores") +
    ggplot2::theme(legend.position = "none")

  graphics::plot(p)

  if (save) {
    .hc_export_ggplot_file(
      file = .hc_output_file("Module_scores.pdf"),
      plot = p,
      width = 8,
      height = 7
    )
  }

  gtc$label <- NULL
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "module_scores"), list(scores_per_gene = gtc, plot = p))
}


