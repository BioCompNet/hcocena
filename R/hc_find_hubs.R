#' Find Hub Genes
#'
#' Hub genes are determined per cluster using a combined ranking based on weighted degree centrality, weighted closeness centrality and weighted betweenness centrality.
#'  A table of hub genes per cluster is returned as an excel file and a heatmap of the hub genes expression values is plotted per cluster and per dataset.
#' @param top An integer. All genes are ranked based on their hub potential, this parameter defines the number of the top ranked genes to be considered a hub gene. Default value is 10.
#' @param save A Boolean. Whether or not a labelled hub-network per cluster and the expression heatmap are to be save to PDF. Default is FALSE.
#' @param tree_layout A Boolean. Whether or not to depict the network witht tree layout, implying a sort of hierarchical structure to the network.
#' @param TF_only Either FALSE (default, all genes in clsuter are considered for hub genes), or "all" (all genes from transcriptionfactor supplementary file are considered for hub genes),
#'  or any gene category listed in the last column of the provided transcriptionfactor supplementary file (only that subgroup condired for hub genes).
#' @param plot A Boolean. Wheather or not to plot the network (per cluster) with highlighted hub nodes. Default is FALSE.
#' @param clusters Either "all" (default) or a vector of cluster colours for which the hub detection should be performed.
#' @noRd

.hc_find_hubs_driver <- function(clusters = c("all"),
                      top = 10,
                      tree_layout = FALSE,
                      TF_only = FALSE,
                      save = FALSE,
                      plot = FALSE) {
  gtc <- .hc_gene_to_cluster_impl()
  if (clusters[1] == "all") {
    clusters <- base::unique(gtc$color[!gtc$color == "white"])
  }

  cluster_labels <- .hc_hub_display_labels(clusters)

  hubs <- base::lapply(clusters, function(x) {
    tmp <- hub_node_detection(
      cluster = x, top = top, save = save, tree_layout = tree_layout,
      TF_only = TF_only, plot = plot, label = cluster_labels[[x]]
    )
    return(tmp$hub_nodes)
  })


  hubs_df <- base::lapply(hubs, function(x) {
    if (base::is.null(x)) {
      x <- base::rep(" ", top)
    }
    if (base::length(x) < top) {
      x <- base::c(x, base::rep(" ", top - base::length(x)))
    }
    return(x)
  }) %>%
    rlist::list.cbind() %>%
    base::as.data.frame()

  base::colnames(hubs_df) <- base::unname(cluster_labels[clusters])

  hubs_df[base::is.na(hubs_df)] <- " "

  for (col in base::colnames(hubs_df)) {
    if (base::file.exists(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Hub_genes.xlsx"))) {
      wb <- openxlsx::loadWorkbook(base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Hub_genes.xlsx"))
      if (col %in% wb$sheet_names) {
        openxlsx::removeWorksheet(wb, col)
      }
      openxlsx::addWorksheet(wb, col)
      openxlsx::writeData(wb, sheet = col, dplyr::select(hubs_df, tidyselect::all_of(col)), colNames = TRUE)
      .hc_save_workbook_atomic(
        wb = wb,
        file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Hub_genes.xlsx"),
        overwrite = TRUE
      )
    } else {
      .hc_write_xlsx_atomic(dplyr::select(hubs_df, tidyselect::all_of(col)),
        file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/Hub_genes.xlsx"),
        sheetName = col,
        colNames = TRUE, rowNames = FALSE, append = FALSE, overwrite = TRUE
      )
    }

    if ("hub_out" %in% base::names(hcobject[["satellite_outputs"]])) {
      .hc_set_bridge_hcobject_slot(c("satellite_outputs", "hub_out", col), dplyr::select(hubs_df, tidyselect::all_of(col)))
    } else {
      .hc_set_bridge_hcobject_slot(c("satellite_outputs", "hub_out"), list())
      .hc_set_bridge_hcobject_slot(c("satellite_outputs", "hub_out", col), dplyr::select(hubs_df, tidyselect::all_of(col)))
    }
  }


  for (col in base::colnames(hubs_df)) {
    hub_genes <- base::trimws(base::as.character(dplyr::pull(hubs_df, col)))
    hub_genes <- base::unique(hub_genes[!base::is.na(hub_genes) & base::nzchar(hub_genes)])
    if (base::length(hub_genes) > 0) {
      .hc_visualize_gene_expression_driver(
        genes = hub_genes,
        name = base::paste0("Hub_genes_", col, "_module_expression"),
        width = 10,
        height = 10,
        save = save,
        label_map = cluster_labels
      )
    }
  }
}


#' Resolve module display labels for a vector of cluster colours
#'
#' Clusters are addressed internally by colour, but every user-facing output
#' uses the module label (`M1`, `M2.1`, ...) held in
#' `cluster_calc[["module_label_map"]]`. After [hc_split_modules()] the colours
#' of new submodules are generated hex codes, so labelling by colour produces
#' names such as `#30A89C` that cannot be matched to any other output.
#' Falls back to the colour itself when no label is available.
#' @noRd
.hc_hub_display_labels <- function(clusters) {
  clusters <- base::as.character(clusters)
  cluster_calc <- hcobject[["integrated_output"]][["cluster_calc"]]
  label_map <- if (!base::is.null(cluster_calc) &&
    "module_label_map" %in% base::names(cluster_calc)) {
    cluster_calc[["module_label_map"]]
  } else {
    NULL
  }
  label_map <- .hc_resolve_module_label_map_for_colors(
    label_map = label_map,
    module_colors = clusters
  )
  if (base::is.null(label_map) || base::length(label_map) == 0) {
    return(stats::setNames(clusters, clusters))
  }
  out <- base::vapply(clusters, function(cl) {
    lbl <- label_map[[cl]]
    if (base::is.null(lbl) || base::length(lbl) == 0 || base::is.na(lbl[[1]]) ||
      !base::nzchar(base::as.character(lbl[[1]]))) {
      cl
    } else {
      base::as.character(lbl[[1]])
    }
  }, FUN.VALUE = base::character(1))
  # Duplicate labels would collide as Excel sheet names and file names.
  if (base::anyDuplicated(out) > 0) {
    return(stats::setNames(clusters, clusters))
  }
  stats::setNames(out, clusters)
}
