#' Fail early and legibly when Cytoscape cannot be reached
#'
#' RCy3 talks to Cytoscape over CyREST. With Cytoscape closed it receives an
#' empty response and dies inside its own parsing with "$ operator is invalid
#' for atomic vectors", which gives the user nothing to act on.
#' @noRd
.hc_require_cytoscape <- function(what) {
  if (!base::requireNamespace("RCy3", quietly = TRUE)) {
    stop("Package `RCy3` is required for ", what,
         ". Install it with BiocManager::install(\"RCy3\").", call. = FALSE)
  }
  ok <- tryCatch(
    {
      RCy3::cytoscapePing()
      TRUE
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
  if (!isTRUE(ok)) {
    stop(
      "Cannot reach Cytoscape, which ", what, " needs. Start the Cytoscape ",
      "desktop application and leave it open, then try again. If it is ",
      "already running, check that the CyREST port is reachable with ",
      "`RCy3::cytoscapePing()`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Export Integrated Network to local R session
#' Due to a lacking possibility of communication with Cytoscape from within Docker container, you need to export all necessary information for import into a local R session.
#' @param file Path to the folder where the network information should be saved. Default path is set to the save folder defined in the global settings.
#' @noRd

.hc_export_to_local_folder_driver <- function(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]])) {

  network_df <- igraph::as_data_frame(network_filt())
  network_df$weight <- NULL

  readr::write_delim(network_df,
    file = paste0(file, "/network.txt"),
    delim = "\t",
    col_names = TRUE
  )
}


#' Export Integrated Network To Cytoscape
#'
#' Attention: Cytoscape Software must be open.
#' Due to difficulties in the communication between R/RCy3 and Cytoscape, you need to manually stop the function in R as soon as the table of nodes and edges appears in Cytoscape.(by pressing the little stop sign above the console).
#' @param name A string. The name given to the graph in Cytoscape. Default is "my igraph".
#' @param docker_container Deprecated legacy flag kept for backward compatibility.
#' @noRd

.hc_export_to_cytoscape_driver <- function(name = "my igraph", docker_container = FALSE) {
  .hc_require_cytoscape("`hc_export_to_cytoscape()`")

  RCy3::createNetworkFromIgraph(network_filt(), name)
}


#' Import Layout From Cytoscape
#'
#' Imports the layout of a network currently open in Cytoscape.
#' @noRd

.hc_import_layout_from_cytoscape_driver <- function() {
  .hc_require_cytoscape("`hc_import_layout_from_cytoscape()`")

  l <- RCy3::getNodePosition() %>% as.matrix()
  rn <- rownames(l)
  l <- apply(l, 2, as.numeric)
  rownames(l) <- rn

  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "layout"), l)
}


#' Import network layout from file to your Docker container
#' Imports the layout generated within a local R session and Cytoscape back into the Docker container
#' @param file Exact path and file name containing the Cytoscape layout. Default path is set to the save folder defined in the global settings.
#' @noRd

.hc_import_layout_from_local_folder_driver <- function(file = base::paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/network_layout.csv")) {

  if (!base::file.exists(file)) {
    stop(
      "No layout file at `", file, "`. `hc_import_layout_from_local_folder()` ",
      "reads the layout that `hc_import_layout_from_cytoscape()` wrote in a ",
      "local R session; point `file` at that CSV, and keep its name unchanged.",
      call. = FALSE
    )
  }
  l <- utils::read.csv(
    file = file,
    row.names = 1
  )

  rn <- rownames(l)
  l <- apply(l, 2, as.numeric)
  rownames(l) <- rn

  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "layout"), l)
}





