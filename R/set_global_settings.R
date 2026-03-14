#' Define Global Settings
#' 
#' Receives all settings that are globally valid, i.e., that are not dataset specific.
#' @param organism Specification of the organism needed for all TF-related analyses. Set to "human" for human data and "mouse" for mouse data. Currently, hCoCena only supports human and mouse data, however, all analysis steps not related to TFs can be applied to other organisms.
#' @param control_keyword Either 'none' if no controls are present (only possible when analysing one single dataset) or a string contained in the control sample descriptor of all annotation files, e.g. "healthy".
#' 	The string must only be contained in the descriptor, e.g. "healthy" would work for "rhinovirusSetHealthy" and "influenzaSetHealthy", it does not have to match it perfectly.
#' @param variable_of_interest The name of the column that mus be rpesent in all annotation files and which will be used for grouping samples, e.g., "condition"".
#' @param min_nodes_number_for_network An integer. The minimum number of nodes in the subsequently created network that can define a graph component. Graph components with less nodes will be discarded. Default is 50.
#' @param min_nodes_number_for_cluster An integer. The minimum number of nodes that constitute a module/cluster when detecting community structures in the network. Default is 50.
#' @param range_GFC A float. Defines the maximum value the group fold changes (GFCs) can acquire, all values above this value or beneath its negative will be truncated. Default is 2.0.
#' @param layout_algorithm Layout algorithm used for the network. Supported values are:
#'  `"layout_with_stress"` (default), `"layout_with_sparse_stress"`, `"layout_with_fr"`,
#'  `"layout_with_drl"`, `"layout_with_kk"` and `"cytoscape"`.
#'  The stress-based layouts require the optional `graphlayouts` package and are typically
#'  more stable/readable than Fruchterman-Reingold for medium/large graphs.
#'  `"cytoscape"` uses an externally calculated Cytoscape layout.
#' @param data_in_log Boolean. Whether or not the provided gene expression data is logged to the base of 2.
#' @export



set_global_settings <- function(organism = "human",
								control_keyword = "none",
								variable_of_interest = "merged",
								min_nodes_number_for_network = 50,
								min_nodes_number_for_cluster = 50,
								range_GFC = 2.0,
								layout_algorithm = "layout_with_stress",
								data_in_log = TRUE){
	.hc_legacy_warning("set_global_settings")

	invisible(.hc_run_modern_from_legacy(
		.hc_set_global_settings_impl,
		organism = organism,
		control_keyword = control_keyword,
		variable_of_interest = variable_of_interest,
		min_nodes_number_for_network = min_nodes_number_for_network,
		min_nodes_number_for_cluster = min_nodes_number_for_cluster,
		range_GFC = range_GFC,
		layout_algorithm = layout_algorithm,
		data_in_log = data_in_log
	))
}

.hc_normalize_layout_algorithm <- function(layout_algorithm) {
  if (!is.character(layout_algorithm) || length(layout_algorithm) != 1 || is.na(layout_algorithm)) {
    stop("`layout_algorithm` must be a single character string.")
  }
  x <- tolower(trimws(layout_algorithm))
  aliases <- list(
    layout_with_fr = c("layout_with_fr", "fr", "fruchterman", "fruchterman_reingold", "fruchterman-rheingold"),
    layout_with_stress = c("layout_with_stress", "stress"),
    layout_with_sparse_stress = c("layout_with_sparse_stress", "sparse_stress", "sparsestress"),
    layout_with_drl = c("layout_with_drl", "drl"),
    layout_with_kk = c("layout_with_kk", "kk", "kamada_kawai", "kamada-kawai"),
    cytoscape = c("cytoscape")
  )
  for (canonical in names(aliases)) {
    if (x %in% aliases[[canonical]]) {
      return(canonical)
    }
  }
  stop(
    "`layout_algorithm` not supported. Use one of: ",
    paste(
      c(
        "\"layout_with_stress\"",
        "\"layout_with_sparse_stress\"",
        "\"layout_with_fr\"",
        "\"layout_with_drl\"",
        "\"layout_with_kk\"",
        "\"cytoscape\""
      ),
      collapse = ", "
    ),
    "."
  )
}
