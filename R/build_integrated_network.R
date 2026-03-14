#' Network integration
#' 
#' The previously constructed layer-specific networks are being integrated. 
#' @param mode A string, either "u", if network integration is to be done by union (default), or "i", if integration is to be done by intersection.
#'  For details please refer to the information pages provided in the repository's Wiki.
#' @param multi_edges One of "min", "mean" or "max" resulting in the simplification of the multigraph by using the minimum, the mean or the maximum edge weight among the multiple edges, respectively.
#'  Multiple edges occur when the edge was present in more than one datset.
#' @param GFC_when_missing The value to substitute missing data in the case where some genes were not measured in all but only some of the datasets.
#' @param with Either an integer giving the number of the dataset to be used as reference (e.g., 1) or the name given to the layer. 
#'  Can be ignored when integration is done by union.
#' @export

build_integrated_network <- function(mode = "u", 
                                      with = NULL,
                                      multi_edges = "min",        
                                      GFC_when_missing = -hcobject[["global_settings"]][["range_GFC"]]){
  .hc_legacy_warning("build_integrated_network")
  invisible(.hc_run_modern_from_legacy(
    hc_build_integrated_network,
    mode = mode,
    with = with,
    multi_edges = multi_edges,
    GFC_when_missing = GFC_when_missing
  ))
}



