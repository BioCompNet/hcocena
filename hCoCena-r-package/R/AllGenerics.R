#' Access `MultiAssayExperiment` from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @return A `MultiAssayExperiment::MultiAssayExperiment` object.
#' @export
setGeneric("hc_mae", function(x) standardGeneric("hc_mae"))

#' @rdname hc_mae
#' @export
setMethod("hc_mae", "HCoCenaExperiment", function(x) x@mae)

#' Access configuration from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @return A `HCoCenaConfig` object.
#' @export
setGeneric("hc_config", function(x) standardGeneric("hc_config"))

#' @rdname hc_config
#' @export
setMethod("hc_config", "HCoCenaExperiment", function(x) x@config)

#' Access layer results from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @return A `S4Vectors::SimpleList`.
#' @export
setGeneric("hc_layer_results", function(x) standardGeneric("hc_layer_results"))

#' @rdname hc_layer_results
#' @export
setMethod("hc_layer_results", "HCoCenaExperiment", function(x) x@layer_results)

#' Access integration payload from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @return A `HCoCenaIntegration` object.
#' @export
setGeneric("hc_integration", function(x) standardGeneric("hc_integration"))

#' @rdname hc_integration
#' @export
setMethod("hc_integration", "HCoCenaExperiment", function(x) x@integration)

#' Access cluster information from an `HCoCenaExperiment`
#' @param x An `HCoCenaExperiment`.
#' @return Cluster information or `NULL` if missing.
#' @export
setGeneric("hc_clusters", function(x) standardGeneric("hc_clusters"))

#' @rdname hc_clusters
#' @export
setMethod("hc_clusters", "HCoCenaExperiment", function(x) {
  if ("cluster_information" %in% base::names(x@integration@cluster)) {
    x@integration@cluster[["cluster_information"]]
  } else {
    NULL
  }
})

