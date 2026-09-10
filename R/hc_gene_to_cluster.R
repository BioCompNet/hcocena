#' Gene To Cluster Dictionary
#'
#' The function maps the gene names to their corresponding cluster.
#' @param cluster_information Cluster table, typically
#'   `hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]`.
#' @return A data frame with two columns, the first containing gene names as strings, the second containing cluster colours as strings.
#' @noRd

.hc_gene_to_cluster_impl <- function(cluster_information = hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]) {
  gtc <- base::do.call(rbind, base::apply(cluster_information, 1, function(x) {
    tmp <- x["gene_n"] %>%
      base::strsplit(., split = ",") %>%
      base::unlist(.)
    base::data.frame(gene = tmp, color = base::rep(x["color"], base::length(tmp)))
  }))
  return(gtc)
}

#' Gene-to-module table (S4 API)
#'
#' Returns the mapping of every network gene to the module it was assigned to.
#'
#' @param hc A `HCoCenaExperiment`.
#' @return A data frame with two columns: `gene` (gene symbol) and `color`
#'   (module colour). Genes that were not assigned to any module carry the
#'   colour `"white"`.
#' @export
hc_gene_to_cluster <- function(hc) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment`.")
  }
  cluster_info <- as.list(hc@integration@cluster)[["cluster_information"]]
  if (base::is.null(cluster_info) || base::nrow(cluster_info) == 0) {
    stop("No cluster information found. Run `hc_cluster_calculation()` first.")
  }
  .hc_gene_to_cluster_impl(base::as.data.frame(cluster_info, stringsAsFactors = FALSE))
}
