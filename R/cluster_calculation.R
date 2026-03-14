#' Cluster Calculation
#' 
#' The function offers several community detection algorithms to identify dense regions in the co-expression network.
#'  These dense regions represent collections of highly co-expressed genes likely to form a functional group.
#' @param cluster_algo The clustering algorithm to be used. The choice is between "cluster_leiden" (default), "cluster_louvain", "cluster_label_prop", "cluster_fast_greedy",
#'  "cluster_infomap", "cluster_walktrap" and "auto" (in which case all are tested and the one with the highest modularity is chosen).
#' @param no_of_iterations Some of the algorithms are iterative (e.g. Leiden Algorithm). Set here, how many iterations should be performed. 
#'  For information on which other algorithms are iterative, please refer to their documentation in the igraph or leidenbase package. Default is 2.
#' @param max_cluster_count_per_gene The maximum number of different clusters a
#'  gene is allowed to be associated with during the different iterations
#'  before it is marked as indecisive and removed. Default is 1.
#' @param resolution The cluster resolution if the cluster algorithm is set to "cluster_leiden". Default is 0.1. Higher values result in more clusters and vice versa.
#' @param partition_type Name of the partition type. Select from 'CPMVertexPartition', 'ModularityVertexPartition', 'RBConfigurationVertexPartition' and 'RBERVertexPartition'. Default is 'RBConfigurationVertexPartition'.
#' @param return_result Logical. If `TRUE`, return the cluster table instead of
#'  storing it in `hcobject`.

#' @export 

cluster_calculation <- function(cluster_algo = "cluster_leiden",
                                no_of_iterations = 2,
                                resolution = 0.1,
                                partition_type = "RBConfigurationVertexPartition",
                                max_cluster_count_per_gene = 1,
                                return_result = FALSE) {
  .hc_legacy_warning("cluster_calculation")

  out <- .hc_run_modern_from_legacy_capture(
    .hc_cluster_calculation_impl,
    cluster_algo = cluster_algo,
    no_of_iterations = no_of_iterations,
    resolution = resolution,
    partition_type = partition_type,
    max_cluster_count_per_gene = max_cluster_count_per_gene,
    return_result = return_result
  )

  if (isTRUE(return_result) && !base::is.null(out$result)) {
    return(out$result)
  }

  invisible(NULL)
}
