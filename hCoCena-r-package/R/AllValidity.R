# Validity functions for hCoCena S4 classes.

.hc_unique_genes_from_mae <- function(mae) {
  exps <- MultiAssayExperiment::experiments(mae)
  if (base::length(exps) == 0) {
    return(character())
  }

  base::unique(base::unlist(base::lapply(exps, function(se) {
    base::rownames(SummarizedExperiment::assay(se, "counts"))
  })))
}

setValidity("HCoCenaLayerResult", function(object) {
  if (!inherits(object@part1, "SimpleList") || !inherits(object@part2, "SimpleList")) {
    return("`part1` and `part2` must be `S4Vectors::SimpleList` objects.")
  }
  TRUE
})

setValidity("HCoCenaExperiment", function(object) {
  errors <- character()

  exp_names <- base::names(MultiAssayExperiment::experiments(object@mae))
  layer_names <- base::names(object@layer_results)
  cfg_layers <- character()

  if (base::nrow(object@config@layer) > 0 &&
      "layer_id" %in% base::colnames(object@config@layer)) {
    cfg_layers <- as.character(object@config@layer$layer_id)
  }

  non_empty <- list(
    experiments = exp_names,
    layer_results = layer_names,
    config_layer = cfg_layers
  )
  non_empty <- non_empty[base::vapply(non_empty, base::length, integer(1)) > 0]

  if (base::length(non_empty) > 1) {
    ref <- non_empty[[1]]
    for (i in 2:base::length(non_empty)) {
      if (!base::setequal(ref, non_empty[[i]])) {
        errors <- c(
          errors,
          "Layer names in `mae`, `layer_results` and `config@layer$layer_id` must match."
        )
        break
      }
    }
  }

  voi <- NULL
  if (base::nrow(object@config@global) > 0 &&
      "voi" %in% base::colnames(object@config@global)) {
    voi <- as.character(object@config@global$voi[[1]])
  }

  if (!is.null(voi) && base::nzchar(voi) && base::length(exp_names) > 0) {
    for (nm in exp_names) {
      cd <- SummarizedExperiment::colData(MultiAssayExperiment::experiments(object@mae)[[nm]])
      if (!(voi %in% base::colnames(cd))) {
        errors <- c(errors, paste0(
          "Missing variable_of_interest (`", voi, "`) in colData of layer `", nm, "`."
        ))
      }
    }
  }

  if (base::nrow(object@config@layer) > 0 &&
      "cutoff" %in% base::colnames(object@config@layer) &&
      "layer_id" %in% base::colnames(object@config@layer)) {
    if (base::length(object@config@layer$cutoff) !=
        base::length(object@config@layer$layer_id)) {
      errors <- c(errors, "`cutoff` length in `config@layer` must match number of layers.")
    }
  }

  if (!is.null(object@integration@graph) && inherits(object@integration@graph, "igraph")) {
    graph_nodes <- igraph::V(object@integration@graph)$name
    all_genes <- .hc_unique_genes_from_mae(object@mae)
    if (base::length(base::setdiff(graph_nodes, all_genes)) > 0) {
      errors <- c(errors, "All graph nodes must be contained in the MAE layer gene sets.")
    }

    if ("cluster_information" %in% base::names(object@integration@cluster)) {
      ci <- object@integration@cluster[["cluster_information"]]
      cluster_genes <- character()

      if ((base::is.data.frame(ci) || inherits(ci, "DataFrame")) &&
          "genes" %in% base::colnames(ci)) {
        cluster_genes <- as.character(ci$genes)
      } else if ((base::is.data.frame(ci) || inherits(ci, "DataFrame")) &&
                 "gene" %in% base::colnames(ci)) {
        cluster_genes <- as.character(ci$gene)
      } else if ((base::is.data.frame(ci) || inherits(ci, "DataFrame")) &&
                 "gene_n" %in% base::colnames(ci)) {
        cluster_genes <- base::unlist(base::strsplit(as.character(ci$gene_n), split = ","))
      }

      cluster_genes <- base::unique(cluster_genes[!base::is.na(cluster_genes)])
      if (base::length(base::setdiff(cluster_genes, graph_nodes)) > 0) {
        errors <- c(errors, "Cluster assignments must reference nodes in the integrated graph.")
      }
    }
  }

  if (base::length(errors) == 0) TRUE else errors
})

