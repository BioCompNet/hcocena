# Conversion helpers between legacy `hcobject` and S4 `HCoCenaExperiment`.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.hc_to_data_frame <- function(x) {
  if (is.null(x) || base::length(x) == 0) {
    return(S4Vectors::DataFrame())
  }
  if (inherits(x, "DataFrame")) {
    return(x)
  }
  if (base::is.data.frame(x)) {
    return(S4Vectors::DataFrame(x, check.names = FALSE))
  }
  S4Vectors::DataFrame(as.list(x), check.names = FALSE)
}

.hc_rows_to_data_frame <- function(rows) {
  if (base::length(rows) == 0) {
    return(S4Vectors::DataFrame())
  }

  all_names <- base::unique(base::unlist(base::lapply(rows, base::names)))
  cols <- stats::setNames(vector("list", base::length(all_names)), all_names)

  for (nm in all_names) {
    vals <- base::lapply(rows, function(r) {
      if (nm %in% base::names(r)) r[[nm]] else NA
    })
    scalar <- base::all(base::vapply(vals, function(v) {
      base::length(v) <= 1 && !base::is.list(v)
    }, logical(1)))
    cols[[nm]] <- if (scalar) base::unlist(vals) else I(vals)
  }

  S4Vectors::DataFrame(cols, check.names = FALSE)
}

.hc_build_layer_config <- function(hcobject) {
  layers <- hcobject[["layers"]] %||% list()
  if (base::length(layers) == 0) {
    return(S4Vectors::DataFrame())
  }

  layer_ids <- base::names(layers)
  if (is.null(layer_ids)) {
    layer_ids <- base::paste0("set", base::seq_along(layers))
  }

  layer_names <- hcobject[["layers_names"]] %||% layer_ids
  layer_settings <- hcobject[["layer_settings"]] %||% list()
  cutoffs <- hcobject[["cutoff_vec"]] %||% NULL

  rows <- vector("list", base::length(layers))
  for (i in base::seq_along(layers)) {
    lid <- layer_ids[[i]]
    pair <- layers[[i]]
    set_cfg <- layer_settings[[lid]]
    if (is.null(set_cfg) && base::length(layer_settings) >= i) {
      set_cfg <- layer_settings[[i]]
    }

    row <- list(
      layer_id = lid,
      layer_name = as.character(layer_names[[i]]),
      count_source = if (base::length(pair) >= 1) as.character(pair[[1]]) else NA_character_,
      annotation_source = if (base::length(pair) >= 2) as.character(pair[[2]]) else NA_character_
    )

    if (!is.null(set_cfg) && base::length(set_cfg) > 0) {
      for (nm in base::names(set_cfg)) {
        row[[nm]] <- set_cfg[[nm]]
      }
    }

    if (!is.null(cutoffs) && base::length(cutoffs) >= i) {
      row[["cutoff"]] <- cutoffs[[i]]
    }

    rows[[i]] <- row
  }

  .hc_rows_to_data_frame(rows)
}

.hc_build_mae <- function(hcobject, layer_config) {
  data <- hcobject[["data"]] %||% list()
  if (base::length(data) == 0) {
    exps <- S4Vectors::SimpleList()
    base::names(exps) <- character(0)
    return(MultiAssayExperiment::MultiAssayExperiment(
      experiments = exps
    ))
  }

  if (base::nrow(layer_config) > 0 && "layer_id" %in% base::colnames(layer_config)) {
    layer_ids <- as.character(layer_config$layer_id)
  } else {
    layer_ids <- base::gsub("_counts$", "", base::grep("_counts$", base::names(data), value = TRUE))
  }

  exps <- S4Vectors::SimpleList()
  for (lid in layer_ids) {
    ckey <- base::paste0(lid, "_counts")
    akey <- base::paste0(lid, "_anno")
    if (!(ckey %in% base::names(data) && akey %in% base::names(data))) {
      next
    }

    counts <- data[[ckey]]
    anno <- data[[akey]]
    if (is.null(counts) || is.null(anno)) {
      next
    }

    if (base::is.data.frame(counts)) {
      counts <- base::as.matrix(counts)
    }
    if (!base::is.matrix(counts)) {
      next
    }

    anno_df <- S4Vectors::DataFrame(anno, check.names = FALSE)
    if (is.null(base::rownames(anno_df)) && !is.null(base::colnames(counts))) {
      base::rownames(anno_df) <- base::colnames(counts)
    }
    if (!is.null(base::colnames(counts)) && !is.null(base::rownames(anno_df))) {
      shared <- base::intersect(base::colnames(counts), base::rownames(anno_df))
      counts <- counts[, shared, drop = FALSE]
      anno_df <- anno_df[shared, , drop = FALSE]
    }

    exps[[lid]] <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts),
      colData = anno_df
    )
  }

  if (is.null(base::names(exps))) {
    base::names(exps) <- character(0)
  }
  MultiAssayExperiment::MultiAssayExperiment(experiments = exps)
}

.hc_build_layer_results <- function(layer_specific_outputs) {
  layer_specific_outputs <- layer_specific_outputs %||% list()
  if (base::length(layer_specific_outputs) == 0) {
    return(S4Vectors::SimpleList())
  }

  out <- S4Vectors::SimpleList()
  for (nm in base::names(layer_specific_outputs)) {
    x <- layer_specific_outputs[[nm]]
    out[[nm]] <- new(
      "HCoCenaLayerResult",
      part1 = S4Vectors::SimpleList(x[["part1"]] %||% list()),
      part2 = S4Vectors::SimpleList(x[["part2"]] %||% list())
    )
  }
  out
}

.hc_reference_registry <- function(supplement) {
  supplement <- supplement %||% list()
  if (base::length(supplement) == 0) {
    return(S4Vectors::DataFrame())
  }

  nm <- base::names(supplement)
  if (is.null(nm)) {
    nm <- base::paste0("supplement_", base::seq_along(supplement))
  }

  rows <- vector("list", 0)
  for (i in base::seq_along(supplement)) {
    key <- nm[[i]]
    if (is.null(key) || !base::nzchar(key)) {
      key <- base::paste0("supplement_", i)
    }
    val <- supplement[[i]]

    if (is.null(val) || base::length(val) == 0) {
      next
    }
    if (base::is.list(val)) {
      val <- base::unlist(val, use.names = FALSE)
    }
    src <- as.character(val)
    src <- src[!base::is.na(src) & base::nzchar(src)]
    if (base::length(src) == 0) {
      next
    }

    rows[[base::length(rows) + 1]] <- base::data.frame(
      name = base::rep(key, base::length(src)),
      source = src,
      stringsAsFactors = FALSE
    )
  }

  if (base::length(rows) == 0) {
    return(S4Vectors::DataFrame())
  }
  out <- base::do.call(base::rbind, rows)
  base::rownames(out) <- NULL
  S4Vectors::DataFrame(out, check.names = FALSE)
}

.hc_default_object <- function() {
  list(
    working_directory = list(),
    data = list(),
    supplementary_data = list(),
    global_settings = list(),
    layer_settings = list(),
    layers = list(),
    supplement = list(),
    layers_names = list(),
    layer_specific_outputs = list(),
    integrated_output = list(),
    cutoff_vec = NULL,
    satellite_outputs = list()
  )
}

#' Convert legacy `hcobject` to `HCoCenaExperiment`
#'
#' @param hcobject Legacy list-style analysis object.
#' @return A valid `HCoCenaExperiment`.
#' @export
as_hcocena <- function(hcobject) {
  if (inherits(hcobject, "HCoCenaExperiment")) {
    return(hcobject)
  }
  if (!base::is.list(hcobject)) {
    stop("`hcobject` must be a list or `HCoCenaExperiment`.")
  }

  global_df <- .hc_to_data_frame(hcobject[["global_settings"]])
  paths_df <- .hc_to_data_frame(hcobject[["working_directory"]])
  layer_df <- .hc_build_layer_config(hcobject)

  cfg <- new("HCoCenaConfig", global = global_df, layer = layer_df, paths = paths_df)
  refs <- new(
    "HCoCenaReferences",
    registry = .hc_reference_registry(hcobject[["supplement"]]),
    data = S4Vectors::SimpleList(hcobject[["supplementary_data"]] %||% list())
  )
  layer_results <- .hc_build_layer_results(hcobject[["layer_specific_outputs"]])
  integration <- new(
    "HCoCenaIntegration",
    combined_edgelist = .hc_to_data_frame(hcobject[["integrated_output"]][["combined_edgelist"]]),
    graph = hcobject[["integrated_output"]][["merged_net"]],
    gfc = .hc_to_data_frame(hcobject[["integrated_output"]][["GFC_all_layers"]]),
    cluster = S4Vectors::SimpleList(hcobject[["integrated_output"]][["cluster_calc"]] %||% list())
  )

  hc <- new(
    "HCoCenaExperiment",
    mae = .hc_build_mae(hcobject, layer_df),
    config = cfg,
    references = refs,
    layer_results = layer_results,
    integration = integration,
    satellite = S4Vectors::SimpleList(hcobject[["satellite_outputs"]] %||% list()),
    provenance = S4Vectors::DataFrame(
      created_at = as.character(Sys.time()),
      source = "legacy_hcobject"
    )
  )
  methods::validObject(hc)
  hc
}

.hc_row_to_list <- function(df) {
  if (base::nrow(df) == 0 || base::ncol(df) == 0) {
    return(list())
  }
  out <- as.list(df[1, , drop = FALSE])
  base::lapply(out, function(v) {
    if (base::length(v) == 1) v[[1]] else v
  })
}

#' Convert `HCoCenaExperiment` to legacy `hcobject`
#'
#' @param hc A `HCoCenaExperiment` (or list, returned unchanged).
#' @return Legacy list-style `hcobject`.
#' @export
as_hcobject <- function(hc) {
  if (base::is.list(hc) && !inherits(hc, "HCoCenaExperiment")) {
    return(hc)
  }
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be `HCoCenaExperiment` or legacy list `hcobject`.")
  }

  out <- .hc_default_object()
  out[["working_directory"]] <- .hc_row_to_list(hc@config@paths)
  out[["global_settings"]] <- .hc_row_to_list(hc@config@global)

  if (base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer)) {
    for (i in base::seq_len(base::nrow(hc@config@layer))) {
      lid <- as.character(hc@config@layer$layer_id[[i]])
      count_src <- if ("count_source" %in% base::colnames(hc@config@layer)) as.character(hc@config@layer$count_source[[i]]) else NA_character_
      anno_src <- if ("annotation_source" %in% base::colnames(hc@config@layer)) as.character(hc@config@layer$annotation_source[[i]]) else NA_character_
      out[["layers"]][[lid]] <- c(count_src, anno_src)
    }
    if ("layer_name" %in% base::colnames(hc@config@layer)) {
      out[["layers_names"]] <- as.character(hc@config@layer$layer_name)
    } else {
      out[["layers_names"]] <- as.character(hc@config@layer$layer_id)
    }

    setting_cols <- base::setdiff(
      base::colnames(hc@config@layer),
      c("layer_id", "layer_name", "count_source", "annotation_source", "cutoff")
    )
    for (i in base::seq_len(base::nrow(hc@config@layer))) {
      lid <- as.character(hc@config@layer$layer_id[[i]])
      set_cfg <- list()
      for (nm in setting_cols) {
        set_cfg[[nm]] <- hc@config@layer[[nm]][[i]]
      }
      out[["layer_settings"]][[lid]] <- set_cfg
    }
    if ("cutoff" %in% base::colnames(hc@config@layer)) {
      out[["cutoff_vec"]] <- as.numeric(hc@config@layer$cutoff)
    }
  }

  exps <- MultiAssayExperiment::experiments(hc@mae)
  if (base::length(exps) > 0) {
    for (nm in base::names(exps)) {
      se <- exps[[nm]]
      out[["data"]][[base::paste0(nm, "_counts")]] <- SummarizedExperiment::assay(se, "counts")
      out[["data"]][[base::paste0(nm, "_anno")]] <- base::as.data.frame(SummarizedExperiment::colData(se))
    }
  }

  out[["supplementary_data"]] <- as.list(hc@references@data)
  if (base::nrow(hc@references@registry) > 0 &&
      all(c("name", "source") %in% base::colnames(hc@references@registry))) {
    out[["supplement"]] <- stats::setNames(
      as.list(as.character(hc@references@registry$source)),
      as.character(hc@references@registry$name)
    )
  }

  if (base::length(hc@layer_results) > 0) {
    for (nm in base::names(hc@layer_results)) {
      lr <- hc@layer_results[[nm]]
      if (inherits(lr, "HCoCenaLayerResult")) {
        out[["layer_specific_outputs"]][[nm]] <- list(
          part1 = as.list(lr@part1),
          part2 = as.list(lr@part2)
        )
      } else {
        out[["layer_specific_outputs"]][[nm]] <- as.list(lr)
      }
    }
  }

  out[["integrated_output"]][["combined_edgelist"]] <- base::as.data.frame(hc@integration@combined_edgelist)
  out[["integrated_output"]][["merged_net"]] <- hc@integration@graph
  out[["integrated_output"]][["GFC_all_layers"]] <- base::as.data.frame(hc@integration@gfc)
  out[["integrated_output"]][["cluster_calc"]] <- as.list(hc@integration@cluster)
  out[["satellite_outputs"]] <- as.list(hc@satellite)

  out
}
