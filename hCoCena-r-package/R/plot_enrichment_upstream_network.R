#' Plot a module knowledge network from enrichment + upstream inference
#'
#' Builds a three-column network:
#' modules -> enrichment terms and modules -> upstream terms (TF/Pathway).
#' The output is shown together with the hCoCena module heatmap to preserve the
#' direct connection to module-level expression patterns. The heatmap panel
#' automatically follows the last upstream inference settings (e.g. GFC vs FC).
#'
#' @param enrichment_mode Character scalar. One of `"selected"` (default) or
#'   `"significant"`.
#' @param upstream_mode Character scalar. One of `"selected"` (default) or
#'   `"significant"`.
#' @param clusters Either `"all"` (default) or a character vector of module
#'   colors to include.
#' @param max_enrichment_per_module Optional positive integer. If set, keeps at
#'   most this many enrichment edges per module (best q-values first).
#' @param max_upstream_per_module Optional positive integer. If set, keeps at
#'   most this many upstream edges per module (best q-values first).
#' @param label_mode Character scalar controlling term label density:
#'   `"both"` (default), `"upstream_only"`, or `"focus_only"`.
#' @param show_plot Logical; if `TRUE` (default), prints the combined overview
#'   (heatmap + network) in the active graphics device.
#' @param save_pdf Logical; if `TRUE` (default), writes a multi-page PDF to the
#'   current hCoCena save folder (overview + per-module focus pages).
#' @param pdf_name Output PDF filename.
#' @param gfc_scale_limits Optional numeric vector controlling the left module
#'   heatmap color scale limits. Provide one positive number (`x` -> `c(-x, x)`)
#'   or two numbers (`c(min, max)`). If NULL, uses upstream inference settings
#'   first (if available), then main heatmap settings, then `c(-range_GFC, range_GFC)`.
#' @param pdf_width Optional numeric width (inches) for network PDFs.
#'   If NULL (default), width is auto-estimated from content.
#' @param pdf_height Optional numeric height (inches) for network PDFs.
#'   If NULL (default), height is auto-estimated from content.
#' @param pdf_pointsize Numeric base pointsize used for PDF export devices.
#'   Default is 11.
#' @param overall_plot_scale Numeric scaling factor for text and line sizes.
#'
#' @return A list with overview/focus plots, nodes, edges, and output file.
#' @export
plot_enrichment_upstream_network <- function(enrichment_mode = "selected",
                                             upstream_mode = "selected",
                                             clusters = c("all"),
                                             max_enrichment_per_module = NULL,
                                             max_upstream_per_module = NULL,
                                             label_mode = "both",
                                             show_plot = TRUE,
                                             save_pdf = TRUE,
                                             pdf_name = "Module_Knowledge_Network.pdf",
                                             gfc_scale_limits = NULL,
                                             pdf_width = NULL,
                                             pdf_height = NULL,
                                             pdf_pointsize = 11,
                                             overall_plot_scale = 1) {
  .hc_legacy_warning("plot_enrichment_upstream_network")

  enrichment_mode <- base::match.arg(
    base::tolower(base::as.character(enrichment_mode)),
    choices = c("selected", "significant")
  )
  upstream_mode <- base::match.arg(
    base::tolower(base::as.character(upstream_mode)),
    choices = c("selected", "significant")
  )
  label_mode <- base::match.arg(
    base::tolower(base::as.character(label_mode)),
    choices = c("both", "upstream_only", "focus_only")
  )
  if (!base::is.logical(show_plot) || base::length(show_plot) != 1) {
    stop("`show_plot` must be TRUE or FALSE.")
  }
  if (!base::is.logical(save_pdf) || base::length(save_pdf) != 1) {
    stop("`save_pdf` must be TRUE or FALSE.")
  }
  if (!base::is.character(pdf_name) || base::length(pdf_name) != 1 || base::is.na(pdf_name) || pdf_name == "") {
    stop("`pdf_name` must be a non-empty filename.")
  }
  if (!base::is.null(pdf_width) &&
      (!base::is.numeric(pdf_width) || base::length(pdf_width) != 1 || base::is.na(pdf_width) || pdf_width <= 0)) {
    stop("`pdf_width` must be NULL or a single positive number.")
  }
  if (!base::is.null(pdf_height) &&
      (!base::is.numeric(pdf_height) || base::length(pdf_height) != 1 || base::is.na(pdf_height) || pdf_height <= 0)) {
    stop("`pdf_height` must be NULL or a single positive number.")
  }
  if (!base::is.numeric(pdf_pointsize) ||
      base::length(pdf_pointsize) != 1 ||
      base::is.na(pdf_pointsize) ||
      pdf_pointsize <= 0) {
    stop("`pdf_pointsize` must be a single positive number.")
  }
  if (!base::is.numeric(overall_plot_scale) ||
      base::length(overall_plot_scale) != 1 ||
      base::is.na(overall_plot_scale) ||
      overall_plot_scale <= 0) {
    stop("`overall_plot_scale` must be a positive numeric scalar.")
  }
  overall_plot_scale <- base::max(0.5, base::min(3, overall_plot_scale))

  normalize_scale_limits <- function(x) {
    if (base::is.null(x)) {
      return(NULL)
    }
    x <- suppressWarnings(base::as.numeric(x))
    if (base::length(x) == 1) {
      if (!base::is.finite(x) || x <= 0) {
        stop("`gfc_scale_limits` as single value must be finite and > 0.")
      }
      return(c(-base::abs(x), base::abs(x)))
    }
    if (base::length(x) != 2 || any(!base::is.finite(x))) {
      stop("`gfc_scale_limits` must be NULL, one positive number, or a numeric vector of length 2.")
    }
    x <- base::sort(x)
    if (x[1] == x[2]) {
      stop("`gfc_scale_limits` must have different min/max values.")
    }
    x
  }
  compute_scale_ticks <- function(lims) {
    br <- pretty(lims, n = 5)
    br <- br[
      br >= lims[1] - .Machine$double.eps^0.5 &
        br <= lims[2] + .Machine$double.eps^0.5
    ]
    if (!any(base::abs(br) < .Machine$double.eps^0.5)) {
      br <- base::sort(base::unique(base::c(br, 0)))
    }
    if (base::length(br) < 3) {
      br <- base::seq(lims[1], lims[2], length.out = 5)
    }
    list(
      breaks = br,
      labels = base::formatC(br, format = "fg", digits = 3)
    )
  }

  check_optional_limit <- function(x, arg) {
    if (base::is.null(x)) {
      return(NULL)
    }
    if (!base::is.numeric(x) || base::length(x) != 1 || base::is.na(x) || x < 1) {
      stop("`", arg, "` must be NULL or a positive integer.")
    }
    base::as.integer(x)
  }
  max_enrichment_per_module <- check_optional_limit(max_enrichment_per_module, "max_enrichment_per_module")
  max_upstream_per_module <- check_optional_limit(max_upstream_per_module, "max_upstream_per_module")

  file_prefix <- base::paste0(
    hcobject[["working_directory"]][["dir_output"]],
    hcobject[["global_settings"]][["save_folder"]]
  )

  cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]
  if (base::is.null(cluster_info) || base::nrow(cluster_info) == 0) {
    stop("No cluster information found. Run `cluster_calculation()` first.")
  }
  all_clusters <- base::unique(base::as.character(cluster_info$color))
  all_clusters <- all_clusters[all_clusters != "white" & !base::is.na(all_clusters)]
  if (base::length(all_clusters) == 0) {
    stop("No non-white modules available for plotting.")
  }
  if (clusters[1] == "all") {
    clusters <- all_clusters
  }
  clusters <- base::as.character(clusters)
  clusters <- clusters[clusters %in% all_clusters]
  if (base::length(clusters) == 0) {
    stop("No valid clusters selected.")
  }

  cluster_calc <- hcobject[["integrated_output"]][["cluster_calc"]]
  stored_hm <- cluster_calc[["heatmap_cluster"]]

  cluster_order <- all_clusters
  if (!base::is.null(stored_hm)) {
    row_order <- tryCatch(ComplexHeatmap::row_order(stored_hm), error = function(e) NULL)
    if (base::is.list(row_order) && base::length(row_order) > 0) {
      row_order <- row_order[[1]]
    }
    hm_rows <- tryCatch(base::rownames(stored_hm@ht_list[[1]]@matrix), error = function(e) NULL)
    if (!base::is.null(hm_rows)) {
      if (base::is.numeric(row_order) && base::length(row_order) == base::length(hm_rows)) {
        cluster_order <- hm_rows[row_order]
      } else if (base::is.character(row_order)) {
        cluster_order <- row_order
      }
      cluster_order <- cluster_order[cluster_order %in% all_clusters]
      cluster_order <- base::c(cluster_order, base::setdiff(all_clusters, cluster_order))
    }
  }
  cluster_order <- cluster_order[cluster_order %in% clusters]
  if (base::length(cluster_order) == 0) {
    stop("No modules remain after applying order and `clusters` filter.")
  }

  module_label_map <- cluster_calc[["module_label_map"]]
  if (!base::is.null(module_label_map) && base::length(module_label_map) > 0) {
    map_names <- base::names(cluster_calc[["module_label_map"]])
    module_label_map <- base::as.character(module_label_map)
    if (!base::is.null(map_names) && base::length(map_names) == base::length(module_label_map)) {
      base::names(module_label_map) <- base::as.character(map_names)
    }
    missing_before <- base::setdiff(cluster_order, base::names(module_label_map))
    if (base::length(missing_before) > 0) {
      inverse_map <- stats::setNames(base::names(module_label_map), base::as.character(module_label_map))
      if (base::all(cluster_order %in% base::names(inverse_map))) {
        module_label_map <- inverse_map
      }
    }
  }
  module_prefix <- cluster_calc[["module_prefix"]]
  if (base::is.null(module_prefix) || !base::is.character(module_prefix) || base::length(module_prefix) != 1) {
    module_prefix <- "M"
  }
  if (base::is.null(module_label_map) || base::length(module_label_map) == 0) {
    module_label_map <- stats::setNames(
      base::paste0(module_prefix, base::seq_along(cluster_order)),
      cluster_order
    )
  }
  missing_map <- base::setdiff(cluster_order, base::names(module_label_map))
  if (base::length(missing_map) > 0) {
    module_label_map <- base::c(
      module_label_map,
      stats::setNames(
        base::paste0(module_prefix, base::seq.int(base::length(module_label_map) + 1, length.out = base::length(missing_map))),
        missing_map
      )
    )
  }
  module_label_map <- module_label_map[cluster_order]

  read_xlsx_sheet <- function(file, sheet) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      return(NULL)
    }
    if (!base::file.exists(file)) {
      return(NULL)
    }
    out <- tryCatch(
      openxlsx::read.xlsx(xlsxFile = file, sheet = sheet),
      error = function(e) NULL
    )
    if (base::is.null(out) || !base::is.data.frame(out) || base::nrow(out) == 0) {
      return(NULL)
    }
    out
  }

  normalize_enrichment_df <- function(df) {
    if (base::is.null(df) || !base::is.data.frame(df) || base::nrow(df) == 0) {
      return(NULL)
    }
    out <- base::as.data.frame(df, stringsAsFactors = FALSE)
    if (!("database" %in% base::colnames(out))) {
      if ("DB" %in% base::colnames(out)) {
        out$database <- base::as.character(out$DB)
      } else {
        out$database <- NA_character_
      }
    }
    if (!("cluster" %in% base::colnames(out))) {
      return(NULL)
    }
    if (!("term" %in% base::colnames(out))) {
      if ("Description" %in% base::colnames(out)) {
        out$term <- base::as.character(out$Description)
      } else {
        return(NULL)
      }
    }
    if (!("qvalue" %in% base::colnames(out))) {
      if ("p.adjust" %in% base::colnames(out)) {
        out$qvalue <- suppressWarnings(base::as.numeric(out$p.adjust))
      } else if ("pvalue" %in% base::colnames(out)) {
        out$qvalue <- suppressWarnings(base::as.numeric(out$pvalue))
      } else {
        out$qvalue <- NA_real_
      }
    }
    out$database <- base::as.character(out$database)
    out$cluster <- base::as.character(out$cluster)
    out$term <- base::as.character(out$term)
    out$qvalue <- suppressWarnings(base::as.numeric(out$qvalue))
    out <- out[!base::is.na(out$cluster) & !base::is.na(out$term) & out$term != "", , drop = FALSE]
    if (base::nrow(out) == 0) {
      return(NULL)
    }
    out
  }

  normalize_upstream_df <- function(df) {
    if (base::is.null(df) || !base::is.data.frame(df) || base::nrow(df) == 0) {
      return(NULL)
    }
    out <- base::as.data.frame(df, stringsAsFactors = FALSE)
    if (!("resource" %in% base::colnames(out))) {
      if ("Resource" %in% base::colnames(out)) {
        out$resource <- base::as.character(out$Resource)
      } else {
        out$resource <- NA_character_
      }
    }
    if (!("cluster" %in% base::colnames(out))) {
      return(NULL)
    }
    if (!("term" %in% base::colnames(out))) {
      if ("source" %in% base::colnames(out)) {
        out$term <- base::as.character(out$source)
      } else {
        return(NULL)
      }
    }
    if (!("qvalue" %in% base::colnames(out))) {
      if ("p.adjust" %in% base::colnames(out)) {
        out$qvalue <- suppressWarnings(base::as.numeric(out$p.adjust))
      } else if ("pvalue" %in% base::colnames(out)) {
        out$qvalue <- suppressWarnings(base::as.numeric(out$pvalue))
      } else {
        out$qvalue <- NA_real_
      }
    }
    if (!("direction" %in% base::colnames(out))) {
      out$direction <- NA_character_
    }
    out$resource <- base::as.character(out$resource)
    out$resource[base::tolower(out$resource) == "pathway"] <- "Pathway"
    out$resource[base::tolower(out$resource) == "tf"] <- "TF"
    out$cluster <- base::as.character(out$cluster)
    out$term <- base::as.character(out$term)
    out$qvalue <- suppressWarnings(base::as.numeric(out$qvalue))
    out$direction <- base::as.character(out$direction)
    out <- out[!base::is.na(out$cluster) & !base::is.na(out$term) & out$term != "", , drop = FALSE]
    if (base::nrow(out) == 0) {
      return(NULL)
    }
    out
  }

  resolve_enrichment_df <- function(store, mode) {
    key <- if (identical(mode, "selected")) "selected_enrichments_all_dbs" else "significant_enrichments_all_dbs"
    direct <- normalize_enrichment_df(store[[key]])
    if (!base::is.null(direct) && base::nrow(direct) > 0) {
      return(direct)
    }
    if (identical(mode, "selected") && base::is.list(store[["top_all_dbs"]])) {
      cand <- normalize_enrichment_df(store[["top_all_dbs"]][["result"]])
      if (!base::is.null(cand) && base::nrow(cand) > 0) {
        return(cand)
      }
    }
    per_db_keys <- base::grep("^top_", base::names(store), value = TRUE)
    per_db_keys <- per_db_keys[!per_db_keys %in% c("top_all_dbs", "top_all_dbs_mixed")]
    if (base::length(per_db_keys) > 0) {
      collected <- base::lapply(per_db_keys, function(k) {
        entry <- store[[k]]
        if (!base::is.list(entry)) {
          return(NULL)
        }
        subkey <- if (identical(mode, "selected")) "selected_enrichments" else "significant_enrichments"
        normalize_enrichment_df(entry[[subkey]])
      })
      collected <- collected[!base::vapply(collected, base::is.null, FUN.VALUE = base::logical(1))]
      if (base::length(collected) > 0) {
        out <- base::do.call(base::rbind, collected)
        base::rownames(out) <- NULL
        return(out)
      }
    }
    NULL
  }

  resolve_upstream_df <- function(store, mode) {
    key <- if (identical(mode, "selected")) "selected_upstream_all" else "significant_upstream_all"
    normalize_upstream_df(store[[key]])
  }

  enrich_store <- hcobject[["integrated_output"]][["enrichments"]]
  if (base::is.null(enrich_store)) {
    enrich_store <- hcobject[["satellite_outputs"]][["enrichments"]]
  }
  if (base::is.null(enrich_store)) {
    enrich_store <- list()
  }
  enrich_key <- if (identical(enrichment_mode, "selected")) "selected_enrichments_all_dbs" else "significant_enrichments_all_dbs"
  enrich_df <- resolve_enrichment_df(enrich_store, enrichment_mode)
  if (base::is.null(enrich_df) || base::nrow(enrich_df) == 0) {
    enrich_file <- base::paste0(file_prefix, "/Enrichment_Selected_All_DBs.xlsx")
    enrich_df <- normalize_enrichment_df(read_xlsx_sheet(enrich_file, enrich_key))
  }
  if (base::is.null(enrich_df) || base::nrow(enrich_df) == 0) {
    stop("Missing enrichment results. Run `functional_enrichment()` first.")
  }

  upstream_store <- hcobject[["integrated_output"]][["upstream_inference"]]
  if (base::is.null(upstream_store)) {
    upstream_store <- hcobject[["satellite_outputs"]][["upstream_inference"]]
  }
  if (base::is.null(upstream_store)) {
    upstream_store <- list()
  }
  upstream_key <- if (identical(upstream_mode, "selected")) "selected_upstream_all" else "significant_upstream_all"
  upstream_df <- resolve_upstream_df(upstream_store, upstream_mode)
  if (base::is.null(upstream_df) || base::nrow(upstream_df) == 0) {
    upstream_file <- base::paste0(file_prefix, "/Upstream_Inference.xlsx")
    upstream_df <- normalize_upstream_df(read_xlsx_sheet(upstream_file, upstream_key))
  }
  if (base::is.null(upstream_df) || base::nrow(upstream_df) == 0) {
    stop("Missing upstream results. Run `upstream_inference()` first.")
  }
  upstream_settings <- upstream_store[["settings"]]
  upstream_activity_input <- "gfc"
  if (base::is.list(upstream_settings) &&
      "activity_input" %in% base::names(upstream_settings) &&
      base::length(upstream_settings$activity_input) > 0 &&
      !base::is.na(upstream_settings$activity_input[[1]])) {
    upstream_activity_input <- base::tolower(base::as.character(upstream_settings$activity_input[[1]]))
  }
  if (!upstream_activity_input %in% c("gfc", "fc", "expression")) {
    upstream_activity_input <- "gfc"
  }
  upstream_module_heatmap_mat <- upstream_store[["module_heatmap_matrix"]]
  if (!base::is.null(upstream_module_heatmap_mat)) {
    upstream_module_heatmap_mat <- upstream_module_heatmap_mat %>% base::as.matrix()
    if (base::nrow(upstream_module_heatmap_mat) == 0 || base::ncol(upstream_module_heatmap_mat) == 0) {
      upstream_module_heatmap_mat <- NULL
    }
  }
  upstream_module_heatmap_col_order <- upstream_store[["module_heatmap_col_order"]]
  if (base::is.null(upstream_module_heatmap_col_order)) {
    upstream_module_heatmap_col_order <- base::character(0)
  } else {
    upstream_module_heatmap_col_order <- base::as.character(upstream_module_heatmap_col_order)
  }
  upstream_module_heatmap_name <- upstream_store[["module_heatmap_name"]]
  if (base::is.null(upstream_module_heatmap_name) || base::length(upstream_module_heatmap_name) == 0 || base::is.na(upstream_module_heatmap_name[[1]])) {
    upstream_module_heatmap_name <- if (identical(upstream_activity_input, "fc")) "FC" else "GFC"
  }
  upstream_module_heatmap_name <- base::as.character(upstream_module_heatmap_name[[1]])
  resolved_gfc_scale_limits <- normalize_scale_limits(gfc_scale_limits)
  if (base::is.null(resolved_gfc_scale_limits) &&
      base::is.list(upstream_settings) &&
      "gfc_scale_limits" %in% base::names(upstream_settings)) {
    resolved_gfc_scale_limits <- tryCatch(
      normalize_scale_limits(upstream_settings[["gfc_scale_limits"]]),
      error = function(e) NULL
    )
  }
  if (base::is.null(resolved_gfc_scale_limits)) {
    resolved_gfc_scale_limits <- tryCatch(
      normalize_scale_limits(cluster_calc[["gfc_scale_limits"]]),
      error = function(e) NULL
    )
  }
  if (base::is.null(resolved_gfc_scale_limits)) {
    fallback_lim <- suppressWarnings(base::as.numeric(hcobject[["global_settings"]][["range_GFC"]]))
    if (!base::is.finite(fallback_lim) || fallback_lim <= 0) {
      fallback_lim <- 2
    }
    resolved_gfc_scale_limits <- c(-base::abs(fallback_lim), base::abs(fallback_lim))
  }
  resolved_gfc_scale_ticks <- compute_scale_ticks(resolved_gfc_scale_limits)

  build_fc_module_heatmap_from_settings <- function(fc_comparisons, cluster_order, cluster_info) {
    if (base::is.null(fc_comparisons) || base::length(fc_comparisons) == 0) {
      stop(
        "Upstream was run with `activity_input = 'fc'`, but no `fc_comparisons` are stored. ",
        "Please re-run `hc_upstream_inference(..., activity_input = 'fc', fc_comparisons = ...)`."
      )
    }
    merged_net <- hcobject[["integrated_output"]][["merged_net"]]
    if (base::is.null(merged_net)) {
      stop("Missing integrated network for FC heatmap reconstruction.")
    }
    net_genes <- igraph::get.vertex.attribute(merged_net, "name")
    net_genes <- base::unique(base::as.character(net_genes))
    net_genes <- net_genes[!base::is.na(net_genes) & net_genes != ""]
    if (base::length(net_genes) == 0) {
      stop("Could not extract genes from integrated network for FC heatmap reconstruction.")
    }
    parse_fc <- function(x) {
      x <- base::trimws(base::as.character(x))
      parts <- base::strsplit(x, "\\s*_vs_\\s*", perl = TRUE)[[1]]
      if (base::length(parts) != 2) {
        parts <- base::strsplit(x, "\\s+vs\\s+", perl = TRUE)[[1]]
      }
      if (base::length(parts) != 2) {
        stop("Invalid FC comparison format: '", x, "'. Use 'groupA_vs_groupB'.")
      }
      base::data.frame(
        comparison = base::paste0(base::trimws(parts[[1]]), "_vs_", base::trimws(parts[[2]])),
        numerator = base::trimws(parts[[1]]),
        denominator = base::trimws(parts[[2]]),
        stringsAsFactors = FALSE
      )
    }
    comparisons <- base::lapply(fc_comparisons, parse_fc)
    comparisons <- base::do.call(base::rbind, comparisons)
    comparisons <- comparisons[!duplicated(comparisons$comparison), , drop = FALSE]
    rownames(comparisons) <- NULL
    collapse_rows <- function(mat) {
      if (!base::anyDuplicated(base::rownames(mat))) {
        return(mat)
      }
      idx <- base::split(base::seq_len(base::nrow(mat)), base::rownames(mat))
      out <- base::t(base::vapply(
        idx,
        function(i) {
          vals <- base::colMeans(mat[i, , drop = FALSE], na.rm = TRUE)
          vals[base::is.nan(vals)] <- NA_real_
          vals
        },
        FUN.VALUE = base::numeric(base::ncol(mat))
      ))
      out <- out %>% base::as.matrix()
      out
    }
    collapse_cols <- function(mat) {
      if (!base::anyDuplicated(base::colnames(mat))) {
        return(mat)
      }
      idx <- base::split(base::seq_len(base::ncol(mat)), base::colnames(mat))
      out <- base::vapply(
        idx,
        function(i) {
          if (base::length(i) == 1) {
            return(mat[, i])
          }
          vals <- base::rowMeans(mat[, i, drop = FALSE], na.rm = TRUE)
          vals[base::is.nan(vals)] <- NA_real_
          vals
        },
        FUN.VALUE = base::numeric(base::nrow(mat))
      )
      out <- out %>% base::as.matrix()
      base::rownames(out) <- base::rownames(mat)
      out
    }
    resolve_group_labels <- function(info_dataset) {
      info_dataset <- info_dataset %>% base::as.data.frame(stringsAsFactors = FALSE)
      if (base::nrow(info_dataset) == 0) {
        return(base::character(0))
      }
      voi <- hcobject[["global_settings"]][["voi"]]
      voi <- base::intersect(voi, base::colnames(info_dataset))
      if (base::length(voi) > 0) {
        return(purrr::pmap(info_dataset[, voi, drop = FALSE], paste, sep = "-") %>% base::unlist())
      }
      base::as.character(info_dataset[[1]])
    }
    fc_limit <- suppressWarnings(base::as.numeric(hcobject[["global_settings"]][["range_GFC"]]))
    if (!base::is.finite(fc_limit) || base::is.na(fc_limit) || fc_limit <= 0) {
      fc_limit <- 2
    }
    pseudo_count <- 1e-08
    set_indices <- base::seq_len(base::length(hcobject[["layer_specific_outputs"]]))
    set_mats <- list()
    found <- stats::setNames(base::rep(FALSE, base::nrow(comparisons)), comparisons$comparison)
    for (z in set_indices) {
      set_name <- base::paste0("set", z)
      expr_mat <- hcobject[["layer_specific_outputs"]][[set_name]][["part1"]][["topvar"]]
      if (base::is.null(expr_mat)) {
        expr_mat <- hcobject[["data"]][[base::paste0(set_name, "_counts")]]
      }
      anno <- hcobject[["data"]][[base::paste0(set_name, "_anno")]]
      if (base::is.null(expr_mat) || base::is.null(anno)) {
        next
      }
      expr_mat <- expr_mat %>% base::as.matrix()
      if (base::nrow(expr_mat) == 0 || base::ncol(expr_mat) == 0) {
        next
      }
      mode(expr_mat) <- "numeric"
      samples <- base::intersect(base::colnames(expr_mat), base::rownames(anno))
      if (base::length(samples) == 0) {
        next
      }
      expr_mat <- expr_mat[, samples, drop = FALSE]
      anno <- anno[samples, , drop = FALSE]
      grpvar <- resolve_group_labels(anno)
      grpvar <- base::as.character(grpvar)
      if (base::length(grpvar) != base::length(samples)) {
        next
      }
      bad <- base::is.na(grpvar) | grpvar == ""
      if (base::any(bad)) {
        grpvar[bad] <- samples[bad]
      }
      if (isTRUE(hcobject[["global_settings"]][["data_in_log"]])) {
        expr_mat <- antilog(expr_mat, 2)
      }
      grp_levels <- base::unique(grpvar)
      set_mean <- base::vapply(
        grp_levels,
        function(grp) {
          idx <- base::which(grpvar == grp)
          if (base::length(idx) == 1) {
            return(expr_mat[, idx])
          }
          base::rowMeans(expr_mat[, idx, drop = FALSE], na.rm = TRUE)
        },
        FUN.VALUE = base::numeric(base::nrow(expr_mat))
      )
      set_mean <- set_mean %>% base::as.matrix()
      if (base::ncol(set_mean) == 0) {
        next
      }
      base::colnames(set_mean) <- grp_levels
      base::rownames(set_mean) <- base::rownames(expr_mat)
      overlap <- base::intersect(net_genes, base::rownames(set_mean))
      if (base::length(overlap) == 0) {
        next
      }
      for (i in base::seq_len(base::nrow(comparisons))) {
        num_grp <- comparisons$numerator[[i]]
        den_grp <- comparisons$denominator[[i]]
        cmp <- comparisons$comparison[[i]]
        if (!(num_grp %in% base::colnames(set_mean) && den_grp %in% base::colnames(set_mean))) {
          next
        }
        found[[cmp]] <- TRUE
        num_vals <- suppressWarnings(base::as.numeric(set_mean[overlap, num_grp]))
        den_vals <- suppressWarnings(base::as.numeric(set_mean[overlap, den_grp]))
        fc_vals <- suppressWarnings(base::log2((num_vals + pseudo_count) / (den_vals + pseudo_count)))
        fc_vals[!base::is.finite(fc_vals)] <- NA_real_
        fc_vals[fc_vals > fc_limit] <- fc_limit
        fc_vals[fc_vals < (-fc_limit)] <- -fc_limit
        set_full <- base::matrix(NA_real_, nrow = base::length(net_genes), ncol = 1, dimnames = list(net_genes, cmp))
        set_full[overlap, 1] <- fc_vals
        set_mats[[base::length(set_mats) + 1]] <- set_full
      }
    }
    if (base::length(set_mats) == 0 || !base::any(found)) {
      stop(
        "Could not reconstruct FC heatmap for knowledge network. ",
        "Please re-run `hc_upstream_inference(..., activity_input = 'fc', fc_comparisons = ...)` on this object."
      )
    }
    gene_fc <- base::do.call(base::cbind, set_mats)
    gene_fc <- collapse_rows(gene_fc)
    gene_fc <- collapse_cols(gene_fc)
    keep_rows <- base::rowSums(!base::is.na(gene_fc)) > 0
    gene_fc <- gene_fc[keep_rows, , drop = FALSE]
    requested_order <- base::as.character(comparisons$comparison)
    col_order <- requested_order[requested_order %in% base::colnames(gene_fc)]
    col_order <- base::c(col_order, base::setdiff(base::colnames(gene_fc), col_order))
    gene_fc <- gene_fc[, col_order, drop = FALSE]

    module_out <- list()
    for (cl in cluster_order) {
      genes <- dplyr::filter(cluster_info, color == cl) %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(split = ",") %>%
        base::unlist()
      genes <- base::unique(base::as.character(genes))
      genes <- base::intersect(genes, base::rownames(gene_fc))
      if (base::length(genes) == 0) {
        next
      }
      vals <- gene_fc[genes, , drop = FALSE]
      m <- base::colMeans(vals, na.rm = TRUE)
      m[base::is.nan(m)] <- NA_real_
      module_out[[cl]] <- m
    }
    if (base::length(module_out) == 0) {
      stop("No module-level FC means could be computed for the knowledge network heatmap.")
    }
    module_mat <- base::do.call(base::rbind, module_out)
    module_mat <- module_mat %>% base::as.matrix()
    list(mat = module_mat, col_order = col_order)
  }

  if (identical(upstream_activity_input, "fc")) {
    if (base::is.null(upstream_module_heatmap_name) || base::length(upstream_module_heatmap_name) == 0) {
      upstream_module_heatmap_name <- "FC"
    }
    upstream_module_heatmap_name <- "FC"
    if (base::is.null(upstream_module_heatmap_mat) || base::nrow(upstream_module_heatmap_mat) == 0 || base::ncol(upstream_module_heatmap_mat) == 0) {
      rebuilt_fc <- build_fc_module_heatmap_from_settings(
        fc_comparisons = if (base::is.list(upstream_settings)) upstream_settings[["fc_comparisons"]] else NULL,
        cluster_order = cluster_order,
        cluster_info = cluster_info
      )
      upstream_module_heatmap_mat <- rebuilt_fc$mat
      upstream_module_heatmap_col_order <- rebuilt_fc$col_order
    }
  }

  enrich_df <- enrich_df[base::as.character(enrich_df$cluster) %in% cluster_order, , drop = FALSE]
  upstream_df <- upstream_df[base::as.character(upstream_df$cluster) %in% cluster_order, , drop = FALSE]
  if (base::nrow(enrich_df) == 0 && base::nrow(upstream_df) == 0) {
    stop("No enrichment/upstream rows remain after applying cluster filter.")
  }

  clamp_q <- function(x) {
    x <- suppressWarnings(base::as.numeric(x))
    x[base::is.na(x) | x <= 0] <- 1e-300
    x
  }

  limit_per_module <- function(df, cluster_col, q_col, n_max) {
    if (base::is.null(n_max) || base::nrow(df) == 0) {
      return(df)
    }
    split_idx <- base::split(base::seq_len(base::nrow(df)), base::as.character(df[[cluster_col]]))
    kept <- base::lapply(split_idx, function(idx) {
      sub <- df[idx, , drop = FALSE]
      qv <- clamp_q(sub[[q_col]])
      ord <- base::order(qv)
      sub[ord[base::seq_len(base::min(base::length(ord), n_max))], , drop = FALSE]
    })
    out <- base::do.call(base::rbind, kept)
    base::rownames(out) <- NULL
    out
  }

  enrich_df <- limit_per_module(enrich_df, "cluster", "qvalue", max_enrichment_per_module)
  upstream_df <- limit_per_module(upstream_df, "cluster", "qvalue", max_upstream_per_module)
  if (base::nrow(enrich_df) == 0 && base::nrow(upstream_df) == 0) {
    stop("No rows remain after applying per-module limits.")
  }

  cluster_idx_map <- stats::setNames(base::seq_along(cluster_order), cluster_order)
  x_module <- 0
  x_enrichment <- 1.05
  x_upstream <- 2.75
  x_limits <- c(-0.25, 3.75)
  module_nodes <- base::data.frame(
    id = base::paste0("M::", cluster_order),
    key = cluster_order,
    label = base::as.character(module_label_map[cluster_order]),
    node_type = "module",
    x = x_module,
    y = base::rev(base::seq_len(base::length(cluster_order))),
    node_color = cluster_order,
    stringsAsFactors = FALSE
  )
  module_y_map <- stats::setNames(module_nodes$y, module_nodes$key)

  interleaved_order <- function(meta_df, group_col) {
    if (base::is.null(meta_df) || base::nrow(meta_df) == 0) {
      return(base::character(0))
    }
    group_levels <- base::unique(base::as.character(meta_df[[group_col]]))
    bucket_levels <- base::sort(base::unique(meta_df$first_cluster_idx))
    ordered <- base::character(0)
    for (b in bucket_levels) {
      bucket <- meta_df[meta_df$first_cluster_idx == b, , drop = FALSE]
      group_lists <- base::lapply(group_levels, function(g) {
        x <- bucket[base::as.character(bucket[[group_col]]) == g, , drop = FALSE]
        if (base::nrow(x) == 0) {
          return(base::character(0))
        }
        x <- x[base::order(x$best_q, -x$hit_count, x$mean_cluster_idx, x$term), , drop = FALSE]
        x$node_key
      })
      base::names(group_lists) <- group_levels
      max_len <- base::max(base::lengths(group_lists))
      for (k in base::seq_len(max_len)) {
        for (g in group_levels) {
          if (base::length(group_lists[[g]]) >= k) {
            ordered <- base::c(ordered, group_lists[[g]][k])
          }
        }
      }
    }
    ordered <- base::unique(ordered)
    base::c(ordered, base::setdiff(meta_df$node_key, ordered))
  }

  enrich_nodes <- base::data.frame(
    id = base::character(0),
    key = base::character(0),
    label = base::character(0),
    node_type = base::character(0),
    node_subtype = base::character(0),
    x = base::numeric(0),
    y = base::numeric(0),
    stringsAsFactors = FALSE
  )
  enrich_edges <- base::data.frame(
    from_id = base::character(0),
    to_id = base::character(0),
    from_key = base::character(0),
    to_key = base::character(0),
    edge_type = base::character(0),
    edge_subtype = base::character(0),
    qvalue = base::numeric(0),
    direction = base::character(0),
    x = base::numeric(0),
    y = base::numeric(0),
    xend = base::numeric(0),
    yend = base::numeric(0),
    stringsAsFactors = FALSE
  )
  if (base::nrow(enrich_df) > 0) {
    enrich_df$database <- base::as.character(enrich_df$database)
    enrich_df$cluster <- base::as.character(enrich_df$cluster)
    enrich_df$term <- base::as.character(enrich_df$term)
    enrich_df$qvalue <- clamp_q(enrich_df$qvalue)
    enrich_df$cluster_idx <- cluster_idx_map[enrich_df$cluster]
    enrich_df$cluster_idx[base::is.na(enrich_df$cluster_idx)] <- base::length(cluster_order) + 1
    enrich_df$node_key <- base::paste(enrich_df$database, enrich_df$term, sep = "||")

    split_en <- base::split(base::seq_len(base::nrow(enrich_df)), enrich_df$node_key)
    enrich_meta <- base::do.call(base::rbind, base::lapply(base::names(split_en), function(k) {
      idx <- split_en[[k]]
      sub <- enrich_df[idx, , drop = FALSE]
      base::data.frame(
        node_key = k,
        database = sub$database[[1]],
        term = sub$term[[1]],
        first_cluster_idx = base::min(sub$cluster_idx, na.rm = TRUE),
        mean_cluster_idx = base::mean(sub$cluster_idx, na.rm = TRUE),
        best_q = base::min(sub$qvalue, na.rm = TRUE),
        hit_count = base::length(base::unique(sub$cluster)),
        stringsAsFactors = FALSE
      )
    }))
    enrich_meta <- enrich_meta[base::order(enrich_meta$first_cluster_idx, enrich_meta$best_q, -enrich_meta$hit_count, enrich_meta$term), , drop = FALSE]
    en_order <- interleaved_order(enrich_meta, "database")
    enrich_meta <- enrich_meta[base::match(en_order, enrich_meta$node_key), , drop = FALSE]
    n_en <- base::nrow(enrich_meta)
    y_en <- if (n_en <= 1) {
      base::mean(module_nodes$y)
    } else {
      base::seq(base::max(module_nodes$y), base::min(module_nodes$y), length.out = n_en)
    }
    enrich_nodes <- base::data.frame(
      id = base::paste0("E::", enrich_meta$node_key),
      key = enrich_meta$node_key,
      label = base::paste0("[", enrich_meta$database, "] ", enrich_meta$term),
      node_type = "enrichment",
      node_subtype = enrich_meta$database,
      x = x_enrichment,
      y = y_en,
      stringsAsFactors = FALSE
    )
    enrich_y_map <- stats::setNames(enrich_nodes$y, enrich_nodes$key)
    enrich_edges <- base::data.frame(
      from_id = base::paste0("M::", enrich_df$cluster),
      to_id = base::paste0("E::", enrich_df$node_key),
      from_key = enrich_df$cluster,
      to_key = enrich_df$node_key,
      edge_type = "enrichment",
      edge_subtype = enrich_df$database,
      qvalue = enrich_df$qvalue,
      direction = NA_character_,
      stringsAsFactors = FALSE
    )
    enrich_edges$x <- x_module
    enrich_edges$y <- module_y_map[enrich_edges$from_key]
    enrich_edges$xend <- x_enrichment
    enrich_edges$yend <- enrich_y_map[enrich_edges$to_key]
  }

  upstream_nodes <- base::data.frame(
    id = base::character(0),
    key = base::character(0),
    label = base::character(0),
    node_type = base::character(0),
    node_subtype = base::character(0),
    x = base::numeric(0),
    y = base::numeric(0),
    stringsAsFactors = FALSE
  )
  upstream_edges <- base::data.frame(
    from_id = base::character(0),
    to_id = base::character(0),
    from_key = base::character(0),
    to_key = base::character(0),
    edge_type = base::character(0),
    edge_subtype = base::character(0),
    qvalue = base::numeric(0),
    direction = base::character(0),
    x = base::numeric(0),
    y = base::numeric(0),
    xend = base::numeric(0),
    yend = base::numeric(0),
    stringsAsFactors = FALSE
  )
  if (base::nrow(upstream_df) > 0) {
    upstream_df$resource <- base::as.character(upstream_df$resource)
    upstream_df$cluster <- base::as.character(upstream_df$cluster)
    upstream_df$term <- base::as.character(upstream_df$term)
    upstream_df$direction <- base::as.character(upstream_df$direction)
    upstream_df$qvalue <- clamp_q(upstream_df$qvalue)
    upstream_df$cluster_idx <- cluster_idx_map[upstream_df$cluster]
    upstream_df$cluster_idx[base::is.na(upstream_df$cluster_idx)] <- base::length(cluster_order) + 1
    upstream_df$node_key <- base::paste(upstream_df$resource, upstream_df$term, sep = "||")

    split_up <- base::split(base::seq_len(base::nrow(upstream_df)), upstream_df$node_key)
    upstream_meta <- base::do.call(base::rbind, base::lapply(base::names(split_up), function(k) {
      idx <- split_up[[k]]
      sub <- upstream_df[idx, , drop = FALSE]
      base::data.frame(
        node_key = k,
        resource = sub$resource[[1]],
        term = sub$term[[1]],
        first_cluster_idx = base::min(sub$cluster_idx, na.rm = TRUE),
        mean_cluster_idx = base::mean(sub$cluster_idx, na.rm = TRUE),
        best_q = base::min(sub$qvalue, na.rm = TRUE),
        hit_count = base::length(base::unique(sub$cluster)),
        stringsAsFactors = FALSE
      )
    }))
    upstream_meta <- upstream_meta[base::order(upstream_meta$first_cluster_idx, upstream_meta$best_q, -upstream_meta$hit_count, upstream_meta$term), , drop = FALSE]
    up_order <- interleaved_order(upstream_meta, "resource")
    upstream_meta <- upstream_meta[base::match(up_order, upstream_meta$node_key), , drop = FALSE]
    n_up <- base::nrow(upstream_meta)
    y_up <- if (n_up <= 1) {
      base::mean(module_nodes$y)
    } else {
      base::seq(base::max(module_nodes$y), base::min(module_nodes$y), length.out = n_up)
    }
    upstream_nodes <- base::data.frame(
      id = base::paste0("U::", upstream_meta$node_key),
      key = upstream_meta$node_key,
      label = base::paste0(upstream_meta$term, " [", upstream_meta$resource, "]"),
      node_type = "upstream",
      node_subtype = upstream_meta$resource,
      x = x_upstream,
      y = y_up,
      stringsAsFactors = FALSE
    )
    upstream_y_map <- stats::setNames(upstream_nodes$y, upstream_nodes$key)
    upstream_edges <- base::data.frame(
      from_id = base::paste0("M::", upstream_df$cluster),
      to_id = base::paste0("U::", upstream_df$node_key),
      from_key = upstream_df$cluster,
      to_key = upstream_df$node_key,
      edge_type = "upstream",
      edge_subtype = upstream_df$resource,
      qvalue = upstream_df$qvalue,
      direction = upstream_df$direction,
      stringsAsFactors = FALSE
    )
    upstream_edges$x <- x_module
    upstream_edges$y <- module_y_map[upstream_edges$from_key]
    upstream_edges$xend <- x_upstream
    upstream_edges$yend <- upstream_y_map[upstream_edges$to_key]
  }

  edges <- base::rbind(enrich_edges, upstream_edges)
  if (base::nrow(edges) == 0) {
    stop("No edges to plot after filtering.")
  }

  edges$qvalue <- clamp_q(edges$qvalue)
  edges$neglog10_q <- -base::log10(base::pmax(edges$qvalue, 1e-300))
  edges$neglog10_q[!base::is.finite(edges$neglog10_q)] <- 0
  edges$neglog10_q <- base::pmin(12, edges$neglog10_q)
  edges$edge_alpha <- 0.24 + (0.56 * base::pmin(1, edges$neglog10_q / 6))

  normalize_db <- function(x) {
    x <- base::as.character(x)
    low <- base::tolower(x)
    out <- x
    out[low %in% c("go", "gobp", "go_bp")] <- "Go"
    out[low %in% c("kegg")] <- "Kegg"
    out[low %in% c("hallmark")] <- "Hallmark"
    out[low %in% c("reactome")] <- "Reactome"
    out
  }
  edges$edge_subtype <- normalize_db(edges$edge_subtype)

  to_dir <- function(x) {
    x <- base::tolower(base::as.character(x))
    if (x %in% c("activated", "up", "positive")) {
      return("activated")
    }
    if (x %in% c("inhibited", "down", "negative")) {
      return("inhibited")
    }
    ""
  }

  edges$direction_class <- ifelse(
    edges$edge_type == "upstream",
    base::vapply(edges$direction, to_dir, FUN.VALUE = base::character(1)),
    ""
  )

  edges$line_group <- base::vapply(base::seq_len(base::nrow(edges)), function(i) {
    if (identical(edges$edge_type[[i]], "enrichment")) {
      return(base::paste0(edges$edge_subtype[[i]], " enrichment"))
    }
    d <- edges$direction_class[[i]]
    if (d == "activated") {
      return(base::paste0(edges$edge_subtype[[i]], " activated"))
    }
    if (d == "inhibited") {
      return(base::paste0(edges$edge_subtype[[i]], " inhibited"))
    }
    base::paste0(edges$edge_subtype[[i]], " (no direction)")
  }, FUN.VALUE = base::character(1))

  line_color_map <- c(
    "Go enrichment" = "#4E79A7",
    "Kegg enrichment" = "#F28E2B",
    "Hallmark enrichment" = "#59A14F",
    "Reactome enrichment" = "#E15759",
    "TF activated" = "#D62728",
    "TF inhibited" = "#1F77B4",
    "Pathway activated" = "#C03D3E",
    "Pathway inhibited" = "#2C73B8",
    "TF (no direction)" = "#3B7EA1",
    "Pathway (no direction)" = "#B26B2C"
  )
  missing_groups <- base::setdiff(base::unique(edges$line_group), base::names(line_color_map))
  if (base::length(missing_groups) > 0) {
    line_color_map <- base::c(
      line_color_map,
      stats::setNames(base::rep("grey65", base::length(missing_groups)), missing_groups)
    )
  }
  line_group_order <- c(
    "Go enrichment",
    "Kegg enrichment",
    "Hallmark enrichment",
    "Reactome enrichment",
    "TF activated",
    "TF inhibited",
    "Pathway activated",
    "Pathway inhibited",
    "TF (no direction)",
    "Pathway (no direction)"
  )
  line_group_order <- base::c(
    line_group_order[line_group_order %in% base::names(line_color_map)],
    base::setdiff(base::names(line_color_map), line_group_order)
  )

  build_hc_heatmap_data <- function(cluster_order,
                                    module_label_map,
                                    cluster_info,
                                    stored_hm,
                                    override_mat = NULL,
                                    override_col_order = base::character(0),
                                    value_name = "GFC",
                                    gfc_scale_limits = c(-2, 2),
                                    gfc_scale_ticks = NULL) {
    if (base::is.null(gfc_scale_ticks) ||
        !base::is.list(gfc_scale_ticks) ||
        !all(c("breaks", "labels") %in% base::names(gfc_scale_ticks))) {
      gfc_scale_ticks <- compute_scale_ticks(gfc_scale_limits)
    }

    extract_heatmap <- function(stored_hm, cluster_info, cluster_order, module_label_map, override_mat) {
      if (!base::is.null(override_mat) && base::nrow(override_mat) > 0 && base::ncol(override_mat) > 0) {
        hm_mat <- override_mat %>% base::as.matrix()
        mode(hm_mat) <- "numeric"
        return(hm_mat)
      }
      hm_mat <- tryCatch(stored_hm@ht_list[[1]]@matrix, error = function(e) NULL)
      if (!base::is.null(hm_mat) && base::nrow(hm_mat) > 0) {
        hm_mat <- hm_mat %>% base::as.matrix()
        inv_map <- stats::setNames(base::names(module_label_map), base::as.character(module_label_map))
        mapped_rows <- ifelse(base::rownames(hm_mat) %in% base::names(inv_map), inv_map[base::rownames(hm_mat)], base::rownames(hm_mat))
        base::rownames(hm_mat) <- mapped_rows
        return(hm_mat)
      }
      gfc_all <- hcobject[["integrated_output"]][["GFC_all_layers"]]
      if (base::is.null(gfc_all) || base::nrow(gfc_all) == 0 || base::ncol(gfc_all) < 2) {
        return(NULL)
      }
      gene_col <- if ("Gene" %in% base::colnames(gfc_all)) "Gene" else base::colnames(gfc_all)[base::ncol(gfc_all)]
      value_cols <- base::setdiff(base::colnames(gfc_all), gene_col)
      out <- base::lapply(cluster_order, function(cl) {
        genes <- dplyr::filter(cluster_info, color == cl) %>%
          dplyr::pull(., "gene_n") %>%
          base::strsplit(split = ",") %>%
          base::unlist()
        genes <- base::unique(base::as.character(genes))
        if (base::length(genes) == 0) {
          return(NULL)
        }
        tmp <- gfc_all[gfc_all[[gene_col]] %in% genes, value_cols, drop = FALSE]
        if (base::nrow(tmp) == 0) {
          return(NULL)
        }
        vals <- tmp %>% base::as.matrix()
        mode(vals) <- "numeric"
        base::colMeans(vals, na.rm = TRUE)
      })
      names(out) <- cluster_order
      out <- out[!base::vapply(out, base::is.null, FUN.VALUE = base::logical(1))]
      if (base::length(out) == 0) {
        return(NULL)
      }
      hm <- base::do.call(base::rbind, out)
      hm %>% base::as.matrix()
    }

    heatmap_mat <- extract_heatmap(stored_hm, cluster_info, cluster_order, module_label_map, override_mat)
    if (base::is.null(heatmap_mat) || base::nrow(heatmap_mat) == 0) {
      return(NULL)
    }

    if (base::length(override_col_order) > 0) {
      hm_cols <- base::colnames(heatmap_mat)
      ordered_cols <- base::as.character(override_col_order)
      ordered_cols <- ordered_cols[ordered_cols %in% hm_cols]
      ordered_cols <- base::c(ordered_cols, base::setdiff(hm_cols, ordered_cols))
      if (base::length(ordered_cols) > 0) {
        heatmap_mat <- heatmap_mat[, ordered_cols, drop = FALSE]
      }
    } else if (!base::is.null(stored_hm) && base::ncol(heatmap_mat) > 0) {
      hm_cols <- base::colnames(heatmap_mat)
      col_order_ref <- tryCatch(ComplexHeatmap::column_order(stored_hm), error = function(e) NULL)
      if (base::is.list(col_order_ref) && base::length(col_order_ref) > 0) {
        col_order_ref <- col_order_ref[[1]]
      }
      ordered_cols <- NULL
      if (base::is.numeric(col_order_ref) && base::length(col_order_ref) == base::length(hm_cols)) {
        ordered_cols <- hm_cols[col_order_ref]
      } else if (base::is.character(col_order_ref) && base::length(col_order_ref) > 0) {
        ordered_cols <- col_order_ref
      } else {
        stored_cols <- tryCatch(base::colnames(stored_hm@ht_list[[1]]@matrix), error = function(e) NULL)
        if (!base::is.null(stored_cols) && base::length(stored_cols) > 0) {
          ordered_cols <- stored_cols
        }
      }
      if (!base::is.null(ordered_cols)) {
        ordered_cols <- base::as.character(ordered_cols)
        ordered_cols <- ordered_cols[ordered_cols %in% hm_cols]
        ordered_cols <- base::c(ordered_cols, base::setdiff(hm_cols, ordered_cols))
        heatmap_mat <- heatmap_mat[, ordered_cols, drop = FALSE]
      }
    }

    keep_clusters <- cluster_order[cluster_order %in% base::rownames(heatmap_mat)]
    if (base::length(keep_clusters) == 0) {
      return(NULL)
    }
    heatmap_mat <- heatmap_mat[keep_clusters, , drop = FALSE]
    module_labels <- base::as.character(module_label_map[keep_clusters])
    module_labels[is.na(module_labels) | module_labels == ""] <- keep_clusters[is.na(module_labels) | module_labels == ""]
    base::rownames(heatmap_mat) <- module_labels

    cdf <- dplyr::filter(cluster_info, color %in% keep_clusters) %>%
      dplyr::distinct(color, .keep_all = TRUE)
    if ("gene_no" %in% base::colnames(cdf)) {
      gene_counts <- suppressWarnings(base::as.numeric(cdf$gene_no))
    } else if ("gene_n" %in% base::colnames(cdf)) {
      gene_counts <- base::vapply(
        base::as.character(cdf$gene_n),
        function(x) {
          if (base::is.na(x) || x == "") {
            return(0)
          }
          parts <- base::trimws(base::unlist(base::strsplit(x, ",", fixed = TRUE)))
          base::sum(parts != "")
        },
        FUN.VALUE = base::numeric(1)
      )
    } else {
      gene_counts <- base::rep(0, base::nrow(cdf))
    }
    gene_counts <- stats::setNames(gene_counts, base::as.character(cdf$color))[keep_clusters]
    gene_counts[base::is.na(gene_counts)] <- 0
    value_name <- base::as.character(value_name[[1]])
    if (base::is.na(value_name) || value_name == "") {
      value_name <- "GFC"
    }
    value_name_upper <- base::toupper(value_name)
    if (!identical(value_name_upper, "FC")) {
      value_name <- "GFC"
    }
    stored_gfc_colors <- cluster_calc[["gfc_colors"]]
    if (!base::is.character(stored_gfc_colors) ||
        base::length(stored_gfc_colors) < 2 ||
        any(base::is.na(stored_gfc_colors)) ||
        any(stored_gfc_colors == "")) {
      stored_gfc_colors <- NULL
    } else {
      stored_gfc_colors <- base::as.character(stored_gfc_colors)
    }
    if (base::is.null(stored_gfc_colors)) {
      stored_gfc_colors <- base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
    }

    list(
      mat = heatmap_mat,
      keep_clusters = keep_clusters,
      module_labels = module_labels,
      module_colors = keep_clusters,
      gene_counts = as.numeric(gene_counts),
      value_name = value_name,
      gfc_colors = stored_gfc_colors,
      scale_limits = gfc_scale_limits,
      scale_breaks = gfc_scale_ticks$breaks,
      scale_labels = gfc_scale_ticks$labels
    )
  }

  build_hc_heatmap_plot <- function(hm_data, overall_plot_scale) {
    mat <- hm_data$mat
    n_r <- base::nrow(mat)
    n_c <- base::ncol(mat)
    if (n_r == 0 || n_c == 0) {
      return(NULL)
    }

    module_y <- base::rev(base::seq_len(n_r))
    row_levels <- base::rownames(mat)
    col_levels <- base::colnames(mat)

    hm_long <- base::as.data.frame(base::as.table(mat), stringsAsFactors = FALSE)
    base::colnames(hm_long) <- c("module_label", "condition", "value")
    hm_long$module_label <- base::as.character(hm_long$module_label)
    hm_long$condition <- base::as.character(hm_long$condition)
    hm_long$row_idx <- base::match(hm_long$module_label, row_levels)
    hm_long$col_idx <- base::match(hm_long$condition, col_levels)
    hm_long$y <- module_y[hm_long$row_idx]
    hm_long$x <- hm_long$col_idx
    hm_long$value <- suppressWarnings(base::as.numeric(hm_long$value))

    module_df <- base::data.frame(
      cluster = hm_data$keep_clusters,
      module_label = hm_data$module_labels,
      module_color = hm_data$module_colors,
      y = module_y,
      stringsAsFactors = FALSE
    )

    gfc_pal <- grDevices::colorRampPalette(hm_data$gfc_colors)(51)
    stored_module_label_fontsize <- cluster_calc[["module_label_fontsize"]]
    if (!base::is.numeric(stored_module_label_fontsize) ||
        base::length(stored_module_label_fontsize) != 1 ||
        base::is.na(stored_module_label_fontsize) ||
        stored_module_label_fontsize <= 0) {
      stored_module_label_fontsize <- NULL
    }
    stored_module_box_to_cell_ratio <- cluster_calc[["module_box_to_cell_ratio"]]
    if (!base::is.numeric(stored_module_box_to_cell_ratio) ||
        base::length(stored_module_box_to_cell_ratio) != 1 ||
        base::is.na(stored_module_box_to_cell_ratio) ||
        stored_module_box_to_cell_ratio <= 0) {
      stored_module_box_to_cell_ratio <- NULL
    }
    module_box_width_units <- if (!base::is.null(stored_module_box_to_cell_ratio)) {
      base::max(0.68, base::min(1.05, stored_module_box_to_cell_ratio))
    } else {
      0.82
    }
    module_box_height_units <- 0.96
    module_box_gap_units <- 0.16
    heatmap_right_edge <- n_c + 0.5
    x_mod <- heatmap_right_edge + module_box_gap_units + (module_box_width_units / 2)
    x_right <- heatmap_right_edge + module_box_gap_units + module_box_width_units + 0.26
    font_module <- if (!base::is.null(stored_module_label_fontsize)) {
      stored_module_label_fontsize
    } else {
      base::max(6.8, base::min(12, 9.4 * overall_plot_scale))
    }
    module_text_size <- base::max(2.2, base::min(5.2, font_module / 2.845276))

    fill_scale <- ggplot2::scale_fill_gradientn(
      colors = gfc_pal,
      limits = hm_data$scale_limits,
      breaks = hm_data$scale_breaks,
      labels = hm_data$scale_labels,
      oob = scales::squish,
      name = hm_data$value_name,
      guide = "none"
    )

    ggplot2::ggplot(hm_long, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_tile(color = "black", linewidth = 0.3) +
      ggplot2::geom_tile(
        data = module_df,
        ggplot2::aes(x = x_mod, y = y),
        inherit.aes = FALSE,
        fill = module_df$module_color,
        color = "black",
        width = module_box_width_units,
        height = module_box_height_units
      ) +
      ggplot2::geom_text(
        data = module_df,
        ggplot2::aes(x = x_mod, y = y, label = module_label),
        inherit.aes = FALSE,
        color = "white",
        fontface = "bold",
        size = module_text_size
      ) +
      fill_scale +
      ggplot2::scale_x_continuous(
        limits = c(0.5, x_right),
        breaks = base::seq_len(n_c),
        labels = col_levels,
        expand = ggplot2::expansion(mult = 0, add = 0)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0.5, n_r + 0.5),
        breaks = NULL,
        expand = ggplot2::expansion(mult = 0, add = 0)
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme_void(base_size = 11 * overall_plot_scale) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.margin = grid::unit(c(1.2, 0.25, 1.2, 2.6) * overall_plot_scale, "mm")
      )
  }

  hc_heatmap_data <- build_hc_heatmap_data(
    cluster_order = cluster_order,
    module_label_map = module_label_map,
    cluster_info = cluster_info,
    stored_hm = stored_hm,
    override_mat = upstream_module_heatmap_mat,
    override_col_order = upstream_module_heatmap_col_order,
    value_name = upstream_module_heatmap_name,
    gfc_scale_limits = resolved_gfc_scale_limits,
    gfc_scale_ticks = resolved_gfc_scale_ticks
  )
  if (base::is.null(hc_heatmap_data)) {
    stop("Unable to build hCoCena heatmap panel for the network plot.")
  }
  hc_heatmap_plot <- build_hc_heatmap_plot(hc_heatmap_data, overall_plot_scale = overall_plot_scale)
  if (base::is.null(hc_heatmap_plot)) {
    stop("Unable to build hCoCena heatmap plot for the network plot.")
  }
  # Force exact y-alignment between heatmap rows and module rows in the network.
  cluster_order <- base::as.character(hc_heatmap_data$keep_clusters)
  module_label_map <- stats::setNames(
    base::as.character(hc_heatmap_data$module_labels),
    cluster_order
  )
  module_nodes <- module_nodes[module_nodes$key %in% cluster_order, , drop = FALSE]
  module_nodes <- module_nodes[base::match(cluster_order, module_nodes$key), , drop = FALSE]
  module_nodes$label <- base::as.character(module_label_map[module_nodes$key])
  module_nodes$y <- base::rev(base::seq_len(base::length(cluster_order)))
  module_y_map <- stats::setNames(module_nodes$y, module_nodes$key)
  if (base::nrow(enrich_edges) > 0) {
    enrich_edges$y <- module_y_map[base::as.character(enrich_edges$from_key)]
  }
  if (base::nrow(upstream_edges) > 0) {
    upstream_edges$y <- module_y_map[base::as.character(upstream_edges$from_key)]
  }
  edges$y <- module_y_map[base::as.character(edges$from_key)]
  missing_edge_y <- base::is.na(edges$y)
  if (base::any(missing_edge_y)) {
    edges <- edges[!missing_edge_y, , drop = FALSE]
    if (base::nrow(edges) == 0) {
      stop("No edges remain after aligning module rows with heatmap rows.")
    }
  }

  trim_label <- function(x, max_chars = 52) {
    x <- base::as.character(x)
    too_long <- base::nchar(x) > max_chars
    x[too_long] <- base::paste0(base::substr(x[too_long], 1, max_chars - 1), "...")
    x
  }

  line_breaks <- c(1.3, 2, 3, 5, 8)
  max_line_metric <- base::max(edges$neglog10_q, na.rm = TRUE)
  line_breaks <- line_breaks[line_breaks <= (max_line_metric + 1e-8)]
  if (base::length(line_breaks) == 0) {
    line_breaks <- pretty(c(0, max_line_metric), n = 3)
    line_breaks <- line_breaks[line_breaks > 0]
  }
  if (base::length(line_breaks) == 0) {
    line_breaks <- base::c(1)
  }

  overview_title <- base::paste0(
    "Module knowledge network (enrichment: ", enrichment_mode,
    ", upstream: ", upstream_mode, ")"
  )
  network_subtitle <- "Squares = modules, circles = terms; line width reflects significance; arrows = activation; T-end = inhibition"
  page_title_for_focus <- function(focus_cluster = NULL) {
    if (base::is.null(focus_cluster)) {
      return(overview_title)
    }
    cl <- base::as.character(focus_cluster[[1]])
    lbl <- module_label_map[[cl]]
    if (base::is.null(lbl) || base::is.na(lbl) || lbl == "") {
      lbl <- cl
    }
    base::paste0(overview_title, " - focus on ", lbl)
  }

  build_network_plot <- function(focus_cluster = NULL,
                                 show_headers = FALSE,
                                 show_titles = FALSE) {
    edges_plot <- edges
    modules_plot <- module_nodes
    enrich_plot <- enrich_nodes
    upstream_plot <- upstream_nodes

    focus_label <- NULL
    if (!base::is.null(focus_cluster)) {
      focus_cluster <- base::as.character(focus_cluster[[1]])
      if (focus_cluster %in% cluster_order) {
        focus_label <- module_label_map[[focus_cluster]]
      } else {
        focus_cluster <- NULL
      }
    }

    if (base::is.null(focus_cluster)) {
      edges_plot$is_focus <- TRUE
      edges_plot$line_group_plot <- edges_plot$line_group
      edges_plot$alpha_plot <- edges_plot$edge_alpha
      modules_plot$is_focus <- TRUE
      enrich_plot$is_focus <- TRUE
      upstream_plot$is_focus <- TRUE
    } else {
      edges_plot$is_focus <- base::as.character(edges_plot$from_key) == focus_cluster
      edges_plot$line_group_plot <- ifelse(edges_plot$is_focus, edges_plot$line_group, ".greyed")
      edges_plot$alpha_plot <- ifelse(edges_plot$is_focus, edges_plot$edge_alpha, 0.08)
      modules_plot$is_focus <- base::as.character(modules_plot$key) == focus_cluster
      focus_enrichment <- base::unique(edges_plot$to_id[edges_plot$is_focus & edges_plot$edge_type == "enrichment"])
      focus_upstream <- base::unique(edges_plot$to_id[edges_plot$is_focus & edges_plot$edge_type == "upstream"])
      enrich_plot$is_focus <- enrich_plot$id %in% focus_enrichment
      upstream_plot$is_focus <- upstream_plot$id %in% focus_upstream
    }

    cmap <- line_color_map
    cmap[[".greyed"]] <- "grey83"
    line_group_labels <- stats::setNames(base::names(cmap), base::names(cmap))
    line_group_labels[["TF activated"]] <- "TF activated  ->"
    line_group_labels[["TF inhibited"]] <- "TF inhibited  -|- (T-end)"
    line_group_labels[["Pathway activated"]] <- "Pathway activated  ->"
    line_group_labels[["Pathway inhibited"]] <- "Pathway inhibited  -|- (T-end)"
    color_breaks <- base::unique(edges_plot$line_group_plot[edges_plot$line_group_plot != ".greyed"])
    color_breaks <- line_group_order[line_group_order %in% color_breaks]
    if (base::length(color_breaks) == 0) {
      color_breaks <- base::setdiff(base::unique(edges_plot$line_group_plot), ".greyed")
    }

    modules_plot$fill_plot <- ifelse(modules_plot$is_focus, modules_plot$node_color, "grey85")
    modules_plot$stroke_plot <- ifelse(modules_plot$is_focus, "black", "grey70")
    modules_plot$text_plot <- ifelse(modules_plot$is_focus, "black", "grey65")
    enrich_plot$stroke_plot <- ifelse(enrich_plot$is_focus, "#666666", "grey82")
    enrich_plot$text_plot <- ifelse(enrich_plot$is_focus, "black", "grey72")
    upstream_plot$stroke_plot <- ifelse(upstream_plot$is_focus, "#666666", "grey82")
    upstream_plot$text_plot <- ifelse(upstream_plot$is_focus, "black", "grey72")
    enrich_plot$label_plot <- trim_label(enrich_plot$label, max_chars = 42)
    upstream_plot$label_plot <- trim_label(upstream_plot$label, max_chars = 34)
    if (identical(label_mode, "upstream_only")) {
      enrich_plot$label_plot[] <- ""
    }
    if (identical(label_mode, "focus_only") && base::is.null(focus_cluster)) {
      enrich_plot$label_plot[] <- ""
      upstream_plot$label_plot[] <- ""
    }
    if (!base::is.null(focus_cluster)) {
      enrich_plot$label_plot[!enrich_plot$is_focus] <- ""
      upstream_plot$label_plot[!upstream_plot$is_focus] <- ""
    }

    module_label_df <- modules_plot[, c("x", "y", "label", "text_plot"), drop = FALSE]
    enrich_label_df <- enrich_plot[, c("x", "y", "label_plot", "text_plot"), drop = FALSE]
    upstream_label_df <- upstream_plot[, c("x", "y", "label_plot", "text_plot"), drop = FALSE]
    node_legend_df <- base::data.frame(x = c(0, 0), y = c(0, 0), node_type = c("Module node", "Term node"), stringsAsFactors = FALSE)
    label_nudge_y <- 0.22
    term_label_size <- if (base::is.null(focus_cluster)) {
      (7.5 * overall_plot_scale) / 3.2
    } else {
      (9.2 * overall_plot_scale) / 3.2
    }
    use_check_overlap <- base::is.null(focus_cluster)
    font_base <- 9.4 * overall_plot_scale
    edges_plot_non_activated <- edges_plot[
      !(edges_plot$edge_type == "upstream" & edges_plot$direction_class %in% c("activated", "inhibited")),
      ,
      drop = FALSE
    ]
    edges_plot_activated <- edges_plot[
      edges_plot$edge_type == "upstream" & edges_plot$direction_class == "activated",
      ,
      drop = FALSE
    ]
    edges_plot_inhibited <- edges_plot[
      edges_plot$edge_type == "upstream" & edges_plot$direction_class == "inhibited",
      ,
      drop = FALSE
    ]
    inhibited_stems <- edges_plot_inhibited[0, , drop = FALSE]
    inhibited_tbars <- base::data.frame(
      x = base::numeric(0),
      y = base::numeric(0),
      xend = base::numeric(0),
      yend = base::numeric(0),
      line_group_plot = base::character(0),
      alpha_plot = base::numeric(0),
      stringsAsFactors = FALSE
    )
    if (base::nrow(edges_plot_inhibited) > 0) {
      # Correct perpendicular direction for non-square plotting panels so T-ends
      # are visually orthogonal to the line on screen.
      dev_size <- tryCatch(grDevices::dev.size("in"), error = function(e) base::c(10, 7))
      if (base::length(dev_size) < 2 || base::any(!base::is.finite(dev_size)) || base::any(dev_size <= 0)) {
        dev_size <- base::c(10, 7)
      }
      panel_aspect <- (dev_size[[2]] * 0.905) / (dev_size[[1]] * 0.65)
      xrange <- base::max(1e-6, diff(x_limits))
      yrange <- base::max(1, base::length(cluster_order))
      perp_corr <- panel_aspect * (xrange / yrange)
      perp_corr <- base::max(1e-6, perp_corr)

      dx <- edges_plot_inhibited$xend - edges_plot_inhibited$x
      dy <- edges_plot_inhibited$yend - edges_plot_inhibited$y
      seg_len <- base::sqrt(dx^2 + dy^2)
      seg_len[seg_len == 0] <- 1
      ux <- dx / seg_len
      uy <- dy / seg_len
      t_offset <- 0.085
      t_center_x <- edges_plot_inhibited$xend - (ux * t_offset)
      t_center_y <- edges_plot_inhibited$yend - (uy * t_offset)
      inhibited_stems <- edges_plot_inhibited
      inhibited_stems$xend <- t_center_x
      inhibited_stems$yend <- t_center_y
      px_raw <- -dy * perp_corr
      py_raw <- dx / perp_corr
      p_len <- base::sqrt(px_raw^2 + py_raw^2)
      p_len[p_len == 0] <- 1
      px <- px_raw / p_len
      py <- py_raw / p_len
      bar_half <- 0.126
      inhibited_tbars <- base::data.frame(
        x = t_center_x + (px * bar_half),
        y = t_center_y + (py * bar_half),
        xend = t_center_x - (px * bar_half),
        yend = t_center_y - (py * bar_half),
        line_group_plot = edges_plot_inhibited$line_group_plot,
        alpha_plot = ifelse(edges_plot_inhibited$is_focus, 1, edges_plot_inhibited$alpha_plot),
        stringsAsFactors = FALSE
      )
    }
    inhibited_tbar_linewidth <- base::max(1.15, 1.7 * overall_plot_scale)

    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = edges_plot_non_activated,
        ggplot2::aes(
          x = x, y = y, xend = xend, yend = yend,
          color = line_group_plot,
          linewidth = neglog10_q,
          alpha = alpha_plot
        ),
        lineend = "round"
      ) +
      ggplot2::geom_segment(
        data = edges_plot_activated,
        ggplot2::aes(
          x = x, y = y, xend = xend, yend = yend,
          color = line_group_plot,
          linewidth = neglog10_q,
          alpha = alpha_plot
        ),
        lineend = "round",
        arrow = grid::arrow(
          type = "closed",
          length = grid::unit(2.2 * overall_plot_scale, "mm")
        )
      ) +
      ggplot2::geom_segment(
        data = inhibited_stems,
        ggplot2::aes(
          x = x, y = y, xend = xend, yend = yend,
          color = line_group_plot,
          linewidth = neglog10_q,
          alpha = alpha_plot
        ),
        lineend = "round"
      ) +
      ggplot2::scale_alpha_identity(guide = "none") +
      ggplot2::scale_color_manual(
        values = cmap,
        breaks = color_breaks,
        labels = line_group_labels[color_breaks],
        name = "Line color"
      ) +
      ggplot2::scale_linewidth_continuous(name = "-log10(qvalue)", breaks = line_breaks, range = c(0.25, 1.9)) +
      ggplot2::geom_point(
        data = modules_plot,
        ggplot2::aes(x = x, y = y),
        shape = 22,
        size = 3.2 * overall_plot_scale,
        fill = modules_plot$fill_plot,
        color = modules_plot$stroke_plot,
        stroke = 0.28
      ) +
      ggplot2::geom_point(
        data = enrich_plot,
        ggplot2::aes(x = x, y = y),
        shape = 21,
        size = 2.55 * overall_plot_scale,
        fill = "white",
        color = enrich_plot$stroke_plot,
        stroke = 0.3
      ) +
      ggplot2::geom_point(
        data = upstream_plot,
        ggplot2::aes(x = x, y = y),
        shape = 21,
        size = 2.55 * overall_plot_scale,
        fill = "white",
        color = upstream_plot$stroke_plot,
        stroke = 0.3
      ) +
      ggplot2::geom_segment(
        data = inhibited_tbars,
        ggplot2::aes(
          x = x, y = y, xend = xend, yend = yend,
          color = line_group_plot,
          alpha = alpha_plot
        ),
        linewidth = inhibited_tbar_linewidth,
        lineend = "butt",
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = module_label_df,
        ggplot2::aes(x = x - 0.05, y = y, label = label),
        hjust = 1,
        size = font_base / 3.2,
        fontface = "bold",
        color = module_label_df$text_plot
      ) +
      ggplot2::geom_text(
        data = enrich_label_df,
        ggplot2::aes(x = x, y = y + label_nudge_y, label = label_plot),
        hjust = 0.5,
        vjust = 0,
        size = term_label_size,
        color = enrich_label_df$text_plot,
        check_overlap = use_check_overlap
      ) +
      ggplot2::geom_text(
        data = upstream_label_df,
        ggplot2::aes(x = x, y = y + label_nudge_y, label = label_plot),
        hjust = 0.5,
        vjust = 0,
        size = term_label_size,
        color = upstream_label_df$text_plot,
        check_overlap = use_check_overlap
      ) +
      ggplot2::geom_point(
        data = node_legend_df,
        ggplot2::aes(x = x, y = y, shape = node_type),
        alpha = 0
      ) +
      ggplot2::scale_shape_manual(name = "Nodes", values = c("Module node" = 22, "Term node" = 21)) +
      ggplot2::guides(
        color = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1, linewidth = 1.8)),
        linewidth = ggplot2::guide_legend(order = 2),
        shape = ggplot2::guide_legend(
          order = 3,
          override.aes = list(
            alpha = 1,
            fill = c("grey60", "white"),
            color = c("black", "#666666"),
            size = c(4, 3.6),
            stroke = c(0.35, 0.35)
          )
        )
      ) +
      ggplot2::scale_x_continuous(
        limits = x_limits,
        breaks = NULL,
        expand = ggplot2::expansion(mult = 0, add = 0)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0.5, base::length(cluster_order) + 0.5),
        breaks = NULL,
        expand = ggplot2::expansion(mult = 0, add = 0)
      ) +
      ggplot2::coord_cartesian(
        clip = "off"
      ) +
      ggplot2::theme_void(base_size = 11 * overall_plot_scale) +
      ggplot2::theme(
        plot.title = if (isTRUE(show_titles)) {
          ggplot2::element_text(face = "bold", size = 13 * overall_plot_scale, hjust = 0.5)
        } else {
          ggplot2::element_blank()
        },
        plot.subtitle = if (isTRUE(show_titles)) {
          ggplot2::element_text(size = 9.8 * overall_plot_scale, hjust = 0.5, color = "grey30")
        } else {
          ggplot2::element_blank()
        },
        legend.title = ggplot2::element_text(size = 10 * overall_plot_scale, face = "bold"),
        legend.text = ggplot2::element_text(size = 8.7 * overall_plot_scale),
        legend.box = "vertical",
        legend.spacing.y = grid::unit(1.2 * overall_plot_scale, "mm"),
        plot.margin = grid::unit(
          c(
            if (isTRUE(show_titles)) 8 else 1.5,
            22,
            2,
            1.1
          ) * overall_plot_scale,
          "mm"
        )
      )

    if (isTRUE(show_titles)) {
      p <- p + ggplot2::labs(
        title = if (base::is.null(focus_label)) {
          overview_title
        } else {
          base::paste0(overview_title, " - focus on ", focus_label)
        },
        subtitle = network_subtitle
      )
    }
    p
  }

  heatmap_grob <- ggplot2::ggplotGrob(hc_heatmap_plot)
  header_x_norm <- (c(x_module, x_enrichment, x_upstream) - x_limits[[1]]) / (x_limits[[2]] - x_limits[[1]])
  build_network_column_header_grob <- function(network_grob) {
    panel_layout_idx <- which(network_grob$layout$name == "panel")
    if (base::length(panel_layout_idx) == 0) {
      return(grid::nullGrob())
    }
    panel_layout_idx <- panel_layout_idx[[1]]
    panel_l <- network_grob$layout$l[[panel_layout_idx]]
    panel_r <- network_grob$layout$r[[panel_layout_idx]]
    header_tbl <- gtable::gtable(
      widths = network_grob$widths,
      heights = grid::unit(1, "null")
    )
    header_text <- grid::grobTree(
      grid::textGrob(
        label = "Modules",
        x = header_x_norm[[1]],
        y = 0.5,
        just = "center",
        gp = grid::gpar(fontsize = 10.5 * overall_plot_scale, fontface = "bold")
      ),
      grid::textGrob(
        label = "Enriched terms/pathways",
        x = header_x_norm[[2]],
        y = 0.5,
        just = "center",
        gp = grid::gpar(fontsize = 10.5 * overall_plot_scale, fontface = "bold")
      ),
      grid::textGrob(
        label = "Inferred TF/Pathway",
        x = header_x_norm[[3]],
        y = 0.5,
        just = "center",
        gp = grid::gpar(fontsize = 10.5 * overall_plot_scale, fontface = "bold")
      )
    )
    gtable::gtable_add_grob(
      x = header_tbl,
      grobs = header_text,
      t = 1,
      l = panel_l,
      r = panel_r,
      clip = "off",
      name = "network_column_headers"
    )
  }

  draw_combined_page <- function(network_plot_obj, page_title_text) {
    top_grob <- grid::textGrob(
      label = page_title_text,
      x = 0.5,
      y = 0.5,
      just = "center",
      gp = grid::gpar(fontsize = 13 * overall_plot_scale, fontface = "bold")
    )
    network_grob <- ggplot2::ggplotGrob(network_plot_obj)
    column_header_row <- gridExtra::arrangeGrob(
      grobs = list(grid::nullGrob(), build_network_column_header_grob(network_grob)),
      ncol = 2,
      widths = c(0.35, 0.65)
    )
    main_grob <- gridExtra::arrangeGrob(
      grobs = list(heatmap_grob, network_grob),
      ncol = 2,
      widths = c(0.35, 0.65)
    )
    gridExtra::grid.arrange(
      grobs = list(top_grob, column_header_row, main_grob),
      ncol = 1,
      heights = c(0.06, 0.04, 0.90)
    )
  }

  overview_plot <- build_network_plot(show_headers = FALSE, show_titles = FALSE)
  focus_plots <- stats::setNames(
    lapply(cluster_order, function(cl) {
      build_network_plot(focus_cluster = cl, show_headers = FALSE, show_titles = FALSE)
    }),
    base::as.character(module_label_map[cluster_order])
  )
  focus_titles <- stats::setNames(
    lapply(cluster_order, function(cl) page_title_for_focus(cl)),
    base::as.character(module_label_map[cluster_order])
  )

  max_term_count <- base::max(base::nrow(enrich_nodes), base::nrow(upstream_nodes), base::length(cluster_order))
  pdf_width_auto <- base::max(14, base::min(27, 12 + 0.08 * max_term_count)) * overall_plot_scale
  pdf_height_auto <- base::max(9, base::min(18, 8 + 0.11 * max_term_count)) * overall_plot_scale
  pdf_width_use <- if (base::is.null(pdf_width)) pdf_width_auto else as.numeric(pdf_width)
  pdf_height_use <- if (base::is.null(pdf_height)) pdf_height_auto else as.numeric(pdf_height)

  out_file <- NULL
  focus_files <- stats::setNames(base::character(0), base::character(0))
  if (isTRUE(save_pdf)) {
    out_file <- base::paste0(file_prefix, "/", pdf_name)
    if (requireNamespace("Cairo", quietly = TRUE)) {
      Cairo::CairoPDF(file = out_file, width = pdf_width_use, height = pdf_height_use, pointsize = pdf_pointsize)
    } else {
      grDevices::pdf(file = out_file, width = pdf_width_use, height = pdf_height_use, pointsize = pdf_pointsize)
    }
    draw_combined_page(overview_plot, page_title_for_focus())
    for (i in base::seq_along(focus_plots)) {
      draw_combined_page(focus_plots[[i]], focus_titles[[i]])
    }
    grDevices::dev.off()

    focus_dir <- base::paste0(file_prefix, "/Module_Knowledge_Network_by_module")
    if (!base::dir.exists(focus_dir)) {
      base::dir.create(focus_dir, recursive = TRUE, showWarnings = FALSE)
    }
    sanitize_file_stem <- function(x) {
      x <- base::as.character(x)
      x <- gsub("[^A-Za-z0-9_-]+", "_", x)
      x <- gsub("_+", "_", x)
      x <- gsub("^_|_$", "", x)
      if (base::nchar(x) == 0) {
        x <- "module"
      }
      x
    }
    focus_names <- base::names(focus_plots)
    focus_files <- stats::setNames(base::character(base::length(focus_names)), focus_names)
    for (i in base::seq_along(focus_plots)) {
      mod_nm <- focus_names[[i]]
      mod_file <- base::paste0(
        focus_dir,
        "/Module_Knowledge_Network_",
        sanitize_file_stem(mod_nm),
        ".pdf"
      )
      if (requireNamespace("Cairo", quietly = TRUE)) {
        Cairo::CairoPDF(file = mod_file, width = pdf_width_use, height = pdf_height_use, pointsize = pdf_pointsize)
      } else {
        grDevices::pdf(file = mod_file, width = pdf_width_use, height = pdf_height_use, pointsize = pdf_pointsize)
      }
      draw_combined_page(focus_plots[[i]], focus_titles[[i]])
      grDevices::dev.off()
      focus_files[[i]] <- mod_file
    }
  }

  if (isTRUE(show_plot)) {
    draw_combined_page(overview_plot, page_title_for_focus())
    for (i in base::seq_along(focus_plots)) {
      draw_combined_page(focus_plots[[i]], focus_titles[[i]])
    }
  }

  nodes <- base::rbind(
    module_nodes[, c("id", "key", "label", "node_type", "x", "y"), drop = FALSE],
    if (base::nrow(enrich_nodes) > 0) {
      enrich_nodes[, c("id", "key", "label", "node_type", "x", "y"), drop = FALSE]
    } else {
      base::data.frame()
    },
    if (base::nrow(upstream_nodes) > 0) {
      upstream_nodes[, c("id", "key", "label", "node_type", "x", "y"), drop = FALSE]
    } else {
      base::data.frame()
    }
  )

  output <- list(
    plot = overview_plot,
    focus_plots = focus_plots,
    focus_files = focus_files,
    nodes = nodes,
    edges = edges,
    file = out_file,
    settings = list(
      enrichment_mode = enrichment_mode,
      upstream_mode = upstream_mode,
      upstream_activity_input = upstream_activity_input,
      heatmap_value_name = hc_heatmap_data$value_name,
      heatmap_scale_limits = hc_heatmap_data$scale_limits,
      label_mode = label_mode,
      clusters = cluster_order,
      max_enrichment_per_module = max_enrichment_per_module,
      max_upstream_per_module = max_upstream_per_module
    )
  )
  hcobject[["integrated_output"]][["knowledge_network"]] <<- output
  hcobject[["satellite_outputs"]][["knowledge_network"]] <<- output
  output
}
