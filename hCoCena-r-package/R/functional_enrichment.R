#' Functional enrichment analysis
#' 
#' Returns a cluster-wise Gene Ontology Database enrichment.
#' @param gene_sets A vector. The names of databases enrichment should be performed for. Choose one or multiple of "Go", "Kegg", "Hallmark", and/or "Reactome".
#'  Available databases depend on supplement files previously set
#'  Default is "Hallmark".
#' @param top Integer. The number of most strongly enriched terms to return per cluster. Default is 5.
#' @param clusters Either "all" (default) or a vector of clusters as strings. Defines for which clusters to perform the enrichment.
#' @param padj Method to use for multiple testing correction. Can be one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  Default is "BH" (Benjamini-Hochberg).
#' @param qval Upper threshold for the adjusted p-value. Default is 0.05.
#' @param heatmap_side Position of the hCoCena heatmap in the combined output.
#'  Choose one of "left" (default) or "right".
#' @param heatmap_cluster_rows A Boolean whether or not to cluster rows in the hCoCena heatmap.
#' @param heatmap_cluster_columns A Boolean whether or not to cluster columns in the hCoCena heatmap.
#' @param heatmap_show_row_dend A Boolean whether to show the row dendrogram when `heatmap_cluster_rows = TRUE`.
#' @param heatmap_show_column_dend A Boolean whether to show the column dendrogram when `heatmap_cluster_columns = TRUE`.
#' @param heatmap_order Optional character vector specifying module order in the hCoCena heatmap.
#'  Entries can be module colors (e.g. "turquoise") or module labels from the main heatmap
#'  (e.g. "M1", "M2", ...). Modules not listed are appended afterwards.
#' @param heatmap_module_label_mode Controls labels inside module color boxes of the hCoCena heatmap.
#'  One of "same" or "none". "same" reuses the prefix from `plot_cluster_heatmap()`
#'  and reindexes modules consecutively (M1, M2, ...) after final module ordering.
#' @param heatmap_show_gene_counts A Boolean whether or not to show gene counts per module
#'  in the right annotation of the hCoCena heatmap. Default is FALSE.
#' @param enrichment_vertical_line_mode Controls vertical helper lines in the enrichment panel.
#'  One of "to_term" (default, current behavior), "full" (full-height lines), or "none".
#' @param enrichment_row_height_scale Numeric scaling factor for module/enrichment row height
#'  in the combined plot. Values < 1 make rows slimmer; values > 1 make rows taller.
#'  Default is 0.9.
#' @param enrichment_column_spacing_scale Numeric scaling factor for spacing between
#'  enrichment terms (x direction). Values < 1 reduce spacing and make the enrichment
#'  panel narrower; values > 1 increase spacing. Default is 1.
#' @param enrichment_line_width_scale Numeric scaling factor for line widths in the
#'  enrichment panel (module guide lines and vertical helper lines). Values < 1 make
#'  lines thinner; values > 1 make lines thicker. Default is 1.
#' @param heatmap_column_label_fontsize Optional numeric font size for hCoCena
#'  heatmap column labels in enrichment plots. If NULL (default), uses automatic sizing.
#' @param heatmap_module_label_fontsize Optional numeric font size for module labels
#'  inside module color boxes (`M1`, `M2`, ...). If NULL (default), uses automatic sizing.
#' @param legend_fontsize Optional numeric base font size for enrichment-related legends
#'  (GFC legend and enrichment significance legend). If NULL (default), uses automatic sizing.
#' @param enrichment_label_fontsize Optional numeric font size for enrichment term labels.
#'  If NULL (default), uses automatic sizing.
#' @param enrichment_label_wrap Logical. If TRUE, wraps enrichment term labels using
#'  `enrichment_label_wrap_width`. Default is FALSE.
#' @param enrichment_label_wrap_width Integer wrap width used when
#'  `enrichment_label_wrap = TRUE`. Default is 30.
#' @param gfc_colors Optional character vector of colors for the hCoCena heatmap
#'  GFC scale. If NULL, uses the legacy default
#'  `rev(RColorBrewer::brewer.pal(11, "RdBu"))`.
#' @param gfc_scale_limits Optional numeric vector controlling the module-heatmap
#'  color scale limits used in enrichment plots. Provide one positive number
#'  (`x` -> `c(-x, x)`) or two numbers (`c(min, max)`). If NULL, uses stored
#'  limits from the latest main module heatmap when available, otherwise
#'  falls back to `c(-range_GFC, range_GFC)`.
#' @param pdf_width Optional numeric width (inches) for enrichment PDFs.
#'  If NULL (default), width is auto-estimated from content.
#' @param pdf_height Optional numeric height (inches) for enrichment PDFs.
#'  If NULL (default), height is auto-estimated from content.
#' @param pdf_pointsize Numeric base pointsize used for PDF devices.
#'  Default is 11.
#' @param overall_plot_scale Numeric scaling factor for the entire combined output
#'  (hCoCena heatmap + enrichment panel + title/labels/legend). Values > 1 enlarge
#'  the full plot, values < 1 make it smaller. Default is 1.
#' @export

functional_enrichment <- function(gene_sets = "Hallmark",
                                  top = 5,
                                  clusters = c("all"),
                                  padj = "BH",
                                  qval = 0.05,
                                  heatmap_side = "left",
                                  heatmap_cluster_rows = FALSE,
                                  heatmap_cluster_columns = FALSE,
                                  heatmap_show_row_dend = FALSE,
                                  heatmap_show_column_dend = FALSE,
                                  heatmap_order = NULL,
                                  heatmap_module_label_mode = "same",
                                  heatmap_show_gene_counts = FALSE,
                                  enrichment_vertical_line_mode = "to_term",
                                  enrichment_row_height_scale = 0.9,
                                  enrichment_column_spacing_scale = 1,
                                  enrichment_line_width_scale = 1,
                                  heatmap_column_label_fontsize = NULL,
                                  heatmap_module_label_fontsize = NULL,
                                  legend_fontsize = NULL,
                                  enrichment_label_fontsize = NULL,
                                  enrichment_label_wrap = FALSE,
                                  enrichment_label_wrap_width = 30,
                                  gfc_colors = NULL,
                                  gfc_scale_limits = NULL,
                                  pdf_width = NULL,
                                  pdf_height = NULL,
                                  pdf_pointsize = 11,
                                  overall_plot_scale = 1) {
  .hc_legacy_warning("functional_enrichment")
  gfc_colors_was_missing <- missing(gfc_colors)

  heatmap_side <- base::match.arg(heatmap_side, choices = c("left", "right"))
  heatmap_module_label_mode <- base::match.arg(
    heatmap_module_label_mode,
    choices = c("same", "none")
  )
  enrichment_vertical_line_mode <- base::match.arg(
    enrichment_vertical_line_mode,
    choices = c("to_term", "full", "none")
  )
  if (!base::is.logical(heatmap_cluster_rows) || base::length(heatmap_cluster_rows) != 1) {
    stop("`heatmap_cluster_rows` must be TRUE or FALSE.")
  }
  if (!base::is.logical(heatmap_cluster_columns) || base::length(heatmap_cluster_columns) != 1) {
    stop("`heatmap_cluster_columns` must be TRUE or FALSE.")
  }
  if (!base::is.logical(heatmap_show_row_dend) || base::length(heatmap_show_row_dend) != 1) {
    stop("`heatmap_show_row_dend` must be TRUE or FALSE.")
  }
  if (!base::is.logical(heatmap_show_column_dend) || base::length(heatmap_show_column_dend) != 1) {
    stop("`heatmap_show_column_dend` must be TRUE or FALSE.")
  }
  if (!base::is.logical(heatmap_show_gene_counts) || base::length(heatmap_show_gene_counts) != 1) {
    stop("`heatmap_show_gene_counts` must be TRUE or FALSE.")
  }
  if (!base::is.numeric(enrichment_row_height_scale) ||
      base::length(enrichment_row_height_scale) != 1 ||
      base::is.na(enrichment_row_height_scale) ||
      enrichment_row_height_scale <= 0) {
    stop("`enrichment_row_height_scale` must be a positive numeric scalar.")
  }
  if (!base::is.numeric(enrichment_column_spacing_scale) ||
      base::length(enrichment_column_spacing_scale) != 1 ||
      base::is.na(enrichment_column_spacing_scale) ||
      enrichment_column_spacing_scale <= 0) {
    stop("`enrichment_column_spacing_scale` must be a positive numeric scalar.")
  }
  if (!base::is.numeric(enrichment_line_width_scale) ||
      base::length(enrichment_line_width_scale) != 1 ||
      base::is.na(enrichment_line_width_scale) ||
      enrichment_line_width_scale <= 0) {
    stop("`enrichment_line_width_scale` must be a positive numeric scalar.")
  }
  if (!base::is.numeric(overall_plot_scale) ||
      base::length(overall_plot_scale) != 1 ||
      base::is.na(overall_plot_scale) ||
      overall_plot_scale <= 0) {
    stop("`overall_plot_scale` must be a positive numeric scalar.")
  }
  if (!base::is.null(heatmap_order) && !base::is.character(heatmap_order)) {
    stop("`heatmap_order` must be NULL or a character vector.")
  }
  overall_plot_scale <- base::max(0.5, base::min(3, overall_plot_scale))
  if (!base::is.null(gfc_colors)) {
    if (!base::is.character(gfc_colors) || base::length(gfc_colors) < 2) {
      stop("`gfc_colors` must be NULL or a character vector with at least two colors.")
    }
    if (any(base::is.na(gfc_colors)) || any(gfc_colors == "")) {
      stop("`gfc_colors` must not contain NA or empty strings.")
    }
    gfc_colors <- base::as.character(gfc_colors)
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
  if (!base::is.null(heatmap_column_label_fontsize) &&
      (!base::is.numeric(heatmap_column_label_fontsize) ||
       base::length(heatmap_column_label_fontsize) != 1 ||
       base::is.na(heatmap_column_label_fontsize) ||
       heatmap_column_label_fontsize <= 0)) {
    stop("`heatmap_column_label_fontsize` must be NULL or a single positive number.")
  }
  if (!base::is.null(heatmap_module_label_fontsize) &&
      (!base::is.numeric(heatmap_module_label_fontsize) ||
       base::length(heatmap_module_label_fontsize) != 1 ||
       base::is.na(heatmap_module_label_fontsize) ||
       heatmap_module_label_fontsize <= 0)) {
    stop("`heatmap_module_label_fontsize` must be NULL or a single positive number.")
  }
  if (!base::is.null(legend_fontsize) &&
      (!base::is.numeric(legend_fontsize) ||
       base::length(legend_fontsize) != 1 ||
       base::is.na(legend_fontsize) ||
       legend_fontsize <= 0)) {
    stop("`legend_fontsize` must be NULL or a single positive number.")
  }
  if (!base::is.null(enrichment_label_fontsize) &&
      (!base::is.numeric(enrichment_label_fontsize) ||
       base::length(enrichment_label_fontsize) != 1 ||
       base::is.na(enrichment_label_fontsize) ||
       enrichment_label_fontsize <= 0)) {
    stop("`enrichment_label_fontsize` must be NULL or a single positive number.")
  }
  if (!base::is.logical(enrichment_label_wrap) ||
      base::length(enrichment_label_wrap) != 1 ||
      base::is.na(enrichment_label_wrap)) {
    stop("`enrichment_label_wrap` must be TRUE or FALSE.")
  }
  if (!base::is.numeric(enrichment_label_wrap_width) ||
      base::length(enrichment_label_wrap_width) != 1 ||
      base::is.na(enrichment_label_wrap_width) ||
      enrichment_label_wrap_width <= 0) {
    stop("`enrichment_label_wrap_width` must be a single positive number.")
  }
  enrichment_label_wrap_width <- base::as.integer(base::round(enrichment_label_wrap_width))
  if (enrichment_label_wrap_width < 5) {
    stop("`enrichment_label_wrap_width` must be >= 5.")
  }

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

  resolve_gfc_scale_limits <- function(input_limits) {
    lims <- normalize_scale_limits(input_limits)
    if (!base::is.null(lims)) {
      return(lims)
    }
    stored <- tryCatch(
      normalize_scale_limits(hcobject[["integrated_output"]][["cluster_calc"]][["gfc_scale_limits"]]),
      error = function(e) NULL
    )
    if (!base::is.null(stored)) {
      return(stored)
    }
    fallback_lim <- suppressWarnings(base::as.numeric(hcobject[["global_settings"]][["range_GFC"]]))
    if (!base::is.finite(fallback_lim) || fallback_lim <= 0) {
      fallback_lim <- 2
    }
    c(-base::abs(fallback_lim), base::abs(fallback_lim))
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

  gfc_scale_limits <- resolve_gfc_scale_limits(gfc_scale_limits)
  gfc_scale_ticks <- compute_scale_ticks(gfc_scale_limits)
  pdf_width_input <- pdf_width
  pdf_height_input <- pdf_height

  # Publication-oriented typography presets
  font_title <- 16 * overall_plot_scale
  font_axis <- 10 * overall_plot_scale
  font_annotation <- 11 * overall_plot_scale
  font_module <- 8.2 * overall_plot_scale
  # Keep enrichment panels closer to the visual footprint of the main heatmap.
  enrichment_panel_compact_scale <- 0.82

  estimate_rotated_label_height_cm <- function(labels,
                                                font_size,
                                                min_cm = 4.8,
                                                max_cm = 11.2,
                                                extra_cm = 0.12) {
    if (base::is.null(labels) || base::length(labels) == 0) {
      return(min_cm)
    }
    labels <- base::as.character(labels)
    labels <- labels[!base::is.na(labels) & labels != ""]
    if (base::length(labels) == 0) {
      return(min_cm)
    }
    label_width_cm <- base::vapply(labels, function(lbl) {
      grid::convertWidth(
        grid::grobWidth(
          grid::textGrob(lbl, gp = grid::gpar(fontsize = font_size))
        ),
        "cm",
        valueOnly = TRUE
      )
    }, FUN.VALUE = base::numeric(1))
    est_cm <- base::max(label_width_cm, na.rm = TRUE) + extra_cm
    base::max(min_cm, base::min(max_cm, est_cm))
  }

  compute_term_label_fontsize <- function(max_chars, base_font = font_axis) {
    if (!base::is.finite(max_chars) || max_chars <= 0) {
      return(base_font)
    }
    adaptive <- 280 / max_chars
    base::max(7.0, base::min(base_font, adaptive))
  }

  maybe_wrap_labels <- function(labels) {
    labels <- base::as.character(labels)
    if (!isTRUE(enrichment_label_wrap) || base::length(labels) == 0) {
      return(labels)
    }
    base::vapply(labels, function(lbl) {
      if (base::is.na(lbl) || lbl == "") {
        return(lbl)
      }
      wrapped <- base::strwrap(lbl, width = enrichment_label_wrap_width, simplify = FALSE)
      if (base::length(wrapped) == 0) {
        return(lbl)
      }
      base::paste(wrapped[[1]], collapse = "\n")
    }, FUN.VALUE = base::character(1))
  }

  visual_max_label_chars <- function(labels) {
    if (base::length(labels) == 0) {
      return(10)
    }
    labels <- base::as.character(labels)
    labels <- labels[!base::is.na(labels) & labels != ""]
    if (base::length(labels) == 0) {
      return(10)
    }
    base::max(
      base::vapply(
        base::strsplit(labels, "\n", fixed = TRUE),
        function(parts) {
          if (base::length(parts) == 0) {
            return(0)
          }
          base::max(base::nchar(parts), na.rm = TRUE)
        },
        FUN.VALUE = base::numeric(1)
      ),
      na.rm = TRUE
    )
  }

  q_to_sig_score <- function(q) {
    if (base::is.na(q) || q <= 0) {
      return(0)
    }
    base::min(1, base::max(0, (-base::log10(q) - 1) / 4))
  }
  interleaved_mixed_term_order <- function(plot_df, cluster_levels) {
    if (base::is.null(plot_df) || base::nrow(plot_df) == 0) {
      return(base::character(0))
    }
    needed_cols <- c("term_with_db", "database", "cluster", "qvalue", "rank")
    if (!all(needed_cols %in% base::colnames(plot_df))) {
      return(base::unique(base::as.character(plot_df$term_with_db)))
    }

    df <- plot_df
    df$term_with_db <- base::as.character(df$term_with_db)
    df$database <- base::as.character(df$database)
    df$cluster <- base::as.character(df$cluster)
    df$qvalue <- suppressWarnings(base::as.numeric(df$qvalue))
    df$rank <- suppressWarnings(base::as.numeric(df$rank))

    cluster_idx_map <- stats::setNames(base::seq_along(cluster_levels), cluster_levels)
    df$cluster_idx <- cluster_idx_map[df$cluster]
    df$cluster_idx[base::is.na(df$cluster_idx)] <- base::length(cluster_levels) + 1

    idx_by_term <- base::split(base::seq_len(base::nrow(df)), df$term_with_db)
    term_meta_list <- base::lapply(base::names(idx_by_term), function(term_nm) {
      idx <- idx_by_term[[term_nm]]
      sub <- df[idx, , drop = FALSE]
      qv <- sub$qvalue[!base::is.na(sub$qvalue)]
      rk <- sub$rank[!base::is.na(sub$rank)]
      base::data.frame(
        term_with_db = term_nm,
        database = sub$database[1],
        first_cluster_idx = base::min(sub$cluster_idx, na.rm = TRUE),
        mean_cluster_idx = base::mean(sub$cluster_idx, na.rm = TRUE),
        hit_count = base::length(base::unique(sub$cluster)),
        best_q = if (base::length(qv) == 0) {
          Inf
        } else {
          base::min(qv)
        },
        best_rank = if (base::length(rk) == 0) {
          Inf
        } else {
          base::min(rk)
        },
        stringsAsFactors = FALSE
      )
    })
    term_meta <- base::do.call(base::rbind, term_meta_list)
    if (base::is.null(term_meta) || base::nrow(term_meta) == 0) {
      return(base::character(0))
    }
    term_meta$first_cluster_idx[!is.finite(term_meta$first_cluster_idx)] <- base::length(cluster_levels) + 1
    term_meta$mean_cluster_idx[!is.finite(term_meta$mean_cluster_idx)] <- base::length(cluster_levels) + 1
    term_meta$best_q[!is.finite(term_meta$best_q)] <- Inf
    term_meta$best_rank[!is.finite(term_meta$best_rank)] <- Inf

    db_order <- base::unique(df$database)
    cluster_buckets <- base::sort(base::unique(term_meta$first_cluster_idx))
    ordered_terms <- base::character(0)

    for (bucket_idx in cluster_buckets) {
      bucket_meta <- term_meta[term_meta$first_cluster_idx == bucket_idx, , drop = FALSE]
      if (base::nrow(bucket_meta) == 0) {
        next
      }
      db_lists <- base::lapply(db_order, function(db_nm) {
        x <- bucket_meta[bucket_meta$database == db_nm, , drop = FALSE]
        if (base::nrow(x) == 0) {
          return(base::character(0))
        }
        x <- x[base::order(
          x$best_rank,
          x$best_q,
          -x$hit_count,
          x$mean_cluster_idx,
          x$term_with_db
        ), , drop = FALSE]
        x$term_with_db
      })
      base::names(db_lists) <- db_order

      max_len <- base::max(base::lengths(db_lists))
      for (k in base::seq_len(max_len)) {
        for (db_nm in db_order) {
          cur_terms <- db_lists[[db_nm]]
          if (base::length(cur_terms) >= k) {
            ordered_terms <- base::c(ordered_terms, cur_terms[k])
          }
        }
      }
    }

    ordered_terms <- base::unique(ordered_terms)
    missing_terms <- base::setdiff(term_meta$term_with_db, ordered_terms)
    if (base::length(missing_terms) > 0) {
      missing_meta <- term_meta[base::match(missing_terms, term_meta$term_with_db), , drop = FALSE]
      missing_meta <- missing_meta[base::order(
        missing_meta$best_rank,
        missing_meta$best_q,
        -missing_meta$hit_count,
        missing_meta$mean_cluster_idx,
        missing_meta$term_with_db
      ), , drop = FALSE]
      ordered_terms <- base::c(ordered_terms, missing_meta$term_with_db)
    }
    ordered_terms
  }

  # Define general stats
  universe <- base::lapply(1:base::length(hcobject[["layers"]]),
                           function(x) {
                             return(base::rownames(hcobject[["data"]][[base::paste0("set", x, "_counts")]]))
                           }) %>% base::unlist() %>% base::unique()

  cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]

  all_clusters <- base::unique(cluster_info$color)
  all_clusters <- all_clusters[all_clusters != "white"]

  if (clusters[1] == "all") {
    clusters <- all_clusters
  }
  clusters <- base::as.character(clusters)
  missing_clusters <- base::setdiff(clusters, all_clusters)
  if (base::length(missing_clusters) > 0) {
    warning(
      "Ignoring unknown cluster(s): ",
      base::paste(missing_clusters, collapse = ", ")
    )
    clusters <- clusters[clusters %in% all_clusters]
  }
  if (base::length(clusters) == 0) {
    stop("No valid clusters available for enrichment.")
  }

  # Build a robust hCoCena heatmap matrix (rows = modules, columns = conditions).
  cluster_calc <- hcobject[["integrated_output"]][["cluster_calc"]]
  stored_hm <- cluster_calc[["heatmap_cluster"]]
  m <- tryCatch(
    {
      stored_hm@ht_list[[1]]@matrix
    },
    error = function(e) NULL
  )

  if (base::is.null(m)) {
    gfc_all <- hcobject[["integrated_output"]][["GFC_all_layers"]]
    if (base::is.null(gfc_all) || base::nrow(gfc_all) == 0 || base::ncol(gfc_all) < 2) {
      stop("Unable to build hCoCena heatmap matrix: missing `GFC_all_layers`.")
    }
    value_cols <- base::colnames(gfc_all)[1:(base::ncol(gfc_all) - 1)]
    m_list <- list()
    for (c in all_clusters) {
      genes <- dplyr::filter(cluster_info, color == c) %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(split = ",") %>%
        base::unlist()
      if (base::length(genes) == 0) {
        next
      }
      c_gfcs <- dplyr::filter(gfc_all, Gene %in% genes)
      if (base::nrow(c_gfcs) == 0) {
        next
      }
      c_means <- base::apply(
        c_gfcs[, value_cols, drop = FALSE] %>% base::as.data.frame(),
        2,
        base::mean
      )
      m_list[[c]] <- c_means
    }
    if (base::length(m_list) == 0) {
      stop("Unable to build hCoCena heatmap matrix: no module GFC means available.")
    }
    m <- base::do.call(base::rbind, m_list)
    m <- m %>% base::as.matrix()
  }

  module_colors <- base::rownames(m)
  if (base::is.null(module_colors) || base::length(module_colors) == 0) {
    stop("Unable to build hCoCena heatmap: no module rows available.")
  }

  # Keep condition order from the previously drawn main heatmap when available.
  if (!base::is.null(stored_hm)) {
    col_order <- tryCatch(ComplexHeatmap::column_order(stored_hm), error = function(e) NULL)
    if (base::is.list(col_order) && base::length(col_order) > 0) {
      col_order <- col_order[[1]]
    }
    if (base::is.numeric(col_order) && base::length(col_order) == base::ncol(m)) {
      m <- m[, col_order, drop = FALSE]
    } else if (base::is.character(col_order)) {
      keep_cols <- col_order[col_order %in% base::colnames(m)]
      if (base::length(keep_cols) > 0) {
        m <- m[, keep_cols, drop = FALSE]
      }
    }
  }

  module_order <- module_colors
  if (!base::is.null(stored_hm)) {
    row_order <- tryCatch(ComplexHeatmap::row_order(stored_hm), error = function(e) NULL)
    if (base::is.list(row_order) && base::length(row_order) > 0) {
      row_order <- row_order[[1]]
    }
    if (base::is.numeric(row_order) && base::length(row_order) == base::length(module_colors)) {
      module_order <- module_colors[row_order]
    } else if (base::is.character(row_order)) {
      keep_rows <- row_order[row_order %in% module_colors]
      if (base::length(keep_rows) > 0) {
        module_order <- base::c(keep_rows, base::setdiff(module_colors, keep_rows))
      }
    }
  }

  module_label_map <- cluster_calc[["module_label_map"]]
  module_prefix <- cluster_calc[["module_prefix"]]
  if (!base::is.null(module_label_map) && base::length(module_label_map) > 0) {
    map_names <- base::names(cluster_calc[["module_label_map"]])
    module_label_map <- base::as.character(module_label_map)
    if (!base::is.null(map_names) && base::length(map_names) == base::length(module_label_map)) {
      base::names(module_label_map) <- base::as.character(map_names)
    }
    missing_before <- base::setdiff(module_order, base::names(module_label_map))
    if (base::length(missing_before) > 0) {
      inverse_map <- stats::setNames(base::names(module_label_map), base::as.character(module_label_map))
      if (base::all(module_order %in% base::names(inverse_map))) {
        module_label_map <- inverse_map
      }
    }
  }
  stored_module_label_fontsize <- cluster_calc[["module_label_fontsize"]]
  if (!base::is.numeric(stored_module_label_fontsize) ||
      base::length(stored_module_label_fontsize) != 1 ||
      base::is.na(stored_module_label_fontsize) ||
      stored_module_label_fontsize <= 0) {
    stored_module_label_fontsize <- NULL
  }
  stored_module_label_pt_size <- cluster_calc[["module_label_pt_size"]]
  if (!base::is.numeric(stored_module_label_pt_size) ||
      base::length(stored_module_label_pt_size) != 1 ||
      base::is.na(stored_module_label_pt_size) ||
      stored_module_label_pt_size <= 0) {
    stored_module_label_pt_size <- NULL
  }
  stored_module_box_width_cm <- cluster_calc[["module_box_width_cm"]]
  if (!base::is.numeric(stored_module_box_width_cm) ||
      base::length(stored_module_box_width_cm) != 1 ||
      base::is.na(stored_module_box_width_cm) ||
      stored_module_box_width_cm <= 0) {
    stored_module_box_width_cm <- NULL
  }
  stored_module_box_to_cell_ratio <- cluster_calc[["module_box_to_cell_ratio"]]
  if (!base::is.numeric(stored_module_box_to_cell_ratio) ||
      base::length(stored_module_box_to_cell_ratio) != 1 ||
      base::is.na(stored_module_box_to_cell_ratio) ||
      stored_module_box_to_cell_ratio <= 0) {
    stored_module_box_to_cell_ratio <- NULL
  }
  stored_heatmap_cell_size_mm <- cluster_calc[["heatmap_cell_size_mm"]]
  if (!base::is.numeric(stored_heatmap_cell_size_mm) ||
      base::length(stored_heatmap_cell_size_mm) != 1 ||
      base::is.na(stored_heatmap_cell_size_mm) ||
      stored_heatmap_cell_size_mm <= 0) {
    stored_heatmap_cell_size_mm <- NULL
  }
  if (base::is.null(stored_module_box_to_cell_ratio) &&
      !base::is.null(stored_module_box_width_cm) &&
      !base::is.null(stored_heatmap_cell_size_mm)) {
    stored_module_box_to_cell_ratio <- (stored_module_box_width_cm * 10) / stored_heatmap_cell_size_mm
    if (!base::is.finite(stored_module_box_to_cell_ratio) ||
        stored_module_box_to_cell_ratio <= 0) {
      stored_module_box_to_cell_ratio <- NULL
    }
  }
  stored_gene_count_fontsize <- cluster_calc[["gene_count_fontsize"]]
  if (!base::is.numeric(stored_gene_count_fontsize) ||
      base::length(stored_gene_count_fontsize) != 1 ||
      base::is.na(stored_gene_count_fontsize) ||
      stored_gene_count_fontsize <= 0) {
    stored_gene_count_fontsize <- NULL
  }
  stored_gene_count_renderer <- cluster_calc[["gene_count_renderer"]]
  if (base::is.null(stored_gene_count_renderer) ||
      !base::is.character(stored_gene_count_renderer) ||
      base::length(stored_gene_count_renderer) != 1 ||
      !(stored_gene_count_renderer %in% c("pch", "text"))) {
    stored_gene_count_renderer <- "pch"
  }
  stored_gene_count_pt_size <- cluster_calc[["gene_count_pt_size"]]
  if (!base::is.numeric(stored_gene_count_pt_size) ||
      base::length(stored_gene_count_pt_size) != 1 ||
      base::is.na(stored_gene_count_pt_size) ||
      stored_gene_count_pt_size <= 0) {
    stored_gene_count_pt_size <- NULL
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
  if ((gfc_colors_was_missing || base::is.null(gfc_colors)) && !base::is.null(stored_gfc_colors)) {
    gfc_colors <- stored_gfc_colors
  }
  if (base::is.null(gfc_colors)) {
    gfc_colors <- base::rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  }
  if (!base::is.character(gfc_colors) || base::length(gfc_colors) < 2) {
    stop("`gfc_colors` must be NULL or a character vector with at least two colors.")
  }
  if (any(base::is.na(gfc_colors)) || any(gfc_colors == "")) {
    stop("`gfc_colors` must not contain NA or empty strings.")
  }
  gfc_colors <- base::as.character(gfc_colors)
  if (base::is.null(module_prefix) || !base::is.character(module_prefix) || base::length(module_prefix) != 1) {
    module_prefix <- "M"
  }
  if (base::is.null(module_label_map) || base::length(module_label_map) == 0) {
    module_label_map <- stats::setNames(
      base::paste0(module_prefix, base::seq_along(module_order)),
      module_order
    )
  }
  missing_map <- base::setdiff(module_order, base::names(module_label_map))
  if (base::length(missing_map) > 0) {
    start_idx <- base::length(module_label_map) + 1
    module_label_map <- base::c(
      module_label_map,
      stats::setNames(
        base::paste0(module_prefix, base::seq.int(start_idx, length.out = base::length(missing_map))),
        missing_map
      )
    )
  }

  if (!base::is.null(heatmap_order)) {
    inverse_map <- stats::setNames(base::names(module_label_map), base::as.character(module_label_map))
    mapped_order <- base::vapply(
      heatmap_order,
      function(x) {
        if (x %in% module_order) {
          return(x)
        }
        if (x %in% base::names(inverse_map)) {
          return(inverse_map[[x]])
        }
        return(NA_character_)
      },
      FUN.VALUE = base::character(1)
    )
    unknown_entries <- base::unique(heatmap_order[base::is.na(mapped_order)])
    if (base::length(unknown_entries) > 0) {
      stop(
        "Unknown entries in `heatmap_order`: ",
        base::paste(unknown_entries, collapse = ", ")
      )
    }
    mapped_order <- base::unique(mapped_order)
    module_order <- base::c(mapped_order, base::setdiff(module_order, mapped_order))
    heatmap_cluster_rows <- FALSE
  }

  cluster_order <- module_order[module_order %in% clusters]
  if (base::length(cluster_order) == 0) {
    stop("No valid module rows available after applying cluster selection/order.")
  }
  heatmap_mat <- m[cluster_order, , drop = FALSE]

  cluster_counts_df <- dplyr::filter(cluster_info, color %in% cluster_order) %>%
    dplyr::select(color, gene_no) %>%
    dplyr::distinct(color, .keep_all = TRUE)
  gene_count_map <- stats::setNames(cluster_counts_df$gene_no, cluster_counts_df$color)
  gene_counts <- gene_count_map[cluster_order]
  gene_counts[base::is.na(gene_counts)] <- 0
  module_labels <- base::as.character(module_label_map[cluster_order])
  missing_label_mask <- base::is.na(module_labels) | module_labels == ""
  if (base::any(missing_label_mask)) {
    module_labels[missing_label_mask] <- base::paste0(module_prefix, base::which(missing_label_mask))
  }
  module_label_map_current <- stats::setNames(module_labels, cluster_order)

  heatmap_row_names <- if (heatmap_module_label_mode == "same") {
    module_labels
  } else {
    cluster_order
  }
  base::rownames(heatmap_mat) <- heatmap_row_names

  n_modules <- base::length(cluster_order)
  label_fontsize <- if (n_modules <= 8) {
    font_module + 1.0
  } else if (n_modules <= 14) {
    font_module + 0.4
  } else if (n_modules <= 20) {
    font_module - 0.2
  } else if (n_modules <= 30) {
    font_module - 0.9
  } else {
    font_module - 1.4
  }
  max_label_chars <- if (heatmap_module_label_mode == "same") {
    base::max(base::nchar(module_labels), na.rm = TRUE)
  } else {
    1
  }
  label_fontsize <- base::max(
    5.5,
    base::min(
      10.5,
      base::min(label_fontsize, 17 / base::max(1.5, max_label_chars))
    )
  )
  if (heatmap_module_label_mode == "same" && !base::is.null(stored_module_label_fontsize)) {
    label_fontsize <- stored_module_label_fontsize
  }
  if (!base::is.null(heatmap_module_label_fontsize)) {
    label_fontsize <- base::as.numeric(heatmap_module_label_fontsize)
  }

  # Match module box size to the hCoCena heatmap cell geometry used in this
  # enrichment figure, while preserving proportions from the main heatmap.
  n_hc_rows <- base::nrow(heatmap_mat)
  n_hc_cols <- base::ncol(heatmap_mat)
  if (heatmap_module_label_mode == "same" && !base::is.null(stored_heatmap_cell_size_mm)) {
    # Match the main hCoCena heatmap geometry exactly when label mode is "same".
    hc_cell_mm <- stored_heatmap_cell_size_mm
  } else {
    hc_cell_mm_base <- if (n_hc_rows <= 10) {
      7.0
    } else if (n_hc_rows <= 18) {
      6.0
    } else if (n_hc_rows <= 30) {
      5.0
    } else {
      4.2
    }
    hc_cell_mm <- hc_cell_mm_base * enrichment_row_height_scale
    hc_cell_mm <- base::max(3.6, base::min(10.5, hc_cell_mm))
    if (n_hc_cols > 8) {
      hc_cell_mm <- base::min(hc_cell_mm, 5.4)
    }
    if (n_hc_cols > 12) {
      hc_cell_mm <- base::min(hc_cell_mm, 4.6)
    }
  }
  hc_cell_mm <- hc_cell_mm * overall_plot_scale

  use_stored_module_box_width <- (
    heatmap_module_label_mode == "same" &&
      !base::is.null(stored_module_box_width_cm)
  )
  module_box_width_cm <- if (heatmap_module_label_mode == "same") {
    if (use_stored_module_box_width) {
      stored_module_box_width_cm
    } else {
      base::max(
        0.65,
        base::min(
          3.2,
          0.34 + (0.14 * max_label_chars) + (0.035 * label_fontsize)
        )
      )
    }
  } else {
    0.35
  }
  if (!use_stored_module_box_width) {
    module_box_width_cm <- module_box_width_cm * overall_plot_scale
  }
  module_label_pt_size <- if (heatmap_module_label_mode == "same" && !base::is.null(stored_module_label_pt_size)) {
    stored_module_label_pt_size
  } else {
    base::max(
      0.16,
      base::min(
        0.72,
        0.06 * label_fontsize
      )
    )
  }
  if (!base::is.null(heatmap_module_label_fontsize)) {
    module_label_pt_size <- base::max(
      0.16,
      base::min(
        0.72,
        0.06 * label_fontsize
      )
    )
  }

  module_color_map <- stats::setNames(cluster_order, cluster_order)
  module_box_anno <- ComplexHeatmap::anno_simple(
    cluster_order,
    col = module_color_map,
    pch = if (heatmap_module_label_mode == "same") module_labels else NULL,
    pt_gp = grid::gpar(col = "white", fontsize = label_fontsize, fontface = "bold"),
    pt_size = grid::unit(module_label_pt_size, "snpc"),
    simple_anno_size = grid::unit(module_box_width_cm, "cm"),
    gp = grid::gpar(col = "black"),
    which = "row"
  )

  gene_count_fontsize <- if (!base::is.null(stored_gene_count_fontsize)) {
    stored_gene_count_fontsize
  } else {
    font_annotation
  }
  gene_count_pt_size <- if (!base::is.null(stored_gene_count_pt_size)) {
    stored_gene_count_pt_size
  } else if (!base::is.null(stored_module_label_pt_size)) {
    stored_module_label_pt_size
  } else {
    0.5
  }
  gene_count_anno <- if (heatmap_show_gene_counts) {
    if (stored_gene_count_renderer == "pch") {
      ComplexHeatmap::anno_simple(
        x = base::rep("count_text", base::length(gene_counts)),
        col = c(count_text = "transparent"),
        pch = base::as.character(gene_counts),
        pt_gp = grid::gpar(col = "black", fontsize = gene_count_fontsize, fontface = "plain"),
        pt_size = grid::unit(gene_count_pt_size, "snpc"),
        gp = grid::gpar(col = NA),
        simple_anno_size = grid::unit(1.2 * overall_plot_scale, "cm"),
        which = "row"
      )
    } else {
      ComplexHeatmap::anno_text(
        x = base::as.character(gene_counts),
        width = grid::unit(1.2 * overall_plot_scale, "cm"),
        just = "left",
        gp = grid::gpar(fontsize = gene_count_fontsize),
        which = "row"
      )
    }
  } else {
    ComplexHeatmap::anno_empty(width = grid::unit(0, "mm"), which = "row", border = FALSE)
  }

  right_anno <- ComplexHeatmap::HeatmapAnnotation(
    modules = module_box_anno,
    `# genes` = gene_count_anno,
    which = "row",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    annotation_name_side = "top",
    annotation_name_gp = grid::gpar(fontsize = font_annotation),
    gap = grid::unit(2 * overall_plot_scale, "mm")
  )

  # Keep module-expression tiles square-like in the combined enrichment figure.
  hc_body_w_mm <- base::max(18, n_hc_cols * hc_cell_mm)
  hc_body_h_mm <- base::max(20, n_hc_rows * hc_cell_mm)
  heatmap_column_fontsize <- if (base::is.null(heatmap_column_label_fontsize)) {
    font_axis
  } else {
    base::as.numeric(heatmap_column_label_fontsize)
  }
  gfc_palette <- grDevices::colorRampPalette(gfc_colors)(51)
  gfc_col_fun <- circlize::colorRamp2(
    seq(gfc_scale_limits[1], gfc_scale_limits[2], length.out = base::length(gfc_palette)),
    gfc_palette
  )

  gfc_legend_param <- list(
    title = "GFC",
    at = gfc_scale_ticks$breaks,
    labels = gfc_scale_ticks$labels
  )
  if (!base::is.null(legend_fontsize)) {
    gfc_legend_param$labels_gp <- grid::gpar(fontsize = base::as.numeric(legend_fontsize))
    gfc_legend_param$title_gp <- grid::gpar(
      fontsize = base::as.numeric(legend_fontsize) + 0.8,
      fontface = "bold"
    )
  }

  hc_ht <- ComplexHeatmap::Heatmap(
    heatmap_mat,
    name = "GFC",
    right_annotation = right_anno,
    col = gfc_col_fun,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    cluster_rows = heatmap_cluster_rows,
    cluster_columns = heatmap_cluster_columns,
    show_row_dend = heatmap_cluster_rows && heatmap_show_row_dend,
    show_column_dend = heatmap_cluster_columns && heatmap_show_column_dend,
    column_names_rot = 90,
    show_row_names = FALSE,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = font_axis),
    column_names_gp = grid::gpar(fontsize = heatmap_column_fontsize),
    width = grid::unit(hc_body_w_mm, "mm"),
    height = grid::unit(hc_body_h_mm, "mm"),
    rect_gp = grid::gpar(col = "black"),
    show_heatmap_legend = TRUE,
    heatmap_legend_param = gfc_legend_param
  )

  build_enrichment_export <- function(res_list, selected_summary_tbl, significant_summary_tbl = NULL) {
    export_tables <- list(selected_enrichments_summary = selected_summary_tbl)
    if (!base::is.null(significant_summary_tbl)) {
      export_tables <- base::c(
        export_tables,
        list(significant_enrichments_summary = significant_summary_tbl)
      )
    }
    if (base::length(res_list) > 0) {
      export_tables <- base::c(export_tables, res_list)
    } else {
      export_tables <- base::c(export_tables, list(no_enrichment = base::data.frame()))
    }
    export_tables
  }
  enrichment_summary_columns <- c(
    "database",
    "cluster",
    "module_label",
    "rank",
    "term",
    "qvalue",
    "p.adjust",
    "pvalue",
    "GeneRatio",
    "BgRatio",
    "Count",
    "geneID"
  )
  empty_enrichment_summary <- function() {
    base::data.frame(
      database = base::character(0),
      cluster = base::character(0),
      module_label = base::character(0),
      rank = base::integer(0),
      term = base::character(0),
      qvalue = base::numeric(0),
      p.adjust = base::numeric(0),
      pvalue = base::numeric(0),
      GeneRatio = base::character(0),
      BgRatio = base::character(0),
      Count = base::numeric(0),
      geneID = base::character(0),
      stringsAsFactors = FALSE
    )
  }
  normalize_enrichment_summary <- function(rows_list, database_name, cluster_levels) {
    if (base::length(rows_list) == 0) {
      return(empty_enrichment_summary())
    }
    ordered_clusters <- cluster_levels[cluster_levels %in% base::names(rows_list)]
    if (base::length(ordered_clusters) == 0) {
      return(empty_enrichment_summary())
    }
    rows_list <- rows_list[ordered_clusters]
    out <- base::do.call(base::rbind, rows_list)
    out <- out[, base::intersect(
      c(
        "cluster",
        "module_label",
        "rank",
        "Description",
        "qvalue",
        "p.adjust",
        "pvalue",
        "GeneRatio",
        "BgRatio",
        "Count",
        "geneID"
      ),
      base::colnames(out)
    ), drop = FALSE]
    if ("Description" %in% base::colnames(out)) {
      base::colnames(out)[base::colnames(out) == "Description"] <- "term"
    }
    out$database <- database_name
    if ("cluster" %in% base::colnames(out)) {
      out$cluster <- base::factor(out$cluster, levels = cluster_levels)
      out <- out[base::order(out$cluster, out$rank), , drop = FALSE]
      out$cluster <- base::as.character(out$cluster)
    }
    missing_cols <- base::setdiff(enrichment_summary_columns, base::colnames(out))
    if (base::length(missing_cols) > 0) {
      for (mc in missing_cols) {
        out[[mc]] <- NA
      }
    }
    out <- out[, enrichment_summary_columns, drop = FALSE]
    base::rownames(out) <- NULL
    out
  }
  combine_summary_tables <- function(summary_list) {
    combined <- empty_enrichment_summary()
    if (base::length(summary_list) > 0) {
      non_empty <- base::lapply(summary_list, function(x) {
        if (base::is.null(x) || base::nrow(x) == 0) {
          return(NULL)
        }
        x
      })
      non_empty <- non_empty[!base::vapply(non_empty, base::is.null, FUN.VALUE = base::logical(1))]
      if (base::length(non_empty) > 0) {
        combined <- base::do.call(base::rbind, non_empty)
        combined$cluster <- base::as.character(combined$cluster)
        base::rownames(combined) <- NULL
      }
    }
    combined
  }
  selected_summary_all_dbs <- list()
  significant_summary_all_dbs <- list()

  # Perform the analysis for each of the selected gene sets
  for (i in gene_sets) {
    i_title <- stringr::str_to_title(i)
    if (i_title %in% base::names(hcobject[["supplementary_data"]])) {
      top_enr <- list()
      res <- list()
      selected_rows <- list()
      significant_rows <- list()

      # Perform the enrichment for each of the clusters
      for (c in clusters) {
        genes <- dplyr::filter(cluster_info, color == c) %>%
          dplyr::pull(., "gene_n") %>%
          base::strsplit(split = ",") %>%
          base::unlist()

        enrich <- clusterProfiler::enricher(
          gene = genes,
          TERM2GENE = hcobject[["supplementary_data"]][[i_title]],
          qvalueCutoff = qval,
          pAdjustMethod = padj,
          universe = universe
        )

        if (!base::is.null(enrich)) {
          tmp <- enrich@result
          tmp <- dplyr::filter(tmp, qvalue <= qval)
          if (base::nrow(tmp) == 0) {
            next
          } else {
            tmp <- tmp[base::order(tmp$qvalue, decreasing = FALSE), , drop = FALSE]
            tmp_significant <- tmp
            tmp_significant$cluster <- c
            tmp_significant$module_label <- module_label_map_current[[c]]
            tmp_significant$rank <- base::seq_len(base::nrow(tmp_significant))
            significant_rows[[c]] <- tmp_significant
            pick_n <- base::min(top, base::nrow(tmp))
            tmp_selected <- tmp[base::seq_len(pick_n), , drop = FALSE]
            top_enr[[c]] <- tmp_selected$Description
            res[[c]] <- enrich@result
            tmp_selected$cluster <- c
            tmp_selected$module_label <- module_label_map_current[[c]]
            tmp_selected$rank <- base::seq_len(base::nrow(tmp_selected))
            selected_rows[[c]] <- tmp_selected
          }
        }
      }

      plot_cluster_order <- cluster_order
      selected_summary <- normalize_enrichment_summary(
        rows_list = selected_rows,
        database_name = i_title,
        cluster_levels = plot_cluster_order
      )
      significant_summary <- normalize_enrichment_summary(
        rows_list = significant_rows,
        database_name = i_title,
        cluster_levels = plot_cluster_order
      )
      selected_summary_all_dbs[[i_title]] <- selected_summary
      significant_summary_all_dbs[[i_title]] <- significant_summary

      if (base::length(top_enr) == 0) {
        output <- list(
          p = NULL,
          result = base::data.frame(),
          enrichment = res,
          selected_enrichments = selected_summary,
          significant_enrichments = significant_summary
        )
        hcobject[["integrated_output"]][["enrichments"]][[base::paste0("top_", i)]] <<- output
        openxlsx::write.xlsx(
          x = build_enrichment_export(res, selected_summary, significant_summary),
          file = base::paste0(
            hcobject[["working_directory"]][["dir_output"]],
            hcobject[["global_settings"]][["save_folder"]],
            "/Enrichment_",
            i_title,
            ".xlsx"
          ),
          overwrite = TRUE
        )
        warning("No enriched terms found for database: ", i_title)
        next
      }

      max_terms <- base::max(base::lengths(top_enr))
      top_enr_matrix <- base::do.call(
        base::cbind,
        base::lapply(plot_cluster_order, function(cl) {
          x <- top_enr[[cl]]
          if (base::length(x) < max_terms) {
            x <- base::c(x, base::rep(NA_character_, max_terms - base::length(x)))
          }
          x
        })
      )
      if (base::is.null(base::dim(top_enr_matrix))) {
        top_enr_matrix <- base::matrix(top_enr_matrix, ncol = 1)
      }
      base::colnames(top_enr_matrix) <- plot_cluster_order
      top_enr_out <- top_enr_matrix %>% base::as.data.frame(stringsAsFactors = FALSE)

      ggplot_df <- base::do.call(
        base::rbind,
        base::lapply(plot_cluster_order, function(cl) {
          terms <- top_enr_matrix[, cl]
          terms <- terms[!base::is.na(terms)]
          if (base::length(terms) == 0) {
            return(NULL)
          }
          base::data.frame(
            terms = terms,
            cluster = cl,
            stringsAsFactors = FALSE
          )
        })
      )

      if (base::is.null(ggplot_df) || base::nrow(ggplot_df) == 0) {
        output <- list(
          p = NULL,
          result = top_enr_out,
          enrichment = res,
          selected_enrichments = selected_summary,
          significant_enrichments = significant_summary
        )
        hcobject[["integrated_output"]][["enrichments"]][[base::paste0("top_", i)]] <<- output
        openxlsx::write.xlsx(
          x = build_enrichment_export(res, selected_summary, significant_summary),
          file = base::paste0(
            hcobject[["working_directory"]][["dir_output"]],
            hcobject[["global_settings"]][["save_folder"]],
            "/Enrichment_",
            i_title,
            ".xlsx"
          ),
          overwrite = TRUE
        )
        next
      }

      ggplot_df$cluster <- base::factor(ggplot_df$cluster, levels = plot_cluster_order)
      ggplot_df <- ggplot_df[stats::complete.cases(ggplot_df), ] %>% dplyr::arrange(cluster)
      term_levels <- base::unique(base::as.character(ggplot_df$terms))
      # Apply the same smart ordering logic used in mixed all-DB plots.
      # For single DB plots, this still gives a module-aware, significance-aware term order.
      term_order_df <- selected_summary
      if (base::nrow(term_order_df) > 0 &&
          all(c("database", "cluster", "term", "qvalue", "rank") %in% base::colnames(term_order_df))) {
        term_order_df <- term_order_df[
          !base::is.na(term_order_df$term) & term_order_df$term != "",
          ,
          drop = FALSE
        ]
        term_order_df$term_with_db <- base::as.character(term_order_df$term)
        ordered_terms <- interleaved_mixed_term_order(
          plot_df = term_order_df,
          cluster_levels = cluster_order
        )
        ordered_terms <- ordered_terms[ordered_terms %in% term_levels]
        term_levels <- base::c(ordered_terms, base::setdiff(term_levels, ordered_terms))
      }
      term_levels_display <- maybe_wrap_labels(term_levels)
      max_term_chars <- visual_max_label_chars(term_levels_display)
      term_label_fontsize <- compute_term_label_fontsize(max_term_chars)
      if (!base::is.null(enrichment_label_fontsize)) {
        term_label_fontsize <- base::as.numeric(enrichment_label_fontsize)
      }
      enrichment_colname_max_cm <- estimate_rotated_label_height_cm(
        labels = term_levels_display,
        font_size = term_label_fontsize,
        min_cm = 4.8,
        max_cm = 11.0,
        extra_cm = 0.10
      )
      row_names_for_enrichment <- base::rownames(heatmap_mat)
      row_name_by_cluster <- stats::setNames(row_names_for_enrichment, cluster_order)
      row_color_by_name <- stats::setNames(cluster_order, row_names_for_enrichment)

      enrich_mat <- base::matrix(
        NA_real_,
        nrow = base::length(row_names_for_enrichment),
        ncol = base::length(term_levels),
        dimnames = list(row_names_for_enrichment, term_levels)
      )
      enrich_q_mat <- base::matrix(
        NA_real_,
        nrow = base::length(row_names_for_enrichment),
        ncol = base::length(term_levels),
        dimnames = list(row_names_for_enrichment, term_levels)
      )
      q_lookup <- list()
      for (cl in plot_cluster_order) {
        cl_res <- res[[cl]]
        if (base::is.null(cl_res) || base::nrow(cl_res) == 0) {
          next
        }
        cl_res <- dplyr::filter(cl_res, qvalue <= qval)
        if (base::nrow(cl_res) == 0) {
          next
        }
        cl_res <- cl_res[base::order(cl_res$qvalue, decreasing = FALSE), , drop = FALSE]
        keep_n <- base::min(top, base::nrow(cl_res))
        cl_res <- cl_res[base::seq_len(keep_n), c("Description", "qvalue"), drop = FALSE]
        q_lookup[[cl]] <- stats::setNames(cl_res$qvalue, cl_res$Description)
      }
      for (r in base::seq_len(base::nrow(ggplot_df))) {
        cl <- base::as.character(ggplot_df$cluster[r])
        tr <- base::as.character(ggplot_df$terms[r])
        rn <- row_name_by_cluster[[cl]]
        if (!base::is.null(rn) && !base::is.na(rn) && rn %in% base::rownames(enrich_mat) && tr %in% base::colnames(enrich_mat)) {
          enrich_mat[rn, tr] <- 1
          qv <- NA_real_
          if (!base::is.null(q_lookup[[cl]]) && tr %in% base::names(q_lookup[[cl]])) {
            qv <- q_lookup[[cl]][[tr]]
          }
          enrich_q_mat[rn, tr] <- qv
        }
      }
      row_has_hits <- base::rowSums(!base::is.na(enrich_mat)) > 0
      row_has_hits_map <- stats::setNames(row_has_hits, base::rownames(enrich_mat))
      col_has_hits <- base::colSums(!base::is.na(enrich_mat)) > 0
      term_line_color <- base::rep("grey70", base::ncol(enrich_mat))
      base::names(term_line_color) <- base::colnames(enrich_mat)
      if (base::ncol(enrich_mat) > 0) {
        for (jj in base::seq_len(base::ncol(enrich_mat))) {
          hit_rows <- which(!base::is.na(enrich_mat[, jj]))
          if (base::length(hit_rows) == 0) {
            next
          }
          first_hit_row <- base::rownames(enrich_mat)[hit_rows[1]]
          row_col <- row_color_by_name[[first_hit_row]]
          if (!base::is.null(row_col) && !base::is.na(row_col)) {
            term_line_color[jj] <- row_col
          }
        }
      }
      q_legend_breaks <- base::c(qval, qval / 5, qval / 25)
      q_legend_breaks[q_legend_breaks <= 0] <- qval
      # Keep legend in low -> high significance direction.
      q_legend_breaks <- base::sort(base::unique(q_legend_breaks), decreasing = TRUE)
      q_legend_breaks <- q_legend_breaks[q_legend_breaks <= qval]
      if (base::length(q_legend_breaks) == 0) {
        q_legend_breaks <- qval
      }
      sig_scores <- base::vapply(q_legend_breaks, q_to_sig_score, FUN.VALUE = base::numeric(1))
      sig_sizes_mm <- (2.8 + (2.6 * sig_scores)) * overall_plot_scale
      sig_level_names <- if (base::length(q_legend_breaks) == 1) {
        "Low"
      } else if (base::length(q_legend_breaks) == 2) {
        c("Low", "High")
      } else if (base::length(q_legend_breaks) == 3) {
        c("Low", "Medium", "High")
      } else {
        base::c("Low", "Medium", "High", base::rep("Very high", base::length(q_legend_breaks) - 3))
      }
      sig_labels <- base::paste0(
        sig_level_names,
        "\n(q <= ",
        format(q_legend_breaks, digits = 2, scientific = TRUE),
        ")"
      )
      sig_legend_graphics <- base::lapply(sig_sizes_mm, function(sz) {
        force(sz)
        function(x, y, w, h) {
          grid::grid.points(
            x = x,
            y = y,
            pch = 16,
            size = grid::unit(sz * 0.66, "mm"),
            gp = grid::gpar(col = "grey30")
          )
        }
      })
      # Match enrichment-significance legend text size to module labels (M1, M2, ...).
      sig_legend_label_fontsize <- base::max(5.2, label_fontsize * 0.88)
      sig_legend_title_fontsize <- base::max(sig_legend_label_fontsize, label_fontsize * 0.94 + 0.6)
      if (!base::is.null(legend_fontsize)) {
        sig_legend_label_fontsize <- base::as.numeric(legend_fontsize)
        sig_legend_title_fontsize <- base::as.numeric(legend_fontsize) + 0.8
      }
      sig_legend <- ComplexHeatmap::Legend(
        title = "Enrichment\nsignificance",
        labels = sig_labels,
        labels_gp = grid::gpar(fontsize = sig_legend_label_fontsize, col = "grey20"),
        title_gp = grid::gpar(fontsize = sig_legend_title_fontsize, fontface = "bold", col = "grey20"),
        graphics = sig_legend_graphics,
        gap = grid::unit(0.8, "mm")
      )
      n_enrich_cols <- base::ncol(enrich_mat)
      enrich_cell_w_mm_base <- if (n_enrich_cols <= 8) {
        8.2
      } else if (n_enrich_cols <= 14) {
        7.0
      } else if (n_enrich_cols <= 22) {
        6.0
      } else {
        5.2
      }
      enrich_cell_w_mm <- enrich_cell_w_mm_base * enrichment_column_spacing_scale
      enrich_cell_w_mm <- base::max(2.8, base::min(9.0, enrich_cell_w_mm))
      enrich_cell_w_mm <- enrich_cell_w_mm * overall_plot_scale
      enrich_cell_w_mm <- enrich_cell_w_mm * enrichment_panel_compact_scale
      enrichment_body_w_mm <- base::max(24, n_enrich_cols * enrich_cell_w_mm)
      line_width_scale <- base::max(0.2, base::min(4, enrichment_line_width_scale))
      lwd_row <- 1.1 * line_width_scale
      lwd_vertical <- 0.9 * line_width_scale

      enrichment_ht <- ComplexHeatmap::Heatmap(
        enrich_mat,
        name = "enrichment",
        col = c("1" = "white"),
        na_col = "white",
        show_heatmap_legend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        border = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_names_rot = 90,
        column_labels = term_levels_display,
        column_names_gp = grid::gpar(fontsize = term_label_fontsize),
        column_names_max_height = grid::unit(enrichment_colname_max_cm, "cm"),
        width = grid::unit(enrichment_body_w_mm, "mm"),
        height = grid::unit(hc_body_h_mm, "mm"),
        rect_gp = grid::gpar(col = NA),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (enrichment_vertical_line_mode == "full" && i == 1 && isTRUE(col_has_hits[j])) {
            grid::grid.lines(
              x = base::c(x, x),
              y = base::c(grid::unit(0, "npc"), grid::unit(1, "npc")),
              gp = grid::gpar(col = term_line_color[j], lwd = lwd_vertical, alpha = 0.4)
            )
          }
          row_name <- base::rownames(enrich_mat)[i]
          row_col <- row_color_by_name[[row_name]]
          row_has_hit <- row_has_hits_map[[row_name]]
          if (base::is.null(row_has_hit) || isFALSE(row_has_hit)) {
            row_col <- "grey75"
            line_alpha <- 0.35
          } else {
            if (base::is.null(row_col) || base::is.na(row_col)) {
              row_col <- "grey70"
            }
            line_alpha <- 0.5
          }
          grid::grid.lines(
            x = base::c(x - w * 0.5, x + w * 0.5),
            y = base::c(y, y),
            gp = grid::gpar(col = row_col, lwd = lwd_row, alpha = line_alpha)
          )
          if (!base::is.na(enrich_mat[i, j])) {
            qv <- enrich_q_mat[i, j]
            sig_score <- q_to_sig_score(qv)
            point_size_mm <- (2.8 + (2.6 * sig_score)) * overall_plot_scale
            if (enrichment_vertical_line_mode == "to_term") {
              grid::grid.lines(
                x = base::c(x, x),
                y = base::c(y, grid::unit(0, "npc")),
                gp = grid::gpar(col = row_col, lwd = lwd_vertical, alpha = 0.42)
              )
            }
            grid::grid.points(
              x = x,
              y = y,
              pch = 16,
              size = grid::unit(point_size_mm, "mm"),
              gp = grid::gpar(col = row_col, alpha = 0.92)
            )
          }
        }
      )

      cp <- if (heatmap_side == "left") {
        hc_ht + enrichment_ht
      } else {
        enrichment_ht + hc_ht
      }
      panel_border_gp <- grid::gpar(col = "black", fill = NA, lwd = 1)
      sig_legend_width_mm <- grid::convertWidth(
        grid::grobWidth(sig_legend@grob),
        "mm",
        valueOnly = TRUE
      )
      draw_padding_mm <- c(
        8,
        base::max(22, sig_legend_width_mm + 20),
        6,
        6
      ) * overall_plot_scale
      draw_padding <- grid::unit(draw_padding_mm, "mm")
      draw_with_body_borders <- function(ht_obj, ...) {
        drawn_ht <- ComplexHeatmap::draw(
          ht_obj,
          newpage = TRUE,
          merge_legends = TRUE,
          annotation_legend_list = list(sig_legend),
          annotation_legend_side = "right",
          heatmap_legend_side = "right",
          show_annotation_legend = TRUE,
          show_heatmap_legend = TRUE,
          ...
        )
        for (nm in c("GFC", "enrichment")) {
          try(
            ComplexHeatmap::decorate_heatmap_body(nm, {
              # Draw border slightly inset so PDF viewport clipping does not trim the top edge.
              grid::grid.rect(
                x = grid::unit(0.5, "npc"),
                y = grid::unit(0.5, "npc"),
                width = grid::unit(0.996, "npc"),
                height = grid::unit(0.996, "npc"),
                gp = panel_border_gp
              )
            }),
            silent = TRUE
          )
        }
        invisible(drawn_ht)
      }
      pdf_height <- base::max(
        8.0,
        base::min(
          24.0,
          (
            hc_body_h_mm +
              base::max(20, enrichment_colname_max_cm * 10) +
              34
          ) / 25.4
        )
      )
      hc_annotation_extra_mm <- (module_box_width_cm * 10) +
        (if (heatmap_show_gene_counts) {
          14 * overall_plot_scale
        } else {
          2
        }) +
        (2 * overall_plot_scale)
      pdf_width <- base::max(
        10.0,
        base::min(
          38.0,
          (
            hc_body_w_mm +
              hc_annotation_extra_mm +
              enrichment_body_w_mm +
              sig_legend_width_mm +
              42
          ) / 25.4
        )
      )

      Cairo::CairoPDF(
        file = base::paste0(
          hcobject[["working_directory"]][["dir_output"]],
          hcobject[["global_settings"]][["save_folder"]],
          "/Enrichment_",
          i_title,
          "_top_",
          top,
          ".pdf"
        ),
        width = if (base::is.null(pdf_width_input)) pdf_width else as.numeric(pdf_width_input),
        height = if (base::is.null(pdf_height_input)) pdf_height else as.numeric(pdf_height_input),
        pointsize = pdf_pointsize
      )
      draw_with_body_borders(
        cp,
        column_title = base::paste0(i_title, " enrichment"),
        column_title_gp = grid::gpar(fontsize = font_title, fontface = "bold"),
        padding = draw_padding
      )
      grDevices::dev.off()
      cp_w_lgd <- draw_with_body_borders(
        cp,
        column_title = base::paste0(i_title, " enrichment"),
        column_title_gp = grid::gpar(fontsize = font_title, fontface = "bold"),
        padding = draw_padding
      )

      openxlsx::write.xlsx(
        x = build_enrichment_export(res, selected_summary, significant_summary),
        file = base::paste0(
          hcobject[["working_directory"]][["dir_output"]],
          hcobject[["global_settings"]][["save_folder"]],
          "/Enrichment_",
          i_title,
          ".xlsx"
        ),
        overwrite = TRUE
      )

      output <- list()
      output[["p"]] <- cp_w_lgd
      output[["enrichment_plot"]] <- enrichment_ht
      output[["hc_heatmap"]] <- hc_ht
      output[["module_label_map"]] <- module_label_map_current
      output[["result"]] <- top_enr_out
      output[["enrichment"]] <- res
      output[["selected_enrichments"]] <- selected_summary
      output[["significant_enrichments"]] <- significant_summary
      hcobject[["integrated_output"]][["enrichments"]][[base::paste0("top_", i)]] <<- output
    } else {
      print(base::paste0("invalid database: ", i))
    }
  }

  combined_selected_summary <- combine_summary_tables(selected_summary_all_dbs)
  combined_significant_summary <- combine_summary_tables(significant_summary_all_dbs)
  hcobject[["integrated_output"]][["enrichments"]][["top_all_dbs"]] <<- NULL
  hcobject[["integrated_output"]][["enrichments"]][["top_all_dbs_mixed"]] <<- NULL

  if (base::nrow(combined_selected_summary) > 0) {
    combined_plot_df <- combined_selected_summary
    combined_plot_df <- combined_plot_df[
      !base::is.na(combined_plot_df$cluster) &
      !base::is.na(combined_plot_df$term) &
      combined_plot_df$term != "" &
      combined_plot_df$cluster %in% cluster_order,
      ,
      drop = FALSE
    ]
    if (base::nrow(combined_plot_df) > 0) {
      combined_plot_df$database <- base::as.character(combined_plot_df$database)
      combined_plot_df$cluster <- base::as.character(combined_plot_df$cluster)
      combined_plot_df$term <- base::as.character(combined_plot_df$term)
      combined_plot_df$rank <- suppressWarnings(base::as.numeric(combined_plot_df$rank))
      combined_plot_df$qvalue <- suppressWarnings(base::as.numeric(combined_plot_df$qvalue))

      combined_plot_df <- combined_plot_df[
        base::order(
          combined_plot_df$database,
          combined_plot_df$rank,
          combined_plot_df$cluster,
          combined_plot_df$term
        ),
        ,
        drop = FALSE
      ]
      combined_plot_df$term_with_db <- base::paste0(
        "[",
        combined_plot_df$database,
        "] ",
        combined_plot_df$term
      )

      term_levels_all <- interleaved_mixed_term_order(
        plot_df = combined_plot_df,
        cluster_levels = cluster_order
      )
      raw_term_levels_all <- base::unique(combined_plot_df$term_with_db)
      term_levels_all <- term_levels_all[term_levels_all %in% raw_term_levels_all]
      term_levels_all <- base::c(
        term_levels_all,
        base::setdiff(raw_term_levels_all, term_levels_all)
      )
      row_names_for_enrichment <- base::rownames(heatmap_mat)
      row_name_by_cluster <- stats::setNames(row_names_for_enrichment, cluster_order)
      row_color_by_name <- stats::setNames(cluster_order, row_names_for_enrichment)
      enrich_mat_all <- base::matrix(
        NA_real_,
        nrow = base::length(row_names_for_enrichment),
        ncol = base::length(term_levels_all),
        dimnames = list(row_names_for_enrichment, term_levels_all)
      )
      enrich_q_mat_all <- base::matrix(
        NA_real_,
        nrow = base::length(row_names_for_enrichment),
        ncol = base::length(term_levels_all),
        dimnames = list(row_names_for_enrichment, term_levels_all)
      )

      for (r in base::seq_len(base::nrow(combined_plot_df))) {
        cl <- combined_plot_df$cluster[r]
        tr <- combined_plot_df$term_with_db[r]
        rn <- row_name_by_cluster[[cl]]
        if (!base::is.null(rn) &&
            !base::is.na(rn) &&
            rn %in% base::rownames(enrich_mat_all) &&
            tr %in% base::colnames(enrich_mat_all)) {
          enrich_mat_all[rn, tr] <- 1
          qv <- combined_plot_df$qvalue[r]
          cur_qv <- enrich_q_mat_all[rn, tr]
          if (base::is.na(cur_qv) || (!base::is.na(qv) && qv < cur_qv)) {
            enrich_q_mat_all[rn, tr] <- qv
          }
        }
      }

      row_has_hits_all <- base::rowSums(!base::is.na(enrich_mat_all)) > 0
      row_has_hits_map_all <- stats::setNames(row_has_hits_all, base::rownames(enrich_mat_all))
      col_has_hits_all <- base::colSums(!base::is.na(enrich_mat_all)) > 0
      term_line_color_all <- base::rep("grey70", base::ncol(enrich_mat_all))
      base::names(term_line_color_all) <- base::colnames(enrich_mat_all)
      if (base::ncol(enrich_mat_all) > 0) {
        for (jj in base::seq_len(base::ncol(enrich_mat_all))) {
          hit_rows <- which(!base::is.na(enrich_mat_all[, jj]))
          if (base::length(hit_rows) == 0) {
            next
          }
          first_hit_row <- base::rownames(enrich_mat_all)[hit_rows[1]]
          row_col <- row_color_by_name[[first_hit_row]]
          if (!base::is.null(row_col) && !base::is.na(row_col)) {
            term_line_color_all[jj] <- row_col
          }
        }
      }

      term_levels_all_display <- maybe_wrap_labels(term_levels_all)
      max_term_chars_all <- visual_max_label_chars(term_levels_all_display)
      term_label_fontsize_all <- compute_term_label_fontsize(max_term_chars_all)
      if (!base::is.null(enrichment_label_fontsize)) {
        term_label_fontsize_all <- base::as.numeric(enrichment_label_fontsize)
      }
      enrichment_colname_max_cm_all <- estimate_rotated_label_height_cm(
        labels = term_levels_all_display,
        font_size = term_label_fontsize_all,
        min_cm = 4.8,
        max_cm = 11.2,
        extra_cm = 0.10
      )

      q_legend_breaks_all <- base::c(qval, qval / 5, qval / 25)
      q_legend_breaks_all[q_legend_breaks_all <= 0] <- qval
      q_legend_breaks_all <- base::sort(base::unique(q_legend_breaks_all), decreasing = TRUE)
      q_legend_breaks_all <- q_legend_breaks_all[q_legend_breaks_all <= qval]
      if (base::length(q_legend_breaks_all) == 0) {
        q_legend_breaks_all <- qval
      }
      sig_scores_all <- base::vapply(q_legend_breaks_all, q_to_sig_score, FUN.VALUE = base::numeric(1))
      sig_sizes_mm_all <- (2.8 + (2.6 * sig_scores_all)) * overall_plot_scale
      sig_level_names_all <- if (base::length(q_legend_breaks_all) == 1) {
        "Low"
      } else if (base::length(q_legend_breaks_all) == 2) {
        c("Low", "High")
      } else if (base::length(q_legend_breaks_all) == 3) {
        c("Low", "Medium", "High")
      } else {
        base::c("Low", "Medium", "High", base::rep("Very high", base::length(q_legend_breaks_all) - 3))
      }
      sig_labels_all <- base::paste0(
        sig_level_names_all,
        "\n(q <= ",
        format(q_legend_breaks_all, digits = 2, scientific = TRUE),
        ")"
      )
      sig_legend_graphics_all <- base::lapply(sig_sizes_mm_all, function(sz) {
        force(sz)
        function(x, y, w, h) {
          grid::grid.points(
            x = x,
            y = y,
            pch = 16,
            size = grid::unit(sz * 0.66, "mm"),
            gp = grid::gpar(col = "grey30")
          )
        }
      })
      sig_legend_label_fontsize_all <- base::max(5.2, label_fontsize * 0.88)
      sig_legend_title_fontsize_all <- base::max(sig_legend_label_fontsize_all, label_fontsize * 0.94 + 0.6)
      if (!base::is.null(legend_fontsize)) {
        sig_legend_label_fontsize_all <- base::as.numeric(legend_fontsize)
        sig_legend_title_fontsize_all <- base::as.numeric(legend_fontsize) + 0.8
      }
      sig_legend_all <- ComplexHeatmap::Legend(
        title = "Enrichment\nsignificance",
        labels = sig_labels_all,
        labels_gp = grid::gpar(fontsize = sig_legend_label_fontsize_all, col = "grey20"),
        title_gp = grid::gpar(fontsize = sig_legend_title_fontsize_all, fontface = "bold", col = "grey20"),
        graphics = sig_legend_graphics_all,
        gap = grid::unit(0.8, "mm")
      )

      n_enrich_cols_all <- base::ncol(enrich_mat_all)
      enrich_cell_w_mm_base_all <- if (n_enrich_cols_all <= 8) {
        8.2
      } else if (n_enrich_cols_all <= 14) {
        7.0
      } else if (n_enrich_cols_all <= 22) {
        6.0
      } else {
        5.2
      }
      enrich_cell_w_mm_all <- enrich_cell_w_mm_base_all * enrichment_column_spacing_scale
      enrich_cell_w_mm_all <- base::max(2.8, base::min(9.0, enrich_cell_w_mm_all))
      enrich_cell_w_mm_all <- enrich_cell_w_mm_all * overall_plot_scale
      enrich_cell_w_mm_all <- enrich_cell_w_mm_all * enrichment_panel_compact_scale
      enrichment_body_w_mm_all <- base::max(24, n_enrich_cols_all * enrich_cell_w_mm_all)
      line_width_scale_all <- base::max(0.2, base::min(4, enrichment_line_width_scale))
      lwd_row_all <- 1.1 * line_width_scale_all
      lwd_vertical_all <- 0.9 * line_width_scale_all

      term_db_map <- stats::setNames(combined_plot_df$database, combined_plot_df$term_with_db)
      term_db_levels <- base::as.character(term_db_map[base::colnames(enrich_mat_all)])
      db_order <- base::unique(base::as.character(combined_plot_df$database))
      term_db_levels <- base::factor(term_db_levels, levels = db_order)

      enrichment_ht_all <- ComplexHeatmap::Heatmap(
        enrich_mat_all,
        name = "enrichment_all",
        col = c("1" = "white"),
        na_col = "white",
        show_heatmap_legend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        border = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_names_rot = 90,
        column_labels = term_levels_all_display,
        column_names_gp = grid::gpar(fontsize = term_label_fontsize_all),
        column_names_max_height = grid::unit(enrichment_colname_max_cm_all, "cm"),
        width = grid::unit(enrichment_body_w_mm_all, "mm"),
        height = grid::unit(hc_body_h_mm, "mm"),
        column_split = term_db_levels,
        column_title_gp = grid::gpar(fontsize = font_axis, fontface = "bold"),
        rect_gp = grid::gpar(col = NA),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (enrichment_vertical_line_mode == "full" && i == 1 && isTRUE(col_has_hits_all[j])) {
            grid::grid.lines(
              x = base::c(x, x),
              y = base::c(grid::unit(0, "npc"), grid::unit(1, "npc")),
              gp = grid::gpar(col = term_line_color_all[j], lwd = lwd_vertical_all, alpha = 0.4)
            )
          }
          row_name <- base::rownames(enrich_mat_all)[i]
          row_col <- row_color_by_name[[row_name]]
          row_has_hit <- row_has_hits_map_all[[row_name]]
          if (base::is.null(row_has_hit) || isFALSE(row_has_hit)) {
            row_col <- "grey75"
            line_alpha <- 0.35
          } else {
            if (base::is.null(row_col) || base::is.na(row_col)) {
              row_col <- "grey70"
            }
            line_alpha <- 0.5
          }
          grid::grid.lines(
            x = base::c(x - w * 0.5, x + w * 0.5),
            y = base::c(y, y),
            gp = grid::gpar(col = row_col, lwd = lwd_row_all, alpha = line_alpha)
          )
          if (!base::is.na(enrich_mat_all[i, j])) {
            qv <- enrich_q_mat_all[i, j]
            sig_score <- q_to_sig_score(qv)
            point_size_mm <- (2.8 + (2.6 * sig_score)) * overall_plot_scale
            if (enrichment_vertical_line_mode == "to_term") {
              grid::grid.lines(
                x = base::c(x, x),
                y = base::c(y, grid::unit(0, "npc")),
                gp = grid::gpar(col = row_col, lwd = lwd_vertical_all, alpha = 0.42)
              )
            }
            grid::grid.points(
              x = x,
              y = y,
              pch = 16,
              size = grid::unit(point_size_mm, "mm"),
              gp = grid::gpar(col = row_col, alpha = 0.92)
            )
          }
        }
      )

      cp_all <- if (heatmap_side == "left") {
        hc_ht + enrichment_ht_all
      } else {
        enrichment_ht_all + hc_ht
      }
      panel_border_gp_all <- grid::gpar(col = "black", fill = NA, lwd = 1)
      sig_legend_width_mm_all <- grid::convertWidth(
        grid::grobWidth(sig_legend_all@grob),
        "mm",
        valueOnly = TRUE
      )
      draw_padding_mm_all <- c(
        8,
        base::max(22, sig_legend_width_mm_all + 20),
        6,
        6
      ) * overall_plot_scale
      draw_padding_all <- grid::unit(draw_padding_mm_all, "mm")
      draw_with_body_borders_all <- function(ht_obj, ...) {
        drawn_ht <- ComplexHeatmap::draw(
          ht_obj,
          newpage = TRUE,
          merge_legends = TRUE,
          annotation_legend_list = list(sig_legend_all),
          annotation_legend_side = "right",
          heatmap_legend_side = "right",
          show_annotation_legend = TRUE,
          show_heatmap_legend = TRUE,
          ...
        )
        for (nm in c("GFC", "enrichment_all")) {
          try(
            ComplexHeatmap::decorate_heatmap_body(nm, {
              grid::grid.rect(
                x = grid::unit(0.5, "npc"),
                y = grid::unit(0.5, "npc"),
                width = grid::unit(0.996, "npc"),
                height = grid::unit(0.996, "npc"),
                gp = panel_border_gp_all
              )
            }),
            silent = TRUE
          )
        }
        invisible(drawn_ht)
      }

      pdf_height_all <- base::max(
        8.2,
        base::min(
          24.0,
          (
            hc_body_h_mm +
              base::max(20, enrichment_colname_max_cm_all * 10) +
              36
          ) / 25.4
        )
      )
      hc_annotation_extra_mm <- (module_box_width_cm * 10) +
        (if (heatmap_show_gene_counts) {
          14 * overall_plot_scale
        } else {
          2
        }) +
        (2 * overall_plot_scale)
      pdf_width_all <- base::max(
        11.0,
        base::min(
          42.0,
          (
            hc_body_w_mm +
              hc_annotation_extra_mm +
              enrichment_body_w_mm_all +
              sig_legend_width_mm_all +
              46
          ) / 25.4
        )
      )

      Cairo::CairoPDF(
        file = base::paste0(
          hcobject[["working_directory"]][["dir_output"]],
          hcobject[["global_settings"]][["save_folder"]],
          "/Enrichment_All_DBs_top_",
          top,
          ".pdf"
        ),
        width = if (base::is.null(pdf_width_input)) pdf_width_all else as.numeric(pdf_width_input),
        height = if (base::is.null(pdf_height_input)) pdf_height_all else as.numeric(pdf_height_input),
        pointsize = pdf_pointsize
      )
      draw_with_body_borders_all(
        cp_all,
        column_title = "Combined enrichment (all DBs)",
        column_title_gp = grid::gpar(fontsize = font_title, fontface = "bold"),
        padding = draw_padding_all
      )
      grDevices::dev.off()
      cp_all_w_lgd <- draw_with_body_borders_all(
        cp_all,
        column_title = "Combined enrichment (all DBs)",
        column_title_gp = grid::gpar(fontsize = font_title, fontface = "bold"),
        padding = draw_padding_all
      )

      hcobject[["integrated_output"]][["enrichments"]][["top_all_dbs"]] <<- list(
        p = cp_all_w_lgd,
        enrichment_plot = enrichment_ht_all,
        hc_heatmap = hc_ht,
        module_label_map = module_label_map_current,
        result = combined_selected_summary
      )

      # Additional mixed all-DB plot: keep terms interleaved across DBs.
      # This ordering lets GO/KEGG/Hallmark terms sit next to each other.
      term_levels_mixed <- interleaved_mixed_term_order(
        plot_df = combined_plot_df,
        cluster_levels = cluster_order
      )
      term_levels_mixed <- term_levels_mixed[term_levels_mixed %in% base::colnames(enrich_mat_all)]
      term_levels_mixed <- base::c(
        term_levels_mixed,
        base::setdiff(base::colnames(enrich_mat_all), term_levels_mixed)
      )
      enrich_mat_all_mixed <- enrich_mat_all[, term_levels_mixed, drop = FALSE]
      enrich_q_mat_all_mixed <- enrich_q_mat_all[, term_levels_mixed, drop = FALSE]
      col_has_hits_all_mixed <- base::colSums(!base::is.na(enrich_mat_all_mixed)) > 0
      term_line_color_all_mixed <- term_line_color_all[term_levels_mixed]

      term_levels_mixed_display <- maybe_wrap_labels(term_levels_mixed)
      max_term_chars_all_mixed <- visual_max_label_chars(term_levels_mixed_display)
      term_label_fontsize_all_mixed <- compute_term_label_fontsize(max_term_chars_all_mixed)
      if (!base::is.null(enrichment_label_fontsize)) {
        term_label_fontsize_all_mixed <- base::as.numeric(enrichment_label_fontsize)
      }
      enrichment_colname_max_cm_all_mixed <- estimate_rotated_label_height_cm(
        labels = term_levels_mixed_display,
        font_size = term_label_fontsize_all_mixed,
        min_cm = 4.8,
        max_cm = 11.2,
        extra_cm = 0.10
      )
      n_enrich_cols_all_mixed <- base::ncol(enrich_mat_all_mixed)
      enrich_cell_w_mm_base_all_mixed <- if (n_enrich_cols_all_mixed <= 8) {
        8.2
      } else if (n_enrich_cols_all_mixed <= 14) {
        7.0
      } else if (n_enrich_cols_all_mixed <= 22) {
        6.0
      } else {
        5.2
      }
      enrich_cell_w_mm_all_mixed <- enrich_cell_w_mm_base_all_mixed * enrichment_column_spacing_scale
      enrich_cell_w_mm_all_mixed <- base::max(2.8, base::min(9.0, enrich_cell_w_mm_all_mixed))
      enrich_cell_w_mm_all_mixed <- enrich_cell_w_mm_all_mixed * overall_plot_scale
      enrich_cell_w_mm_all_mixed <- enrich_cell_w_mm_all_mixed * enrichment_panel_compact_scale
      enrichment_body_w_mm_all_mixed <- base::max(24, n_enrich_cols_all_mixed * enrich_cell_w_mm_all_mixed)

      enrichment_ht_all_mixed <- ComplexHeatmap::Heatmap(
        enrich_mat_all_mixed,
        name = "enrichment_all_mixed",
        col = c("1" = "white"),
        na_col = "white",
        show_heatmap_legend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        border = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_names_rot = 90,
        column_labels = term_levels_mixed_display,
        column_names_gp = grid::gpar(fontsize = term_label_fontsize_all_mixed),
        column_names_max_height = grid::unit(enrichment_colname_max_cm_all_mixed, "cm"),
        width = grid::unit(enrichment_body_w_mm_all_mixed, "mm"),
        height = grid::unit(hc_body_h_mm, "mm"),
        rect_gp = grid::gpar(col = NA),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (enrichment_vertical_line_mode == "full" && i == 1 && isTRUE(col_has_hits_all_mixed[j])) {
            grid::grid.lines(
              x = base::c(x, x),
              y = base::c(grid::unit(0, "npc"), grid::unit(1, "npc")),
              gp = grid::gpar(col = term_line_color_all_mixed[j], lwd = lwd_vertical_all, alpha = 0.4)
            )
          }
          row_name <- base::rownames(enrich_mat_all_mixed)[i]
          row_col <- row_color_by_name[[row_name]]
          row_has_hit <- row_has_hits_map_all[[row_name]]
          if (base::is.null(row_has_hit) || isFALSE(row_has_hit)) {
            row_col <- "grey75"
            line_alpha <- 0.35
          } else {
            if (base::is.null(row_col) || base::is.na(row_col)) {
              row_col <- "grey70"
            }
            line_alpha <- 0.5
          }
          grid::grid.lines(
            x = base::c(x - w * 0.5, x + w * 0.5),
            y = base::c(y, y),
            gp = grid::gpar(col = row_col, lwd = lwd_row_all, alpha = line_alpha)
          )
          if (!base::is.na(enrich_mat_all_mixed[i, j])) {
            qv <- enrich_q_mat_all_mixed[i, j]
            sig_score <- q_to_sig_score(qv)
            point_size_mm <- (2.8 + (2.6 * sig_score)) * overall_plot_scale
            if (enrichment_vertical_line_mode == "to_term") {
              grid::grid.lines(
                x = base::c(x, x),
                y = base::c(y, grid::unit(0, "npc")),
                gp = grid::gpar(col = row_col, lwd = lwd_vertical_all, alpha = 0.42)
              )
            }
            grid::grid.points(
              x = x,
              y = y,
              pch = 16,
              size = grid::unit(point_size_mm, "mm"),
              gp = grid::gpar(col = row_col, alpha = 0.92)
            )
          }
        }
      )

      cp_all_mixed <- if (heatmap_side == "left") {
        hc_ht + enrichment_ht_all_mixed
      } else {
        enrichment_ht_all_mixed + hc_ht
      }
      pdf_height_all_mixed <- base::max(
        8.2,
        base::min(
          24.0,
          (
            hc_body_h_mm +
              base::max(20, enrichment_colname_max_cm_all_mixed * 10) +
              36
          ) / 25.4
        )
      )
      hc_annotation_extra_mm <- (module_box_width_cm * 10) +
        (if (heatmap_show_gene_counts) {
          14 * overall_plot_scale
        } else {
          2
        }) +
        (2 * overall_plot_scale)
      pdf_width_all_mixed <- base::max(
        11.5,
        base::min(
          44.0,
          (
            hc_body_w_mm +
              hc_annotation_extra_mm +
              enrichment_body_w_mm_all_mixed +
              sig_legend_width_mm_all +
              48
          ) / 25.4
        )
      )

      Cairo::CairoPDF(
        file = base::paste0(
          hcobject[["working_directory"]][["dir_output"]],
          hcobject[["global_settings"]][["save_folder"]],
          "/Enrichment_All_DBs_mixed_top_",
          top,
          ".pdf"
        ),
        width = if (base::is.null(pdf_width_input)) pdf_width_all_mixed else as.numeric(pdf_width_input),
        height = if (base::is.null(pdf_height_input)) pdf_height_all_mixed else as.numeric(pdf_height_input),
        pointsize = pdf_pointsize
      )
      draw_with_body_borders_all(
        cp_all_mixed,
        column_title = "Combined enrichment (all DBs, mixed)",
        column_title_gp = grid::gpar(fontsize = font_title, fontface = "bold"),
        padding = draw_padding_all
      )
      grDevices::dev.off()
      cp_all_mixed_w_lgd <- draw_with_body_borders_all(
        cp_all_mixed,
        column_title = "Combined enrichment (all DBs, mixed)",
        column_title_gp = grid::gpar(fontsize = font_title, fontface = "bold"),
        padding = draw_padding_all
      )

      hcobject[["integrated_output"]][["enrichments"]][["top_all_dbs_mixed"]] <<- list(
        p = cp_all_mixed_w_lgd,
        enrichment_plot = enrichment_ht_all_mixed,
        hc_heatmap = hc_ht,
        module_label_map = module_label_map_current,
        result = combined_selected_summary
      )
    }
  }

  openxlsx::write.xlsx(
    x = list(
      selected_enrichments_all_dbs = combined_selected_summary,
      significant_enrichments_all_dbs = combined_significant_summary
    ),
    file = base::paste0(
      hcobject[["working_directory"]][["dir_output"]],
      hcobject[["global_settings"]][["save_folder"]],
      "/Enrichment_Selected_All_DBs.xlsx"
    ),
    overwrite = TRUE
  )
  hcobject[["integrated_output"]][["enrichments"]][["selected_enrichments_all_dbs"]] <<- combined_selected_summary
  hcobject[["integrated_output"]][["enrichments"]][["significant_enrichments_all_dbs"]] <<- combined_significant_summary
  # Mirror enrichment outputs into satellite storage so S4 <-> legacy conversion
  # keeps them available for downstream plotting helpers.
  hcobject[["satellite_outputs"]][["enrichments"]] <<- hcobject[["integrated_output"]][["enrichments"]]
  hcobject[["satellite_outputs"]][["selected_enrichments_all_dbs"]] <<- combined_selected_summary
  hcobject[["satellite_outputs"]][["significant_enrichments_all_dbs"]] <<- combined_significant_summary
}
