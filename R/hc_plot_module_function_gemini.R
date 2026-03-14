#' Plot AI-assisted module function summaries
#'
#' Draws a compact overview of stored AI-assisted module function summaries with
#' a colored module box on the left and the inferred overarching function on the
#' right.
#'
#' @param hc An `HCoCenaExperiment` with stored results from
#'   `hc_module_function_llm()`.
#' @param slot_name Satellite slot name used for storage. Default is
#'   `"llm_module_function"`.
#' @param modules Optional character vector to subset modules.
#' @param fields Character vector selecting which LLM fields to plot. Supported
#'   values are `"general_processes"`, `"contextual_state"`, and
#'   `"key_regulators"`. Defaults to all three.
#' @param max_chars Maximum number of characters shown per term. Default is
#'   `90`.
#' @param text_size Numeric text size passed to `ggplot2::geom_text()`.
#' @param with_heatmap Logical. If `TRUE`, reuse the stored hCoCena cluster
#'   heatmap and place the LLM interpretations as a right-side annotation.
#'   Falls back to a text-only plot if no heatmap cache is available.
#' @param heatmap_col_order Optional character vector overriding the hCoCena
#'   heatmap column order for this LLM plot only. If `NULL` (default), the
#'   column order from the main module heatmap is reused when available.
#' @param heatmap_cluster_columns Logical. If `FALSE` (default), reuse the
#'   column order from the main hCoCena heatmap when available. If `TRUE`,
#'   cluster the columns for this LLM plot instead.
#' @param heatmap_rel_width Relative width of the heatmap panel when
#'   `with_heatmap = TRUE`.
#' @param text_rel_width Relative width of the text panel when
#'   `with_heatmap = TRUE`.
#' @param title Plot title. Default is `"AI-assisted module interpretation"`.
#' @param ... Used by the backward-compatible wrapper alias
#'   `hc_plot_module_function_gemini()`.
#'
#' @return A single plot object when one field is selected, otherwise a named
#'   list of plot objects.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- hc_plot_module_function_llm(hc)
#' print(p)
#' }
hc_plot_module_function_llm <- function(hc,
                                        slot_name = "llm_module_function",
                                        modules = NULL,
                                        fields = c("general_processes", "contextual_state", "key_regulators"),
                                        max_chars = 90,
                                        text_size = 4,
                                        with_heatmap = TRUE,
                                        heatmap_col_order = NULL,
                                        heatmap_cluster_columns = FALSE,
                                        heatmap_rel_width = 1.7,
                                        text_rel_width = 1.25,
                                        title = "AI-assisted module interpretation") {
  if (missing(hc) || is.null(hc)) {
    stop("`hc` must be provided.")
  }
  if (!is.numeric(max_chars) || length(max_chars) != 1 || is.na(max_chars) || max_chars < 20) {
    stop("`max_chars` must be a single number >= 20.")
  }
  if (!is.numeric(text_size) || length(text_size) != 1 || is.na(text_size) || text_size <= 0) {
    stop("`text_size` must be a single positive number.")
  }
  if (!is.logical(with_heatmap) || length(with_heatmap) != 1 || is.na(with_heatmap)) {
    stop("`with_heatmap` must be TRUE or FALSE.")
  }
  if (!is.null(heatmap_col_order)) {
    heatmap_col_order <- as.character(heatmap_col_order)
  }
  if (!is.logical(heatmap_cluster_columns) || length(heatmap_cluster_columns) != 1 || is.na(heatmap_cluster_columns)) {
    stop("`heatmap_cluster_columns` must be TRUE or FALSE.")
  }
  if (!is.numeric(heatmap_rel_width) || length(heatmap_rel_width) != 1 || is.na(heatmap_rel_width) || heatmap_rel_width <= 0) {
    stop("`heatmap_rel_width` must be a single positive number.")
  }
  if (!is.numeric(text_rel_width) || length(text_rel_width) != 1 || is.na(text_rel_width) || text_rel_width <= 0) {
    stop("`text_rel_width` must be a single positive number.")
  }
  fields <- unique(as.character(fields))
  fields <- fields[!is.na(fields) & fields != ""]
  valid_fields <- c("general_processes", "contextual_state", "key_regulators")
  if (length(fields) == 0 || !all(fields %in% valid_fields)) {
    stop("`fields` must contain one or more of: ", paste(valid_fields, collapse = ", "), ".")
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  stored_results <- sat[[slot_name]]
  if (is.null(stored_results) || !is.list(stored_results) || length(stored_results) == 0) {
    stop(
      "No stored module interpretations found in `hc@satellite$", slot_name,
      "`. Run `hc_module_function_llm()` first."
    )
  }

  summary_tbl <- sat[[paste0(slot_name, "_summary")]]
  if (is.null(summary_tbl) || !is.data.frame(summary_tbl) || nrow(summary_tbl) == 0) {
    summary_tbl <- .hc_llm_summary_from_results(results = stored_results, hc = hc)
  }
  if (nrow(summary_tbl) == 0) {
    stop("No module interpretation summaries are available for plotting.")
  }

  if (!is.null(modules)) {
    modules <- as.character(modules)
    summary_tbl <- summary_tbl[summary_tbl$module %in% modules, , drop = FALSE]
  }
  if (nrow(summary_tbl) == 0) {
    stop("No matching modules found for plotting.")
  }

  heatmap_info <- if (isTRUE(with_heatmap)) .hc_llm_heatmap_info(hc) else NULL
  summary_tbl <- .hc_llm_reorder_summary_for_display(summary_tbl = summary_tbl, heatmap_info = heatmap_info)
  summary_tbl$module <- as.character(summary_tbl$module)
  summary_tbl$module_color <- as.character(summary_tbl$module_color)
  summary_tbl$module_color[is.na(summary_tbl$module_color) | summary_tbl$module_color == ""] <- "grey70"
  for (nm in valid_fields) {
    if (!nm %in% colnames(summary_tbl)) {
      summary_tbl[[nm]] <- NA_character_
    }
    summary_tbl[[nm]] <- as.character(summary_tbl[[nm]])
  }
  summary_tbl$text_color <- vapply(summary_tbl$module_color, .hc_llm_darken_color, FUN.VALUE = character(1))
  summary_tbl$label_color <- vapply(summary_tbl$module_color, .hc_llm_contrast_text_color, FUN.VALUE = character(1))
  summary_tbl$module_factor <- factor(summary_tbl$module, levels = rev(summary_tbl$module))

  title_map <- c(
    general_processes = "General processes",
    contextual_state = "Contextual state",
    key_regulators = "Key regulators"
  )

  out <- lapply(fields, function(field_nm) {
    field_tbl <- summary_tbl
    field_tbl$term_plot <- vapply(
      ifelse(
        is.na(field_tbl[[field_nm]]) | field_tbl[[field_nm]] == "",
        "No interpretation available.",
        field_tbl[[field_nm]]
      ),
      .hc_llm_prepare_display_title,
      FUN.VALUE = character(1),
      max_chars = as.integer(max_chars[[1]])
    )

    plot_title <- paste0(title, ": ", title_map[[field_nm]])

    if (isTRUE(with_heatmap) && !is.null(heatmap_info) && isTRUE(heatmap_info$draw_supported)) {
      combined_grob <- .hc_llm_capture_combined_heatmap_grob(
        heatmap_info = heatmap_info,
        summary_tbl = field_tbl,
        max_chars = max_chars,
        text_size = text_size,
        heatmap_col_order = heatmap_col_order,
        heatmap_cluster_columns = heatmap_cluster_columns
      )
      return(
        patchwork::wrap_elements(full = combined_grob) +
          patchwork::plot_annotation(
            title = plot_title,
            theme = ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = text_size * 4)
            )
          )
      )
    }

    ggplot2::ggplot(field_tbl, ggplot2::aes(y = module_factor)) +
      ggplot2::geom_tile(
        ggplot2::aes(x = 1, fill = module_color),
        width = 0.34,
        height = 0.78,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = 1, label = module, color = label_color),
        fontface = "bold",
        size = text_size * 0.85,
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 1.22, xend = 1.34, yend = module_factor, color = text_color),
        linewidth = 0.6,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = 1.38, label = term_plot, color = text_color),
        hjust = 0,
        size = text_size,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_color_identity() +
      ggplot2::coord_cartesian(xlim = c(0.75, 3.6), clip = "off") +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
      ggplot2::labs(title = plot_title) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = text_size * 4),
        plot.margin = ggplot2::margin(t = 10, r = 180, b = 10, l = 4)
      )
  })
  names(out) <- fields

  if (length(out) == 1) {
    return(out[[1]])
  }
  out
}

#' @rdname hc_plot_module_function_llm
#' @export
hc_plot_module_function_gemini <- function(...) {
  hc_plot_module_function_llm(...)
}

.hc_llm_darken_color <- function(col, factor = 0.35) {
  rgb_mat <- tryCatch(grDevices::col2rgb(col) / 255, error = function(e) NULL)
  if (is.null(rgb_mat)) {
    return("#222222")
  }
  rgb_new <- pmax(0, rgb_mat[, 1] * (1 - factor))
  grDevices::rgb(rgb_new[[1]], rgb_new[[2]], rgb_new[[3]])
}

.hc_llm_contrast_text_color <- function(col) {
  rgb_mat <- tryCatch(grDevices::col2rgb(col) / 255, error = function(e) NULL)
  if (is.null(rgb_mat)) {
    return("white")
  }
  luminance <- 0.299 * rgb_mat[1, 1] + 0.587 * rgb_mat[2, 1] + 0.114 * rgb_mat[3, 1]
  if (luminance > 0.65) "black" else "white"
}

.hc_llm_reorder_summary_for_display <- function(summary_tbl, heatmap_info = NULL) {
  if (is.null(heatmap_info) || is.null(heatmap_info$module_order) || length(heatmap_info$module_order) == 0) {
    return(summary_tbl)
  }
  keep <- heatmap_info$module_order[heatmap_info$module_order %in% summary_tbl$module]
  extras <- base::setdiff(summary_tbl$module, keep)
  ord <- base::match(base::c(keep, extras), summary_tbl$module)
  ord <- ord[!base::is.na(ord)]
  summary_tbl[ord, , drop = FALSE]
}

.hc_llm_heatmap_info <- function(hc) {
  cache_info <- .hc_heatmap_cache_info(hc@integration@cluster)
  raw_heatmap_obj <- cache_info$raw_heatmap_obj
  heatmap_obj <- cache_info$heatmap_obj
  mat <- cache_info$matrix
  if (is.null(mat) || is.null(rownames(mat))) {
    return(NULL)
  }

  row_ids <- cache_info$row_order
  if (is.null(row_ids) || length(row_ids) == 0) {
    row_ids <- rownames(mat)
  }

  module_order <- row_ids
  label_map <- tryCatch(hc@integration@cluster[["module_label_map"]], error = function(e) NULL)
  if (!is.null(label_map) && length(label_map) > 0) {
    label_map <- as.character(label_map)
    map_names <- tryCatch(names(hc@integration@cluster[["module_label_map"]]), error = function(e) NULL)
    if (!is.null(map_names) && length(map_names) == length(label_map)) {
      names(label_map) <- as.character(map_names)
    }
    hit <- label_map[row_ids]
    if (length(hit) == length(row_ids) && any(!is.na(hit))) {
      module_order <- as.character(hit)
      module_order[is.na(module_order) | module_order == ""] <- row_ids[is.na(module_order) | module_order == ""]
    }
  }

  list(
    matrix = mat,
    raw_heatmap_obj = raw_heatmap_obj,
    heatmap_obj = heatmap_obj,
    draw_supported = !is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0,
    col_order = cache_info$col_order,
    row_ids = row_ids,
    module_order = module_order,
    module_by_row = module_order,
    stored_module_label_fontsize = tryCatch(as.numeric(hc@integration@cluster[["module_label_fontsize"]]), error = function(e) NA_real_),
    stored_module_label_pt_size = tryCatch(as.numeric(hc@integration@cluster[["module_label_pt_size"]]), error = function(e) NA_real_),
    stored_module_box_width_cm = tryCatch(as.numeric(hc@integration@cluster[["module_box_width_cm"]]), error = function(e) NA_real_)
  )
}

.hc_llm_capture_heatmap_grob <- function(heatmap_obj) {
  clone <- .hc_safe_deep_clone(heatmap_obj, context = "LLM module function heatmap")
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      clone,
      newpage = FALSE,
      merge_legends = TRUE,
      show_annotation_legend = TRUE,
      show_heatmap_legend = TRUE
    )
  )
}

.hc_llm_extract_heatmap_matrix <- function(heatmap_obj) {
  if (inherits(heatmap_obj, "Heatmap")) {
    return(tryCatch(heatmap_obj@matrix, error = function(e) NULL))
  }
  if (inherits(heatmap_obj, "HeatmapList")) {
    return(tryCatch(heatmap_obj@ht_list[[1]]@matrix, error = function(e) NULL))
  }
  NULL
}

.hc_llm_capture_combined_heatmap_grob <- function(heatmap_info,
                                                  summary_tbl,
                                                  max_chars,
                                                  text_size,
                                                  heatmap_col_order = NULL,
                                                  heatmap_cluster_columns = FALSE) {
  source_obj <- if (!is.null(heatmap_info$heatmap_obj)) heatmap_info$heatmap_obj else heatmap_info$raw_heatmap_obj
  mat <- heatmap_info$matrix
  if (is.null(mat) && !is.null(source_obj)) {
    mat <- .hc_llm_extract_heatmap_matrix(source_obj)
  }
  if (is.null(mat) || is.null(rownames(mat)) || nrow(mat) == 0) {
    stop("No valid heatmap matrix available for combined LLM visualization.")
  }

  row_ids <- as.character(heatmap_info$row_ids)
  module_by_row <- as.character(heatmap_info$module_by_row)
  keep_row <- module_by_row %in% as.character(summary_tbl$module)
  row_ids <- row_ids[keep_row]
  module_by_row <- module_by_row[keep_row]
  if (length(row_ids) == 0) {
    stop("No overlapping modules between LLM summaries and stored heatmap.")
  }

  summary_tbl$module <- as.character(summary_tbl$module)
  if (!"term_plot" %in% colnames(summary_tbl)) {
    if (!"short_title" %in% colnames(summary_tbl)) {
      summary_tbl$short_title <- NA_character_
    }
    summary_tbl$term_plot <- vapply(
      ifelse(
        is.na(summary_tbl$short_title) | summary_tbl$short_title == "",
        as.character(summary_tbl$overarching_function),
        as.character(summary_tbl$short_title)
      ),
      .hc_llm_prepare_display_title,
      FUN.VALUE = character(1),
      max_chars = max_chars
    )
  } else {
    summary_tbl$term_plot <- vapply(
      as.character(summary_tbl$term_plot),
      .hc_llm_prepare_display_title,
      FUN.VALUE = character(1),
      max_chars = max_chars
    )
  }
  summary_tbl$text_color <- vapply(
    as.character(summary_tbl$module_color),
    .hc_llm_darken_color,
    FUN.VALUE = character(1)
  )
  summary_tbl$module_color <- as.character(summary_tbl$module_color)
  summary_tbl$module_color[is.na(summary_tbl$module_color) | summary_tbl$module_color == ""] <- "grey70"
  summary_tbl$label_color <- vapply(summary_tbl$module_color, .hc_llm_contrast_text_color, FUN.VALUE = character(1))

  fallback_col_ids <- if (!is.null(source_obj)) .hc_llm_extract_column_ids(source_obj, mat) else base::colnames(mat)
  prepared_cols <- .hc_prepare_plot_heatmap_columns(
    mat = mat[row_ids, , drop = FALSE],
    cluster_columns = heatmap_cluster_columns,
    plot_order = heatmap_col_order,
    main_order = heatmap_info$col_order,
    fallback_order = fallback_col_ids,
    context = "LLM module function heatmap"
  )
  mat_use <- prepared_cols$mat

  module_lookup <- stats::setNames(summary_tbl$term_plot, summary_tbl$module)
  color_lookup <- stats::setNames(summary_tbl$text_color, summary_tbl$module)
  fill_lookup <- stats::setNames(summary_tbl$module_color, summary_tbl$module)
  label_lookup <- stats::setNames(summary_tbl$label_color, summary_tbl$module)

  term_labels <- unname(module_lookup[module_by_row])
  term_cols <- unname(color_lookup[module_by_row])
  fill_cols <- unname(fill_lookup[module_by_row])
  label_cols <- unname(label_lookup[module_by_row])

  term_labels[is.na(term_labels)] <- ""
  term_cols[is.na(term_cols) | term_labels == ""] <- "#00000000"
  fill_cols[is.na(fill_cols) | fill_cols == ""] <- "grey70"
  label_cols[is.na(label_cols) | label_cols == ""] <- "black"

  max_width_chars <- if (length(term_labels) == 0) 0 else max(nchar(term_labels), na.rm = TRUE)
  text_width_cm <- max(7.5, min(18, (max_width_chars * 0.14) + 0.8))
  module_levels <- unique(module_by_row)
  module_col_map <- stats::setNames(fill_lookup[module_levels], module_levels)
  module_col_map[is.na(module_col_map) | module_col_map == ""] <- "grey70"
  stored_module_label_fontsize <- suppressWarnings(as.numeric(heatmap_info$stored_module_label_fontsize[[1]]))
  if (!is.finite(stored_module_label_fontsize) || stored_module_label_fontsize <= 0) {
    stored_module_label_fontsize <- NULL
  }
  stored_module_label_pt_size <- suppressWarnings(as.numeric(heatmap_info$stored_module_label_pt_size[[1]]))
  if (!is.finite(stored_module_label_pt_size) || stored_module_label_pt_size <= 0) {
    stored_module_label_pt_size <- NULL
  }
  stored_module_box_width_cm <- suppressWarnings(as.numeric(heatmap_info$stored_module_box_width_cm[[1]]))
  if (!is.finite(stored_module_box_width_cm) || stored_module_box_width_cm <= 0) {
    stored_module_box_width_cm <- NULL
  }

  label_fontsize <- max(8, min(10.5, 17 / max(1.5, max(nchar(module_by_row), na.rm = TRUE))))
  if (!is.null(stored_module_label_fontsize)) {
    label_fontsize <- stored_module_label_fontsize
  }
  module_box_width_cm <- if (!is.null(stored_module_box_width_cm)) {
    stored_module_box_width_cm
  } else {
    max(0.65, min(3.2, 0.34 + (0.14 * max(nchar(module_by_row), na.rm = TRUE)) + (0.035 * label_fontsize)))
  }
  module_pt_size <- if (!is.null(stored_module_label_pt_size)) {
    stored_module_label_pt_size
  } else {
    max(0.16, min(0.72, 0.06 * label_fontsize))
  }
  cell_mm <- 6.5
  heatmap_body_w_mm <- max(18, ncol(mat_use) * cell_mm)
  heatmap_body_h_mm <- max(20, nrow(mat_use) * cell_mm)

  module_box_anno <- ComplexHeatmap::anno_simple(
    module_by_row,
    col = module_col_map,
    pch = module_by_row,
    pt_gp = grid::gpar(col = "white", fontsize = label_fontsize, fontface = "bold"),
    pt_size = grid::unit(module_pt_size, "snpc"),
    simple_anno_size = grid::unit(module_box_width_cm, "cm"),
    gp = grid::gpar(col = "black"),
    which = "row"
  )

  right_anno <- ComplexHeatmap::HeatmapAnnotation(
    modules = module_box_anno,
    llm = ComplexHeatmap::anno_text(
      term_labels,
      which = "row",
      just = "left",
      location = 0,
      gp = grid::gpar(col = term_cols, fontsize = max(7.2, text_size * 2.25)),
      width = grid::unit(text_width_cm, "cm")
    ),
    which = "row",
    show_legend = FALSE,
    show_annotation_name = FALSE,
    gap = grid::unit(2, "mm")
  )

  row_dend <- NULL
  col_dend <- prepared_cols$col_dend
  gfc_lim <- suppressWarnings(max(abs(mat_use), na.rm = TRUE))
  if (!is.finite(gfc_lim) || gfc_lim <= 0) {
    gfc_lim <- 2
  }
  gfc_colors <- .hc_default_gfc_colors()
  gfc_breaks <- seq(-gfc_lim, gfc_lim, length.out = length(gfc_colors))
  gfc_col_fun <- circlize::colorRamp2(gfc_breaks, gfc_colors)

  combined_ht <- ComplexHeatmap::Heatmap(
    mat_use,
    name = "GFC",
    col = gfc_col_fun,
    cluster_rows = if (is.null(row_dend)) FALSE else row_dend,
    cluster_columns = if (is.null(col_dend)) FALSE else col_dend,
    show_row_names = FALSE,
    show_heatmap_legend = FALSE,
    show_row_dend = !is.null(row_dend),
    show_column_dend = !is.null(col_dend),
    right_annotation = right_anno,
    width = grid::unit(heatmap_body_w_mm, "mm"),
    height = grid::unit(heatmap_body_h_mm, "mm"),
    rect_gp = grid::gpar(col = "grey35", lwd = 0.7),
    row_names_side = "right",
    column_names_side = "bottom",
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = max(8, text_size * 2.4))
  )
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      combined_ht,
      newpage = FALSE,
      merge_legends = TRUE,
      show_annotation_legend = FALSE,
      show_heatmap_legend = FALSE,
      padding = grid::unit(c(5, 8, 8, 5), "mm")
    )
  )
}

.hc_llm_extract_column_ids <- function(heatmap_obj, mat) {
  col_ids <- colnames(mat)
  ord <- tryCatch(ComplexHeatmap::column_order(heatmap_obj), error = function(e) NULL)
  if (is.list(ord) && length(ord) > 0) {
    ord <- ord[[1]]
  }
  if (is.numeric(ord) && length(ord) == ncol(mat)) {
    col_ids <- col_ids[ord]
  } else if (is.character(ord) && length(ord) > 0) {
    col_ids <- ord[ord %in% col_ids]
  }
  col_ids
}

.hc_llm_extract_dendrogram <- function(heatmap_obj, which = c("row", "column")) {
  which <- match.arg(which)
  fn <- if (identical(which, "row")) ComplexHeatmap::row_dend else ComplexHeatmap::column_dend
  out <- tryCatch(fn(heatmap_obj), error = function(e) NULL)
  if (is.list(out) && length(out) > 0) {
    out <- out[[1]]
  }
  if (inherits(out, "dendrogram") || inherits(out, "hclust")) {
    return(out)
  }
  NULL
}

.hc_llm_prepare_display_term <- function(term, max_chars = 90) {
  term <- as.character(term[[1]])
  if (!nzchar(term) || is.na(term)) {
    return("No interpretation available.")
  }
  term <- .hc_llm_clean_text(term)
  term <- stringr::str_squish(term)
  term <- sub("^This module primarily reflects\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module reflects\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module primarily captures\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module captures\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module is characterized by\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module represents\\s+", "", term, ignore.case = TRUE)
  term <- sub("\\.$", "", term)
  if (is.infinite(max_chars) || is.na(max_chars) || max_chars <= 0) {
    return(term)
  }
  stringr::str_trunc(term, width = max_chars)
}

.hc_llm_prepare_display_title <- function(term, max_chars = 60) {
  term <- .hc_llm_short_title_fallback(term = term, max_chars = max_chars)
  term <- .hc_llm_clean_text(term)
  term <- stringr::str_squish(term)
  if (!nzchar(term) || is.na(term)) {
    return("No interpretation available")
  }
  stringr::str_trunc(term, width = max_chars)
}
