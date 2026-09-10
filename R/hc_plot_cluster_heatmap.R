#' Plot Cluster Heatmap
#'
#' Plots a heatmap with sample groups as columns and gene clusters as rows. The cells are coloured according to the mean GFC of a given cluster in the respective sample group.
#'  If categorical or numerical metadata annotations or user-defined enrichments have been created with satellite functions, they will be incorporated into the heatmap as column/row annotations.
#'  A module-gene export (columns: `genes`, `module`) is written to the save
#'  folder based on the currently displayed module labels. Unsplit analyses use
#'  `Module_Gene_List.xlsx`; analyses containing split-module labels use
#'  `Module_Gene_splitted_List.xlsx` so the original export is preserved.
#'  For the available options check out the satellite functions.
#' @param col_order Defines the order in which the sample groups (conditions) appear in the heatmap.
#'  Accepts a vector of strings giving the conditions in their desired order.
#'  If `NULL` and `cluster_columns = FALSE`, the column order from the previous
#'  main hCoCena heatmap is reused when available.
#'  If `cluster_columns = TRUE`, this order is overwritten by clustering.
#' @param row_order Like col_order but with module colors, module labels (for
#'   example `"M1"`), or numeric module indices in the current heatmap order.
#' @param cluster_columns A Boolean, whether or not to cluster the columns of
#'  the heatmap. Default is FALSE so the main hCoCena column order is preserved.
#' @param cluster_rows Like cluster_columns but for rows.
#' @param k The resulting cluster_columns tree is cut into k groups. Default is 0 (no cutting).
#' @param return_HM A Boolean whether of not to return the ComplexHeatmap object to hcobject$integrated_output$cluster_calc$heatmap_cluster in addition to plotting it. Default is FALSE.
#' @param cat_as_bp  A vector of Booleans which length is equivalent to the number of categorical meta data annotations created with the satellite functions.
#'  Each Boolean states whether or not the categorical variable should be annotated as a bar plot (TRUE) or as a line plot (FALSE).
#'  If you did not perform any meta data annotation, ignore this parameter.
#' @param file_name A string giving the name of the file (with .pdf ending) to which the heatmap should be written.
#'  Use `NULL` or `FALSE` to skip PDF export. Default is `"Heatmap_modules.pdf"` in the legacy API.
#' @param module_label_mode Controls labels in module color boxes. One of "legacy", "prefix", "color", "none".
#'  "legacy" keeps the current behavior. "prefix" writes indexed labels like "M1", "M2", ...
#'  Default is "prefix".
#' @param module_prefix Prefix used when `module_label_mode = "prefix"`. Default is "M".
#' @param module_label_numbering Controls when prefix labels are numbered.
#'  One of "before_clustering", "after_clustering", or "preserve_existing".
#'  With "after_clustering", modules are numbered along the clustered row order (M1, M2, ...).
#'  With "preserve_existing", existing labels from `module_label_map` are reused
#'  (useful after module splitting, e.g. M1.1, M1.2, ...).
#' @param show_module_color_names A Boolean. If FALSE, row names with legacy module color names are hidden.
#' @param gene_count_mode Controls how gene counts per module are shown. One of "legacy", "bar_and_text", "bar", "text", "none".
#'  Default is "text".
#' @param module_label_preset Layout preset for module labels. One of "auto", "compact", "balanced", "presentation".
#'  Presets adjust automatic sizing only; explicitly set sizing parameters still
#'  define the requested style, with a final fit guard applied at draw time.
#'  Default is "balanced".
#' @param module_label_color Text color for labels drawn in module boxes.
#' @param module_label_fontsize Optional numeric fontsize for module box labels. If NULL, a data-driven size is chosen.
#' @param module_label_pt_size Optional numeric scale for label glyph size inside module boxes.
#'  Uses `snpc` units in `anno_simple()`. If NULL, a data-driven size is chosen.
#'  The drawn size is automatically capped to one common safe size per heatmap
#'  so all module names fit inside their boxes without mixed label sizes.
#' @param module_box_width_cm Optional numeric width (cm) for module boxes. If NULL, width adapts to label length and fontsize.
#' @param gene_count_fontsize Optional numeric fontsize for textual gene counts. If NULL,
#'  it defaults to the module label fontsize to keep both visually consistent.
#' @param gene_count_fontface Font face for textual gene counts. One of "plain", "bold", "italic", "bold.italic".
#' @param gene_count_renderer Renderer for textual gene counts. One of "pch" or "text".
#'  Use "pch" to match module-label rendering as closely as possible.
#'  Default is "pch".
#' @param gene_count_pt_size Optional numeric scale for gene-count glyph size when `gene_count_renderer = "pch"`.
#'  Uses `snpc` units. If NULL, it follows the module-label glyph size.
#' @param gfc_colors Optional character vector of colors for the GFC color scale.
#'  If NULL, uses the hCoCena default GFC palette (RdBu-based with darker extremes).
#' @param gfc_scale_limits Optional numeric vector controlling the module-heatmap
#'  color scale limits. Provide either one positive number (`x` -> `c(-x, x)`) or
#'  two numbers (`c(min, max)`). If NULL, uses stored limits from the previous
#'  main heatmap run, otherwise falls back to `c(-range_GFC, range_GFC)`.
#' @param gfc_legend_side Side of the GFC heatmap legend. One of `"left"`,
#'  `"right"`, or `"bottom"`. Default is `"right"`.
#' @param pdf_width Numeric PDF width in inches for saved heatmap output.
#'  Default is 50.
#' @param pdf_height Numeric PDF height in inches for saved heatmap output.
#'  Default is 30.
#' @param pdf_pointsize Numeric base pointsize for saved heatmap PDF.
#'  Default is 11.
#' @param pdf_dpi Numeric raster DPI passed to the Cairo PDF backend.
#'  Default is 300.
#' @param overall_plot_scale Numeric scaling factor for the overall heatmap output.
#'  Values > 1 enlarge the plot, values < 1 shrink it. Default is 1.
#' @param smart_column_gaps Logical. If `TRUE`, insert subtle column gaps at
#'  automatically detected condition/layer/prefix group boundaries. Default is
#'  `FALSE`.
#' @param column_gap_by Optional metadata column name used to split heatmap
#'  columns. If supplied, gaps are enabled and this explicit metadata split
#'  takes precedence over automatic smart gap detection. The metadata value must
#'  be constant within each displayed heatmap column.
#' @param column_gap_mm Numeric gap size in millimeters used when column gaps
#'  are enabled. Default is `0.6`.
#' @param include_dynamic_enrichment_slots Logical. If `TRUE`, also include
#'  dynamically named enrichment slots (e.g. `enriched_per_cluster_<db>`).
#'  Default is `FALSE` to keep the standard heatmap behavior unchanged.
#' @param celltype_bar_top_n Integer. Keep only the globally top-N cell types
#'  (by summed annotation weights) in stacked row barplots; optionally collapse
#'  the rest into `celltype_bar_other_label`. Set `NULL` to keep all.
#'  Default is `3`.
#' @param celltype_bar_include_other Logical. If `TRUE`, collapsed categories
#'  are shown as one additional segment (`celltype_bar_other_label`).
#' @param celltype_bar_other_label Label used for the collapsed segment.
#' @param celltype_bar_show_dominant Logical. If `TRUE`, add a text column next
#'  to each stacked bar with the dominant cell type per module.
#' @param celltype_bar_dominant_width_cm Width in cm of the dominant cell-type
#'  text column.
#' @param celltype_bar_mode Controls display of dynamic cell-type annotations.
#'  One of `"bar_and_text"`, `"bar"`, or `"text"`.
#'  Default is `"bar_and_text"`.
#' @param include_module_significance Logical. If `TRUE`, add a right-side row
#'  annotation with module-significance labels from
#'  `satellite_outputs[[module_significance_slot]]`.
#' @param module_significance_slot Character scalar naming the satellite slot
#'  produced by `.hc_module_condition_significance_driver()`. Default is
#'  `"module_condition_significance"`.
#' @param module_significance_method Which method to visualize. One of
#'  `"auto"`, `"wilcox"`, `"limma"`, `"lmm"`. `"auto"` prefers `wilcox`,
#'  then `limma`, then `lmm`, then summary fallback columns.
#' @param module_significance_show_qvalue Logical. If `TRUE`, print q-values
#'  next to significance stars.
#' @param module_significance_width_cm Width (cm) of the significance
#'  annotation column.
#' @param module_significance_p_cutoffs Numeric length-3 vector used for star
#'  labels (`***`, `**`, `*`).
#' @param module_significance_annotation_name Column name displayed above the
#'  significance annotation.


.hc_plot_cluster_heatmap_driver <- function(col_order = NULL,
                                            row_order = NULL,
                                            cluster_columns = FALSE,
                                            cluster_rows = TRUE,
                                            k = 0,
                                            return_HM = FALSE,
                                            cat_as_bp = NULL,
                                            file_name = "Heatmap_modules.pdf",
                                            module_label_mode = "prefix",
                                            module_prefix = "M",
                                            module_label_numbering = "after_clustering",
                                            show_module_color_names = FALSE,
                                            gene_count_mode = "text",
                                            module_label_preset = "balanced",
                                            module_label_color = "white",
                                            module_label_fontsize = NULL,
                                            module_label_pt_size = NULL,
                                            module_box_width_cm = NULL,
                                            gene_count_fontsize = NULL,
                                            gene_count_fontface = "plain",
                                            gene_count_renderer = "pch",
                                            gene_count_pt_size = NULL,
                                            gfc_colors = NULL,
                                            gfc_scale_limits = NULL,
                                            gfc_legend_side = "right",
                                            pdf_width = 50,
                                            pdf_height = 30,
                                            pdf_pointsize = 11,
                                            pdf_dpi = 300,
                                            overall_plot_scale = 1,
                                            smart_column_gaps = FALSE,
                                            column_gap_by = NULL,
                                            column_gap_mm = 0.6,
                                            include_dynamic_enrichment_slots = FALSE,
                                            celltype_bar_top_n = 3,
                                            celltype_bar_include_other = TRUE,
                                            celltype_bar_other_label = "Other",
                                            celltype_bar_show_dominant = TRUE,
                                            celltype_bar_dominant_width_cm = 2.8,
                                            celltype_bar_mode = "bar_and_text",
                                            include_module_significance = FALSE,
                                            module_significance_slot = "module_condition_significance",
                                            module_significance_method = "auto",
                                            module_significance_show_qvalue = FALSE,
                                            module_significance_width_cm = 1.6,
                                            module_significance_p_cutoffs = c(0.001, 0.01, 0.05),
                                            module_significance_annotation_name = "sig",
                                            write_module_tables = TRUE) {
  plot_cluster_heatmap_new(
    col_order = col_order,
    row_order = row_order,
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    k = k,
    return_HM = return_HM,
    cat_as_bp = cat_as_bp,
    file_name = file_name,
    module_label_mode = module_label_mode,
    module_prefix = module_prefix,
    module_label_numbering = module_label_numbering,
    show_module_color_names = show_module_color_names,
    gene_count_mode = gene_count_mode,
    module_label_preset = module_label_preset,
    module_label_color = module_label_color,
    module_label_fontsize = module_label_fontsize,
    module_label_pt_size = module_label_pt_size,
    module_box_width_cm = module_box_width_cm,
    gene_count_fontsize = gene_count_fontsize,
    gene_count_fontface = gene_count_fontface,
    gene_count_renderer = gene_count_renderer,
    gene_count_pt_size = gene_count_pt_size,
    gfc_colors = gfc_colors,
    gfc_scale_limits = gfc_scale_limits,
    gfc_legend_side = gfc_legend_side,
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    pdf_pointsize = pdf_pointsize,
    pdf_dpi = pdf_dpi,
    overall_plot_scale = overall_plot_scale,
    smart_column_gaps = smart_column_gaps,
    column_gap_by = column_gap_by,
    column_gap_mm = column_gap_mm,
    include_dynamic_enrichment_slots = include_dynamic_enrichment_slots,
    celltype_bar_top_n = celltype_bar_top_n,
    celltype_bar_include_other = celltype_bar_include_other,
    celltype_bar_other_label = celltype_bar_other_label,
    celltype_bar_show_dominant = celltype_bar_show_dominant,
    celltype_bar_dominant_width_cm = celltype_bar_dominant_width_cm,
    celltype_bar_mode = celltype_bar_mode,
    include_module_significance = include_module_significance,
    module_significance_slot = module_significance_slot,
    module_significance_method = module_significance_method,
    module_significance_show_qvalue = module_significance_show_qvalue,
    module_significance_width_cm = module_significance_width_cm,
    module_significance_p_cutoffs = module_significance_p_cutoffs,
    module_significance_annotation_name = module_significance_annotation_name,
    write_module_tables = write_module_tables
  )
}


.hc_resolve_cluster_heatmap_row_order <- function(row_order,
                                                  cluster_calc,
                                                  available_colors,
                                                  module_prefix = "M") {
  if (base::is.null(row_order)) {
    return(NULL)
  }

  available_colors <- base::unique(base::as.character(available_colors))
  available_colors <- available_colors[!base::is.na(available_colors) & base::nzchar(available_colors)]
  if (base::length(available_colors) == 0) {
    stop("No valid module rows available before applying `row_order`.")
  }

  stored_module_prefix <- tryCatch(
    cluster_calc[["module_prefix"]],
    error = function(e) NULL
  )
  if (!base::is.null(stored_module_prefix) &&
    base::length(stored_module_prefix) == 1 &&
    !base::is.na(stored_module_prefix) &&
    base::nzchar(base::as.character(stored_module_prefix[[1]]))) {
    module_prefix <- base::as.character(stored_module_prefix[[1]])
  }

  module_label_map <- .hc_normalize_module_label_map_for_split(
    module_label_map = tryCatch(cluster_calc[["module_label_map"]], error = function(e) NULL),
    available_colors = available_colors,
    module_prefix = module_prefix
  )
  module_order <- .hc_split_module_order(
    cluster_calc = cluster_calc,
    available_colors = available_colors
  )
  resolved <- .hc_resolve_modules_for_split(
    modules = row_order,
    available_colors = available_colors,
    module_label_map = module_label_map,
    module_order = module_order
  )
  unresolved <- resolved$resolution_table[
    base::as.character(resolved$resolution_table$status) != "ok", ,
    drop = FALSE
  ]
  if (base::nrow(unresolved) > 0) {
    unresolved_inputs <- base::unique(base::as.character(unresolved$input))
    available_labels <- base::unique(base::as.character(module_label_map[module_order]))
    available_labels <- available_labels[!base::is.na(available_labels) & base::nzchar(available_labels)]
    available_preview <- if (base::length(available_labels) > 0) {
      out <- base::paste(utils::head(available_labels, 12L), collapse = ", ")
      if (base::length(available_labels) > 12L) {
        out <- base::paste0(out, ", ...")
      }
      base::paste0("\nAvailable module labels: ", out)
    } else {
      ""
    }
    stop(
      "Unknown entries in `row_order`: ",
      base::paste(unresolved_inputs, collapse = ", "),
      "\nUse module labels, module colors, or numeric indices.",
      available_preview
    )
  }

  out <- base::as.character(resolved$target_colors)
  out <- out[!base::is.na(out) & base::nzchar(out)]
  if (base::length(out) == 0) {
    stop("No valid module rows available after applying `row_order`.")
  }
  out
}

.hc_module_label_map_has_split_labels <- function(module_label_map) {
  if (base::is.null(module_label_map) || base::length(module_label_map) == 0) {
    return(FALSE)
  }
  labels <- base::as.character(module_label_map)
  labels <- labels[!base::is.na(labels) & base::nzchar(labels)]
  # A split child carries a numeric suffix (e.g. "M3.1"); repeated splits append
  # further suffixes ("M3.1.2"). Anchor on the trailing ".<number>" so nested
  # splits are still detected.
  base::any(base::grepl("\\.[0-9]+$", labels))
}

.hc_module_gene_list_filename <- function(module_label_map = NULL,
                                          split_history = NULL) {
  has_split_history <- base::is.list(split_history) &&
    base::length(split_history) > 0
  if (has_split_history ||
    .hc_module_label_map_has_split_labels(module_label_map)) {
    return("Module_Gene_splitted_List.xlsx")
  }
  "Module_Gene_List.xlsx"
}

.hc_module_label_draw_width_cm <- function(module_box_width_cm,
                                           module_labels_display = NULL,
                                           user_set_module_box_width_cm = FALSE,
                                           module_sig_integrated = FALSE,
                                           max_sig_stars = 0) {
  draw_width <- module_box_width_cm
  if (isTRUE(user_set_module_box_width_cm)) {
    return(draw_width)
  }

  labels_chr <- base::as.character(module_labels_display)
  labels_chr <- labels_chr[!base::is.na(labels_chr) & base::nzchar(labels_chr)]
  max_label_chars <- if (base::length(labels_chr) == 0) {
    0
  } else {
    base::max(base::nchar(labels_chr), na.rm = TRUE)
  }
  has_split_like_labels <- base::length(labels_chr) > 0 &&
    base::any(base::grepl("\\.[0-9]+", labels_chr))
  has_sig_suffix <- isTRUE(module_sig_integrated) && max_sig_stars > 0
  label_width_step_cm <- 0.16
  split_sig_extra_cm <- 0.03 * base::min(3, max_sig_stars)

  if ((isTRUE(has_split_like_labels) || isTRUE(has_sig_suffix)) && max_label_chars > 2) {
    required_width <- 0.62 + (label_width_step_cm * (max_label_chars - 2))
    if (isTRUE(has_split_like_labels) && isTRUE(has_sig_suffix)) {
      required_width <- required_width + split_sig_extra_cm
    }
    draw_width <- base::min(4.8, base::max(draw_width, required_width))
  } else if (isTRUE(has_sig_suffix)) {
    required_width <- 0.62 + (label_width_step_cm * max_sig_stars)
    draw_width <- base::min(4.8, base::max(draw_width, required_width))
  }

  draw_width
}

.hc_cluster_heatmap_cell_size_mm <- function(n_heat_rows,
                                             n_heat_cols,
                                             duplicate_condition_width_scale = 1,
                                             module_box_width_cm_draw = 0,
                                             overall_plot_scale = 1) {
  cell_size_mm <- 5.4
  if (n_heat_rows > 20) {
    cell_size_mm <- 4.9
  }
  if (n_heat_rows > 30) {
    cell_size_mm <- 4.3
  }
  if (n_heat_rows > 45) {
    cell_size_mm <- 3.8
  }
  if (n_heat_cols > 10) {
    cell_size_mm <- base::min(cell_size_mm, 4.6)
  }
  if (n_heat_cols <= 4 && n_heat_rows <= 24) {
    min_body_w_mm <- if (n_heat_cols <= 3) 30 else 34
    min_body_h_mm <- if (n_heat_rows <= 12) 90 else 108
    max_cell_mm <- if (n_heat_rows <= 12) 10 else 8
    boosted_cell_mm <- base::max(
      min_body_w_mm / base::max(1, n_heat_cols),
      min_body_h_mm / base::max(1, n_heat_rows)
    )
    cell_size_mm <- base::max(cell_size_mm, base::min(max_cell_mm, boosted_cell_mm))
  }
  if (duplicate_condition_width_scale > 1) {
    target_module_box_to_cell_ratio <- 0.68
    min_cell_mm_from_box_ratio <- (module_box_width_cm_draw * 10) / target_module_box_to_cell_ratio
    min_cell_mm_from_box_ratio <- base::min(10, min_cell_mm_from_box_ratio)
    cell_size_mm <- base::max(cell_size_mm, min_cell_mm_from_box_ratio)
  }
  cell_size_mm * overall_plot_scale
}

.hc_module_label_text_width_cm <- function(labels,
                                           fontsize_pt,
                                           fontface = "bold") {
  labels_chr <- base::as.character(labels)
  labels_chr[base::is.na(labels_chr)] <- ""
  fontsize_pt <- base::as.numeric(fontsize_pt)
  if (base::length(fontsize_pt) == 1) {
    fontsize_pt <- base::rep(fontsize_pt, base::length(labels_chr))
  }
  base::mapply(
    FUN = function(label, font_pt) {
      if (!base::nzchar(label) || !base::is.finite(font_pt) || font_pt <= 0) {
        return(0)
      }
      tryCatch(
        grid::convertWidth(
          grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = font_pt, fontface = fontface))),
          "cm",
          valueOnly = TRUE
        ),
        error = function(e) {
          base::nchar(label) * font_pt * 0.021
        }
      )
    },
    labels_chr,
    fontsize_pt,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
}

.hc_module_label_fit_pt <- function(module_label_pt_size,
                                    module_box_width_cm,
                                    module_labels_display = NULL,
                                    n_heat_rows = 1,
                                    cell_size_mm = 5.4,
                                    module_label_fontsize = NULL,
                                    use_fontsize_request = FALSE,
                                    fontface = "bold") {
  labels_chr <- base::as.character(module_labels_display)
  labels_chr[base::is.na(labels_chr)] <- ""
  n_labels <- base::length(labels_chr)
  if (n_labels == 0) {
    return(list(
      pt_size = base::numeric(0),
      base_pt_size = base::numeric(0),
      available_width_cm = NA_real_,
      available_height_pt = NA_real_,
      width_limited = FALSE,
      height_limited = FALSE,
      shrunk = FALSE
    ))
  }

  module_label_pt_size <- .hc_first_numeric_value(module_label_pt_size)
  module_label_fontsize <- .hc_first_numeric_value(module_label_fontsize)
  module_box_width_cm <- .hc_first_numeric_value(module_box_width_cm)
  n_heat_rows <- base::max(1L, base::as.integer(n_heat_rows))
  cell_size_mm <- .hc_first_numeric_value(cell_size_mm)
  if (!base::is.finite(module_label_pt_size) || module_label_pt_size <= 0 ||
    !base::is.finite(module_box_width_cm) || module_box_width_cm <= 0 ||
    !base::is.finite(cell_size_mm) || cell_size_mm <= 0) {
    return(list(
      pt_size = base::rep(NA_real_, n_labels),
      base_pt_size = base::rep(NA_real_, n_labels),
      available_width_cm = NA_real_,
      available_height_pt = NA_real_,
      width_limited = FALSE,
      height_limited = FALSE,
      shrunk = FALSE
    ))
  }

  annotation_height_cm <- (n_heat_rows * cell_size_mm) / 10
  snpc_cm <- base::min(module_box_width_cm, annotation_height_cm)
  base_pt <- if (isTRUE(use_fontsize_request) &&
    base::is.finite(module_label_fontsize) &&
    module_label_fontsize > 0) {
    module_label_fontsize
  } else {
    module_label_pt_size * snpc_cm * 72 / 2.54
  }
  base_pt_vec <- base::rep(base_pt, n_labels)
  width_fill <- 0.95
  height_fill <- 0.95
  available_width_cm <- base::max(0.02, module_box_width_cm * width_fill)
  available_height_pt <- (cell_size_mm * 72 / 25.4) * height_fill

  # `anno_simple(pch = "M1")` draws text-like symbols a little wider than a
  # bare textGrob on some devices, so keep a tiny render guard while targeting
  # 95% of the actual module box width.
  symbol_width_guard <- 1.32
  width_at_10pt <- .hc_module_label_text_width_cm(labels_chr, fontsize_pt = 10, fontface = fontface) *
    symbol_width_guard
  max_pt_by_width <- base::ifelse(
    width_at_10pt > 0,
    10 * available_width_cm / width_at_10pt,
    Inf
  )
  preferred_min_pt <- 2.8
  absolute_min_pt <- 0.2
  finite_width_pt <- max_pt_by_width[base::is.finite(max_pt_by_width)]
  width_limit_pt <- if (base::length(finite_width_pt) > 0) {
    base::min(finite_width_pt)
  } else {
    Inf
  }

  common_pt <- base::min(base_pt, width_limit_pt, available_height_pt, na.rm = TRUE)
  if (!base::is.finite(common_pt) || common_pt <= 0) {
    common_pt <- base_pt
  }
  common_pt <- base::max(absolute_min_pt, common_pt)
  if (common_pt >= preferred_min_pt ||
    (width_limit_pt >= preferred_min_pt && available_height_pt >= preferred_min_pt && base_pt >= preferred_min_pt)) {
    common_pt <- base::max(preferred_min_pt, common_pt)
  }

  fitted_pt <- base::rep(common_pt, n_labels)
  fitted_width_cm <- .hc_module_label_text_width_cm(labels_chr, fontsize_pt = fitted_pt, fontface = fontface) *
    symbol_width_guard
  for (i in base::seq_len(5L)) {
    max_width_cm <- base::max(fitted_width_cm, na.rm = TRUE)
    if (!base::is.finite(max_width_cm) ||
      max_width_cm <= available_width_cm ||
      common_pt <= absolute_min_pt) {
      break
    }
    common_pt <- base::pmax(
      absolute_min_pt,
      common_pt * (available_width_cm / max_width_cm) * 0.995
    )
    fitted_pt <- base::rep(common_pt, n_labels)
    fitted_width_cm <- .hc_module_label_text_width_cm(labels_chr, fontsize_pt = fitted_pt, fontface = fontface) *
      symbol_width_guard
  }

  list(
    pt_size = fitted_pt,
    base_pt_size = base_pt_vec,
    fitted_width_cm = fitted_width_cm,
    available_width_cm = available_width_cm,
    available_height_pt = available_height_pt,
    width_limited = base::any(max_pt_by_width < base_pt_vec - 1e-8, na.rm = TRUE),
    height_limited = base::any(available_height_pt < base_pt_vec - 1e-8, na.rm = TRUE),
    below_preferred_min = base::is.finite(common_pt) && common_pt < preferred_min_pt,
    shrunk = base::any(fitted_pt < base_pt_vec - 1e-8, na.rm = TRUE)
  )
}

.hc_module_label_effective_fontsize <- function(module_label_fit,
                                                fallback_fontsize = NULL,
                                                fallback_pt_size = NULL,
                                                min_pt = 0.2) {
  fit_pt <- tryCatch(
    .hc_as_numeric_safely(module_label_fit$pt_size),
    error = function(e) numeric(0)
  )
  fit_pt <- fit_pt[base::is.finite(fit_pt) & fit_pt > min_pt]
  if (base::length(fit_pt) > 0) {
    return(base::min(fit_pt))
  }
  fallback <- .hc_first_numeric_value(fallback_fontsize)
  if (base::is.finite(fallback) && fallback > min_pt) {
    return(fallback)
  }
  fallback <- .hc_first_numeric_value(fallback_pt_size)
  if (base::is.finite(fallback) && fallback > min_pt) {
    return(fallback)
  }
  5
}

.hc_module_label_box_annotation <- function(values,
                                            colors,
                                            labels = NULL,
                                            label_color = "white",
                                            label_fontsize_pt = NULL,
                                            fontface = "bold",
                                            width_cm,
                                            border_gp = grid::gpar(col = "black", lwd = 0.5),
                                            which = "row",
                                            width_fill = 0.95,
                                            height_fill = 0.90,
                                            font_size_width_fill = 0.78,
                                            min_visible_font_pt = 5.5,
                                            minimum_sizing_label = "M4.2") {
  values_chr <- base::as.character(values)
  n_values <- base::length(values_chr)
  labels_chr <- if (base::is.null(labels)) {
    base::rep("", n_values)
  } else {
    label_names <- base::names(labels)
    labels_raw <- base::as.character(labels)
    if (base::length(labels_raw) == n_values) {
      labels_raw
    } else if (!base::is.null(label_names) &&
      base::length(label_names) == base::length(labels_raw) &&
      base::any(values_chr %in% base::as.character(label_names))) {
      mapped <- base::as.character(labels_raw[base::match(values_chr, base::as.character(label_names))])
      missing_mapped <- base::is.na(mapped) | !base::nzchar(mapped)
      mapped[missing_mapped] <- values_chr[missing_mapped]
      mapped
    } else if (base::length(labels_raw) > 0) {
      base::rep_len(labels_raw, n_values)
    } else {
      values_chr
    }
  }
  labels_chr[base::is.na(labels_chr)] <- ""
  if (base::length(labels_chr) != n_values && n_values > 0) {
    labels_chr <- base::rep_len(labels_chr, n_values)
  }

  color_names <- base::names(colors)
  colors_chr <- base::as.character(colors)
  if (!base::is.null(color_names)) {
    base::names(colors_chr) <- color_names
  }
  if (base::is.null(base::names(colors_chr))) {
    base::names(colors_chr) <- base::as.character(colors_chr)
  }
  width_cm <- .hc_first_numeric_value(width_cm)
  if (!base::is.finite(width_cm) || width_cm <= 0) {
    width_cm <- 0.8
  }
  label_fontsize_pt <- .hc_as_numeric_safely(label_fontsize_pt)
  label_fontsize_pt <- label_fontsize_pt[base::is.finite(label_fontsize_pt) & label_fontsize_pt > 0]
  max_font_pt <- if (base::length(label_fontsize_pt) > 0) {
    base::min(label_fontsize_pt)
  } else {
    Inf
  }
  border_col <- if (!base::is.null(border_gp$col)) border_gp$col else "black"
  border_lwd <- if (!base::is.null(border_gp$lwd)) border_gp$lwd else 0.5

  fills_chr <- base::as.character(base::unname(colors_chr[values_chr]))
  missing_fill <- base::is.na(fills_chr) | !base::nzchar(fills_chr)
  fills_chr[missing_fill] <- values_chr[missing_fill]
  valid_fill <- base::vapply(fills_chr, function(cl) {
    base::isTRUE(tryCatch({
      grDevices::col2rgb(cl)
      TRUE
    }, error = function(e) FALSE))
  }, FUN.VALUE = base::logical(1))
  fills_chr[!valid_fill] <- "#d9d9d9"

  draw_fun <- function(index, ...) {
    n <- base::length(index)
    if (n == 0) {
      return(invisible(NULL))
    }
    y <- (n - base::seq_len(n) + 0.5) / n
    fill <- fills_chr[index]
    for (i in base::seq_len(n)) {
      grid::grid.rect(
        x = grid::unit(0.5, "npc"),
        y = grid::unit(y[[i]], "npc"),
        width = grid::unit(1, "npc"),
        height = grid::unit(1 / n, "npc"),
        gp = grid::gpar(fill = fill[[i]], col = border_col, lwd = border_lwd)
      )
    }

    draw_labels <- labels_chr[index]
    has_label <- !base::is.na(draw_labels) & base::nzchar(draw_labels)
    if (!base::any(has_label)) {
      return(invisible(NULL))
    }
    all_labels <- labels_chr[!base::is.na(labels_chr) & base::nzchar(labels_chr)]
    minimum_sizing_label <- base::as.character(minimum_sizing_label)
    minimum_sizing_label <- minimum_sizing_label[!base::is.na(minimum_sizing_label) &
      base::nzchar(minimum_sizing_label)]
    sizing_labels <- base::unique(base::c(all_labels, minimum_sizing_label))
    target_width_cm <- grid::convertWidth(grid::unit(width_fill, "npc"), "cm", valueOnly = TRUE)
    row_height_pt <- grid::convertHeight(grid::unit(1 / n, "npc"), "pt", valueOnly = TRUE)
    target_height_pt <- base::max(row_height_pt * height_fill, min_visible_font_pt)
    max_font_from_box_width_pt <- grid::convertWidth(
      grid::unit(font_size_width_fill, "npc"),
      "pt",
      valueOnly = TRUE
    )
    width_at_10pt <- .hc_module_label_text_width_cm(sizing_labels, fontsize_pt = 10, fontface = fontface)
    max_width_10pt <- base::max(width_at_10pt, na.rm = TRUE)
    width_fit_pt <- if (base::is.finite(max_width_10pt) && max_width_10pt > 0) {
      10 * target_width_cm / max_width_10pt
    } else {
      Inf
    }
    font_pt <- base::min(
      max_font_pt,
      width_fit_pt,
      target_height_pt,
      max_font_from_box_width_pt,
      na.rm = TRUE
    )
    if (!base::is.finite(font_pt) || font_pt <= 0) {
      font_pt <- base::min(width_fit_pt, target_height_pt, max_font_from_box_width_pt, na.rm = TRUE)
    }
    if (!base::is.finite(font_pt) || font_pt <= 0) {
      font_pt <- 5
    }
    grid::grid.text(
      draw_labels[has_label],
      x = grid::unit(base::rep(0.5, base::sum(has_label)), "npc"),
      y = grid::unit(y[has_label], "npc"),
      gp = grid::gpar(col = label_color, fontsize = font_pt, fontface = fontface)
    )
    invisible(NULL)
  }

  ComplexHeatmap::AnnotationFunction(
    fun = draw_fun,
    fun_name = "module_label_boxes",
    which = match.arg(which, c("column", "row")),
    width = grid::unit(width_cm, "cm"),
    n = n_values,
    data_scale = c(0.5, 1.5),
    var_import = list(
      labels_chr = labels_chr,
      fills_chr = fills_chr,
      max_font_pt = max_font_pt,
      label_color = label_color,
      fontface = fontface,
      width_fill = width_fill,
      height_fill = height_fill,
      font_size_width_fill = font_size_width_fill,
      min_visible_font_pt = min_visible_font_pt,
      minimum_sizing_label = minimum_sizing_label,
      border_col = border_col,
      border_lwd = border_lwd,
      # `AnnotationFunction(var_import = ...)` re-homes `draw_fun` into an
      # isolated environment whose parent is not the hcocena namespace, so the
      # internal label-sizing helper must be imported explicitly or it cannot be
      # found at draw time.
      .hc_module_label_text_width_cm = .hc_module_label_text_width_cm
    )
  )
}

.hc_heatmap_screen_fit_scale <- function(total_width_mm,
                                         total_height_mm,
                                         device_size_mm = NULL,
                                         width_fill = 0.95,
                                         height_fill = 0.90) {
  total_width_mm <- .hc_first_numeric_value(total_width_mm[[1]])
  total_height_mm <- .hc_first_numeric_value(total_height_mm[[1]])
  if (!base::is.finite(total_width_mm) || total_width_mm <= 0 ||
    !base::is.finite(total_height_mm) || total_height_mm <= 0) {
    return(1)
  }

  if (base::is.null(device_size_mm)) {
    device_size_mm <- tryCatch(
      grDevices::dev.size("in"),
      error = function(e) c(NA_real_, NA_real_)
    )
    device_size_mm <- .hc_as_numeric_safely(device_size_mm) * 25.4
  } else {
    device_size_mm <- .hc_as_numeric_safely(device_size_mm)
  }

  if (base::length(device_size_mm) != 2 ||
    any(!base::is.finite(device_size_mm)) ||
    any(device_size_mm <= 0)) {
    return(1)
  }

  screen_scale <- base::min(
    1,
    (device_size_mm[[1]] * width_fill) / total_width_mm,
    (device_size_mm[[2]] * height_fill) / total_height_mm
  )
  if (!base::is.finite(screen_scale) || screen_scale <= 0) {
    return(1)
  }

  screen_scale
}


plot_cluster_heatmap_new <- function(col_order = NULL,
                                     row_order = NULL,
                                     cluster_columns = FALSE,
                                     cluster_rows = TRUE,
                                     k = 0,
                                     return_HM = FALSE,
                                     cat_as_bp = NULL,
                                     file_name = "Heatmap_modules.pdf",
                                     module_label_mode = "prefix",
                                     module_prefix = "M",
                                     module_label_numbering = "after_clustering",
                                     show_module_color_names = FALSE,
                                     gene_count_mode = "text",
                                     module_label_preset = "balanced",
                                     module_label_color = "white",
                                     module_label_fontsize = NULL,
                                     module_label_pt_size = NULL,
                                     module_box_width_cm = NULL,
                                     gene_count_fontsize = NULL,
                                     gene_count_fontface = "plain",
                                     gene_count_renderer = "pch",
                                     gene_count_pt_size = NULL,
                                     gfc_colors = NULL,
                                     gfc_scale_limits = NULL,
                                     gfc_legend_side = "right",
                                     pdf_width = 50,
                                     pdf_height = 30,
                                     pdf_pointsize = 11,
                                     pdf_dpi = 300,
                                     overall_plot_scale = 1,
                                     smart_column_gaps = FALSE,
                                     column_gap_by = NULL,
                                     column_gap_mm = 0.6,
                                     include_dynamic_enrichment_slots = FALSE,
                                     celltype_bar_top_n = 3,
                                     celltype_bar_include_other = TRUE,
                                     celltype_bar_other_label = "Other",
                                     celltype_bar_show_dominant = TRUE,
                                     celltype_bar_dominant_width_cm = 2.8,
                                     celltype_bar_mode = "bar_and_text",
                                     include_module_significance = FALSE,
                                     module_significance_slot = "module_condition_significance",
                                     module_significance_method = "auto",
                                     module_significance_show_qvalue = FALSE,
                                     module_significance_width_cm = 1.6,
                                     module_significance_p_cutoffs = c(0.001, 0.01, 0.05),
                                     module_significance_annotation_name = "sig",
                                     write_module_tables = TRUE) {
  pdf_width_is_default <- isTRUE(all.equal(as.numeric(pdf_width), 50))
  pdf_height_is_default <- isTRUE(all.equal(as.numeric(pdf_height), 30))

  module_label_mode <- base::match.arg(
    module_label_mode,
    choices = c("legacy", "prefix", "color", "none")
  )
  module_label_numbering <- base::match.arg(
    module_label_numbering,
    choices = c("before_clustering", "after_clustering", "preserve_existing")
  )
  gene_count_mode <- base::match.arg(
    gene_count_mode,
    choices = c("legacy", "bar_and_text", "bar", "text", "none")
  )
  module_label_preset <- base::match.arg(
    module_label_preset,
    choices = c("auto", "compact", "balanced", "presentation")
  )
  gene_count_renderer <- base::match.arg(
    gene_count_renderer,
    choices = c("pch", "text")
  )
  module_significance_method <- base::match.arg(
    module_significance_method,
    choices = c("auto", "wilcox", "limma", "lmm")
  )
  gfc_legend_side <- base::match.arg(
    gfc_legend_side,
    choices = c("left", "right", "bottom")
  )
  celltype_bar_mode <- base::match.arg(
    celltype_bar_mode,
    choices = c("bar_and_text", "bar", "text")
  )
  if (!base::is.logical(smart_column_gaps) || base::length(smart_column_gaps) != 1 ||
    base::is.na(smart_column_gaps)) {
    stop("`smart_column_gaps` must be TRUE or FALSE.")
  }
  if (!base::is.null(column_gap_by)) {
    if (!base::is.character(column_gap_by) || base::length(column_gap_by) != 1 ||
      base::is.na(column_gap_by) || !base::nzchar(base::trimws(column_gap_by))) {
      stop("`column_gap_by` must be NULL or a non-empty metadata column name.")
    }
    column_gap_by <- base::trimws(column_gap_by)
  }
  if (!base::is.numeric(column_gap_mm) || base::length(column_gap_mm) != 1 ||
    base::is.na(column_gap_mm) || !base::is.finite(column_gap_mm) || column_gap_mm <= 0) {
    stop("`column_gap_mm` must be a single positive number.")
  }
  show_celltype_bars <- celltype_bar_mode %in% c("bar_and_text", "bar")
  show_celltype_text <- (celltype_bar_mode %in% c("bar_and_text", "text")) && isTRUE(celltype_bar_show_dominant)
  if (!base::is.null(celltype_bar_top_n)) {
    if (!base::is.numeric(celltype_bar_top_n) || base::length(celltype_bar_top_n) != 1 ||
      base::is.na(celltype_bar_top_n) || celltype_bar_top_n < 1) {
      stop("`celltype_bar_top_n` must be NULL or a single integer >= 1.")
    }
    celltype_bar_top_n <- base::as.integer(celltype_bar_top_n)
  }
  if (!base::is.logical(celltype_bar_include_other) || base::length(celltype_bar_include_other) != 1 ||
    base::is.na(celltype_bar_include_other)) {
    stop("`celltype_bar_include_other` must be TRUE or FALSE.")
  }
  if (!base::is.character(celltype_bar_other_label) || base::length(celltype_bar_other_label) != 1 ||
    base::is.na(celltype_bar_other_label) || !base::nzchar(base::trimws(celltype_bar_other_label))) {
    stop("`celltype_bar_other_label` must be a non-empty character scalar.")
  }
  celltype_bar_other_label <- base::trimws(celltype_bar_other_label)
  if (!base::is.logical(celltype_bar_show_dominant) || base::length(celltype_bar_show_dominant) != 1 ||
    base::is.na(celltype_bar_show_dominant)) {
    stop("`celltype_bar_show_dominant` must be TRUE or FALSE.")
  }
  if (!base::is.numeric(celltype_bar_dominant_width_cm) || base::length(celltype_bar_dominant_width_cm) != 1 ||
    base::is.na(celltype_bar_dominant_width_cm) || celltype_bar_dominant_width_cm <= 0) {
    stop("`celltype_bar_dominant_width_cm` must be a single positive number.")
  }

  user_set_module_label_fontsize <- !base::is.null(module_label_fontsize)
  user_set_module_label_pt_size <- !base::is.null(module_label_pt_size)
  user_set_module_box_width_cm <- !base::is.null(module_box_width_cm)
  user_set_gene_count_fontsize <- !base::is.null(gene_count_fontsize)
  user_set_gene_count_pt_size <- !base::is.null(gene_count_pt_size)
  if (!base::is.null(module_label_fontsize) &&
    (!base::is.numeric(module_label_fontsize) || base::length(module_label_fontsize) != 1 || module_label_fontsize <= 0)) {
    stop("`module_label_fontsize` must be NULL or a single positive number.")
  }
  if (!base::is.null(module_label_pt_size) &&
    (!base::is.numeric(module_label_pt_size) || base::length(module_label_pt_size) != 1 || module_label_pt_size <= 0)) {
    stop("`module_label_pt_size` must be NULL or a single positive number.")
  }
  if (!base::is.null(module_box_width_cm) &&
    (!base::is.numeric(module_box_width_cm) || base::length(module_box_width_cm) != 1 || module_box_width_cm <= 0)) {
    stop("`module_box_width_cm` must be NULL or a single positive number.")
  }
  if (!base::is.null(gene_count_fontsize) &&
    (!base::is.numeric(gene_count_fontsize) || base::length(gene_count_fontsize) != 1 || gene_count_fontsize <= 0)) {
    stop("`gene_count_fontsize` must be NULL or a single positive number.")
  }
  if (!base::is.null(gene_count_pt_size) &&
    (!base::is.numeric(gene_count_pt_size) || base::length(gene_count_pt_size) != 1 || gene_count_pt_size <= 0)) {
    stop("`gene_count_pt_size` must be NULL or a single positive number.")
  }
  if (!is.character(gene_count_fontface) || base::length(gene_count_fontface) != 1 ||
    !(gene_count_fontface %in% c("plain", "bold", "italic", "bold.italic"))) {
    stop("`gene_count_fontface` must be one of: \"plain\", \"bold\", \"italic\", \"bold.italic\".")
  }
  if (!base::is.logical(show_module_color_names) || base::length(show_module_color_names) != 1) {
    stop("`show_module_color_names` must be TRUE or FALSE.")
  }
  if (!base::is.numeric(overall_plot_scale) ||
    base::length(overall_plot_scale) != 1 ||
    base::is.na(overall_plot_scale) ||
    overall_plot_scale <= 0) {
    stop("`overall_plot_scale` must be a positive numeric scalar.")
  }
  write_pdf <- FALSE
  if (base::is.null(file_name) || (base::is.logical(file_name) && base::identical(file_name, FALSE))) {
    write_pdf <- FALSE
  } else if (base::is.character(file_name) && base::length(file_name) == 1 && base::nzchar(file_name)) {
    write_pdf <- TRUE
  } else {
    stop("`file_name` must be a non-empty string, NULL, or FALSE.")
  }
  if (!base::is.logical(include_module_significance) ||
    base::length(include_module_significance) != 1 ||
    base::is.na(include_module_significance)) {
    stop("`include_module_significance` must be TRUE or FALSE.")
  }
  if (!base::is.character(module_significance_slot) ||
    base::length(module_significance_slot) != 1 ||
    !base::nzchar(module_significance_slot)) {
    stop("`module_significance_slot` must be a non-empty string.")
  }
  if (!base::is.logical(module_significance_show_qvalue) ||
    base::length(module_significance_show_qvalue) != 1 ||
    base::is.na(module_significance_show_qvalue)) {
    stop("`module_significance_show_qvalue` must be TRUE or FALSE.")
  }
  if (!base::is.numeric(module_significance_width_cm) ||
    base::length(module_significance_width_cm) != 1 ||
    base::is.na(module_significance_width_cm) ||
    module_significance_width_cm <= 0) {
    stop("`module_significance_width_cm` must be a positive number.")
  }
  module_significance_width_cm <- base::max(0.8, base::min(4.0, module_significance_width_cm))
  module_significance_p_cutoffs <- .hc_as_numeric_safely(module_significance_p_cutoffs)
  if (base::length(module_significance_p_cutoffs) != 3 ||
    base::any(!base::is.finite(module_significance_p_cutoffs))) {
    stop("`module_significance_p_cutoffs` must be a numeric vector of length 3.")
  }
  module_significance_p_cutoffs <- base::sort(module_significance_p_cutoffs)
  if (module_significance_p_cutoffs[[1]] <= 0 ||
    module_significance_p_cutoffs[[3]] > 1) {
    stop("`module_significance_p_cutoffs` must be within (0, 1].")
  }
  if (!base::is.character(module_significance_annotation_name) ||
    base::length(module_significance_annotation_name) != 1 ||
    !base::nzchar(module_significance_annotation_name)) {
    stop("`module_significance_annotation_name` must be a non-empty string.")
  }
  if (!base::is.logical(write_module_tables) ||
    base::length(write_module_tables) != 1 ||
    base::is.na(write_module_tables)) {
    stop("`write_module_tables` must be TRUE or FALSE.")
  }
  overall_plot_scale <- base::max(0.5, base::min(3, overall_plot_scale))
  if (base::is.null(gfc_colors)) {
    gfc_colors <- .hc_default_gfc_colors()
  }
  if (!base::is.character(gfc_colors) || base::length(gfc_colors) < 2) {
    stop("`gfc_colors` must be NULL or a character vector with at least two colors.")
  }
  if (any(base::is.na(gfc_colors)) || any(gfc_colors == "")) {
    stop("`gfc_colors` must not contain NA or empty strings.")
  }
  gfc_colors <- base::as.character(gfc_colors)
  if (!base::is.numeric(pdf_width) ||
    base::length(pdf_width) != 1 ||
    base::is.na(pdf_width) ||
    pdf_width <= 0) {
    stop("`pdf_width` must be a single positive number.")
  }
  if (!base::is.numeric(pdf_height) ||
    base::length(pdf_height) != 1 ||
    base::is.na(pdf_height) ||
    pdf_height <= 0) {
    stop("`pdf_height` must be a single positive number.")
  }
  if (!base::is.numeric(pdf_pointsize) ||
    base::length(pdf_pointsize) != 1 ||
    base::is.na(pdf_pointsize) ||
    pdf_pointsize <= 0) {
    stop("`pdf_pointsize` must be a single positive number.")
  }
  if (!base::is.numeric(pdf_dpi) ||
    base::length(pdf_dpi) != 1 ||
    base::is.na(pdf_dpi) ||
    pdf_dpi <= 0) {
    stop("`pdf_dpi` must be a single positive number.")
  }

  normalize_scale_limits <- function(x) {
    if (base::is.null(x)) {
      return(NULL)
    }
    x <- .hc_as_numeric_safely(x)
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

  gfc_scale_limits <- normalize_scale_limits(gfc_scale_limits)
  if (base::is.null(gfc_scale_limits)) {
    stored_limits <- tryCatch(
      normalize_scale_limits(hcobject[["integrated_output"]][["cluster_calc"]][["gfc_scale_limits"]]),
      error = function(e) NULL
    )
    if (!base::is.null(stored_limits)) {
      gfc_scale_limits <- stored_limits
    } else {
      fallback_lim <- .hc_first_numeric_value(hcobject[["global_settings"]][["range_GFC"]])
      if (base::length(fallback_lim) != 1 || any(!base::is.finite(fallback_lim)) || any(fallback_lim <= 0)) {
        fallback_lim <- 2
      }
      gfc_scale_limits <- c(-base::abs(fallback_lim), base::abs(fallback_lim))
    }
  }
  gfc_scale_breaks <- pretty(gfc_scale_limits, n = 5)
  gfc_scale_breaks <- gfc_scale_breaks[
    gfc_scale_breaks >= gfc_scale_limits[1] - .Machine$double.eps^0.5 &
      gfc_scale_breaks <= gfc_scale_limits[2] + .Machine$double.eps^0.5
  ]
  if (!any(base::abs(gfc_scale_breaks) < .Machine$double.eps^0.5)) {
    gfc_scale_breaks <- base::sort(base::unique(base::c(gfc_scale_breaks, 0)))
  }
  if (base::length(gfc_scale_breaks) < 3) {
    gfc_scale_breaks <- base::seq(gfc_scale_limits[1], gfc_scale_limits[2], length.out = 5)
  }
  gfc_scale_labels <- base::formatC(gfc_scale_breaks, format = "fg", digits = 3)
  gfc_scale_labels <- base::trimws(gfc_scale_labels)
  gfc_label_width <- base::max(base::nchar(gfc_scale_labels), na.rm = TRUE)
  gfc_scale_labels <- base::format(gfc_scale_labels, width = gfc_label_width, justify = "right")

  # --- 1. Load Satellite Outputs (Enrichments & Metadata) ---

  # Collect available row-enrichment slots (legacy + dynamic per-database slots).
  satellite_outputs <- hcobject[["satellite_outputs"]]
  preferred_slots <- base::character(0)
  dynamic_slots <- base::character(0)
  if (isTRUE(include_dynamic_enrichment_slots)) {
    if ("celltype_annotation" %in% base::names(satellite_outputs) &&
      !base::is.null(satellite_outputs[["celltype_annotation"]]) &&
      base::is.list(satellite_outputs[["celltype_annotation"]]) &&
      "annotation_slots" %in% base::names(satellite_outputs[["celltype_annotation"]])) {
      preferred_slots <- base::as.character(satellite_outputs[["celltype_annotation"]][["annotation_slots"]])
    }
    dynamic_slots <- base::names(satellite_outputs)[base::grepl("^enriched_per_cluster_", base::names(satellite_outputs))]
  }
  enrichment_slot_names <- base::unique(base::c(
    preferred_slots,
    "enriched_per_cluster",
    "enriched_per_cluster2",
    dynamic_slots
  ))
  enrichment_slot_names <- enrichment_slot_names[enrichment_slot_names %in% base::names(satellite_outputs)]

  user_enrichment_slots <- list()
  for (slot_nm in enrichment_slot_names) {
    slot_obj <- satellite_outputs[[slot_nm]]
    if (base::is.null(slot_obj) || !base::is.list(slot_obj) || !("categories_per_cluster" %in% base::names(slot_obj))) {
      next
    }
    slot_df <- slot_obj[["categories_per_cluster"]]
    if (!base::is.data.frame(slot_df) || base::nrow(slot_df) == 0) {
      next
    }
    if (!all(c("cluster", "cell_type", "count") %in% base::colnames(slot_df))) {
      next
    }
    slot_df$cluster <- base::as.character(slot_df$cluster)
    slot_df$cell_type <- base::as.character(slot_df$cell_type)
    slot_df$count <- .hc_as_numeric_safely(slot_df$count)
    slot_df$count[!base::is.finite(slot_df$count)] <- 0

    slot_hidden <- base::attr(slot_df, "hidden")
    slot_label <- slot_nm
    if (base::is.list(slot_hidden)) {
      if ("database" %in% base::names(slot_hidden) &&
        base::length(slot_hidden$database) > 0 &&
        base::nzchar(base::as.character(slot_hidden$database[[1]]))) {
        slot_label <- base::as.character(slot_hidden$database[[1]])
      }
      if ("label" %in% base::names(slot_hidden) &&
        base::length(slot_hidden$label) > 0 &&
        base::nzchar(base::as.character(slot_hidden$label[[1]]))) {
        slot_label <- base::as.character(slot_hidden$label[[1]])
      }
    }

    user_enrichment_slots[[slot_nm]] <- list(
      slot = slot_nm,
      label = slot_label,
      data = slot_df
    )
  }

  # Get categorical and numerical sample group annotations if they exist:
  if ("column_annos_categorical" %in% base::names(hcobject[["satellite_outputs"]])) {
    column_anno_categorical <- hcobject[["satellite_outputs"]][["column_annos_categorical"]]
  } else {
    column_anno_categorical <- NULL
  }

  if ("column_annos_numerical" %in% base::names(hcobject[["satellite_outputs"]])) {
    column_anno_numerical <- hcobject[["satellite_outputs"]][["column_annos_numerical"]]
  } else {
    column_anno_numerical <- NULL
  }

  if (base::is.null(cat_as_bp)) {
    if (!base::is.null(column_anno_categorical)) {
      cat_as_bp <- base::rep(FALSE, base::length(column_anno_categorical))
    }
  }

  base::gc()

  parse_gene_string <- function(x) {
    if (base::is.null(x) || base::length(x) == 0 || base::all(base::is.na(x))) {
      return(base::character(0))
    }
    genes <- base::trimws(
      base::unlist(
        base::strsplit(base::as.character(x), ",", fixed = TRUE),
        use.names = FALSE
      )
    )
    genes <- genes[!base::is.na(genes) & genes != ""]
    base::unique(genes)
  }

  build_cluster_gene_map <- function(cluster_df, target_clusters) {
    gene_strings <- base::split(base::as.character(cluster_df$gene_n), base::as.character(cluster_df$color))
    out <- stats::setNames(base::vector("list", base::length(target_clusters)), target_clusters)
    for (cl in target_clusters) {
      out[[cl]] <- parse_gene_string(gene_strings[[cl]])
    }
    out
  }

  build_mat_heatmap_fast <- function(gfc_df, cluster_gene_map, target_clusters) {
    if (!("Gene" %in% base::colnames(gfc_df))) {
      stop("`GFC_all_layers` must contain a `Gene` column.")
    }

    gfc_cols <- base::seq_len(base::ncol(gfc_df) - 1)
    gfc_mat <- base::as.matrix(gfc_df[, gfc_cols, drop = FALSE])
    storage.mode(gfc_mat) <- "numeric"
    gene_ids <- base::as.character(gfc_df$Gene)
    valid_gene_ids <- !base::is.na(gene_ids) & base::nzchar(gene_ids)
    gfc_mat <- gfc_mat[valid_gene_ids, , drop = FALSE]
    gene_ids <- gene_ids[valid_gene_ids]
    base::rownames(gfc_mat) <- gene_ids
    gene_index <- base::split(base::seq_along(gene_ids), gene_ids)

    out <- base::vapply(
      target_clusters,
      function(cl) {
        genes <- cluster_gene_map[[cl]]
        if (base::length(genes) == 0) {
          return(base::rep(NA_real_, base::ncol(gfc_mat)))
        }
        idx <- base::unlist(gene_index[genes], use.names = FALSE)
        if (base::length(idx) == 0) {
          return(base::rep(NA_real_, base::ncol(gfc_mat)))
        }
        base::colMeans(gfc_mat[idx, , drop = FALSE], na.rm = TRUE)
      },
      FUN.VALUE = base::numeric(base::ncol(gfc_mat))
    )

    out <- base::t(out)
    base::rownames(out) <- target_clusters
    base::colnames(out) <- base::colnames(gfc_df)[gfc_cols]
    out
  }

  # Filter for included clusters (non-white)
  c_df <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], cluster_included == "yes")
  available_row_colors <- base::unique(base::as.character(c_df$color))
  if (!base::is.null(row_order)) {
    row_order <- .hc_resolve_cluster_heatmap_row_order(
      row_order = row_order,
      cluster_calc = hcobject[["integrated_output"]][["cluster_calc"]],
      available_colors = available_row_colors,
      module_prefix = module_prefix
    )
  }

  # --- 2. Build Heatmap Matrix (GFC Values) ---
  target_clusters <- if (!base::is.null(row_order)) row_order else base::unique(c_df$color)
  target_clusters <- base::as.character(target_clusters)
  cluster_gene_map <- build_cluster_gene_map(c_df, target_clusters)
  mat_heatmap <- build_mat_heatmap_fast(
    gfc_df = hcobject[["integrated_output"]][["GFC_all_layers"]],
    cluster_gene_map = cluster_gene_map,
    target_clusters = target_clusters
  )

  existing_heatmap_col_order <- .hc_heatmap_cache_info(
    hcobject[["integrated_output"]][["cluster_calc"]]
  )$col_order
  if (!isTRUE(cluster_columns)) {
    selected_col_order <- .hc_select_heatmap_col_order(
      available_cols = base::colnames(mat_heatmap),
      plot_order = col_order,
      main_order = existing_heatmap_col_order,
      fallback_order = base::colnames(mat_heatmap),
      context = "cluster heatmap"
    )
    if (base::length(selected_col_order) > 0) {
      mat_heatmap <- .hc_subset_matrix_cols_with_duplicates(
        mat_heatmap,
        selected_col_order
      ) %>% base::as.matrix()
    }
  }

  # --- 3. Prepare Enrichment Data ---

  enrichment_entries <- list()
  .truncate_for_label <- function(x, max_chars = 34) {
    x <- base::as.character(x)
    x[base::is.na(x)] <- ""
    need_cut <- base::nchar(x) > max_chars
    x[need_cut] <- base::paste0(base::substr(x[need_cut], 1, max_chars - 3), "...")
    x
  }
  .clean_celltype_label <- function(x) {
    x <- base::as.character(x)
    x[base::is.na(x)] <- ""
    x <- base::trimws(x)
    x <- gsub("^\\[[^\\]]+\\]\\s*", "", x, perl = TRUE)
    x[x %in% c("", "NA", "N/A", "na", "n/a")] <- ""
    x
  }

  # Helper logic to extract enrichment data
  if (base::length(user_enrichment_slots) > 0) {
    base_palette <- ggsci::pal_d3(palette = "category20")(20)
    for (slot_idx in base::seq_along(user_enrichment_slots)) {
      slot_item <- user_enrichment_slots[[slot_idx]]
      slot_df <- slot_item$data
      slot_df <- slot_df[
        !base::is.na(slot_df$cluster) & base::nzchar(base::trimws(base::as.character(slot_df$cluster))) &
          !base::is.na(slot_df$cell_type), ,
        drop = FALSE
      ]
      slot_df$cell_type <- .clean_celltype_label(slot_df$cell_type)
      slot_df <- slot_df[base::nzchar(slot_df$cell_type), , drop = FALSE]
      if (base::nrow(slot_df) == 0) {
        next
      }
      cell_types <- base::unique(base::as.character(slot_df$cell_type))
      if (base::length(cell_types) == 0) {
        next
      }

      slot_mat <- base::matrix(
        0,
        nrow = base::length(target_clusters),
        ncol = base::length(cell_types),
        dimnames = list(target_clusters, cell_types)
      )
      slot_hits <- stats::setNames(base::rep(0, base::length(target_clusters)), target_clusters)

      for (x in target_clusters) {
        tmp <- dplyr::filter(slot_df, cluster == x)
        if (base::nrow(tmp) == 0) {
          next
        }
        tmp_agg <- stats::aggregate(count ~ cell_type, data = tmp, FUN = base::sum)
        idx_all <- base::match(tmp_agg$cell_type, cell_types)
        keep_idx <- base::is.finite(idx_all)
        if (base::any(keep_idx)) {
          slot_mat[x, idx_all[keep_idx]] <- as.numeric(tmp_agg$count[keep_idx])
        }
        if ("hits" %in% base::colnames(tmp)) {
          h <- .hc_first_numeric_value(dplyr::first(tmp$hits))
          if (base::is.finite(h)) {
            slot_hits[[x]] <- h
          } else {
            slot_hits[[x]] <- base::sum(tmp_agg$count, na.rm = TRUE)
          }
        } else {
          slot_hits[[x]] <- base::sum(tmp_agg$count, na.rm = TRUE)
        }
      }

      if (base::all(slot_mat == 0)) {
        next
      }

      slot_mat_full <- slot_mat
      if (!base::is.null(celltype_bar_top_n) && base::ncol(slot_mat) > celltype_bar_top_n) {
        col_totals <- base::colSums(slot_mat, na.rm = TRUE)
        keep_idx <- base::order(col_totals, decreasing = TRUE)
        keep_idx <- keep_idx[base::seq_len(base::min(base::length(keep_idx), celltype_bar_top_n))]
        slot_mat <- slot_mat[, keep_idx, drop = FALSE]
        other_idx <- base::setdiff(base::seq_len(base::ncol(slot_mat_full)), keep_idx)
        if (isTRUE(celltype_bar_include_other) && base::length(other_idx) > 0) {
          other_vals <- base::rowSums(slot_mat_full[, other_idx, drop = FALSE], na.rm = TRUE)
          if (celltype_bar_other_label %in% base::colnames(slot_mat)) {
            slot_mat[, celltype_bar_other_label] <- slot_mat[, celltype_bar_other_label] + other_vals
          } else {
            slot_mat <- base::cbind(slot_mat, .other = other_vals)
            base::colnames(slot_mat)[base::ncol(slot_mat)] <- celltype_bar_other_label
          }
        }
      }

      dominant_text <- stats::setNames(base::rep("", base::nrow(slot_mat_full)), base::rownames(slot_mat_full))
      for (cl in base::rownames(slot_mat_full)) {
        v <- slot_mat_full[cl, , drop = TRUE]
        if (!base::any(base::is.finite(v) & v > 0)) {
          dominant_text[[cl]] <- ""
          next
        }
        ct <- .clean_celltype_label(base::colnames(slot_mat_full))
        is_other_like <- base::tolower(ct) == base::tolower(celltype_bar_other_label)
        idx_pool <- base::which(!is_other_like & base::is.finite(v) & v > 0)
        if (base::length(idx_pool) == 0) {
          idx_pool <- base::which(base::is.finite(v) & v > 0)
        }
        if (base::length(idx_pool) == 0) {
          dominant_text[[cl]] <- ""
          next
        }
        idx <- idx_pool[[base::which.max(v[idx_pool])]]
        dom_ct <- .clean_celltype_label(base::colnames(slot_mat_full)[[idx]])
        dom_pct <- .hc_first_numeric_value(v[[idx]])
        if (!base::nzchar(dom_ct) || !base::is.finite(dom_pct)) {
          dominant_text[[cl]] <- ""
        } else {
          dominant_text[[cl]] <- base::paste0(dom_ct, " (", base::format(round(dom_pct, 1), trim = TRUE, nsmall = 0), "%)")
        }
      }
      dominant_text <- .truncate_for_label(dominant_text, max_chars = 34)

      pal_offset <- ((slot_idx - 1) * 3) %% base::length(base_palette)
      ct_names <- base::colnames(slot_mat)
      is_other <- ct_names == celltype_bar_other_label
      slot_cols <- base::character(base::length(ct_names))
      n_main <- base::sum(!is_other)
      if (n_main > 0) {
        slot_cols[!is_other] <- base_palette[((base::seq_len(n_main) + pal_offset - 1) %% base::length(base_palette)) + 1]
      }
      if (base::any(is_other)) {
        slot_cols[is_other] <- "#bdbdbd"
      }
      enrichment_entries[[base::length(enrichment_entries) + 1]] <- list(
        slot = slot_item$slot,
        label = slot_item$label,
        mat = slot_mat,
        hits = slot_hits,
        cell_types = base::colnames(slot_mat),
        colors = slot_cols,
        dominant_text = base::unname(dominant_text[base::as.character(target_clusters)])
      )
    }
  }

  if (base::length(enrichment_entries) > 0) {
    raw_labels <- base::vapply(
      enrichment_entries,
      function(x) base::as.character(x$label)[1],
      FUN.VALUE = base::character(1)
    )
    uniq_labels <- base::make.unique(raw_labels, sep = "_")
    for (i in base::seq_along(enrichment_entries)) {
      enrichment_entries[[i]]$label <- uniq_labels[[i]]
    }
  }


  # --- 4. Setup Row Annotations (Colors, Labels, Gene Counts) ---

  # Re-sort c_df to match row_order if necessary
  if (!base::is.null(row_order)) {
    cluster_colors <- base::factor(row_order)
    base::names(cluster_colors) <- row_order
    c_df <- c_df[base::match(row_order, c_df$color), ]
  } else {
    cluster_colors <- base::factor(c_df$color)
    base::names(cluster_colors) <- c_df$color
    row_order <- base::unique(c_df$color)
  }

  # Optional module-significance labels from `.hc_module_condition_significance_driver()`.
  module_sig_labels <- NULL
  module_sig_q <- NULL
  module_sig_method_used <- NULL
  module_sig_table <- NULL
  if (isTRUE(include_module_significance)) {
    sig_slot <- satellite_outputs[[module_significance_slot]]
    if (!base::is.null(sig_slot) && base::is.list(sig_slot)) {
      .sig_star_fun <- function(qvals, cuts) {
        out <- base::rep("", base::length(qvals))
        ok <- base::is.finite(qvals)
        out[ok & qvals <= cuts[[1]]] <- "***"
        out[ok & qvals > cuts[[1]] & qvals <= cuts[[2]]] <- "**"
        out[ok & qvals > cuts[[2]] & qvals <= cuts[[3]]] <- "*"
        out
      }

      candidate_names <- if (base::identical(module_significance_method, "auto")) {
        base::c("wilcox", "limma_summary", "lmm")
      } else if (base::identical(module_significance_method, "limma")) {
        "limma_summary"
      } else {
        module_significance_method
      }

      for (nm in candidate_names) {
        if (!(nm %in% base::names(sig_slot))) {
          next
        }
        sig_df <- sig_slot[[nm]]
        if (base::is.null(sig_df) || !base::is.data.frame(sig_df) || base::nrow(sig_df) == 0) {
          next
        }
        if (!("cluster" %in% base::colnames(sig_df)) && "module" %in% base::colnames(sig_df)) {
          sig_df$cluster <- sig_df$module
        }
        q_col <- base::intersect(
          base::c("p_adj", "adj.P.Val", "qvalue", "q", "best_q"),
          base::colnames(sig_df)
        )
        if (!("cluster" %in% base::colnames(sig_df)) || base::length(q_col) == 0) {
          next
        }
        sig_df$cluster <- base::as.character(sig_df$cluster)
        sig_df$q_tmp <- .hc_as_numeric_safely(sig_df[[q_col[[1]]]])
        sig_df <- sig_df[!base::is.na(sig_df$cluster) & base::nzchar(sig_df$cluster), , drop = FALSE]
        if (base::nrow(sig_df) == 0) {
          next
        }
        split_sig <- base::split(sig_df, sig_df$cluster)
        best_sig <- base::lapply(split_sig, function(x) {
          qv <- .hc_as_numeric_safely(x$q_tmp)
          idx <- if (base::any(base::is.finite(qv))) {
            base::which.min(base::ifelse(base::is.finite(qv), qv, Inf))
          } else {
            1
          }
          x[idx, , drop = FALSE]
        })
        best_sig <- base::do.call(base::rbind, best_sig)
        hit <- base::match(row_order, best_sig$cluster)
        module_sig_q <- best_sig$q_tmp[hit]
        base::names(module_sig_q) <- row_order
        stars <- .sig_star_fun(module_sig_q, module_significance_p_cutoffs)
        module_sig_labels <- if (isTRUE(module_significance_show_qvalue)) {
          base::ifelse(
            base::is.finite(module_sig_q),
            base::ifelse(
              stars == "",
              base::paste0("q=", base::formatC(module_sig_q, format = "fg", digits = 2)),
              base::paste0(stars, " q=", base::formatC(module_sig_q, format = "fg", digits = 2))
            ),
            ""
          )
        } else {
          stars
        }
        module_sig_method_used <- if (base::identical(nm, "limma_summary")) "limma" else nm
        break
      }

      if (base::is.null(module_sig_labels) &&
        "summary" %in% base::names(sig_slot) &&
        base::is.data.frame(sig_slot[["summary"]]) &&
        base::nrow(sig_slot[["summary"]]) > 0) {
        sum_df <- sig_slot[["summary"]]
        if (!("cluster" %in% base::colnames(sum_df)) && "module" %in% base::colnames(sum_df)) {
          sum_df$cluster <- sum_df$module
        }
        if ("cluster" %in% base::colnames(sum_df)) {
          sum_df$cluster <- base::as.character(sum_df$cluster)
          hit <- base::match(row_order, sum_df$cluster)
          q_cols <- if (base::identical(module_significance_method, "auto")) {
            base::c("best_q", "wilcox_q", "limma_q", "lmm_q", "qvalue")
          } else {
            base::c(
              base::paste0(module_significance_method, "_q"),
              "best_q",
              "qvalue"
            )
          }
          sig_cols <- if (base::identical(module_significance_method, "auto")) {
            base::c("best_sig", "wilcox_sig", "limma_sig", "lmm_sig", "significance")
          } else {
            base::c(
              base::paste0(module_significance_method, "_sig"),
              "best_sig",
              "significance"
            )
          }
          q_col <- base::intersect(q_cols, base::colnames(sum_df))
          sig_col <- base::intersect(sig_cols, base::colnames(sum_df))
          if (base::length(q_col) > 0) {
            module_sig_q <- .hc_as_numeric_safely(sum_df[[q_col[[1]]]])[hit]
            base::names(module_sig_q) <- row_order
          }
          if (base::length(sig_col) > 0) {
            module_sig_labels <- base::as.character(sum_df[[sig_col[[1]]]])[hit]
          } else if (!base::is.null(module_sig_q)) {
            module_sig_labels <- .sig_star_fun(module_sig_q, module_significance_p_cutoffs)
          }
          if ("best_method" %in% base::colnames(sum_df)) {
            method_vals <- base::as.character(sum_df$best_method)[hit]
            method_vals <- method_vals[!base::is.na(method_vals) & base::nzchar(method_vals)]
            if (base::length(method_vals) > 0) {
              module_sig_method_used <- base::paste(base::unique(method_vals), collapse = ", ")
            }
          }
        }
      }
    }

    if (!base::is.null(module_sig_labels)) {
      module_sig_labels <- base::as.character(module_sig_labels)
      module_sig_labels[base::is.na(module_sig_labels)] <- ""
      if (base::is.null(module_sig_q) || base::length(module_sig_q) != base::length(row_order)) {
        module_sig_q <- base::rep(NA_real_, base::length(row_order))
      }
      module_sig_table <- base::data.frame(
        cluster = row_order,
        method = if (base::is.null(module_sig_method_used)) "unknown" else module_sig_method_used,
        qvalue = .hc_as_numeric_safely(module_sig_q),
        significance = module_sig_labels,
        stringsAsFactors = FALSE
      )
      .hc_set_bridge_hcobject_slot(c("satellite_outputs", "module_significance_last_heatmap"), module_sig_table)
      message(
        "Module-significance summary table (method: ",
        if (base::is.null(module_sig_method_used)) "unknown" else module_sig_method_used,
        "):"
      )
      .hc_display_object(module_sig_table, row.names = FALSE)
    } else {
      warning(
        "No module-significance labels could be drawn. Run `hc_module_condition_significance()` first ",
        "and check slot `", module_significance_slot, "`.",
        call. = FALSE
      )
    }
  }

  has_enrichment <- base::length(enrichment_entries) > 0
  show_gene_bar <- FALSE
  show_gene_text <- FALSE
  if (gene_count_mode == "legacy") {
    show_gene_bar <- TRUE
    show_gene_text <- !has_enrichment
  } else if (gene_count_mode == "bar_and_text") {
    show_gene_bar <- TRUE
    show_gene_text <- TRUE
  } else if (gene_count_mode == "bar") {
    show_gene_bar <- TRUE
  } else if (gene_count_mode == "text") {
    show_gene_text <- TRUE
  }

  clustered_row_order <- row_order
  cluster_rows_for_heatmap <- cluster_rows
  if (isTRUE(cluster_rows) && base::nrow(c_df) > 1) {
    row_cluster_info <- tryCatch(
      {
        hc_rows <- stats::hclust(stats::dist(mat_heatmap), method = "complete")
        row_dend <- stats::as.dendrogram(hc_rows)
        # Match ComplexHeatmap default dendrogram reordering (row_dend_reorder = TRUE).
        row_weights <- base::rowMeans(mat_heatmap, na.rm = TRUE)
        row_weights <- row_weights[base::rownames(mat_heatmap)]
        row_dend <- stats::reorder(row_dend, row_weights, agglo.FUN = base::mean)
        # Use an explicitly flipped dendrogram so the plotted top-to-bottom order
        # follows `order.dendrogram(dend_flipped)` (M1 at the top).
        dend_flipped <- base::rev(row_dend)
        list(
          order = row_order[stats::order.dendrogram(dend_flipped)],
          dend = dend_flipped
        )
      },
      error = function(e) {
        hc_rows_fallback <- tryCatch(
          stats::hclust(stats::dist(mat_heatmap), method = "complete"),
          error = function(e2) NULL
        )
        if (base::is.null(hc_rows_fallback)) {
          list(
            order = row_order,
            dend = cluster_rows
          )
        } else {
          dend_fallback <- base::rev(stats::as.dendrogram(hc_rows_fallback))
          list(
            order = row_order[stats::order.dendrogram(dend_fallback)],
            dend = dend_fallback
          )
        }
      }
    )
    clustered_row_order <- row_cluster_info$order
    cluster_rows_for_heatmap <- row_cluster_info$dend
  }
  displayed_row_order <- if (inherits(cluster_rows_for_heatmap, "dendrogram")) {
    row_order[stats::order.dendrogram(cluster_rows_for_heatmap)]
  } else {
    clustered_row_order
  }
  module_labels <- NULL
  if (module_label_mode == "prefix") {
    existing_label_map <- hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]]
    map_complete <- FALSE
    preserved_labels <- NULL
    if (!base::is.null(existing_label_map) && base::length(existing_label_map) > 0) {
      existing_label_map <- base::as.character(existing_label_map)
      map_names <- base::names(hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]])
      if (!base::is.null(map_names) && base::length(map_names) == base::length(existing_label_map)) {
        base::names(existing_label_map) <- base::as.character(map_names)
      }
      # Recover from accidentally inverted map orientation (labels->colors).
      missing_before <- base::setdiff(row_order, base::names(existing_label_map))
      if (base::length(missing_before) > 0) {
        inverse_map <- stats::setNames(base::names(existing_label_map), base::as.character(existing_label_map))
        if (base::all(row_order %in% base::names(inverse_map))) {
          existing_label_map <- inverse_map
        }
      }
      preserved_labels <- base::as.character(existing_label_map[row_order])
      map_complete <- base::length(preserved_labels) == base::length(row_order) &&
        base::all(!base::is.na(preserved_labels) & base::nzchar(preserved_labels))
    }

    base_labels <- if (isTRUE(map_complete)) {
      preserved_labels
    } else {
      base::paste0(module_prefix, base::seq_len(base::nrow(c_df)))
    }

    if (module_label_numbering == "preserve_existing") {
      module_labels <- base_labels
    } else {
      numbering_order <- if (module_label_numbering == "after_clustering" &&
        isTRUE(cluster_rows) &&
        base::nrow(c_df) > 1) {
        displayed_row_order
      } else {
        row_order
      }
      label_by_cluster <- stats::setNames(base::as.character(base_labels), row_order)
      labels_in_order <- base::as.character(label_by_cluster[numbering_order])
      has_split_like_labels <- base::any(base::grepl("^[^.]+\\.[0-9]+$", labels_in_order))

      if (isTRUE(has_split_like_labels)) {
        split_mask <- base::grepl("^[^.]+\\.[0-9]+$", labels_in_order)
        parent_labels_raw <- base::vapply(
          labels_in_order,
          function(lbl) {
            if (base::grepl("^[^.]+\\.[0-9]+$", lbl)) {
              base::sub("\\.[0-9]+$", "", lbl)
            } else {
              lbl
            }
          },
          FUN.VALUE = base::character(1)
        )
        parent_labels_in_order <- parent_labels_raw
        parent_labels_in_order <- base::unique(parent_labels_in_order)
        parent_remap <- stats::setNames(
          base::paste0(module_prefix, base::seq_along(parent_labels_in_order)),
          parent_labels_in_order
        )
        counters <- base::list()
        for (ii in base::seq_along(labels_in_order)) {
          current_label <- labels_in_order[[ii]]
          parent_label_old <- parent_labels_raw[[ii]]
          parent_label_new <- if (parent_label_old %in% base::names(parent_remap)) {
            base::as.character(parent_remap[[parent_label_old]])
          } else {
            parent_label_old
          }
          if (split_mask[[ii]]) {
            if (!(parent_label_new %in% base::names(counters))) {
              counters[[parent_label_new]] <- 0L
            }
            counters[[parent_label_new]] <- counters[[parent_label_new]] + 1L
            labels_in_order[[ii]] <- base::paste0(parent_label_new, ".", counters[[parent_label_new]])
          } else {
            labels_in_order[[ii]] <- parent_label_new
          }
        }
        remap <- stats::setNames(labels_in_order, numbering_order)
        module_labels <- base::as.character(remap[row_order])
      } else {
        remap <- stats::setNames(
          base::paste0(module_prefix, base::seq_along(numbering_order)),
          numbering_order
        )
        module_labels <- base::as.character(remap[row_order])
      }
    }
  } else if (module_label_mode == "color") {
    module_labels <- base::as.character(row_order)
  }

  module_sig_integrated <- (
    isTRUE(include_module_significance) &&
      !base::is.null(module_sig_labels) &&
      !isTRUE(module_significance_show_qvalue) &&
      !base::is.null(module_labels) &&
      base::identical(module_label_mode, "prefix")
  )
  sig_suffix <- NULL
  module_labels_display <- if (isTRUE(module_sig_integrated)) {
    sig_suffix <- base::ifelse(
      base::is.na(module_sig_labels) | !base::nzchar(module_sig_labels),
      "",
      base::as.character(module_sig_labels)
    )
    base::paste0(module_labels, sig_suffix)
  } else {
    module_labels
  }
  module_labels_have_split_suffix <- !base::is.null(module_labels_display) &&
    base::length(module_labels_display) > 0 &&
    base::any(base::grepl("\\.[0-9]+", module_labels_display))
  max_sig_stars <- if (isTRUE(module_sig_integrated) && !base::is.null(sig_suffix)) {
    base::max(base::nchar(sig_suffix), na.rm = TRUE)
  } else {
    0
  }

  max_chars <- if (base::is.null(module_labels_display)) {
    1
  } else {
    base::max(base::nchar(module_labels_display), na.rm = TRUE)
  }

  n_heat_rows <- base::nrow(mat_heatmap)
  n_heat_cols <- base::ncol(mat_heatmap)
  duplicate_condition_width_scale <- .hc_gfc_duplicate_condition_width_scale(
    hcobject,
    base::colnames(mat_heatmap)
  )
  n_rows <- base::nrow(c_df)
  preset_scale_font <- switch(module_label_preset,
    compact = 0.82,
    balanced = 1.0,
    presentation = 1.2,
    auto = 1.0
  )
  preset_scale_width <- switch(module_label_preset,
    compact = 0.95,
    balanced = 1.0,
    presentation = 1.15,
    auto = 1.0
  )
  preset_scale_pt <- switch(module_label_preset,
    compact = 0.85,
    balanced = 1.0,
    presentation = 1.15,
    auto = 1.0
  )

  if (!user_set_module_label_fontsize) {
    # Set a fixed, readable font size for module labels.
    module_label_fontsize <- 5
  }

  if (!user_set_module_box_width_cm) {
    if (base::is.null(module_labels_display) || module_label_mode == "legacy") {
      base_width <- 0.5
    } else {
      # Keep the default layout unchanged; split labels and significance
      # suffixes are widened later in `.hc_module_label_draw_width_cm()`.
      base_width <- (max_chars * 0.09) + 0.15
    }
    module_box_width_cm <- base::max(
      0.80,
      base::min(
        4.8,
        base_width * preset_scale_width
      )
    )
  }
  module_box_width_cm_draw <- .hc_module_label_draw_width_cm(
    module_box_width_cm = module_box_width_cm,
    module_labels_display = module_labels_display,
    user_set_module_box_width_cm = user_set_module_box_width_cm,
    module_sig_integrated = module_sig_integrated,
    max_sig_stars = max_sig_stars
  )
  cell_size_mm <- .hc_cluster_heatmap_cell_size_mm(
    n_heat_rows = n_heat_rows,
    n_heat_cols = n_heat_cols,
    duplicate_condition_width_scale = duplicate_condition_width_scale,
    module_box_width_cm_draw = module_box_width_cm_draw,
    overall_plot_scale = overall_plot_scale
  )

  if (!user_set_module_label_pt_size) {
    base_pt <- if (n_rows <= 10) {
      0.95
    } else if (n_rows <= 15) {
      0.90
    } else if (n_rows <= 25) {
      0.82
    } else if (n_rows <= 40) {
      0.72
    } else {
      0.62
    }
    # Request a large default glyph; the box-fit guard chooses the largest
    # common rendered size that still fits every module label.
    char_penalty <- 1
    module_label_pt_size <- base::max(
      0.14,
      base::min(
        1.05,
        (base_pt * preset_scale_pt) / char_penalty
      )
    )
  }
  module_label_pt_size_draw <- module_label_pt_size
  if (!user_set_module_label_pt_size && isTRUE(module_sig_integrated) && max_sig_stars > 0) {
    sig_pt_scale <- base::max(0.64, 1 - (0.11 * max_sig_stars))
    module_label_pt_size_draw <- base::max(
      0.14,
      base::min(0.55, module_label_pt_size * sig_pt_scale)
    )
  }
  module_label_fit <- .hc_module_label_fit_pt(
    module_label_pt_size = module_label_pt_size_draw,
    module_box_width_cm = module_box_width_cm_draw,
    module_labels_display = module_labels_display,
    n_heat_rows = n_heat_rows,
    cell_size_mm = cell_size_mm,
    module_label_fontsize = module_label_fontsize,
    use_fontsize_request = user_set_module_label_fontsize && !user_set_module_label_pt_size,
    fontface = "bold"
  )
  module_label_pt_size_unit <- if (base::length(module_labels_display) > 0 &&
    base::length(module_label_fit$pt_size) == base::length(module_labels_display)) {
    grid::unit(module_label_fit$pt_size, "pt")
  } else {
    grid::unit(module_label_pt_size_draw, "snpc")
  }
  module_label_fontsize_draw <- .hc_module_label_effective_fontsize(
    module_label_fit = module_label_fit,
    fallback_fontsize = module_label_fontsize,
    fallback_pt_size = module_label_pt_size_draw
  )

  if (!user_set_gene_count_fontsize) {
    gene_count_fontsize <- base::max(
      5,
      base::min(10, module_label_fontsize_draw * 0.72)
    )
  }

  if (!user_set_gene_count_pt_size) {
    if (gene_count_renderer == "pch") {
      # Independent default for gene counts: smaller and stable across module label changes
      gene_count_pt_size <- if (n_rows <= 10) {
        0.25
      } else if (n_rows <= 15) {
        0.22
      } else if (n_rows <= 25) {
        0.20
      } else if (n_rows <= 40) {
        0.18
      } else {
        0.15
      }
    } else {
      gene_count_pt_size <- 0.5
    }
  }

  shared_heatmap_line_lwd <- 0.5
  module_box_border_gp <- grid::gpar(col = "black", lwd = shared_heatmap_line_lwd)
  heatmap_cell_border_gp <- grid::gpar(col = "black", lwd = shared_heatmap_line_lwd)
  dendrogram_line_gp <- grid::gpar(col = "black", lwd = shared_heatmap_line_lwd)

  show_module_sig <- isTRUE(include_module_significance) &&
    !base::is.null(module_sig_labels) &&
    !isTRUE(module_sig_integrated)
  module_sig_width_cm_use <- if (isTRUE(show_module_sig)) {
    if (isTRUE(module_significance_show_qvalue)) {
      base::max(2.4, module_significance_width_cm)
    } else {
      base::max(1.8, module_significance_width_cm)
    }
  } else {
    module_significance_width_cm
  }
  module_box_anno <- .hc_module_label_box_annotation(
    values = row_order,
    colors = cluster_colors,
    labels = module_labels_display,
    label_color = module_label_color,
    label_fontsize_pt = module_label_fit$base_pt_size,
    fontface = "bold",
    width_cm = module_box_width_cm_draw,
    border_gp = module_box_border_gp,
    which = "row"
  )

  gene_bar_anno <- if (show_gene_bar) {
    ComplexHeatmap::anno_barplot(c_df$gene_no, width = grid::unit(2.5, "cm"), which = "row")
  } else {
    ComplexHeatmap::anno_empty(width = grid::unit(0, "mm"), which = "row", border = FALSE)
  }

  gene_text_anno <- if (show_gene_text) {
    if (gene_count_renderer == "pch") {
      ComplexHeatmap::anno_simple(
        x = base::rep("count_text", base::nrow(c_df)),
        col = c(count_text = "transparent"),
        pch = base::as.character(c_df$gene_no),
        pt_gp = grid::gpar(col = "black", fontsize = gene_count_fontsize, fontface = gene_count_fontface),
        pt_size = grid::unit(gene_count_pt_size, "snpc"),
        gp = grid::gpar(col = NA),
        simple_anno_size = grid::unit(1.2, "cm"),
        which = "row"
      )
    } else {
      ComplexHeatmap::anno_text(
        base::as.character(c_df$gene_no),
        width = grid::unit(1.2, "cm"),
        just = "left",
        gp = grid::gpar(fontsize = gene_count_fontsize, fontface = gene_count_fontface),
        which = "row"
      )
    }
  } else {
    ComplexHeatmap::anno_empty(width = grid::unit(0, "mm"), which = "row", border = FALSE)
  }

  module_sig_anno <- if (isTRUE(show_module_sig)) {
    sig_cols <- if (!base::is.null(module_sig_q) && base::length(module_sig_q) == base::length(module_sig_labels)) {
      base::ifelse(
        base::is.finite(module_sig_q) & module_sig_q <= module_significance_p_cutoffs[[3]],
        "black",
        "#8a8a8a"
      )
    } else {
      base::ifelse(base::nzchar(module_sig_labels), "black", "#8a8a8a")
    }
    if (isTRUE(module_significance_show_qvalue)) {
      ComplexHeatmap::anno_text(
        module_sig_labels,
        width = grid::unit(module_sig_width_cm_use, "cm"),
        just = "center",
        gp = grid::gpar(
          col = sig_cols,
          fontsize = base::max(9.5, gene_count_fontsize),
          fontface = "bold"
        ),
        which = "row"
      )
    } else {
      ComplexHeatmap::anno_simple(
        x = base::rep("sig_bg", base::length(module_sig_labels)),
        col = c(sig_bg = "#f5f5f5"),
        pch = module_sig_labels,
        pt_gp = grid::gpar(
          col = sig_cols,
          fontsize = base::max(11, gene_count_fontsize + 1),
          fontface = "bold"
        ),
        pt_size = grid::unit(0.9, "snpc"),
        simple_anno_size = grid::unit(module_sig_width_cm_use, "cm"),
        gp = grid::gpar(col = "#d0d0d0"),
        which = "row"
      )
    }
  } else {
    ComplexHeatmap::anno_empty(width = grid::unit(0, "mm"), which = "row", border = FALSE)
  }

  base_row_width_cm <- module_box_width_cm_draw +
    if (show_gene_text) {
      1.2
    } else {
      0 +
        if (show_gene_bar) {
          2.5
        } else {
          0 +
            if (show_module_sig) {
              module_sig_width_cm_use
            } else {
              0 +
                0.8
            }
        }
    }
  enrichment_slot_bar_width_cm <- if (show_celltype_bars) 2.2 else 0
  enrichment_slot_text_width_cm <- if (show_celltype_text) celltype_bar_dominant_width_cm else 0
  enrichment_slot_gap_cm <- 0.15
  n_enrichment_slots <- base::length(enrichment_entries)
  enrichment_total_width_cm <- if (n_enrichment_slots > 0) {
    (n_enrichment_slots * (enrichment_slot_bar_width_cm + enrichment_slot_text_width_cm)) +
      ((n_enrichment_slots - 1) * enrichment_slot_gap_cm)
  } else {
    0
  }
  row_annotation_total_width_cm <- base_row_width_cm + enrichment_total_width_cm

  lgd_list <- list()
  ha_annos <- list(
    modules = module_box_anno,
    `# genes` = gene_text_anno,
    genes = gene_bar_anno
  )
  if (isTRUE(show_module_sig)) {
    sig_anno_name <- module_significance_annotation_name
    if (sig_anno_name %in% base::names(ha_annos)) {
      sig_anno_name <- base::paste0(sig_anno_name, "_sig")
    }
    ha_annos[[sig_anno_name]] <- module_sig_anno
  }
  if (n_enrichment_slots > 0) {
    for (i in base::seq_along(enrichment_entries)) {
      entry <- enrichment_entries[[i]]
      anno_name <- entry$label
      if (show_celltype_bars) {
        ha_annos[[anno_name]] <- ComplexHeatmap::anno_barplot(
          entry$mat,
          width = grid::unit(enrichment_slot_bar_width_cm, "cm"),
          gp = grid::gpar(fill = entry$colors, col = entry$colors),
          baseline = 0,
          which = "row"
        )
      }
      if (show_celltype_text) {
        dominant_name <- base::paste0(anno_name, "_dominant")
        dom_mat <- entry$mat
        dom_labels <- base::character(base::nrow(dom_mat))
        for (ri in base::seq_len(base::nrow(dom_mat))) {
          vv <- .hc_as_numeric_safely(dom_mat[ri, , drop = TRUE])
          if (!base::any(base::is.finite(vv) & vv > 0)) {
            dom_labels[[ri]] <- ""
            next
          }
          ctn <- .clean_celltype_label(base::colnames(dom_mat))
          is_other_like <- base::tolower(ctn) == base::tolower(celltype_bar_other_label)
          idx_pool <- base::which(!is_other_like & base::is.finite(vv) & vv > 0)
          if (base::length(idx_pool) == 0) {
            idx_pool <- base::which(base::is.finite(vv) & vv > 0)
          }
          if (base::length(idx_pool) == 0) {
            dom_labels[[ri]] <- ""
            next
          }
          j <- idx_pool[[base::which.max(vv[idx_pool])]]
          ct_lab <- ctn[[j]]
          ct_val <- vv[[j]]
          if (!base::nzchar(ct_lab) || !base::is.finite(ct_val)) {
            dom_labels[[ri]] <- ""
          } else {
            dom_labels[[ri]] <- base::paste0(ct_lab, " (", base::format(round(ct_val, 1), trim = TRUE, nsmall = 0), "%)")
          }
        }
        dominant_labels <- .truncate_for_label(dom_labels, max_chars = 34)
        dominant_labels[base::is.na(dominant_labels)] <- ""
        ha_annos[[dominant_name]] <- ComplexHeatmap::anno_text(
          dominant_labels,
          which = "row",
          just = "left",
          location = 0,
          width = grid::unit(celltype_bar_dominant_width_cm, "cm"),
          gp = grid::gpar(fontsize = 7.2, col = "#333333")
        )
      }
      if (show_celltype_bars) {
        lgd_list[[base::length(lgd_list) + 1]] <- ComplexHeatmap::Legend(
          labels = .clean_celltype_label(entry$cell_types),
          title = if (n_enrichment_slots == 1) "Cell type" else anno_name,
          legend_gp = grid::gpar(col = entry$colors),
          type = "points",
          pch = 15
        )
      }
    }
  }
  ha <- base::do.call(
    ComplexHeatmap::HeatmapAnnotation,
    base::c(
      ha_annos,
      list(
        which = "row",
        width = grid::unit(row_annotation_total_width_cm, "cm"),
        annotation_name_side = "top",
        gap = grid::unit(2, "mm"),
        annotation_name_gp = grid::gpar(fontsize = 8)
      )
    )
  )
  build_row_annotation <- function(anno_scale = 1, build_legends = TRUE) {
    anno_scale <- .hc_first_numeric_value(anno_scale[[1]])
    if (!base::is.finite(anno_scale) || anno_scale <= 0) {
      anno_scale <- 1
    }
    anno_scale <- base::max(0.55, base::min(1, anno_scale))

    module_box_width_cm_use <- module_box_width_cm_draw * anno_scale
    gene_bar_width_cm_use <- if (show_gene_bar) 2.5 * anno_scale else 0
    gene_text_width_cm_use <- if (show_gene_text) 1.2 * anno_scale else 0
    module_sig_width_cm_scaled <- if (show_module_sig) module_sig_width_cm_use * anno_scale else 0
    enrichment_slot_bar_width_cm <- if (show_celltype_bars) 2.2 * anno_scale else 0
    enrichment_slot_text_width_cm <- if (show_celltype_text) celltype_bar_dominant_width_cm * anno_scale else 0
    enrichment_slot_gap_cm <- 0.15 * anno_scale
    n_enrichment_slots <- base::length(enrichment_entries)
    enrichment_total_width_cm <- if (n_enrichment_slots > 0) {
      (n_enrichment_slots * (enrichment_slot_bar_width_cm + enrichment_slot_text_width_cm)) +
        ((n_enrichment_slots - 1) * enrichment_slot_gap_cm)
    } else {
      0
    }
    base_row_width_cm <- module_box_width_cm_use +
      gene_text_width_cm_use +
      gene_bar_width_cm_use +
      module_sig_width_cm_scaled +
      (0.8 * anno_scale)
    row_annotation_total_width_cm <- base_row_width_cm + enrichment_total_width_cm
    module_label_fit_use <- .hc_module_label_fit_pt(
      module_label_pt_size = module_label_pt_size_draw,
      module_box_width_cm = module_box_width_cm_use,
      module_labels_display = module_labels_display,
      n_heat_rows = n_heat_rows,
      cell_size_mm = cell_size_mm * anno_scale,
      module_label_fontsize = module_label_fontsize,
      use_fontsize_request = user_set_module_label_fontsize && !user_set_module_label_pt_size,
      fontface = "bold"
    )
    module_label_pt_size_unit_use <- if (base::length(module_labels_display) > 0 &&
      base::length(module_label_fit_use$pt_size) == base::length(module_labels_display)) {
      grid::unit(module_label_fit_use$pt_size, "pt")
    } else {
      grid::unit(module_label_pt_size_draw, "snpc")
    }
    module_label_fontsize_draw_use <- .hc_module_label_effective_fontsize(
      module_label_fit = module_label_fit_use,
      fallback_fontsize = module_label_fontsize,
      fallback_pt_size = module_label_pt_size_draw
    )

    module_box_anno <- .hc_module_label_box_annotation(
      values = row_order,
      colors = cluster_colors,
      labels = module_labels_display,
      label_color = module_label_color,
      label_fontsize_pt = module_label_fit_use$base_pt_size,
      fontface = "bold",
      width_cm = module_box_width_cm_use,
      border_gp = module_box_border_gp,
      which = "row"
    )

    gene_bar_anno <- if (show_gene_bar) {
      ComplexHeatmap::anno_barplot(
        c_df$gene_no,
        width = grid::unit(gene_bar_width_cm_use, "cm"),
        which = "row"
      )
    } else {
      ComplexHeatmap::anno_empty(width = grid::unit(0, "mm"), which = "row", border = FALSE)
    }

    gene_text_anno <- if (show_gene_text) {
      if (gene_count_renderer == "pch") {
        ComplexHeatmap::anno_simple(
          x = base::rep("count_text", base::nrow(c_df)),
          col = c(count_text = "transparent"),
          pch = base::as.character(c_df$gene_no),
          pt_gp = grid::gpar(
            col = "black",
            fontsize = gene_count_fontsize,
            fontface = gene_count_fontface
          ),
          pt_size = grid::unit(gene_count_pt_size, "snpc"),
          gp = grid::gpar(col = NA),
          simple_anno_size = grid::unit(gene_text_width_cm_use, "cm"),
          which = "row"
        )
      } else {
        ComplexHeatmap::anno_text(
          base::as.character(c_df$gene_no),
          width = grid::unit(gene_text_width_cm_use, "cm"),
          just = "left",
          gp = grid::gpar(
            fontsize = gene_count_fontsize,
            fontface = gene_count_fontface
          ),
          which = "row"
        )
      }
    } else {
      ComplexHeatmap::anno_empty(width = grid::unit(0, "mm"), which = "row", border = FALSE)
    }

    module_sig_anno <- if (isTRUE(show_module_sig)) {
      sig_cols <- if (!base::is.null(module_sig_q) && base::length(module_sig_q) == base::length(module_sig_labels)) {
        base::ifelse(
          base::is.finite(module_sig_q) & module_sig_q <= module_significance_p_cutoffs[[3]],
          "black",
          "#8a8a8a"
        )
      } else {
        base::ifelse(base::nzchar(module_sig_labels), "black", "#8a8a8a")
      }
      if (isTRUE(module_significance_show_qvalue)) {
        ComplexHeatmap::anno_text(
          module_sig_labels,
          width = grid::unit(module_sig_width_cm_scaled, "cm"),
          just = "center",
          gp = grid::gpar(
            col = sig_cols,
            fontsize = base::max(9.5, gene_count_fontsize),
            fontface = "bold"
          ),
          which = "row"
        )
      } else {
        ComplexHeatmap::anno_simple(
          x = base::rep("sig_bg", base::length(module_sig_labels)),
          col = c(sig_bg = "#f5f5f5"),
          pch = module_sig_labels,
          pt_gp = grid::gpar(
            col = sig_cols,
            fontsize = base::max(11, gene_count_fontsize + 1),
            fontface = "bold"
          ),
          pt_size = grid::unit(0.9, "snpc"),
          simple_anno_size = grid::unit(module_sig_width_cm_scaled, "cm"),
          gp = grid::gpar(col = "#d0d0d0"),
          which = "row"
        )
      }
    } else {
      ComplexHeatmap::anno_empty(width = grid::unit(0, "mm"), which = "row", border = FALSE)
    }

    lgd_list_local <- if (isTRUE(build_legends)) list() else NULL
    ha_annos <- list(
      modules = module_box_anno,
      `# genes` = gene_text_anno,
      genes = gene_bar_anno
    )
    if (isTRUE(show_module_sig)) {
      sig_anno_name <- module_significance_annotation_name
      if (sig_anno_name %in% base::names(ha_annos)) {
        sig_anno_name <- base::paste0(sig_anno_name, "_sig")
      }
      ha_annos[[sig_anno_name]] <- module_sig_anno
    }
    if (n_enrichment_slots > 0) {
      for (i in base::seq_along(enrichment_entries)) {
        entry <- enrichment_entries[[i]]
        anno_name <- entry$label
        if (show_celltype_bars) {
          ha_annos[[anno_name]] <- ComplexHeatmap::anno_barplot(
            entry$mat,
            width = grid::unit(enrichment_slot_bar_width_cm, "cm"),
            gp = grid::gpar(fill = entry$colors, col = entry$colors),
            baseline = 0,
            which = "row"
          )
        }
        if (show_celltype_text) {
          dominant_name <- base::paste0(anno_name, "_dominant")
          dom_mat <- entry$mat
          dom_labels <- base::character(base::nrow(dom_mat))
          for (ri in base::seq_len(base::nrow(dom_mat))) {
            vv <- .hc_as_numeric_safely(dom_mat[ri, , drop = TRUE])
            if (!base::any(base::is.finite(vv) & vv > 0)) {
              dom_labels[[ri]] <- ""
              next
            }
            ctn <- .clean_celltype_label(base::colnames(dom_mat))
            is_other_like <- base::tolower(ctn) == base::tolower(celltype_bar_other_label)
            idx_pool <- base::which(!is_other_like & base::is.finite(vv) & vv > 0)
            if (base::length(idx_pool) == 0) {
              idx_pool <- base::which(base::is.finite(vv) & vv > 0)
            }
            if (base::length(idx_pool) == 0) {
              dom_labels[[ri]] <- ""
              next
            }
            j <- idx_pool[[base::which.max(vv[idx_pool])]]
            ct_lab <- ctn[[j]]
            ct_val <- vv[[j]]
            if (!base::nzchar(ct_lab) || !base::is.finite(ct_val)) {
              dom_labels[[ri]] <- ""
            } else {
              dom_labels[[ri]] <- base::paste0(ct_lab, " (", base::format(round(ct_val, 1), trim = TRUE, nsmall = 0), "%)")
            }
          }
          dominant_labels <- .truncate_for_label(dom_labels, max_chars = 34)
          dominant_labels[base::is.na(dominant_labels)] <- ""
          ha_annos[[dominant_name]] <- ComplexHeatmap::anno_text(
            dominant_labels,
            which = "row",
            just = "left",
            location = 0,
            width = grid::unit(enrichment_slot_text_width_cm, "cm"),
            gp = grid::gpar(fontsize = 7.2, col = "#333333")
          )
        }
        if (isTRUE(show_celltype_bars) && isTRUE(build_legends)) {
          lgd_list_local[[base::length(lgd_list_local) + 1]] <- ComplexHeatmap::Legend(
            labels = .clean_celltype_label(entry$cell_types),
            title = if (n_enrichment_slots == 1) "Cell type" else anno_name,
            legend_gp = grid::gpar(col = entry$colors),
            type = "points",
            pch = 15
          )
        }
      }
    }

    ha <- base::do.call(
      ComplexHeatmap::HeatmapAnnotation,
      base::c(
        ha_annos,
        list(
          which = "row",
          width = grid::unit(row_annotation_total_width_cm, "cm"),
          annotation_name_side = "top",
          gap = grid::unit(base::max(1, 2 * anno_scale), "mm"),
          annotation_name_gp = grid::gpar(fontsize = 8)
        )
      )
    )

    list(
      annotation = ha,
      row_annotation_total_width_cm = row_annotation_total_width_cm,
      legend_list = lgd_list_local
    )
  }
  # --- 6. Setup Column Annotations ---

  anno_list <- NULL

  if (!base::length(column_anno_categorical) == 0) {
    for (a in base::seq_along(column_anno_categorical)) {
      tmp_colour <- grDevices::colorRampPalette(c("#332288", "#117733", "#44aa99", "#88ccee", "#cc6677", "#aa4499", "#882255"))(base::ncol(column_anno_categorical[[a]]))

      if (cat_as_bp[a] == TRUE) {
        column_anno_categorical[[a]][base::is.na(column_anno_categorical[[a]])] <- 0
        current_anno <- ComplexHeatmap::HeatmapAnnotation(
          col_anno = ComplexHeatmap::anno_barplot(column_anno_categorical[[a]] %>% base::as.matrix(),
            width = grid::unit(2, "cm"),
            gp = grid::gpar(fill = tmp_colour, col = tmp_colour)
          ),
          which = "column",
          height = grid::unit(1, "cm"),
          annotation_name_side = "right",
          gap = grid::unit(2, "mm"),
          annotation_name_rot = 0,
          annotation_name_gp = grid::gpar(fontsize = 8),
          annotation_label = base::names(column_anno_categorical)[a]
        )
      } else {
        current_anno <- ComplexHeatmap::HeatmapAnnotation(
          col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),
            width = grid::unit(2, "cm"),
            gp = grid::gpar(col = tmp_colour),
            add_points = TRUE,
            pt_gp = grid::gpar(col = tmp_colour), pch = 16
          ),
          which = "column",
          height = grid::unit(1, "cm"),
          annotation_name_side = "right",
          gap = grid::unit(2, "mm"),
          annotation_name_rot = 0,
          annotation_name_gp = grid::gpar(fontsize = 8),
          annotation_label = base::names(column_anno_categorical)[a]
        )
      }

      if (base::is.null(anno_list)) {
        anno_list <- current_anno
      } else {
        anno_list <- ComplexHeatmap::add_heatmap(anno_list, current_anno, direction = c("vertical"))
      }

      lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(
        labels = base::colnames(column_anno_categorical[[a]] %>% base::as.matrix()),
        title = base::names(column_anno_categorical)[a],
        legend_gp = grid::gpar(col = tmp_colour),
        type = "points", pch = 15
      ))
    }
  }

  if (!base::length(column_anno_numerical) == 0) {
    for (a in base::seq_along(column_anno_numerical)) {
      tmp_col_anno_2 <- column_anno_numerical[[a]]
      tmp_col_anno_2 <- tmp_col_anno_2[base::colnames(mat_heatmap)]

      current_anno <- ComplexHeatmap::HeatmapAnnotation(
        cont_anno = ComplexHeatmap::anno_boxplot(tmp_col_anno_2, height = grid::unit(1, "cm")),
        which = "column",
        annotation_name_side = "right",
        gap = grid::unit(2, "mm"),
        annotation_name_rot = 0,
        annotation_name_gp = grid::gpar(fontsize = 8),
        annotation_label = base::names(column_anno_numerical)[a], show_legend = FALSE
      )

      if (base::is.null(anno_list)) {
        anno_list <- current_anno
      } else {
        anno_list <- ComplexHeatmap::add_heatmap(anno_list, current_anno, direction = c("vertical"))
      }
    }
  }

  column_labels_display <- .hc_gfc_display_col_labels(hcobject, base::colnames(mat_heatmap))
  all_conditions <- .hc_gfc_display_count_labels(hcobject, base::colnames(mat_heatmap))
  column_gap_k_enabled <- !base::is.numeric(k) || base::length(k) == 0 || base::all(k <= 0)
  column_gap_enabled <- (isTRUE(smart_column_gaps) || !base::is.null(column_gap_by)) &&
    isTRUE(column_gap_k_enabled)
  column_gap_spec <- .hc_heatmap_column_gap_spec(
    hcobject = hcobject,
    cols = base::colnames(mat_heatmap),
    cluster_columns = cluster_columns,
    gap_mm = column_gap_mm * overall_plot_scale,
    enabled = column_gap_enabled,
    metadata_column = if (isTRUE(column_gap_enabled)) column_gap_by else NULL
  )

  if (base::is.null(anno_list)) {
    anno_list <- ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions))
  } else {
    anno_list <- ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions)), direction = c("vertical"))
  }


  # --- 7. Plotting and Output ---

  # Keep module-expression tiles square-like for publication consistency.
  cell_size_mm <- .hc_cluster_heatmap_cell_size_mm(
    n_heat_rows = n_heat_rows,
    n_heat_cols = n_heat_cols,
    duplicate_condition_width_scale = duplicate_condition_width_scale,
    module_box_width_cm_draw = module_box_width_cm_draw,
    overall_plot_scale = overall_plot_scale
  )
  hm_width <- grid::unit((n_heat_cols * cell_size_mm) + column_gap_spec$total_gap_mm, "mm")
  hm_height <- grid::unit(n_heat_rows * cell_size_mm, "mm")
  max_col_chars <- if (base::length(column_labels_display) > 0) {
    base::max(base::nchar(column_labels_display), na.rm = TRUE)
  } else {
    10
  }
  shared_column_name_max_cm <- base::max(
    2.8,
    base::min(
      8.0,
      (max_col_chars * (10 * overall_plot_scale) * 0.022) + 0.6
    )
  )
  row_dend_width_mm <- if (isTRUE(cluster_rows)) {
    base::max(8, base::min(14, n_heat_rows * 0.75))
  } else {
    6
  }
  row_dend_width_mm <- row_dend_width_mm * overall_plot_scale
  column_dend_height_mm <- if (isTRUE(cluster_columns)) {
    base::max(10, base::min(18, n_heat_cols * 2.2))
  } else {
    6
  }
  column_dend_height_mm <- column_dend_height_mm * overall_plot_scale
  gfc_palette <- grDevices::colorRampPalette(gfc_colors)(51)
  gfc_col_fun <- circlize::colorRamp2(
    seq(gfc_scale_limits[1], gfc_scale_limits[2], length.out = base::length(gfc_palette)),
    gfc_palette
  )

  # ComplexHeatmap/grid objects can keep draw-time viewport state.
  # Clone objects before repeated draws (PDF + current device) and retry safely.
  deep_clone <- function(x) {
    .hc_safe_deep_clone(x, context = "cluster heatmap object")
  }

  safe_draw <- function(ht_obj,
                        legend_obj = NULL,
                        padding_obj = NULL,
                        heatmap_legend_side = "right",
                        annotation_legend_side = "right",
                        context = "cluster heatmap",
                        capture_current_device = FALSE) {
    draw_once <- function(obj,
                          lgd,
                          show_ann_legend = TRUE,
                          show_heat_legend = TRUE,
                          use_padding = TRUE) {
      args <- list(
        object = obj,
        newpage = TRUE,
        merge_legends = if (isTRUE(show_ann_legend)) TRUE else FALSE,
        show_heatmap_legend = show_heat_legend,
        show_annotation_legend = show_ann_legend,
        heatmap_legend_side = heatmap_legend_side,
        annotation_legend_side = annotation_legend_side
      )
      if (!base::is.null(lgd) && base::length(lgd) > 0) {
        args$annotation_legend_list <- lgd
      }
      if (isTRUE(use_padding) && !base::is.null(padding_obj)) {
        args$padding <- padding_obj
      }
      base::do.call(ComplexHeatmap::draw, args)
    }
    run_attempt <- function(expr_fun) {
      if (!isTRUE(capture_current_device)) {
        return(list(result = expr_fun(), grob = NULL))
      }
      result <- NULL
      grob <- grid::grid.grabExpr(
        {
          result <- expr_fun()
        },
        wrap = TRUE
      )
      list(result = result, grob = grob)
    }
    finish_attempt <- function(attempt) {
      if (isTRUE(capture_current_device) && !base::is.null(attempt$grob)) {
        grid::grid.newpage()
        grid::grid.draw(attempt$grob)
      }
      attempt$result
    }

    first <- try(
      run_attempt(function() draw_once(ht_obj, legend_obj, show_ann_legend = TRUE)),
      silent = TRUE
    )
    if (!inherits(first, "try-error")) {
      return(finish_attempt(first))
    }

    if (!isTRUE(capture_current_device)) {
      try(grid::grid.newpage(), silent = TRUE)
    }
    second <- try(
      run_attempt(
        function() {
          draw_once(
            deep_clone(ht_obj),
            deep_clone(legend_obj),
            show_ann_legend = TRUE
          )
        }
      ),
      silent = TRUE
    )
    if (!inherits(second, "try-error")) {
      warning(
        "Recovered from a transient grid viewport issue while drawing ",
        context,
        ".",
        call. = FALSE
      )
      return(finish_attempt(second))
    }

    if (!isTRUE(capture_current_device)) {
      try(grid::grid.newpage(), silent = TRUE)
    }
    third <- try(
      run_attempt(
        function() {
          draw_once(
            deep_clone(ht_obj),
            NULL,
            show_ann_legend = FALSE
          )
        }
      ),
      silent = TRUE
    )
    if (!inherits(third, "try-error")) {
      warning(
        "Draw fallback for ",
        context,
        ": annotation legends were disabled after viewport issues.",
        call. = FALSE
      )
      return(finish_attempt(third))
    }

    if (!isTRUE(capture_current_device)) {
      try(grid::grid.newpage(), silent = TRUE)
    }
    fourth <- try(
      run_attempt(
        function() {
          draw_once(
            deep_clone(ht_obj),
            NULL,
            show_ann_legend = FALSE,
            show_heat_legend = TRUE,
            use_padding = FALSE
          )
        }
      ),
      silent = TRUE
    )
    if (!inherits(fourth, "try-error")) {
      warning(
        "Draw fallback for ",
        context,
        ": legends/padding were simplified after viewport issues.",
        call. = FALSE
      )
      return(finish_attempt(fourth))
    }

    if (!isTRUE(capture_current_device)) {
      try(grid::grid.newpage(), silent = TRUE)
    }
    fifth <- try(
      run_attempt(
        function() {
          draw_once(
            deep_clone(ht_obj),
            NULL,
            show_ann_legend = FALSE,
            show_heat_legend = FALSE,
            use_padding = FALSE
          )
        }
      ),
      silent = TRUE
    )
    if (!inherits(fifth, "try-error")) {
      warning(
        "Draw fallback for ",
        context,
        ": all legends were disabled after viewport issues.",
        call. = FALSE
      )
      return(finish_attempt(fifth))
    }

    stop(fifth)
  }

  build_hm <- function(row_dend_mm,
                       col_dend_mm,
                       column_name_max_cm,
                       fixed_size = TRUE,
                       right_annotation_obj = ha,
                       body_width_mm = NULL,
                       body_height_mm = NULL) {
    legend_height_mm <- base::max(24, 4.5 * base::max(1, base::length(gfc_scale_breaks))) * overall_plot_scale
    heat_legend_param <- list(
      title = "GFC",
      at = gfc_scale_breaks,
      labels = gfc_scale_labels,
      title_gp = grid::gpar(fontsize = 7.6 * overall_plot_scale, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 6.6 * overall_plot_scale),
      grid_width = grid::unit(3.2 * overall_plot_scale, "mm"),
      legend_height = grid::unit(legend_height_mm, "mm")
    )
    if (identical(heatmap_legend_side_mode, "bottom")) {
      heat_legend_param$direction <- "horizontal"
      heat_legend_param$legend_width <- grid::unit(32 * overall_plot_scale, "mm")
      heat_legend_param$title_position <- "leftcenter"
    }

    use_column_gap <- !base::is.null(column_gap_spec$column_split) &&
      !base::is.null(column_gap_spec$column_gap) &&
      (!base::is.numeric(k) || base::length(k) == 0 || base::all(k <= 0))

    hm_args <- list(
      mat_heatmap,
      name = "GFC",
      right_annotation = deep_clone(right_annotation_obj),
      col = gfc_col_fun,
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      cluster_columns = cluster_columns,
      cluster_rows = cluster_rows_for_heatmap,
      row_dend_reorder = FALSE,
      column_names_rot = 90,
      column_labels = column_labels_display,
      column_names_centered = FALSE,
      row_dend_width = grid::unit(row_dend_mm, "mm"),
      column_dend_height = grid::unit(col_dend_mm, "mm"),
      row_dend_gp = dendrogram_line_gp,
      column_dend_gp = dendrogram_line_gp,
      column_names_max_height = grid::unit(column_name_max_cm, "cm"),
      show_row_names = show_module_color_names,
      row_names_gp = grid::gpar(fontsize = 8 * overall_plot_scale),
      column_names_gp = grid::gpar(fontsize = 10 * overall_plot_scale),
      rect_gp = heatmap_cell_border_gp,
      show_heatmap_legend = TRUE,
      heatmap_legend_param = heat_legend_param,
      column_km = k
    )
    if (isTRUE(use_column_gap)) {
      hm_args$column_split <- column_gap_spec$column_split
      hm_args$column_gap <- column_gap_spec$column_gap
      hm_args$cluster_column_slices <- FALSE
      hm_args$column_title <- column_gap_spec$slice_titles
    }
    if (isTRUE(fixed_size)) {
      if (!base::is.null(body_width_mm) && !base::is.null(body_height_mm)) {
        hm_args$width <- grid::unit(body_width_mm, "mm")
        hm_args$height <- grid::unit(body_height_mm, "mm")
      } else {
        hm_args$width <- hm_width
        hm_args$height <- hm_height
      }
    }
    do.call(ComplexHeatmap::Heatmap, hm_args)
  }

  draw_state <- new.env(parent = emptyenv())
  draw_state$hm_pdf_drawn <- NULL
  heatmap_export_files <- NULL
  heatmap_legend_side_mode <- gfc_legend_side
  annotation_legend_side_mode <- "right"
  right_pad_export_mm <- if (show_celltype_text && show_celltype_bars && base::length(lgd_list) > 0) {
    52
  } else if (show_celltype_text && !show_celltype_bars) {
    26
  } else {
    36
  }
  bottom_pad_export_mm <- if (heatmap_legend_side_mode == "bottom") {
    52
  } else {
    base::max(18, base::min(34, shared_column_name_max_cm * 3.0))
  }
  draw_padding_export_mm <- c(26, 18, bottom_pad_export_mm, right_pad_export_mm) * overall_plot_scale
  draw_padding_export <- grid::unit(draw_padding_export_mm, "mm")
  export_total_width_mm <- base::as.numeric(hm_width) +
    row_dend_width_mm +
    (row_annotation_total_width_cm * 10) +
    (18 * overall_plot_scale) +
    right_pad_export_mm
  export_total_height_mm <- base::as.numeric(hm_height) +
    column_dend_height_mm +
    (shared_column_name_max_cm * 10) +
    (26 * overall_plot_scale) +
    bottom_pad_export_mm
  png_width_in <- if (pdf_width_is_default && pdf_height_is_default) {
    base::max(4.8, export_total_width_mm / 25.4)
  } else {
    pdf_width
  }
  png_height_in <- if (pdf_width_is_default && pdf_height_is_default) {
    base::max(6.0, export_total_height_mm / 25.4)
  } else {
    pdf_height
  }
  if (isTRUE(write_pdf)) {
    export_file <- paste0(
      hcobject[["working_directory"]][["dir_output"]],
      hcobject[["global_settings"]][["save_folder"]],
      "/",
      file_name
    )
    draw_export_heatmap <- function() {
      hm_export <- build_hm(
        row_dend_mm = row_dend_width_mm,
        col_dend_mm = column_dend_height_mm,
        column_name_max_cm = shared_column_name_max_cm,
        fixed_size = TRUE,
        right_annotation_obj = ha,
        body_width_mm = base::as.numeric(hm_width),
        body_height_mm = base::as.numeric(hm_height)
      )

      anno_list_export_src <- deep_clone(anno_list)
      lgd_list_export <- deep_clone(lgd_list)
      anno_list_export <- if (base::is.null(anno_list_export_src)) {
        hm_export
      } else {
        ComplexHeatmap::add_heatmap(hm_export, anno_list_export_src, direction = c("vertical"))
      }

      drawn_export <- tryCatch(
        safe_draw(
          ht_obj = anno_list_export,
          legend_obj = lgd_list_export,
          padding_obj = draw_padding_export,
          heatmap_legend_side = heatmap_legend_side_mode,
          annotation_legend_side = annotation_legend_side_mode,
          context = "cluster heatmap export"
        ),
        error = function(e) {
          warning(
            "Could not fully draw cluster heatmap export: ",
            base::conditionMessage(e),
            call. = FALSE
          )
          NULL
        }
      )
      if (base::is.null(draw_state$hm_pdf_drawn)) {
        draw_state$hm_pdf_drawn <- drawn_export
      }
      invisible(drawn_export)
    }

    heatmap_export_files <- tryCatch(
      .hc_export_single_page_plot(
        file = export_file,
        width = pdf_width,
        height = pdf_height,
        png_width = png_width_in,
        png_height = png_height_in,
        pointsize = pdf_pointsize,
        res = pdf_dpi,
        pdf_dpi = pdf_dpi,
        draw_fun = draw_export_heatmap
      ),
      error = function(e) {
        warning(
          "Could not export cluster heatmap PDF/PNG: ",
          base::conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
  }

  if (isTRUE(module_labels_have_split_suffix) || duplicate_condition_width_scale > 1) {
    screen_fit_scale <- .hc_heatmap_screen_fit_scale(
      total_width_mm = export_total_width_mm,
      total_height_mm = export_total_height_mm
    )
    row_dend_screen_mm <- base::max(6 * overall_plot_scale, row_dend_width_mm * screen_fit_scale)
    col_dend_screen_mm <- base::max(6 * overall_plot_scale, column_dend_height_mm * screen_fit_scale)
    col_label_max_cm_screen <- base::max(2.4, shared_column_name_max_cm * screen_fit_scale)
    draw_padding_screen_mm <- base::pmax(
      c(10, 10, 14, 14) * overall_plot_scale,
      draw_padding_export_mm * screen_fit_scale
    )
    draw_padding_screen <- grid::unit(draw_padding_screen_mm, "mm")
    hm_screen_w_mm <- base::max(20 * overall_plot_scale, base::as.numeric(hm_width) * screen_fit_scale)
    hm_screen_h_mm <- base::max(24 * overall_plot_scale, base::as.numeric(hm_height) * screen_fit_scale)
    row_annotation_screen_scale <- base::max(0.58, screen_fit_scale)
  } else {
    row_dend_screen_mm <- row_dend_width_mm
    col_dend_screen_mm <- column_dend_height_mm
    col_label_max_cm_screen <- shared_column_name_max_cm
    draw_padding_screen <- draw_padding_export
    hm_screen_w_mm <- base::as.numeric(hm_width)
    hm_screen_h_mm <- base::as.numeric(hm_height)
    row_annotation_screen_scale <- 1
  }
  ha_screen <- if (row_annotation_screen_scale < 0.999) {
    build_row_annotation(anno_scale = row_annotation_screen_scale, build_legends = FALSE)$annotation
  } else {
    ha
  }

  hm_screen <- build_hm(
    row_dend_mm = row_dend_screen_mm,
    col_dend_mm = col_dend_screen_mm,
    column_name_max_cm = col_label_max_cm_screen,
    fixed_size = TRUE,
    right_annotation_obj = ha_screen,
    body_width_mm = hm_screen_w_mm,
    body_height_mm = hm_screen_h_mm
  )

  anno_list_screen_src <- deep_clone(anno_list)
  lgd_list_screen <- deep_clone(lgd_list)
  anno_list_screen <- if (base::is.null(anno_list_screen_src)) {
    hm_screen
  } else {
    ComplexHeatmap::add_heatmap(hm_screen, anno_list_screen_src, direction = c("vertical"))
  }

  hm_w_lgd <- tryCatch(
    safe_draw(
      ht_obj = anno_list_screen,
      legend_obj = lgd_list_screen,
      padding_obj = draw_padding_screen,
      heatmap_legend_side = heatmap_legend_side_mode,
      annotation_legend_side = annotation_legend_side_mode,
      context = "cluster heatmap current device",
      capture_current_device = TRUE
    ),
    error = function(e) {
      warning(
        "Could not draw cluster heatmap on current device: ",
        base::conditionMessage(e),
        call. = FALSE
      )
      fallback <- try(
        ComplexHeatmap::draw(
          deep_clone(anno_list_screen),
          newpage = TRUE,
          merge_legends = FALSE,
          show_heatmap_legend = TRUE,
          show_annotation_legend = FALSE,
          heatmap_legend_side = heatmap_legend_side_mode,
          annotation_legend_side = annotation_legend_side_mode
        ),
        silent = TRUE
      )
      if (!inherits(fallback, "try-error")) {
        return(fallback)
      }
      draw_state$hm_pdf_drawn
    }
  )
  module_label_map <- NULL
  if (!base::is.null(module_labels)) {
    module_label_map <- stats::setNames(base::as.character(module_labels), base::as.character(row_order))
  }
  module_export_labels <- if (!base::is.null(module_labels)) {
    base::as.character(module_labels)
  } else {
    base::as.character(row_order)
  }
  if (base::length(module_export_labels) != base::length(row_order)) {
    module_export_labels <- base::as.character(row_order)
  }
  module_export_map <- stats::setNames(module_export_labels, base::as.character(row_order))
  module_gene_rows <- base::lapply(base::as.character(row_order), function(cl) {
    genes <- cluster_gene_map[[cl]]
    if (base::length(genes) == 0) {
      return(NULL)
    }
    base::data.frame(
      genes = genes,
      module = base::as.character(module_export_map[[cl]]),
      stringsAsFactors = FALSE
    )
  })
  module_gene_rows <- module_gene_rows[!base::vapply(module_gene_rows, base::is.null, FUN.VALUE = base::logical(1))]
  module_gene_list_tbl <- if (base::length(module_gene_rows) > 0) {
    out <- base::do.call(base::rbind, module_gene_rows)
    out <- out[, base::c("genes", "module"), drop = FALSE]
    base::rownames(out) <- NULL
    out
  } else {
    base::data.frame(
      genes = base::character(0),
      module = base::character(0),
      stringsAsFactors = FALSE
    )
  }
  if (base::isTRUE(write_module_tables)) {
    stored_module_label_map <- hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]]
    split_history <- hcobject[["satellite_outputs"]][["module_split_history"]]
    module_gene_list_name <- .hc_module_gene_list_filename(
      module_label_map = base::c(stored_module_label_map, module_export_labels),
      split_history = split_history
    )
    module_gene_list_file <- base::paste0(
      hcobject[["working_directory"]][["dir_output"]],
      hcobject[["global_settings"]][["save_folder"]],
      "/",
      module_gene_list_name
    )
    tryCatch(
      {
        .hc_write_xlsx_atomic(
          x = list(module_gene_list = module_gene_list_tbl),
          file = module_gene_list_file,
          overwrite = TRUE
        )
      },
      error = function(e) {
        warning(
          "Could not write ", module_gene_list_name, ": ",
          base::conditionMessage(e)
        )
      }
    )
  }
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_map"), module_label_map)
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "module_gene_list"), module_gene_list_tbl)

  # Mean-GFC-per-module-and-group table: the numeric values behind the heatmap
  # cells (rows = modules, columns = sample groups). `mat_heatmap` is keyed by
  # cluster colour, so relabel rows with the displayed module labels and keep the
  # colour as a reference column.
  module_gfc_means_tbl <- tryCatch(
    {
      ordered_colors <- base::as.character(row_order)
      ordered_colors <- ordered_colors[ordered_colors %in% base::rownames(mat_heatmap)]
      ordered_colors <- base::c(
        ordered_colors,
        base::setdiff(base::rownames(mat_heatmap), ordered_colors)
      )
      gfc_means_mat <- mat_heatmap[ordered_colors, , drop = FALSE]
      row_labels <- base::as.character(module_export_map[ordered_colors])
      missing_lab <- base::is.na(row_labels) | !base::nzchar(row_labels)
      row_labels[missing_lab] <- ordered_colors[missing_lab]
      base::data.frame(
        module = row_labels,
        cluster_color = ordered_colors,
        base::as.data.frame(gfc_means_mat, check.names = FALSE, stringsAsFactors = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    },
    error = function(e) {
      base::warning(
        "Could not build module GFC means table: ",
        base::conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )
  if (!base::is.null(module_gfc_means_tbl)) {
    if (base::isTRUE(write_module_tables)) {
      module_gfc_means_file <- base::paste0(
        hcobject[["working_directory"]][["dir_output"]],
        hcobject[["global_settings"]][["save_folder"]],
        "/Module_GFC_Means.xlsx"
      )
      tryCatch(
        .hc_write_xlsx_atomic(
          x = base::list(module_gfc_means = module_gfc_means_tbl),
          file = module_gfc_means_file,
          overwrite = TRUE
        ),
        error = function(e) {
          base::warning(
            "Could not write Module_GFC_Means.xlsx: ",
            base::conditionMessage(e)
          )
        }
      )
    }
    .hc_set_bridge_hcobject_slot(c("satellite_outputs", "module_gfc_means"), module_gfc_means_tbl)
  }

  module_label_fit_pt <- .hc_as_numeric_safely(module_label_fit$pt_size)
  module_label_fit_base_pt <- .hc_as_numeric_safely(module_label_fit$base_pt_size)
  module_label_fit_pt <- module_label_fit_pt[base::is.finite(module_label_fit_pt)]
  module_label_fit_base_pt <- module_label_fit_base_pt[base::is.finite(module_label_fit_base_pt)]
  module_label_fit_min <- if (base::length(module_label_fit_pt) > 0) base::min(module_label_fit_pt) else NA_real_
  module_label_fit_max <- if (base::length(module_label_fit_pt) > 0) base::max(module_label_fit_pt) else NA_real_
  module_label_fit_requested_max <- if (base::length(module_label_fit_base_pt) > 0) {
    base::max(module_label_fit_base_pt)
  } else {
    NA_real_
  }
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_mode"), module_label_mode)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_numbering"), module_label_numbering)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_fontsize"), module_label_fontsize)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_pt_size"), module_label_pt_size)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_pt_size_effective_pt_min"), module_label_fit_min)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_pt_size_effective_pt_max"), module_label_fit_max)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_pt_size_requested_pt_max"), module_label_fit_requested_max)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_auto_fit_shrunk"), isTRUE(module_label_fit$shrunk))
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_auto_fit_width_limited"), isTRUE(module_label_fit$width_limited))
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_label_auto_fit_height_limited"), isTRUE(module_label_fit$height_limited))
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_box_width_cm"), module_box_width_cm)
  module_box_to_cell_ratio <- (.hc_first_numeric_value(module_box_width_cm) * 10) /
    .hc_first_numeric_value(cell_size_mm)
  if (!base::is.finite(module_box_to_cell_ratio) || module_box_to_cell_ratio <= 0) {
    module_box_to_cell_ratio <- NULL
  }
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_box_to_cell_ratio"), module_box_to_cell_ratio)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_cell_size_mm"), cell_size_mm)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "duplicate_condition_width_scale"), duplicate_condition_width_scale)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "smart_column_gaps"), smart_column_gaps)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "column_gap_by"), column_gap_by)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "column_gap_mm"), column_gap_mm)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_column_gap_mm"), column_gap_spec$total_gap_mm)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_column_gap_source"), column_gap_spec$source)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "gene_count_fontsize"), gene_count_fontsize)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "gene_count_renderer"), gene_count_renderer)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "gene_count_pt_size"), gene_count_pt_size)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "gfc_colors"), gfc_colors)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "gfc_scale_limits"), gfc_scale_limits)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "overall_plot_scale"), overall_plot_scale)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_output_files"), heatmap_export_files)
  final_row_order <- .hc_normalize_heatmap_axis_order(
    tryCatch(ComplexHeatmap::row_order(hm_w_lgd), error = function(e) NULL),
    base::rownames(mat_heatmap)
  )
  if (base::is.null(final_row_order) || base::length(final_row_order) == 0) {
    final_row_order <- .hc_normalize_heatmap_axis_order(displayed_row_order, base::rownames(mat_heatmap))
  }
  final_col_order <- .hc_normalize_heatmap_axis_order(
    tryCatch(ComplexHeatmap::column_order(hm_w_lgd), error = function(e) NULL),
    base::colnames(mat_heatmap)
  )
  if (base::is.null(final_col_order) || base::length(final_col_order) == 0) {
    final_col_order <- .hc_normalize_heatmap_axis_order(base::colnames(mat_heatmap), base::colnames(mat_heatmap))
  }
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_matrix"), mat_heatmap)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_row_order"), final_row_order)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_column_order"), final_col_order)
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_column_labels_display"), column_labels_display)
  if (module_label_mode == "prefix") {
    .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_prefix"), module_prefix)
  } else {
    .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "module_prefix"), NULL)
  }

  # Always cache the reusable raw heatmap object so downstream views such as
  # the LLM module summaries can mirror the main heatmap styling even when the
  # user does not request a heatmap object return value.
  .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_cluster_raw"), deep_clone(anno_list_screen))
  if (return_HM) {
    .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_cluster"), hm_w_lgd)
  } else {
    .hc_set_bridge_hcobject_slot(c("integrated_output", "cluster_calc", "heatmap_cluster"), NULL)
  }
}
