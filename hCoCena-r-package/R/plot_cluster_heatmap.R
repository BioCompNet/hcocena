
#' Plot Cluster Heatmap
#' 
#' Plots a heatmap with sample groups as columns and gene clusters as rows. The cells are coloured according to the mean GFC of a given cluster in the respective sample group.
#'  If categorical or numerical metadata annotations or user-defined enrichments have been created with satellite functions, they will be incorporated into the heatmap as column/row annotations.
#'  A `Module_Gene_List.xlsx` export (columns: `genes`, `module`) is written to the save folder
#'  based on the currently displayed module labels (including split-module labels, if present).
#'  For the available options check out the satellite functions.
#' @param col_order Defines the order in which the sample groups (conditions) appear in the heatmap. 
#'  Accepts a vector of strings giving the conditions in their desired order, default is NULL (automatic order, determined by clustering).
#'  If the parameter cluster_columns is not set to FALSE, this order will be overwritten.
#' @param row_order Like col_order but with cluster names.
#' @param cluster_columns A Boolean, whether or not to cluster the columns of the heatmap. Default is TRUE, overwrites col_order.
#' @param cluster_rows Like cluster_columns but for rows.
#' @param k The resulting cluster_columns tree is cut into k groups. Default is 0 (no cutting).
#' @param return_HM A Boolean whether of not to return the ComplexHeatmap object to hcobject$integrated_output$cluster_calc$heatmap_cluster in addition to plotting it. Default is TRUE.
#' @param cat_as_bp  A vector of Booleans which length is equivalent to the number of categorical meta data annotations created with the satellite functions. 
#'  Each Boolean states whether or not the categorical variable should be annotated as a bar plot (TRUE) or as a line plot (FALSE). 
#'  If you did not perform any meta data annotation, ignore this parameter. 
#' @param file_name A string giving the name of the file (with .pdf ending) to which the heatmap should be written. Default is "Heatmap_modules.pdf".
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
#'  Presets adjust automatic sizing only; explicitly set sizing parameters still take precedence.
#'  Default is "balanced".
#' @param module_label_color Text color for labels drawn in module boxes.
#' @param module_label_fontsize Optional numeric fontsize for module box labels. If NULL, a data-driven size is chosen.
#' @param module_label_pt_size Optional numeric scale for label glyph size inside module boxes.
#'  Uses `snpc` units in `anno_simple()`. If NULL, a data-driven size is chosen.
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
#'  If NULL, uses the legacy default `rev(RColorBrewer::brewer.pal(11, "RdBu"))`.
#' @param gfc_scale_limits Optional numeric vector controlling the module-heatmap
#'  color scale limits. Provide either one positive number (`x` -> `c(-x, x)`) or
#'  two numbers (`c(min, max)`). If NULL, uses stored limits from the previous
#'  main heatmap run, otherwise falls back to `c(-range_GFC, range_GFC)`.
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
#' @param include_dynamic_enrichment_slots Logical. If `TRUE`, also include
#'  dynamically named enrichment slots (e.g. `enriched_per_cluster_<db>`).
#'  Default is `FALSE` to keep the standard heatmap behavior unchanged.
#' @export



plot_cluster_heatmap <- function(col_order = NULL,
                                 row_order = NULL,
                                 cluster_columns = T,
                                 cluster_rows = T,
                                 k = 0,
                                 return_HM = T,
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
                                 pdf_width = 50,
                                 pdf_height = 30,
                                 pdf_pointsize = 11,
                                 pdf_dpi = 300,
                                 overall_plot_scale = 1,
                                 include_dynamic_enrichment_slots = FALSE){

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
    pdf_width = pdf_width,
    pdf_height = pdf_height,
    pdf_pointsize = pdf_pointsize,
    pdf_dpi = pdf_dpi,
    overall_plot_scale = overall_plot_scale,
    include_dynamic_enrichment_slots = include_dynamic_enrichment_slots
  )
}


plot_cluster_heatmap_new <- function(col_order = NULL, 
                                 row_order = NULL, 
                                 cluster_columns = T,
                                 cluster_rows = T, 
                                 k = 0, 
                                 return_HM = T, 
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
                                 pdf_width = 50,
                                 pdf_height = 30,
                                 pdf_pointsize = 11,
                                 pdf_dpi = 300,
                                 overall_plot_scale = 1,
                                 include_dynamic_enrichment_slots = FALSE){

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
  overall_plot_scale <- base::max(0.5, base::min(3, overall_plot_scale))
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

  gfc_scale_limits <- normalize_scale_limits(gfc_scale_limits)
  if (base::is.null(gfc_scale_limits)) {
    stored_limits <- tryCatch(
      normalize_scale_limits(hcobject[["integrated_output"]][["cluster_calc"]][["gfc_scale_limits"]]),
      error = function(e) NULL
    )
    if (!base::is.null(stored_limits)) {
      gfc_scale_limits <- stored_limits
    } else {
      fallback_lim <- suppressWarnings(base::as.numeric(hcobject[["global_settings"]][["range_GFC"]]))
      if (!base::is.finite(fallback_lim) || fallback_lim <= 0) {
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
    slot_df$count <- suppressWarnings(base::as.numeric(slot_df$count))
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
  if("column_annos_categorical" %in% base::names(hcobject[["satellite_outputs"]])){
    column_anno_categorical <- hcobject[["satellite_outputs"]][["column_annos_categorical"]]
  } else {
    column_anno_categorical <- NULL
  }
  
  if("column_annos_numerical" %in% base::names(hcobject[["satellite_outputs"]])){
    column_anno_numerical <- hcobject[["satellite_outputs"]][["column_annos_numerical"]]
  } else {
    column_anno_numerical <- NULL
  }
  
  if(base::is.null(cat_as_bp)){
    if(!base::is.null(column_anno_categorical)){
      cat_as_bp <- base::rep(F, base::length(column_anno_categorical))
    }
  }
  
  base::gc()
  
  # Filter for included clusters (non-white)
  c_df <- dplyr::filter(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], cluster_included == "yes")
  mat_heatmap <- NULL
  
  # --- 2. Build Heatmap Matrix (GFC Values) ---
  
  if(!base::is.null(row_order)){
    for (c in row_order){
      # Get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)
      
      # GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(hcobject[["integrated_output"]][["GFC_all_layers"]], Gene %in% genes)
      c_GFC_means <- base::apply(c_GFCs[, base::c(1:(base::ncol(c_GFCs)-1))], 2, base::mean)
      
      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)
    }
    base::rownames(mat_heatmap) <- row_order
    
  } else {
    for (c in base::unique(c_df$color)){
      # Get genes from the original cluster
      genes <- c_df[c_df$color == c, ] %>%
        dplyr::pull(., "gene_n") %>%
        base::strsplit(., split = ",") %>%
        base::unlist(.)
      
      # GFCs of new data set
      c_GFCs <- dplyr::filter(hcobject[["integrated_output"]][["GFC_all_layers"]], Gene %in% genes)
      
      if(base::is.vector(c_GFCs)){
        c_GFC_means <- cGFCs
      } else {
        c_GFC_means <- base::apply(c_GFCs[, base::c(1:(base::ncol(c_GFCs)-1))] %>% 
                                     base::as.data.frame(), 2, base::mean)
      }
      mat_heatmap <- base::rbind(mat_heatmap, c_GFC_means)
    }
    base::rownames(mat_heatmap) <- c_df$color
  }
  
  base::colnames(mat_heatmap) <- base::colnames(hcobject[["integrated_output"]][["GFC_all_layers"]])[1:(base::ncol(hcobject[["integrated_output"]][["GFC_all_layers"]])-1)]
  
  if(!base::is.null(col_order)){
    mat_heatmap <- mat_heatmap %>% base::as.data.frame()
    mat_heatmap <- mat_heatmap[, col_order] %>% base::as.matrix()
  }
  
  # --- 3. Prepare Enrichment Data ---
  
  enrichment_entries <- list()
  
  # Helper logic to extract enrichment data
  if(!base::is.null(row_order)){
    target_clusters <- row_order
  } else {
    target_clusters <- base::unique(c_df$color) 
  }
  
  if (base::length(user_enrichment_slots) > 0) {
    base_palette <- ggsci::pal_d3(palette = "category20")(20)
    for (slot_idx in base::seq_along(user_enrichment_slots)) {
      slot_item <- user_enrichment_slots[[slot_idx]]
      slot_df <- slot_item$data
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
          h <- suppressWarnings(base::as.numeric(dplyr::first(tmp$hits)))
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

      pal_offset <- ((slot_idx - 1) * 3) %% base::length(base_palette)
      slot_cols <- base_palette[((base::seq_len(base::ncol(slot_mat)) + pal_offset - 1) %% base::length(base_palette)) + 1]
      enrichment_entries[[base::length(enrichment_entries) + 1]] <- list(
        slot = slot_item$slot,
        label = slot_item$label,
        mat = slot_mat,
        hits = slot_hits,
        cell_types = base::colnames(slot_mat),
        colors = slot_cols
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
  if(!base::is.null(row_order)){
    cluster_colors <- base::factor(row_order)
    base::names(cluster_colors) <- row_order
    c_df <- c_df[base::match(row_order, c_df$color),]
  } else {
    cluster_colors <- base::factor(c_df$color)
    base::names(cluster_colors) <- c_df$color
    row_order <- base::unique(c_df$color)
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

  max_chars <- if (base::is.null(module_labels)) {
    1
  } else {
    base::max(base::nchar(module_labels), na.rm = TRUE)
  }

  n_rows <- base::nrow(c_df)
  preset_scale_font <- switch(module_label_preset,
                              compact = 0.82,
                              balanced = 1.0,
                              presentation = 1.2,
                              auto = 1.0)
  preset_scale_width <- switch(module_label_preset,
                               compact = 0.95,
                               balanced = 1.0,
                               presentation = 1.15,
                               auto = 1.0)
  preset_scale_pt <- switch(module_label_preset,
                            compact = 0.85,
                            balanced = 1.0,
                            presentation = 1.15,
                            auto = 1.0)

  if (!user_set_module_label_fontsize) {
    if (base::is.null(module_labels)) {
      base_font <- 8
    } else if (user_set_module_box_width_cm) {
      # Respect user-defined box width and fit text into it.
      base_font <- base::max(
        5,
        base::min(
          10,
          base::floor(((module_box_width_cm * 10) - 1.5) / base::max(1, max_chars) * 1.55)
        )
      )
    } else {
      # Auto-size by number of modules and label length to avoid overlap.
      fontsize_by_rows <- if (n_rows <= 10) {
        9
      } else if (n_rows <= 15) {
        8
      } else if (n_rows <= 25) {
        7.5
      } else if (n_rows <= 40) {
        6.8
      } else {
        6
      }
      fontsize_by_chars <- base::max(5, base::floor(16 / base::max(1, max_chars)))
      base_font <- base::max(5, base::min(10, fontsize_by_rows, fontsize_by_chars))
    }
    module_label_fontsize <- base::max(5, base::min(10, base::round(base_font * preset_scale_font, 1)))
  }

  if (!user_set_module_box_width_cm) {
    if (base::is.null(module_labels) || module_label_mode == "legacy") {
      base_width <- 0.5
    } else {
      # Approximate text width in cm; expand with longer labels automatically.
      base_width <- (max_chars * module_label_fontsize * 0.07) + 0.22
    }
    module_box_width_cm <- base::max(
      0.62,
      base::min(
        4.8,
        base_width * preset_scale_width
      )
    )
  }

  if (!user_set_module_label_pt_size) {
    base_pt <- if (n_rows <= 10) {
      0.34
    } else if (n_rows <= 15) {
      0.30
    } else if (n_rows <= 25) {
      0.26
    } else if (n_rows <= 40) {
      0.22
    } else {
      0.18
    }
    char_penalty <- 1 + base::max(0, max_chars - 2) * 0.12
    module_label_pt_size <- base::max(
      0.14,
      base::min(
        0.55,
        (base_pt * preset_scale_pt) / char_penalty
      )
    )
  }

  if (!user_set_gene_count_fontsize) {
    # Keep numeric gene-count text aligned with module label size by default.
    gene_count_fontsize <- module_label_fontsize
  }

  if (!user_set_gene_count_pt_size) {
    if (gene_count_renderer == "pch") {
      gene_count_pt_size <- module_label_pt_size
    } else {
      gene_count_pt_size <- 0.5
    }
  }

  module_box_anno <- ComplexHeatmap::anno_simple(
    row_order,
    col = cluster_colors,
    pch = module_labels,
    pt_gp = grid::gpar(col = module_label_color, fontsize = module_label_fontsize, fontface = "bold"),
    pt_size = grid::unit(module_label_pt_size, "snpc"),
    simple_anno_size = grid::unit(module_box_width_cm, "cm"),
    gp = grid::gpar(col = "black"),
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

  base_row_width_cm <- module_box_width_cm +
    if (show_gene_text) 1.2 else 0 +
    if (show_gene_bar) 2.5 else 0 +
    0.8
  enrichment_slot_bar_width_cm <- 2.2
  enrichment_slot_gap_cm <- 0.15
  n_enrichment_slots <- base::length(enrichment_entries)
  enrichment_total_width_cm <- if (n_enrichment_slots > 0) {
    (n_enrichment_slots * enrichment_slot_bar_width_cm) +
      ((n_enrichment_slots - 1) * enrichment_slot_gap_cm)
  } else {
    0
  }
  row_annotation_total_width_cm <- base_row_width_cm + enrichment_total_width_cm
  
  # --- 5. Assemble HeatmapAnnotation (Row) ---
  
  lgd_list <- list()
  ha_annos <- list(
    modules = module_box_anno,
    `# genes` = gene_text_anno,
    genes = gene_bar_anno
  )
  if (n_enrichment_slots > 0) {
    for (i in base::seq_along(enrichment_entries)) {
      entry <- enrichment_entries[[i]]
      anno_name <- entry$label
      ha_annos[[anno_name]] <- ComplexHeatmap::anno_barplot(
        entry$mat,
        width = grid::unit(enrichment_slot_bar_width_cm, "cm"),
        gp = grid::gpar(fill = entry$colors, col = entry$colors),
        baseline = 0,
        which = "row"
      )
      lgd_list[[base::length(lgd_list) + 1]] <- ComplexHeatmap::Legend(
        labels = entry$cell_types,
        title = anno_name,
        legend_gp = grid::gpar(col = entry$colors),
        type = "points",
        pch = 15
      )
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
  
  
  # --- 6. Setup Column Annotations ---
  
  anno_list <- NULL
  
  if(!base::length(column_anno_categorical) == 0){
    for(a in 1:base::length(column_anno_categorical)){
      base::set.seed(a)
      tmp_colour <- grDevices::colorRampPalette(c("#332288", "#117733", "#44aa99", "#88ccee", "#cc6677", "#aa4499", "#882255"))(base::ncol(column_anno_categorical[[a]]))
      
      if(cat_as_bp[a] == T){
        column_anno_categorical[[a]][base::is.na(column_anno_categorical[[a]])] <- 0
        current_anno <- ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_barplot(column_anno_categorical[[a]]%>%base::as.matrix(),
                                                                                                  width = grid::unit(2, "cm"),
                                                                                                  gp = grid::gpar(fill = tmp_colour, col = tmp_colour)),
                                                          which = "column",
                                                          height = grid::unit(1, "cm"),
                                                          annotation_name_side = "right",
                                                          gap = grid::unit(2, "mm"),
                                                          annotation_name_rot = 0,
                                                          annotation_name_gp = grid::gpar(fontsize = 8),
                                                          annotation_label = base::names(column_anno_categorical)[a])
      } else {
        current_anno <- ComplexHeatmap::HeatmapAnnotation(col_anno = ComplexHeatmap::anno_lines(column_anno_categorical[[a]] %>% base::as.matrix(),
                                                                                                width = grid::unit(2, "cm"),
                                                                                                gp = grid::gpar(col = tmp_colour),
                                                                                                add_points = TRUE,
                                                                                                pt_gp = grid::gpar(col = tmp_colour), pch = 16),
                                                          which = "column",
                                                          height = grid::unit(1, "cm"),
                                                          annotation_name_side = "right",
                                                          gap = grid::unit(2, "mm"),
                                                          annotation_name_rot = 0,
                                                          annotation_name_gp = grid::gpar(fontsize = 8),
                                                          annotation_label = base::names(column_anno_categorical)[a])
      }
      
      if(base::is.null(anno_list)){
        anno_list <- current_anno
      } else {
        anno_list <- ComplexHeatmap::add_heatmap(anno_list, current_anno, direction = c("vertical"))
      }
      
      lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(labels = base::colnames(column_anno_categorical[[a]]%>%base::as.matrix()), 
                                                                      title = base::names(column_anno_categorical)[a],
                                                                      legend_gp = grid::gpar(col = tmp_colour),
                                                                      type = "points", pch = 15))
    }
  }
  
  if(!base::length(column_anno_numerical) == 0){
    for(a in 1:base::length(column_anno_numerical)){
      tmp_col_anno_2 <- column_anno_numerical[[a]]
      tmp_col_anno_2 <- tmp_col_anno_2[base::colnames(mat_heatmap)]
      
      current_anno <- ComplexHeatmap::HeatmapAnnotation(cont_anno = ComplexHeatmap::anno_boxplot(tmp_col_anno_2, height = grid::unit(1, "cm")),
                                                        which = "column",
                                                        annotation_name_side = "right",
                                                        gap = grid::unit(2, "mm"),
                                                        annotation_name_rot = 0,
                                                        annotation_name_gp = grid::gpar(fontsize = 8),
                                                        annotation_label = base::names(column_anno_numerical)[a], show_legend = F)
      
      if(base::is.null(anno_list)){
        anno_list <- current_anno
      } else {
        anno_list <- ComplexHeatmap::add_heatmap(anno_list, current_anno, direction = c("vertical"))
      }
    }
  }
  
  all_conditions <- NULL
  for(setnum in 1:base::length(hcobject[["layers"]])){
    all_conditions <- base::c(all_conditions, base::as.character(dplyr::pull(hcobject[["data"]][[base::paste0("set", setnum, "_anno")]], hcobject[["global_settings"]][["voi"]])))
  }
  all_conditions <- base::table(all_conditions) %>%
    base::as.data.frame() %>%
    dplyr::filter(., all_conditions %in% base::colnames(mat_heatmap))
  all_conditions <- all_conditions[base::match(base::colnames(mat_heatmap), base::as.character(all_conditions$all_conditions)),]
  all_conditions <- base::paste0(all_conditions$all_conditions, "  [", all_conditions$Freq, "]")
  
  if(base::is.null(anno_list)){
    anno_list <- ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions))
  } else {
    anno_list <- ComplexHeatmap::add_heatmap(anno_list, ComplexHeatmap::columnAnnotation(groups = ComplexHeatmap::anno_text(all_conditions)), direction = c("vertical"))
  }
  
  
  # --- 7. Plotting and Output ---
  
  Cairo::Cairo(file = paste0(hcobject[["working_directory"]][["dir_output"]], hcobject[["global_settings"]][["save_folder"]], "/", file_name),
               width = pdf_width,
               height = pdf_height,
               pointsize = pdf_pointsize,
               dpi = pdf_dpi,
               type = "pdf",
               units = "in")

  # Keep module-expression tiles square-like for publication consistency.
  n_heat_rows <- base::nrow(mat_heatmap)
  n_heat_cols <- base::ncol(mat_heatmap)
  cell_size_mm <- 5
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
  cell_size_mm <- cell_size_mm * overall_plot_scale
  hm_width <- grid::unit(base::max(30 * overall_plot_scale, n_heat_cols * cell_size_mm), "mm")
  hm_height <- grid::unit(base::max(44 * overall_plot_scale, n_heat_rows * cell_size_mm), "mm")
  row_dend_width_mm <- if (isTRUE(cluster_rows)) {
    base::max(20, base::min(36, n_heat_rows * 1.0))
  } else {
    8
  }
  row_dend_width_mm <- row_dend_width_mm * overall_plot_scale
  column_dend_height_mm <- if (isTRUE(cluster_columns)) {
    base::max(24, base::min(44, n_heat_cols * 4.0))
  } else {
    8
  }
  column_dend_height_mm <- column_dend_height_mm * overall_plot_scale
  gfc_palette <- grDevices::colorRampPalette(gfc_colors)(51)
  gfc_col_fun <- circlize::colorRamp2(
    seq(gfc_scale_limits[1], gfc_scale_limits[2], length.out = base::length(gfc_palette)),
    gfc_palette
  )
  
  build_hm <- function(row_dend_mm,
                       col_dend_mm,
                       column_name_max_cm,
                       fixed_size = TRUE,
                       body_width_mm = NULL,
                       body_height_mm = NULL) {
    hm_args <- list(
      mat_heatmap,
      name = "GFC",
      right_annotation = ha,
      col = gfc_col_fun,
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      cluster_columns = cluster_columns,
      cluster_rows = cluster_rows_for_heatmap,
      row_dend_reorder = FALSE,
      column_names_rot = 90,
      column_names_centered = FALSE,
      row_dend_width = grid::unit(row_dend_mm, "mm"),
      column_dend_height = grid::unit(col_dend_mm, "mm"),
      column_names_max_height = grid::unit(column_name_max_cm, "cm"),
      show_row_names = show_module_color_names,
      row_names_gp = grid::gpar(fontsize = 8 * overall_plot_scale),
      column_names_gp = grid::gpar(fontsize = 10 * overall_plot_scale),
      rect_gp = grid::gpar(col = "black"),
      show_heatmap_legend = TRUE,
      heatmap_legend_param = list(
        title = "GFC",
        at = gfc_scale_breaks,
        labels = gfc_scale_labels,
        labels_gp = grid::gpar(fontfamily = "mono")
      ),
      column_km = k
    )
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

  hm_pdf <- build_hm(
    row_dend_mm = row_dend_width_mm,
    col_dend_mm = column_dend_height_mm,
    column_name_max_cm = 22,
    fixed_size = TRUE,
    body_width_mm = base::as.numeric(hm_width),
    body_height_mm = base::as.numeric(hm_height)
  )

  anno_list_pdf <- if (base::is.null(anno_list)) {
    hm_pdf
  } else {
    ComplexHeatmap::add_heatmap(hm_pdf, anno_list, direction = c("vertical"))
  }

  draw_padding_pdf_mm <- c(26, 18, 34, 36) * overall_plot_scale
  draw_padding <- grid::unit(draw_padding_pdf_mm, "mm")

  hm_pdf_drawn <- ComplexHeatmap::draw(
    anno_list_pdf,
    merge_legends = TRUE,
    annotation_legend_list = lgd_list,
    show_heatmap_legend = TRUE,
    show_annotation_legend = TRUE,
    padding = draw_padding
  )

  grDevices::dev.off()

  # For knitr/RStudio chunk devices keep layout compact to avoid clipping.
  dev_in <- tryCatch(grDevices::dev.size("in"), error = function(e) c(8, 6))
  dev_mm <- dev_in * 25.4
  small_device <- (dev_in[1] < 9.5) || (dev_in[2] < 7.0)
  row_dend_screen_mm <- if (small_device) {
    base::max(8, base::min(16, row_dend_width_mm * 0.55))
  } else {
    base::max(10, base::min(24, row_dend_width_mm * 0.8))
  }
  col_dend_screen_mm <- if (small_device) {
    base::max(8, base::min(16, column_dend_height_mm * 0.45))
  } else {
    base::max(10, base::min(28, column_dend_height_mm * 0.75))
  }
  max_col_chars <- if (!base::is.null(base::colnames(mat_heatmap)) && base::length(base::colnames(mat_heatmap)) > 0) {
    base::max(base::nchar(base::colnames(mat_heatmap)), na.rm = TRUE)
  } else {
    10
  }
  # Conservative estimate for vertical x-label space (in cm), capped by device height.
  label_est_cm <- (max_col_chars * (10 * overall_plot_scale) * 0.022) + 0.6
  label_cap_cm <- (dev_mm[2] * if (small_device) 0.24 else 0.28) / 10
  col_label_max_cm_screen <- base::max(
    2.8,
    base::min(
      if (small_device) 6.5 else 8.0,
      label_est_cm,
      label_cap_cm
    )
  )
  bottom_pad_mm <- base::max(
    if (small_device) 12 else 16,
    base::min(34, col_label_max_cm_screen * 3.0)
  )
  draw_padding_screen_mm <- if (small_device) {
    c(8, 10, bottom_pad_mm, 12) * overall_plot_scale
  } else {
    c(12, 12, bottom_pad_mm, 16) * overall_plot_scale
  }
  draw_padding_screen <- grid::unit(draw_padding_screen_mm, "mm")
  # Keep square tiles on knitr/RStudio devices while using more of the available area.
  # This avoids unnecessary shrinking of the main heatmap in notebook/html previews.
  base_body_w_cap_mm <- if (small_device) dev_mm[1] * 0.78 else dev_mm[1] * 0.82
  base_body_h_cap_mm <- if (small_device) dev_mm[2] * 0.80 else dev_mm[2] * 0.84
  row_anno_width_cm <- row_annotation_total_width_cm
  legend_reserve_mm <- if (base::length(lgd_list) == 0) 28 else 40
  available_body_w_mm <- dev_mm[1] -
    (row_anno_width_cm * 10) -
    row_dend_screen_mm -
    draw_padding_screen_mm[2] -
    draw_padding_screen_mm[4] -
    legend_reserve_mm
  available_body_h_mm <- dev_mm[2] -
    (col_label_max_cm_screen * 10) -
    col_dend_screen_mm -
    draw_padding_screen_mm[1] -
    draw_padding_screen_mm[3] -
    8
  body_w_cap_mm <- base::min(
    dev_mm[1] * 0.92,
    base_body_w_cap_mm * overall_plot_scale,
    base::max(22, available_body_w_mm)
  )
  body_h_cap_mm <- base::min(
    dev_mm[2] * 0.92,
    base_body_h_cap_mm * overall_plot_scale,
    base::max(26, available_body_h_mm)
  )
  min_cell_mm <- (if (small_device) 3.0 else 3.4) * overall_plot_scale
  target_cell_mm <- (if (small_device) 8.0 else 8.8) * overall_plot_scale
  cell_size_screen_mm <- base::max(
    min_cell_mm,
    base::min(
      target_cell_mm,
      body_w_cap_mm / base::max(1, n_heat_cols),
      body_h_cap_mm / base::max(1, n_heat_rows)
    )
  )
  hm_screen_w_mm <- base::max(22 * overall_plot_scale, n_heat_cols * cell_size_screen_mm)
  hm_screen_h_mm <- base::max(26 * overall_plot_scale, n_heat_rows * cell_size_screen_mm)

  hm_screen <- build_hm(
    row_dend_mm = row_dend_screen_mm,
    col_dend_mm = col_dend_screen_mm,
    column_name_max_cm = col_label_max_cm_screen,
    fixed_size = TRUE,
    body_width_mm = hm_screen_w_mm,
    body_height_mm = hm_screen_h_mm
  )

  anno_list_screen <- if (base::is.null(anno_list)) {
    hm_screen
  } else {
    ComplexHeatmap::add_heatmap(hm_screen, anno_list, direction = c("vertical"))
  }

  hm_w_lgd <- ComplexHeatmap::draw(
    anno_list_screen,
    merge_legends = TRUE,
    annotation_legend_list = lgd_list,
    show_heatmap_legend = TRUE,
    show_annotation_legend = TRUE,
    padding = draw_padding_screen
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
    idx <- base::which(base::as.character(c_df$color) == cl)
    if (base::length(idx) == 0) {
      return(NULL)
    }
    gene_n_value <- c_df$gene_n[[idx[[1]]]]
    if (base::is.null(gene_n_value) ||
        base::length(gene_n_value) == 0 ||
        base::all(base::is.na(gene_n_value))) {
      return(NULL)
    }
    genes <- base::trimws(
      base::unlist(
        base::strsplit(base::as.character(gene_n_value), ",", fixed = TRUE),
        use.names = FALSE
      )
    )
    genes <- genes[genes != "" & !base::is.na(genes)]
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
  module_gene_list_file <- base::paste0(
    hcobject[["working_directory"]][["dir_output"]],
    hcobject[["global_settings"]][["save_folder"]],
    "/Module_Gene_List.xlsx"
  )
  tryCatch(
    {
      openxlsx::write.xlsx(
        x = list(module_gene_list = module_gene_list_tbl),
        file = module_gene_list_file,
        overwrite = TRUE
      )
    },
    error = function(e) {
      warning(
        "Could not write Module_Gene_List.xlsx: ",
        base::conditionMessage(e)
      )
    }
  )
  hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]] <<- module_label_map
  hcobject[["integrated_output"]][["cluster_calc"]][["module_gene_list"]] <<- module_gene_list_tbl
  hcobject[["satellite_outputs"]][["module_gene_list"]] <<- module_gene_list_tbl
  hcobject[["integrated_output"]][["cluster_calc"]][["module_label_mode"]] <<- module_label_mode
  hcobject[["integrated_output"]][["cluster_calc"]][["module_label_numbering"]] <<- module_label_numbering
  hcobject[["integrated_output"]][["cluster_calc"]][["module_label_fontsize"]] <<- module_label_fontsize
  hcobject[["integrated_output"]][["cluster_calc"]][["module_label_pt_size"]] <<- module_label_pt_size
  hcobject[["integrated_output"]][["cluster_calc"]][["module_box_width_cm"]] <<- module_box_width_cm
  module_box_to_cell_ratio <- suppressWarnings((as.numeric(module_box_width_cm) * 10) / as.numeric(cell_size_screen_mm))
  if (!base::is.finite(module_box_to_cell_ratio) || module_box_to_cell_ratio <= 0) {
    module_box_to_cell_ratio <- NULL
  }
  hcobject[["integrated_output"]][["cluster_calc"]][["module_box_to_cell_ratio"]] <<- module_box_to_cell_ratio
  hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cell_size_mm"]] <<- cell_size_screen_mm
  hcobject[["integrated_output"]][["cluster_calc"]][["gene_count_fontsize"]] <<- gene_count_fontsize
  hcobject[["integrated_output"]][["cluster_calc"]][["gene_count_renderer"]] <<- gene_count_renderer
  hcobject[["integrated_output"]][["cluster_calc"]][["gene_count_pt_size"]] <<- gene_count_pt_size
  hcobject[["integrated_output"]][["cluster_calc"]][["gfc_colors"]] <<- gfc_colors
  hcobject[["integrated_output"]][["cluster_calc"]][["gfc_scale_limits"]] <<- gfc_scale_limits
  hcobject[["integrated_output"]][["cluster_calc"]][["overall_plot_scale"]] <<- overall_plot_scale
  if (module_label_mode == "prefix") {
    hcobject[["integrated_output"]][["cluster_calc"]][["module_prefix"]] <<- module_prefix
  } else {
    hcobject[["integrated_output"]][["cluster_calc"]][["module_prefix"]] <<- NULL
  }

  if(return_HM){
    hcobject[["integrated_output"]][["cluster_calc"]][["heatmap_cluster"]] <<- hm_w_lgd
  }
}
