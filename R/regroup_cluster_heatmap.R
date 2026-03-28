#' Change Grouping Parameter
#' 
#' The variable by which the samples are grouped and based on which the GFCs are calculated can be changed and the cluster heatmap will be replotted. 
#'  This is particularly useful in cases where different variables are potential candidates for driving the genes' expression changes in the data 
#'  and an explorative approach is required to decide on the most suitable one.
#'  Note that any previously generated column annotation will not be plotted, since the grouping will change.
#'  If you eventually decide on another grouping variable, please run the analysis again entirely with the changed "voi" from the very beginning. 
#' @param group_by A string giving the new grouping variable. This must be a name of a column present in all annotation files.
#' @param col_order Defines the order in which the sample groups (conditions) appear in the heatmap. 
#'  Accepts a vector of strings giving the conditions in their desired order.
#'  If `NULL` and `cluster_columns = FALSE`, the column order from the previous
#'  main hCoCena heatmap is reused when available.
#'  If `cluster_columns = TRUE`, this order is overwritten by clustering.
#' @param row_order Like col_order but with cluster names.
#' @param cluster_columns A Boolean, whether or not to cluster the columns of
#'  the heatmap. Default is FALSE so the main hCoCena column order is preserved.
#' @param cluster_rows Like cluster_columns but for rows.
#' @export


change_grouping_parameter <- function(group_by, col_order = NULL, cluster_columns = FALSE, row_order = NULL, cluster_rows = TRUE){
  
  # check if grouping variables are present:
  if(base::length(group_by) == 1){
    for(i in 1:base::length(hcobject[["layers"]])){
      if(!group_by %in% base::colnames(hcobject[["data"]][[base::paste0("set", i, "_anno")]])){
        stop("Grouping variable not present as column name in all annotation tables.")
      }
    }
  }else{
    stop("Please provide only one grouping variable that is a column name present in ALL annotation files.")
  }
 
  
  # to store temporary GFCs per layer:
  sep_GFCs_list <- list()
  tmp_data <- list()
  layer_ids <- base::names(hcobject[["layers"]])
  if (base::is.null(layer_ids) || base::length(layer_ids) != base::length(hcobject[["layers"]])) {
    layer_ids <- base::paste0("set", base::seq_along(hcobject[["layers"]]))
  }
  old_control <- hcobject[["global_settings"]][["control"]]
  old_voi <- hcobject[["global_settings"]][["voi"]]
  old_gfc_all_layers <- hcobject[["integrated_output"]][["GFC_all_layers"]]
  old_data <- hcobject[["data"]]
  old_layer_specific <- hcobject[["layer_specific_outputs"]]
  sat_outputs_present <- "satellite_outputs" %in% base::names(hcobject) &&
    !base::is.null(hcobject[["satellite_outputs"]])
  old_col_annos_categorical <- if (sat_outputs_present) {
    hcobject[["satellite_outputs"]][["column_annos_categorical"]]
  } else {
    NULL
  }
  old_col_annos_numerical <- if (sat_outputs_present) {
    hcobject[["satellite_outputs"]][["column_annos_numerical"]]
  } else {
    NULL
  }
  old_module_label_mode <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["module_label_mode"]], error = function(e) NULL)
  old_module_label_numbering <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["module_label_numbering"]], error = function(e) NULL)
  old_module_label_fontsize <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["module_label_fontsize"]], error = function(e) NULL)
  old_module_label_pt_size <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["module_label_pt_size"]], error = function(e) NULL)
  old_module_box_width_cm <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["module_box_width_cm"]], error = function(e) NULL)
  old_gene_count_fontsize <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["gene_count_fontsize"]], error = function(e) NULL)
  old_gene_count_renderer <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["gene_count_renderer"]], error = function(e) NULL)
  old_gene_count_pt_size <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["gene_count_pt_size"]], error = function(e) NULL)
  old_gfc_colors <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["gfc_colors"]], error = function(e) NULL)
  old_gfc_scale_limits <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["gfc_scale_limits"]], error = function(e) NULL)
  old_overall_plot_scale <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["overall_plot_scale"]], error = function(e) NULL)
  on.exit({
    hcobject[["global_settings"]][["control"]] <<- old_control
    hcobject[["global_settings"]][["voi"]] <<- old_voi
    hcobject[["integrated_output"]][["GFC_all_layers"]] <<- old_gfc_all_layers
    hcobject[["data"]] <<- old_data
    hcobject[["layer_specific_outputs"]] <<- old_layer_specific
    if (sat_outputs_present) {
      hcobject[["satellite_outputs"]][["column_annos_categorical"]] <<- old_col_annos_categorical
      hcobject[["satellite_outputs"]][["column_annos_numerical"]] <<- old_col_annos_numerical
    }
  }, add = TRUE)

  # iterate over data sets:  
  for(i in 1:base::length(hcobject[["layers"]])){
    
    meta_data <- hcobject[["data"]][[base::paste0("set", i, "_anno")]]
    
    meta_data$regrouped <- base::paste0(meta_data[[group_by]],"_", hcobject[["layers_names"]][i])
    
    # mark controls:
    if(!hcobject[["global_settings"]][["control"]] == "none"){
      # temporarily change control to none, since controls in new grouping variable is not known:
      hcobject[["global_settings"]][["control"]] <<- "none"
      print("Regrouped GFCs will be calculated without reference to controls.")
    }
    
    sep_GFCs_list[[i]] <- GFC_calculation(info_dataset = meta_data, grouping_v = "regrouped", x = i)
    tmp_data[[base::paste0("set", i, "_anno")]] <- meta_data
  }
  
  sep_GFCs <- purrr::reduce(sep_GFCs_list, dplyr::full_join, by = "Gene")
  sep_GFCs[base::is.na(sep_GFCs)] <- -hcobject[["global_settings"]][["range_GFC"]]
  gene_col <- sep_GFCs$Gene
  sep_GFCs$Gene <- NULL
  sep_GFCs$Gene <- gene_col
  
  
  # numerical annotation:
  
    # to be added
  
  # categorical annotation:
  
    # to be added

  tmp_data_full <- old_data
  for (nm in base::names(tmp_data)) {
    tmp_data_full[[nm]] <- tmp_data[[nm]]
  }

  tmp_layer_specific <- old_layer_specific
  if (base::is.null(tmp_layer_specific) || !base::is.list(tmp_layer_specific)) {
    tmp_layer_specific <- list()
  }
  for (i in base::seq_along(layer_ids)) {
    lid <- layer_ids[[i]]
    current_slot <- if (!base::is.null(tmp_layer_specific[[lid]])) {
      tmp_layer_specific[[lid]]
    } else if (base::length(tmp_layer_specific) >= i && !base::is.null(tmp_layer_specific[[i]])) {
      tmp_layer_specific[[i]]
    } else {
      list()
    }
    if (!base::is.list(current_slot)) {
      current_slot <- list()
    }
    if (base::is.null(current_slot[["part2"]]) || !base::is.list(current_slot[["part2"]])) {
      current_slot[["part2"]] <- list()
    }
    current_slot[["part2"]][["GFC_all_genes"]] <- sep_GFCs_list[[i]]
    tmp_layer_specific[[lid]] <- current_slot
  }

  if (sat_outputs_present) {
    hcobject[["satellite_outputs"]][["column_annos_categorical"]] <<- NULL
    hcobject[["satellite_outputs"]][["column_annos_numerical"]] <<- NULL
  }
  hcobject[["integrated_output"]][["GFC_all_layers"]] <<- sep_GFCs
  hcobject[["data"]] <<- tmp_data_full
  hcobject[["layer_specific_outputs"]] <<- tmp_layer_specific
  hcobject[["global_settings"]][["voi"]] <<- "regrouped"

  plot_args <- list(
    col_order = col_order,
    row_order = row_order,
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    return_HM = TRUE,
    file_name = base::paste0("module_heatmap_regrouped_", group_by, ".pdf"),
    gene_count_mode = "text"
  )

  if (!base::is.null(old_module_label_mode) && base::nzchar(base::as.character(old_module_label_mode))) {
    plot_args$module_label_mode <- old_module_label_mode
  }
  if (!base::is.null(old_module_label_numbering) && base::nzchar(base::as.character(old_module_label_numbering))) {
    plot_args$module_label_numbering <- old_module_label_numbering
  }
  if (base::is.numeric(old_module_label_fontsize) && base::length(old_module_label_fontsize) == 1 && base::is.finite(old_module_label_fontsize)) {
    plot_args$module_label_fontsize <- old_module_label_fontsize
  }
  if (base::is.numeric(old_module_label_pt_size) && base::length(old_module_label_pt_size) == 1 && base::is.finite(old_module_label_pt_size)) {
    plot_args$module_label_pt_size <- old_module_label_pt_size
  }
  if (base::is.numeric(old_module_box_width_cm) && base::length(old_module_box_width_cm) == 1 && base::is.finite(old_module_box_width_cm)) {
    plot_args$module_box_width_cm <- old_module_box_width_cm
  }
  if (base::is.numeric(old_gene_count_fontsize) && base::length(old_gene_count_fontsize) == 1 && base::is.finite(old_gene_count_fontsize)) {
    plot_args$gene_count_fontsize <- old_gene_count_fontsize
  }
  if (!base::is.null(old_gene_count_renderer) && base::nzchar(base::as.character(old_gene_count_renderer))) {
    plot_args$gene_count_renderer <- old_gene_count_renderer
  }
  if (base::is.numeric(old_gene_count_pt_size) && base::length(old_gene_count_pt_size) == 1 && base::is.finite(old_gene_count_pt_size)) {
    plot_args$gene_count_pt_size <- old_gene_count_pt_size
  }
  if (base::is.character(old_gfc_colors) && base::length(old_gfc_colors) >= 2) {
    plot_args$gfc_colors <- old_gfc_colors
  }
  if (base::is.numeric(old_gfc_scale_limits) && base::length(old_gfc_scale_limits) >= 1 && base::all(base::is.finite(old_gfc_scale_limits))) {
    plot_args$gfc_scale_limits <- old_gfc_scale_limits
  }
  if (base::is.numeric(old_overall_plot_scale) && base::length(old_overall_plot_scale) == 1 && base::is.finite(old_overall_plot_scale)) {
    plot_args$overall_plot_scale <- old_overall_plot_scale
  }

  hm <- base::do.call(plot_cluster_heatmap_new, plot_args)

  invisible(hm)

}

.hc_change_grouping_parameter_legacy_driver <- change_grouping_parameter

change_grouping_parameter <- function(group_by, col_order = NULL, cluster_columns = FALSE, row_order = NULL, cluster_rows = TRUE) {
  .hc_run_legacy_via_modern(
    "change_grouping_parameter",
    hc_change_grouping_parameter,
    group_by = group_by,
    col_order = col_order,
    cluster_columns = cluster_columns,
    row_order = row_order,
    cluster_rows = cluster_rows
  )
}
