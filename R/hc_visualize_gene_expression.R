#' Visualize Gene Expression
#'
#' Plots the mean expression values per condition for the given genes for each dataset as a heatmap. Values are scaled across rows.
#' @param genes A vector of strings giving the genes to be plotted.
#' @param save Bolean to determine whether the plot should be saved as PDF. TRUE by default.
#' @param name A string giving the plot title and the name for the save file (WITHOUT file ending).
#' @param width Width of the PDF, default ist 15, only change if plots overlap in PDF.
#' @param height Height of the PDF, default ist 10, only change if plots overlap in PDF.
#' @param label_map Optional named character vector mapping cluster colours to
#'  module labels (`M1`, `M2.1`, ...). When supplied, row labels show the module
#'  label instead of the raw colour. Defaults to NULL (colour, as before).
#' @noRd

.hc_visualize_gene_expression_driver <- function(genes, name = NULL, width = 15, height = 10, save = TRUE,
                                                 label_map = NULL) {
  # Filter for genes present in the network

  gtc <- .hc_gene_to_cluster_impl() %>%
    dplyr::filter(., gene %in% genes) %>%
    dplyr::filter(., !color == "white")

  if (nrow(gtc) == 0) {
    message("None of these genes are present in the network.")
    stop()
  } else {
    message(base::nrow(gtc), " out of ", base::length(genes), " requested genes are present in the network.")
  }
  genes <- gtc$gene


  # Heatmap

  plotls <- NULL
  cp <- NULL
  heatmap_count <- 0L

  for (x in base::seq_along(hcobject[["layers"]])) {
    counts <- hcobject[["data"]][[base::paste0("set", x, "_counts")]]
    genes_in_layer <- base::intersect(genes, base::rownames(counts))
    missing_genes <- base::setdiff(genes, genes_in_layer)
    if (base::length(missing_genes) > 0) {
      message(
        "Skipping ",
        base::length(missing_genes),
        " requested gene(s) absent from layer `",
        hcobject[["layers_names"]][x],
        "`."
      )
    }
    if (base::length(genes_in_layer) == 0) {
      next
    }
    exp <- counts[genes_in_layer, , drop = FALSE]

    anno <- hcobject[["data"]][[base::paste0("set", x, "_anno")]]
    missing_annotations <- base::setdiff(base::colnames(exp), base::rownames(anno))
    if (base::length(missing_annotations) > 0) {
      stop(
        "Samples from expression matrix not found in annotation for layer ",
        x,
        ": ",
        base::paste(missing_annotations, collapse = ", ")
      )
    }
    hm_anno <- anno[base::colnames(exp), , drop = FALSE] %>%
      dplyr::select(., dplyr::all_of(hcobject[["global_settings"]][["voi"]]))
    base::colnames(hm_anno) <- "voi"

    conditions <- base::sort(base::unique(hm_anno$voi[!base::is.na(hm_anno$voi)]))

    mexp <- base::lapply(conditions, function(y) {
      samples <- base::subset(hm_anno, voi == y) %>% base::rownames()
      missing_samples <- base::setdiff(samples, base::colnames(exp))
      if (base::length(missing_samples) > 0) {
        stop(
          "Samples from annotation not found in expression matrix for layer ",
          x, ": ", base::paste(missing_samples, collapse = ", ")
        )
      }
      tmp_exp <- exp[, samples, drop = FALSE]
      tmp_mexp <- base::data.frame(V1 = base::apply(tmp_exp, 1, base::mean))
      base::colnames(tmp_mexp) <- y
      return(tmp_mexp)
    }) %>% rlist::list.cbind()

    base::rownames(mexp) <- base::lapply(base::rownames(mexp), function(r) {
      color <- dplyr::filter(gtc, gene == r) %>% dplyr::pull(., "color")
      if (!base::is.null(label_map) && base::as.character(color)[1] %in% base::names(label_map)) {
        color <- base::as.character(label_map[[base::as.character(color)[1]]])
      }
      return(paste0(r, " [", color, "]"))
    }) %>% base::unlist()

    hm <- ComplexHeatmap::pheatmap(
      mat = as.matrix(mexp),
      scale = "row",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      cellwidth = 30,
      cellheight = 15,
      angle_col = "90",
      legend = if (x == 1) {
        FALSE
      } else {
        TRUE
      },
      heatmap_legend_param = list(title = "scaled mean expr."),
      main = hcobject[["layers_names"]][x],
      fontsize = 10,
      color = grDevices::colorRampPalette(base::rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG")))(51),
      treeheight_col = 25, treeheight_row = 25
    )
    plotls <- if (base::is.null(plotls)) hm else plotls + hm
    heatmap_count <- heatmap_count + 1L
  }

  if (base::is.null(plotls) || heatmap_count == 0L) {
    message("None of the requested network genes are present in any expression layer; skipping the expression heatmap.")
    return(base::invisible(NULL))
  }

  heatmap_title <- stringr::str_replace_all(string = name, pattern = "_", replacement = " ")

  if (save) {
    if (base::is.null(name)) {
      message("Cannot save the file since no unique file name was provided (See function parameter 'name').")
    } else {
      .hc_export_single_page_plot(
        file = base::paste0(
          hcobject[["working_directory"]][["dir_output"]],
          hcobject[["global_settings"]][["save_folder"]],
          "/",
          name,
          ".pdf"
        ),
        width = width,
        height = height,
        draw_fun = function() {
          ComplexHeatmap::plot.HeatmapList(
            plotls,
            column_title = heatmap_title,
            column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
          )
        }
      )
    }
  }

  graphics::plot.new()
  ComplexHeatmap::plot.HeatmapList(plotls, column_title = heatmap_title, column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"))
}


