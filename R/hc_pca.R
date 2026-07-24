#' Run And Plot PCA
#'
#' Plots one PCA for each dataset.
#' @param which One of "all" (uses all genes), "topvar" (uses top most variant genes, cannot be used before running "Data processing part I" in the main markdown),
#' 	or "network" (uses the genes present in network, cannot be run before finishing "Data processing part II" in the main markdown).
#' @param color_by `NULL` (default) auto-uses the global `variable_of_interest` (e.g. `"merged"`).
#'  Set `"none"` to draw ungrouped points. Alternatively, set this to a vector of
#'  column names (`c("name of col in annotation 1", "name of col in annotation 2", ...)`)
#'  containing the groups by which you want to color.
#' @param ellipses A Boolean. Whether or not to add ellipses to the PCA plot. For details see documentation on factoextra::fviz_pca_ind. Default is FALSE.
#' @param cols Optional color palette passed to `factoextra::fviz_pca_ind()`.
#' @export


.hc_pca_prepare_expression <- function(x,
                                       layer_label,
                                       scale_features = FALSE) {
  if (base::is.null(x)) {
    stop(
      "No expression data are available for PCA in layer `",
      layer_label,
      "`.",
      call. = FALSE
    )
  }

  mat <- tryCatch(
    base::as.matrix(x),
    error = function(e) {
      stop(
        "Could not convert PCA expression data for layer `",
        layer_label,
        "` to a matrix: ",
        base::conditionMessage(e),
        call. = FALSE
      )
    }
  )
  if (base::length(base::dim(mat)) != 2 || !base::is.numeric(mat)) {
    stop(
      "PCA expression data for layer `",
      layer_label,
      "` must be a two-dimensional numeric matrix.",
      call. = FALSE
    )
  }
  if (base::nrow(mat) == 0 || base::ncol(mat) < 2) {
    stop(
      "PCA expression data for layer `",
      layer_label,
      "` must contain at least one gene and two samples.",
      call. = FALSE
    )
  }
  if (base::any(!base::is.finite(mat))) {
    stop(
      "PCA expression data for layer `",
      layer_label,
      "` contain non-finite values.",
      call. = FALSE
    )
  }

  if (base::isTRUE(scale_features)) {
    feature_sd <- base::apply(mat, 1, stats::sd)
    keep <- base::is.finite(feature_sd) & feature_sd > 0
    if (base::any(!keep)) {
      warning(
        "PCA: removed ",
        base::sum(!keep),
        " constant feature(s) from layer `",
        layer_label,
        "` before scaling.",
        call. = FALSE
      )
      mat <- mat[keep, , drop = FALSE]
    }
    if (base::nrow(mat) == 0) {
      stop(
        "PCA expression data for layer `",
        layer_label,
        "` contain no variable genes after filtering.",
        call. = FALSE
      )
    }
  }

  base::t(mat)
}


.hc_pca_individual_plot <- function(res.pca,
                                    groups,
                                    palette,
                                    ellipses,
                                    title,
                                    pointsize = 5) {
  if (base::ncol(res.pca$x) >= 2) {
    return(
      factoextra::fviz_pca_ind(
        res.pca,
        geom = "point",
        addEllipses = ellipses,
        habillage = groups,
        palette = palette,
        label = "none",
        title = title,
        pointsize = pointsize,
        invisible = "quali"
      ) +
        ggplot2::scale_shape_manual(
          values = base::rep(19, base::length(base::unique(groups)))
        ) +
        ggplot2::theme_bw()
    )
  }

  pca_var <- res.pca$sdev^2
  pc1_percent <- if (base::length(pca_var) > 0 &&
    base::sum(pca_var) > 0) {
    base::round(pca_var[[1]] / base::sum(pca_var) * 100, 1)
  } else {
    0
  }
  plot_data <- base::data.frame(
    PC1 = res.pca$x[, 1],
    PC2 = base::rep(0, base::nrow(res.pca$x)),
    Group = groups,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = PC1, y = PC2, colour = Group)
  ) +
    ggplot2::geom_point(size = pointsize, shape = 19) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::labs(
      title = title,
      x = base::paste0("PC1 (", pc1_percent, "%)"),
      y = "PC2 (0%)"
    ) +
    ggplot2::theme_bw()
}


PCA <- function(which = "all", color_by = NULL, ellipses = FALSE, cols = NULL) {
  default_color_by <- hcobject[["global_settings"]][["voi"]]
  if (is.null(default_color_by) || length(default_color_by) == 0 || is.na(default_color_by[[1]])) {
    default_color_by <- "none"
  } else {
    default_color_by <- as.character(default_color_by[[1]])
  }


  out <- base::lapply(base::seq_along(hcobject[["layers"]]), function(x) {
    if (which == "all") {
      pca_input <- .hc_pca_prepare_expression(
        hcobject[["data"]][[base::paste0("set", x, "_counts")]],
        layer_label = hcobject[["layers_names"]][x],
        scale_features = FALSE
      )
      res.pca <- stats::prcomp(pca_input, scale. = FALSE)
    } else if (which == "topvar") {
      pca_input <- .hc_pca_prepare_expression(
        hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part1"]][["topvar"]],
        layer_label = hcobject[["layers_names"]][x],
        scale_features = FALSE
      )
      res.pca <- stats::prcomp(pca_input, scale. = FALSE)
    } else if (which == "network") {
      genes <- igraph::V(hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["part2"]][["heatmap_out"]][["filt_cutoff_graph"]])$name
      dat <- hcobject[["data"]][[base::paste0("set", x, "_counts")]]
      dat <- dat[base::rownames(dat) %in% genes, , drop = FALSE]
      pca_input <- .hc_pca_prepare_expression(
        dat,
        layer_label = hcobject[["layers_names"]][x],
        scale_features = FALSE
      )
      res.pca <- stats::prcomp(pca_input, scale. = FALSE)
    } else {
      stop(
        "Invalid parameter setting: `which = \"",
        which,
        "\"` is not recognized.",
        call. = FALSE
      )
    }

    this_anno <- hcobject[["data"]][[base::paste0("set", x, "_anno")]]
    this_color_by <- color_by
    if (is.null(this_color_by) || (length(this_color_by) == 1 && is.na(this_color_by[[1]]))) {
      this_color_by <- default_color_by
    }

    if (length(this_color_by) == 1) {
      this_color_by <- as.character(this_color_by[[1]])
      if (this_color_by %in% c("auto", "default", "voi")) {
        this_color_by <- default_color_by
      }

      if (this_color_by == "none") {
        groups <- base::rep("all", base::nrow(this_anno))
      } else if (this_color_by %in% base::colnames(this_anno)) {
        groups <- dplyr::pull(this_anno, this_color_by)
      } else {
        warning(
          "PCA: color_by column '", this_color_by,
          "' not found in annotation for layer ", hcobject[["layers_names"]][x],
          ". Falling back to ungrouped points.",
          call. = FALSE
        )
        groups <- base::rep("all", base::nrow(this_anno))
      }
    } else {
      this_color_by <- as.character(this_color_by[[x]])
      if (this_color_by %in% base::colnames(this_anno)) {
        groups <- dplyr::pull(this_anno, this_color_by)
      } else {
        warning(
          "PCA: color_by column '", this_color_by,
          "' not found in annotation for layer ", hcobject[["layers_names"]][x],
          ". Falling back to ungrouped points.",
          call. = FALSE
        )
        groups <- base::rep("all", base::nrow(this_anno))
      }
    }

    if (is.null(cols)) {
      if (base::length(base::unique(groups)) > base::length(ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(8))) {
        my_palette <- grDevices::colorRampPalette(ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(8))(base::length(base::unique(groups)))
      } else {
        my_palette <- ggsci::pal_nejm(palette = base::c("default"), alpha = 1)(base::length(base::unique(groups)))
      }
    } else {
      my_palette <- cols
    }

    g <- .hc_pca_individual_plot(
      res.pca = res.pca,
      groups = groups,
      palette = my_palette,
      ellipses = ellipses,
      title = base::paste0("PCA ", hcobject[["layers_names"]][x], " ", which),
      pointsize = 5
    )

    graphics::plot(g)

    .hc_export_ggplot_file(
      file = .hc_output_file(base::paste0("PCA_", which, "_", hcobject[["layers_names"]][x], ".pdf")),
      plot = g,
      width = 10,
      height = 7
    )

    return(g)
  })


  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "pca"), out)
}

.hc_PCA_driver <- PCA

PCA <- function(which = "all", color_by = NULL, ellipses = FALSE, cols = NULL) {
  .hc_run_alias_via_modern(
    "PCA",
    hc_pca,
    which = which,
    color_by = color_by,
    ellipses = ellipses,
    cols = cols
  )
}
