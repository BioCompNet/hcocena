.hc_display_object <- function(x, row.names = TRUE, col.names = TRUE) {
  if (base::is.null(x)) {
    return(base::invisible(NULL))
  }

  if (inherits(x, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(x, merge_legends = TRUE, newpage = TRUE)
    return(base::invisible(x))
  }

  if (inherits(x, "gtable")) {
    grid::grid.newpage()
    grid::grid.draw(x)
    return(base::invisible(x))
  }

  if (base::is.list(x) && !base::is.null(x$gtable) && inherits(x$gtable, "gtable")) {
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    return(base::invisible(x))
  }

  if (inherits(x, "htmlwidget")) {
    utils::getFromNamespace("print.htmlwidget", "htmlwidgets")(x)
    return(base::invisible(x))
  }

  if (ggplot2::is.ggplot(x) || inherits(x, "patchwork")) {
    graphics::plot(x)
    return(base::invisible(x))
  }

  if (base::is.data.frame(x) || base::is.matrix(x)) {
    utils::write.table(
      x = base::as.data.frame(x),
      file = "",
      sep = "\t",
      quote = FALSE,
      row.names = row.names,
      col.names = col.names
    )
    return(base::invisible(x))
  }

  base::writeLines(utils::capture.output(utils::str(x, max.level = 1)))
  base::invisible(x)
}
