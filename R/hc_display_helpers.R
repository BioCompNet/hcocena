.hc_display_htmlwidget <- function(x) {
  display_mode <- base::tolower(base::as.character(
    getOption("hcocena.htmlwidget_display", "auto")
  )[[1]])
  if (!display_mode %in% c("auto", "inline", "viewer")) {
    display_mode <- "auto"
  }

  if (!identical(display_mode, "viewer") &&
      isTRUE(getOption("knitr.in.progress", FALSE)) &&
      requireNamespace("knitr", quietly = TRUE) &&
      isTRUE(knitr::is_html_output())) {
    rendered <- tryCatch(
      knitr::knit_print(x, options = knitr::opts_current$get()),
      error = function(e) NULL
    )
    if (inherits(rendered, "knit_asis")) {
      meta <- attr(rendered, "knit_meta", exact = TRUE)
      if (!is.null(meta)) {
        knitr::knit_meta_add(meta)
      }
      cat(as.character(rendered), sep = "\n")
      return(invisible(x))
    }
  }

  # Emit a raw text/html representation only for front-ends that consume HTML
  # written to stdout (e.g. some Jupyter/IRkernel setups). This is strictly
  # opt-in: RStudio notebooks and the console must NOT take this path, otherwise
  # the widget is dumped as raw HTML text below the chunk instead of rendering.
  emit_html <- identical(display_mode, "inline") ||
    isTRUE(getOption("hcocena.htmlwidget_emit_html", FALSE))

  if (isTRUE(emit_html) &&
      requireNamespace("repr", quietly = TRUE)) {
    displayed <- tryCatch(
      {
        html <- repr::repr_html(x)
        if (base::is.character(html) &&
            base::length(html) > 0 &&
            base::any(base::nzchar(html))) {
          base::cat(html, sep = "\n")
          TRUE
        } else {
          FALSE
        }
      },
      error = function(e) FALSE
    )
    if (isTRUE(displayed)) {
      return(invisible(x))
    }
  }

  # Default path. Dispatch through the `print` generic exactly like a bare
  # `print(widget)` in a chunk would. RStudio's notebook execution renders
  # widgets inline only when they go through the generic dispatch (this works
  # even from inside a function); calling `print.htmlwidget` directly (e.g. via
  # getFromNamespace) bypasses that hook, so the widget ends up in the Viewer
  # pane instead of inline below the chunk. View defaults to interactive(), the
  # same as a top-level `print(widget)`.
  print(x)
  invisible(x)
}

.hc_display_object <- function(x, row.names = TRUE, col.names = TRUE) {
  if (base::is.null(x)) {
    return(base::invisible(NULL))
  }

  if (inherits(x, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(x, merge_legends = TRUE, newpage = TRUE)
    return(base::invisible(x))
  }

  if (inherits(x, "hc_llm_heatmap_plot")) {
    print(x)
    return(base::invisible(x))
  }

  if (inherits(x, "gtable")) {
    grid::grid.newpage()
    grid::grid.draw(x)
    return(base::invisible(x))
  }

  if (inherits(x, "grob")) {
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
    return(.hc_display_htmlwidget(x))
  }

  if (inherits(x, "ggplot") || inherits(x, "patchwork")) {
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
