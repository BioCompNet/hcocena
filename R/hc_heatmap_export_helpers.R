.hc_export_path_with_ext <- function(path, ext) {
  ext <- base::paste0(".", base::sub("^\\.", "", base::as.character(ext[[1]])))
  base::paste0(tools::file_path_sans_ext(path), ext)
}

.hc_export_sanitize_stem <- function(x, default = "page") {
  x <- base::as.character(x[[1]])
  if (base::is.na(x) || !base::nzchar(x)) {
    return(default)
  }
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!base::nzchar(x)) {
    x <- default
  }
  x
}

.hc_output_root <- function() {
  out_dir <- tryCatch(hcobject[["working_directory"]][["dir_output"]], error = function(e) NULL)
  save_folder <- tryCatch(hcobject[["global_settings"]][["save_folder"]], error = function(e) NULL)

  if (base::is.null(out_dir) || base::length(out_dir) == 0) {
    stop("No output directory configured in `hcobject`.")
  }
  out_dir <- base::as.character(out_dir[[1]])
  if (base::is.na(out_dir) || !base::nzchar(out_dir)) {
    stop("No output directory configured in `hcobject`.")
  }

  save_folder_value <- ""
  if (!base::is.null(save_folder) && base::length(save_folder) > 0) {
    save_folder_value <- base::as.character(save_folder[[1]])
    if (base::is.na(save_folder_value) || save_folder_value %in% c("", "FALSE", "false")) {
      save_folder_value <- ""
    }
  }

  if (!base::nzchar(save_folder_value)) {
    return(out_dir)
  }
  base::file.path(out_dir, save_folder_value)
}

.hc_output_dir <- function(...) {
  parts <- base::list(...)
  path <- .hc_output_root()
  for (part in parts) {
    if (base::is.null(part) || !base::nzchar(base::as.character(part[[1]]))) {
      next
    }
    path <- base::file.path(path, base::as.character(part[[1]]))
  }
  if (!base::dir.exists(path)) {
    base::dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

.hc_output_file <- function(filename, ..., create_dirs = TRUE) {
  dir_path <- .hc_output_root()
  parts <- base::list(...)
  for (part in parts) {
    if (base::is.null(part) || !base::nzchar(base::as.character(part[[1]]))) {
      next
    }
    dir_path <- base::file.path(dir_path, base::as.character(part[[1]]))
  }
  if (isTRUE(create_dirs) && !base::dir.exists(dir_path)) {
    base::dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  base::file.path(dir_path, filename)
}

.hc_open_pdf_device <- function(file,
                                width,
                                height,
                                pointsize = 11,
                                dpi = NULL,
                                onefile = TRUE) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    if (!base::is.null(dpi)) {
      Cairo::Cairo(
        file = file,
        width = width,
        height = height,
        pointsize = pointsize,
        dpi = dpi,
        type = "pdf",
        units = "in"
      )
    } else {
      Cairo::CairoPDF(
        file = file,
        width = width,
        height = height,
        pointsize = pointsize,
        onefile = onefile
      )
    }
  } else {
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      pointsize = pointsize,
      onefile = onefile
    )
  }
}

.hc_open_png_device <- function(file,
                                width,
                                height,
                                res = 300,
                                pointsize = 11,
                                bg = "white") {
  grDevices::png(
    filename = file,
    width = width,
    height = height,
    units = "in",
    res = res,
    pointsize = pointsize,
    bg = bg
  )
}

.hc_export_single_page_plot <- function(file,
                                        width,
                                        height,
                                        png_width = NULL,
                                        png_height = NULL,
                                        pointsize = 11,
                                        res = 300,
                                        pdf_dpi = NULL,
                                        draw_fun) {
  if (!base::is.character(file) || base::length(file) != 1 || !base::nzchar(file)) {
    stop("`file` must be a non-empty file path.")
  }
  if (!base::is.function(draw_fun)) {
    stop("`draw_fun` must be a function.")
  }

  pdf_file <- file
  png_file <- .hc_export_path_with_ext(file, "png")
  if (base::is.null(png_width)) {
    png_width <- width
  }
  if (base::is.null(png_height)) {
    png_height <- height
  }

  render_page <- function(open_device) {
    open_device()
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    draw_fun()
    invisible(NULL)
  }

  render_page(function() {
    .hc_open_pdf_device(
      file = pdf_file,
      width = width,
      height = height,
      pointsize = pointsize,
      dpi = pdf_dpi
    )
  })
  render_page(function() {
    .hc_open_png_device(
      file = png_file,
      width = png_width,
      height = png_height,
      res = res,
      pointsize = pointsize
    )
  })

  list(
    pdf = pdf_file,
    png = png_file
  )
}

.hc_export_multi_page_plot <- function(file,
                                       page_labels,
                                       width,
                                       height,
                                       pointsize = 11,
                                       res = 300,
                                       draw_page_fun,
                                       display = FALSE) {
  if (!base::is.character(file) || base::length(file) != 1 || !base::nzchar(file)) {
    stop("`file` must be a non-empty file path.")
  }
  if (!base::is.function(draw_page_fun)) {
    stop("`draw_page_fun` must be a function.")
  }

  page_labels <- base::as.character(page_labels)
  if (base::length(page_labels) == 0) {
    return(list(pdf = NULL, png = base::character(0)))
  }

  pdf_file <- file
  png_files <- stats::setNames(
    base::character(base::length(page_labels)),
    page_labels
  )

  .hc_open_pdf_device(
    file = pdf_file,
    width = width,
    height = height,
    pointsize = pointsize,
    onefile = TRUE
  )
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  for (idx in base::seq_along(page_labels)) {
    draw_page_fun(idx, page_labels[[idx]])
  }
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)

  stem <- tools::file_path_sans_ext(pdf_file)
  for (idx in base::seq_along(page_labels)) {
    png_file <- base::sprintf(
      "%s_%03d_%s.png",
      stem,
      idx,
      .hc_export_sanitize_stem(page_labels[[idx]], default = base::sprintf("page_%03d", idx))
    )
    .hc_open_png_device(
      file = png_file,
      width = width,
      height = height,
      res = res,
      pointsize = pointsize
    )
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    draw_page_fun(idx, page_labels[[idx]])
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
    png_files[[idx]] <- png_file
  }

  # Optionally replay each page on the active graphics device so the figures
  # also appear inline (e.g. under an R Markdown chunk / in a notebook), in
  # addition to the exported files. Gated to contexts where a display target
  # exists, so batch runs do not spawn a stray Rplots.pdf.
  if (isTRUE(display) &&
    (base::interactive() || isTRUE(base::getOption("knitr.in.progress", FALSE)))) {
    for (idx in base::seq_along(page_labels)) {
      draw_page_fun(idx, page_labels[[idx]])
    }
  }

  list(
    pdf = pdf_file,
    png = png_files
  )
}

.hc_export_ggplot_file <- function(file,
                                   plot,
                                   width,
                                   height,
                                   pointsize = 11,
                                   res = 300) {
  if (!inherits(plot, "ggplot")) {
    stop("`plot` must be a ggplot object.")
  }
  .hc_export_single_page_plot(
    file = file,
    width = width,
    height = height,
    pointsize = pointsize,
    res = res,
    draw_fun = function() {
      .hc_display_object(plot)
    }
  )
}

# Save a ggplot to both a (cairo) PDF and a PNG companion in one call.
#
# Drop-in replacement for `ggplot2::ggsave(filename = "...pdf", ...)`: the PDF
# is written with `cairo_pdf` (matching the rest of the package) and a PNG of
# the same dimensions is written next to it at `res` dpi. Extra `...` arguments
# are forwarded to both saves; any `device` is ignored (PDF forces cairo_pdf,
# PNG infers from the `.png` extension). The PNG failing only warns, so a PDF is
# still produced.
.hc_ggsave_pdf_png <- function(filename,
                               plot,
                               width,
                               height,
                               units = "in",
                               res = 300,
                               ...) {
  if (!base::is.character(filename) || base::length(filename) != 1 || !base::nzchar(filename)) {
    stop("`filename` must be a non-empty file path.")
  }
  dots <- base::list(...)
  dots[["device"]] <- NULL
  if (!base::is.null(dots[["dpi"]])) {
    res <- dots[["dpi"]]
    dots[["dpi"]] <- NULL
  }

  pdf_file <- .hc_export_path_with_ext(filename, "pdf")
  png_file <- .hc_export_path_with_ext(filename, "png")

  pdf_args <- base::c(
    base::list(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height,
      units = units,
      device = grDevices::cairo_pdf
    ),
    dots
  )
  base::do.call(ggplot2::ggsave, pdf_args)

  png_args <- base::c(
    base::list(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      units = units,
      dpi = res
    ),
    dots
  )
  if (base::is.null(png_args[["bg"]])) {
    png_args[["bg"]] <- "white"
  }
  tryCatch(
    base::do.call(ggplot2::ggsave, png_args),
    error = function(e) {
      base::warning(
        "Could not write PNG companion ", png_file, ": ",
        base::conditionMessage(e),
        call. = FALSE
      )
    }
  )

  base::invisible(base::list(pdf = pdf_file, png = png_file))
}
