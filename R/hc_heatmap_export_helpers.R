.hc_export_path_with_ext <- function(path, ext) {
  ext <- base::paste0(".", base::sub("^\\.", "", base::as.character(ext[[1]])))
  base::paste0(tools::file_path_sans_ext(path), ext)
}

# Atomic, verified file write.
#
# Sync clients (Sciebo / ownCloud / OneDrive) can grab a file mid-write when the
# output folder is a synced directory, producing a truncated or empty result
# (e.g. an unreadable PDF or a zero-row .xlsx). `producer(tmp)` writes to a
# sibling temporary file in the *same* directory; on success the temp file is
# renamed onto `final_path` (an atomic operation on the same filesystem), so a
# consumer never observes a half-written final file. If the finished file is
# missing or smaller than `min_bytes`, an error is raised so the failure is
# surfaced loudly instead of leaving a silently-broken output behind.
.hc_output_payload_valid <- function(path, min_bytes = 1) {
  if (!base::file.exists(path)) {
    return(FALSE)
  }
  size <- base::file.info(path)$size
  if (base::is.na(size) || size < min_bytes) {
    return(FALSE)
  }

  ext <- base::tolower(tools::file_ext(path))
  if (base::identical(ext, "pdf")) {
    signature <- tryCatch(
      readBin(path, what = "raw", n = 5L),
      error = function(e) raw(0)
    )
    return(base::identical(signature, charToRaw("%PDF-")))
  }
  if (base::identical(ext, "png")) {
    signature <- tryCatch(
      readBin(path, what = "raw", n = 8L),
      error = function(e) raw(0)
    )
    return(base::identical(
      signature,
      as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
    ))
  }
  if (base::identical(ext, "xlsx")) {
    entries <- tryCatch(
      base::suppressWarnings(utils::unzip(path, list = TRUE)$Name),
      error = function(e) character(0)
    )
    return(
      "xl/workbook.xml" %in% entries &&
        base::any(base::grepl("^xl/worksheets/sheet[^/]*\\.xml$", entries))
    )
  }

  TRUE
}

.hc_write_atomic <- function(final_path, producer, min_bytes = 1) {
  final_path <- base::as.character(final_path[[1]])
  if (base::is.na(final_path) || !base::nzchar(final_path)) {
    stop("`final_path` must be a non-empty file path.")
  }
  if (!base::is.function(producer)) {
    stop("`producer` must be a function of one argument (the temp path).")
  }
  if (!base::is.numeric(min_bytes) || base::length(min_bytes) != 1L ||
    base::is.na(min_bytes) || !base::is.finite(min_bytes) || min_bytes < 1) {
    stop("`min_bytes` must be a positive finite number.")
  }

  dir_path <- base::dirname(final_path)
  if (!base::dir.exists(dir_path)) {
    base::dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  # Keep the original extension on the temp file: graphics devices (Cairo),
  # ggplot2::ggsave and openxlsx all infer the output format from the file
  # extension, so the temp name must end in the same `.pdf`/`.png`/`.xlsx`.
  ext <- tools::file_ext(final_path)
  file_ext <- if (base::nzchar(ext)) base::paste0(".", ext) else ""
  tmp_path <- base::tempfile(
    pattern = base::paste0(
      base::basename(tools::file_path_sans_ext(final_path)),
      ".part-"
    ),
    tmpdir = dir_path,
    fileext = file_ext
  )
  backup_path <- NULL
  backup_created <- FALSE
  replacement_started <- FALSE
  committed <- FALSE
  base::on.exit(
    {
      if (base::file.exists(tmp_path)) {
        base::try(base::file.remove(tmp_path), silent = TRUE)
      }
      if (isTRUE(backup_created) && !base::is.null(backup_path) &&
        base::file.exists(backup_path)) {
        if (isTRUE(committed) || !isTRUE(replacement_started)) {
          base::try(base::file.remove(backup_path), silent = TRUE)
        } else {
          if (base::file.exists(final_path)) {
            base::try(base::file.remove(final_path), silent = TRUE)
          }
          restored <- base::suppressWarnings(base::file.rename(backup_path, final_path))
          if (!isTRUE(restored)) {
            restored <- base::file.copy(backup_path, final_path, overwrite = TRUE)
            if (isTRUE(restored)) {
              base::try(base::file.remove(backup_path), silent = TRUE)
            }
          }
          if (!isTRUE(restored)) {
            base::warning(
              "Could not restore the previous output after a failed write. ",
              "The backup remains at '", backup_path, "'.",
              call. = FALSE
            )
          }
        }
      } else if (!isTRUE(committed) && isTRUE(replacement_started) &&
        base::file.exists(final_path)) {
        base::try(base::file.remove(final_path), silent = TRUE)
      }
    },
    add = TRUE
  )

  producer(tmp_path)

  if (!.hc_output_payload_valid(tmp_path, min_bytes = min_bytes)) {
    stop(base::sprintf(
      "Temporary output '%s' was not written or has an invalid payload.", tmp_path
    ))
  }

  if (base::file.exists(final_path)) {
    backup_path <- base::tempfile(
      pattern = base::paste0(
        base::basename(tools::file_path_sans_ext(final_path)),
        ".backup-"
      ),
      tmpdir = dir_path,
      fileext = file_ext
    )
    backup_created <- base::file.copy(final_path, backup_path, overwrite = FALSE)
    if (!isTRUE(backup_created) ||
      !base::file.exists(backup_path) ||
      base::file.info(backup_path)$size != base::file.info(final_path)$size) {
      stop("Could not create a backup of the existing output '", final_path, "'.")
    }
  } else {
    replacement_started <- TRUE
  }

  moved <- base::suppressWarnings(base::file.rename(tmp_path, final_path))
  if (!isTRUE(moved)) {
    if (base::file.exists(final_path)) {
      removed <- base::suppressWarnings(base::file.remove(final_path))
      if (!isTRUE(removed)) {
        stop("Could not replace the existing output '", final_path, "'.")
      }
      replacement_started <- TRUE
    }
    moved <- base::suppressWarnings(base::file.rename(tmp_path, final_path))
    if (!isTRUE(moved)) {
      moved <- base::file.copy(tmp_path, final_path, overwrite = FALSE)
    }
  } else {
    replacement_started <- TRUE
  }
  if (!isTRUE(moved) || !.hc_output_payload_valid(final_path, min_bytes = min_bytes)) {
    stop(base::sprintf(
      "Failed to finalize output '%s' (write to a synced folder such as Sciebo/OneDrive may have been interrupted).",
      final_path
    ))
  }

  committed <- TRUE
  invisible(final_path)
}

# Verify a file that should already exist is present and non-empty; warn (do not
# stop) so a failed export becomes visible in the log instead of silent.
.hc_verify_output_file <- function(path, label = NULL, min_bytes = 1) {
  if (base::is.null(path) || base::length(path) == 0) {
    return(invisible(FALSE))
  }
  path <- base::as.character(path[[1]])
  if (base::is.na(path) || !base::nzchar(path)) {
    return(invisible(FALSE))
  }
  ok <- .hc_output_payload_valid(path, min_bytes = min_bytes)
  if (!ok) {
    base::warning(
      base::sprintf(
        "Expected output %s'%s' is missing, empty, or invalid after writing. If the save folder is on a synced drive (Sciebo/OneDrive), the sync client may have interrupted the write; try a local output folder.",
        if (!base::is.null(label)) base::paste0(label, " ") else "",
        path
      ),
      call. = FALSE
    )
  }
  invisible(ok)
}

# Atomic drop-in for `openxlsx::write.xlsx(x = ..., file = ..., ...)`: keep the
# call site identical except for the function name. Writes the workbook to a
# temporary sibling and atomically moves it onto `file`, so a synced output
# folder (Sciebo/OneDrive) cannot leave a truncated / zero-row .xlsx behind.
.hc_write_xlsx_atomic <- function(x, file, ...) {
  dots <- base::list(...)
  .hc_write_atomic(
    final_path = file,
    producer = function(tmp) {
      base::do.call(
        openxlsx::write.xlsx,
        base::c(base::list(x = x, file = tmp), dots)
      )
    }
  )
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
                                onefile = TRUE,
                                bg = "white") {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    if (!base::is.null(dpi)) {
      Cairo::Cairo(
        file = file,
        width = width,
        height = height,
        pointsize = pointsize,
        dpi = dpi,
        type = "pdf",
        units = "in",
        bg = bg,
        canvas = bg
      )
    } else {
      Cairo::CairoPDF(
        file = file,
        width = width,
        height = height,
        pointsize = pointsize,
        onefile = onefile,
        bg = bg
      )
    }
  } else {
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      pointsize = pointsize,
      onefile = onefile,
      bg = bg
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

.hc_draw_white_page_background <- function() {
  grid::grid.rect(
    x = 0.5,
    y = 0.5,
    width = 1,
    height = 1,
    gp = grid::gpar(fill = "white", col = NA)
  )
  invisible(NULL)
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

  # Render to a temporary sibling file, then atomically move it onto the final
  # path, so a synced output folder (Sciebo/OneDrive) cannot leave a truncated
  # PDF/PNG behind.
  render_page <- function(target, open_device) {
    .hc_write_atomic(target, function(tmp) {
      open_device(tmp)
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      .hc_draw_white_page_background()
      draw_fun()
      invisible(NULL)
    })
  }

  render_page(pdf_file, function(f) {
    .hc_open_pdf_device(
      file = f,
      width = width,
      height = height,
      pointsize = pointsize,
      dpi = pdf_dpi
    )
  })
  render_page(png_file, function(f) {
    .hc_open_png_device(
      file = f,
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

  # Multi-page PDF: render every page to a temporary sibling, then atomically
  # move it onto the final path (protects synced output folders from truncation).
  .hc_write_atomic(pdf_file, function(tmp) {
    .hc_open_pdf_device(
      file = tmp,
      width = width,
      height = height,
      pointsize = pointsize,
      onefile = TRUE
    )
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    for (idx in base::seq_along(page_labels)) {
      .hc_draw_white_page_background()
      draw_page_fun(idx, page_labels[[idx]])
    }
    invisible(NULL)
  })

  stem <- tools::file_path_sans_ext(pdf_file)
  for (idx in base::seq_along(page_labels)) {
    png_file <- base::sprintf(
      "%s_%03d_%s.png",
      stem,
      idx,
      .hc_export_sanitize_stem(page_labels[[idx]], default = base::sprintf("page_%03d", idx))
    )
    .hc_write_atomic(png_file, function(tmp) {
      .hc_open_png_device(
        file = tmp,
        width = width,
        height = height,
        res = res,
        pointsize = pointsize
      )
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      .hc_draw_white_page_background()
      draw_page_fun(idx, page_labels[[idx]])
      invisible(NULL)
    })
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

  .hc_write_atomic(pdf_file, function(tmp) {
    pdf_args <- base::c(
      base::list(
        filename = tmp,
        plot = plot,
        width = width,
        height = height,
        units = units,
        device = grDevices::cairo_pdf
      ),
      dots
    )
    if (base::is.null(pdf_args[["bg"]])) {
      pdf_args[["bg"]] <- "white"
    }
    base::do.call(ggplot2::ggsave, pdf_args)
  })

  tryCatch(
    .hc_write_atomic(png_file, function(tmp) {
      png_args <- base::c(
        base::list(
          filename = tmp,
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
      base::do.call(ggplot2::ggsave, png_args)
    }),
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
