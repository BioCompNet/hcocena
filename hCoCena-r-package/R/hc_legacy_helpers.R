# Internal helpers for the S4 transition phase.

.hc_legacy_warning <- function(fun) {
  if (isTRUE(getOption("hcocena.suppress_legacy_warning", FALSE))) {
    return(invisible(NULL))
  }

  warning(
    paste0(
      "`", fun, "()` uses the legacy global `hcobject` API. ",
      "Prefer the S4 workflow (`hc_*()` functions) for new analyses."
    ),
    call. = FALSE
  )
}

.hc_run_legacy <- function(hc, fun, ..., envo = .GlobalEnv) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment` object.")
  }

  base::assign("hcobject", as_hcobject(hc), envir = envo)

  old_opt <- getOption("hcocena.suppress_legacy_warning", FALSE)
  options(hcocena.suppress_legacy_warning = TRUE)
  on.exit(options(hcocena.suppress_legacy_warning = old_opt), add = TRUE)

  if (is.character(fun)) {
    if (!base::exists(fun, mode = "function", inherits = TRUE)) {
      stop("Legacy function `", fun, "` was not found.")
    }
    fun <- base::get(fun, mode = "function", inherits = TRUE)
  }

  base::do.call(fun, list(...))

  as_hcocena(base::get("hcobject", envir = envo, inherits = FALSE))
}

