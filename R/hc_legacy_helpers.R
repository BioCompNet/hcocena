# Internal helpers for the S4 transition phase.

hcobject <- NULL

.hc_legacy_state_env <- function() {
  asNamespace("hcocena")
}

.hc_bind_legacy_hcobject <- function(hcobject, envo = .hc_legacy_state_env()) {
  had_existing <- base::exists("hcobject", envir = envo, inherits = FALSE)
  old_hcobject <- NULL
  if (isTRUE(had_existing)) {
    old_hcobject <- base::get("hcobject", envir = envo, inherits = FALSE)
  }

  was_locked <- base::bindingIsLocked("hcobject", env = envo)
  if (isTRUE(was_locked)) {
    base::unlockBinding("hcobject", envo)
  }
  base::assign("hcobject", hcobject, envir = envo)

  list(
    envo = envo,
    had_existing = had_existing,
    old_hcobject = old_hcobject,
    binding_locked = was_locked
  )
}

.hc_restore_legacy_hcobject <- function(state) {
  if (base::is.null(state) || base::is.null(state$envo)) {
    return(invisible(NULL))
  }

  was_locked <- isTRUE(state$binding_locked) &&
    base::bindingIsLocked("hcobject", env = state$envo)
  if (isTRUE(was_locked)) {
    base::unlockBinding("hcobject", state$envo)
  }
  on.exit(
    if (isTRUE(state$binding_locked)) {
      base::lockBinding("hcobject", state$envo)
    },
    add = TRUE
  )

  if (isTRUE(state$had_existing)) {
    base::assign("hcobject", state$old_hcobject, envir = state$envo)
  } else if (base::exists("hcobject", envir = state$envo, inherits = FALSE)) {
    base::rm("hcobject", envir = state$envo)
  }

  invisible(NULL)
}

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

.hc_has_usable_legacy_hcobject <- function(envo = .hc_legacy_state_env()) {
  if (!base::exists("hcobject", envir = envo, inherits = FALSE)) {
    return(FALSE)
  }
  base::is.list(base::get("hcobject", envir = envo, inherits = FALSE))
}

.hc_get_legacy_hcobject <- function(envo = .hc_legacy_state_env()) {
  if (!.hc_has_usable_legacy_hcobject(envo = envo)) {
    return(.hc_default_object())
  }

  hcobject <- base::get("hcobject", envir = envo, inherits = FALSE)

  defaults <- .hc_default_object()
  missing_names <- base::setdiff(base::names(defaults), base::names(hcobject))
  for (nm in missing_names) {
    hcobject[[nm]] <- defaults[[nm]]
  }

  hcobject
}

.hc_assign_legacy_hcobject <- function(hcobject, envo = .hc_legacy_state_env()) {
  was_locked <- base::bindingIsLocked("hcobject", env = envo)
  if (isTRUE(was_locked)) {
    base::unlockBinding("hcobject", envo)
  }
  on.exit(
    if (isTRUE(was_locked)) {
      base::lockBinding("hcobject", envo)
    },
    add = TRUE
  )

  base::assign("hcobject", hcobject, envir = envo)
  invisible(hcobject)
}

.hc_resolve_legacy_bridge_envs <- function(envo = .hc_legacy_state_env(),
                                           mirror_envo = NULL) {
  source_envo <- envo
  if (!.hc_has_usable_legacy_hcobject(envo = source_envo) &&
      .hc_has_usable_legacy_hcobject(envo = .GlobalEnv)) {
    source_envo <- .GlobalEnv
    if (base::is.null(mirror_envo) && !base::identical(envo, .GlobalEnv)) {
      mirror_envo <- envo
    }
  }

  if (base::is.null(mirror_envo) &&
      .hc_has_usable_legacy_hcobject(envo = .GlobalEnv) &&
      !base::identical(source_envo, .GlobalEnv)) {
    mirror_envo <- .GlobalEnv
  }

  list(
    source_envo = source_envo,
    mirror_envo = mirror_envo
  )
}

.hc_normalize_modern_bridge_result <- function(out) {
  if (inherits(out, "HCoCenaExperiment")) {
    return(list(hc = out, result = NULL))
  }

  if (base::is.list(out) &&
      "hc" %in% base::names(out) &&
      inherits(out[["hc"]], "HCoCenaExperiment")) {
    return(list(
      hc = out[["hc"]],
      result = out[["result"]]
    ))
  }

  stop(
    "Modern bridge functions must return a `HCoCenaExperiment` ",
    "or a list with `hc` and optional `result`."
  )
}

.hc_run_modern_from_legacy_capture <- function(fun,
                                               ...,
                                               envo = .hc_legacy_state_env(),
                                               mirror_envo = NULL) {
  bridge_envs <- .hc_resolve_legacy_bridge_envs(
    envo = envo,
    mirror_envo = mirror_envo
  )
  source_envo <- bridge_envs$source_envo
  mirror_envo <- bridge_envs$mirror_envo

  hc <- as_hcocena(.hc_get_legacy_hcobject(envo = source_envo))
  out <- base::do.call(fun, c(list(hc), list(...)))
  normalized <- .hc_normalize_modern_bridge_result(out)
  legacy <- as_hcobject(normalized$hc)
  .hc_assign_legacy_hcobject(legacy, envo = source_envo)

  if (!base::is.null(mirror_envo) && !base::identical(mirror_envo, source_envo)) {
    .hc_assign_legacy_hcobject(legacy, envo = mirror_envo)
  }

  normalized
}

.hc_run_modern_from_legacy <- function(fun,
                                       ...,
                                       envo = .hc_legacy_state_env(),
                                       mirror_envo = NULL) {
  .hc_run_modern_from_legacy_capture(
    fun = fun,
    ...,
    envo = envo,
    mirror_envo = mirror_envo
  )

  invisible(NULL)
}

.hc_run_legacy_via_modern <- function(legacy_name, modern_fun, ...) {
  .hc_legacy_warning(legacy_name)
  invisible(.hc_run_modern_from_legacy(modern_fun, ...))
}

.hc_run_legacy_capture <- function(hc, fun, ..., envo = .hc_legacy_state_env()) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment` object.")
  }

  legacy_state <- .hc_bind_legacy_hcobject(as_hcobject(hc), envo = envo)
  on.exit(.hc_restore_legacy_hcobject(legacy_state), add = TRUE)

  old_opt <- getOption("hcocena.suppress_legacy_warning", FALSE)
  options(hcocena.suppress_legacy_warning = TRUE)
  on.exit(options(hcocena.suppress_legacy_warning = old_opt), add = TRUE)

  fun_name <- NULL
  if (is.character(fun)) {
    fun_name <- fun[[1]]
    if (!base::exists(fun_name, envir = envo, mode = "function", inherits = TRUE)) {
      stop("Legacy function `", fun, "` was not found.")
    }
    fun <- base::get(fun_name, envir = envo, mode = "function", inherits = TRUE)
  }

  dot_args <- list(...)
  formal_names <- base::names(base::formals(fun))
  if (!("..." %in% formal_names)) {
    arg_names <- base::names(dot_args)
    extra_args <- arg_names[base::nzchar(arg_names) & !(arg_names %in% formal_names)]
    if (base::length(extra_args) > 0) {
      stop(
        "Legacy function `", fun_name %||% "<anonymous>", "` does not accept argument(s): ",
        base::paste(extra_args, collapse = ", "),
        ". This usually means an outdated/stale function version is loaded. ",
        "Reload `hcocena` with `devtools::load_all(reset = TRUE)` or restart R.",
        call. = FALSE
      )
    }
  }

  result <- base::do.call(fun, dot_args, envir = envo)

  list(
    hc = as_hcocena(base::get("hcobject", envir = envo, inherits = FALSE)),
    result = result
  )
}

.hc_run_legacy <- function(hc, fun, ..., envo = .hc_legacy_state_env()) {
  .hc_run_legacy_capture(
    hc = hc,
    fun = fun,
    ...,
    envo = envo
  )[["hc"]]
}

.hc_as_hcobject_for_cluster_plot <- function(hc) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment` object.")
  }

  out <- .hc_default_object()
  out[["working_directory"]] <- .hc_row_to_list(hc@config@paths)
  out[["global_settings"]] <- .hc_row_to_list(hc@config@global)
  out[["integrated_output"]][["GFC_all_layers"]] <- base::as.data.frame(hc@integration@gfc)
  out[["integrated_output"]][["cluster_calc"]] <- as.list(hc@integration@cluster)
  out[["satellite_outputs"]] <- as.list(hc@satellite)

  exps <- MultiAssayExperiment::experiments(hc@mae)
  if (base::length(exps) > 0) {
    layer_ids <- base::names(exps)
    if (is.null(layer_ids)) {
      layer_ids <- base::paste0("set", base::seq_along(exps))
    }
    out[["layers"]] <- stats::setNames(
      base::replicate(base::length(layer_ids), base::character(0), simplify = FALSE),
      layer_ids
    )
    for (nm in layer_ids) {
      se <- exps[[nm]]
      out[["data"]][[base::paste0(nm, "_anno")]] <- base::as.data.frame(SummarizedExperiment::colData(se))
    }
  } else if (base::nrow(hc@config@layer) > 0 && "layer_id" %in% base::colnames(hc@config@layer)) {
    layer_ids <- base::as.character(hc@config@layer$layer_id)
    out[["layers"]] <- stats::setNames(
      base::replicate(base::length(layer_ids), base::character(0), simplify = FALSE),
      layer_ids
    )
  }

  out
}

.hc_update_hc_from_cluster_plot <- function(hc, hcobject) {
  if (!inherits(hc, "HCoCenaExperiment")) {
    stop("`hc` must be a `HCoCenaExperiment` object.")
  }

  cluster_calc <- hcobject[["integrated_output"]][["cluster_calc"]] %||% list()
  satellite_outputs <- hcobject[["satellite_outputs"]] %||% list()

  hc@integration@cluster <- S4Vectors::SimpleList(cluster_calc)
  hc@satellite <- S4Vectors::SimpleList(satellite_outputs)
  methods::validObject(hc)
  hc
}
