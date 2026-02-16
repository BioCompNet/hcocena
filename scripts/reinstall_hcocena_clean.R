#!/usr/bin/env Rscript

# Clean reinstall helper for local hcocena development installs on Windows.
# Usage:
#   Rscript scripts/reinstall_hcocena_clean.R [pkg_path] [lib_path]
# Examples:
#   Rscript scripts/reinstall_hcocena_clean.R
#   Rscript scripts/reinstall_hcocena_clean.R "E:/.../hCoCena-r-package"
#   Rscript scripts/reinstall_hcocena_clean.R "E:/.../hCoCena-r-package" "C:/R-lib-4.5"

args <- commandArgs(trailingOnly = TRUE)

pkg_path <- if (length(args) >= 1) args[[1]] else "hCoCena-r-package"
pkg_path <- normalizePath(pkg_path, winslash = "/", mustWork = FALSE)

custom_lib <- if (length(args) >= 2) args[[2]] else NULL
if (!is.null(custom_lib)) {
  custom_lib <- normalizePath(custom_lib, winslash = "/", mustWork = FALSE)
  dir.create(custom_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(custom_lib, .libPaths()))
}

if (!dir.exists(pkg_path)) {
  stop("Package path does not exist: ", pkg_path)
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Unload namespace if currently attached/loaded.
if ("package:hcocena" %in% search()) {
  detach("package:hcocena", unload = TRUE, character.only = TRUE)
}
if ("hcocena" %in% loadedNamespaces()) {
  unloadNamespace("hcocena")
}

lib <- .libPaths()[1]
pkg_dir <- file.path(lib, "hcocena")
lock_dirs <- Sys.glob(file.path(lib, "00LOCK*"))

if (dir.exists(pkg_dir)) {
  unlink(pkg_dir, recursive = TRUE, force = TRUE)
}
if (length(lock_dirs) > 0) {
  unlink(lock_dirs, recursive = TRUE, force = TRUE)
}

# Staged install can occasionally leave corrupted payloads on some Windows setups.
Sys.setenv(R_INSTALL_STAGED = "FALSE")

remotes::install_local(
  pkg_path,
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build_vignettes = FALSE
)

pkg_name <- "hcocena"
if (!requireNamespace(pkg_name, quietly = TRUE)) {
  stop("Package `", pkg_name, "` could not be loaded after installation.")
}
library(pkg_name, character.only = TRUE)
options(hcocena.suppress_legacy_warning = TRUE)
init_object()

message("hcocena clean reinstall completed successfully.")
message("Package source: ", pkg_path)
message("Library path used: ", .libPaths()[1])
