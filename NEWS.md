# hcocena 0.99.4

## Fixes

- Added an explicit early error and documentation note that
  `impute_method = "rfcont"` requires `library(CALIBERrfimpute)` in the
  current session.

# hcocena 0.99.3

## Fixes

- Registered `mice.impute.rfcont` directly in the local longitudinal
  imputation call frame in addition to temporary global exposure for more
  robust `mice` lookup in interactive sessions.

# hcocena 0.99.2

## Fixes

- Hardened `rfcont` longitudinal imputation registration for interactive
  sessions by exposing `mice.impute.rfcont` temporarily during `mice` calls.

# hcocena 0.99.1

## Fixes

- Fixed `rfcont` longitudinal imputation so `CALIBERrfimpute` methods are
  available to `mice` during package use.

# hcocena 0.99.0

## Bioconductor preparation

- Aligned the package version with Bioconductor pre-submission conventions.
- Removed redundant `Author` and `Maintainer` fields from `DESCRIPTION` in favor
  of `Authors@R`.
- Added a package-level `README.md` and `inst/CITATION`.
- Expanded the workflow and migration vignettes to use `BiocStyle` and
  reproducible examples based on `inst/extdata`.
- Added ignore rules for local build artifacts and updated the submission export
  helper.
