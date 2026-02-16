#' Fetch reference data from ExperimentHub
#'
#' This is the transition API for the planned `hcocenaData` companion package.
#' The function expects an ExperimentHub ID (e.g. `"EH12345"`).
#'
#' @param hub_id Character scalar with an ExperimentHub record identifier.
#' @param hub Optional pre-initialized `ExperimentHub::ExperimentHub` object.
#' @return The resource stored under `hub_id`.
#' @export
hc_get_reference_data <- function(hub_id, hub = NULL) {
  if (!is.character(hub_id) || length(hub_id) != 1 || !nzchar(hub_id)) {
    stop("`hub_id` must be a non-empty character scalar.")
  }

  if (is.null(hub)) {
    if (!requireNamespace("ExperimentHub", quietly = TRUE)) {
      stop(
        "Package `ExperimentHub` is required. Install it with ",
        "`BiocManager::install('ExperimentHub')`."
      )
    }
    hub <- ExperimentHub::ExperimentHub()
  }

  hub[[hub_id]]
}

