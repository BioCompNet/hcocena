#' Antilog
#'
#' Reverses a log operation.
#' @param lx The log value.
#' @param base The base to which the log was performed

.hc_antilog_impl <- function(lx, base) {
  lbx <- lx / base::log(base::exp(1), base = base)
  result <- base::exp(lbx)
  return(result)
}

#' Reverse a log transformation
#'
#' Inverse of [base::log()] for an arbitrary base, i.e. `hc_antilog(x, 2)`
#' returns `2^x`.
#'
#' @param lx Numeric vector, matrix or data frame of logged values.
#' @param base The base the logarithm was taken to.
#' @return The de-logged values, in the same shape as `lx`.
#' @export
hc_antilog <- function(lx, base) {
  .hc_antilog_impl(lx = lx, base = base)
}
