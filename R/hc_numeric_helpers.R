.hc_as_numeric_safely <- function(x) {
  if (base::is.null(x)) {
    return(base::numeric())
  }

  if (base::is.data.frame(x)) {
    x <- base::unlist(x, use.names = FALSE)
  } else if (base::is.list(x) && !base::is.object(x)) {
    x <- base::unlist(x, use.names = FALSE)
  }

  nm <- base::names(x)

  if (base::is.numeric(x) || base::is.integer(x)) {
    out <- base::as.numeric(x)
  } else if (base::is.logical(x)) {
    out <- base::as.numeric(x)
  } else {
    x_chr <- base::trimws(base::as.character(x))
    out <- base::rep(NA_real_, base::length(x_chr))
    numeric_like <- !base::is.na(x_chr) &
      base::grepl("^[-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?$", x_chr)
    special_like <- !base::is.na(x_chr) & x_chr %in% base::c("Inf", "+Inf", "-Inf", "NaN")
    if (base::any(numeric_like)) {
      out[numeric_like] <- base::as.numeric(x_chr[numeric_like])
    }
    if (base::any(special_like)) {
      out[special_like] <- base::as.numeric(x_chr[special_like])
    }
  }

  if (!base::is.null(nm) && base::length(nm) == base::length(out)) {
    base::names(out) <- nm
  }
  out
}

.hc_as_integer_safely <- function(x) {
  num <- .hc_as_numeric_safely(x)
  out <- base::rep(NA_integer_, base::length(num))
  ok <- base::is.finite(num) &
    num >= -2147483648 &
    num <= 2147483647 &
    base::abs(num - base::round(num)) < 1e-8
  if (base::any(ok)) {
    out[ok] <- base::as.integer(base::round(num[ok]))
  }
  if (!base::is.null(base::names(num))) {
    base::names(out) <- base::names(num)
  }
  out
}

.hc_first_numeric_value <- function(x, default = NA_real_) {
  out <- .hc_as_numeric_safely(x)
  if (base::length(out) < 1) {
    return(default)
  }
  out[[1]]
}

.hc_max_finite <- function(x) {
  num <- .hc_as_numeric_safely(x)
  num <- num[base::is.finite(num)]
  if (base::length(num) == 0) {
    return(NA_real_)
  }
  base::max(num)
}

.hc_min_finite <- function(x) {
  num <- .hc_as_numeric_safely(x)
  num <- num[base::is.finite(num)]
  if (base::length(num) == 0) {
    return(NA_real_)
  }
  base::min(num)
}

.hc_range_span <- function(x) {
  max_val <- .hc_max_finite(x)
  min_val <- .hc_min_finite(x)
  if (!base::is.finite(max_val) || !base::is.finite(min_val)) {
    return(NA_real_)
  }
  max_val - min_val
}

.hc_log2_safely <- function(x) {
  num <- .hc_as_numeric_safely(x)
  out <- base::rep(NA_real_, base::length(num))
  ok <- base::is.finite(num) & num > 0
  if (base::any(ok)) {
    out[ok] <- base::log2(num[ok])
  }
  if (!base::is.null(base::names(num))) {
    base::names(out) <- base::names(num)
  }
  out
}
