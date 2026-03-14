#' Evaluate an expression with a temporary RNG seed.
#' @noRd
.hc_with_seed <- function(seed, expr) {
  expr <- substitute(expr)
  seed <- suppressWarnings(as.integer(seed[[1]]))
  if (length(seed) == 0L || is.na(seed)) {
    return(eval.parent(expr))
  }

  seed_env <- .GlobalEnv
  had_seed <- exists(".Random.seed", envir = seed_env, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = seed_env, inherits = FALSE) else NULL

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = seed_env)
    } else if (exists(".Random.seed", envir = seed_env, inherits = FALSE)) {
      rm(".Random.seed", envir = seed_env)
    }
  }, add = TRUE)

  base::get("set.seed", mode = "function", envir = baseenv())(seed)
  eval.parent(expr)
}
