#' Fixes Directories 
#' 
#' Iteratively calls fix_dir() on the provided working directory paths to fix them if necessary or throw an error if they are invalid.
#' @param create_output_dir Boolean. If TRUE and `dir_output` does not exist, it will be created.
#' @export

check_dirs <- function(create_output_dir = FALSE){
  invisible(.hc_run_modern_from_legacy(
    .hc_check_dirs_impl,
    create_output_dir = create_output_dir
  ))
}

