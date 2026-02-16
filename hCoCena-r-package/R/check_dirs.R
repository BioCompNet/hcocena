#' Fixes Directories 
#' 
#' Iteratively calls fix_dir() on the provided working directory paths to fix them if necessary or throw an error if they are invalid.
#' @param create_output_dir Boolean. If TRUE and `dir_output` does not exist, it will be created.
#' @export

check_dirs <- function(create_output_dir = FALSE){
  for(x in base::names(hcobject[["working_directory"]])){
    current_dir <- hcobject[["working_directory"]][[x]]
    if (x == "dir_output" && !identical(current_dir, FALSE) && isTRUE(create_output_dir) && !base::dir.exists(current_dir)) {
      base::dir.create(current_dir, recursive = TRUE, showWarnings = FALSE)
      message("Created missing output directory: ", current_dir)
    }
    hcobject[["working_directory"]][[x]] <<- fix_dir(hcobject[["working_directory"]][[x]])
  }
}

