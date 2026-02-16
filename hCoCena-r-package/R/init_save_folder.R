#' Creates a Save Folder
#' 
#' A folder with the given name is created in the output directory. All analysis outputs will be saved to this folder.
#' @param name The name of the folder to be created.
#'  Use `""` to write outputs directly into `dir_output` without creating a subfolder.
#' @export

init_save_folder <- function(name){
  if (is.null(name) || base::length(name) != 1 || is.na(name)) {
    stop("`name` must be a non-NA character scalar. Use \"\" to skip a subfolder.")
  }

  out_dir <- hcobject[["working_directory"]][["dir_output"]]
  if (is.null(out_dir) || identical(out_dir, FALSE) || !base::nzchar(as.character(out_dir))) {
    stop("`dir_output` is not set. Please run `init_wd()`/`hc_set_paths()` first.")
  }

  if (!base::dir.exists(out_dir)) {
    base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("Created output directory: ", out_dir)
  }

  save_folder <- as.character(name)

  if (save_folder == "") {
    message("Using output directory directly (no additional save subfolder): ", out_dir)
  } else {
    target_dir <- base::file.path(out_dir, save_folder)
    if (!base::dir.exists(target_dir)) {
      base::dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
      message("Created save folder: ", target_dir)
    } else {
      message("Using existing save folder: ", target_dir)
    }
  }

  hcobject[["global_settings"]][["save_folder"]] <<- name
}
