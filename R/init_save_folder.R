#' Creates a Save Folder
#' 
#' A folder with the given name is created in the output directory. All analysis outputs will be saved to this folder.
#' @param name The name of the folder to be created.
#'  Use `""` to write outputs directly into `dir_output` without creating a subfolder.
#' @export

init_save_folder <- function(name){
  invisible(.hc_run_modern_from_legacy(
    .hc_init_save_folder_impl,
    name = name
  ))
}
