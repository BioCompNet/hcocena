#' Read Supplementary Data
#' 
#' Function to read and collect the supplementary files set in set_supp_files().
#' @export

read_supplementary <- function (){
  .hc_legacy_warning("read_supplementary")
  invisible(.hc_run_modern_from_legacy(
    .hc_read_supplementary_impl
  ))
}
