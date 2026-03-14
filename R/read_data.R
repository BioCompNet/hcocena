
#' Function To Read All Count And Annotation Files
#' 
#' The function loads the count and annotation data for each layer and saves it in the hCoCena-Object's "data" slot.
#' Genes with a variance of 0 will be automatically be excluded from the analysis.
#' @param sep_counts The separator of the count files. Default is tab separated files. Ignore when loading data from objects instead of files.
#' @param sep_anno The separator of the annotation files. Default is tab separated files. Ignore when loading data from objects instead of files.
#' @param gene_symbol_col A String. Name of the column that contains the gene symbols. Ignore when loading data from objects instead of files.
#' @param sample_col A String. Name of the column that contains the sample IDs. Ignore when loading data from objects instead of files.
#' @param count_has_rn A Boolean. Whether or not the count file has rownames. Default is TRUE. Ignore when loading data from objects instead of files.
#' @param anno_has_rn A Boolean. Whether or not the annotation file has rownames. Default is TRUE. Ignore when loading data from objects instead of files.
#' @export

read_data <- function(sep_counts = "\t", 
						sep_anno = "\t", 
						gene_symbol_col = NULL, 
						sample_col = NULL, 
						count_has_rn = TRUE, 
						anno_has_rn = TRUE){
	.hc_legacy_warning("read_data")

	invisible(.hc_run_modern_from_legacy(
		.hc_read_data_impl,
		sep_counts = sep_counts,
		sep_anno = sep_anno,
		gene_symbol_col = gene_symbol_col,
		sample_col = sample_col,
		count_has_rn = count_has_rn,
		anno_has_rn = anno_has_rn,
		auto_setup_output = FALSE,
		source_env = parent.frame()
	))
}
