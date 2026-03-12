#' Defines the Data Layers
#' 
#' The function collects descriptive names of the datasets as well as the names of the count and annotation files for each dataset.
#' @param data_sets A named list. Each list entry represents a dataset and must have a name and a value. The name should describe the dataset and will be used for plot labels etc., while the entry should be a vector or exactly two strings.
#' 	The first string should be the name of the dataset's count file (incl. file ending), the second the name of the dataset's annotation file (incl. file ending).
#' 	Alternatively to providing the count and annotation data as files, they can also be provided as existing dataframe objects. 
#' 	Give the object names as strings instead of the file names. The count object must be a data frame with sample names as columns and gene names as rows. 
#' 	There must be no additional columns other than those representing the counts per sample. 
#' 	The annotation object must also be a data frame, where row names are sample names that match the column names of the count object, and column names are information categories.
#' @export
#' @examples
#' \dontrun{
#' init_object()
#' define_layers(
#'   data_sets = list(
#'     rhinovirusSet = c("rhinovirusSetCount.txt", "rhinovirusSetAnno.txt"),
#'     influenzaSet = c("influenzaSetCount.txt", "influenzaSetAnno.txt")
#'   )
#' )
#'
#' # or from data frames:
#' define_layers(
#'   data_sets = list(
#'     rhinovirusSet = c("rhinovirusSetCountDf", "rhinovirusSetAnnoDf"),
#'     influenzaSet = c("influenzaSetCountDf", "influenzaSetAnnoDf")
#'   )
#' )
#' }

define_layers <- function(data_sets = list()){
	.hc_legacy_warning("define_layers")

	hcobject[["layers"]] <<- list()
	for(setnum in 1:base::length(data_sets)){
		pair <- data_sets[[setnum]]
		set_name <- base::names(data_sets)[[setnum]]
		if(base::is.null(set_name) || !base::nzchar(set_name)) {
			set_name <- base::paste0("set", setnum)
		}
		if(base::length(pair) < 2){
			stop(
				"Layer `", set_name, "` must provide exactly two entries: ",
				"count source and annotation source."
			)
		}

		count_src <- pair[[1]]
		anno_src <- pair[[2]]
		if(!base::is.character(count_src) || base::length(count_src) != 1 || base::is.na(count_src) || !base::nzchar(count_src)){
			stop(
				"Layer `", set_name, "`: count source must be a single string.\n",
				"If you use objects from the environment, pass quoted object names, e.g. ",
				"`c(\"counts_df\", \"anno_df\")` (not `c(counts_df, anno_df)`)."
			)
		}
		if(!base::is.character(anno_src) || base::length(anno_src) != 1 || base::is.na(anno_src) || !base::nzchar(anno_src)){
			stop(
				"Layer `", set_name, "`: annotation source must be a single string.\n",
				"If you use objects from the environment, pass quoted object names, e.g. ",
				"`c(\"counts_df\", \"anno_df\")` (not `c(counts_df, anno_df)`)."
			)
		}

		hcobject[["layers"]][[base::paste0("set", setnum)]] <<- base::c(base::as.character(count_src), base::as.character(anno_src))
	}

	hcobject[["layers_names"]] <<- base::names(data_sets)
	
}
