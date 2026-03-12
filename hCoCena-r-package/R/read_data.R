
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

	base::options(dplyr.summarise.inform = F)

	normalize_layer_source <- function(src, kind, layer_idx){
		if(!base::is.character(src) || base::length(src) != 1 || base::is.na(src) || !base::nzchar(src)){
			stop(
				"Layer ", layer_idx, ": ", kind, " source must be a single non-empty string.\n",
				"If you provide objects, pass quoted names in `hc_define_layers()`, e.g. ",
				"`c(\"counts_df\", \"anno_df\")`."
			)
		}
		base::as.character(src)
	}

	legacy_dir_is_false <- function(x){
		base::isFALSE(x) || (base::is.character(x) && base::identical(x, "FALSE"))
	}

	resolve_layer_file <- function(dir_path, source_name){
		if(legacy_dir_is_false(dir_path)){
			return(source_name)
		}
		base::paste0(dir_path, source_name)
	}

	normalize_count_object <- function(obj, source_name){
		if(base::is.matrix(obj)){
			obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
		}
		if(!base::is.data.frame(obj)){
			stop("Count source `", source_name, "` exists but is not a data.frame or matrix.")
		}
		obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)

		if(!base::is.null(gene_symbol_col) && gene_symbol_col %in% base::colnames(obj)){
			obj <- make_rownames_unique(counts = obj, gene_symbol_col = gene_symbol_col)
		}else{
			if(base::is.null(base::rownames(obj)) || base::any(!base::nzchar(base::rownames(obj)))){
				stop(
					"Count object `", source_name, "` has no usable rownames.\n",
					"Provide `gene_symbol_col` (matching a column in the count object) ",
					"or set rownames to gene symbols before calling `hc_read_data()`."
				)
			}
			numeric_cols <- base::vapply(obj, base::is.numeric, logical(1))
			if(!base::all(numeric_cols)){
				obj <- obj[, numeric_cols, drop = FALSE]
			}
			if(base::ncol(obj) == 0){
				stop("Count object `", source_name, "` has no numeric sample columns.")
			}
		}

		obj
	}

	normalize_anno_object <- function(obj, source_name){
		if(base::is.matrix(obj)){
			obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
		}
		if(!base::is.data.frame(obj)){
			stop("Annotation source `", source_name, "` exists but is not a data.frame or matrix.")
		}
		obj <- base::as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)

		if(!base::is.null(sample_col) && sample_col %in% base::colnames(obj)){
			base::rownames(obj) <- base::as.character(obj[[sample_col]])
		}
		if(base::is.null(base::rownames(obj)) || base::any(!base::nzchar(base::rownames(obj)))){
			stop(
				"Annotation object `", source_name, "` has no usable rownames.\n",
				"Provide `sample_col` (matching a column in the annotation object) ",
				"or set rownames to sample IDs before calling `hc_read_data()`."
			)
		}

		obj[] <- base::lapply(obj, base::factor)
		obj
	}

	data <- list()

	# iterate over layers:
	for (x in 1:base::length(hcobject[["layers"]])){
		count_source <- normalize_layer_source(
			hcobject[["layers"]][[base::paste0("set",x)]][1],
			"count",
			x
		)
		anno_source <- normalize_layer_source(
			hcobject[["layers"]][[base::paste0("set",x)]][2],
			"annotation",
			x
		)

		# read expression data:

		# check if provided string refers to an existing object, if so load the object instead of reading from file:
		if(base::exists(count_source)){
		  
		  data[[base::paste0("set", x, "_counts")]] <- normalize_count_object(
		  	base::get(count_source),
		  	count_source
		  )
		  
		}else{
		  
		  if(base::is.null(gene_symbol_col)){
		  	stop("You must provide the 'gene_symbol_col' parameter.")
		  }
		  count_file <- resolve_layer_file(
		  	hcobject[["working_directory"]][["dir_count_data"]],
		  	count_source
		  )
		  if(!base::file.exists(count_file)){
		  	stop(
		  		"Count source `", count_source, "` was not found as an object, and file `",
		  		count_file, "` does not exist.\n",
		  		"If you use object input, pass quoted object names in `hc_define_layers()`."
		  	)
		  }
		  data[[base::paste0("set", x, "_counts")]] <- read_expression_data(file = count_file, rown = count_has_rn,
		                                                              sep = sep_counts, gene_symbol_col = gene_symbol_col)
		}
	  # check for zero-variance genes
	    var.df <- rank_variance(data[[base::paste0("set", x, "_counts")]])
	    if(base::all(base::is.na(var.df$variance))){
	    	stop(
	    		"Count data for layer ", x, " contains no valid numeric expression values after preprocessing."
	    	)
	    }
	    if(any(var.df$variance == 0, na.rm = TRUE)){
	    message(base::paste0("Detected genes with 0 variance in dataset ", x, "."))
	    data[[base::paste0("set", x, "_counts_unfiltered")]] <- data[[base::paste0("set", x, "_counts")]]
	    data[[base::paste0("set", x, "_counts")]] <- data[[base::paste0("set", x, "_counts")]] %>%
	      dplyr::filter(!(row.names(.) %in% (var.df %>% dplyr::filter(variance == 0) %>% dplyr::pull(gene))))
	    message(base::paste0((nrow(data[[base::paste0("set", x, "_counts_unfiltered")]])-nrow(data[[base::paste0("set", x, "_counts")]])), " gene(s) were removed from dataset ", x, "."))
	    } 
	  

		# read annotation data:

		# check if provided string refers to an existing object, if so load the object instead of reading from file:
		if(base::exists(anno_source)){
		  data[[base::paste0("set", x, "_anno")]] <- normalize_anno_object(
		  	base::get(anno_source),
		  	anno_source
		  )
		  
		}else{
		  if(base::is.null(sample_col)){
		  	stop("You must provide the 'sample_col' parameter.")
		  }
		  anno_file <- resolve_layer_file(
		  	hcobject[["working_directory"]][["dir_annotation"]],
		  	anno_source
		  )
		  if(!base::file.exists(anno_file)){
		  	stop(
		  		"Annotation source `", anno_source, "` was not found as an object, and file `",
		  		anno_file, "` does not exist.\n",
		  		"If you use object input, pass quoted object names in `hc_define_layers()`."
		  	)
		  }
		  data[[base::paste0("set", x, "_anno")]] <- read_anno(file = anno_file, rown = anno_has_rn,
		                                                 sep = sep_anno, sample_col = sample_col)
		  
		}

		# check if samples match between annotation and counts:
		if(!base::ncol(data[[base::paste0("set", x, "_counts")]]) == base::nrow(data[[base::paste0("set", x, "_anno")]])){
		  stop(base::paste0("The count table has ",  base::ncol(data[[base::paste0("set", x, "_counts")]]), " columns but the annotation has ", 
		              base::nrow(data[[base::paste0("set", x, "_anno")]]), " rows. These values are required to be the same since they
		                 should correspond to the number of samples. THE LOADING OF THE DATA WILL BE TERMINATED."))
		}else if(!base::all(base::as.character(base::colnames(data[[base::paste0("set", x, "_counts")]]))  %in% base::as.character(base::rownames(data[[base::paste0("set", x, "_anno")]])))){
		  stop(base::paste0("The column names of the count file do not all match the rownames of the annotation. Please make sure they contain the same samples."))
		}else{
		  data[[base::paste0("set", x, "_anno")]] <- data[[base::paste0("set", x, "_anno")]][base::colnames(data[[base::paste0("set", x, "_counts")]]),]
		}

	}

	hcobject[["data"]] <<- data
}
