
#' Plot of Cutoff Statistics
#' 
#' Plots the R-squared value, number of edges, number of genes and number of networks for different cut-offs.
#' @param interactive Boolean. If TRUE (default) plot is created using plotly including and interactive cutoff slider. If FALSE, plot is created as a static ggplot.
#' @param hline A list with four slots ("R.squared", "no_edges", "no_nodes", "no_networks") each of which can be set either to NULL (default) or a number to introduce a horizontal line for orientation at that value in the respective plot. Only used in the non-interactive plot.
#' @export

plot_cutoffs <- function(interactive = T, 
                         hline = list("R.squared" = NULL, 
                                      "no_edges" = NULL, 
                                      "no_nodes" = NULL, 
                                      "no_networks" = NULL)){

  .hc_get_cutoff_tuning_marker <- function(layer_id, layer_index) {
    .as_num1 <- function(x) {
      y <- suppressWarnings(base::as.numeric(x))
      if (base::length(y) < 1) {
        return(NA_real_)
      }
      y[[1]]
    }
    .as_chr1 <- function(x) {
      y <- base::as.character(x)
      if (base::length(y) < 1) {
        return(NA_character_)
      }
      y[[1]]
    }
    .normalize_policy <- function(x) {
      pol <- base::tolower(.as_chr1(x))
      if (!base::is.character(pol) || base::is.na(pol) || !base::nzchar(pol)) {
        return(NA_character_)
      }
      base::switch(
        pol,
        tier1 = "strict",
        tier2 = "relaxed",
        fallback = "best_available",
        pol
      )
    }
    .is_missing_chr <- function(x) {
      !base::is.character(x) || base::is.na(x) || !base::nzchar(x)
    }
    .same_cutoff <- function(a, b, tol = 1e-8) {
      base::is.finite(a) && base::is.finite(b) && base::abs(a - b) <= tol
    }
    .extract_from_summary <- function(summary_df, flag_col) {
      if (!base::is.data.frame(summary_df) || !(flag_col %in% base::colnames(summary_df))) {
        return(NA_real_)
      }
      if (!("cutoff" %in% base::colnames(summary_df))) {
        return(NA_real_)
      }
      cuts <- suppressWarnings(base::as.numeric(summary_df$cutoff))
      flags <- base::as.logical(summary_df[[flag_col]])
      idx <- base::which(!base::is.na(flags) & flags)
      if (base::length(idx) < 1) {
        return(NA_real_)
      }
      idx <- idx[base::is.finite(cuts[idx])]
      if (base::length(idx) < 1) {
        return(NA_real_)
      }
      cuts[[idx[[1]]]]
    }

    sat <- hcobject[["satellite_outputs"]]
    if (base::is.null(sat) || !base::is.list(sat)) {
      return(list(
        selected_cutoff = NA_real_,
        strict_cutoff = NA_real_,
        relaxed_cutoff = NA_real_,
        best_available_cutoff = NA_real_,
        simple_cutoff = NA_real_,
        selected_policy = NA_character_,
        tiered_cutoff = NA_real_,
        tier = NA_character_
      ))
    }
    tuning <- sat[["cutoff_tuning"]]
    if (base::is.null(tuning) || !base::is.list(tuning)) {
      tuning <- sat[["auto_tune_cutoff_tiered"]]
    }
    if (base::is.null(tuning) || !base::is.list(tuning)) {
      return(list(
        selected_cutoff = NA_real_,
        strict_cutoff = NA_real_,
        relaxed_cutoff = NA_real_,
        best_available_cutoff = NA_real_,
        simple_cutoff = NA_real_,
        selected_policy = NA_character_,
        tiered_cutoff = NA_real_,
        tier = NA_character_
      ))
    }

    selected_cutoff <- NA_real_
    strict_cutoff <- NA_real_
    relaxed_cutoff <- NA_real_
    best_available_cutoff <- NA_real_
    simple_cutoff <- NA_real_
    selected_policy <- NA_character_

    rec <- tuning[["recommended_cutoff_vector"]]
    if (!base::is.null(rec)) {
      rec_num <- suppressWarnings(base::as.numeric(rec))
      rec_nm <- base::names(rec)
      if (!base::is.null(rec_nm) && layer_id %in% rec_nm) {
        selected_cutoff <- .as_num1(rec[[layer_id]])
      } else if (base::length(rec_num) >= layer_index) {
        selected_cutoff <- rec_num[[layer_index]]
      }
    }

    cmp <- tuning[["comparison_with_simple"]]
    if (base::is.data.frame(cmp) && base::nrow(cmp) > 0) {
      row_idx <- NULL
      if ("layer_id" %in% base::colnames(cmp)) {
        hit <- base::which(base::as.character(cmp$layer_id) == layer_id)
        if (base::length(hit) > 0) {
          row_idx <- hit[[1]]
        }
      }
      if (base::is.null(row_idx) && base::nrow(cmp) >= layer_index) {
        row_idx <- layer_index
      }
      if (!base::is.null(row_idx)) {
        if ("cutoff_simple" %in% base::colnames(cmp)) {
          simple_cutoff <- .as_num1(cmp$cutoff_simple[[row_idx]])
        }
        if ("selected_policy" %in% base::colnames(cmp)) {
          selected_policy <- .normalize_policy(cmp$selected_policy[[row_idx]])
        }
        if (("tier" %in% base::colnames(cmp)) &&
            (!base::is.character(selected_policy) || base::is.na(selected_policy) || !base::nzchar(selected_policy))) {
          selected_policy <- .normalize_policy(cmp$tier[[row_idx]])
        }
        if (!base::is.finite(selected_cutoff) && "cutoff_tiered" %in% base::colnames(cmp)) {
          selected_cutoff <- .as_num1(cmp$cutoff_tiered[[row_idx]])
        }
        if ("cutoff_strict" %in% base::colnames(cmp)) {
          strict_cutoff <- .as_num1(cmp$cutoff_strict[[row_idx]])
        } else if ("strict_cutoff" %in% base::colnames(cmp)) {
          strict_cutoff <- .as_num1(cmp$strict_cutoff[[row_idx]])
        }
        if ("cutoff_relaxed" %in% base::colnames(cmp)) {
          relaxed_cutoff <- .as_num1(cmp$cutoff_relaxed[[row_idx]])
        } else if ("relaxed_cutoff" %in% base::colnames(cmp)) {
          relaxed_cutoff <- .as_num1(cmp$relaxed_cutoff[[row_idx]])
        }
        if ("cutoff_best_available" %in% base::colnames(cmp)) {
          best_available_cutoff <- .as_num1(cmp$cutoff_best_available[[row_idx]])
        } else if ("best_available_cutoff" %in% base::colnames(cmp)) {
          best_available_cutoff <- .as_num1(cmp$best_available_cutoff[[row_idx]])
        }
      }
    }

    ld_all <- tuning[["layer_details"]]
    ld <- NULL
    if (base::is.list(ld_all) && base::length(ld_all) > 0) {
      if (layer_id %in% base::names(ld_all)) {
        ld <- ld_all[[layer_id]]
      } else if (base::length(ld_all) >= layer_index) {
        ld <- ld_all[[layer_index]]
      }
    }
    if (base::is.list(ld)) {
      if (!base::is.finite(selected_cutoff) && "selected_cutoff" %in% base::names(ld)) {
        selected_cutoff <- .as_num1(ld$selected_cutoff)
      }
      if (.is_missing_chr(selected_policy) && "selected_policy" %in% base::names(ld)) {
        selected_policy <- .normalize_policy(ld$selected_policy)
      }
      if (.is_missing_chr(selected_policy) && "selected_tier" %in% base::names(ld)) {
        selected_policy <- .normalize_policy(ld$selected_tier)
      }
      if (!base::is.finite(strict_cutoff) && "strict_cutoff" %in% base::names(ld)) {
        strict_cutoff <- .as_num1(ld$strict_cutoff)
      }
      if (!base::is.finite(relaxed_cutoff) && "relaxed_cutoff" %in% base::names(ld)) {
        relaxed_cutoff <- .as_num1(ld$relaxed_cutoff)
      }
      if (!base::is.finite(best_available_cutoff) && "best_available_cutoff" %in% base::names(ld)) {
        best_available_cutoff <- .as_num1(ld$best_available_cutoff)
      }

      summary_df <- ld[["cutoff_summary"]]
      if (base::is.data.frame(summary_df) && base::nrow(summary_df) > 0) {
        if (!base::is.finite(strict_cutoff)) {
          strict_cutoff <- .extract_from_summary(summary_df, "strict_ok")
        }
        if (!base::is.finite(strict_cutoff)) {
          strict_cutoff <- .extract_from_summary(summary_df, "tier1_ok")
        }
        if (!base::is.finite(relaxed_cutoff)) {
          relaxed_cutoff <- .extract_from_summary(summary_df, "relaxed_ok")
        }
        if (!base::is.finite(relaxed_cutoff)) {
          relaxed_cutoff <- .extract_from_summary(summary_df, "tier2_ok")
        }
      }
    }

    if (!base::is.finite(selected_cutoff)) {
      if (base::is.finite(strict_cutoff)) {
        selected_cutoff <- strict_cutoff
      } else if (base::is.finite(relaxed_cutoff)) {
        selected_cutoff <- relaxed_cutoff
      } else if (base::is.finite(best_available_cutoff)) {
        selected_cutoff <- best_available_cutoff
      }
    }
    if (.is_missing_chr(selected_policy) && base::is.finite(selected_cutoff)) {
      if (.same_cutoff(selected_cutoff, strict_cutoff)) {
        selected_policy <- "strict"
      } else if (.same_cutoff(selected_cutoff, relaxed_cutoff)) {
        selected_policy <- "relaxed"
      } else if (.same_cutoff(selected_cutoff, best_available_cutoff)) {
        selected_policy <- "best_available"
      }
    }

    list(
      selected_cutoff = selected_cutoff,
      strict_cutoff = strict_cutoff,
      relaxed_cutoff = relaxed_cutoff,
      best_available_cutoff = best_available_cutoff,
      simple_cutoff = simple_cutoff,
      selected_policy = selected_policy,
      # compatibility aliases for older plotting internals
      tiered_cutoff = selected_cutoff,
      tier = selected_policy
    )
  }

	
	for(x in 1:base::length(hcobject[["layers"]])){
		cutoff_df <- hcobject[["layer_specific_outputs"]][[base::paste0("set",x)]][["part1"]][["cutoff_calc_out"]][["cutoff_stats_concise"]]
    tuning_marker <- .hc_get_cutoff_tuning_marker(
      layer_id = base::paste0("set", x),
      layer_index = x
    )
		if(interactive == T){
		  p <- plot_cutoffs_internal_interactive(cutoff_stats = cutoff_df,
		                                         x = x,
		                                         tuning_marker = tuning_marker)
		}else{
		  p <- plot_cutoffs_internal_static(cutoff_stats = cutoff_df, 
		                                    hline = hline, 
		                                    x = x,
		                                    tuning_marker = tuning_marker)
		}

	print(p)
	hcobject[["layer_specific_outputs"]][[base::paste0("set", x)]][["cutoff_plot"]] <<- p

	}
}	
