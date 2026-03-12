#' Module-Level Significance Between Conditions
#'
#' Computes module-level statistics on sample-wise module means to identify
#' modules with significant differences between conditions.
#'
#' Methods:
#' - `Wilcoxon` (2 groups; for >2 groups an omnibus Kruskal-Wallis fallback is
#'   used automatically),
#' - `limma` pairwise contrasts,
#' - `LMM` (`nlme::lme`) for repeated measures with donor random intercepts.
#'
#' Results are stored in `hcobject[["satellite_outputs"]][[slot_name]]` and can
#' be added to the main module heatmap with
#' `plot_cluster_heatmap(include_module_significance = TRUE)`.
#'
#' @param set Integer vector of layer indices or `"all"` (default).
#' @param condition_col Column name in annotation files defining conditions.
#'   If `NULL`, uses `hcobject[["global_settings"]][["voi"]]`.
#' @param donor_col Optional donor identifier column in annotation files.
#'   Required if `run_lmm = TRUE`.
#' @param time_col Optional time column in annotation files (used by LMM when
#'   available and `lmm_include_time = TRUE`).
#' @param run_wilcox Logical; run Wilcoxon/Kruskal module-wise test.
#' @param run_limma Logical; run limma pairwise contrasts.
#' @param run_lmm Logical; run linear mixed model (`nlme::lme`) per module.
#' @param lmm_include_time Logical; include `time_col` as fixed effect in LMM
#'   when available.
#' @param limma_reference Optional reference condition for one-vs-reference
#'   limma contrasts. If `NULL`, all pairwise contrasts are tested.
#' @param limma_trend Logical; passed to `limma::eBayes(trend = ...)`.
#' @param padj Multiple-testing correction method (passed to `p.adjust`/limma).
#' @param export_excel Logical; write result tables to Excel in save folder.
#' @param excel_file File name for the Excel export.
#' @param slot_name Satellite slot name for storing results.
#' @return Invisibly returns the stored result list.
#' @export
module_condition_significance <- function(set = "all",
                                          condition_col = NULL,
                                          donor_col = NULL,
                                          time_col = NULL,
                                          run_wilcox = TRUE,
                                          run_limma = TRUE,
                                          run_lmm = FALSE,
                                          lmm_include_time = TRUE,
                                          limma_reference = NULL,
                                          limma_trend = TRUE,
                                          padj = "BH",
                                          export_excel = TRUE,
                                          excel_file = "Module_condition_significance.xlsx",
                                          slot_name = "module_condition_significance") {

  if (isTRUE(run_limma) && !requireNamespace("limma", quietly = TRUE)) {
    warning(
      "Package `limma` is not installed; skipping limma contrasts (`run_limma = FALSE`). ",
      "Install with `BiocManager::install('limma')` to enable this step."
    )
    run_limma <- FALSE
  }
  if (isTRUE(run_lmm) && !requireNamespace("nlme", quietly = TRUE)) {
    warning(
      "Package `nlme` is not installed; skipping mixed-model step (`run_lmm = FALSE`). ",
      "Install with `install.packages('nlme')` to enable this step."
    )
    run_lmm <- FALSE
  }

  if (!isTRUE(run_wilcox) && !isTRUE(run_limma) && !isTRUE(run_lmm)) {
    stop(
      "No analysis step left to run. Enable at least one of ",
      "`run_wilcox`, `run_limma`, or `run_lmm` (and ensure required packages are installed)."
    )
  }
  if (!base::is.null(limma_reference) && base::length(limma_reference) > 1) {
    limma_reference <- base::as.character(limma_reference[[1]])
    warning("More than one value passed to `limma_reference`; using the first value only.")
  }

  if (base::is.null(condition_col) || !base::nzchar(base::as.character(condition_col))) {
    condition_col <- base::as.character(hcobject[["global_settings"]][["voi"]])
  } else {
    condition_col <- base::as.character(condition_col)
  }
  if (!base::nzchar(condition_col)) {
    stop("Could not infer `condition_col`. Set it explicitly.")
  }

  set_indices <- .hc_mc_parse_set_indices(set = set)
  prep <- .hc_mc_collect_module_scores(
    set_indices = set_indices,
    condition_col = condition_col,
    donor_col = donor_col,
    time_col = time_col
  )
  module_mat <- prep$module_mat
  sample_meta <- prep$sample_meta

  out_wilcox <- NULL
  out_wilcox_pairwise <- NULL
  out_limma <- NULL
  out_limma_summary <- NULL
  out_lmm <- NULL

  if (isTRUE(run_wilcox)) {
    wc <- .hc_mc_run_wilcox(module_mat = module_mat, sample_meta = sample_meta, padj = padj)
    out_wilcox <- wc$global
    out_wilcox_pairwise <- wc$pairwise
  }

  if (isTRUE(run_limma)) {
    lm_out <- .hc_mc_run_limma(
      module_mat = module_mat,
      sample_meta = sample_meta,
      padj = padj,
      reference = limma_reference,
      trend = limma_trend
    )
    out_limma <- lm_out$pairwise
    out_limma_summary <- lm_out$summary
  }

  if (isTRUE(run_lmm)) {
    if (base::is.null(donor_col) || !base::nzchar(base::as.character(donor_col))) {
      stop("`donor_col` must be provided when `run_lmm = TRUE`.")
    }
    out_lmm <- .hc_mc_run_lmm(
      module_mat = module_mat,
      sample_meta = sample_meta,
      include_time = lmm_include_time,
      padj = padj
    )
  }

  summary_tbl <- base::data.frame(
    cluster = prep$modules,
    stringsAsFactors = FALSE
  )

  module_label_map <- hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]]
  if (!base::is.null(module_label_map) && base::length(module_label_map) > 0) {
    module_label_map <- base::as.character(module_label_map)
    base::names(module_label_map) <- base::as.character(base::names(hcobject[["integrated_output"]][["cluster_calc"]][["module_label_map"]]))
    summary_tbl$module_label <- base::as.character(module_label_map[summary_tbl$cluster])
  } else {
    summary_tbl$module_label <- summary_tbl$cluster
  }

  if (!base::is.null(out_wilcox) && base::nrow(out_wilcox) > 0) {
    m <- base::match(summary_tbl$cluster, out_wilcox$cluster)
    summary_tbl$wilcox_p <- out_wilcox$p[m]
    summary_tbl$wilcox_q <- out_wilcox$p_adj[m]
    summary_tbl$wilcox_sig <- out_wilcox$significance[m]
  }
  if (!base::is.null(out_limma_summary) && base::nrow(out_limma_summary) > 0) {
    m <- base::match(summary_tbl$cluster, out_limma_summary$cluster)
    summary_tbl$limma_contrast <- out_limma_summary$contrast[m]
    summary_tbl$limma_logFC <- out_limma_summary$logFC[m]
    summary_tbl$limma_p <- out_limma_summary$p[m]
    summary_tbl$limma_q <- out_limma_summary$p_adj[m]
    summary_tbl$limma_sig <- out_limma_summary$significance[m]
  }
  if (!base::is.null(out_lmm) && base::nrow(out_lmm) > 0) {
    m <- base::match(summary_tbl$cluster, out_lmm$cluster)
    summary_tbl$lmm_p <- out_lmm$p[m]
    summary_tbl$lmm_q <- out_lmm$p_adj[m]
    summary_tbl$lmm_sig <- out_lmm$significance[m]
  }

  q_cols <- base::intersect(base::c("wilcox_q", "limma_q", "lmm_q"), base::colnames(summary_tbl))
  if (base::length(q_cols) > 0) {
    q_mat <- base::as.matrix(summary_tbl[, q_cols, drop = FALSE])
    summary_tbl$best_q <- base::apply(q_mat, 1, function(x) {
      x <- suppressWarnings(base::as.numeric(x))
      x <- x[base::is.finite(x)]
      if (base::length(x) == 0) {
        return(NA_real_)
      }
      base::min(x, na.rm = TRUE)
    })
    summary_tbl$best_method <- base::apply(q_mat, 1, function(x) {
      x <- suppressWarnings(base::as.numeric(x))
      if (!base::any(base::is.finite(x))) {
        return(NA_character_)
      }
      q_cols[[base::which.min(base::ifelse(base::is.finite(x), x, Inf))]]
    })
    summary_tbl$best_sig <- .hc_mc_signif_label(summary_tbl$best_q)
  }

  slot_out <- list(
    summary = summary_tbl,
    wilcox = out_wilcox,
    wilcox_pairwise = out_wilcox_pairwise,
    limma_pairwise = out_limma,
    limma_summary = out_limma_summary,
    lmm = out_lmm,
    sample_metadata = sample_meta,
    parameters = list(
      set = set_indices,
      condition_col = condition_col,
      donor_col = donor_col,
      time_col = time_col,
      run_wilcox = run_wilcox,
      run_limma = run_limma,
      run_lmm = run_lmm,
      lmm_include_time = lmm_include_time,
      limma_reference = limma_reference,
      limma_trend = limma_trend,
      padj = padj
    )
  )

  hcobject[["satellite_outputs"]][[slot_name]] <<- slot_out

  if (isTRUE(export_excel)) {
    xlsx_tables <- list(summary = summary_tbl)
    if (!base::is.null(out_wilcox) && base::nrow(out_wilcox) > 0) {
      xlsx_tables[["wilcox_global"]] <- out_wilcox
    }
    if (!base::is.null(out_wilcox_pairwise) && base::nrow(out_wilcox_pairwise) > 0) {
      xlsx_tables[["wilcox_pairwise"]] <- out_wilcox_pairwise
    }
    if (!base::is.null(out_limma) && base::nrow(out_limma) > 0) {
      xlsx_tables[["limma_pairwise"]] <- out_limma
    }
    if (!base::is.null(out_limma_summary) && base::nrow(out_limma_summary) > 0) {
      xlsx_tables[["limma_summary"]] <- out_limma_summary
    }
    if (!base::is.null(out_lmm) && base::nrow(out_lmm) > 0) {
      xlsx_tables[["lmm"]] <- out_lmm
    }

    out_dir <- base::file.path(
      hcobject[["working_directory"]][["dir_output"]],
      hcobject[["global_settings"]][["save_folder"]]
    )
    if (!base::dir.exists(out_dir)) {
      base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    out_file <- base::file.path(out_dir, excel_file)
    tryCatch(
      openxlsx::write.xlsx(x = xlsx_tables, file = out_file, overwrite = TRUE),
      error = function(e) warning("Could not write module-significance Excel: ", conditionMessage(e))
    )
  }

  invisible(slot_out)
}

.hc_mc_signif_label <- function(p) {
  p <- suppressWarnings(base::as.numeric(p))
  out <- base::rep("", base::length(p))
  out[base::is.finite(p) & p <= 0.001] <- "***"
  out[base::is.finite(p) & p > 0.001 & p <= 0.01] <- "**"
  out[base::is.finite(p) & p > 0.01 & p <= 0.05] <- "*"
  out
}

.hc_mc_parse_set_indices <- function(set = "all") {
  n_layers <- base::length(hcobject[["layers_names"]])
  if (n_layers == 0) {
    stop("No layers found in `hcobject`. Run data import first.")
  }
  if (base::is.character(set) && base::length(set) == 1 &&
      base::tolower(base::as.character(set)) == "all") {
    return(base::seq_len(n_layers))
  }
  idx <- suppressWarnings(base::as.integer(set))
  if (base::length(idx) == 0 || base::any(!base::is.finite(idx))) {
    stop("`set` must be \"all\" or integer layer indices.")
  }
  idx <- base::unique(idx)
  if (base::any(idx < 1) || base::any(idx > n_layers)) {
    stop("`set` indices must be in 1..", n_layers, ".")
  }
  idx
}

.hc_mc_collect_module_scores <- function(set_indices,
                                         condition_col,
                                         donor_col = NULL,
                                         time_col = NULL) {
  cluster_info <- hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]]
  if (base::is.null(cluster_info) || base::nrow(cluster_info) == 0) {
    stop("Missing `integrated_output$cluster_calc$cluster_information`.")
  }
  if (!all(base::c("color", "cluster_included") %in% base::colnames(cluster_info))) {
    stop("`cluster_information` is missing required columns (`color`, `cluster_included`).")
  }

  c_df <- dplyr::filter(cluster_info, cluster_included == "yes")
  modules <- base::unique(base::as.character(c_df$color))
  modules <- modules[!base::is.na(modules) & base::nzchar(modules) & modules != "white"]
  if (base::length(modules) == 0) {
    stop("No included non-white modules found.")
  }

  layer_names <- base::as.character(hcobject[["layers_names"]])
  mats <- list()
  metas <- list()

  for (si in set_indices) {
    anno_name <- base::paste0("set", si, "_anno")
    if (!(anno_name %in% base::names(hcobject[["data"]]))) {
      stop("Missing annotation table `", anno_name, "`.")
    }
    anno <- hcobject[["data"]][[anno_name]]
    if (base::is.null(base::rownames(anno)) || base::nrow(anno) == 0) {
      stop("Annotation table `", anno_name, "` has no row names / samples.")
    }
    if (!(condition_col %in% base::colnames(anno))) {
      stop("`condition_col` = '", condition_col, "' is missing in `", anno_name, "`.")
    }
    if (!base::is.null(donor_col) && !(donor_col %in% base::colnames(anno))) {
      stop("`donor_col` = '", donor_col, "' is missing in `", anno_name, "`.")
    }
    if (!base::is.null(time_col) && !(time_col %in% base::colnames(anno))) {
      stop("`time_col` = '", time_col, "' is missing in `", anno_name, "`.")
    }

    m <- sample_wise_cluster_expression(set = si)
    if (base::is.null(m) || base::length(m) == 0) {
      next
    }
    m <- base::as.matrix(m)
    if (base::is.null(base::rownames(m)) || base::is.null(base::colnames(m))) {
      next
    }
    samples <- base::intersect(base::colnames(m), base::rownames(anno))
    if (base::length(samples) == 0) {
      next
    }

    m <- m[, samples, drop = FALSE]
    sample_uid <- base::paste0("set", si, "::", samples)

    full <- base::matrix(
      NA_real_,
      nrow = base::length(modules),
      ncol = base::length(samples),
      dimnames = list(modules, sample_uid)
    )
    row_common <- base::intersect(base::rownames(m), modules)
    if (base::length(row_common) > 0) {
      full[row_common, ] <- base::as.matrix(m[row_common, samples, drop = FALSE])
    }

    donor_vals <- base::rep(NA_character_, base::length(samples))
    if (!base::is.null(donor_col)) {
      donor_vals <- base::as.character(anno[samples, donor_col])
    }

    time_vals <- base::rep(NA_character_, base::length(samples))
    if (!base::is.null(time_col)) {
      time_vals <- base::as.character(anno[samples, time_col])
    }

    meta <- base::data.frame(
      sample_uid = sample_uid,
      sample = samples,
      set = base::paste0("set", si),
      layer = layer_names[[si]],
      condition = base::as.character(anno[samples, condition_col]),
      donor = donor_vals,
      time = time_vals,
      stringsAsFactors = FALSE
    )

    mats[[base::length(mats) + 1]] <- full
    metas[[base::length(metas) + 1]] <- meta
  }

  if (base::length(mats) == 0) {
    stop("No usable samples found for selected `set`.")
  }

  module_mat <- base::do.call(base::cbind, mats)
  sample_meta <- base::do.call(base::rbind, metas)
  sample_meta <- sample_meta[base::match(base::colnames(module_mat), sample_meta$sample_uid), , drop = FALSE]
  base::rownames(sample_meta) <- sample_meta$sample_uid

  list(
    module_mat = module_mat,
    sample_meta = sample_meta,
    modules = modules
  )
}

.hc_mc_run_wilcox <- function(module_mat, sample_meta, padj = "BH") {
  keep <- !base::is.na(sample_meta$condition) & base::nzchar(base::as.character(sample_meta$condition))
  if (base::sum(keep) < 2) {
    stop("Not enough samples with non-missing conditions for Wilcoxon test.")
  }

  cond <- base::factor(base::as.character(sample_meta$condition[keep]))
  mat <- module_mat[, keep, drop = FALSE]
  n_groups <- base::nlevels(cond)
  if (n_groups < 2) {
    stop("At least two condition groups are required.")
  }

  global_test <- if (n_groups == 2) "wilcox" else "kruskal_fallback"
  global_rows <- base::lapply(base::rownames(mat), function(cl) {
    x <- suppressWarnings(base::as.numeric(mat[cl, ]))
    k <- base::is.finite(x) & !base::is.na(cond)
    x <- x[k]
    g <- base::droplevels(cond[k])

    if (base::length(x) < 2 || base::nlevels(g) < 2) {
      return(base::data.frame(
        cluster = cl,
        test = global_test,
        comparison = NA_character_,
        n = base::length(x),
        group_sizes = NA_character_,
        effect = NA_real_,
        p = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    tab <- base::table(g)
    group_sizes <- base::paste(base::paste0(base::names(tab), "=", base::as.integer(tab)), collapse = "; ")
    med <- tapply(x, g, stats::median, na.rm = TRUE)
    effect <- if (base::length(med) >= 2) {
      suppressWarnings(base::as.numeric(base::max(med, na.rm = TRUE) - base::min(med, na.rm = TRUE)))
    } else {
      NA_real_
    }

    if (base::nlevels(g) == 2) {
      lev <- base::levels(g)
      p <- tryCatch(
        stats::wilcox.test(x[g == lev[[1]]], x[g == lev[[2]]], exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
      cmp <- base::paste0(lev[[1]], " vs ", lev[[2]])
      used_test <- "wilcox"
    } else {
      p <- tryCatch(stats::kruskal.test(x ~ g)$p.value, error = function(e) NA_real_)
      cmp <- base::paste(base::levels(g), collapse = " | ")
      used_test <- global_test
    }

    base::data.frame(
      cluster = cl,
      test = used_test,
      comparison = cmp,
      n = base::length(x),
      group_sizes = group_sizes,
      effect = effect,
      p = suppressWarnings(base::as.numeric(p)),
      stringsAsFactors = FALSE
    )
  })
  global <- base::do.call(base::rbind, global_rows)
  global$p_adj <- stats::p.adjust(global$p, method = padj)
  global$significance <- .hc_mc_signif_label(global$p_adj)

  pairwise_rows <- list()
  for (cl in base::rownames(mat)) {
    x <- suppressWarnings(base::as.numeric(mat[cl, ]))
    k <- base::is.finite(x) & !base::is.na(cond)
    x <- x[k]
    g <- base::droplevels(cond[k])
    if (base::length(x) < 2 || base::nlevels(g) < 2) {
      next
    }
    pw <- tryCatch(
      stats::pairwise.wilcox.test(x = x, g = g, p.adjust.method = padj, exact = FALSE),
      error = function(e) NULL
    )
    if (base::is.null(pw) || base::is.null(pw$p.value)) {
      next
    }
    pv <- pw$p.value
    if (base::is.null(base::dim(pv))) {
      next
    }
    for (ri in base::seq_len(base::nrow(pv))) {
      for (ci in base::seq_len(base::ncol(pv))) {
        qv <- suppressWarnings(base::as.numeric(pv[ri, ci]))
        if (!base::is.finite(qv)) {
          next
        }
        g2 <- base::rownames(pv)[ri]
        g1 <- base::colnames(pv)[ci]
        pairwise_rows[[base::length(pairwise_rows) + 1]] <- base::data.frame(
          cluster = cl,
          group1 = g1,
          group2 = g2,
          contrast = base::paste0(g1, " vs ", g2),
          p_adj = qv,
          significance = .hc_mc_signif_label(qv),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  pairwise <- if (base::length(pairwise_rows) > 0) {
    base::do.call(base::rbind, pairwise_rows)
  } else {
    base::data.frame(
      cluster = character(0),
      group1 = character(0),
      group2 = character(0),
      contrast = character(0),
      p_adj = numeric(0),
      significance = character(0),
      stringsAsFactors = FALSE
    )
  }

  list(global = global, pairwise = pairwise)
}

.hc_mc_run_limma <- function(module_mat,
                             sample_meta,
                             padj = "BH",
                             reference = NULL,
                             trend = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package `limma` is required for `run_limma = TRUE`.")
  }
  keep <- !base::is.na(sample_meta$condition) & base::nzchar(base::as.character(sample_meta$condition))
  if (base::sum(keep) < 2) {
    stop("Not enough samples with non-missing conditions for limma.")
  }

  mat <- module_mat[, keep, drop = FALSE]
  cond_raw <- base::as.character(sample_meta$condition[keep])
  cond_levels <- base::unique(cond_raw)
  if (base::length(cond_levels) < 2) {
    stop("At least two condition groups are required for limma.")
  }

  safe_levels <- base::make.unique(base::make.names(cond_levels))
  safe_map <- stats::setNames(safe_levels, cond_levels)
  safe_to_orig <- stats::setNames(cond_levels, safe_levels)
  cond_safe <- base::factor(safe_map[cond_raw], levels = safe_levels)

  design <- stats::model.matrix(~0 + cond_safe)
  base::colnames(design) <- safe_levels

  if (base::is.null(reference)) {
    combos <- utils::combn(safe_levels, 2, simplify = FALSE)
    contrast_defs <- base::vapply(combos, function(v) base::paste0(v[[1]], " - ", v[[2]]), FUN.VALUE = base::character(1))
    contrast_labels <- base::vapply(combos, function(v) {
      base::paste0(safe_to_orig[[v[[1]]]], " vs ", safe_to_orig[[v[[2]]]])
    }, FUN.VALUE = base::character(1))
  } else {
    ref <- base::as.character(reference[[1]])
    if (!(ref %in% cond_levels)) {
      stop("`limma_reference` ('", ref, "') is not present in `condition_col`.")
    }
    ref_safe <- safe_map[[ref]]
    others <- safe_levels[safe_levels != ref_safe]
    if (base::length(others) == 0) {
      stop("No contrasts can be formed for limma with the chosen `limma_reference`.")
    }
    contrast_defs <- base::vapply(others, function(x) base::paste0(x, " - ", ref_safe), FUN.VALUE = base::character(1))
    contrast_labels <- base::vapply(others, function(x) base::paste0(safe_to_orig[[x]], " vs ", ref), FUN.VALUE = base::character(1))
  }

  contr <- limma::makeContrasts(contrasts = contrast_defs, levels = design)
  if (base::is.null(base::dim(contr))) {
    contr <- base::matrix(contr, ncol = 1, dimnames = list(base::rownames(design), base::names(contrast_defs)))
  }

  fit <- limma::lmFit(mat, design = design)
  fit <- limma::contrasts.fit(fit, contrasts = contr)
  fit <- limma::eBayes(fit, trend = isTRUE(trend))

  pairwise_rows <- base::lapply(base::seq_len(base::ncol(contr)), function(i) {
    tt <- limma::topTable(
      fit,
      coef = i,
      number = Inf,
      sort.by = "none",
      adjust.method = padj
    )
    if (base::nrow(tt) == 0) {
      return(NULL)
    }
    tt$cluster <- base::rownames(tt)
    tt$contrast <- contrast_labels[[i]]
    tt$contrast_id <- base::colnames(contr)[i]
    tt
  })
  pairwise_rows <- pairwise_rows[!base::vapply(pairwise_rows, base::is.null, FUN.VALUE = logical(1))]

  pairwise <- if (base::length(pairwise_rows) > 0) {
    out <- base::do.call(base::rbind, pairwise_rows)
    out$significance <- .hc_mc_signif_label(out$adj.P.Val)
    out
  } else {
    base::data.frame(
      cluster = character(0),
      contrast = character(0),
      contrast_id = character(0),
      logFC = numeric(0),
      AveExpr = numeric(0),
      t = numeric(0),
      P.Value = numeric(0),
      adj.P.Val = numeric(0),
      B = numeric(0),
      significance = character(0),
      stringsAsFactors = FALSE
    )
  }

  summary <- if (base::nrow(pairwise) > 0) {
    split_pw <- base::split(pairwise, pairwise$cluster)
    summary_rows <- base::lapply(base::names(split_pw), function(cl) {
      sub <- split_pw[[cl]]
      qv <- suppressWarnings(base::as.numeric(sub$adj.P.Val))
      if (base::any(base::is.finite(qv))) {
        idx <- base::which.min(base::ifelse(base::is.finite(qv), qv, Inf))
      } else {
        idx <- 1
      }
      base::data.frame(
        cluster = cl,
        contrast = base::as.character(sub$contrast[idx]),
        logFC = suppressWarnings(base::as.numeric(sub$logFC[idx])),
        p = suppressWarnings(base::as.numeric(sub$P.Value[idx])),
        p_adj = suppressWarnings(base::as.numeric(sub$adj.P.Val[idx])),
        t = suppressWarnings(base::as.numeric(sub$t[idx])),
        B = suppressWarnings(base::as.numeric(sub$B[idx])),
        stringsAsFactors = FALSE
      )
    })
    sm <- base::do.call(base::rbind, summary_rows)
    sm$significance <- .hc_mc_signif_label(sm$p_adj)
    sm
  } else {
    base::data.frame(
      cluster = character(0),
      contrast = character(0),
      logFC = numeric(0),
      p = numeric(0),
      p_adj = numeric(0),
      t = numeric(0),
      B = numeric(0),
      significance = character(0),
      stringsAsFactors = FALSE
    )
  }

  list(pairwise = pairwise, summary = summary)
}

.hc_mc_run_lmm <- function(module_mat,
                           sample_meta,
                           include_time = TRUE,
                           padj = "BH") {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("Package `nlme` is required for `run_lmm = TRUE`.")
  }
  keep <- !base::is.na(sample_meta$condition) &
    base::nzchar(base::as.character(sample_meta$condition)) &
    !base::is.na(sample_meta$donor) &
    base::nzchar(base::as.character(sample_meta$donor))
  if (base::sum(keep) < 3) {
    stop("Not enough samples with non-missing condition and donor for LMM.")
  }

  mat <- module_mat[, keep, drop = FALSE]
  cond_vec <- base::as.character(sample_meta$condition[keep])
  donor_vec <- base::as.character(sample_meta$donor[keep])
  time_vec <- base::as.character(sample_meta$time[keep])

  rows <- base::lapply(base::rownames(mat), function(cl) {
    x <- suppressWarnings(base::as.numeric(mat[cl, ]))
    df <- base::data.frame(
      value = x,
      cond = cond_vec,
      donor = donor_vec,
      time = time_vec,
      stringsAsFactors = FALSE
    )
    df <- df[base::is.finite(df$value) &
               !base::is.na(df$cond) & base::nzchar(df$cond) &
               !base::is.na(df$donor) & base::nzchar(df$donor), , drop = FALSE]

    if (base::nrow(df) < 3) {
      return(base::data.frame(
        cluster = cl,
        test = "lmm",
        n = base::nrow(df),
        n_condition = base::length(base::unique(df$cond)),
        n_donor = base::length(base::unique(df$donor)),
        effect = NA_real_,
        p = NA_real_,
        status = "insufficient_data",
        stringsAsFactors = FALSE
      ))
    }

    df$cond <- base::factor(df$cond)
    df$donor <- base::factor(df$donor)
    if (base::nlevels(df$cond) < 2 || base::nlevels(df$donor) < 2) {
      return(base::data.frame(
        cluster = cl,
        test = "lmm",
        n = base::nrow(df),
        n_condition = base::nlevels(df$cond),
        n_donor = base::nlevels(df$donor),
        effect = NA_real_,
        p = NA_real_,
        status = "insufficient_levels",
        stringsAsFactors = FALSE
      ))
    }

    use_time <- FALSE
    if (isTRUE(include_time) && "time" %in% base::colnames(df)) {
      has_time <- !base::is.na(df$time) & base::nzchar(base::as.character(df$time))
      use_time <- base::any(has_time)
      if (use_time) {
        if (base::all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", base::as.character(df$time[has_time])))) {
          df$time <- suppressWarnings(base::as.numeric(df$time))
        } else {
          df$time <- base::factor(df$time, levels = base::unique(df$time))
        }
      }
    }

    fixed_formula <- if (isTRUE(use_time)) {
      stats::as.formula("value ~ cond + time")
    } else {
      stats::as.formula("value ~ cond")
    }

    fit <- tryCatch(
      nlme::lme(
        fixed = fixed_formula,
        random = ~1 | donor,
        data = df,
        method = "REML",
        na.action = stats::na.omit,
        control = nlme::lmeControl(msMaxIter = 200, returnObject = TRUE)
      ),
      error = function(e) NULL
    )
    if (base::is.null(fit)) {
      return(base::data.frame(
        cluster = cl,
        test = if (isTRUE(use_time)) "lmm_cond_plus_time" else "lmm_cond",
        n = base::nrow(df),
        n_condition = base::nlevels(df$cond),
        n_donor = base::nlevels(df$donor),
        effect = NA_real_,
        p = NA_real_,
        status = "fit_failed",
        stringsAsFactors = FALSE
      ))
    }

    a <- tryCatch(nlme::anova.lme(fit), error = function(e) NULL)
    p <- NA_real_
    if (!base::is.null(a)) {
      p_col <- base::colnames(a)[base::grepl("p", base::tolower(base::colnames(a)))]
      if (base::length(p_col) > 0 && "cond" %in% base::rownames(a)) {
        p <- suppressWarnings(base::as.numeric(a["cond", p_col[[1]]]))
      }
    }

    m <- tapply(df$value, df$cond, base::mean, na.rm = TRUE)
    effect <- if (base::length(m) >= 2) {
      suppressWarnings(base::as.numeric(base::max(m, na.rm = TRUE) - base::min(m, na.rm = TRUE)))
    } else {
      NA_real_
    }

    base::data.frame(
      cluster = cl,
      test = if (isTRUE(use_time)) "lmm_cond_plus_time" else "lmm_cond",
      n = base::nrow(df),
      n_condition = base::nlevels(df$cond),
      n_donor = base::nlevels(df$donor),
      effect = effect,
      p = p,
      status = if (base::is.finite(p)) "ok" else "no_pvalue",
      stringsAsFactors = FALSE
    )
  })

  out <- base::do.call(base::rbind, rows)
  out$p_adj <- stats::p.adjust(out$p, method = padj)
  out$significance <- .hc_mc_signif_label(out$p_adj)
  out
}
