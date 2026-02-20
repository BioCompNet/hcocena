#' decoupleR-based module cell-type activity from Enrichr marker libraries
#'
#' Computes module-wise activity scores for Enrichr-derived cell-type marker sets
#' using `decoupleR`, stores selected/significant results, and writes stacked-bar
#' annotation slots (`enriched_per_cluster_*`) that can be rendered via
#' `plot_cluster_heatmap(include_dynamic_enrichment_slots = TRUE)`.
#'
#' @param databases Character vector of Enrichr library names.
#' @param custom_gmt_files Optional custom GMT file(s) to include as additional
#' marker resources. Accepts a character vector or named list of file paths.
#' @param clusters Either `"all"` (default) or a character vector of module IDs/colors.
#' @param mode Either `"coarse"` (default) or `"fine"`.
#' @param top Number of selected marker activities per module.
#' @param qval Adjusted p-value cutoff.
#' @param padj Multiple-testing correction method.
#' @param activity_input One of `"gfc"` (default), `"fc"`, `"expression"`.
#' @param fc_comparisons Vector like `c("A_vs_B", "C_vs_B")` when `activity_input = "fc"`.
#' @param method decoupleR method suffix (e.g. `"ulm"`).
#' @param minsize Minimum target size passed to decoupleR.
#' @param min_term_genes Minimum genes per Enrichr term.
#' @param annotation_slot Legacy slot override for single-database runs.
#' @param slot_suffix Optional slot suffix; default `"decoupler"`.
#' @param clear_previous_slots If TRUE, clears old `enriched_per_cluster*` slots first.
#' @param coarse_map Optional named regex map used in coarse mode.
#' @param coarse_include_other Keep `"Other"` class in coarse mode.
#' @param refresh_db Refresh Enrichr cache.
#' @param export_excel Write summary workbook.
#' @param excel_file Excel filename.
#' @param plot_heatmap Replot cluster heatmap with dynamic slots.
#' @param heatmap_file_name Heatmap filename when `plot_heatmap = TRUE`.
#' @param ... Passed to `plot_cluster_heatmap()`.
#' @return Invisibly returns a list with activity tables and slot mapping.
#' @export
celltype_activity_decoupler <- function(
    databases = c("Descartes_Cell_Types_and_Tissue_2021", "Human_Gene_Atlas"),
    custom_gmt_files = NULL,
    clusters = c("all"),
    mode = c("coarse", "fine"),
    top = 3,
    qval = 0.1,
    padj = "BH",
    activity_input = c("gfc", "fc", "expression"),
    fc_comparisons = NULL,
    method = "ulm",
    minsize = 5,
    min_term_genes = 5,
    annotation_slot = c("auto", "enriched_per_cluster", "enriched_per_cluster2"),
    slot_suffix = "decoupler",
    clear_previous_slots = FALSE,
    coarse_map = NULL,
    coarse_include_other = TRUE,
    refresh_db = FALSE,
    export_excel = TRUE,
    excel_file = "Module_Celltype_Activity_decoupler.xlsx",
    plot_heatmap = FALSE,
    heatmap_file_name = "Heatmap_modules_celltype_activity_decoupler.pdf",
    ...) {
  .hc_legacy_warning("celltype_activity_decoupler")

  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package `decoupleR` is required for `celltype_activity_decoupler()`.")
  }
  mode <- base::match.arg(mode)
  annotation_slot <- base::match.arg(annotation_slot)
  activity_input <- base::match.arg(base::tolower(base::as.character(activity_input[[1]])), c("gfc", "fc", "expression"))

  if (base::is.null(databases)) databases <- base::character(0)
  databases <- base::unique(base::trimws(base::as.character(databases)))
  databases <- databases[!base::is.na(databases) & databases != ""]
  if (!base::is.numeric(top) || base::length(top) != 1 || base::is.na(top) || top < 1) stop("`top` must be a positive integer.")
  top <- base::as.integer(top)
  if (!base::is.numeric(qval) || base::length(qval) != 1 || base::is.na(qval) || qval <= 0 || qval > 1) stop("`qval` must be in (0,1].")
  if (!base::is.numeric(minsize) || base::length(minsize) != 1 || base::is.na(minsize) || minsize < 1) stop("`minsize` must be >= 1.")
  minsize <- base::as.integer(minsize)
  if (!base::is.numeric(min_term_genes) || base::length(min_term_genes) != 1 || base::is.na(min_term_genes) || min_term_genes < 1) stop("`min_term_genes` must be >= 1.")
  min_term_genes <- base::as.integer(min_term_genes)
  if (!base::is.logical(coarse_include_other) || base::length(coarse_include_other) != 1 || base::is.na(coarse_include_other)) stop("`coarse_include_other` must be TRUE/FALSE.")
  if (!base::is.logical(refresh_db) || base::length(refresh_db) != 1 || base::is.na(refresh_db)) stop("`refresh_db` must be TRUE/FALSE.")
  if (!base::is.logical(export_excel) || base::length(export_excel) != 1 || base::is.na(export_excel)) stop("`export_excel` must be TRUE/FALSE.")
  if (!base::is.logical(plot_heatmap) || base::length(plot_heatmap) != 1 || base::is.na(plot_heatmap)) stop("`plot_heatmap` must be TRUE/FALSE.")
  if (!base::is.logical(clear_previous_slots) || base::length(clear_previous_slots) != 1 || base::is.na(clear_previous_slots)) stop("`clear_previous_slots` must be TRUE/FALSE.")
  if (!base::is.null(slot_suffix)) {
    slot_suffix <- .hc_sanitize_db_token(slot_suffix)
    if (!base::nzchar(slot_suffix)) stop("`slot_suffix` is invalid.")
  }

  if (mode == "coarse") {
    if (base::is.null(coarse_map)) coarse_map <- .hc_default_coarse_celltype_map()
    if (!base::is.character(coarse_map) || base::is.null(base::names(coarse_map)) || any(base::names(coarse_map) == "")) {
      stop("`coarse_map` must be a named character vector.")
    }
  }

  cluster_info <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]][["cluster_information"]], error = function(e) NULL)
  if (base::is.null(cluster_info) || base::nrow(cluster_info) == 0) stop("No cluster information found.")
  if (!all(c("color", "gene_n") %in% base::colnames(cluster_info))) stop("`cluster_information` is missing columns.")
  all_clusters <- base::unique(base::as.character(cluster_info$color))
  all_clusters <- all_clusters[!base::is.na(all_clusters) & all_clusters != "white"]
  if (base::length(all_clusters) == 0) stop("No non-white clusters found.")
  if (base::length(clusters) == 1 && identical(clusters[[1]], "all")) {
    clusters <- all_clusters
  } else {
    clusters <- base::as.character(clusters)
    clusters <- clusters[clusters %in% all_clusters]
    if (base::length(clusters) == 0) stop("No valid clusters selected.")
  }

  cluster_calc <- tryCatch(hcobject[["integrated_output"]][["cluster_calc"]], error = function(e) NULL)
  module_prefix <- "M"
  if (!base::is.null(cluster_calc) && "module_prefix" %in% base::names(cluster_calc)) {
    tmp <- base::as.character(cluster_calc[["module_prefix"]])
    if (base::length(tmp) > 0 && !base::is.na(tmp[[1]]) && base::nzchar(tmp[[1]])) module_prefix <- tmp[[1]]
  }
  module_label_map <- if (!base::is.null(cluster_calc) && "module_label_map" %in% base::names(cluster_calc)) cluster_calc[["module_label_map"]] else NULL
  if (!base::is.null(module_label_map) && base::length(module_label_map) > 0) {
    nm <- base::names(module_label_map)
    module_label_map <- base::as.character(module_label_map)
    if (!base::is.null(nm) && base::length(nm) == base::length(module_label_map)) base::names(module_label_map) <- base::as.character(nm)
  }
  if (base::is.null(module_label_map) || base::length(module_label_map) == 0) module_label_map <- stats::setNames(base::paste0(module_prefix, base::seq_along(all_clusters)), all_clusters)
  missing_map <- base::setdiff(all_clusters, base::names(module_label_map))
  if (base::length(missing_map) > 0) {
    module_label_map <- base::c(module_label_map, stats::setNames(base::paste0(module_prefix, base::seq.int(base::length(module_label_map) + 1, length.out = base::length(missing_map))), missing_map))
  }

  collapse_dup_rows <- function(mat) {
    if (!base::anyDuplicated(base::rownames(mat))) return(mat)
    idx <- base::split(base::seq_len(base::nrow(mat)), base::rownames(mat))
    out <- base::t(base::vapply(idx, function(i) {
      v <- base::colMeans(mat[i, , drop = FALSE], na.rm = TRUE)
      v[base::is.nan(v)] <- NA_real_
      v
    }, FUN.VALUE = base::numeric(base::ncol(mat))))
    out %>% base::as.matrix()
  }
  collapse_dup_cols <- function(mat) {
    if (!base::anyDuplicated(base::colnames(mat))) return(mat)
    idx <- base::split(base::seq_len(base::ncol(mat)), base::colnames(mat))
    out <- base::vapply(idx, function(i) {
      if (base::length(i) == 1) return(mat[, i])
      v <- base::rowMeans(mat[, i, drop = FALSE], na.rm = TRUE)
      v[base::is.nan(v)] <- NA_real_
      v
    }, FUN.VALUE = base::numeric(base::nrow(mat)))
    out <- out %>% base::as.matrix()
    base::rownames(out) <- base::rownames(mat)
    out
  }
  resolve_group_labels <- function(anno_df) {
    voi <- hcobject[["global_settings"]][["voi"]]
    voi <- base::intersect(voi, base::colnames(anno_df))
    if (base::length(voi) == 0) return(base::as.character(anno_df[[1]]))
    vals <- anno_df[, voi, drop = FALSE]
    if (base::ncol(vals) == 1) return(base::as.character(vals[[1]]))
    base::apply(vals, 1, function(x) base::paste(x, collapse = "-"))
  }
  get_net_genes <- function() {
    merged_net <- hcobject[["integrated_output"]][["merged_net"]]
    if (base::is.null(merged_net)) stop("Missing integrated network.")
    g <- igraph::get.vertex.attribute(merged_net, "name")
    g <- base::unique(base::toupper(base::trimws(base::as.character(g))))
    g[!base::is.na(g) & g != ""]
  }
  prepare_group_mean_matrix <- function(set_name) {
    expr <- hcobject[["layer_specific_outputs"]][[set_name]][["part1"]][["topvar"]]
    if (base::is.null(expr)) expr <- hcobject[["data"]][[base::paste0(set_name, "_counts")]]
    anno <- hcobject[["data"]][[base::paste0(set_name, "_anno")]]
    if (base::is.null(expr) || base::is.null(anno)) return(NULL)
    expr <- expr %>% base::as.matrix()
    if (base::nrow(expr) == 0 || base::ncol(expr) == 0) return(NULL)
    mode(expr) <- "numeric"
    samples <- base::intersect(base::colnames(expr), base::rownames(anno))
    if (base::length(samples) == 0) return(NULL)
    expr <- expr[, samples, drop = FALSE]
    anno <- anno[samples, , drop = FALSE]
    grp <- resolve_group_labels(anno)
    bad <- base::is.na(grp) | grp == ""
    if (base::any(bad)) grp[bad] <- samples[bad]
    if (isTRUE(hcobject[["global_settings"]][["data_in_log"]])) expr <- antilog(expr, 2)
    lv <- base::unique(grp)
    mat <- base::vapply(lv, function(v) {
      i <- base::which(grp == v)
      if (base::length(i) == 1) return(expr[, i])
      base::rowMeans(expr[, i, drop = FALSE], na.rm = TRUE)
    }, FUN.VALUE = base::numeric(base::nrow(expr)))
    mat <- mat %>% base::as.matrix()
    base::colnames(mat) <- lv
    base::rownames(mat) <- base::toupper(base::trimws(base::rownames(expr)))
    mat <- mat[!base::is.na(base::rownames(mat)) & base::rownames(mat) != "", , drop = FALSE]
    collapse_dup_rows(mat)
  }
  parse_fc <- function(fc) {
    if (base::is.null(fc) || base::length(fc) == 0) stop("`fc_comparisons` must be provided when `activity_input = 'fc'`.")
    out <- base::lapply(fc, function(x) {
      x <- base::trimws(base::as.character(x))
      p <- base::strsplit(x, "\\s*_vs_\\s*", perl = TRUE)[[1]]
      if (base::length(p) != 2) p <- base::strsplit(x, "\\s+vs\\s+", perl = TRUE)[[1]]
      if (base::length(p) != 2) stop("Invalid `fc_comparisons`: ", x)
      base::data.frame(comparison = base::paste0(base::trimws(p[[1]]), "_vs_", base::trimws(p[[2]])), numerator = base::trimws(p[[1]]), denominator = base::trimws(p[[2]]), stringsAsFactors = FALSE)
    })
    out <- base::do.call(base::rbind, out)
    out <- out[!duplicated(out$comparison), , drop = FALSE]
    base::rownames(out) <- NULL
    out
  }

  fc_summary <- NULL
  if (identical(activity_input, "gfc")) {
    gfc_df <- hcobject[["integrated_output"]][["GFC_all_layers"]]
    if (base::is.null(gfc_df) || base::nrow(gfc_df) == 0 || base::ncol(gfc_df) < 2) stop("Missing `GFC_all_layers`.")
    gene_col <- if ("Gene" %in% base::colnames(gfc_df)) "Gene" else base::colnames(gfc_df)[base::ncol(gfc_df)]
    val_cols <- base::setdiff(base::colnames(gfc_df), gene_col)
    activity_mat <- gfc_df[, val_cols, drop = FALSE] %>% base::as.matrix()
    mode(activity_mat) <- "numeric"
    base::rownames(activity_mat) <- base::toupper(base::trimws(base::as.character(gfc_df[[gene_col]])))
    activity_mat <- activity_mat[!base::is.na(base::rownames(activity_mat)) & base::rownames(activity_mat) != "", , drop = FALSE]
    activity_mat <- collapse_dup_rows(activity_mat)
  } else if (identical(activity_input, "expression")) {
    net_genes <- get_net_genes()
    set_mats <- list()
    for (z in base::seq_len(base::length(hcobject[["layer_specific_outputs"]]))) {
      set_name <- base::paste0("set", z)
      sm <- prepare_group_mean_matrix(set_name)
      if (base::is.null(sm)) next
      full <- base::matrix(NA_real_, nrow = base::length(net_genes), ncol = base::ncol(sm), dimnames = list(net_genes, base::colnames(sm)))
      ov <- base::intersect(net_genes, base::rownames(sm))
      if (base::length(ov) > 0) full[ov, ] <- sm[ov, , drop = FALSE]
      set_mats[[base::length(set_mats) + 1]] <- full
    }
    if (base::length(set_mats) == 0) stop("Could not build expression activity matrix.")
    activity_mat <- base::do.call(base::cbind, set_mats)
    activity_mat <- collapse_dup_rows(activity_mat)
    activity_mat <- collapse_dup_cols(activity_mat)
  } else {
    net_genes <- get_net_genes()
    fc_def <- parse_fc(fc_comparisons)
    comparison_found <- stats::setNames(base::rep(FALSE, base::nrow(fc_def)), base::as.character(fc_def$comparison))
    fc_limit <- suppressWarnings(base::as.numeric(hcobject[["global_settings"]][["range_GFC"]]))
    if (!base::is.finite(fc_limit) || fc_limit <= 0) fc_limit <- 2
    pseudo <- 1e-08
    set_mats <- list()
    for (z in base::seq_len(base::length(hcobject[["layer_specific_outputs"]]))) {
      set_name <- base::paste0("set", z)
      sm <- prepare_group_mean_matrix(set_name)
      if (base::is.null(sm)) next
      ov <- base::intersect(net_genes, base::rownames(sm))
      if (base::length(ov) == 0) next
      for (i in base::seq_len(base::nrow(fc_def))) {
        num <- base::as.character(fc_def$numerator[i]); den <- base::as.character(fc_def$denominator[i]); cmp <- base::as.character(fc_def$comparison[i])
        if (!(num %in% base::colnames(sm) && den %in% base::colnames(sm))) next
        comparison_found[[cmp]] <- TRUE
        fc_vals <- suppressWarnings(base::log2((as.numeric(sm[ov, num]) + pseudo) / (as.numeric(sm[ov, den]) + pseudo)))
        fc_vals[!base::is.finite(fc_vals)] <- NA_real_
        fc_vals[fc_vals > fc_limit] <- fc_limit
        fc_vals[fc_vals < (-fc_limit)] <- -fc_limit
        full <- base::matrix(NA_real_, nrow = base::length(net_genes), ncol = 1, dimnames = list(net_genes, cmp))
        full[ov, 1] <- fc_vals
        set_mats[[base::length(set_mats) + 1]] <- full
      }
    }
    if (base::length(set_mats) == 0 || !base::any(unlist(comparison_found))) stop("None of the requested `fc_comparisons` could be computed.")
    activity_mat <- base::do.call(base::cbind, set_mats)
    activity_mat <- collapse_dup_rows(activity_mat)
    activity_mat <- collapse_dup_cols(activity_mat)
    used <- base::as.character(fc_def$comparison)[base::as.character(fc_def$comparison) %in% base::colnames(activity_mat)]
    used <- base::c(used, base::setdiff(base::colnames(activity_mat), used))
    activity_mat <- activity_mat[, used, drop = FALSE]
    fc_summary <- list(used = used, missing = base::names(comparison_found)[!comparison_found])
  }

  keep_rows <- base::rowSums(!base::is.na(activity_mat)) > 0
  activity_mat <- activity_mat[keep_rows, , drop = FALSE]
  if (base::nrow(activity_mat) == 0 || base::ncol(activity_mat) == 0) stop("Activity matrix is empty.")
  condition_levels <- base::colnames(activity_mat)

  empty_df <- function(with_condition = FALSE) {
    cols <- base::c("resource", "database", "cluster", "module_label", if (with_condition) "condition", "rank", "cell_type", "score", "abs_score", "pvalue", "qvalue", "direction", "n_conditions", "n_genes", "n_targets")
    out <- base::as.data.frame(stats::setNames(base::replicate(base::length(cols), base::vector(mode = "list", length = 0), simplify = FALSE), cols), stringsAsFactors = FALSE)
    out
  }
  combine_non_empty <- function(lst, with_condition = FALSE) {
    lst <- lst[!base::vapply(lst, function(x) base::is.null(x) || base::nrow(x) == 0, FUN.VALUE = base::logical(1))]
    if (base::length(lst) == 0) return(empty_df(with_condition = with_condition))
    out <- base::do.call(base::rbind, lst)
    base::rownames(out) <- NULL
    out
  }

  custom_t2g <- .hc_load_custom_gmt_term2gene(
    custom_gmt_files = custom_gmt_files,
    default_prefix = "CustomCellType",
    uppercase_genes = TRUE,
    add_database_prefix = TRUE,
    title_case_names = FALSE
  )
  custom_gmt_paths <- base::attr(custom_t2g, "paths")
  if (base::length(databases) > 0 && base::length(custom_t2g) > 0) {
    overlap_db <- base::intersect(databases, base::names(custom_t2g))
    if (base::length(overlap_db) > 0) {
      stop(
        "Overlapping Enrichr and custom GMT database names: ",
        base::paste(overlap_db, collapse = ", "),
        ". Rename custom GMT entries (names(custom_gmt_files))."
      )
    }
  }

  db_t2g <- list()
  if (base::length(databases) > 0) {
    .hc_enrichr_db_table(refresh = isTRUE(refresh_db))
    db_t2g <- lapply(
      databases,
      function(db) .hc_enrichr_fetch_library_term2gene(
        database = db,
        refresh = isTRUE(refresh_db),
        min_genes_per_term = min_term_genes
      )$term2gene
    )
    base::names(db_t2g) <- databases
  }
  if (base::length(custom_t2g) > 0) {
    db_t2g <- base::c(db_t2g, custom_t2g)
  }
  db_names <- base::names(db_t2g)
  if (base::length(db_names) == 0) stop("No databases available. Provide `databases` and/or `custom_gmt_files`.")

  selected_by_db <- stats::setNames(vector("list", base::length(db_names)), db_names)
  significant_by_db <- stats::setNames(vector("list", base::length(db_names)), db_names)
  selected_by_condition_by_db <- stats::setNames(vector("list", base::length(db_names)), db_names)
  significant_by_condition_by_db <- stats::setNames(vector("list", base::length(db_names)), db_names)

  for (db_nm in db_names) {
    t2g <- db_t2g[[db_nm]]
    if (base::is.null(t2g) || base::nrow(t2g) == 0) next
    src <- if (mode == "fine") base::as.character(t2g$term) else .hc_assign_coarse_celltype(base::as.character(t2g$term_raw), coarse_map)
    net_db <- base::data.frame(source = src, target = base::toupper(base::trimws(base::as.character(t2g$gene))), mor = 1, stringsAsFactors = FALSE)
    net_db <- net_db[!base::is.na(net_db$source) & net_db$source != "" & !base::is.na(net_db$target) & net_db$target != "", , drop = FALSE]
    if (mode == "coarse" && !isTRUE(coarse_include_other)) net_db <- net_db[net_db$source != "Other", , drop = FALSE]
    net_db <- base::unique(net_db)
    net_db <- net_db[net_db$target %in% base::rownames(activity_mat), , drop = FALSE]
    if (base::nrow(net_db) == 0) next

    sel_rows <- list(); sig_rows <- list(); selc_rows <- list(); sigc_rows <- list()
    for (cl in clusters) {
      genes <- .hc_extract_cluster_genes(cluster_info, cl)
      genes <- genes[genes %in% base::rownames(activity_mat)]
      if (base::length(genes) == 0) next
      net_mod <- net_db[net_db$target %in% genes, , drop = FALSE]
      if (base::nrow(net_mod) == 0) next
      raw <- .hc_ui_run_decouple(mat = activity_mat, network = net_mod, method = method, minsize = minsize)
      if (base::is.null(raw) || base::nrow(raw) == 0) next

      agg <- .hc_ui_summarize_decouple_result(raw, padj = padj)
      if (base::nrow(agg) > 0) {
        agg$resource <- "CellType"; agg$database <- db_nm; agg$cluster <- cl; agg$module_label <- if (cl %in% base::names(module_label_map)) base::as.character(module_label_map[[cl]]) else cl
        agg$cell_type <- base::as.character(agg$term); agg$n_genes <- base::length(base::unique(genes)); agg$n_targets <- base::length(base::unique(net_mod$target))
        agg <- agg[base::order(agg$qvalue, -agg$abs_score, agg$cell_type), , drop = FALSE]
        sig <- agg[!base::is.na(agg$qvalue) & agg$qvalue <= qval, , drop = FALSE]
        if (base::nrow(sig) > 0) {
          sel <- sig[base::seq_len(base::min(top, base::nrow(sig))), , drop = FALSE]
          sel$rank <- base::seq_len(base::nrow(sel)); sig$rank <- base::seq_len(base::nrow(sig))
          sel_rows[[base::length(sel_rows) + 1]] <- sel
          sig_rows[[base::length(sig_rows) + 1]] <- sig
        }
      }

      aggc <- .hc_ui_summarize_decouple_by_condition(raw, padj = padj)
      if (base::nrow(aggc) > 0) {
        aggc$resource <- "CellType"; aggc$database <- db_nm; aggc$cluster <- cl; aggc$module_label <- if (cl %in% base::names(module_label_map)) base::as.character(module_label_map[[cl]]) else cl
        aggc$cell_type <- base::as.character(aggc$term); aggc$n_genes <- base::length(base::unique(genes)); aggc$n_targets <- base::length(base::unique(net_mod$target))
        aggc <- aggc[base::order(aggc$condition, aggc$qvalue, -aggc$abs_score, aggc$cell_type), , drop = FALSE]
        sigc <- aggc[!base::is.na(aggc$qvalue) & aggc$qvalue <= qval, , drop = FALSE]
        if (base::nrow(sigc) > 0) {
          selc <- base::do.call(base::rbind, base::lapply(base::split(sigc, base::as.character(sigc$condition)), function(x) {
            x <- x[base::order(x$qvalue, -x$abs_score, x$cell_type), , drop = FALSE]
            x <- x[base::seq_len(base::min(top, base::nrow(x))), , drop = FALSE]
            x$rank <- base::seq_len(base::nrow(x))
            x
          }))
          if (!base::is.null(selc) && base::nrow(selc) > 0) selc_rows[[base::length(selc_rows) + 1]] <- selc
          sigc$rank <- stats::ave(sigc$qvalue, sigc$condition, FUN = function(v) base::rank(v, ties.method = "first"))
          sigc_rows[[base::length(sigc_rows) + 1]] <- sigc
        }
      }
    }
    selected_by_db[[db_nm]] <- combine_non_empty(sel_rows, with_condition = FALSE)
    significant_by_db[[db_nm]] <- combine_non_empty(sig_rows, with_condition = FALSE)
    selected_by_condition_by_db[[db_nm]] <- combine_non_empty(selc_rows, with_condition = TRUE)
    significant_by_condition_by_db[[db_nm]] <- combine_non_empty(sigc_rows, with_condition = TRUE)
  }

  selected_all <- combine_non_empty(selected_by_db, with_condition = FALSE)
  significant_all <- combine_non_empty(significant_by_db, with_condition = FALSE)
  selected_by_condition_all <- combine_non_empty(selected_by_condition_by_db, with_condition = TRUE)
  significant_by_condition_all <- combine_non_empty(significant_by_condition_by_db, with_condition = TRUE)

  categories_by_db <- stats::setNames(vector("list", base::length(db_names)), db_names)
  celltype_order_by_db <- stats::setNames(vector("list", base::length(db_names)), db_names)
  for (db_nm in db_names) {
    sel_db <- selected_by_db[[db_nm]]
    if (base::is.null(sel_db) || base::nrow(sel_db) == 0) {
      categories_by_db[[db_nm]] <- base::data.frame(count = base::numeric(0), cell_type = base::character(0), cluster = base::character(0), hits = base::numeric(0), stringsAsFactors = FALSE)
      celltype_order_by_db[[db_nm]] <- base::data.frame(cell_type = base::character(0), total_count = base::numeric(0), best_q = base::numeric(0), n_clusters = base::numeric(0), stringsAsFactors = FALSE)
      next
    }
    sel_db$count <- suppressWarnings(base::as.numeric(sel_db$abs_score)); sel_db$count[!base::is.finite(sel_db$count)] <- 0
    tot <- stats::aggregate(count ~ cell_type, data = sel_db, FUN = base::sum)
    bq <- stats::aggregate(qvalue ~ cell_type, data = sel_db, FUN = base::min)
    nc <- stats::aggregate(cluster ~ cell_type, data = sel_db, FUN = function(x) base::length(base::unique(x)))
    ord <- dplyr::left_join(tot, bq, by = "cell_type")
    base::colnames(ord)[base::colnames(ord) == "count"] <- "total_count"
    base::colnames(ord)[base::colnames(ord) == "qvalue"] <- "best_q"
    ord <- dplyr::left_join(ord, nc, by = "cell_type")
    base::colnames(ord)[base::colnames(ord) == "cluster"] <- "n_clusters"
    ord <- ord[base::order(ord$best_q, -ord$n_clusters, -ord$total_count, ord$cell_type), , drop = FALSE]
    cell_types <- base::as.character(ord$cell_type)
    pct_rows <- list()
    for (cl in clusters) {
      tmp <- sel_db[sel_db$cluster == cl, c("cell_type", "count"), drop = FALSE]
      if (base::nrow(tmp) == 0) next
      tmp <- stats::aggregate(count ~ cell_type, data = tmp, FUN = base::sum)
      hits <- base::sum(tmp$count, na.rm = TRUE)
      tmp$hits <- hits
      tmp$pct <- if (hits > 0) (tmp$count / hits) * 100 else 0
      tmp$cluster <- cl
      pct_rows[[base::length(pct_rows) + 1]] <- tmp[, c("cluster", "cell_type", "pct", "hits"), drop = FALSE]
    }
    pct_df <- if (base::length(pct_rows) == 0) base::data.frame(cluster = base::character(0), cell_type = base::character(0), pct = base::numeric(0), hits = base::numeric(0), stringsAsFactors = FALSE) else base::do.call(base::rbind, pct_rows)
    grid_db <- base::expand.grid(cluster = clusters, cell_type = cell_types, stringsAsFactors = FALSE)
    cat_db <- dplyr::left_join(grid_db, pct_df, by = c("cluster", "cell_type"))
    cat_db$pct[base::is.na(cat_db$pct)] <- 0
    cat_db$hits[base::is.na(cat_db$hits)] <- 0
    cat_db$count <- cat_db$pct
    cat_db <- cat_db[, c("count", "cell_type", "cluster", "hits"), drop = FALSE]
    categories_by_db[[db_nm]] <- cat_db
    celltype_order_by_db[[db_nm]] <- ord
  }

  categories <- base::do.call(base::rbind, base::lapply(base::names(categories_by_db), function(db_nm) {
    x <- categories_by_db[[db_nm]]
    if (base::is.null(x) || base::nrow(x) == 0) return(NULL)
    x$database <- db_nm
    x[, c("database", "count", "cell_type", "cluster", "hits"), drop = FALSE]
  }))
  if (base::is.null(categories)) categories <- base::data.frame(database = base::character(0), count = base::numeric(0), cell_type = base::character(0), cluster = base::character(0), hits = base::numeric(0), stringsAsFactors = FALSE)

  db_tokens <- base::make.unique(base::vapply(db_names, .hc_sanitize_db_token, FUN.VALUE = base::character(1)), sep = "_")
  db_slot_names <- base::paste0("enriched_per_cluster_", db_tokens)
  if (!base::is.null(slot_suffix) && base::nzchar(slot_suffix)) db_slot_names <- base::paste0(db_slot_names, "_", slot_suffix)
  base::names(db_slot_names) <- db_names
  if (base::length(db_names) == 1 && annotation_slot != "auto") {
    db_slot_names[[1]] <- annotation_slot
    if (!base::is.null(slot_suffix) && base::nzchar(slot_suffix)) db_slot_names[[1]] <- base::paste0(db_slot_names[[1]], "_", slot_suffix)
  }

  sat_out <- hcobject[["satellite_outputs"]]
  if (base::is.null(sat_out) || !base::is.list(sat_out)) sat_out <- list()
  if (isTRUE(clear_previous_slots)) {
    sat_names <- base::names(sat_out)
    drop_nm <- sat_names[base::grepl("^enriched_per_cluster_", sat_names) | sat_names %in% c("enriched_per_cluster", "enriched_per_cluster2")]
    for (nm in drop_nm) sat_out[[nm]] <- NULL
  }

  hidden_common <- list(source = "celltype_activity_decoupler", databases = db_names, custom_gmt_files = custom_gmt_paths, clusters = clusters, mode = mode, top = top, qval = qval, padj = padj, activity_input = activity_input, slot_suffix = slot_suffix, clear_previous_slots = clear_previous_slots)
  base::attr(categories, "hidden") <- hidden_common
  for (db_nm in db_names) {
    slot_nm <- db_slot_names[[db_nm]]
    cat_db <- categories_by_db[[db_nm]]
    if (base::is.null(cat_db)) cat_db <- base::data.frame(count = base::numeric(0), cell_type = base::character(0), cluster = base::character(0), hits = base::numeric(0), stringsAsFactors = FALSE)
    hidden_db <- base::c(hidden_common, list(database = db_nm, label = base::paste0(db_nm, " (decoupleR)"), slot = slot_nm))
    base::attr(cat_db, "hidden") <- hidden_db
    sat_out[[slot_nm]] <- list(categories_per_cluster = cat_db)
    base::attr(sat_out[[slot_nm]][["categories_per_cluster"]], "hidden") <- hidden_db
  }

  out <- list(
    categories_per_cluster = categories,
    categories_per_cluster_by_db = categories_by_db,
    selected_celltypes = selected_all,
    selected_celltypes_by_db = selected_by_db,
    selected_celltypes_by_condition = selected_by_condition_all,
    selected_celltypes_by_condition_by_db = selected_by_condition_by_db,
    significant_terms = significant_all,
    significant_terms_by_db = significant_by_db,
    significant_terms_by_condition = significant_by_condition_all,
    significant_terms_by_condition_by_db = significant_by_condition_by_db,
    celltype_order_by_db = celltype_order_by_db,
    annotation_slots = base::unname(db_slot_names),
    annotation_slot_map = base::data.frame(database = db_names, slot = base::unname(db_slot_names), stringsAsFactors = FALSE),
    parameters = list(databases = db_names, custom_gmt_files = custom_gmt_paths, clusters = clusters, mode = mode, top = top, qval = qval, padj = padj, activity_input = activity_input, fc_comparisons = if (!base::is.null(fc_summary)) fc_summary$used else NULL, fc_comparisons_missing = if (!base::is.null(fc_summary)) fc_summary$missing else NULL, method = method, minsize = minsize, min_term_genes = min_term_genes, annotation_slot = annotation_slot, annotation_slots = base::unname(db_slot_names), slot_suffix = slot_suffix, clear_previous_slots = clear_previous_slots, coarse_include_other = coarse_include_other, condition_levels = condition_levels)
  )
  base::attr(out$categories_per_cluster, "hidden") <- hidden_common

  sat_out[["celltype_activity_decoupler"]] <- out
  hcobject[["satellite_outputs"]] <<- sat_out

  if (isTRUE(export_excel)) {
    out_dir <- .hc_find_output_dir()
    if (!base::is.null(out_dir)) {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      xlsx_tables <- list(
        selected_celltypes = selected_all,
        selected_celltypes_by_condition = selected_by_condition_all,
        significant_terms = significant_all,
        significant_terms_by_condition = significant_by_condition_all,
        categories_per_cluster = categories,
        annotation_slot_map = out$annotation_slot_map
      )
      tryCatch(openxlsx::write.xlsx(x = xlsx_tables, file = file.path(out_dir, excel_file), overwrite = TRUE), error = function(e) warning("Could not write Excel: ", conditionMessage(e)))
    }
  }

  if (isTRUE(plot_heatmap)) {
    hm_args <- list(file_name = heatmap_file_name, gene_count_mode = "none", module_label_preset = "compact", include_dynamic_enrichment_slots = TRUE)
    extra <- list(...)
    if (base::length(extra) > 0) hm_args <- utils::modifyList(hm_args, extra)
    base::do.call(plot_cluster_heatmap, hm_args)
  }

  invisible(out)
}
