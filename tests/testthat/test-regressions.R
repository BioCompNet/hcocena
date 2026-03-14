test_that("regression: utils no longer references working_director", {
  src <- paste(deparse(get("run_expression_analysis_1_body", asNamespace("hcocena"))), collapse = "\n")
  expect_false(grepl("working_director", src, fixed = TRUE))
})


test_that("regression: rho path now guards optional propr dependency", {
  src <- paste(deparse(get("pwcorr", asNamespace("hcocena"))), collapse = "\n")
  expect_true(grepl("requireNamespace(\"propr\"", src, fixed = TRUE))
})


test_that("regression: safe deep clone falls back instead of aborting", {
  clone_fun <- get(".hc_safe_deep_clone", asNamespace("hcocena"))

  x <- list(a = 1, b = list(2))
  y <- clone_fun(x, context = "test object")
  y$b[[1]] <- 3
  expect_equal(x$b[[1]], 2)

  env <- new.env(parent = emptyenv())
  env$a <- 1
  testthat::local_mocked_bindings(
    serialize = function(...) stop("mock oom"),
    .package = "base"
  )
  expect_warning(
    out <- clone_fun(env, context = "test object"),
    "Could not deep-clone test object"
  )
  expect_identical(out, env)

  expect_warning(
    out_big <- clone_fun(x, context = "large test object", max_bytes = 1),
    "Skipping deep clone of large test object"
  )
  expect_identical(out_big, x)
})


test_that("regression: cutoff stats match legacy igraph behavior on small graphs", {
  rsquaredfun <- get("rsquaredfun", asNamespace("hcocena"))
  union_find <- get(".hc_union_find_components", asNamespace("hcocena"))

  graph_df <- data.frame(
    V1 = c("g1", "g1", "g2", "g4"),
    V2 = c("g2", "g3", "g3", "g5"),
    rval = c(0.91, 0.88, 0.9, 0.95),
    stringsAsFactors = FALSE
  )

  ref_graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE, vertices = NULL)
  ref_components <- igraph::components(ref_graph)
  keep_components <- which(ref_components$csize >= 3)
  gene_to_comp <- data.frame(
    gene = names(ref_components$membership),
    component = ref_components$membership,
    stringsAsFactors = FALSE
  )
  nodes_to_remove <- dplyr::filter(gene_to_comp, !component %in% keep_components) %>%
    dplyr::pull(gene)
  ref_graph <- igraph::delete_vertices(ref_graph, nodes_to_remove)
  ref_degree <- igraph::degree(ref_graph, mode = "all")
  ref_dd <- igraph::degree_distribution(ref_graph, mode = "all", cumulative = FALSE)
  ref_prob <- ref_dd[-1]
  ref_deg_vals <- seq_len(max(ref_degree))
  nonzero <- which(ref_prob != 0)
  ref_prob <- ref_prob[nonzero]
  ref_deg_vals <- ref_deg_vals[nonzero]
  ref_r2 <- summary(stats::lm(log(ref_prob) ~ log(ref_deg_vals)))$r.squared

  uf <- union_find(from_idx = c(1L, 1L, 2L, 4L), to_idx = c(2L, 3L, 3L, 5L), n_vertices = 5L)
  expect_equal(sort(uf$csize), c(2L, 3L))
  expect_equal(uf$csize[uf$membership], c(3L, 3L, 3L, 2L, 2L))

  out <- rsquaredfun(
    graph_df = graph_df,
    cutoff = 0.88,
    print.all.plots = FALSE,
    min_nodes = 3,
    x = 1
  )

  expect_equal(out$no_of_networks, 1)
  expect_equal(out$no_nodes, 3)
  expect_equal(out$no_edges, 3)
  expect_equal(out$degree[[1]], ref_deg_vals)
  expect_equal(out$Probs[[1]], ref_prob)
  expect_equal(out$R.squared[[1]], ref_r2)
})


test_that("regression: longitudinal top terms are compacted after filtering", {
  rank_lookup <- get(".hc_longitudinal_rank_lookup", asNamespace("hcocena"))
  panel_label <- get(".hc_longitudinal_panel_label", asNamespace("hcocena"))

  top_df <- data.frame(
    module = c("M1", "M1", "M1", "M2", "M2"),
    term = c("KEGG_T1", "KEGG_T2", "KEGG_T3", "HALLMARK_U1", "HALLMARK_U2"),
    module_color = c("red", "red", "red", "blue", "blue"),
    rank = c(1, 2, 3, 1, 2),
    qvalue = c(0.01, 0.02, 0.03, 0.01, 0.02),
    stringsAsFactors = FALSE
  )
  valid_terms <- data.frame(
    module = c("M1", "M1", "M2"),
    term = c("KEGG_T1", "KEGG_T3", "HALLMARK_U2"),
    stringsAsFactors = FALSE
  )

  out <- rank_lookup(
    top_df = top_df,
    wrap_term = function(x) base::gsub("_", " ", x, fixed = TRUE),
    valid_terms = valid_terms
  )

  expect_equal(out[out$module == "M1", "term"], c("KEGG_T1", "KEGG_T3"))
  expect_equal(out[out$module == "M1", "term_rank"], c(1L, 2L))
  expect_equal(out[out$module == "M1", "term_rank_label"], c("Top 1", "Top 2"))
  expect_false(any(base::is.na(out$term_label)))
  expect_false(any(!base::nzchar(base::trimws(out$term_label))))

  fallback <- panel_label(module = NA_character_, rank_label = NA_character_, term_display = NA_character_)
  expect_identical(fallback, "Module: Top\nTerm unavailable")
})


test_that("regression: longitudinal enrichment uses observed module label/color pairs", {
  harmonize <- get(".hc_harmonize_module_lookup_with_enrichment", asNamespace("hcocena"))

  module_lookup <- data.frame(
    module = c("M1", "M2", "M3"),
    module_color = c("gold", "lightblue", "wheat"),
    stringsAsFactors = FALSE
  )
  enrich_df <- data.frame(
    cluster = c("gold", "darkgreen", "cyan", "wheat"),
    module_label = c("M1", "M7", "M10", ""),
    stringsAsFactors = FALSE
  )

  out <- harmonize(module_lookup = module_lookup, enrich_df = enrich_df)

  expect_equal(out$module_color[match("M1", out$module)], "gold")
  expect_equal(out$module_color[match("M7", out$module)], "darkgreen")
  expect_equal(out$module_color[match("M10", out$module)], "cyan")
  expect_equal(out$module[match("wheat", out$module_color)], "M3")
})


test_that("regression: longitudinal plotting drops panel labels outside final facet levels", {
  align_panels <- get(".hc_align_longitudinal_panel_factors", asNamespace("hcocena"))

  df <- data.frame(
    panel_label = c("M1: Top 1\nA", "M1: Top 2\nB", NA, "NA"),
    term_display = c("A", "B", "A", "B"),
    value = 1:4,
    stringsAsFactors = FALSE
  )

  out <- align_panels(
    df = df,
    panel_levels = c("M1: Top 1\nA"),
    term_display_levels = c("A")
  )

  expect_equal(base::as.character(out$panel_label), "M1: Top 1\nA")
  expect_equal(base::as.character(out$term_display), "A")
  expect_equal(out$value, 1L)
})


test_that("regression: emmeans summaries are coerced to numeric safely", {
  normalize_emm <- get(".hc_normalize_emmeans_summary", asNamespace("hcocena"))

  df <- data.frame(
    response = c("0.4", "0.7"),
    lower.CL = c("nonEst", "0.2"),
    upper.CL = c("0.9", "1.1"),
    stringsAsFactors = FALSE
  )

  out <- normalize_emm(df)

  expect_equal(out$emmean, c(0.4, 0.7))
  expect_true(is.na(out$lower.CL[[1]]))
  expect_equal(out$lower.CL[[2]], 0.2)
  expect_equal(out$upper.CL, c(0.9, 1.1))

  out$lower.CL <- suppressWarnings(as.numeric(out$lower.CL))
  out$upper.CL <- suppressWarnings(as.numeric(out$upper.CL))
  ci_missing <- !is.finite(out$lower.CL) | !is.finite(out$upper.CL)
  out$lower.CL[ci_missing] <- out$emmean[ci_missing]
  out$upper.CL[ci_missing] <- out$emmean[ci_missing]

  expect_equal(out$lower.CL[[1]], out$emmean[[1]])
  expect_equal(out$upper.CL[[1]], out$emmean[[1]])
})


test_that("regression: heatmap col_order ignores missing columns safely", {
  resolve_cols <- get(".hc_resolve_heatmap_col_order", asNamespace("hcocena"))

  expect_warning(
    out <- resolve_cols(
      mat_cols = c("MC1_T1", "MC1_T2", "MC2_T1"),
      requested_order = c("old_group", "MC2_T1"),
      context = "regrouped cluster heatmap"
    ),
    "Ignoring 1 `col_order` entries"
  )

  expect_equal(out, c("MC2_T1", "MC1_T1", "MC1_T2"))
})


test_that("regression: enrichment panel storage defaults to memory-saving for multi-db runs", {
  resolve_mode <- get(".hc_resolve_panel_storage_mode", asNamespace("hcocena"))
  has_draw_obj <- get(".hc_functional_enrichment_has_draw_object", asNamespace("hcocena"))

  expect_identical(resolve_mode("auto", n_databases = 1L), "always")
  expect_identical(resolve_mode("auto", n_databases = 3L), "never")
  expect_identical(resolve_mode("always", n_databases = 3L), "always")

  expect_true(has_draw_obj(list(p = structure(list(), class = "HeatmapList"))))
  expect_false(has_draw_obj(list(
    p = NULL,
    hc_heatmap = NULL,
    enrichment_plot = NULL,
    panel_objects_stored = FALSE
  )))
})


test_that("regression: large result stores are no longer mirrored across legacy slots", {
  fun_enrich_src <- paste(deparse(get("functional_enrichment", asNamespace("hcocena"))), collapse = "\n")
  heatmap_src <- paste(deparse(get("plot_cluster_heatmap", asNamespace("hcocena"))), collapse = "\n")
  heatmap_new_src <- paste(deparse(get("plot_cluster_heatmap_new", asNamespace("hcocena"))), collapse = "\n")
  upstream_src <- paste(deparse(get("upstream_inference", asNamespace("hcocena"))), collapse = "\n")
  knowledge_src <- paste(deparse(get("plot_enrichment_upstream_network", asNamespace("hcocena"))), collapse = "\n")
  network_plot_src <- paste(
    deparse(get(".hc_plot_integrated_network_legacy_driver", asNamespace("hcocena"))),
    collapse = "\n"
  )
  gfc_network_src <- paste(
    deparse(get(".hc_plot_GFC_network_legacy_driver", asNamespace("hcocena"))),
    collapse = "\n"
  )

  expect_false(grepl(
    'hcobject[["satellite_outputs"]][["enrichments"]] <<- hcobject[["integrated_output"]][["enrichments"]]',
    fun_enrich_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["enrichments"]][["all_enrichments_all_dbs"]] <<-',
    fun_enrich_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["cluster_calc"]][["module_gene_list"]] <<-',
    heatmap_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    '[["heatmap_matrix"]] <<- mat_heatmap',
    heatmap_new_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    '[["heatmap_row_order"]] <<- final_row_order',
    heatmap_new_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    '[["heatmap_column_order"]] <<- final_col_order',
    heatmap_new_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["upstream_inference"]] <<- output',
    upstream_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["knowledge_network"]] <<- output',
    knowledge_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["cluster_calc"]][["network_col_by_module"]] <<- network',
    network_plot_src,
    fixed = TRUE
  ))
  expect_false(grepl(
    'hcobject[["integrated_output"]][["cluster_calc"]][["labelled_network"]] <<- network2',
    network_plot_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    'if (isTRUE(store_plot))',
    network_plot_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    'sat[["network_col_by_module"]]',
    gfc_network_src,
    fixed = TRUE
  ))
  expect_true(grepl(
    'network <- hcobject[["integrated_output"]][["merged_net"]]',
    gfc_network_src,
    fixed = TRUE
  ))
})


test_that("regression: lightweight heatmap cache works without ComplexHeatmap object", {
  cache_info <- get(".hc_heatmap_cache_info", asNamespace("hcocena"))
  select_col_order <- get(".hc_select_heatmap_col_order", asNamespace("hcocena"))
  llm_heatmap_info <- get(".hc_llm_heatmap_info", asNamespace("hcocena"))
  llm_capture <- get(".hc_llm_capture_combined_heatmap_grob", asNamespace("hcocena"))
  plot_heatmap <- get("plot_cluster_heatmap", asNamespace("hcocena"))
  plot_heatmap_new <- get("plot_cluster_heatmap_new", asNamespace("hcocena"))
  plot_network <- get("plot_integrated_network", asNamespace("hcocena"))
  fun_enrich <- get("functional_enrichment", asNamespace("hcocena"))
  up_inf <- get("upstream_inference", asNamespace("hcocena"))
  knowledge_plot <- get("plot_enrichment_upstream_network", asNamespace("hcocena"))
  llm_plot <- get("hc_plot_module_function_llm", asNamespace("hcocena"))

  cluster_calc <- list(
    heatmap_matrix = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("M1", "M2"), c("T1", "T2"))
    ),
    heatmap_row_order = c("M2", "M1"),
    heatmap_column_order = c("T2", "T1")
  )

  info <- cache_info(cluster_calc)

  expect_equal(info$row_order, c("M2", "M1"))
  expect_equal(info$col_order, c("T2", "T1"))
  expect_equal(base::rownames(info$matrix), c("M1", "M2"))
  expect_equal(
    select_col_order(
      available_cols = c("A", "B", "C"),
      plot_order = c("C"),
      main_order = c("B", "A"),
      fallback_order = c("A", "C")
    ),
    c("C", "A", "B")
  )
  expect_equal(
    select_col_order(
      available_cols = c("A", "B", "C"),
      main_order = c("B", "A"),
      fallback_order = c("C")
    ),
    c("B", "A", "C")
  )
  expect_equal(
    select_col_order(
      available_cols = c("A", "B", "C"),
      main_order = c("Z"),
      fallback_order = c("C", "A")
    ),
    c("C", "A", "B")
  )
  expect_identical(formals(plot_heatmap)$return_HM, FALSE)
  expect_identical(formals(plot_heatmap_new)$return_HM, FALSE)
  expect_identical(formals(plot_network)$store_plot, FALSE)
  expect_true("heatmap_col_order" %in% names(formals(fun_enrich)))
  expect_true("heatmap_col_order" %in% names(formals(up_inf)))
  expect_true("heatmap_col_order" %in% names(formals(knowledge_plot)))
  expect_true("heatmap_col_order" %in% names(formals(llm_plot)))

  hc <- methods::new("HCoCenaExperiment")
  hc@integration@cluster <- S4Vectors::SimpleList(
    heatmap_matrix = matrix(
      c(-1, 0.5, 1, -0.25),
      nrow = 2,
      dimnames = list(c("red", "blue"), c("T1", "T2"))
    ),
    heatmap_row_order = c("red", "blue"),
    heatmap_column_order = c("T2", "T1"),
    module_label_map = c(red = "M1", blue = "M2"),
    module_label_fontsize = 9,
    module_label_pt_size = 0.3,
    module_box_width_cm = 0.9
  )
  info2 <- llm_heatmap_info(hc)
  expect_true(isTRUE(info2$draw_supported))
  expect_equal(info2$col_order, c("T2", "T1"))
  expect_equal(info2$module_order, c("M1", "M2"))

  summary_tbl <- data.frame(
    module = c("M1", "M2"),
    module_color = c("red", "blue"),
    term_plot = c("alpha process", "beta process"),
    text_color = c("#111111", "#222222"),
    label_color = c("white", "white"),
    stringsAsFactors = FALSE
  )
  grob <- llm_capture(
    heatmap_info = info2,
    summary_tbl = summary_tbl,
    max_chars = 90,
    text_size = 4
  )
  expect_s3_class(grob, "grob")
})


test_that("regression: .hc_run_legacy resolves functions from the legacy env", {
  run_legacy <- get(".hc_run_legacy", asNamespace("hcocena"))
  hc <- methods::new("HCoCenaExperiment")
  legacy_env <- new.env(parent = baseenv())
  legacy_env$target_env <- legacy_env

  fun <- function(gene_sets = "Hallmark", heatmap_col_order = NULL) {
    base::assign("hcobject", list(used = heatmap_col_order), envir = target_env)
    invisible(NULL)
  }
  environment(fun) <- legacy_env
  legacy_env$functional_enrichment <- fun
  legacy_env$hcobject <- list(initial = TRUE)

  global_fun <- function(gene_sets = "Hallmark") {
    stop("wrong global function selected")
  }
  base::assign("functional_enrichment", global_fun, envir = .GlobalEnv)
  on.exit(base::rm("functional_enrichment", envir = .GlobalEnv), add = TRUE)

  testthat::local_mocked_bindings(
    as_hcobject = function(hc) list(initial = TRUE),
    as_hcocena = function(x) x,
    .hc_bind_legacy_hcobject = function(hcobject, envo = legacy_env) {
      base::assign("hcobject", hcobject, envir = envo)
      list(envo = envo, had_existing = FALSE, old_hcobject = NULL, binding_locked = FALSE)
    },
    .hc_restore_legacy_hcobject = function(state) invisible(NULL),
    .package = "hcocena"
  )

  out <- run_legacy(
    hc,
    fun = "functional_enrichment",
    gene_sets = "Hallmark",
    heatmap_col_order = c("T1", "T2"),
    envo = legacy_env
  )

  expect_equal(out$used, c("T1", "T2"))
})


test_that("regression: vllm think blocks are stripped before JSON parsing", {
  strip_think <- get(".hc_llm_strip_think_blocks", asNamespace("hcocena"))
  strip_fences <- get(".hc_llm_strip_json_fences", asNamespace("hcocena"))

  raw_text <- paste(
    "<think>",
    "First inspect the genes and reason privately.",
    "</think>",
    "```json",
    '{"general_processes":"interferon signaling","contextual_state":"interferon-high inflammatory state","key_regulators":"STAT1 / IRF7 / IRF9"}',
    "```",
    sep = "\n"
  )

  cleaned <- strip_fences(strip_think(raw_text))

  expect_false(grepl("<think>", cleaned, fixed = TRUE))
  expect_false(grepl("</think>", cleaned, fixed = TRUE))
  expect_true(grepl('"general_processes":"interferon signaling"', cleaned, fixed = TRUE))
})


test_that("regression: vllm gets a longer default timeout", {
  resolve_timeout <- get(".hc_llm_resolve_timeout", asNamespace("hcocena"))

  expect_equal(
    resolve_timeout(timeout_sec = 60, llm = "openai", timeout_was_missing = TRUE),
    60
  )
  expect_equal(
    resolve_timeout(timeout_sec = 60, llm = "vllm", timeout_was_missing = TRUE),
    300
  )
  expect_null(
    resolve_timeout(timeout_sec = 0, llm = "vllm", timeout_was_missing = FALSE)
  )
})


test_that("regression: llm request helpers use ellmer backends", {
  gemini_src <- paste(deparse(get(".hc_llm_request_gemini", asNamespace("hcocena"))), collapse = "\n")
  openai_src <- paste(deparse(get(".hc_llm_request_openai", asNamespace("hcocena"))), collapse = "\n")
  vllm_src <- paste(deparse(get(".hc_llm_request_vllm", asNamespace("hcocena"))), collapse = "\n")

  expect_true(grepl("ellmer::chat_google_gemini", gemini_src, fixed = TRUE))
  expect_true(grepl("ellmer::chat_openai", openai_src, fixed = TRUE))
  expect_true(grepl("ellmer::chat_vllm", vllm_src, fixed = TRUE))
  expect_true(grepl("credentials = .hc_llm_api_key_credentials", gemini_src, fixed = TRUE))
  expect_true(grepl("credentials = .hc_llm_api_key_credentials", openai_src, fixed = TRUE))
  expect_true(grepl("credentials = .hc_llm_api_key_credentials", vllm_src, fixed = TRUE))
  expect_true(grepl("run_request\\(include_temperature = FALSE\\)", gemini_src))
  expect_true(grepl("max_tokens = 8000", vllm_src, fixed = TRUE))
  expect_true(grepl("chat_template_kwargs = list", vllm_src, fixed = TRUE))
  expect_true(grepl("enable_thinking = FALSE", vllm_src, fixed = TRUE))
})


test_that("regression: llm summary builder is robust for error-only results", {
  summary_fun <- get(".hc_llm_summary_from_results", asNamespace("hcocena"))
  err_res <- list(
    M1 = list(
      label = "M1",
      module = "M1",
      llm = "gemini",
      model = "gemini-3.1-pro-preview",
      gene_count_input = 10L,
      gene_count_sent = 0L,
      truncated = FALSE,
      status = "error",
      error_message = "HTTP 400",
      response = list(
        general_processes = NA_character_,
        contextual_state = NA_character_,
        key_regulators = "HTTP 400"
      ),
      timestamp = "2026-03-12 12:00:00"
    )
  )

  out <- summary_fun(err_res, hc = NULL)

  expect_equal(nrow(out), 1)
  expect_equal(out$module, "M1")
  expect_equal(out$status, "error")
  expect_equal(out$error_message, "HTTP 400")
  expect_equal(out$gene_count_sent, 0L)
})


test_that("regression: llm module order uses natural ordering for module='all'", {
  natural_order <- get(".hc_llm_natural_module_order", asNamespace("hcocena"))
  module_order <- get(".hc_llm_module_order", asNamespace("hcocena"))
  resolve_modules <- get(".hc_llm_resolve_modules", asNamespace("hcocena"))

  expect_equal(
    natural_order(c("M1", "M10", "M2", "M3")),
    c("M1", "M2", "M3", "M10")
  )

  hc <- methods::new("HCoCenaExperiment")
  hc@satellite <- S4Vectors::SimpleList(list(
    module_gene_list = data.frame(
      module = c("M1", "M10", "M2", "M3"),
      genes = c("A", "B", "C", "D"),
      stringsAsFactors = FALSE
    )
  ))

  expect_equal(
    module_order(hc),
    c("M1", "M2", "M3", "M10")
  )
  expect_equal(
    resolve_modules(hc, module = "all"),
    c("M1", "M2", "M3", "M10")
  )
})


test_that("regression: plot heatmaps default to main order unless clustering is enabled", {
  prepare_cols <- get(".hc_prepare_plot_heatmap_columns", asNamespace("hcocena"))

  mat <- matrix(
    c(1, 0, 2,
      2, 1, 0,
      3, 2, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("M1", "M2", "M3"), c("T3", "T1", "T2"))
  )

  ordered <- prepare_cols(
    mat = mat,
    cluster_columns = FALSE,
    plot_order = NULL,
    main_order = c("T1", "T2", "T3"),
    fallback_order = c("T3", "T2", "T1"),
    context = "test heatmap"
  )
  expect_equal(colnames(ordered$mat), c("T1", "T2", "T3"))
  expect_null(ordered$col_dend)

  clustered <- prepare_cols(
    mat = mat,
    cluster_columns = TRUE,
    plot_order = c("T1", "T2", "T3"),
    main_order = c("T1", "T2", "T3"),
    fallback_order = c("T3", "T2", "T1"),
    context = "test heatmap"
  )
  expect_setequal(colnames(clustered$mat), c("T1", "T2", "T3"))
  expect_true(inherits(clustered$col_dend, "dendrogram"))

  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena:::upstream_inference)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena:::plot_enrichment_upstream_network)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena::hc_upstream_inference)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena::hc_plot_enrichment_upstream_network)))
  expect_true("heatmap_cluster_columns" %in% names(formals(hcocena::hc_plot_module_function_llm)))
  expect_identical(formals(hcocena:::plot_cluster_heatmap)$cluster_columns, FALSE)
  expect_identical(formals(hcocena:::plot_cluster_heatmap_new)$cluster_columns, FALSE)
  expect_identical(formals(hcocena:::replot_cluster_heatmap)$cluster_columns, FALSE)
  expect_identical(formals(hcocena:::change_grouping_parameter)$cluster_columns, FALSE)
  expect_identical(formals(hcocena::hc_change_grouping_parameter)$cluster_columns, FALSE)
})


test_that("regression: celltype annotation matrix shows module labels instead of raw colors", {
  hc <- methods::new("HCoCenaExperiment")
  hc@integration@cluster <- S4Vectors::SimpleList()
  hc@integration@cluster[["module_prefix"]] <- "M"
  hc@integration@cluster[["module_label_map"]] <- c(yellow = "M1", blue = "M2")
  hc@satellite <- S4Vectors::SimpleList()
  hc@satellite[["celltype_annotation"]] <- list(
    selected_celltypes = base::data.frame(
      cluster = c("yellow", "blue"),
      cell_type = c("B cell", "T cell"),
      pct = c(55, 72),
      count = c(11, 9),
      stringsAsFactors = FALSE
    )
  )

  out <- hcocena::hc_plot_celltype_annotation_matrix(hc, return_data = TRUE)

  expect_equal(base::rownames(out$matrix), c("M1", "M2"))
  expect_setequal(base::unique(base::as.character(out$long_data$module_color)), c("yellow", "blue"))
  expect_true(inherits(out$plot, "ggplot"))
})
