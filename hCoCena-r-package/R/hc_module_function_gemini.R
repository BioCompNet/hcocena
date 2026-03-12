#' AI-assisted module function summary via Gemini, ChatGPT, or vLLM
#'
#' Uses a large language model API to summarize the likely overarching
#' biological function of a module or free gene set. This is an AI-assisted
#' interpretation step, not a statistical enrichment test.
#'
#' The function can use genes from an `HCoCenaExperiment` module or a manually
#' supplied gene vector. Optional free-text context can be added to steer the
#' interpretation towards a disease, tissue, perturbation, timepoint, or other
#' experimental background.
#'
#' By default, results are stored in `hc@satellite[[slot_name]][[label]]` when
#' an `hc` object is provided. A compact table is also written to
#' `hc@satellite[[paste0(slot_name, "_summary")]]`.
#'
#' The function name is retained for backward compatibility, but it can now use
#' Gemini, OpenAI-compatible APIs, and local vLLM servers via `llm`.
#'
#' @param hc Optional `HCoCenaExperiment`. Required when `module` is used or
#'   when `save_to_hc = TRUE`.
#' @param module Optional character vector naming one or more modules. Use
#'   `"all"` to summarize all modules present in
#'   `hc@satellite$module_gene_list`.
#' @param genes Optional character vector of gene symbols. Use this instead of
#'   `module`.
#' @param context Optional character string or vector with biological context
#'   for backward compatibility.
#' @param biological_context Optional character string or vector with biological
#'   context for the interpretation. Preferred over `context`.
#' @param label Optional label used for storing the result. Defaults to the
#'   module name or `"custom_geneset"`. Can only be used for a single module or
#'   a free gene set.
#' @param llm LLM provider. Supported values are `"gemini"`, `"openai"`,
#'   `"chatgpt"` (alias for `"openai"`), or `"vllm"` for a local
#'   OpenAI-compatible vLLM server.
#' @param api_key API key for the selected provider. If omitted, the function
#'   uses `GEMINI_API_KEY` for Gemini, `OPENAI_API_KEY` for OpenAI, and
#'   `VLLM_API_KEY` for vLLM. For local vLLM, the fallback is `"EMPTY"`.
#' @param gemini_model Optional Gemini model code. Only used when
#'   `llm = "gemini"`. Defaults to `"gemini-2.5-pro"`.
#' @param vllm_model Optional vLLM model code. Only used when `llm = "vllm"`.
#'   Defaults to `"Qwen/Qwen2.5-VL-32B-Instruct"`.
#' @param vllm_base_url Base URL for the local OpenAI-compatible vLLM server.
#'   Only used when `llm = "vllm"`. Defaults to `"http://localhost:8000/v1"`.
#' @param model Optional generic model code. Mainly useful for OpenAI, or as a
#'   provider-agnostic override.
#' @param max_genes Maximum number of genes sent to the API. Use `Inf` or
#'   `NULL` to send all genes. Default is `Inf`.
#' @param temperature Generation temperature. Default is `0.2`.
#' @param timeout_sec Request timeout in seconds. Use `0` or `NULL` to disable
#'   the timeout. Default is `60` for hosted APIs and `300` for local `vllm`
#'   requests unless overridden explicitly.
#' @param pause_sec Pause in seconds between module requests. Useful for
#'   `module = "all"`. Use `0` or `NULL` for no pause. Default is `0`.
#' @param continue_on_error Logical. If `TRUE`, continue with the next module
#'   when one request fails and store the error in the result summary.
#' @param save_to_hc Logical. If `TRUE`, store the result in `hc@satellite` and
#'   return updated `hc`. Default is `TRUE` when `hc` is provided.
#' @param slot_name Satellite slot name used for storage. Default is
#'   `"llm_module_function"`.
#' @param system_instruction Optional custom system instruction for the model.
#' @param verbose Logical. If `TRUE`, print short progress messages.
#' @param ... Used by the backward-compatible wrapper aliases
#'   `hc_module_function_gemini()` and `hc_module_function_vllm()`.
#'
#' @return Updated `HCoCenaExperiment` if `save_to_hc = TRUE`; otherwise either
#'   a single result list or, for multiple modules, a list with `results` and
#'   `summary`.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- hc_module_function_llm(
#'   genes = c("STAT1", "IRF7", "CXCL10", "GBP1", "IFI44L"),
#'   context = "Interferon-driven blood module in acute viral infection",
#'   llm = "openai",
#'   api_key = Sys.getenv("OPENAI_API_KEY"),
#'   save_to_hc = FALSE
#' )
#'
#' hc <- hc_module_function_llm(
#'   hc,
#'   module = "all",
#'   llm = "gemini",
#'   api_key = Sys.getenv("GEMINI_API_KEY")
#' )
#'
#' hc <- hc_module_function_llm(
#'   hc,
#'   module = "all",
#'   biological_context = "maturation of monocytes in preterm infants over the first year of life",
#'   llm = "vllm",
#'   vllm_model = "Qwen/Qwen2.5-VL-32B-Instruct",
#'   vllm_base_url = "http://localhost:8000/v1"
#' )
#'
#' hc@satellite$llm_module_function_summary
#' }
hc_module_function_llm <- function(hc = NULL,
                                   module = NULL,
                                   genes = NULL,
                                   context = NULL,
                                   biological_context = NULL,
                                   label = NULL,
                                   llm = c("gemini", "openai", "chatgpt", "vllm"),
                                   api_key = NULL,
                                   gemini_model = NULL,
                                   vllm_model = NULL,
                                   vllm_base_url = NULL,
                                   model = NULL,
                                   max_genes = Inf,
                                   temperature = 0.2,
                                   timeout_sec = 60,
                                   pause_sec = 0,
                                   continue_on_error = FALSE,
                                   save_to_hc = !base::is.null(hc),
                                   slot_name = "llm_module_function",
                                   system_instruction = NULL,
                                   verbose = TRUE) {
  timeout_was_missing <- missing(timeout_sec)

  if (!base::is.null(module) && !base::is.null(genes)) {
    stop("Use either `module` or `genes`, not both.")
  }
  if (base::is.null(module) && base::is.null(genes)) {
    stop("Provide either `module` or `genes`.")
  }
  if (isTRUE(save_to_hc) && base::is.null(hc)) {
    stop("`hc` must be provided when `save_to_hc = TRUE`.")
  }

  llm <- base::tolower(base::as.character(llm[[1]]))
  if (!(llm %in% c("gemini", "openai", "chatgpt", "vllm"))) {
    stop("`llm` must be one of `gemini`, `openai`, `chatgpt`, or `vllm`.")
  }
  if (llm == "chatgpt") {
    llm <- "openai"
  }

  if (llm == "gemini" &&
      !base::is.null(gemini_model) &&
      base::nzchar(base::as.character(gemini_model[[1]]))) {
    model <- base::as.character(gemini_model[[1]])
  } else if (llm == "vllm" &&
             !base::is.null(vllm_model) &&
             base::nzchar(base::as.character(vllm_model[[1]]))) {
    model <- base::as.character(vllm_model[[1]])
  } else if (base::is.null(model) || !base::nzchar(base::as.character(model[[1]]))) {
    model <- if (llm == "gemini") {
      "gemini-2.5-pro"
    } else if (llm == "vllm") {
      "Qwen/Qwen2.5-VL-32B-Instruct"
    } else {
      "gpt-4o-mini"
    }
  } else {
    model <- base::as.character(model[[1]])
  }
  if (!base::is.character(model) || base::length(model) != 1 || !base::nzchar(model)) {
    stop("`model` must be a non-empty character scalar.")
  }
  if (!base::is.null(max_genes)) {
    if (!base::is.numeric(max_genes) || base::length(max_genes) != 1 || base::is.na(max_genes[[1]]) || max_genes[[1]] <= 0) {
      stop("`max_genes` must be NULL, Inf, or a single positive number.")
    }
  }
  if (!base::is.numeric(temperature) || base::length(temperature) != 1 || !base::is.finite(temperature)) {
    stop("`temperature` must be a single finite numeric value.")
  }
  if (!base::is.null(timeout_sec)) {
    if (!base::is.numeric(timeout_sec) || base::length(timeout_sec) != 1 || !base::is.finite(timeout_sec)) {
      stop("`timeout_sec` must be NULL or a single finite numeric value >= 0.")
    }
  }
  if (!base::is.null(pause_sec)) {
    if (!base::is.numeric(pause_sec) || base::length(pause_sec) != 1 || !base::is.finite(pause_sec) || pause_sec < 0) {
      stop("`pause_sec` must be NULL or a single non-negative number.")
    }
  }
  if (!base::is.logical(continue_on_error) || base::length(continue_on_error) != 1 || base::is.na(continue_on_error)) {
    stop("`continue_on_error` must be TRUE or FALSE.")
  }
  if (!base::is.logical(verbose) || base::length(verbose) != 1 || base::is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  api_key <- .hc_llm_resolve_api_key(api_key = api_key, llm = llm)
  vllm_base_url <- .hc_llm_resolve_vllm_base_url(vllm_base_url = vllm_base_url, llm = llm)
  timeout_sec <- .hc_llm_resolve_timeout(
    timeout_sec = timeout_sec,
    llm = llm,
    timeout_was_missing = timeout_was_missing
  )
  pause_sec <- .hc_llm_normalize_pause(pause_sec)

  if (base::is.null(system_instruction) || !base::nzchar(base::as.character(system_instruction[[1]]))) {
    system_instruction <- if (llm == "vllm") {
      paste(
        "You are a careful transcriptomics analyst.",
        "Infer the main biological program of the provided gene module.",
        "Use the biological context as background information, but do not simply repeat it.",
        "Be concise and module-specific.",
        "For contextual_state, describe the immediate transcriptional or cellular state, not the study design or cohort.",
        "Avoid repeating generic context words such as immature, maturation, preterm, infant, or monocyte unless they are essential and specifically supported by the genes.",
        "Do not write sentence starters like 'this module' or 'the genes'.",
        "Use the supplied genes and optional context only as evidence.",
        "Do not claim statistical enrichment was performed.",
        "Return valid JSON only."
      )
    } else {
      paste(
        "You are a careful transcriptomics analyst.",
        "Infer the main biological program of the provided gene module.",
        "Use the biological context only as an interpretation frame, not as text to repeat.",
        "Prioritize module-specific biology over generic cohort-level wording.",
        "Avoid generic answers such as monocyte maturation, immune differentiation, or immune activation unless the genes strongly support nothing more specific.",
        "Prefer concrete programs such as interferon signaling, antigen presentation, cell cycle, erythroid/megakaryocytic bias, platelet program, mitochondrial metabolism, phagolysosome, inflammatory signaling, ribosome biogenesis, stress response, chemotaxis, or tissue contamination if clearly supported.",
        "Make different modules distinguishable from each other.",
        "Return compact labels and compact descriptions only.",
        "Do not write sentence starters like 'this module', 'this program', or 'is best described as'.",
        "Use the supplied genes and optional context only as evidence.",
        "Do not claim statistical enrichment was performed.",
        "State uncertainty when the signal is weak or mixed.",
        "Return valid JSON only."
      )
    }
  } else {
    system_instruction <- base::as.character(system_instruction[[1]])
  }

  if (base::is.null(module)) {
    gene_infos <- list(
      list(
        label = "custom_geneset",
        module = NULL,
        genes = .hc_gemini_normalize_genes(genes)
      )
    )
  } else {
    if (base::is.null(hc)) {
      stop("`hc` must be provided when `module` is used.")
    }
    modules_to_run <- .hc_llm_resolve_modules(hc = hc, module = module)
    gene_infos <- lapply(modules_to_run, function(x) .hc_gemini_get_module_genes(hc = hc, module = x))
  }

  gene_infos <- gene_infos[!base::vapply(gene_infos, base::is.null, FUN.VALUE = base::logical(1))]
  if (base::length(gene_infos) == 0) {
    stop("No valid module or gene inputs were resolved.")
  }
  if (base::length(gene_infos) > 1 && !base::is.null(label) && base::nzchar(base::as.character(label[[1]]))) {
    stop("`label` can only be used for a single module or gene set.")
  }

  biological_context <- if (base::is.null(biological_context) || !base::length(biological_context)) {
    context
  } else {
    biological_context
  }
  context_text <- .hc_gemini_normalize_context(biological_context)
  results <- .hc_llm_run_sequential(
    gene_infos = gene_infos,
    label = label,
    context_text = context_text,
    llm = llm,
    api_key = api_key,
    model = model,
    vllm_base_url = vllm_base_url,
    max_genes = max_genes,
    temperature = temperature,
    timeout_sec = timeout_sec,
    pause_sec = pause_sec,
    continue_on_error = continue_on_error,
    system_instruction = system_instruction,
    response_schema = .hc_llm_response_schema(),
    verbose = verbose
  )

  if (!isTRUE(save_to_hc)) {
    if (base::length(results) == 1) {
      return(results[[1]])
    }
    return(list(
      results = stats::setNames(results, base::vapply(results, function(x) x$label, FUN.VALUE = base::character(1))),
      summary = .hc_llm_summary_from_results(results = results, hc = hc)
    ))
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  slot_obj <- sat[[slot_name]]
  if (base::is.null(slot_obj) || !base::is.list(slot_obj)) {
    slot_obj <- list()
  }
  for (res in results) {
    slot_obj[[res$label]] <- res
  }
  sat[[slot_name]] <- slot_obj
  summary_tbl <- .hc_llm_summary_from_results(results = slot_obj, hc = hc)
  sat[[base::paste0(slot_name, "_summary")]] <- summary_tbl
  hc@satellite <- S4Vectors::SimpleList(sat)
  .hc_llm_export_results_excel(
    hc = hc,
    results = slot_obj,
    summary_tbl = summary_tbl,
    slot_name = slot_name
  )
  hc
}

#' @rdname hc_module_function_llm
#' @export
hc_module_function_gemini <- function(...) {
  hc_module_function_llm(...)
}

#' @rdname hc_module_function_llm
#' @export
hc_module_function_vllm <- function(...) {
  hc_module_function_llm(llm = "vllm", ...)
}

.hc_llm_run_sequential <- function(gene_infos,
                                   label,
                                   context_text,
                                   llm,
                                   api_key,
                                   model,
                                   vllm_base_url,
                                   max_genes,
                                   temperature,
                                   timeout_sec,
                                   pause_sec,
                                   continue_on_error,
                                   system_instruction,
                                   response_schema,
                                   verbose) {
  results <- vector("list", length = base::length(gene_infos))
  for (i in base::seq_along(gene_infos)) {
    gene_info <- gene_infos[[i]]
    this_label <- if (base::is.null(label) || !base::nzchar(base::as.character(label[[1]]))) {
      gene_info$label
    } else {
      base::as.character(label[[1]])
    }
    results[[i]] <- tryCatch(
      .hc_llm_summarize_one(
        gene_info = gene_info,
        label = this_label,
          context_text = context_text,
          llm = llm,
          api_key = api_key,
          model = model,
          vllm_base_url = vllm_base_url,
          max_genes = max_genes,
        temperature = temperature,
        timeout_sec = timeout_sec,
        system_instruction = system_instruction,
        response_schema = response_schema,
        verbose = verbose
      ),
      error = function(e) {
        if (!isTRUE(continue_on_error)) {
          stop(e)
        }
        if (isTRUE(verbose)) {
          message("LLM module summary failed for `", this_label, "`: ", base::conditionMessage(e))
        }
        .hc_llm_error_result(
          gene_info = gene_info,
          label = this_label,
          context_text = context_text,
          llm = llm,
          model = model,
          error_message = base::conditionMessage(e)
        )
      }
    )
    if (pause_sec > 0 && i < base::length(gene_infos)) {
      if (isTRUE(verbose)) {
        message("Waiting ", pause_sec, " seconds before next module request.")
      }
      Sys.sleep(pause_sec)
    }
  }
  results
}

.hc_llm_summarize_many <- function(gene_infos,
                                   context_text,
                                   llm,
                                   api_key,
                                   model,
                                   vllm_base_url,
                                   max_genes,
                                   temperature,
                                   timeout_sec,
                                   system_instruction,
                                   verbose) {
  batch_inputs <- lapply(gene_infos, function(gene_info) {
    genes_all <- .hc_gemini_normalize_genes(gene_info$genes)
    genes_use <- .hc_llm_limit_genes(genes = genes_all, max_genes = max_genes)
    list(
      label = gene_info$label,
      module = gene_info$module,
      genes_all = genes_all,
      genes_use = genes_use,
      truncated = base::length(genes_use) < base::length(genes_all)
    )
  })

  if (isTRUE(verbose)) {
    message(
      "LLM module summary: sending ",
      base::length(batch_inputs),
      " modules in one combined request using provider `",
      llm,
      "` and model `",
      model,
      "`."
    )
  }

  prompt <- .hc_llm_build_batch_prompt(batch_inputs = batch_inputs, context_text = context_text)
  response_schema <- .hc_llm_batch_response_schema(n_modules = base::length(batch_inputs))

  req_result <- if (llm == "gemini") {
    .hc_llm_request_gemini(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      response_schema = response_schema,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  } else if (llm == "vllm") {
    .hc_llm_request_vllm(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      temperature = temperature,
      timeout_sec = timeout_sec,
      base_url = vllm_base_url
    )
  } else {
    .hc_llm_request_openai(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      response_schema = response_schema,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  }

  parsed_result <- tryCatch(
    jsonlite::fromJSON(req_result$result_text, simplifyVector = FALSE),
    error = function(e) {
      stop("Could not parse combined LLM response JSON: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  module_entries <- parsed_result[["modules"]]
  if (base::is.null(module_entries) || !base::is.list(module_entries) || base::length(module_entries) == 0) {
    stop("Combined LLM response did not contain a `modules` list.", call. = FALSE)
  }

  entry_labels <- base::vapply(
    module_entries,
    function(x) {
      lbl <- x[["label"]]
      if (base::is.null(lbl)) "" else base::as.character(lbl[[1]])
    },
    FUN.VALUE = base::character(1)
  )

  out <- lapply(batch_inputs, function(inp) {
    idx <- base::match(inp$label, entry_labels)
    if (base::is.na(idx)) {
      return(.hc_llm_error_result(
        gene_info = list(module = inp$module, genes = inp$genes_all),
        label = inp$label,
        context_text = context_text,
        llm = llm,
        model = model,
        error_message = "Combined LLM response did not return an entry for this module."
      ))
    }
    entry <- module_entries[[idx]]
    list(
      label = inp$label,
      module = inp$module,
      genes_input = inp$genes_all,
      genes_sent = inp$genes_use,
      gene_count_input = base::length(inp$genes_all),
      gene_count_sent = base::length(inp$genes_use),
      truncated = inp$truncated,
      context = if (base::nzchar(context_text)) context_text else NULL,
      llm = llm,
      model = model,
      status = "ok",
      error_message = NA_character_,
      prompt = prompt,
      response = entry,
      response_text = jsonlite::toJSON(entry, auto_unbox = TRUE),
      raw_response_text = req_result$raw_response_text,
      timestamp = base::as.character(Sys.time())
    )
  })

  out
}

.hc_llm_summarize_one <- function(gene_info,
                                  label,
                                  context_text,
                                  llm,
                                  api_key,
                                  model,
                                  vllm_base_url,
                                  max_genes,
                                  temperature,
                                  timeout_sec,
                                  system_instruction,
                                  response_schema,
                                  verbose) {
  genes_all <- .hc_gemini_normalize_genes(gene_info$genes)
  if (base::length(genes_all) == 0) {
    stop("No valid genes available for LLM interpretation.")
  }

  genes_use <- .hc_llm_limit_genes(genes = genes_all, max_genes = max_genes)
  truncated <- base::length(genes_use) < base::length(genes_all)

  if (isTRUE(verbose)) {
    message(
      "LLM module summary: sending ",
      base::length(genes_use),
      if (truncated) base::paste0(" of ", base::length(genes_all)) else "",
      " genes for `", label, "` using provider `", llm,
      "` and model `", model, "`."
    )
  }

  prompt <- .hc_gemini_build_prompt(
    label = label,
    genes = genes_use,
    total_gene_count = base::length(genes_all),
    context_text = context_text,
    truncated = truncated,
    llm = llm
  )

  req_result <- if (llm == "gemini") {
    .hc_llm_request_gemini(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      response_schema = response_schema,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  } else if (llm == "vllm") {
    .hc_llm_request_vllm(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      temperature = temperature,
      timeout_sec = timeout_sec,
      base_url = vllm_base_url
    )
  } else {
    .hc_llm_request_openai(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      response_schema = response_schema,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  }

  result_text <- req_result$result_text
  parsed_result <- tryCatch(
    jsonlite::fromJSON(result_text, simplifyVector = TRUE),
    error = function(e) {
      list(
        general_processes = NA_character_,
        contextual_state = NA_character_,
        key_regulators = result_text,
        parse_error = base::conditionMessage(e)
      )
    }
  )

  list(
    label = label,
    module = gene_info$module,
    genes_input = genes_all,
    genes_sent = genes_use,
    gene_count_input = base::length(genes_all),
    gene_count_sent = base::length(genes_use),
    truncated = truncated,
    context = if (base::nzchar(context_text)) context_text else NULL,
    llm = llm,
    model = model,
    status = "ok",
    error_message = NA_character_,
    prompt = prompt,
    response = parsed_result,
    response_text = result_text,
    raw_response_text = req_result$raw_response_text,
    timestamp = base::as.character(Sys.time())
  )
}

.hc_llm_limit_genes <- function(genes, max_genes) {
  if (base::is.null(max_genes) || !base::is.finite(max_genes[[1]])) {
    return(genes)
  }
  utils::head(genes, base::as.integer(max_genes[[1]]))
}

.hc_llm_error_result <- function(gene_info,
                                 label,
                                 context_text,
                                 llm,
                                 model,
                                 error_message) {
  genes_all <- .hc_gemini_normalize_genes(gene_info$genes)
  list(
    label = label,
    module = gene_info$module,
    genes_input = genes_all,
    genes_sent = base::character(0),
    gene_count_input = base::length(genes_all),
    gene_count_sent = 0L,
    truncated = FALSE,
    context = if (base::nzchar(context_text)) context_text else NULL,
    llm = llm,
    model = model,
    status = "error",
    error_message = error_message,
    prompt = NULL,
    response = list(
      general_processes = NA_character_,
      contextual_state = NA_character_,
      key_regulators = error_message
    ),
    response_text = NA_character_,
    raw_response_text = NA_character_,
    timestamp = base::as.character(Sys.time())
  )
}

.hc_llm_natural_module_order <- function(modules) {
  modules <- base::unique(base::as.character(modules))
  modules <- modules[!base::is.na(modules) & base::nzchar(modules)]
  if (base::length(modules) <= 1) {
    return(modules)
  }

  parsed <- lapply(modules, function(x) {
    hit <- regexec("^([^0-9]*?)([0-9]+)([^0-9]*)$", x, perl = TRUE)
    parts <- regmatches(x, hit)[[1]]
    if (base::length(parts) == 4) {
      return(list(
        prefix = base::tolower(parts[[2]]),
        number = suppressWarnings(base::as.integer(parts[[3]])),
        suffix = base::tolower(parts[[4]])
      ))
    }
    list(
      prefix = base::tolower(x),
      number = Inf,
      suffix = ""
    )
  })

  ord <- base::order(
    base::vapply(parsed, function(x) x$prefix, FUN.VALUE = base::character(1)),
    base::vapply(parsed, function(x) x$number, FUN.VALUE = base::numeric(1)),
    base::vapply(parsed, function(x) x$suffix, FUN.VALUE = base::character(1)),
    base::tolower(modules),
    modules
  )
  modules[ord]
}

.hc_llm_is_quota_error <- function(msg) {
  msg <- base::tolower(base::as.character(msg[[1]]))
  any(base::grepl(
    pattern = c("quota exceeded", "resource_exhausted", "rate limit", "too many requests", "retry in"),
    x = msg,
    fixed = TRUE
  ))
}

.hc_llm_resolve_modules <- function(hc, module) {
  module <- base::as.character(module)
  module <- base::trimws(module)
  module <- module[!base::is.na(module) & module != ""]
  if (base::length(module) == 0) {
    stop("No valid `module` values provided.")
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  tbl <- sat[["module_gene_list"]]
  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !"module" %in% base::colnames(tbl)) {
    stop(
      "No `hc@satellite$module_gene_list` found. Run `hc_plot_cluster_heatmap()` once first ",
      "to create the module-to-gene table."
    )
  }

  available_modules <- unique(base::as.character(tbl$module))
  available_modules <- available_modules[!base::is.na(available_modules) & available_modules != ""]
  if (base::length(available_modules) == 0) {
    stop("No modules found in `hc@satellite$module_gene_list`.")
  }

  if (base::length(module) == 1 && base::tolower(module[[1]]) == "all") {
    return(.hc_llm_module_order(hc = hc))
  }

  module
}

.hc_llm_resolve_api_key <- function(api_key, llm) {
  if (base::is.null(api_key) || !base::nzchar(base::as.character(api_key[[1]]))) {
    api_key <- if (llm == "gemini") {
      Sys.getenv("GEMINI_API_KEY", unset = "")
    } else if (llm == "vllm") {
      Sys.getenv("VLLM_API_KEY", unset = "EMPTY")
    } else {
      Sys.getenv("OPENAI_API_KEY", unset = "")
    }
  } else {
    api_key <- base::as.character(api_key[[1]])
  }

  if (!base::nzchar(api_key)) {
    if (llm == "gemini") {
      stop("No Gemini API key found. Pass `api_key` or set `GEMINI_API_KEY`.")
    } else if (llm == "vllm") {
      return("EMPTY")
    }
    stop("No OpenAI API key found. Pass `api_key` or set `OPENAI_API_KEY`.")
  }

  api_key
}

.hc_llm_resolve_vllm_base_url <- function(vllm_base_url, llm) {
  if (llm != "vllm") {
    return(NULL)
  }
  if (base::is.null(vllm_base_url) || !base::nzchar(base::as.character(vllm_base_url[[1]]))) {
    vllm_base_url <- Sys.getenv("VLLM_BASE_URL", unset = "http://localhost:8000/v1")
  } else {
    vllm_base_url <- base::as.character(vllm_base_url[[1]])
  }
  vllm_base_url <- sub("/+$", "", vllm_base_url)
  if (!base::nzchar(vllm_base_url)) {
    stop("No vLLM base URL found. Pass `vllm_base_url` or set `VLLM_BASE_URL`.")
  }
  vllm_base_url
}

.hc_llm_normalize_timeout <- function(timeout_sec) {
  if (base::is.null(timeout_sec)) {
    return(NULL)
  }
  timeout_sec <- base::as.numeric(timeout_sec[[1]])
  if (!base::is.finite(timeout_sec) || timeout_sec <= 0) {
    return(NULL)
  }
  timeout_sec
}

.hc_llm_resolve_timeout <- function(timeout_sec, llm, timeout_was_missing = FALSE) {
  if (isTRUE(timeout_was_missing) && identical(llm, "vllm")) {
    timeout_sec <- Sys.getenv("VLLM_TIMEOUT_SEC", unset = "300")
  }
  .hc_llm_normalize_timeout(timeout_sec)
}

.hc_llm_normalize_pause <- function(pause_sec) {
  if (base::is.null(pause_sec)) {
    return(0)
  }
  pause_sec <- base::as.numeric(pause_sec[[1]])
  if (!base::is.finite(pause_sec) || pause_sec <= 0) {
    return(0)
  }
  pause_sec
}

.hc_llm_require_ellmer <- function() {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    stop(
      "Package `ellmer` is required for LLM requests. Install it first.",
      call. = FALSE
    )
  }
}

.hc_llm_api_key_credentials <- function(api_key) {
  force(api_key)
  function() api_key
}

.hc_llm_with_ellmer_timeout <- function(timeout_sec, expr) {
  if (base::is.null(timeout_sec)) {
    return(force(expr))
  }

  old_timeout <- getOption("ellmer_timeout_s")
  options(ellmer_timeout_s = base::as.numeric(timeout_sec))
  on.exit(options(ellmer_timeout_s = old_timeout), add = TRUE)
  force(expr)
}

.hc_llm_gemini_temperature_unsupported_error <- function(msg) {
  msg <- base::tolower(base::as.character(msg[[1]]))
  base::grepl("unknown name \"temperature\"", msg, fixed = TRUE) ||
    base::grepl("cannot find field", msg, fixed = TRUE)
}

.hc_llm_request_gemini <- function(api_key,
                                   model,
                                   prompt,
                                   system_instruction,
                                   response_schema,
                                   temperature,
                                   timeout_sec) {
  .hc_llm_require_ellmer()

  run_request <- function(include_temperature = TRUE) {
    api_args <- list()
    if (isTRUE(include_temperature) &&
        !base::is.null(temperature) &&
        base::is.finite(temperature)) {
      api_args$temperature <- temperature
    }

    chat <- tryCatch(
      do.call(
        ellmer::chat_google_gemini,
        list(
          system_prompt = system_instruction,
          base_url = "https://generativelanguage.googleapis.com/v1beta/",
          credentials = .hc_llm_api_key_credentials(api_key),
          model = model,
          api_args = api_args,
          api_headers = c("x-goog-api-client" = "hcocena/1.9"),
          echo = "none"
        )
      ),
      error = function(e) {
        stop("Could not initialize Gemini chat via ellmer: ", base::conditionMessage(e), call. = FALSE)
      }
    )

    .hc_llm_with_ellmer_timeout(timeout_sec, chat$chat(prompt))
  }

  result_text <- tryCatch(
    run_request(include_temperature = TRUE),
    error = function(e) {
      msg <- base::conditionMessage(e)
      if (.hc_llm_gemini_temperature_unsupported_error(msg)) {
        return(
          tryCatch(
            run_request(include_temperature = FALSE),
            error = function(e2) {
              stop("Gemini request failed via ellmer: ", base::conditionMessage(e2), call. = FALSE)
            }
          )
        )
      }
      stop("Gemini request failed via ellmer: ", msg, call. = FALSE)
    }
  )

  result_text <- .hc_llm_strip_json_fences(result_text)
  if (!base::nzchar(result_text)) {
    stop("Gemini response contained no text payload.", call. = FALSE)
  }

  list(
    result_text = result_text,
    raw_response_text = result_text
  )
}

.hc_llm_request_openai <- function(api_key,
                                   model,
                                   prompt,
                                   system_instruction,
                                   response_schema,
                                   temperature,
                                   timeout_sec) {
  .hc_llm_require_ellmer()

  chat <- tryCatch(
    ellmer::chat_openai(
      system_prompt = system_instruction,
      base_url = "https://api.openai.com/v1",
      credentials = .hc_llm_api_key_credentials(api_key),
      model = model,
      api_args = list(
        temperature = temperature
      ),
      echo = "none"
    ),
    error = function(e) {
      stop("Could not initialize OpenAI chat via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- tryCatch(
    .hc_llm_with_ellmer_timeout(timeout_sec, chat$chat(prompt)),
    error = function(e) {
      stop("OpenAI request failed via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- .hc_llm_strip_json_fences(result_text)
  if (!base::nzchar(result_text)) {
    stop("OpenAI response contained no text payload.", call. = FALSE)
  }

  list(
    result_text = result_text,
    raw_response_text = result_text
  )
}

.hc_llm_request_vllm <- function(api_key,
                                 model,
                                 prompt,
                                 system_instruction,
                                 temperature,
                                 timeout_sec,
                                 base_url) {
  .hc_llm_require_ellmer()

  chat <- tryCatch(
    ellmer::chat_vllm(
      base_url = base_url,
      model = model,
      credentials = .hc_llm_api_key_credentials(api_key),
      system_prompt = system_instruction,
      api_args = list(
        temperature = temperature,
        max_tokens = 8000,
        chat_template_kwargs = list(enable_thinking = FALSE)
      )
    ),
    error = function(e) {
      stop("Could not initialize vLLM chat via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- tryCatch(
    .hc_llm_with_ellmer_timeout(timeout_sec, chat$chat(prompt)),
    error = function(e) {
      stop("vLLM request failed via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- .hc_llm_strip_think_blocks(result_text)
  result_text <- .hc_llm_strip_json_fences(result_text)
  if (!base::nzchar(result_text)) {
    stop("vLLM response contained no text payload.", call. = FALSE)
  }

  list(
    result_text = result_text,
    raw_response_text = result_text
  )
}

.hc_gemini_normalize_genes <- function(genes) {
  genes <- base::as.character(genes)
  genes <- base::trimws(genes)
  genes <- genes[!base::is.na(genes) & genes != ""]
  genes <- genes[!duplicated(genes)]
  genes
}

.hc_gemini_normalize_context <- function(context) {
  if (base::is.null(context)) {
    return("")
  }
  context <- base::as.character(context)
  context <- base::trimws(context)
  context <- context[!base::is.na(context) & context != ""]
  if (base::length(context) == 0) {
    return("")
  }
  base::paste(context, collapse = "\n")
}

.hc_gemini_get_module_genes <- function(hc, module) {
  module <- base::as.character(module[[1]])
  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  tbl <- sat[["module_gene_list"]]

  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !"genes" %in% base::colnames(tbl) || !"module" %in% base::colnames(tbl)) {
    stop(
      "No `hc@satellite$module_gene_list` found. Run `hc_plot_cluster_heatmap()` once first ",
      "to create the module-to-gene table."
    )
  }

  tbl$genes <- base::as.character(tbl$genes)
  tbl$module <- base::as.character(tbl$module)

  module_lookup <- module
  if (!(module_lookup %in% tbl$module)) {
    label_map <- tryCatch(hc@integration@cluster[["module_label_map"]], error = function(e) NULL)
    if (!base::is.null(label_map) && base::length(label_map) > 0) {
      label_map <- base::as.character(label_map)
      map_names <- base::names(label_map)
      if (!base::is.null(map_names) && base::length(map_names) == base::length(label_map)) {
        base::names(label_map) <- base::as.character(map_names)
        if (module %in% base::names(label_map)) {
          module_lookup <- base::as.character(label_map[[module]])
        }
      }
    }
  }

  genes <- tbl$genes[tbl$module == module_lookup]
  genes <- .hc_gemini_normalize_genes(genes)
  if (base::length(genes) == 0) {
    stop("Could not resolve genes for module `", module, "`.")
  }

  list(
    label = module_lookup,
    module = module_lookup,
    genes = genes
  )
}

.hc_gemini_build_prompt <- function(label,
                                    genes,
                                    total_gene_count,
                                    context_text,
                                    truncated,
                                    llm = "gemini") {
  trunc_note <- if (isTRUE(truncated)) {
    base::paste0(
      "Only the first ", base::length(genes),
      " genes are sent here out of ", total_gene_count,
      " total input genes due to `max_genes`."
    )
  } else {
    base::paste0("All ", total_gene_count, " input genes are included.")
  }

  biological_context <- if (base::nzchar(context_text)) context_text else "none provided"

  prompt_instructions <- if (llm == "vllm") {
    base::paste(
      "Infer the main biological program of this transcriptomic module.",
      base::paste0("The biological context is: ", biological_context, "."),
      "Return JSON matching the provided schema.",
      "Provide exactly three compact but informative text fields.",
      "Do not start with phrases like this module or the genes.",
      "Use the biological context as framing information, but do not simply repeat it.",
      "Make the answer module-specific and distinguishable from other modules from the same study.",
      "Avoid one-word or overly generic labels if the genes support something more specific.",
      "For `general_processes`, list 2 to 4 specific biological themes separated by ' / '.",
      "For `contextual_state`, give one compact but informative phrase of about 4 to 10 words describing the immediate module state.",
      "Do not use `contextual_state` to restate the cohort or timeline.",
      "Avoid generic phrases such as immature monocyte maturation, preterm infant maturation, or developing monocytes unless the module is truly nonspecific.",
      "Prefer specific states such as interferon-high inflammatory state, ribosome-high proliferative state, phagolysosomal activated state, antigen-presenting inflammatory state, platelet-like metabolic state, erythroid-skewed progenitor-like state, or macrophage-like transition state when supported.",
      "For `key_regulators`, list 2 to 5 likely transcription factors or signaling regulators separated by ' / '.",
      "If regulator evidence is weak, provide the most plausible regulators briefly rather than repeating the process.",
      "Example style only:",
      '{"general_processes":"interferon signaling / antiviral innate immunity / antigen presentation","contextual_state":"activated interferon-high inflammatory monocyte state","key_regulators":"STAT1 / IRF7 / IRF9 / NFKB1"}'
    )
  } else {
    base::paste(
      "Infer the main biological program of this transcriptomic module.",
      base::paste0("The biological context is: ", biological_context, "."),
      "Return JSON matching the provided schema.",
      "Be extremely concise and provide exactly three short text elements.",
      "Do not start with phrases like this module or the genes.",
      "Treat the biological context as framing information, not as wording to repeat.",
      "Do not simply restate broad phrases from the context such as monocyte maturation unless the module is truly nonspecific.",
      "Make the answer module-specific and distinguishable from other modules from the same study.",
      "If the genes support a more specific program such as interferon response, phagolysosome, antigen presentation, cell cycle, platelet-like program, erythroid bias, mitochondrial metabolism, ribosome biogenesis, glycolysis, chemotaxis, inflammatory signaling, or tissue contamination, prefer that over generic immune wording.",
      "Provide `general_processes` as a short noun phrase listing the main biological program.",
      "Provide `contextual_state` as a short phrase describing the specific monocyte or transcriptional state in this study context.",
      "Provide `key_regulators` as a short phrase naming likely driving transcription factors or signaling regulators.",
      "If regulator evidence is weak, state the most plausible regulators briefly rather than repeating the biological process."
    )
  }

  base::paste(
    prompt_instructions,
    "",
    base::paste0("Label: ", label),
    base::paste0("Gene-count note: ", trunc_note),
    base::paste0("Genes:\n", base::paste(genes, collapse = ", ")),
    sep = "\n"
  )
}

.hc_llm_build_batch_prompt <- function(batch_inputs, context_text) {
  context_block <- if (base::nzchar(context_text)) {
    base::paste0("Context:\n", context_text, "\n")
  } else {
    "Context:\nNone provided.\n"
  }

  module_blocks <- base::vapply(
    batch_inputs,
    function(inp) {
      trunc_note <- if (isTRUE(inp$truncated)) {
        base::paste0(
          "Only the first ", base::length(inp$genes_use),
          " genes are listed here out of ", base::length(inp$genes_all),
          " total input genes due to `max_genes`."
        )
      } else {
        base::paste0("All ", base::length(inp$genes_all), " input genes are included.")
      }
      base::paste(
        base::paste0("Module label: ", inp$label),
        base::paste0("Gene-count note: ", trunc_note),
        base::paste0("Genes:\n", base::paste(inp$genes_use, collapse = ", ")),
        sep = "\n"
      )
    },
    FUN.VALUE = base::character(1)
  )

  base::paste(
    "Summarize the likely overarching biological function of each transcriptomic module separately.",
    "Return JSON matching the provided schema with one entry per module.",
    "For each module, provide `short_title` as a short plot-ready label with about 3 to 8 words.",
    "The provided context is the primary interpretation frame and must strongly constrain the answer.",
    "Prefer explanations that are compatible with the given biological context, cell type, cohort, and tissue.",
    "Do not assign unrelated tissue programs such as neuronal, epithelial, ciliary, muscular, or organ-specific identities unless the evidence is overwhelming and no context-compatible explanation fits.",
    "If genes look context-mismatched, keep the interpretation context-aware and mention uncertainty in `caveats` instead of drifting to an unrelated lineage.",
    "This is an interpretation task, not a statistical enrichment test.",
    "",
    context_block,
    base::paste(module_blocks, collapse = "\n\n---\n\n"),
    sep = "\n"
  )
}

.hc_llm_response_schema <- function() {
  json_schema <- list(
    type = "object",
    properties = list(
      general_processes = list(
        type = "string",
        description = "Two to four specific biological themes separated by ' / '."
      ),
      contextual_state = list(
        type = "string",
        description = "A compact but informative phrase of about four to ten words describing the specific transcriptional or cellular state."
      ),
      key_regulators = list(
        type = "string",
        description = "Two to five likely transcription factors or signaling regulators separated by ' / '."
      )
    ),
    required = c("general_processes", "contextual_state", "key_regulators")
  )
  json_schema
}

.hc_llm_batch_response_schema <- function(n_modules) {
  list(
    type = "object",
    properties = list(
      modules = list(
        type = "array",
        minItems = base::as.integer(n_modules),
        maxItems = base::as.integer(n_modules),
        items = .hc_llm_response_schema()
      )
    ),
    required = base::as.list("modules"),
    additionalProperties = FALSE
  )
}

.hc_gemini_extract_response_text <- function(resp_obj) {
  candidates <- resp_obj[["candidates"]]
  if (base::is.null(candidates) || base::length(candidates) == 0) {
    return("")
  }
  first_candidate <- candidates[[1]]
  content <- first_candidate[["content"]]
  if (base::is.null(content)) {
    return("")
  }
  parts <- content[["parts"]]
  if (base::is.null(parts) || base::length(parts) == 0) {
    return("")
  }
  text_parts <- base::vapply(
    parts,
    function(x) {
      txt <- x[["text"]]
      if (base::is.null(txt)) "" else base::as.character(txt)
    },
    FUN.VALUE = base::character(1)
  )
  text_parts <- text_parts[text_parts != ""]
  if (base::length(text_parts) == 0) {
    return("")
  }
  base::paste(text_parts, collapse = "\n")
}

.hc_openai_extract_response_text <- function(resp_obj) {
  message_obj <- tryCatch(resp_obj$choices[[1]]$message, error = function(e) NULL)
  if (base::is.null(message_obj)) {
    return("")
  }

  content <- message_obj$content
  if (base::is.null(content)) {
    return("")
  }

  if (base::is.character(content) && base::length(content) >= 1) {
    return(base::as.character(content[[1]]))
  }

  if (base::is.list(content) && base::length(content) > 0) {
    text_parts <- base::character(0)
    for (part in content) {
      part_text <- part$text
      if (base::is.character(part_text) && base::length(part_text) >= 1) {
        text_parts <- c(text_parts, base::as.character(part_text[[1]]))
      } else if (base::is.list(part_text) && !base::is.null(part_text$value)) {
        text_parts <- c(text_parts, base::as.character(part_text$value[[1]]))
      }
    }
    text_parts <- text_parts[!base::is.na(text_parts) & text_parts != ""]
    if (base::length(text_parts) > 0) {
      return(base::paste(text_parts, collapse = "\n"))
    }
  }

  ""
}

.hc_llm_strip_json_fences <- function(x) {
  x <- base::as.character(x[[1]])
  if (!base::nzchar(x)) {
    return("")
  }
  x <- sub("^\\s*```(?:json)?\\s*", "", x, perl = TRUE)
  x <- sub("\\s*```\\s*$", "", x, perl = TRUE)
  stringr::str_squish(x)
}

.hc_llm_strip_think_blocks <- function(x) {
  x <- base::as.character(x[[1]])
  if (!base::nzchar(x)) {
    return("")
  }
  # Local Qwen/vLLM deployments may emit reasoning in <think>...</think>.
  x <- gsub("(?is)<think>.*?</think>", " ", x, perl = TRUE)
  x <- sub("(?is)<think>.*$", " ", x, perl = TRUE)
  stringr::str_squish(x)
}

.hc_llm_summary_from_results <- function(results, hc = NULL) {
  if (base::is.null(results) || base::length(results) == 0) {
    return(base::data.frame(
      module = base::character(0),
      module_color = base::character(0),
      llm = base::character(0),
      model = base::character(0),
      general_processes = base::character(0),
      contextual_state = base::character(0),
      key_regulators = base::character(0),
      llm_long_output = base::character(0),
      response_json = base::character(0),
      short_title = base::character(0),
      overarching_function = base::character(0),
      confidence = base::character(0),
      gene_count_input = base::integer(0),
      gene_count_sent = base::integer(0),
      truncated = base::logical(0),
      status = base::character(0),
      error_message = base::character(0),
      timestamp = base::character(0),
      stringsAsFactors = FALSE
    ))
  }

  results <- results[base::vapply(results, base::is.list, FUN.VALUE = base::logical(1))]
  if (base::length(results) == 0) {
    return(base::data.frame(
      module = base::character(0),
      module_color = base::character(0),
      llm = base::character(0),
      model = base::character(0),
      general_processes = base::character(0),
      contextual_state = base::character(0),
      key_regulators = base::character(0),
      llm_long_output = base::character(0),
      response_json = base::character(0),
      short_title = base::character(0),
      overarching_function = base::character(0),
      confidence = base::character(0),
      gene_count_input = base::integer(0),
      gene_count_sent = base::integer(0),
      truncated = base::logical(0),
      status = base::character(0),
      error_message = base::character(0),
      timestamp = base::character(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- lapply(results, function(res) {
    base::data.frame(
      module = .hc_llm_result_scalar(res, c("label"), default = .hc_llm_result_scalar(res, c("module"), default = NA_character_)),
      module_color = .hc_llm_module_color(label = .hc_llm_result_label(res), hc = hc),
      llm = .hc_llm_result_scalar(res, c("llm")),
      model = .hc_llm_result_scalar(res, c("model")),
      general_processes = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "general_processes"))),
      contextual_state = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "contextual_state"))),
      key_regulators = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "key_regulators"))),
      llm_long_output = .hc_llm_result_long_output(res),
      response_json = .hc_llm_result_json(res),
      short_title = .hc_llm_result_short_title(res),
      overarching_function = .hc_llm_result_overarching(res),
      confidence = NA_character_,
      gene_count_input = .hc_llm_result_scalar(res, c("gene_count_input"), default = NA_integer_, mode = "integer"),
      gene_count_sent = .hc_llm_result_scalar(res, c("gene_count_sent"), default = NA_integer_, mode = "integer"),
      truncated = .hc_llm_result_scalar(res, c("truncated"), default = FALSE, mode = "logical"),
      status = .hc_llm_result_scalar(res, c("status")),
      error_message = .hc_llm_result_scalar(res, c("error_message")),
      timestamp = .hc_llm_result_scalar(res, c("timestamp")),
      stringsAsFactors = FALSE
    )
  })
  out <- base::do.call(base::rbind, out)
  base::rownames(out) <- NULL

  if (!base::is.null(hc)) {
    module_order <- .hc_llm_module_order(hc = hc)
    if (base::length(module_order) > 0) {
      keep <- module_order[module_order %in% out$module]
      extras <- base::setdiff(out$module, keep)
      ord <- base::match(base::c(keep, extras), out$module)
      ord <- ord[!base::is.na(ord)]
      out <- out[ord, , drop = FALSE]
      base::rownames(out) <- NULL
    }
  }

  out
}

.hc_llm_result_label <- function(res) {
  lbl <- .hc_llm_result_scalar(res, c("label"), default = NA_character_)
  if (base::is.na(lbl) || !base::nzchar(lbl)) {
    lbl <- .hc_llm_result_scalar(res, c("module"), default = NA_character_)
  }
  lbl
}

.hc_llm_result_short_title <- function(res) {
  short_title <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "contextual_state")))
  if (base::is.null(short_title) || base::length(short_title) == 0) {
    short_title <- NA_character_
  } else {
    short_title <- base::as.character(short_title[[1]])
  }
  if (!base::nzchar(short_title) || base::is.na(short_title)) {
    short_title <- .hc_llm_short_title_fallback(
      .hc_llm_result_field(res, c("response", "general_processes"))
    )
  }
  short_title
}

.hc_llm_result_overarching <- function(res) {
  val <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "general_processes")))
  if (base::is.null(val) || base::length(val) == 0) {
    return(NA_character_)
  }
  base::as.character(val[[1]])
}

.hc_llm_result_long_output <- function(res) {
  gp <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "general_processes")))
  cs <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "contextual_state")))
  kr <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "key_regulators")))
  parts <- c(
    if (!base::is.null(gp) && base::nzchar(base::as.character(gp[[1]]))) base::paste0("General processes: ", base::as.character(gp[[1]])) else NULL,
    if (!base::is.null(cs) && base::nzchar(base::as.character(cs[[1]]))) base::paste0("Contextual state: ", base::as.character(cs[[1]])) else NULL,
    if (!base::is.null(kr) && base::nzchar(base::as.character(kr[[1]]))) base::paste0("Key regulators: ", base::as.character(kr[[1]])) else NULL
  )
  if (base::length(parts) == 0) {
    return(NA_character_)
  }
  base::paste(parts, collapse = " | ")
}

.hc_llm_result_json <- function(res) {
  resp <- .hc_llm_result_field(res, c("response"))
  if (base::is.null(resp)) {
    return(NA_character_)
  }
  out <- tryCatch(
    jsonlite::toJSON(resp, auto_unbox = TRUE, null = "null"),
    error = function(e) NA_character_
  )
  base::as.character(out[[1]])
}

.hc_llm_result_field <- function(x, path) {
  cur <- x
  for (nm in path) {
    if (base::is.null(cur) || !(nm %in% base::names(cur))) {
      return(NULL)
    }
    cur <- cur[[nm]]
  }
  cur
}

.hc_llm_result_scalar <- function(res, path, default = NA_character_, mode = c("character", "integer", "logical")) {
  mode <- base::match.arg(mode)
  val <- .hc_llm_result_field(res, path)
  if (base::is.null(val) || base::length(val) == 0) {
    return(default)
  }

  if (identical(mode, "character")) {
    val <- base::as.character(val[[1]])
    if (base::length(val) == 0) {
      return(base::as.character(default[[1]]))
    }
    return(val)
  }
  if (identical(mode, "integer")) {
    val <- suppressWarnings(base::as.integer(val[[1]]))
    if (base::length(val) == 0 || base::is.na(val)) {
      return(suppressWarnings(base::as.integer(default[[1]])))
    }
    return(val)
  }

  val <- base::as.logical(val[[1]])
  if (base::length(val) == 0 || base::is.na(val)) {
    return(base::as.logical(default[[1]]))
  }
  val
}

.hc_llm_short_title_fallback <- function(term, max_chars = 64) {
  if (base::is.null(term) || base::length(term) == 0) {
    return("No interpretation available")
  }
  term <- .hc_llm_clean_text(base::as.character(term[[1]]))
  if (!base::nzchar(term) || base::is.na(term)) {
    return("No interpretation available")
  }

  term <- stringr::str_squish(term)
  term <- sub("^This module is best described as\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module can be summarized as\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module primarily reflects\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module reflects\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module primarily captures\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module captures\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module is characterized by\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module represents\\s+", "", term, ignore.case = TRUE)
  term <- sub("^processes related to\\s+", "", term, ignore.case = TRUE)
  term <- sub("^processes involving\\s+", "", term, ignore.case = TRUE)
  term <- sub("^the coordinated regulation of\\s+", "", term, ignore.case = TRUE)
  term <- sub("^the orchestration of\\s+", "", term, ignore.case = TRUE)
  term <- sub("^the activation of\\s+", "", term, ignore.case = TRUE)
  term <- sub("^the diverse functions of\\s+", "", term, ignore.case = TRUE)
  term <- sub("\\.$", "", term)

  pieces <- unlist(strsplit(
    term,
    "\\s*,\\s*|\\s*;\\s*|\\s+and\\s+|\\s+with\\s+|\\s+linked to\\s+|\\s+coupled to\\s+",
    perl = TRUE
  ))
  pieces <- trimws(pieces)
  pieces <- pieces[nzchar(pieces)]
  if (length(pieces) == 0) {
    pieces <- term
  }

  chosen <- character(0)
  for (piece in pieces) {
    candidate <- paste(c(chosen, piece), collapse = " / ")
    if (nchar(candidate) > max_chars && length(chosen) > 0) {
      break
    }
    chosen <- c(chosen, piece)
    if (nchar(candidate) >= (max_chars - 8)) {
      break
    }
  }

  if (length(chosen) == 0) {
    chosen <- pieces[[1]]
  }

  short_title <- paste(chosen, collapse = " / ")
  short_title <- .hc_llm_clean_text(stringr::str_squish(short_title))
  stringr::str_trunc(short_title, width = max_chars)
}

.hc_llm_clean_text <- function(x) {
  if (base::is.null(x) || base::length(x) == 0) {
    return(x)
  }
  x <- base::as.character(x)
  x <- stringr::str_squish(x)
  x <- gsub("\\s*/\\s*/+\\s*", " / ", x, perl = TRUE)
  x <- gsub("\\s*\\|\\s*\\|+\\s*", " | ", x, perl = TRUE)
  x <- gsub("\\s{2,}", " ", x, perl = TRUE)
  x <- base::trimws(x)
  x
}

.hc_llm_export_results_excel <- function(hc,
                                         results,
                                         summary_tbl,
                                         slot_name) {
  if (base::is.null(hc) || !requireNamespace("openxlsx", quietly = TRUE)) {
    return(invisible(NULL))
  }

  details_tbl <- base::do.call(
    base::rbind,
    lapply(results, function(res) {
      base::data.frame(
        module = .hc_llm_result_label(res),
        llm = .hc_llm_result_scalar(res, c("llm")),
        model = .hc_llm_result_scalar(res, c("model")),
        general_processes = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "general_processes"))),
        contextual_state = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "contextual_state"))),
        key_regulators = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "key_regulators"))),
        llm_long_output = .hc_llm_result_long_output(res),
        response_json = .hc_llm_result_json(res),
        gene_count_input = .hc_llm_result_scalar(res, c("gene_count_input"), default = NA_integer_, mode = "integer"),
        gene_count_sent = .hc_llm_result_scalar(res, c("gene_count_sent"), default = NA_integer_, mode = "integer"),
        truncated = .hc_llm_result_scalar(res, c("truncated"), default = FALSE, mode = "logical"),
        status = .hc_llm_result_scalar(res, c("status")),
        error_message = .hc_llm_result_scalar(res, c("error_message")),
        prompt = .hc_llm_result_scalar(res, c("prompt")),
        timestamp = .hc_llm_result_scalar(res, c("timestamp")),
        stringsAsFactors = FALSE
      )
    })
  )
  base::rownames(details_tbl) <- NULL

  out_dir <- .hc_resolve_output_dir(hc)
  file <- base::file.path(out_dir, base::paste0(slot_name, "_summary.xlsx"))
  tryCatch(
    openxlsx::write.xlsx(
      x = list(
        summary = summary_tbl,
        details = details_tbl
      ),
      file = file,
      overwrite = TRUE
    ),
    error = function(e) warning("Could not write LLM Excel summary: ", base::conditionMessage(e))
  )
  invisible(file)
}

.hc_llm_module_order <- function(hc) {
  heatmap_info <- tryCatch(.hc_heatmap_cache_info(hc@integration@cluster), error = function(e) NULL)
  label_map <- tryCatch(hc@integration@cluster[["module_label_map"]], error = function(e) NULL)
  if (!base::is.null(heatmap_info) &&
      !base::is.null(heatmap_info$row_order) &&
      base::length(heatmap_info$row_order) > 0 &&
      !base::is.null(label_map) &&
      base::length(label_map) > 0) {
    label_map <- base::as.character(label_map)
    map_names <- base::names(hc@integration@cluster[["module_label_map"]])
    if (!base::is.null(map_names) && base::length(map_names) == base::length(label_map)) {
      base::names(label_map) <- base::as.character(map_names)
    }
    mapped <- base::as.character(label_map[base::as.character(heatmap_info$row_order)])
    mapped <- mapped[!base::is.na(mapped) & base::nzchar(mapped)]
    if (base::length(mapped) > 0) {
      return(base::unique(mapped))
    }
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  tbl <- sat[["module_gene_list"]]
  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !"module" %in% base::colnames(tbl)) {
    return(base::character(0))
  }
  .hc_llm_natural_module_order(tbl$module)
}

.hc_llm_module_color <- function(label, hc = NULL) {
  if (base::is.null(hc)) {
    return("grey70")
  }
  label <- base::as.character(label[[1]])
  label_map <- tryCatch(hc@integration@cluster[["module_label_map"]], error = function(e) NULL)
  if (base::is.null(label_map) || base::length(label_map) == 0) {
    return("grey70")
  }
  label_map <- base::as.character(label_map)
  map_names <- base::names(hc@integration@cluster[["module_label_map"]])
  if (!base::is.null(map_names) && base::length(map_names) == base::length(label_map)) {
    base::names(label_map) <- base::as.character(map_names)
  }
  inv_map <- stats::setNames(base::names(label_map), label_map)
  col <- inv_map[[label]]
  if (base::is.null(col) || !base::nzchar(base::as.character(col))) {
    return("grey70")
  }
  base::as.character(col)
}
