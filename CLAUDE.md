# CLAUDE.md

Guidance for Claude Code / other agents working in this repository.

## Project

`hcocena` — Bioconductor-style R package (v0.99.6) for horizontal integration and
downstream analysis of transcriptomics datasets. Modern S4 workflow around
`HCoCenaExperiment` (built on `MultiAssayExperiment`/`SummarizedExperiment`).
Package source lives at the repository root. Docker setup is in `docker/`,
GitHub walkthrough notebooks in `github_workflows/`.

## Running R (important — read before running anything)

- Use the **native Windows Rscript**, invoked from PowerShell:
  `& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' ...`
- Do **not** run R through MSYS / Git-bash. The bash R wrapper segfaults on the
  htmlwidget viewer code path.
- Load the package in dev mode: `pkgload::load_all(".")`
  (`pkgload` + `devtools` are installed; `hcocena` is also installed).
- Regenerate docs after roxygen changes: `devtools::document()`.

## Testing

- `testthat` edition 3 under `tests/testthat/`. Run all: `Rscript -e "devtools::test()"`.
- Real-data regression is **opt-in** (STAR-protocol data, see README
  "Real-data regression checks"):
  `HCOCENA_RUN_REALDATA=true HCOCENA_REALDATA_MODE=quick Rscript -e "testthat::test_file('tests/testthat/test-realdata-regression.R')"`
  or `Rscript scripts/run_realdata_regression.R --mode quick`.

## LLM-assisted module interpretation (active area of work)

This is an AI interpretation step, **not** a statistical enrichment test.

Files:
- `R/hc_module_function_gemini.R` — core `hc_module_function_llm()`; thin wrappers
  `hc_module_function_gemini(...)` and `hc_module_function_vllm(...)`.
- `R/hc_plot_module_function_gemini.R` — plotting of LLM results.
- `R/hc_list_llm_models.R` — `hc_list_llm_models()` lists available models per
  provider (newer file; currently untracked).

Key facts:
- **Providers** via `llm=` (or alias `provider=`): `"gemini"`, `"claude"`/`"anthropic"`,
  `"openai"`/`"chatgpt"`, `"vllm"` (local OpenAI-compatible). Requests go through
  the `ellmer` package (in Suggests).
- **API keys** from env vars: `GEMINI_API_KEY`, `ANTHROPIC_API_KEY`,
  `OPENAI_API_KEY`, `VLLM_API_KEY`. Local server URL via `base_url=` /
  `vllm_base_url=` / `VLLM_BASE_URL`.
- **Default models**: `gemini-2.5-pro`, `claude-sonnet-4-6`, `gpt-4o-mini`,
  `Qwen/Qwen2.5-VL-32B-Instruct`. Override with `claude_model=`, `model=`, etc.
- **Testing with a custom gene list** (no `hc` object needed): pass
  `genes = c(...)` together with `save_to_hc = FALSE`. Use `module=` instead for
  real modules (needs an `hc` with `module_gene_list`, created by
  `hc_plot_cluster_heatmap()`). `module` and `genes` are mutually exclusive.
- Result object: `res$response` holds the parsed JSON
  (`general_processes` / `contextual_state` / `key_regulators`); also
  `res$model`, `res$genes_sent`, `res$prompt`, `res$raw_response_text`.

History note (so agents don't re-derive it): the `genes=` custom-list capability
existed since the first LLM release (v1.9). **Claude/Anthropic support** and the
`provider=` / `base_url=` aliases plus `hc_list_llm_models()` were added later.
In v1.9 `llm` only accepted gemini/openai/chatgpt/vllm.

Minimal Claude smoke test:
```r
pkgload::load_all(".")
res <- hc_module_function_llm(
  genes      = c("STAT1", "IRF7", "CXCL10", "GBP1", "IFI44L"),
  context    = "Interferon-driven blood module",   # optional
  llm        = "claude",
  api_key    = Sys.getenv("ANTHROPIC_API_KEY"),
  save_to_hc = FALSE
)
res$response
```

A ready-to-edit script lives at repo root: **`test_llm_claude.R`** (set
`ANTHROPIC_API_KEY`, edit `my_genes`, run with native Rscript). It is a local
scratch script, not part of the package.

## Other recent context

- **htmlwidgets**: in RStudio notebooks widgets render via the viewer option;
  `repr`/`cat` output is Jupyter-only and opt-in (`R/hc_display_helpers.R`).
  Several recent commits route htmlwidgets through the print generic for inline
  notebook rendering.
- **Heatmaps**: module-label auto-fit to box size and split-label ordering were
  recently reworked (`R/hc_plot_cluster_heatmap.R`,
  `R/hc_heatmap_export_helpers.R`).

## Conventions

- Functions use explicit `base::` / `pkg::` namespace prefixes throughout.
- Roxygen2 with `markdown = TRUE`.
- The public API is **`hc_*` only**. The legacy entry points (`TF_overrep_module`,
  `find_hubs`, `read_data`, …) were removed; the old names still exist as
  unexported implementations that the `hc_*` wrappers call through
  `.hc_run_driver()`, but they are not part of the API and must not gain
  `@export` tags again.
- Every exported name must start with `hc_` — `tests/testthat/test-api-config.R`
  enforces this. If `devtools::document()` ever adds non-`hc_` exports to
  NAMESPACE, an internal function has picked up a stray `#' @export`; remove the
  tag rather than editing NAMESPACE by hand.
- S3 methods (`print.hc_llm_heatmap_plot`, `plot.hc_sample_regrouping`) are the
  one exception: roxygen registers them *via* `@export`, so those tags must stay.
