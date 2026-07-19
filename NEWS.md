# hcocena 0.99.6

## Enrichment defaults

- Made module-heatmap column gaps opt-in via `smart_column_gaps`, with
  `column_gap_by` for explicit metadata-based splits and `column_gap_mm` for
  gap size control.
- Switched functional-enrichment defaults to consistent term selection across
  modules and wrappers.
- Added optional DoRAG retrieval support to `hc_module_function_llm()` so LLM
  module summaries can be grounded in retrieved passages and stored citations.
- RAG runs now preserve the normal context-aware interpretation in `response`
  and store a separate literature-supported interpretation in `rag_response`.
- Added `rag_connect_timeout_sec` and `rag_continue_on_error` to make DoRAG
  retrieval robust to unreachable or slow RAG servers.
- Extended `hc_plot_module_function_llm()` to plot separate RAG fields such as
  `rag_contextual_state` (also available as `contextual_state_rag`).
- Added common RNA-seq differential-expression packages to the Docker image,
  including `DESeq2`, `limma`, `sva`, `edgeR`, and supporting visualization,
  shrinkage, import, and organism annotation packages.
- Allowed `hc_split_modules()` to use one Leiden `resolution` value per
  selected module.
- Clarified that `hc_read_data()` now removes zero-variance genes and drops
  non-numeric helper columns from object-based count inputs, which can shift
  `hc_suggest_topvar()` inflection points slightly compared with older
  releases.
- Refreshed the Docker release metadata for the next public image tag.
- Rendered interactive cutoff plots inline during HTML/R Markdown knitting
  instead of opening them in the RStudio Viewer.

# hcocena 0.99.5

## Bioconductor readiness and Docker release

- Finalized the S4/legacy bridge cleanup, including removal of remaining
  package-level `<<-` usage from the active R sources.
- Hardened regression coverage for the updated heatmap and auto-tuning paths.
- Refined package formatting and documentation metadata ahead of submission.
- Refreshed the public Docker release metadata for the next image tag.

# hcocena 0.99.1

## Bioconductor preparation

- Aligned the package version with Bioconductor pre-submission conventions.
- Simplified `DESCRIPTION` metadata for the first Bioconductor submission.
- Added a package-level `README.md` and `inst/CITATION`.
- Expanded the workflow and migration vignettes to use `BiocStyle` and
  reproducible examples based on `inst/extdata`.
- Added ignore rules for local build artifacts and reduced the default branch to
  package source for submission.
- Added compatibility fixes for longitudinal `rfcont` imputation and clarified
  that this workflow requires `library(CALIBERrfimpute)` in the active session.
- Made LLM-related examples safe for package checks and non-interactive builds.
