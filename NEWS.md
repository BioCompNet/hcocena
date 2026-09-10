# hcocena 0.99.7

## Examples and fixtures

- Rebuilt the bundled example fixtures with a 16-donor, two-timepoint
  longitudinal design, so the meta-clustering functions can be demonstrated on
  real output rather than only described. `inst/scripts/make-fixtures.R`
  documents how they are generated.
- Dropped two stored-but-never-read plot objects (`dd_plot_calculated_optimal`
  and the layer heatmap) from the fixtures. Together with the new design this
  takes the source tarball from 2.8 MB to 1.2 MB.
- Added runnable examples to 17 further help pages, covering the longitudinal
  meta-clustering chain, the direct workflow, `hc_meta_correlation_num()`, the
  cell-type database helpers, and the `hc_sample_regrouping` plot method.

## Fixes

- `hc_longitudinal_workflow_direct()` failed immediately with `'arg' must be of
  length 1`. It forwards `cap_na_impute` explicitly, but the receiving formal
  defaults to the already-matched `na_impute`, so `match.arg()` saw a single
  choice. The three affected call sites now pass `choices` explicitly.
- `hc_meta_correlation_cat()` correlates across groups of the variable of
  interest, so it needs at least three of them. With two it fell through to a
  cryptic `cor.test()` error; it now reports the requirement and the groups it
  found.

## Export reliability

- Build and validate XLSX workbooks on R's local temporary filesystem before
  publishing them to synchronized or bind-mounted output directories.
- Verify staged XLSX transfers byte-for-byte and atomically replace existing
  outputs without exposing partially written workbooks.
- Route both table-based exports and updated `Hub_genes.xlsx` workbooks through
  the same local-staging path.

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
