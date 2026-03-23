## GitHub Workflows

This branch keeps extended hCoCena workflow documents available on GitHub
without adding them to the Bioconductor submission branch.

The files in [`github_workflows/`](github_workflows/) keep the familiar
`main` / `satellite` workflow structure while using the current object-based
API:

- `hcocena_main.Rmd`
- `hcocena_satellite.Rmd`
- `hcocena_main_seq_only.Rmd`

These documents are useful as:

- a reference for the GitHub-only walkthrough structure
- a starting point for GitHub-only walkthroughs
- longer examples that stay outside the Bioconductor submission branch

They are not treated as package vignettes and are intentionally kept outside
the Bioconductor submission flow.

See [`github_workflows/workflow_notes.md`](github_workflows/workflow_notes.md)
for a short map of the current workflow set on this branch.
