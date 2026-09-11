## Workflow Notes

These R Markdown files are kept on the `workflows` branch as GitHub-only
workflow templates. They are not package vignettes.

### Best starting points

- `hcocena_main.Rmd`: closest to the current object-passing workflow
- `hcocena_satellite.Rmd`: useful for optional analysis ideas

### Other files in git history

The Git history still contains older STAR protocol and earlier workflow
variants, but they are intentionally not part of this GitHub workflow set.

### Workflow direction

When updating these workflow files, prefer:

- explicit object-passing with `hc <- ...`
- current exported package functions
- current `hc@satellite` / `satellite_outputs` storage conventions
- examples that match the current S4/container workflow

### Wrapper status

The restored workflows now use current `hc_*` wrappers consistently, including
for helpers that previously required manual object bridging:

- `hc_import_clusters()`
- `hc_user_specific_cluster_profiling()`
- `hc_col_anno_numerical()`
- `hc_meta_correlation_num()`
- `hc_cut_hclust()`

### Practical recommendation

If these documents are updated further, keep them on this GitHub-only
`workflows` branch or move them into a separate companion repository such as
`hcocena-workflows`. Avoid adding them back to the Bioconductor submission
branch unless they are intentionally converted into maintained package
vignettes.
