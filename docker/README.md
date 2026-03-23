# Docker

Build the image from the repository root:

```bash
docker build -f docker/Dockerfile -t hcocena .
```

The image includes:

- the local `hcocena` package installation
- a ready-to-use workspace at `/home/rstudio/hcocena`
- bundled reference files at `/home/rstudio/hcocena/reference_files`
- visible workflow notebooks at `/home/rstudio/hcocena/github_workflows`
- direct copies in the workspace root:
  `01_hcocena_main.Rmd`, `02_hcocena_satellite.Rmd`,
  `03_hcocena_main_seq_only.Rmd`
- preinstalled optional packages for common workflows, including
  `CALIBERrfimpute`, `RCy3`, `SpatialExperiment`, and `GSVA`
- empty `count_data`, `annotation_data`, and `output` directories

RStudio Server is exposed on port `8787` by the base image.
