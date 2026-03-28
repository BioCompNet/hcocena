# Docker

Build the image from the repository root:

```bash
docker build -f docker/Dockerfile -t hcocena .
```

The image includes:

- the local `hcocena` package installation
- a ready-to-use workspace directly at `/home/rstudio`
- TB-neo-specific `01_hcocena_main.Rmd` and `02_hcocena_satellite.Rmd`
- bundled reference files at `/home/rstudio/reference_files`
- an empty `/home/rstudio/project_data` folder with instructions for
  adding the project-specific `.rds` inputs yourself
- workflow notebooks directly in the workspace root:
  `01_hcocena_main.Rmd`, `02_hcocena_satellite.Rmd`, `03_hcocena_main_seq_only.Rmd`
- preinstalled optional packages for common workflows, including
  `graphlayouts`, `ellmer`, `CALIBERrfimpute`, `RCy3`, `SpatialExperiment`, and `GSVA`
- empty `count_data`, `annotation_data`, and `output` directories

RStudio Server is exposed on port `8787` by the base image.
