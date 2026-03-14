# Docker

Build the image from the repository root:

```bash
docker build -f docker/Dockerfile -t hcocena .
```

The image includes:

- the local `hcocena` package installation
- a ready-to-use workspace at `/home/rstudio/hcocena`
- bundled reference files at `/home/rstudio/hcocena/reference_files`
- preinstalled optional packages for common workflows, including
  `CALIBERrfimpute`, `RCy3`, `SpatialExperiment`, and `GSVA`
- empty `count_data`, `annotation_data`, and `output` directories

RStudio Server is exposed on port `8787` by the base image.
