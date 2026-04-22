# Docker

Build the image from the repository root:

```bash
docker build -f docker/Dockerfile -t hcocena .
```

The image includes:

- the local `hcocena` package installation
- a ready-to-use workspace directly at `/home/rstudio`
- TB-neo-specific `01_hcocena_main.Rmd` and `02_hcocena_satellite.Rmd`
- a prepared `/home/rstudio/reference_files` folder with a README for adding
  your reference files
- an empty `/home/rstudio/project_data` folder with instructions for
  adding the project-specific `.rds` inputs yourself
- workflow notebooks directly in the workspace root:
  `01_hcocena_main.Rmd`, `02_hcocena_satellite.Rmd`, `03_hcocena_main_seq_only.Rmd`
- preinstalled optional packages for common workflows, including
  `graphlayouts`, `ellmer`, `CALIBERrfimpute`, `RCy3`, `SpatialExperiment`, and `GSVA`
- empty `count_data`, `annotation_data`, and `output` directories

RStudio Server is exposed on port `8787` by the base image.

## Updating the Docker Hub description

When publishing a new public Docker tag, update the overview file from the
repository root with:

```bash
Rscript docker/update_dockerhub_overview.R <new-tag>
```

This updates [`DOCKERHUB_OVERVIEW.md`](DOCKERHUB_OVERVIEW.md) so the new tag
becomes the recommended pinned version, moves the previous recommended tag into
the older reproducibility list, and refreshes the `docker pull` / `docker run`
examples.

If Docker Hub is not syncing this file automatically, copy the updated markdown
to the repository Overview field there right after pushing the new tag.
