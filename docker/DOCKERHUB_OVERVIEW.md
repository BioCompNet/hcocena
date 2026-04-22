# hcocena Docker image

Public RStudio-based image for `hcocena`, the R package for horizontal
integration and downstream analysis of transcriptomics datasets.

## What is included

- the `hcocena` package preinstalled in the container
- an RStudio workspace at `/home/rstudio/hcocena`
- bundled `reference_files/`
- visible workflow notebooks:
  - `01_hcocena_main.Rmd`
  - `02_hcocena_satellite.Rmd`
  - `03_hcocena_main_seq_only.Rmd`
- additional notebooks under `/home/rstudio/hcocena/github_workflows/`

## Recommended tags

- `1.96` for a pinned, reproducible setup
- `latest` if you prefer the moving convenience tag

Older tags still kept for older reproducible runs:

- `1.95`
- `1.94`
- `1.9`
- `1.28`
- `1.1.2`

## Quick start

Pull the current image:

```bash
docker pull therealtomek/hcocena:1.96
```

Run RStudio Server:

```bash
docker run --rm -p 8787:8787 -e PASSWORD=hcocena therealtomek/hcocena:1.96
```

Then open:

```text
http://localhost:8787
```

Login:

- user: `rstudio`
- password: the value passed in `PASSWORD`

## Project links

- GitHub repository: https://github.com/BioCompNet/hcocena
- Docker Hub tags: https://hub.docker.com/r/therealtomek/hcocena/tags

## Notes

The image is intended as a ready-to-use environment for exploring the package
and running the main and satellite workflows in RStudio.
