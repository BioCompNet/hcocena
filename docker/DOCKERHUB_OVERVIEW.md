# hcocena Docker image

Public RStudio-based image for `hcocena`, the R package for horizontal
integration and downstream analysis of transcriptomics datasets.

## What is included

- the `hcocena` package preinstalled in the container
- an RStudio workspace at `/home/rstudio/hcocena`
- bundled `reference_files/` with pathway, GO, hallmark, TF, and immune helper references
- visible workflow notebooks:
  - `01_hcocena_main.Rmd`
  - `02_hcocena_satellite.Rmd`
  - `03_hcocena_main_seq_only.Rmd`
- additional notebooks under `/home/rstudio/hcocena/github_workflows/`
- RNA-seq and differential-expression packages including `DESeq2`, `limma`,
  `sva`, `edgeR`, `tximport`, `apeglm`, `ashr`, `EnhancedVolcano`,
  `pheatmap`, `org.Hs.eg.db`, and `org.Mm.eg.db`

## Recommended tags

- `latest` for the current image and quick-start commands
- `1.100` for a pinned, reproducible setup

Older tags still kept for older reproducible runs:

- `1.99`
- `1.98`
- `1.97`
- `1.96`
- `1.95`
- `1.94`
- `1.9`
- `1.28`
- `1.1.2`

## Quick start

Pull the current image:

```bash
docker pull therealtomek/hcocena:latest
```

Run RStudio Server:

```bash
docker run --rm -p 8787:8787 -e PASSWORD=hcocena therealtomek/hcocena:latest
```

For reproducible runs, replace `latest` with a pinned tag such as `1.100`.

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
