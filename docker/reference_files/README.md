Reference files are not versioned in this public repository snapshot.

The Docker image still creates `/home/rstudio/reference_files/` so workflows
can use the standard path. Add the required supplementary files to that folder
before running analyses that depend on them.

Typical examples include transcription factor lists and pathway gene-set files
used by `hc_read_supplementary()`.
