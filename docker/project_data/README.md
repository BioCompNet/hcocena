Project-specific TB-neo input files are intentionally not bundled into the
public Docker image.

Before running `01_hcocena_main.Rmd`, place these files into
`/home/rstudio/project_data/`:

- `sample_table_230530.rds`
- `vst_anno_log_230530.rds`
- `exclusion_vector_lm_full_240625.rds`

The Docker notebooks already point to this folder.
