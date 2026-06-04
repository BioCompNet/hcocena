test_that("real-data regression runner works when explicitly enabled", {
  run_realdata <- tolower(Sys.getenv("HCOCENA_RUN_REALDATA", unset = "false")) %in%
    c("1", "true", "yes", "y", "on")
  skip_if_not(run_realdata, "Set HCOCENA_RUN_REALDATA=true to run local real-data regression tests.")

  repo_root <- normalizePath(test_path("..", ".."), winslash = "/", mustWork = TRUE)
  runner <- file.path(repo_root, "scripts", "run_realdata_regression.R")
  skip_if_not(file.exists(runner), "Real-data regression runner is unavailable.")
  source(runner)

  mode <- Sys.getenv("HCOCENA_REALDATA_MODE", unset = "quick")
  data_dir <- Sys.getenv(
    "HCOCENA_REALDATA_DIR",
    unset = file.path(dirname(repo_root), "data")
  )
  skip_if_not(dir.exists(data_dir), paste0("Real-data directory is unavailable: ", data_dir))
  reference_parent <- Sys.getenv("HCOCENA_REALDATA_REFERENCE_DIR", unset = "")
  reference_dir <- if (nzchar(reference_parent)) {
    file.path(reference_parent, mode)
  } else {
    file.path(repo_root, "realdata-reference", mode)
  }

  out <- run_hcocena_realdata_regression(
    mode = mode,
    data_dir = data_dir,
    output_dir = file.path(repo_root, "realdata-output", paste0("test-", mode)),
    reference_dir = reference_dir,
    compare_reference = dir.exists(reference_dir),
    update_reference = FALSE,
    load_source = TRUE
  )

  expect_s3_class(out, "hcocena_realdata_regression")
  expect_true(file.exists(file.path(out$output_dir, "manifest.json")))
  expect_true(file.exists(file.path(out$output_dir, "font_probe.csv")))
  expect_true(file.exists(file.path(out$output_dir, "plot_variants.csv")))
  expect_true(file.exists(file.path(out$output_dir, "plot_variant_parameters.csv")))
  expect_true(file.exists(file.path(out$output_dir, paste0("visual_check_report_", mode, ".pdf"))))

  cfg <- .hcr_mode_config(mode)
  plot_variants <- utils::read.csv(
    file.path(out$output_dir, "plot_variants.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_true(any(plot_variants$plot == "Cluster heatmap" & plot_variants$variant == "baseline before module split"))
  if (isTRUE(cfg$run_module_split)) {
    expect_true(any(plot_variants$plot == "Cluster heatmap" & grepl("after split", plot_variants$variant)))
  }
  if (isTRUE(cfg$run_enrichment)) {
    expect_true(any(plot_variants$plot == "Functional enrichment" & plot_variants$variant == "combined all DBs"))
  }
  png_files <- plot_variants$png_file[nzchar(plot_variants$png_file)]
  expect_true(all(file.exists(file.path(out$output_dir, png_files))))

  if (!is.null(out$comparison)) {
    expect_true(all(out$comparison$status == "ok"))
  }
})
