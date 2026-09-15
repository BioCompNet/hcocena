## An S4 call must not touch a variable in the user's workspace that happens to
## be called `hcobject`. The bridge used to mirror its whole legacy object into
## `.GlobalEnv` whenever one was found there.

test_that("S4 calls leave a global `hcobject` alone", {
  path <- system.file("extdata", "hc_clustered.rds", package = "hcocena")
  skip_if(!nzchar(path), "The clustered fixture is unavailable.")

  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  grDevices::pdf(file.path(dir, "plots.pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  hc <- readRDS(path)
  marker <- list(marker = "untouched", data = list())

  had_before <- exists("hcobject", envir = globalenv(), inherits = FALSE)
  previous <- if (had_before) get("hcobject", envir = globalenv(), inherits = FALSE) else NULL
  withr::defer({
    if (had_before) {
      assign("hcobject", previous, envir = globalenv())
    } else {
      suppressWarnings(rm("hcobject", envir = globalenv()))
    }
  })

  calls <- list(
    hc_get_module_scores = function() hc_get_module_scores(hc),
    hc_find_hubs = function() hc_find_hubs(hc, top = 3),
    hc_gene_to_cluster = function() hc_gene_to_cluster(hc),
    hc_build_integrated_network = function() hc_build_integrated_network(hc, mode = "u")
  )

  for (nm in names(calls)) {
    assign("hcobject", marker, envir = globalenv())
    suppressMessages(suppressWarnings(try(calls[[nm]](), silent = TRUE)))
    expect_identical(
      get("hcobject", envir = globalenv(), inherits = FALSE), marker,
      info = paste("global hcobject changed by", nm)
    )
  }
})
