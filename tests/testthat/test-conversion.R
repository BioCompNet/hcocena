test_that("legacy <-> S4 conversion roundtrip preserves key content", {
  legacy <- list(
    working_directory = list(
      dir_count_data = FALSE,
      dir_annotation = FALSE,
      dir_reference_files = tempdir(),
      dir_output = tempdir()
    ),
    data = list(
      set1_counts = matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("g1", "g2"), c("s1", "s2"))
      ),
      set1_anno = data.frame(group = c("A", "B"), row.names = c("s1", "s2"))
    ),
    supplementary_data = list(),
    global_settings = list(voi = "group", range_GFC = 2),
    layer_settings = list(set1 = list(top_var = "all", min_corr = 0.9, range_cutoff_length = 100)),
    layers = list(set1 = c("set1_counts", "set1_anno")),
    supplement = list(),
    layers_names = c("Layer1"),
    layer_specific_outputs = list(),
    integrated_output = list(),
    cutoff_vec = c(0.95),
    satellite_outputs = list()
  )

  hc <- as_hcocena(legacy)
  expect_s4_class(hc, "HCoCenaExperiment")

  legacy2 <- as_hcobject(hc)
  expect_true("set1_counts" %in% names(legacy2$data))
  expect_equal(legacy2$global_settings$voi, "group")
  expect_equal(legacy2$cutoff_vec[[1]], 0.95)
})

test_that("conversion handles NULL entries in supplement registry", {
  legacy <- list(
    working_directory = list(),
    data = list(),
    supplementary_data = list(),
    global_settings = list(),
    layer_settings = list(),
    layers = list(),
    supplement = list(
      Tf = "TFcat.txt",
      Hallmark = "hallmark.gmt",
      Reactome = NULL,
      CustomGoCC = "custom_go_cc.gmt"
    ),
    layers_names = list(),
    layer_specific_outputs = list(),
    integrated_output = list(),
    cutoff_vec = NULL,
    satellite_outputs = list()
  )

  hc <- as_hcocena(legacy)
  reg <- as.data.frame(hc@references@registry)

  expect_s4_class(hc, "HCoCenaExperiment")
  expect_equal(nrow(reg), 3)
  expect_true(all(c("Tf", "Hallmark", "CustomGoCC") %in% reg$name))
  expect_false("Reactome" %in% reg$name)
})
