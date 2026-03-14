test_that("empty HCoCenaExperiment is valid", {
  hc <- hc_init()
  expect_s4_class(hc, "HCoCenaExperiment")
  expect_true(methods::validObject(hc))
})

test_that("layer-name mismatches are invalid", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "s1"))),
    colData = S4Vectors::DataFrame(group = "A", row.names = "s1")
  )

  hc <- hc_init()
  hc@mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(set2 = se)
  )
  hc@config@layer <- S4Vectors::DataFrame(layer_id = "set1")
  hc@layer_results <- S4Vectors::SimpleList(set1 = new("HCoCenaLayerResult"))

  expect_error(methods::validObject(hc), "Layer names")
})
