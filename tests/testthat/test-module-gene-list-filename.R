test_that("module gene-list filenames preserve the unsplit export", {
  filename <- get(".hc_module_gene_list_filename", asNamespace("hcocena"))

  expect_identical(
    filename(c(red = "M1", blue = "M2")),
    "Module_Gene_List.xlsx"
  )
  expect_identical(
    filename(c(red_1 = "M1.1", red_2 = "M1.2", blue = "M2")),
    "Module_Gene_splitted_List.xlsx"
  )
  expect_identical(
    filename(c(red_1_2 = "M1.1.2", blue = "M2")),
    "Module_Gene_splitted_List.xlsx"
  )
  expect_identical(
    filename(
      c(red = "M1", blue = "M2"),
      split_history = list(list(timestamp = "2026-07-29"))
    ),
    "Module_Gene_splitted_List.xlsx"
  )
})

test_that("module gene-list tables reflect current split labels", {
  build_table <- get(
    ".hc_module_gene_list_from_cluster_info",
    asNamespace("hcocena")
  )
  cluster_info <- data.frame(
    color = c("red_1", "red_2", "white"),
    gene_n = c("A,B", "C,D", "E"),
    cluster_included = c("yes", "yes", "no"),
    stringsAsFactors = FALSE
  )

  result <- build_table(
    cluster_info = cluster_info,
    module_label_map = c(red_1 = "M1.1", red_2 = "M1.2")
  )

  expect_identical(result$genes, c("A", "B", "C", "D"))
  expect_identical(result$module, c("M1.1", "M1.1", "M1.2", "M1.2"))
})
