\name{hc_s4_api}
\alias{hc_init}
\alias{hc_set_paths}
\alias{hc_define_layers}
\alias{hc_read_data}
\alias{hc_set_global_settings}
\alias{hc_set_layer_settings}
\alias{hc_set_supp_files}
\alias{hc_read_supplementary}
\alias{hc_run_expression_analysis_1}
\alias{hc_set_cutoff}
\alias{hc_run_expression_analysis_2}
\alias{hc_build_integrated_network}
\alias{hc_cluster_calculation}
\alias{hc_functional_enrichment}
\alias{hc_get_reference_data}
\alias{as_hcocena}
\alias{as_hcobject}
\title{S4-first hCoCena API}
\description{Object-oriented entry points for hCoCena analyses based on \code{HCoCenaExperiment}.}
\usage{
hc_init()

hc_set_paths(hc, dir_count_data, dir_annotation, dir_reference_files, dir_output)

hc_define_layers(hc, data_sets = list())

hc_read_data(hc, sep_counts = "\t", sep_anno = "\t", gene_symbol_col = NULL,
  sample_col = NULL, count_has_rn = TRUE, anno_has_rn = TRUE,
  auto_setup_output = TRUE, create_output_dir = TRUE,
  project_folder = NULL)

hc_set_global_settings(hc, organism, control_keyword, variable_of_interest,
  min_nodes_number_for_network = 15, min_nodes_number_for_cluster = 15,
  range_GFC = 2, layout_algorithm = "layout_with_fr", data_in_log)

hc_set_layer_settings(hc, top_var, min_corr = 0.7, range_cutoff_length,
  print_distribution_plots = FALSE)

hc_set_supp_files(hc, Tf = NULL, Hallmark = NULL, Go = NULL, Kegg = NULL,
  Reactome = NULL, ...)

hc_read_supplementary(hc)

hc_run_expression_analysis_1(hc, padj = "none", export = FALSE, import = NULL,
  bayes = FALSE, prior = 2, alpha = 0.5, corr_method = "pearson")

hc_set_cutoff(hc, cutoff_vector = c(), auto = FALSE, fallback_cutoff = 0.982,
  verbose = TRUE)

hc_run_expression_analysis_2(hc, grouping_v = NULL, plot_HM = TRUE,
  method = "complete", additional_anno = NULL, cols = NULL)

hc_build_integrated_network(hc, mode = "u", with = NULL, multi_edges = "min",
  GFC_when_missing = NULL)

hc_cluster_calculation(hc, cluster_algo = "cluster_leiden", no_of_iterations = 2,
  resolution = 0.1, partition_type = "RBConfigurationVertexPartition",
  max_cluster_count_per_gene = 1, return_result = FALSE)

hc_functional_enrichment(hc, gene_sets = c("Go", "Kegg", "Hallmark", "Reactome"),
  top = 5, clusters = c("all"), padj = "BH", qval = 0.05,
  heatmap_side = "left", heatmap_cluster_rows = FALSE,
  heatmap_cluster_columns = FALSE, heatmap_show_row_dend = FALSE,
  heatmap_show_column_dend = FALSE, heatmap_order = NULL,
  heatmap_module_label_mode = "same", heatmap_show_gene_counts = FALSE)

hc_get_reference_data(hub_id, hub = NULL)

as_hcocena(hcobject)

as_hcobject(hc)
}
\arguments{
  \item{hc}{An \code{HCoCenaExperiment} object.}
  \item{data_sets}{Named list of layers; each entry stores count and annotation source.}
  \item{auto_setup_output}{Boolean. If \code{TRUE}, run output-directory checks and initialize the save folder before importing data (if \code{dir_output} is set).}
  \item{create_output_dir}{Boolean passed to \code{hc_check_dirs()} when \code{auto_setup_output = TRUE}.}
  \item{project_folder}{Optional save subfolder name. Use \code{""} for direct output into \code{dir_output}. If \code{NULL}, existing \code{save_folder} is reused or defaults to \code{""}.}
  \item{hub_id}{ExperimentHub record identifier.}
  \item{hcobject}{Legacy list-based object.}
}
\value{Updated \code{HCoCenaExperiment} object or converted legacy object, depending on function.}
\details{Most \code{hc_*} functions currently orchestrate the legacy internals and therefore retain feature parity while moving users to an S4 API.}
