% Manual entry added because full roxygen regeneration is currently blocked by
% an unrelated existing template issue elsewhere in the package.
\name{hc_longitudinal_workflow_direct}
\alias{hc_longitudinal_workflow_direct}
\title{Convenience longitudinal workflow using the direct step-1 backend}
\usage{
hc_longitudinal_workflow_direct(
  hc,
  donor_col = "Subject",
  time_col = "Time_token",
  layer = NULL,
  time_levels,
  k = 2:6,
  method = c("kmeans", "hclust"),
  nstart = 50,
  cap_runs = 50,
  impute = TRUE,
  impute_method = "rfcont",
  ntree = 10,
  min_cluster_fraction = 0.1,
  score_method = c("calinski_harabasz", "tot_withinss"),
  scale_features = TRUE,
  na_impute = c("median", "zero"),
  cap_scale_features = scale_features,
  cap_na_impute = na_impute,
  cap_jitter_sd = 0,
  seed = 42,
  means_slot = "longitudinal_module_means_direct",
  output_slot = "longitudinal_endotypes_direct"
)
}
\arguments{
\item{hc}{A \code{HCoCenaExperiment}.}

\item{donor_col}{Annotation column containing donor IDs.}

\item{time_col}{Annotation column containing timepoint labels.}

\item{layer}{Optional layer index, layer id, layer name, or \code{"all"}.}

\item{time_levels}{Required explicit timepoint order vector.}

\item{k}{Candidate donor-trajectory cluster numbers per module.}

\item{method}{One of \code{"kmeans"} or \code{"hclust"} for direct clustering
on the donor-by-time matrix.}

\item{nstart}{Number of random starts per direct clustering run.}

\item{cap_runs}{Number of repeated direct clustering runs used to build CAP.}

\item{impute}{Logical. If \code{TRUE}, impute missing donor-time values before
direct clustering.}

\item{impute_method}{Missing-value method passed to \code{mice}.}

\item{ntree}{Number of trees for \code{rfcont} imputation.}

\item{min_cluster_fraction}{Minimum allowed fraction of donors in the smallest
module cluster.}

\item{score_method}{One of \code{"calinski_harabasz"} or
\code{"tot_withinss"}.}

\item{scale_features}{Logical. If \code{TRUE}, z-score timepoint columns before
direct clustering.}

\item{na_impute}{One of \code{"median"} or \code{"zero"}.}

\item{cap_scale_features}{Logical. If \code{TRUE}, z-score features for CAP
reruns.}

\item{cap_na_impute}{One of \code{"median"} or \code{"zero"} for CAP reruns.}

\item{cap_jitter_sd}{Optional Gaussian jitter SD for direct CAP reruns.}

\item{seed}{Random seed used across the three quick steps.}

\item{means_slot}{Satellite slot for module means.}

\item{output_slot}{Satellite slot for direct endotype outputs.}
}
\value{
A list with updated \code{hc}, nested \code{plots}, and per-step diagnostics.
}
\description{
Runs the direct donor-by-time step 1 backend and then reuses the existing quick
step 2 and step 3 wrappers on the resulting output slot family.
}
\examples{
exists("hc_longitudinal_workflow_direct", mode = "function")

\donttest{
TRUE
}
}
