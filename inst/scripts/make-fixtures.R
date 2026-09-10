# Provenance for the example fixtures in inst/extdata.
#
# The fixtures are small HCoCenaExperiment objects captured at four stages of
# the workflow, used by the examples on the man pages and by a regression test.
# This script documents how the annotation carried by those objects is shaped;
# run it from the package root with the native Rscript.
#
#   Rscript inst/scripts/make-fixtures.R
#
# It augments the existing fixtures in place rather than regenerating the
# expression data, so the numeric content of the objects (and therefore the
# regression test that reads hc_after_part1.rds) stays untouched.

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})
pkgload::load_all(".", quiet = TRUE)

# saveRDS() runs DelayedArray's .S4_object_contains_out_of_memory_data hook,
# which walks the S4 slots recursively and overflows the node stack on the
# larger fixtures. serialize() to a gzfile writes the same format without that
# hook, and readRDS() reads the result back unchanged.
write_rds <- function(x, path) {
  con <- gzfile(path, "wb")
  on.exit(close(con), add = TRUE)
  serialize(x, con)
}

extdata <- file.path("inst", "extdata")
fixtures <- c("hc_prepared.rds", "hc_after_part1.rds",
              "hc_after_part2.rds", "hc_clustered.rds")

# A longitudinal design on top of the existing samples: the six samples per
# layer become three donors followed across two timepoints. This is what the
# hc_longitudinal_* examples need; without donor/time columns they cannot run.
# Three donors is the minimum the endotype and meta-clustering steps accept,
# which is why the design is 3x2 rather than 2x3.
add_longitudinal_columns <- function(hc) {
  exps <- MultiAssayExperiment::experiments(hc@mae)
  for (nm in names(exps)) {
    se <- exps[[nm]]
    cd <- SummarizedExperiment::colData(se)
    n <- nrow(cd)
    cd$donor <- rep(c("D1", "D2", "D3"), length.out = n)
    cd$timepoint <- rep(c("T1", "T2"), each = 3, length.out = n)
    SummarizedExperiment::colData(se) <- cd
    exps[[nm]] <- se
  }
  MultiAssayExperiment::experiments(hc@mae) <- exps
  hc
}

for (f in fixtures) {
  p <- file.path(extdata, f)
  if (!file.exists(p)) {
    message("skipping missing fixture: ", f)
    next
  }
  hc <- readRDS(p)
  hc <- add_longitudinal_columns(hc)
  write_rds(hc, p)
  # validObject() recurses too deeply on the larger fixtures, so verify by
  # reading the file back instead
  chk <- readRDS(p)
  stopifnot(inherits(chk, "HCoCenaExperiment"))
  cd <- SummarizedExperiment::colData(
    MultiAssayExperiment::experiments(chk@mae)[[1]]
  )
  stopifnot(all(c("donor", "timepoint") %in% colnames(cd)))
  message(f, ": ", paste(colnames(cd), collapse = ", "))
}
