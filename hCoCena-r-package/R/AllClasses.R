#' hCoCena S4 Classes
#'
#' S4 container classes used by the object-oriented hCoCena API.
#'
#' @name HCoCenaClasses
#' @docType class
#' @aliases
#'   HCoCenaConfig-class
#'   HCoCenaReferences-class
#'   HCoCenaLayerResult-class
#'   HCoCenaIntegration-class
#'   HCoCenaExperiment-class
#' @keywords classes
#' @importClassesFrom S4Vectors DataFrame SimpleList
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
NULL

.hc_empty_mae <- function() {
  MultiAssayExperiment::MultiAssayExperiment()
}

setClass(
  "HCoCenaConfig",
  slots = c(
    global = "DataFrame",
    layer = "DataFrame",
    paths = "DataFrame"
  ),
  prototype = prototype(
    global = S4Vectors::DataFrame(),
    layer = S4Vectors::DataFrame(),
    paths = S4Vectors::DataFrame()
  )
)

setClass(
  "HCoCenaReferences",
  slots = c(
    registry = "DataFrame",
    data = "SimpleList"
  ),
  prototype = prototype(
    registry = S4Vectors::DataFrame(),
    data = S4Vectors::SimpleList()
  )
)

setClass(
  "HCoCenaLayerResult",
  slots = c(
    part1 = "SimpleList",
    part2 = "SimpleList"
  ),
  prototype = prototype(
    part1 = S4Vectors::SimpleList(),
    part2 = S4Vectors::SimpleList()
  )
)

setClass(
  "HCoCenaIntegration",
  slots = c(
    combined_edgelist = "DataFrame",
    graph = "ANY",
    gfc = "DataFrame",
    cluster = "SimpleList"
  ),
  prototype = prototype(
    combined_edgelist = S4Vectors::DataFrame(),
    graph = NULL,
    gfc = S4Vectors::DataFrame(),
    cluster = S4Vectors::SimpleList()
  )
)

setClass(
  "HCoCenaExperiment",
  slots = c(
    mae = "MultiAssayExperiment",
    config = "HCoCenaConfig",
    references = "HCoCenaReferences",
    layer_results = "SimpleList",
    integration = "HCoCenaIntegration",
    satellite = "SimpleList",
    provenance = "DataFrame"
  ),
  prototype = prototype(
    mae = .hc_empty_mae(),
    config = new("HCoCenaConfig"),
    references = new("HCoCenaReferences"),
    layer_results = S4Vectors::SimpleList(),
    integration = new("HCoCenaIntegration"),
    satellite = S4Vectors::SimpleList(),
    provenance = S4Vectors::DataFrame()
  )
)
