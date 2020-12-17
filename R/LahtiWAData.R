#' Obtain the Lahti Western Adtult microbiome data
#'
#' Obtain the microbiome data from Lahti et al. (2014).
#'
#' @details
#' The LahtiWA dataset contains high-throughput genus-level microbiota profiling
#' with HITChip for 1006 western adults with no reported health complications,
#' reported in Lahti et al. (2014).
#'
#' Column metadata includes the age, sex, nationality, DNA extraction method,
#' project identifier, precomputed diversity measurement, BMI group,
#' subject identifier, time information and sample identifier.
#'
#' Row metadata contains taxonomic information on the Phylum, Family and Genus
#' level.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#'
#' @return a \linkS4class{TreeSummarizedExperiment}
#'
#' @author Felix G.M. Ernst
#'
#' @references
#' Lahti L et al. (2014): Tipping elements in the human intestinal ecosystem.
#' \emph{Nature Communications} 5:4344 (2014)
#' \url{https://doi.org/10.1038/ncomms5344}
#'
#' @name LahtiWAData
#' @export
#'
#' @examples
#' tse <- LahtiWAData()
LahtiWAData <- function() {
    dataset <- "lahti-wa"
    tse <- .create_tse(dataset,
                       assays = c("counts"),
                       has.rowdata = TRUE,
                       has.coldata = TRUE)
    tse
}
