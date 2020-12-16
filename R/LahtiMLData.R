#' Obtain the Lahti Microbiome and Lipid data
#'
#' Obtain the microbiome and lipid data from Lahti et al. (2013).
#'
#' @details
#' The LahtiML dataset contains high-through put profiling data from 389 human
#' blood serum lipids and 130 intestinal genus-level bacteria from 44 samples
#' (22 subjects from 2 time points; before and after probiotic/placebo
#' intervention). The data set can be used to investigate associations between
#' intestinal bacteria and host lipid metabolism
#'
#' Column metadata includes the time points, sex, subject identifier, sample
#' identifier and treatment group.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#'
#' @return A \linkS4class{MultiAssayExperiment} object with
#'   \linkS4class{TreeSummarizedExperiment} and a
#'   \linkS4class{SummarizedExperiment}
#'
#' @author Felix G.M. Ernst
#'
#' @references
#' Lahti L et al. (2013).
#' Associations between the human intestinal microbiota, Lactobacillus rhamnosus
#' GG and serum lipids indicated by integrated analysis of high-throughput
#' profiling data.
#' \emph{PeerJ} 1:e32 \url{https://doi.org/10.7717/peerj.32}
#'
#' @export
#'
#' @examples
#' mae <- LahtiMLData()
LahtiMLData <- function() {
    mae <- .create_mae("lahti-ml",
                       types = list(microbiome = list("TSE" = c("counts")),
                                    lipids = list("SE" = c("counts"))),
                       coldata = TRUE,
                       samplemap = FALSE,
                       has.rowdata = list(microbiome = TRUE,
                                          lipids = FALSE),
                       has.coldata = list(microbiome = FALSE,
                                          lipids = FALSE))
    mae
}
