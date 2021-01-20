#' Obtain the Lahti Western Adult microbiome data
#'
#' Obtain the microbiome data from Lahti et al. (2014).
#'
#' @details
#' The LahtiWA dataset contains high-throughput genus-level microbiota profiling
#' with HITChip for 1006 western adults with no reported health complications,
#' reported in Lahti et al. (2014).
#'
#' This data set is based on the Human Intestinal Tract (HIT)Chip
#' phylogenetic 16S microarray (Rajilić‐Stojanović _et al._ 2009.
#' This profiling technology differs from the more widely used 16S rRNA
#' amplicon sequencing.
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
#' @author Felix G.M. Ernst and Leo Lahti
#'
#' @references
#' Lahti L et al. (2014): Tipping elements in the human intestinal ecosystem.
#' \emph{Nature Communications} 5:4344 (2014)
#' \url{https://doi.org/10.1038/ncomms5344}
#'
#' Rajilić-Stojanović M, Heilig HG, Molenaar D, Kajander K, Surakka A, Smidt H,
#' de Vos WM (2009). Development and application of the human intestinal tract chip, a
#' phylogenetic microarray: analysis of universally conserved phylotypes in the
#' abundant microbiota of young and elderly adults.
#' \emph{Environ Microbiol.} 11(7):1736-51
#' \url{https://doi.org/10.1111/j.1462-2920.2009.01900.x}
#'
#' @name LahtiWAData
#' @export
#'
#' @examples
#' tse <- LahtiWAData()
LahtiWAData <- function() {
    dataset <- "3.13/lahti-wa"
    tse <- .create_tse(dataset,
                       assays = c("counts"),
                       has.rowdata = TRUE,
                       has.coldata = TRUE)
    tse
}

#' @rdname LahtiWAData
#' @export
atlas1006 <- LahtiWAData
