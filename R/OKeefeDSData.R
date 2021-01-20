#' Obtain the O'Keefe diet swap microbiome data
#'
#' Obtain the microbiome data from O'Keefe et al. (2015).
#'
#' @details
#' The OKeefeDS dataset contains microbiome data from a study with African and
#' African American groups undergoing a two-week diet swap.
#'
#' This data set is based on the Human Intestinal Tract (HIT)Chip
#' phylogenetic 16S microarray (Rajilić‐Stojanović _et al._ 2009.
#' This profiling technology differs from the more widely used 16S rRNA
#' amplicon sequencing.
#'
#' Column metadata includes the subject identifier, sex, nationality, group
#' information, sample identifier, time point information, time point
#' information within group and BMI group.
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
#' O'Keefe S et al. (2015): Fat, fibre and cancer risk in African Americans and
#' rural Africans. \emph{Nature Communications} 6:6342 (2015)
#' \url{https://dx.doi.org/10.1038/ncomms7342}
#'
#' Rajilić-Stojanović M, Heilig HG, Molenaar D, Kajander K, Surakka A, Smidt H,
#' de Vos WM (2009). Development and application of the human intestinal tract chip, a
#' phylogenetic microarray: analysis of universally conserved phylotypes in the
#' abundant microbiota of young and elderly adults.
#' \emph{Environ Microbiol.} 11(7):1736-51
#' \url{https://doi.org/10.1111/j.1462-2920.2009.01900.x}
#'
#' @name OKeefeDSData
#' @export
#'
#' @examples
#' tse <- OKeefeDSData()
OKeefeDSData <- function() {
    dataset <- "3.13/okeefe-ds"
    tse <- .create_tse(dataset,
                       assays = c("counts"),
                       has.rowdata = TRUE,
                       has.coldata = TRUE)
    tse
}

#' @rdname OKeefeDSData
#' @export
dietswap <- OKeefeDSData
