#' Obtain the O'Keefe diet swap microbiome data
#'
#' Obtain the microbiome data from O'Keefe et al. (2015).
#'
#' @details
#' The OKeefeDS dataset contains microbiome data from a study with African and
#' African American groups undergoing a two-week diet swap.
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
#' @author Felix G.M. Ernst
#'
#' @references
#' O'Keefe S et al. (2015): Fat, fibre and cancer risk in African Americans and
#' rural Africans. \emph{Nature Communications} 6:6342 (2015)
#' \url{https://dx.doi.org/10.1038/ncomms7342}
#'
#' @name OKeefeDSData
#' @export
#'
#' @examples
#' tse <- OKeefeDSData()
OKeefeDSData <- function() {
    dataset <- "okeefe-ds"
    tse <- .create_tse(dataset,
                       assays = c("counts"),
                       has.rowdata = TRUE,
                       has.coldata = TRUE)
    tse
}
