#' Obtain the Human Microbiome Project Consortium dataset
#'
#' Obtain the selected microbiome datasets from Human Microbiome Project 
#' Consortium publication (2012).
#'
#' @details
#' The HMPC dataset comes in two different flavours, both based on 16s rRNA
#' sequencing using pyrosequencing. HMPCV13Data contains information generated 
#' from sequencing the variable regions 1-3, whereas HMPCV35Data contains
#' the information extracted from sequencing region 3-5.
#'
#' Column metadata includes among others information on body sites and 
#' extraction methods used.
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
#' Human Microbiome Project Consortium. Structure, function and diversity of the
#' healthy human microbiome. Nature. 2012 Jun 13;486(7402):207-14. doi: 
#' 10.1038/nature11234. PMID: 22699609; PMCID: PMC3564958.
#'
#' @name HMPCData
#' @export
#'
#' @examples
#' tse <- HMPCV13Data()
HMPCV13Data <- function() {
    dataset <- "hmpc-v13"
    tse <- .create_tse(dataset,
                       assays = c("counts"),
                       has.rowdata = TRUE,
                       has.coldata = TRUE,
                       has.rowtree = TRUE)
    tse
}

#' @name HMPCData
#' @export
HMPCV35Data <- function() {
    dataset <- "hmpc-v35"
    tse <- .create_tse(dataset,
                       assays = c("counts"),
                       has.rowdata = TRUE,
                       has.coldata = TRUE,
                       has.rowtree = TRUE)
    tse
}
