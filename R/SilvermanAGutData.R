#' Obtain the Silverman Artificial Gut data
#'
#' Obtain the microbiome data from Silverman et al. (2018).
#'
#' @details
#' The SilvermanAGutData dataset contains 16S rRNA gene based high-throughput
#' profiling of 4 in vitro artificial gut models. The sampling was done hourly
#' and daily to capture sub-daily dynamics of microbial community originating
#' from human feces. The data consists of 413 taxa from 639 samples. The data
#' set can be used to investigate longitudinal dynamics of microbial community
#' in a controlled environment.
#'
#' Column metadata includes the days of sampling, vessel identifier, sampling 
#' frequency pre-post challenge with Bacteroides ovatus.
#'
#' The wow metadata of the microbiome data contains taxonomic information on the
#' Kingdom, Phylum, Class, Order, Family and Genus and Species level.  
#' 
#' The row tree consists of a phylogenetic tree build using sequence information
#' of 413 taxa.  
#' 
#' As reference sequences the ASV are provided.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#'
#' @return a \linkS4class{TreeSummarizedExperiment}
#'
#' @author Sudarshan A. Shetty and Felix G.M. Ernst
#'
#' @references
#' Silveman J.D et al. (2018): Dynamic linear models guide design and 
#' analysis of microbiota studies within artificial human guts. 
#' \emph{Microbiome} 6:202 \url{https://doi.org/10.1186/s40168-018-0584-3}
#'
#' @name SilvermanAGutData
#' @export
#'
#' @examples
#' tse <- SilvermanAGutData()
SilvermanAGutData <- function() {
    dataset <- "3.13/silverman-ag"
    tse <- .create_tse(dataset,
                    assays = c("counts"),
                    has.rowdata = TRUE,
                    has.coldata = TRUE,
                    has.rowtree = TRUE,
                    has.refseq = TRUE)
    tse
}

#' @rdname SilvermanAGutData
#' @export
artificialgut <- SilvermanAGutData
