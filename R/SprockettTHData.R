#' Obtain the Sprockett Tsimane Horticulturalists data
#'
#' Obtain the microbiome data from Sprockett et al. (2020).
#'
#' @details
#' The SprockettTHData dataset contains 16S rRNA gene based high-throughput
#' profiling of 1966 Feces, 120 Saliva samples from 319 participants from 
#' Bolivia, Finland and Bangladesh. These include samples from adults, 
#' children, and infants. Several participants have longitudinal samples. 
#' The data consists of 2319 taxa from 2086 samples. The data
#' set can be used to investigate assembly, structure, and dynamics as well 
#' as associations between several host related parameters with microbiota.
#'
#' Column metadata includes the sex, age, feeding status, delivery mode, 
#' country, and other information.
#'
#' The row metadata of the microbiome data contains taxonomic information 
#' on the Kingdom, Phylum, Class, Order, Family and Genus, Species and 
#' lowest taxonomic rank.  
#' 
#' The row tree consists of a phylogenetic tree build using sequence 
#' information of 2319 taxa.
#' 
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
#' Sprockett, D.D., Martin, M., Costello, E.K. et al. (2020)
#' Microbiota assembly, structure, and dynamics among Tsimane
#' horticulturalists of the Bolivian Amazon. 
#' \emph{Nat Commun} 11, 3772 \url{https://doi.org/10.1038/s41467-020-17541-6}
#'
#' Subramanian, S., Huq, S., Yatsunenko, T., et al. (2014) Persistent gut
#' microbiota immaturity in malnourished Bangladeshi children. 
#' \emph{Nature} 510, 417-421 \url{https://doi.org/10.1038/nature13421}
#' 
#' Vatanen, T., Kostic A.D., d'Hennezel E., et al. (2016) Variation in
#' microbiome LPS immunogenicity contributes to autoimmunity in humans. 
#' \emph{Cell} 165, 842-853 \url{https://doi.org/10.1016/j.cell.2016.04.007}
#' 
#' @name SprockettTHData
#' @export
#'
#' @examples
#' tse <- SprockettTHData()
SprockettTHData <- function() {
    dataset <- "3.13/sprockett-th"
    tse <- .create_tse(dataset,
                    assays = c("counts"),
                    has.rowdata = TRUE,
                    has.coldata = TRUE,
                    has.rowtree = TRUE,
                    has.refseq = TRUE)
    tse
}
