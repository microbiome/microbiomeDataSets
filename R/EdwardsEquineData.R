#' Obtain the Edwards Equine Data 
#'
#' Obtain the microbiome data from Edwards et al. (2020).
#'
#' @details
#' The EdwardsEquineData dataset contains 16S rRNA gene based high-throughput
#' profiling of bacterial and archaeal composition and anaerobic fungal 
#' profiling using 18S rRNA gene (~130bp), full ITS1 region and partial 5.8S 
#' rRNA gene (~31bp) from fecal samples of equines (n = 70) that included five 
#' extant species (i.e. Equus ferus, Equus africanus, Equus quagga, Equus zebra 
#' and Equus greyvii), as well as mules and hinnies (i.e. horse × donkey). The 
#' raw sequencing data were processed using NG-Tax.
#'
#' Column metadata includes the sampled host, country, location, Breeding,
#' warm or cold.blood and weight group and other technical information related 
#' to sequencing experiment. For bacterial data, number of equines is 70 while 
#' for fungal profiling, data from 64 equines is available.
#'
#' The row metadata for both bacteria and fungi of the microbiome data contains 
#' taxonomic information on the Kingdom, Phylum, Class, Order, Family and 
#' Genus level.  
#' 
#' The bacterial row tree consists of a phylogenetic tree build using sequence 
#' information of 2118 taxa. The fungal row tree consists of a phylogenetic 
#' tree build using sequence information of 358 taxa  
#' 
#' Reference sequences are not available.  
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#'
#' @return a \linkS4class{TreeSummarizedExperiment}
#'
#' @author Sudarshan A. Shetty 
#'
#' @references
#' Edwards J.D et al. (2020): Multi-kingdom characterization of the core equine 
#' fecal microbiota based on multiple equine (sub)species. 
#' \emph{Animal Microbiome} 6:202 \url{https://doi.org/10.1186/s42523-020-0023-1}
#'
#' @name EdwardsEquineData
#' @export
#'
#' @examples
#' mae <- EdwardsEquineData()
EdwardsEquineData <- function() {
  mae <- .create_mae("3.14/edwards-equine",
                     types = list(bacteria = list("TSE" = c("counts")),
                                  fungi = list("TSE" = c("counts"))),
                     coldata = TRUE,
                     samplemap = FALSE,
                     has.rowdata = list(bacteria = TRUE,
                                        fungi = TRUE),
                     has.coldata = list(bacteria = TRUE,
                                        fungi = TRUE),
                     has.rowtree = list(bacteria = TRUE,
                                        fungi = TRUE),
                     has.refseq = list(bacteria = TRUE,
                                        fungi = TRUE))
  mae
}

#' @rdname EdwardsEquineData
#' @export
equinegut <- EdwardsEquineData
