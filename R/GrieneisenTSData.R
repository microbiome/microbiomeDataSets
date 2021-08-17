#' Retrieve GrieneisenTS data
#'
#' Obtain gut microbiome in Papio cynocephalus from Grieneisen et al. (2021).
#'
#' @details
#' The GrieneisenTS dataset contains 16,234 16S rRNA gene sequencing-based microbiome profiles from 585 baboon
#' samples collected over 14 years.Each baboon had an average of 28 samples collected over 4.5 years to 
#' determine whether the heritability of the gut microbiome was dependent on various environmental factors such as diet, age, season.
#' 
#' Column metadata includes the sample ID, baboon ID, collection date, sex, age, social group,
#' group size, rain per month(mm), season(wet/dry), hydro year, month, readcount, plate, post PCR DNA(ng), 
#' diets, ASV shannon H, diet shannon h, ASV richness, principal coordinates(PC).
#' 
#' Row metadata of the microbiome data contains taxonomic information on the
#' Domain, Phylum, Class, Order, Family, Genus, and ASV levels.
#'
#' The row tree consists of a phylogenetic tree build using sequence information of 613 taxa. 
#'    
#' As reference sequences the ASV are provided.
#'
#' @return for  \code{GrieneisenTS} a \linkS4class{TreeSummarizedExperiment}
#'   is provided with the TreeSummarizedExperiment object.
#'
#' @author Laura Grieneisen
#'
#' @references
#' Grieneisen et al. (2021) : Gut microbiome heritability is nearly universal but
#' environmentally contingent 
#' \emph{Science}
#' 373:6551 \url{https://science.sciencemag.org/content/373/6551/181.full}
#'
#' @name GrieneisenTSData
#' @importFrom TreeSummarizedExperiment  TreeSummarizedExperiment<-
#' @export
#'
#' @example
#' # tse <- GrieneisenTSData()
GrieneisenTSData <- function() {
  dataset <- "3.13/grienesen-ts"
  tse <- .create_tse(dataset,
                     assays = "counts",
                     has.rowdata = TRUE,
                     has.coldata = TRUE,
                     has.rowtree = TRUE,
                     has.refseq = TRUE,
                     prefix = NULL)
  tse
}