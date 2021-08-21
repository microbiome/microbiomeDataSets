#' Retrieve GrieneisenTS data
#'
#' Obtain gut microbiome in Papio cynocephalus from Grieneisen et al. (2021).
#'
#' @details
#' The GrieneisenTS dataset contains 16,234 16S rRNA gene
#' sequencing-based microbiome profiles from 585 baboon samples
#' collected over 14 years to determine the heritability of the 
#' gut microbiome on various environmental factors such as 
#' diet(high/low diversity), age, season(wet/dry).
#' Each baboon had an average of 28 samples collected over 4.5 years.
#' The data set can be used to investigate significance 
#' of longitudinal sampling at large sample sizes.
#'
#' This data set contains the 613 most prevalent taxa, including the
#' phylogenetic tree.
#' 
#' Column metadata includes the samples, baboon ID, collection date of 
#' the samples, sex(F/M), age, social group, social group size, 
#' rain per month(mm), season(wet/dry), hydro year, month of the year,  
#' readcount information, plate information, post PCR DNA(ng) information,  
#' diets, ASV Shannon's H index, Dietary Shannon's H, Amplicon sequence
#' variant (ASV) richness, the five first principal coordinates(PCs) of a
#' Bray-Curtis dissimilarity matrix.
#' 
#' Row metadata of the microbiome data contains taxonomic information on the
#' Domain, Phylum, Class, Order, Family, Genus, and ASV levels.
#'
#' The row tree consists of a phylogenetic tree build using sequence
#' information of 613 taxa. 
#'    
#' As reference sequences the ASV are provided.
#'
#' @return for  \code{GrieneisenTS} a \linkS4class{TreeSummarizedExperiment}
#'   is provided with the TreeSummarizedExperiment object.
#'
#' @author Yagmur Simsek and Leo Lahti
#'
#' @references
#' Grieneisen et al. (2021):
#' Gut microbiome heritability is nearly universal but
#' environmentally contingent 
#' \emph{Science}
#' 373:6551 \url{https://science.sciencemag.org/content/373/6551/181.full}
#'
#' @name GrieneisenTSData
#' @export
#'
#' @example
#' \donttest{
#'   tse <- GrieneisenTSData()
#' }
GrieneisenTSData <- function() {
    dataset <- "3.14/grienesen-ts"
    tse <- .create_tse(dataset,
                    assays = "counts",
                    has.rowdata = TRUE,
                    has.coldata = TRUE,
                    has.rowtree = TRUE,
                    has.refseq = TRUE,
                    prefix = NULL)
    tse
}
