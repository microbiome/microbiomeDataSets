#' Retrieve GrieneisenTS data
#'
#' Obtain longitudinal gut microbiome data in wild baboons from
#' Grieneisen et al. (2021).
#'
#' @details
#' The GrieneisenTS dataset contains 16,234 16S rRNA gene
#' sequencing-based microbiome profiles from 585 baboon samples
#' collected over 14 years to determine the heritability of the
#' gut microbiome on various environmental factors such as
#' diet, age, season. Each baboon had an average of 28 samples
#' collected over 4.5 years. The data set can be used to
#' investigate significance of longitudinal sampling at large
#' sample sizes.
#'
#' This data set contains the 613 most prevalent taxa with a 
#' phylogenetic tree.
#'
#' Column metadata includes the following fields:
#'\describe{
#'  \item{sample}{Sample ID (character)}
#'  \item{baboon_id}{Baboon ID (factor)}
#'  \item{collection_date}{Sample collection date (date; YYYY-MM-DD)}
#'  \item{sex}{Sex (factor; F/M)}
#'  \item{age}{Age (numeric)}
#'  \item{social_group}{Social group ID (factor)}
#'  \item{group_size}{Social group size (integer)}
#'  \item{rain_month_mm}{Rain per month(mm) (numeric)}
#'  \item{season}{Season (factor; dry/wet)}
#'  \item{hydro_year}{Hydro year (integer)}
#'  \item{month}{Month (integer)}
#'  \item{readcount}{Read count (numeric)}
#'  \item{plate}{Plate (factor)}
#'  \item{post_pcr_dna_ng}{Post PCR DNA(ng) (numeric)}
#'  \item{diet_PC1}{Diet Principal coordinate 1 (numeric)}
#'  \item{diet_PC2}{Diet Principal coordinate 2 (numeric)}
#'  \item{diet_PC3}{Diet Principal coordinate 3 (numeric)}
#'  \item{diet_PC4}{Diet Principal coordinate 4 (numeric)}
#'  \item{diet_PC5}{Diet Principal coordinate 5 (numeric)}
#'  \item{diet_PC6}{Diet Principal coordinate 6 (numeric)}
#'  \item{diet_PC7}{Diet Principal coordinate 7 (numeric)}
#'  \item{diet_PC8}{Diet Principal coordinate 8 (numeric)}
#'  \item{diet_PC9}{Diet Principal coordinate 9 (numeric)}
#'  \item{diet_PC10}{Diet Principal coordinate 10 (numeric)}
#'  \item{diet_PC11}{Diet Principal coordinate 11 (numeric)}
#'  \item{diet_PC12}{Diet Principal coordinate 12 (numeric)}
#'  \item{diet_PC13}{Diet Principal coordinate 13 (numeric)}
#'  \item{diet_shannon_h}{Dietary Shannon's H index (numeric)}
#'  \item{asv_richness}{Amplicon sequence variant (ASV) richness (integer)}
#'  \item{asv_shannon_h}{ASV Shannon's H index (numeric)}
#'  \item{pc1_bc}{Principal coordinate 1 Bray-Curtis dissimilarity (numeric)}
#'  \item{pc2_bc}{Principal coordinate 2 Bray-Curtis dissimilarity (numeric)}
#'  \item{pc3_bc}{Principal coordinate 3 Bray-Curtis dissimilarity (numeric)}
#'  \item{pc4_bc}{Principal coordinate 4 Bray-Curtis dissimilarity (numeric)}
#'  \item{pc5_bc}{Principal coordinate 5 Bray-Curtis dissimilarity (numeric)}
#'}
#'
#' Row metadata of the microbiome data contains taxonomic information on the
#' Domain, Phylum, Class, Order, Family, Genus, and ASV levels.
#'
#' The row tree consists of a phylogenetic tree build using sequence
#' information of 613 taxa.
#'
#' As reference sequences the ASV are provided.
#'
#' @return A \linkS4class{TreeSummarizedExperiment} object.
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
#' @aliases baboongut
#' @export
#'
#' @examples
#' tse <- GrieneisenTSData()
#'
GrieneisenTSData <- function() {
    dataset <- "3.14/grieneisen-ts"
    tse <- .create_tse(dataset,
                    assays = "counts",
                    has.rowdata = TRUE,
                    has.coldata = TRUE,
                    has.rowtree = TRUE,
                    has.refseq = TRUE,
                    prefix = NULL)
    tse
}

#' @rdname GrieneisenTSData
#' @export
baboongut <- GrieneisenTSData
