#' Retrieve HintikkaXO data
#'
#' Obtain microbiome and lipid data from Hintikka et al. (2021).
#'
#' @details
#' The HintikkaXO dataset contains high-throughput profiling data from 40 rat
#' samples, including 39 biomarkers, 38 metabolites (NMR), and 12706 OTUs from
#' 318 species, measured from Cecum. This is diet comparison study with
#' High/Low fat diet and xylo-oligosaccaride supplementation.
#'
#' Column metadata includes the SampleID, RatID, Measurement site, Diet group,
#' Diet Fat (Low/High),  Diet Supplement (XOS 0/1)
#' 
#' Row metadata of the microbiome data contains taxonomic information on the
#' Phylum, Class, Order, Family, Genus, Species, and OTU levels.
#'
#' Biomarker and NMR metabolite data are provided as altExp (alternative
#' Experiments).
#' 
#' Biomarker data contains 39 biomarkers.
#'
#' NMR data contains 38 metabolites.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#'
#' @return for  \code{HintikkaXOData} a \linkS4class{SummarizedExperiment}
#'   is provided with the altExp objects.
#'
#' @author Leo Lahti
#'
#' @references
#' Hintikka L et al. (2021): Xylo-oligosaccharides in prevention of hepatic
#' steatosis and adipose tissue inflammation: associating taxonomic and
#' metabolomic patterns in fecal microbiomes with biclustering.  
#' \emph{International Journal of Environmental Research and Public Health}
#' 18(8):4049 \url{https://doi.org/10.3390/ijerph18084049}
#'
#' @name HintikkaXOData
#' @importFrom SingleCellExperiment altExp
#' @importFrom SingleCellExperiment altExp<-
#' @export
#'
#' @examples
#' # tse <- HintikkaXOData()
#'
HintikkaXOData <- function() {

    # Circumvent warnings
    meta_cecum <- otu_cecum <- tax <- NULL

    dataset <- "3.14/hintikka-xo"
    se <- .create_se(dataset,
                    assays = c("counts"),
                    has.rowdata = TRUE,
                    has.coldata = TRUE)

    tse <- SummarizedExperiment(assays = list(counts = otu_cecum),
                        colData = meta_cecum,
                        rowData = tax)

    # Create SE objects
    bm  <- SummarizedExperiment(bm)
    nmr <- SummarizedExperiment(nmr)

    # Add Biomarkers as "alternative experiment":
    altExp(se, "Biomarkers") <- bm

    # Add NMR metabolite abundances as "alternative experiment":
    altExp(se, "NMR") <- nmr

    se

}


