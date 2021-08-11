#' Retrieve HintikkaXO data
#'
#' Obtain microbiome, metabolite and biomarker data from Hintikka et al. (2021).
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
#' Biomarker and NMR metabolite data are provided through the
#' \linkS4class{MultiAssayExperiment} mechanism.
#'
#' Biomarker data contains 39 biomarkers.
#'
#' NMR data contains 38 metabolites.
#' 
#' All data are downloaded from ExperimentHub and cached for local re-use.
#'
#' @return for  \code{HintikkaXOData} a \linkS4class{MultiAssayExperiment}
#' is provided.
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
#' @export
#'
#' @examples
#' # tse <- HintikkaXOData()
#'
HintikkaXOData <- function() {

    mae <- .create_mae("3.14/hintikka-xo",
    
        types = list(microbiome  = list("TSE" = c("counts")),
                     metabolites = list("SE"  = NULL),
                     biomarkers  = list("SE"  = NULL)
             	),
 
        coldata = TRUE,
        samplemap = FALSE,
		    
        has.rowdata = list(microbiome = TRUE,
                           metabolites = FALSE,
	                   biomarkers  = FALSE),
	
        has.coldata = list(microbiome = FALSE,
                          metabolites = FALSE,
		          biomarkers  = FALSE))
    mae
    
}




HintikkaXODataTSE <- function() {

    dataset <- "3.14/hintikka-xo"
    hub <- ExperimentHub()
    tse <- .create_tse(dataset,
                    hub = hub,
                    assays = c("counts"),
                    has.rowdata = TRUE,
                    has.coldata = FALSE,
                    prefix = "microbiome")
		    
    args <- .get_col_row_map_data(dataset,
                                hub = hub,
                                has.rowdata = FALSE,
                                has.coldata = TRUE)

    colData(tse) <- DataFrame(args$colData)
    tse

}


HintikkaXODataSE <- function() {

    dataset <- "3.14/hintikka-xo"
    hub <- ExperimentHub()
    tse <- .create_se(dataset,
                    hub = hub,
                    assays = c("counts"),
                    has.rowdata = TRUE,
                    has.coldata = FALSE,
                    prefix = "microbiome")
		    
    args <- .get_col_row_map_data(dataset,
                                hub = hub,
                                has.rowdata = FALSE,
                                has.coldata = TRUE)

    colData(se) <- DataFrame(args$colData)
    se

}