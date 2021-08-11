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
#' Column metadata is common for all experiments (microbiome, metabolites,
#' biomarkers) and includes the Sample (ID), Rat (ID), Site of measurement,
#' Diet group, Fat in Diet (Low/High), XOS Diet Supplement (0/1)
#' 
#' Row metadata of the microbiome data contains taxonomic information on the
#' Phylum, Class, Order, Family, Genus, Species, and OTU levels.
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
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment rowData<-
#' 
#' @name HintikkaXOData
#' @export
#'
#' @examples
#' # Retrieve the MAE data
#' mae <- HintikkaXOData()
#'
#' # List the experiments in this MultiAssayExperiment object
#' print(experiments(mae))
#'
#' # colData for this MAE data object
#' colData(mae)
#'
#' # metabolite assay data
#' # assays(mae[["metabolites"]])$assay
#'
#' # biomarker assay data
#' # assays(mae[["biomarkers"]])$assay
# '
#' # microbiome assay counts
#' # assays(mae[["microbiome"]])$counts
#'
#' # microbiome rowData
#' rowData(mae[["microbiome"]])
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


    # Fix order in rowData
    rowData(mae[["microbiome"]]) <- rowData(mae[["microbiome"]])[, c(2:7, 1)]
    # Add sample names to rownames
    rownames(mae[["microbiome"]]) <- rowData(mae[["microbiome"]])$OTU

    mae
    
}

