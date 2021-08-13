#' Retrieve HintikkaXO data
#'
#' Obtain microbiota, metabolite and biomarker data from Hintikka et al. (2021).
#'
#' @details
#' The HintikkaXO dataset contains high-throughput profiling data from 40 rat
#' samples, including 39 biomarkers, 38 metabolites (NMR), and 12706 OTUs from
#' 318 species, measured from Cecum. This is diet comparison study with
#' High/Low fat diet and xylo-oligosaccaride supplementation.
#'
#' Column metadata is common for all experiments (microbiota, metabolites,
#' biomarkers) and includes the following fields:
#'
#' \itemize{
#'   \item{Sample: } {Sample ID (character)}
#'   \item{Rat: } {Rat ID (factor)}
#'   \item{Site: } {Site of measurement ("Cecum"); single value}
#'   \item{Diet: } {Diet group (factor; combination of the Fat and XOS fields)}
#'   \item{Fat: } {Fat in Diet (factor; Low/High)}
#'   \item{XOS: } {XOS Diet Supplement (numeric; 0/1)}
#' }
#' 
#' Row metadata of the microbiota data contains taxonomic information on the
#' Phylum, Class, Order, Family, Genus, Species, and OTU levels.
#'
#' Biomarker data contains 39 biomarkers.
#'
#' Metabolite data contains 38 NMR metabolites.
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
#' metabolomic patterns in fecal microbiotas with biclustering.  
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
#' # colData for this MAE data object (print first few rows)
#' print(head(colData(mae)))
#'
#' # metabolite assay data
#' nmr <- assays(mae[["metabolites"]])$nmr
#'
#' # biomarker assay data
#' bm <- assays(mae[["biomarkers"]])$signals
# '
#' # microbiota assay counts
#' counts <- assays(mae[["microbiota"]])$counts
#'
#' # microbiota rowData
#' taxtab <- rowData(mae[["microbiota"]])
#' 
HintikkaXOData <- function() {

    mae <- .create_mae("3.14/hintikka-xo",
    
        types = list(microbiota  = list("SE" = c("counts")),
                    metabolites = list("SE"  = nmr),
                    biomarkers  = list("SE"  = signals)
                ),
 
        coldata = TRUE,
        samplemap = FALSE,

        has.rowdata = list(microbiota = TRUE,
                        metabolites = FALSE,
                        biomarkers  = FALSE),

        has.coldata = list(microbiota = FALSE,
                        metabolites = FALSE,
                        biomarkers  = FALSE))

    mae
    
}

