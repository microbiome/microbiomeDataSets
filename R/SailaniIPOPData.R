#' Obtain the Sailani deep longitudinal multiomics data
#'
#' Obtain the deep longitudinal multiomics data from Sailani et al. (2020).
#'
#' @details
#' The SailaniIPOP dataset contains longitudinal multiomics data from profiling 
#' of 105 generally healthy individuals. This dataset consists of 16S rRNA gene 
#' based profiling of gut and nasal microbiota. Host transcriptome, proteome and 
#' metabolome from plasma, cytokine profile, and associated clinical data. The
#' samples were collected over a period of 4 years. The data set can be used 
#' to investigate longitudinal dynamics and interaction between microbiota and 
#' host associated factors.  
#'
#' This microbiota data as currently available as relative abundances. Hierarchical
#' taxonomic information is unavailable for microbiota data.
#'
#' Column metadata includes the time points, sex, subject identifier, sample
#' identifier and treatment group.
#'
#'
#' All data are downloaded from ExperimentHub and cached for local re-use.
#'
#' @return for  \code{SailaniIPOPData} a \linkS4class{MultiAssayExperiment} object
#'   consisting of two \linkS4class{TreeSummarizedExperiment} objects and 
#'   five \linkS4class{SummarizedExperiment} objects.
#'
#' @author Sudarshan Shetty and Felix G.M. Ernst
#'
#' @references
#' Sailani, M.R., Metwally, A.A., Zhou, W. et al. (2020) Deep longitudinal multiomics 
#' profiling reveals two biological seasonal patterns in California.
#' \emph{Nat Commun.} 11, 4933 \url{https://doi.org/10.1038/s41467-020-18758-1}
#'
#' Zhou, W., Sailani, M.R., Contrepois, K. et al. (2019) Longitudinal multi-omics of 
#' host-microbe dynamics in prediabetes.
#' \emph{Nature} 569, 663-671 \url{https://doi.org/10.1038/s41586-019-1236-x}
#'
#' @name SailaniIPOPData
#' @export
#'
#' @examples
#' mae <- SailaniIPOPData()
#' #tse <- SailaniIPOPData()
SailaniIPOPData <- function() {
    mae <- .create_mae("sailani-ipop",
                       types = list(gutMicrobiota = list("TSE" = c("relabundance")),
                                    nasalMicrobiota = list("TSE" = c("relabundance")),
                                    hostTranscriptome = list("SE" = c("counts")),
                                    hostProteome = list("SE" = c("counts")),
                                    hostMetabolome = list("SE" = c("counts")),
                                    hostCytokines = list("SE" = c("counts")),
                                    hostClinical = list("SE" = c("counts"))),
                       coldata = FALSE,
                       samplemap = FALSE,
                       has.rowdata = list(gutMicrobiota = FALSE,
                                          nasalMicrobiota = FALSE,
                                          hostTranscriptome = FALSE,
                                          hostProteome = FALSE,
                                          hostMetabolome = FALSE,
                                          hostCytokines = FALSE,
                                          hostClinical = FALSE),
                       has.coldata = list(gutMicrobiota = TRUE,
                                          nasalMicrobiota = TRUE,
                                          hostTranscriptome = TRUE,
                                          hostProteome = TRUE,
                                          hostMetabolome = TRUE,
                                          hostCytokines = TRUE,
                                          hostClinical = TRUE))
    mae
}


# #' @rdname SailaniIPOPData
# #' @importFrom SummarizedExperiment colData<-
# #' @export

SailaniIPOPData <- function() {
    dataset <- "sailani-ipop"
    hub <- ExperimentHub()
    tse <- .create_tse(dataset,
                       hub = hub,
                       assays = c("relabundance"),
                       has.rowdata = FALSE,
                       has.coldata = TRUE) #"$Microbiota"
    args <- .get_col_row_map_data(dataset,
                                  hub = hub,
                                  has.rowdata = FALSE,
                                  has.coldata = TRUE)
    colData(tse) <- args$colData
    tse
}

#' @rdname SailaniIPOPData
#' @export
sailaniipop <- SailaniIPOPData
