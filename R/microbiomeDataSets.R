#' Microbiome data sets 
#' 
#' \code{microbiomeDataSets} is a collation of data from microbiome and 
#' associated studies, which are publicaly available. 
#' 
#' The data is made available through the ExperimentHub resources of the
#' Bioconductor project. It is loaded as \code{TreeSummarizedExperiment} 
#' object or a \code{MultiAssayExperiment} objects.
#' 
#' @section Deprecated functions:
#' The following functions are deprecated and will be removed in a future 
#' release:
#' \itemize{
#'   \item \code{LahtiWAdata}: Please use the \code{data()} function from the 
#'   \code{mia} package to load this dataset.
#'   \item \code{HintikkaXOData}: Please use the \code{data()} function from 
#'   the \code{miaViz} package to load this dataset.
#' }
#'
#' @name microbiomeDataSets-package
NULL

#' @importFrom BiocGenerics updateObject
#' @importFrom methods as
NULL

# Wrapper for LahtiWAdata
LahtiWAdata <- function() {
    .Deprecated("This function is deprecated and will be removed in a future 
                release. Please load the dataset directly using the `data()` 
                function from the mia/miaViz/mitTime packages.")
    data_env <- new.env()
    data("LahtiWAdata", package = "mia", envir = data_env)
    return(data_env$LahtiWAdata)
}

# Wrapper for HintikkaXOData
HintikkaXOData <- function() {
    .Deprecated("This function is deprecated and will be removed in a future 
                release. Please load the dataset directly using the `data()` 
                function from the mia/miaViz/mitTime packages.")
    data_env <- new.env()
    data("HintikkaXOData", package = "miaViz", envir = data_env)
    return(data_env$HintikkaXOData)
}

