
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom MultiAssayExperiment MultiAssayExperiment
.create_mae <- function(dataset,
                        types = list(),
                        hub = ExperimentHub(),
                        coldata = TRUE,
                        samplemap = TRUE,
                        has.rowdata = list(),
                        has.coldata = list()){

    el <- .get_experiment_list(dataset, hub, types, has.rowdata, has.coldata)
    args <- .get_col_row_map_data(dataset, hub,
                                  has.rowdata = FALSE,
                                  has.coldata = coldata,
                                  has.samplemap = samplemap)
    do.call(MultiAssayExperiment, c(list(el), args))
}

#' @importFrom MultiAssayExperiment ExperimentList
.get_experiment_list <- function(dataset, hub,
                                 types = list(),
                                 has.rowdata = list(),
                                 has.coldata = list()){
    experiments <- mapply(
        function(prefix, se, hr, hc){
            stopifnot(length(se) == 1L)
            FUN <- switch(names(se),
                          TSE = .create_tse,
                          SE = .create_se)
            do.call(FUN, list(dataset = dataset,
                              hub = hub,
                              assays = se[[1L]],
                              has.rowdata = hr,
                              has.coldata = hc,
                              prefix = prefix))
        },
        names(types),
        types,
        has.rowdata,
        has.coldata)
    names(experiments) <- names(types)
    do.call(ExperimentList, experiments)
}

#' @importFrom ExperimentHub ExperimentHub
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
.create_tse <- function(dataset,
                        hub = ExperimentHub(),
                        assays = "counts",
                        has.rowdata = TRUE,
                        has.coldata = TRUE,
                        prefix = NULL) {
    assays <- .get_assays(dataset, hub, assays, prefix)
    args <- .get_col_row_map_data(dataset, hub, prefix,
                                  has.rowdata = has.rowdata,
                                  has.coldata = has.coldata)
    do.call(TreeSummarizedExperiment, c(list(assays=assays), args))
}
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom SummarizedExperiment SummarizedExperiment
.create_se <- function(dataset,
                       hub = ExperimentHub(),
                       assays = "counts",
                       has.rowdata = TRUE,
                       has.coldata = TRUE,
                       prefix = NULL) {
    assays <- .get_assays(dataset, hub, assays, prefix)
    args <- .get_col_row_map_data(dataset, hub, prefix,
                                  has.rowdata = has.rowdata,
                                  has.coldata = has.coldata)
    do.call(SummarizedExperiment, c(list(assays=assays), args))
}

.get_base_path <- function(dataset){
    paste0(file.path(system.file(package = "microbiomeDataSets",mustWork = FALSE),
                     "extdata","hub",dataset))
    # file.path("microbiomeDataSet", dataset)
}

.norm_prefix <- function(prefix){
    if (is.null(prefix)) {
        prefix <- ""
    } else {
        prefix <- paste0(prefix, "_")
    }
    prefix
}

.get_assays <- function(dataset, hub, assays, prefix = NULL){
    base <- .get_base_path(dataset)
    prefix <- .norm_prefix(prefix)
    assay_list <- list()
    for (a in assays) {
        assay_list[[a]] <- readRDS(file.path(base, sprintf("%s%s.rds", prefix, a)))
        # assay_list[[a]] <- hub[hub$rdatapath==file.path(base, sprintf("%s%s.rds", prefix, a))][[1]]
    }
    names(assay_list) <- assays
    assay_list
}

.get_col_row_map_data <- function(dataset, hub,
                                  prefix = NULL,
                                  has.rowdata = TRUE,
                                  has.coldata = TRUE,
                                  has.samplemap = FALSE){
    base <- .get_base_path(dataset)
    prefix <- .norm_prefix(prefix)
    args <- list()
    if (has.coldata) {
        args$colData <- readRDS(file.path(base, sprintf("%scoldata.rds", prefix)))
        # args$colData <- hub[hub$rdatapath==file.path(base, sprintf("%scoldata.rds", prefix))][[1]]
    }
    if (has.rowdata) {
        args$rowData <- readRDS(file.path(base, sprintf("%srowdata.rds", prefix)))
        # args$rowData <- hub[hub$rdatapath==file.path(base, sprintf("%srowdata.rds", prefix))][[1]]
    }
    if (has.samplemap) {
        args$sampleMap <- readRDS(file.path(base, sprintf("%ssamplemap.rds", prefix)))
        # args$sampleMap <- hub[hub$rdatapath==file.path(base, sprintf("%ssamplemap.rds", prefix))][[1]]
    }
    args
}
