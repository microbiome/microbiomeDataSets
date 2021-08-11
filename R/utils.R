#' Load available microbiome data sets
#' 
#' To list the available datasets in \code{microbiomeDataSets},
#' run \code{availableDataSets}.
#' 
#' For information visit the individual man pages.
#' 
#' @return 
#' A \code{data.frame} containing the following columns:
#' \itemize{
#'   \item{Dataset:} {the name of the function to load a dataset}
#' }
#' 
#' @export
#' 
#' @importFrom utils read.csv
#' 
#' @examples
#' availableDataSets()
availableDataSets <- function(){
    file <- system.file("extdata","datasets.csv",package = "microbiomeDataSets")
    read.csv(file, stringsAsFactors = FALSE)
}

####################################################################
# ExperimentHub data loading

.get_res_by_path <- function(hub, path, failOnError = TRUE){
    res <- hub[hub$rdatapath==path]
    if(length(res) > 1L){
        stop("Multiple resources found.")
    }
    if(length(res) == 0L && failOnError){
        stop("No resources found for ",path,".")
    } else if(length(res) == 0L) {
        return(NULL)
    }
    res[[1L]]
}

###################################################################
# object creation 

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
                        has.rowtree = FALSE,
                        has.coltree = FALSE,
                        has.refseq = FALSE,
                        prefix = NULL) {
    assays <- .get_assays(dataset, hub, assays, prefix)
    args <- .get_col_row_map_data(dataset, hub, prefix,
                                  has.rowdata = has.rowdata,
                                  has.coldata = has.coldata)
    tse <- do.call(TreeSummarizedExperiment, c(list(assays=assays), args))
    tse <- .add_refseq(dataset, hub, prefix = prefix,
                       tse, has.refseq)
    tse <- .add_trees(dataset, hub, prefix = prefix,
                      tse, has.rowtree, has.coltree)
    tse
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

#########################################################################
# path utils

.get_base_path <- function(dataset){
    # paste0(file.path(system.file(package = "microbiomeDataSets",
    #     mustWork = FALSE),
    #                  "extdata","hub",dataset))
    file.path("microbiomeDataSets", dataset)
}

.norm_prefix <- function(prefix){

    if (is.null(prefix)) {
        prefix <- ""
    } else {
        prefix <- paste0(prefix, "_")
    }
    prefix
}

################################################################################
# assay data

.get_assays <- function(dataset, hub, assays, prefix = NULL){
    base <- .get_base_path(dataset)
    # With multiple assays, add "_" to prefix
    if (!is.null(assays)) {
        prefix <- .norm_prefix(prefix)
    } else {
        prefix <- prefix
    }
    assay_list <- list()
    if (!is.null(assays)) {
        for (a in assays) {
            # assay_list[[a]] <- readRDS(file.path(base, sprintf("%s%s.rds",
            #  prefix, a)))
            path <- file.path(base, sprintf("%s%s.rds", prefix, a))
	    assay_list[[a]] <- .get_res_by_path(hub, path)	    
        }
	names(assay_list) <- assays
     } else { # Needed if no assay has been specified in data resource name
              # e.g. metabolites.rds (vs. microbes_counts.rds)
         f <- sprintf("%s%s.rds", prefix, "")
         path <- file.path(base, f)
	 assay_list[[1]] <- .get_res_by_path(hub, path)
	 names(assay_list) <- "assay"
     }         
    assay_list
}

################################################################################
# row, column and sample map data

.get_col_row_map_data <- function(dataset, hub,
                                  prefix = NULL,
                                  has.rowdata = TRUE,
                                  has.coldata = TRUE,
                                  has.samplemap = FALSE){
    base <- .get_base_path(dataset)
    prefix <- .norm_prefix(prefix)
    args <- list()
    if (has.coldata) {
        # args$colData <- readRDS(file.path(base, sprintf("%scoldata.rds",
    #   prefix)))
        path <- file.path(base, sprintf("%scoldata.rds", prefix))
        args$colData <- .get_res_by_path(hub, path)
    }
    if (has.rowdata) {
        # args$rowData <- readRDS(file.path(base, sprintf("%srowdata.rds",
    #   prefix)))
        path <- file.path(base, sprintf("%srowdata.rds", prefix))
        args$rowData <- .get_res_by_path(hub, path)
    }
    if (has.samplemap) {
        # args$sampleMap <- readRDS(file.path(base, sprintf("%ssamplemap.rds",
    #  prefix)))
        path <- file.path(base, sprintf("%ssamplemap.rds", prefix))
        args$sampleMap <- .get_res_by_path(hub, path)
    }
    args
}

################################################################################
# tree data loading

#' @importFrom ape read.tree
#' @importFrom TreeSummarizedExperiment LinkDataFrame
#' @importClassesFrom TreeSummarizedExperiment LinkDataFrame
.get_tree_data <- function(dataset, hub,
                           prefix = NULL,
                           type = c("row","column")){
    base <- .get_base_path(dataset)
    prefix <- .norm_prefix(prefix)
    name <- switch(type,
                   row = "rowtree",
                   column = "coltree")
    # tree <- read.tree(file.path(base, sprintf("%s%s.tre.gz", prefix, name)))
    # links <- list()
    # if(file.exists(file.path(base, sprintf("%s%s_links.rds", prefix, name)))){
    # links <- readRDS(file.path(base, sprintf("%s%s_links.rds", prefix, name)))
    # }
    tree <- file.path(base, sprintf("%s%s.tre.gz", prefix, name))
    links <- file.path(base, sprintf("%s%s_links.rds", prefix, name))
    tree <- read.tree(.get_res_by_path(hub, tree))
    links <- .get_res_by_path(hub, links, failOnError = FALSE)
    if(!is.null(links)){
        links <- as(links,"LinkDataFrame")
    }
    
    list(tree = tree, links = links)
}

#' @importFrom TreeSummarizedExperiment changeTree
.add_tree <- function(tse, tree_data, type = c("row","column")){
    if(length(tree_data$links) == 0L){
        if(type == "row"){
            tse <- changeTree(tse, rowTree = tree_data$tree)
        } else {
            tse <- changeTree(tse, colTree = tree_data$tree)
        }
    } else {
        if(type == "row"){
            tse@rowLinks <- tree_data$links
            tse@rowTree$phylo <- tree_data$tree
        } else {
            tse@colinks <- tree_data$links
            tse@colTree$phylo <- tree_data$tree
        }
    }
    updateObject(tse)
}

.add_trees <- function(dataset, hub,
                       prefix = NULL,
                       tse,
                       has.rowtree = FALSE,
                       has.coltree = FALSE){
    
    if(has.rowtree){
        tree_data <- .get_tree_data(dataset, hub, prefix = prefix, type = "row")
        tse <- .add_tree(tse, tree_data = tree_data, type = "row")
    }
    if(has.coltree){
        tree_data <- .get_tree_data(dataset, hub, prefix = prefix,
                                    type = "column") 
        tse <- .add_tree(tse, tree_data = tree_data, type = "column")
    }
    tse
}

################################################################################
# reference sequence

#' @importFrom Biostrings readDNAStringSet
#' @importFrom TreeSummarizedExperiment referenceSeq<- 
.add_refseq <- function(dataset, hub, prefix = NULL,
                        tse,
                        has.refseq = FALSE){
    if(has.refseq){
        base <- .get_base_path(dataset)
        prefix <- .norm_prefix(prefix)
        # refSeq <- readDNAStringSet(file.path(base,
        # sprintf("%srefseq.fasta.gz", prefix)))
        refSeq <- file.path(base, sprintf("%srefseq.fasta.gz", prefix))
        refSeq <- .get_res_by_path(hub,refSeq)
        refSeq <- readDNAStringSet(refSeq)
        names <- names(refSeq)
        if(!is.null(names) && all(grepl("_ \\|\\|_",names))){
            groups <- regmatches(names,regexec("(.+)_\\|\\|_.*",names))
            groups <- vapply(groups,"[[",character(1),2L)
            names <- regmatches(names,regexec(".*_\\|\\|_(.+)",names))
            names <- vapply(names,"[[",character(1),2L)
            names(refSeq) <- names
            refSeq <- split(refSeq, factor(groups,unique(groups)))
            names(refSeq) <- unique(groups)
        }
        referenceSeq(tse) <- refSeq
    }
    tse
}
