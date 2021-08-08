
library(tidyverse)
library(remotes)
library(mia)
library(BiocManager)
library(SummarizedExperiment)
library(biomformat)
library(TreeSummarizedExperiment)


#Sample phenotype data
samples <- read.csv("git_metadata_and_community_phenotypes.csv")

#Abundance table
counts <- readRDS("git_ASV_table.RDS")

#Taxonomic mapping table
tax <- counts %>% select(3:8)

# In order to make sure that the samples match between the abundance and phenodata tables
counts_trimmed <- counts[, samples$sample]


se <- SummarizedExperiment(assays = list(counts = counts_trimmed),
                           colData = samples,
                           rowData = tax)

tse <- as(se, "TreeSummarizedExperiment")

#Phylogenetic tree
tree <- ape::read.tree("philr_tree_139_asv.nwk")


tse <- TreeSummarizedExperiment(assays = list(counts = counts_trimmed),
                           colData = samples,
                           rowData = tax,
                           rowTree = tree)

#sequence information
refSeq <- DNAStringSet(counts$asv_sequence, start=NA, end=NA, width=NA, use.names=TRUE)

referenceSeq(tse) <- refSeq



