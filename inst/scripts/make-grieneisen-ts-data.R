library(tidyverse)
library(remotes)
library(mia)
library(BiocManager)
library(SummarizedExperiment)
library(biomformat)
library(TreeSummarizedExperiment)
library(dplyr)
library(utils)
library(tidyr)
library(readxl)

#Sample phenotype data
samples <- read.csv("git_metadata_and_community_phenotypes.csv")

#Abundance table
counts <- readRDS("git_ASV_table.RDS")

#Taxonomic mapping table
tax <- select(counts, domain, phylum, class, order, family, genus)

# In order to make sure that the samples match between the abundance and phenodata tables
counts_trimmed <- counts[, samples$sample]

#Phylogenetic tree
tree <- ape::read.tree("philr_tree_139_asv.nwk")

tse <- TreeSummarizedExperiment(assays = list(counts = counts_trimmed),
                                colData = samples,
                                rowData = tax,
                                rowTree = tree)
#sequence information
refSeq <- DNAStringSet(counts$asv_sequence, start=NA, end=NA, width=NA, use.names=TRUE)

referenceSeq(tse) <- refSeq



#Abundance Table
counts <- readRDS("git_ASV_table.RDS")
counts_1 <- counts  %>% dplyr::rename(c(
  "ASV ID" = "asv_id",
  "ASV Sequence" = "asv_sequence"))

#converting csv file to rds file
saveRDS(samples, "samples-names.rds")

#Metadata Table
samples_1 <- readRDS("samples-names.rds")
samples_2 <- samples_1 %>% dplyr::rename(c(
  "Number" = "X",
  "Sample ID" = "sample",
  "baboon ID" = "baboon_id",
  "date" = "collection_date",
  "group" = "social_group" ,
  "group size" = "group_size",
  "rain/month (mm) " = "rain_month_mm",
  "hydro year" = "hydro_year",
  "post PCR DNA (ng)" = "post_pcr_dna_ng",
  "diet Shannon h " = "diet_shannon_h",
  "ASV richness" = "asv_richness",
  "ASV Shannon h" = "asv_shannon_h"))


data <- "../MicrobiomeDataSets/3.14/"


saveRDS(samples_2, file = paste0(data, "samples.rds"))
saveRDS(counts_1, file = paste0(data, "counts.rds"))
saveRDS(tree, file = paste0(data, "tree.rds"))
saveRDS(tax, file = paste0(data,"tax.rds"))



