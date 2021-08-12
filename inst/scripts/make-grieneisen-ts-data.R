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
rownames(samples) <- samples$sample
samples$sample <- NULL

#Abundance table
counts <- readRDS("git_ASV_table.RDS")

#sequence information
ASV_Sequence <- counts$asv_sequence
ref_seq <- data.frame(ASV_Sequence)
rownames(ref_seq) <- counts$asv_id


#Taxonomic mapping table
tax_1 <- select(counts, domain, phylum, class, order, family, genus, asv_id)
rownames(tax_1) <- counts$asv_id

# In order to make sure that the samples match between the abundance and phenodata tables
counts_trimmed <- counts[, samples$sample]

#Phylogenetic tree
tree <- ape::read.tree("philr_tree_139_asv.nwk")


#converting csv file into rds file
saveRDS(samples, "samples.rds")

#Metadata Table
samples_1 <- samples %>% dplyr::rename(c(
  Number = "X",
  Baboon = "baboon_id",
  Date = "collection_date",
  Sex = "sex",
  Age = "age",
  Group = "social_group" ,
  "Group_Size" = "group_size",
  "Rain_Month_mm" = "rain_month_mm",
  Season = "season",
  Month = "month",
  Readcount = "readcount",
  Plate = "plate",
  "Hydro_Year" = "hydro_year",
  "Post_PCR_DNA_ng" = "post_pcr_dna_ng",
  "Diet_Shannon_h" = "diet_shannon_h",
  "ASV_Richness" = "asv_richness",
  "ASV_Shannon_h" = "asv_shannon_h",
  PC1_bc = "pc1_bc",
  PC2_bc = "pc2_bc",
  PC3_bc = "pc3_bc",
  PC4_bc = "pc4_bc",
  PC5_bc = "pc5_bc",
  ))

samples <- data.frame(samples_1)

#Taxonomy Table
saveRDS(tax_1, "tax.rds")
tax <- tax_1  %>% dplyr::rename(c(
  Domain = "domain",
  Phylum = "phylum",
  Class = "class",
  Order = "order",
  Family = "family",
  Genus = "genus",
  ASV = "asv_id"
  ))

#Arranging counts
rownames(counts) <- counts$asv_id
counts[, c("asv_id", "asv_sequence", "domain", "phylum", "class", "order", "family", "genus")]  <- NULL

#sequence information
refSeq <- DNAStringSet(counts$asv_sequence, start=NA, end=NA, width=NA, use.names=TRUE)


tse <- TreeSummarizedExperiment(assays = list(counts = counts_trimmed),
                                colData = samples_1,
                                rowData = tax,
                                rowTree = tree)
referenceSeq(tse) <- refSeq

data <- "/MicrobiomeDataSets/3.14/"

saveRDS(samples_1, file = paste0(data, "samples.rds"))
saveRDS(counts, file = paste0(data, "counts.rds"))
saveRDS(tree,  file = paste0(data,"tree.rds"))
saveRDS(tax, file = paste0(data, "tax.rds"))
saveRDS(ref_seq, file = paste0(data,"sequence.rds"))





