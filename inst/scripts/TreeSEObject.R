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


data <- "/Users/Yagmur Simsek/Documents/microbiomeDataSets/inst/extdata/3.14/grieneisen-ts/"

saveRDS(tax, file = paste0(data, "tax.RDS"))
saveRDS(counts, file = paste0(data, "counts.rds"))
saveRDS(samples, file = paste0(data, "samples.rds"))
saveRDS(counts_trimmed, file = paste0(data, "counts_trimmed.rds"))
saveRDS(tree, file = paste0(data, "tree.rds"))

#Metadata Table
samples <- readRDS(paste0(data, "samples.RDS")) %>%
  rename( number = "X") %>%
  rename( " baboon ID" = "baboon_id") %>%
  rename( date = "collection_date") %>%
  rename( group = "social_group") %>%
  rename( "group size" = "group_size") %>%
  rename( "rain/month (mm) " = "rain_month_mm") %>%
  rename( "hydro year" = "hydro_year") %>%
  rename( "post PCR DNA (ng)" = "post_pcr_dna_ng") %>%
  rename( "diet Shannon h " = "diet_shannon_h") %>%
  rename( "ASV richness" = "asv_richness") %>%
  rename( "ASV Shannon h" = "asv_shannon_h" )

# Abundance Table
counts <- readRDS(paste0(data, "counts.rds")) %>%
  rename( "ASV ID" = "asv_id") %>%
  rename( "ASV Sequence" = "asv_sequence" )


