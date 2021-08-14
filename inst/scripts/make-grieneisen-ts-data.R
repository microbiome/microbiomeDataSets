library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)

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
counts_trimmed <- counts[, rownames(samples)]

#Phylogenetic tree
pseq <- readRDS("asv_10prev_phyloseq_for_lahti.RDS")
tree <- phy_tree(pseq)

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
refSeq <- DNAStringSet(ASV_Sequence, start=NA, end=NA, width=NA, use.names=TRUE)

tse <- TreeSummarizedExperiment(assays = list(counts = counts_trimmed),
                                colData = samples_1,
                                rowData = tax,
                                referenceSeq = refSeq)

#Subsetting
tse <- tse[tree$tip.label, ]

rowTree(tse) <- tree
tse 

rowdata <- rowData(tse)
coldata <- colData(tse)
counts <- assay(tse, "counts")
rowtree <- rowTree(tse)
refseq <- referenceSeq(tse)

unlink(c("samples.rds","tax.rds"))

data <- "/MicrobiomeDataSets/3.14/"

saveRDS(counts, file = paste0(data,"counts.rds"))
saveRDS(coldata, file = paste0(data,"coldata.rds"))
saveRDS(rowdata , file = paste0(data, "rowdata.rds"))

ape::write.tree(rowtree, file= paste0(data,"rowtree.tre"))
gz_tree<- gzfile("rowtree.tre.gz", "w")
writeLines(readLines("rowtree.tre"), con=gz_tree)
close(gz_tree)
unlink("rowtree.tre")

Biostrings::writeXStringSet(refseq, "refseq.fasta.gz", compress = TRUE)
