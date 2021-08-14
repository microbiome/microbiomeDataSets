library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ape)

#Sample phenotype data
samples <- read.csv("git_metadata_and_community_phenotypes.csv")
rownames(samples) <- samples$sample
as.factor(c("baboon_id", "sex", "social_group", "season" ,"plate"))
samples$collection_date <- as.Date(samples$collection_date, format="%Y-%m-%d")

#Abundance table
counts <- readRDS("git_ASV_table.RDS")

#sequence information
ASV_Sequence <- counts$asv_sequence
ref_seq <- data.frame(ASV_Sequence)
rownames(ref_seq) <- counts$asv_id

#Taxonomic mapping table
tax <- select(counts, domain, phylum, class, order, family, genus, asv_id)
rownames(tax) <- counts$asv_id

# In order to make sure that the samples match between the abundance and phenodata tables
counts_trimmed <- counts[, rownames(samples)]

#Phylogenetic tree
pseq <- readRDS("asv_10prev_phyloseq_for_lahti.RDS")
tree <- phy_tree(pseq)
tree <- makeNodeLabel(tree, method="number", prefix='n')

#Taxonomy Table
tax <- tax  %>% dplyr::rename(c(
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
counts <- as.matrix(counts)

#sequence information
refSeq <- DNAStringSet(ASV_Sequence, start=NA, end=NA, width=NA, use.names=TRUE)

tse <- TreeSummarizedExperiment(assays = list(counts = counts_trimmed),
                                colData = samples,
                                rowData = tax,
                                referenceSeq = refSeq)

#Subsetting
tse <- tse[tree$tip.label, ]

rowTree(tse) <- tree

# Collapse the tree
tree <- ape::keep.tip(phy = rowTree(tse), tip = rowLinks(tse)$nodeNum)

rowdata <- rowData(tse)
coldata <- colData(tse)
counts <- assay(tse, "counts")
rowtree <- rowTree(tse)
refseq <- referenceSeq(tse)

data <- "/microbiomeDataSets/3.14/"

saveRDS(counts, file = paste0(data,"counts.rds"))
saveRDS(coldata, file = paste0(data,"coldata.rds"))
saveRDS(rowdata , file = paste0(data, "rowdata.rds"))

ape::write.tree(rowtree, file= paste0(data,"rowtree.tre"))
gz_tree<- gzfile("rowtree.tre.gz", "w")
writeLines(readLines("rowtree.tre"), con=gz_tree)
close(gz_tree)
unlink("rowtree.tre")

Biostrings::writeXStringSet(refseq, "refseq.fasta.gz", compress = TRUE)
