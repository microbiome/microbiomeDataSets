library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ape)
library(forcats)

#Sample phenotype data
orig.data <- "grieneisen-ts-original/"

#Phylogenetic tree
pseq <- readRDS(paste0(orig.data, "asv_10prev_phyloseq_for_lahti.RDS"))
tree <- phy_tree(pseq)
tree <- makeNodeLabel(tree, method="number", prefix='n')
leaves <- tree$tip.label

samples <- read.csv(paste0(orig.data, "git_metadata_and_community_phenotypes.csv"))
rownames(samples) <- samples$sample
samples$X <- NULL
samples[,c("baboon_id", "sex", "social_group", "season" ,"plate")] <- lapply(samples[,c("baboon_id", "sex", "social_group", "season" ,"plate")],as.factor)
samples$collection_date <- as.Date(samples$collection_date, format="%Y-%m-%d")

#Abundance table
counts <- readRDS(paste0(orig.data, "git_ASV_table.RDS"))
rownames(counts) <- counts$asv_id
# Pick the subset of taxa for which we have tree
counts <- counts[leaves,]

#sequence information
ASV_Sequence <- counts$asv_sequence
ref_seq <- data.frame(ASV_Sequence)
rownames(ref_seq) <- rownames(counts)

#Taxonomic mapping table
tax <- select(counts, domain, phylum, class, order, family, genus, asv_id)
rownames(tax) <- counts$asv_id

# In order to make sure that the samples match between the abundance and phenodata tables
counts_trimmed <- counts[, rownames(samples)]

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
tax$ASV<- as.factor(tax$ASV)

# sequence information
refSeq <- DNAStringSet(ASV_Sequence, start=NA, end=NA, width=NA, use.names=TRUE)

tse <- TreeSummarizedExperiment(assays = list(counts = counts_trimmed[leaves, ]),
                                colData = samples,
                                rowData = tax[leaves, ],
				rowTree = tree,
                                referenceSeq = refSeq)

#Subsetting- not needed now; already subsetted above
# subsetByLeaf(tse = tse, rowLeaf = tree$tip.label)

rowdata <- rowData(tse)
coldata <- colData(tse)
counts <- assay(tse, "counts")
rowtree <- rowTree(tse)
refseq <- referenceSeq(tse)

counts <- as.matrix(counts)

data <- "microbiomeDataSets/3.14/"

saveRDS(counts, file = paste0(data,"counts.rds"))
saveRDS(coldata, file = paste0(data,"coldata.rds"))
saveRDS(rowdata , file = paste0(data, "rowdata.rds"))

treefile <- paste0(data, "rowtree.tre")
write.tree(rowtree, file=treefile)
gz_tree<- gzfile(treefile, "w")
writeLines(readLines(treefile), con=gz_tree)
close(gz_tree)
unlink(treefile)

Biostrings::writeXStringSet(refseq, "refseq.fasta.gz", compress = TRUE)
