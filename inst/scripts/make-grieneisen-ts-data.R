library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ape)
library(forcats)

#Sample phenotype data
samples <- read.csv("git_metadata_and_community_phenotypes.csv")
rownames(samples) <- samples$sample
samples$X <- NULL
samples[,c("baboon_id", "sex", "social_group", "season" ,"plate")] <- lapply(samples[,c("baboon_id", "sex", "social_group", "season" ,"plate")],as.factor)
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

tax$ASV<- as.factor(tax$ASV)

#Arranging counts
rownames(counts) <- counts$asv_id
counts[, c("asv_id", "asv_sequence", "domain", "phylum", "class", "order", "family", "genus")]  <- NULL

#sequence information
refSeq <- DNAStringSet(ASV_Sequence, start=NA, end=NA, width=NA, use.names=TRUE)

tse <- TreeSummarizedExperiment(assays = list(counts = counts_trimmed),
                                colData = samples,
                                rowData = tax,
                                referenceSeq = refSeq)

#Subsetting
subsetByLeaf <- function(tse, rowLeaf) {
  # if rowLeaf is provided as node labels, convert them to node numbers
  if (is.character(rowLeaf)) {
    rowLeaf <- convertNode(tree = rowTree(tse), node = rowLeaf)
  }
  # subset data by leaves
  sse <- subsetByNode(tse, rowNode = rowLeaf)
  # update the row tree
  ## -------------- new tree: drop leaves ----------
  oldTree <- rowTree(sse)
  newTree <- ape::keep.tip(phy = oldTree, tip = rowLeaf)
  ## -------------- update the row link ----------
  # track the tree
  track <- trackNode(oldTree)
  track <- ape::keep.tip(phy = track, tip = rowLeaf)
  # row links
  rowL <- rowLinks(sse)
  rowL <- DataFrame(rowL)
  # update the row links:
  #   1. use the alias label to track and updates the nodeNum
  #   2. the nodeLab should be updated based on the new tree using the new
  #      nodeNum
  #   3. lastly, update the nodeLab_alias
  rowL$nodeNum <- convertNode(tree = track, node = rowL$nodeLab_alias,
                              message = FALSE)
  rowL$nodeLab <- convertNode(tree = newTree, node = rowL$nodeNum,
                              use.alias = FALSE, message = FALSE)
  rowL$nodeLab_alias <- convertNode(tree = newTree, node = rowL$nodeNum,
                                    use.alias = TRUE, message = FALSE)
  rowL$isLeaf <- isLeaf(tree = newTree, node = rowL$nodeNum)
  rowNL <- new("LinkDataFrame", rowL)
  ## update the row tree and links
  BiocGenerics:::replaceSlots(sse,
                              rowLinks = rowNL,
                              rowTree = list(phylo = newTree))
}

# Choose leaves
leaves <- tree$tip.label[1:613]
subsetByLeaf(tse = tse, rowLeaf = leaves )

rowdata <- rowData(tse)
coldata <- colData(tse)
counts <- assay(tse, "counts")
rowtree <- rowTree(tse)
refseq <- referenceSeq(tse)

counts <- as.matrix(counts)

data <- "/microbiomeDataSets/3.14/"

saveRDS(counts, file = paste0(data,"counts.rds"))
saveRDS(coldata, file = paste0(data,"coldata.rds"))
saveRDS(rowdata , file = paste0(data, "rowdata.rds"))

write.tree(rowtree, file= paste0(data,"rowtree.tre"))
gz_tree<- gzfile("rowtree.tre.gz", "w")
writeLines(readLines("rowtree.tre"), con=gz_tree)
close(gz_tree)
unlink("rowtree.tre")

Biostrings::writeXStringSet(refseq, "refseq.fasta.gz", compress = TRUE)
