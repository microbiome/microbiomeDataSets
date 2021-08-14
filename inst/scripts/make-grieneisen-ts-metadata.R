
library(dplyr)
library(R.utils)
library(stringr)

BiocVersion <- "3.14"
path <- paste0("MicrobiomeDataSets/",BiocVersion,"/")

df_Base <- data.frame(
  BiocVersion = BiocVersion,
  SourceVersion = NA,
  Coordinate_1_based = TRUE,
  Domain = "Bacteria",
  TaxonomyId = "2",
  SourceVersion=Sys.time(),
  Genome = NA,
  SourceUrl = "https://science.sciencemag.org/content/373/6551/181",
  DataProvider = NA
)

df <- rbind(
  cbind(df_Base,
       data.frame(
                  Title = "Counts",
                  Description = paste0("Count matrix for the bacteria microbiome dataset"),
                  SourceType = "Rds",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"counts.rds"),
                  Tags = "Microbiome" )),
  cbind(df_Base,
        data.frame(
                  Title = "Row data",
                  Description = paste0("Taxonomy table for the bacteria microbiome dataset"),
                  SourceType = "Rds",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"rowdata.rds"),
                  Tags = "Microbiome")),
  cbind(df_Base,
        data.frame(
                  Title = "Phylogenetic tree",
                  Description = paste0("Phylogenetic tree for the bacteria microbiome dataset"),
                  SourceType = "Rds",
                  RDataClass = "Phylo",
                  DispatchClass = "Tre.gz",
                  RDataPath = paste0(path,"rowtree.tre.gz"),
                  Tags = "Microbiome")),
  cbind(df_Base,
        data.frame(
                  Title = "Sample data",
                  Description = paste0("Sample Matrix for the bacteria microbiome dataset"),
                  SourceType = "Csv",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"coldata.rds"),
                  Tags = "Microbiome")),
  cbind(df_Base,
        data.frame(
                  Title = "Sequence data",
                  Description = paste0("Sequence information for the bacteria microbiome dataset"),
                  SourceType = "Rds",
                  RDataClass = "DFrame",
                  DispatchClass = "Fasta.gz",
                  RDataPath = paste0(path,"refseq.fasta.gz"),
                  Tags = "Microbiome"))
)

write.csv(df, file = paste0("../extdata/",BiocVersion,"/metadata-grieneisen-ts.csv"),
          row.names = FALSE)
