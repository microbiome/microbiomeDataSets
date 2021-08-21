
library(dplyr)
library(R.utils)
library(stringr)

BiocVersion <- "3.14"
path <- paste0("microbiomeDataSets/",BiocVersion,"/")

df_Base <- data.frame(
  BiocVersion = BiocVersion,
  SourceVersion = NA,
  Coordinate_1_based = TRUE,
  Species = "Papio cynocephalus",
  TaxonomyId = "9556",
  SourceVersion=Sys.time(),
  Genome = NA,
  SourceUrl = "https://science.sciencemag.org/content/373/6551/181",
  DataProvider = "University of Minnesota"
)

df <- rbind(
  cbind(df_Base,
       data.frame(
                  Title = "Grieneisen Baboon counts data set",
                  Description = paste0("Count matrix for the Grieneisen Baboon dataset"),
                  SourceType = "Rds",
                  RDataClass = "matrix",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"baboon-ts/counts.rds"),
                  Tags = "NA")),
  cbind(df_Base,
        data.frame(
                  Title = "Grieneisen Baboon row data set",
                  Description = paste0("Taxonomy table for the Grieneisen Baboon dataset"),
                  SourceType = "Rds",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"baboon-ts/rowdata.rds"),
                  Tags = "NA")),
  cbind(df_Base,
        data.frame(
                  Title = "Grieneisen Baboon phylogenetic tree data set",
                  Description = paste0("Phylogenetic tree for the Grieneisen Baboon dataset"),
                  SourceType = "TXT",
                  RDataClass = "character",
                  DispatchClass = "FilePath",
                  RDataPath = paste0(path,"baboon-ts/rowtree.tre.gz"),
                  Tags = "NA")),
  cbind(df_Base,
        data.frame(
                  Title = "Grieneisen Baboon sample data set",
                  Description = paste0("Sample information for the Grieneisen Baboon dataset"),
                  SourceType = "CSV",
                  RDataClass = "Dframe",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"baboon-ts/coldata.rds"),
                  Tags = "NA")),
  cbind(df_Base,
        data.frame(
                  Title = "Grieneisen Baboon sequence data set",
                  Description = paste0("Sequence data for the Grieneisen Baboon dataset"),
                  SourceType = "TXT",
                  RDataClass = "character",
                  DispatchClass = "FilePath",
                  RDataPath = paste0(path,"baboon-ts/refseq.fasta.gz"),
                  Tags = "NA"))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="") 

write.csv(df, file = paste0("../extdata/",BiocVersion,"/metadata-grieneisen-ts.csv"),
          row.names = FALSE)
