
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
                  Description = paste0("Row data for the bacteria microbiome dataset"),
                  SourceType = "Rds",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"tax.rds"),
                  Tags = "Microbiome")),
  cbind(df_Base,
        data.frame(
                  Title = "Phylogenetic tree",
                  Description = paste0("Bacteria microbiome phylogenetic tree dataset"),
                  SourceType = "Rds",
                  RDataClass = "Phylo",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"tree.rds"),
                  Tags = "Microbiome")),
  cbind(df_Base,
        data.frame(
                  Title = "Sample data",
                  Description = paste0("Bacteria microbiome sample dataset"),
                  SourceType = "Csv",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"samples.rds"),
                  Tags = "Microbiome")),
  cbind(df_Base,
        data.frame(
                  Title = "Sequence data",
                  Description = paste0("Bacteria microbiome sequence dataset"),
                  SourceType = "Rds",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"sequence.rds"),
                  Tags = "Microbiome"))
)

write.csv(df, file = paste0("../extdata/",BiocVersion,"/metadata-grieneisen-ts.csv"),
          row.names = FALSE)
