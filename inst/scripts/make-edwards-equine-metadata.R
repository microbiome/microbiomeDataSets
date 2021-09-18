
# Base data for all data sets --------------------------------------------------
library(dplyr)
library(R.utils)
library(stringr)

BiocVersion <- "3.14"
path <- paste0("microbiomeDataSets/",BiocVersion,"/")

df_Base <- DataFrame(
    BiocVersion = "3.14",
    SourceVersion = NA,
    Coordinate_1_based = TRUE,
    Species = "Equus",
    TaxonomyId = "9789",
    SourceVersion = Sys.time(),
    Genome = NA,
    SourceUrl = "https://doi.org/10.1186/s42523-020-0023-1",
    DataProvider = "https://github.com/mibwurrepo/EdwardsJ_2019_EquineCoreMicrobiome",
    Maintainer = "Sudarshan Shetty <sudarshanshetty9@gmail.com>"
)

df <- rbind(
    cbind(df_Base,
          data.frame(
              Title = "Edwards equine gut bacteria counts",
              Description = paste0("Bacterial count matrix for the Edwards equine gut dataset"),
              SourceType = "RDS",
              RDataClass = "matrix",
              DispatchClass = "RDS",
              RDataPath = paste0(path,"edwards-equine/bacteria_counts.rds"),
              Tags = NA)),
    cbind(df_Base,
          data.frame(
              Title = "Edwards equine gut row data set",
              Description = paste0("Bacterial taxonomy table for the Edwards equine gut dataset"),
              SourceType = "RDS",
              RDataClass = "DFrame",
              DispatchClass = "RDS",
              RDataPath = paste0(path,"edwards-equine/bacteria_rowdata.rds"),
              Tags = NA)),
    cbind(df_Base,
          data.frame(
              Title = "Edwards equine gut phylogenetic tree data set",
              Description = paste0("Bacterial phylogenetic tree for the Edwards equine gut dataset"),
              SourceType = "TXT",
              RDataClass = "character",
              DispatchClass = "FilePath",
              RDataPath = paste0(path,"edwards-equine/bacteria_rowtree.tre.gz"),
              Tags = NA)),
    cbind(df_Base,
          data.frame(
              Title = "Edwards equine gut sample data set",
              Description = paste0("Sample information for the Edwards equine gut dataset"),
              SourceType = "CSV",
              RDataClass = "Dframe",
              DispatchClass = "RDS",
              RDataPath = paste0(path,"edwards-equine/bacteria_coldata.rds"),
              Tags = NA)),
    cbind(df_Base,
          data.frame(
              Title = "Edwards equine gut bacteria counts",
              Description = paste0("Fungal count matrix for the Edwards equine gut dataset"),
              SourceType = "RDS",
              RDataClass = "matrix",
              DispatchClass = "RDS",
              RDataPath = paste0(path,"edwards-equine/fungi_counts.rds"),
              Tags = NA)),
    cbind(df_Base,
          data.frame(
              Title = "Edwards equine gut row data set",
              Description = paste0("Fungal taxonomy table for the Edwards equine gut dataset"),
              SourceType = "RDS",
              RDataClass = "DFrame",
              DispatchClass = "RDS",
              RDataPath = paste0(path,"edwards-equine/fungi_rowdata.rds"),
              Tags = NA)),
    cbind(df_Base,
          data.frame(
              Title = "Edwards equine gut phylogenetic tree data set",
              Description = paste0("Fungal phylogenetic tree for the Edwards equine gut dataset"),
              SourceType = "TXT",
              RDataClass = "character",
              DispatchClass = "FilePath",
              RDataPath = paste0(path,"edwards-equine/fungi_rowtree.tre.gz"),
              Tags = NA)),
    cbind(df_Base,
          data.frame(
              Title = "Edwards equine gut sample data set",
              Description = paste0("Sample information for the Edwards equine gut dataset"),
              SourceType = "CSV",
              RDataClass = "Dframe",
              DispatchClass = "RDS",
              RDataPath = paste0(path,"edwards-equine/fungi_coldata.rds"),
              Tags = NA))
)


df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="") 


write.csv(df, file = paste0("../extdata/",BiocVersion,"/metadata-edwards-equine.csv"),
          row.names = FALSE)