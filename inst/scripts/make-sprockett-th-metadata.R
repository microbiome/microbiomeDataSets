
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

# BiocVersion <- "3.13"
path <- paste0("microbiomeDataSets/",BiocVersion,"/")

df_Base <- DataFrame(
  BiocVersion = BiocVersion,
  SourceVersion = NA,
  Coordinate_1_based = TRUE,
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  SourceVersion = Sys.time(),
  Genome = NA,
  SourceUrl = "https://www.nature.com/articles/s41467-020-17541-6#Abs1",
  DataProvider = NA,
  Maintainer = "Sudarshan Shetty <sudarshanshetty9@gmail.com>"
)

#

df <- rbind(
  cbind(df_Base,
        DataFrame(Title = "Sprockett Tsimane Horticulturalists microbiome counts",
                  Description = paste0("Count matrix for the Sprockett microbiome dataset"),
                  SourceType = "CSV",
                  RDataClass = "matrix",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"sprockett-th/counts.rds"),
                  Tags = NA)),
  cbind(df_Base,
        DataFrame(Title = "Sprockett Tsimane Horticulturalists row data",
                  Description = paste0("Row data for the Sprockett microbiome dataset"),
                  SourceType = "CSV",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"sprockett-th/rowdata.rds"),
                  Tags = NA)),
  cbind(df_Base,
        DataFrame(Title = "Sprockett Tsimane Horticulturalists sample data",
                  Description = paste0("Sample data for the Sprockett dataset"),
                  SourceType = "CSV",
                  RDataClass = "DFrame",
                  DispatchClass = "Rds",
                  RDataPath = paste0(path,"sprockett-th/coldata.rds"),
                  Tags = NA)),
  cbind(df_Base,
        DataFrame(Title = "Sprockett Tsimane Horticulturalists tree data",
                  Description = paste0("Feature tree data for the Sprockett dataset"),
                  SourceType = "TXT",
                  RDataClass = "character",
                  DispatchClass = "FilePath",
                  RDataPath = paste0(path,"sprockett-th/rowtree.tre.gz"),
                  Tags = NA)),
  cbind(df_Base,
        DataFrame(Title = "Sprockett Tsimane Horticulturalists reference sequence data",
                  Description = paste0("Reference sequences data for the Sprockett dataset"),
                  SourceType = "TXT",
                  RDataClass = "character",
                  DispatchClass = "FilePath",
                  RDataPath = paste0(path,"sprockett-th/refseq.fasta.gz"),
                  Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = paste0("../extdata/",BiocVersion,"/metadata-sprockett-th.csv"),
          row.names = FALSE)
