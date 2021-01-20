
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

BiocVersion <- "3.13"
path <- paste0("microbiomeDataSets/",BiocVersion,"/")

df_Base <- DataFrame(
    BiocVersion = BiocVersion,
    SourceVersion = NA,
    Coordinate_1_based = TRUE,
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    SourceVersion = Sys.time(),
    Genome = NA,
    SourceUrl = "https://doi.org/10.1186/s40168-018-0584-3",
    DataProvider = NA,
    Maintainer = "Sudarshan Shetty <sudarshanshetty9@gmail.com>"
)

#

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Silverman Artificial Gut microbiome counts",
                    Description = paste0("Count matrix for the Silverman microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"silverman-ag/counts.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Silverman Artificial Gut row data",
                    Description = paste0("Row data for the Silverman microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"silverman-ag/rowdata.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Silverman Artificial Gut sample data",
                    Description = paste0("Sample data for the Silverman dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"silverman-ag/coldata.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Silverman Artificial Gut tree data",
                    Description = paste0("Feature tree data for the Silverman dataset"),
                    SourceType = "TXT",
                    RDataClass = NA,
                    DispatchClass = "FilePath",
                    RDataPath = paste0(path,"silverman-ag/rowtree.tre.gz"),
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = paste0("inst/extdata/",BiocVersion,"/metadata-silverman-ag.csv"),
          row.names = FALSE)