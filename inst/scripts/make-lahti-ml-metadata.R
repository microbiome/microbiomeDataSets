
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

df_Base <- DataFrame(
    BiocVersion = "3.13",
    SourceVersion = NA,
    Coordinate_1_based = TRUE
)

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Lahti ML microbiome counts",
                    Description = paste0("Count matrix for the Lahti microbiome dataset"),
                    SourceType = "CSV",
                    SourceUrl = "http://dx.doi.org/10.7717/peerj.32",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <leo.lahti@iki.fi>",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-ml/microbiome_counts.rds",
                    Tags = NA),
          DataFrame(Title = "Lahti ML microbiome row data",
                    Description = paste0("Row data for the Lahti microbiome dataset"),
                    SourceType = "CSV",
                    SourceUrl = "http://dx.doi.org/10.7717/peerj.32",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <leo.lahti@iki.fi>",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-ml/microbiome_rowdata.rds",
                    Tags = NA),
          DataFrame(Title = "Lahti ML lipids counts",
                    Description = paste0("Count matrix for the Lahti lipids dataset"),
                    SourceType = "CSV",
                    SourceUrl = "http://dx.doi.org/10.7717/peerj.32",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <leo.lahti@iki.fi>",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-ml/lipids_counts.rds",
                    Tags = NA),
          DataFrame(Title = "Lahti ML sample data",
                    Description = paste0("Sample data for the Lahti dataset"),
                    SourceType = "CSV",
                    SourceUrl = "http://dx.doi.org/10.7717/peerj.32",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <leo.lahti@iki.fi>",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-ml/coldata.rds",
                    Tags = NA))
)

df$Species <- "Homo sapiens"
df$TaxonomyId <- "9606"
df$SourceVersion <- Sys.time()
df$Genome <- NA
df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = "inst/extdata/metadata-lahti-ml.csv", row.names = FALSE)
