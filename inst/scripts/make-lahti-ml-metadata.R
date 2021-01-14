
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

df_Base <- DataFrame(
    BiocVersion = "3.13",
    SourceVersion = NA,
    Coordinate_1_based = TRUE,
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    SourceVersion = Sys.time(),
    Genome = NA,
    SourceUrl = "http://dx.doi.org/10.7717/peerj.32",
    DataProvider = NA,
    Maintainer = "Leo Lahti <leo.lahti@iki.fi>"
)

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Lahti ML microbiome counts",
                    Description = paste0("Count matrix for the Lahti microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-ml/microbiome_counts.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Lahti ML microbiome row data",
                    Description = paste0("Row data for the Lahti microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-ml/microbiome_rowdata.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Lahti ML lipids counts",
                    Description = paste0("Count matrix for the Lahti lipids dataset"),
                    SourceType = "CSV",
                   RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-ml/lipids_counts.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Lahti ML sample data",
                    Description = paste0("Sample data for the Lahti dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-ml/coldata.rds",
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = "inst/extdata/3.13/metadata-lahti-ml.csv", row.names = FALSE)
