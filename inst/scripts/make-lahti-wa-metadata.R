
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
    SourceUrl = "https://doi.org/10.1038/ncomms5344",
    DataProvider = NA,
    Maintainer = "Leo Lahti <leo.lahti@iki.fi>"
)

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Lahti Western Adult microbiome counts",
                    Description = paste0("Count matrix for the Lahti Western Adult microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-wa/counts.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Lahti Western Adult microbiome row data",
                    Description = paste0("Row data for the Lahti Western Adult microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-wa/rowdata.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Lahti Western Adult sample data",
                    Description = paste0("Sample data for the Lahti Western Adult dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-wa/coldata.rds",
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = "inst/extdata/3.13/metadata-lahti-wa.csv", row.names = FALSE)
