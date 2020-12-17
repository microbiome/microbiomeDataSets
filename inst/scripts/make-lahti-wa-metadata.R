
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

df_Base <- DataFrame(
    BiocVersion = "3.13",
    SourceVersion = NA,
    Coordinate_1_based = TRUE
)

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Lahti Western Adult microbiome counts",
                    Description = paste0("Count matrix for the Lahti Western Adult microbiome dataset"),
                    SourceType = "CSV",
                    SourceUrl = "https://doi.org/10.1038/ncomms5344",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <microbiome-admin@googlegroups.com>",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-wa/counts.rds",
                    Tags = NA),
          DataFrame(Title = "Lahti Western Adult microbiome row data",
                    Description = paste0("Row data for the Lahti Western Adult microbiome dataset"),
                    SourceType = "CSV",
                    SourceUrl = "https://doi.org/10.1038/ncomms5344",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <microbiome-admin@googlegroups.com>",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-wa/rowdata.rds",
                    Tags = NA),
          DataFrame(Title = "Lahti Western Adult sample data",
                    Description = paste0("Sample data for the Lahti Western Adult dataset"),
                    SourceType = "CSV",
                    SourceUrl = "https://doi.org/10.1038/ncomms5344",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <microbiome-admin@googlegroups.com>",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/lahti-wa/coldata.rds",
                    Tags = NA))
)

df$Species <- "Homo sapiens"
df$TaxonomyId <- "9606"
df$SourceVersion <- Sys.time()
df$Genome <- NA
df$Tags <- paste0(df$Tags,"Microbiome")

write.csv(df, file = "inst/extdata/metadata-lahti-wa.csv", row.names = FALSE)
