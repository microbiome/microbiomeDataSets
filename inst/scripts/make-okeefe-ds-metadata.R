
# Base data for all data sets --------------------------------------------------
library(S4Vectors)

df_Base <- DataFrame(
    BiocVersion = "3.13",
    SourceVersion = NA,
    Coordinate_1_based = TRUE
)

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "O'Keefe diet swap microbiome counts",
                    Description = paste0("Count matrix for the O'Keefe diet swap microbiome dataset"),
                    SourceType = "CSV",
                    SourceUrl = "https://dx.doi.org/10.1038/ncomms7342",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <leo.lahti@iki.fi>",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/okeefe-ds/counts.rds",
                    Tags = NA),
          DataFrame(Title = "O'Keefe diet swap microbiome row data",
                    Description = paste0("Row data for the O'Keefe diet swap microbiome dataset"),
                    SourceType = "CSV",
                    SourceUrl = "https://dx.doi.org/10.1038/ncomms7342",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <leo.lahti@iki.fi>",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/okeefe-ds/rowdata.rds",
                    Tags = NA),
          DataFrame(Title = "O'Keefe diet swap sample data",
                    Description = paste0("Sample data for the O'Keefe diet swap microbiome dataset"),
                    SourceType = "CSV",
                    SourceUrl = "https://dx.doi.org/10.1038/ncomms7342",
                    DataProvider = NA,
                    Maintainer = "Leo Lahti <leo.lahti@iki.fi>",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/okeefe-ds/coldata.rds",
                    Tags = NA))
)

df$Species <- "Homo sapiens"
df$TaxonomyId <- "9606"
df$SourceVersion <- Sys.time()
df$Genome <- NA
df$Tags <- paste0(df$Tags,"Microbiome")

write.csv(df, file = "inst/extdata/metadata-okeefe-ds.csv", row.names = FALSE)
