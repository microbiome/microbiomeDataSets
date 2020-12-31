
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
    SourceUrl = "https://dx.doi.org/10.1038/ncomms7342",
    DataProvider = NA,
    Maintainer = "Leo Lahti <leo.lahti@iki.fi>"
)

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "O'Keefe diet swap microbiome counts",
                    Description = paste0("Count matrix for the O'Keefe diet swap microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/okeefe-ds/counts.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "O'Keefe diet swap microbiome row data",
                    Description = paste0("Row data for the O'Keefe diet swap microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/okeefe-ds/rowdata.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "O'Keefe diet swap sample data",
                    Description = paste0("Sample data for the O'Keefe diet swap microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/okeefe-ds/coldata.rds",
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = "inst/extdata/metadata-okeefe-ds.csv", row.names = FALSE)
