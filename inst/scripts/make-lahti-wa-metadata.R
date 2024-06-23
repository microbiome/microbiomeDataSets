
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
    SourceUrl = "https://doi.org/10.1038/ncomms5344",
    DataProvider = NA,
    Maintainer = "Leo Lahti <leo.lahti@iki.fi>"
)

#

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Lahti Western Adult microbiome counts",
                    Description = paste0("Count matrix for the Lahti Western Adult microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"lahti-wa/counts.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Lahti Western Adult microbiome row data",
                    Description = paste0("Row data for the Lahti Western Adult microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"lahti-wa/rowdata.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Lahti Western Adult sample data",
                    Description = paste0("Sample data for the Lahti Western Adult dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"lahti-wa/coldata.rds"),
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = paste0("../extdata/",BiocVersion,"/metadata-lahti-wa.csv"),
          row.names = FALSE)
