# Base data for all data sets --------------------------------------------------
library(S4Vectors)

BiocVersion <- "3.13"
path <- paste0("microbiomeDataSets/",BiocVersion,"/")

df_Base <- DataFrame(
    BiocVersion = BiocVersion,
    SourceVersion = NA,
    Coordinate_1_based = TRUE,
    Species = "Mus musculus",
    TaxonomyId = "10090",
    SourceVersion = Sys.time(),
    Genome = NA,
    SourceUrl = "https://doi.org/10.3390/ijerph18084049",
    DataProvider = NA,
    Maintainer = "Leo Lahti <leo.lahti@iki.fi>"
)


df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Hintikka XO microbiome counts",
                    Description = paste0("Count matrix for the Hintikka microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/microbiome_counts.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Hintikka XO microbiome row data",
                    Description = paste0("Row data for the Hintikka microbiome dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/microbiome_rowdata.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Hintikka XO metabolites counts",
                    Description = paste0("Count matrix for the Hintikka metabolites dataset"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/metabolite_counts.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Hintikka XO sample data",
                    Description = paste0("Sample data for the Hintikka dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/coldata.rds"),
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = paste0("inst/extdata/",BiocVersion,"/metadata-hintikka-xo.csv"),
          row.names = FALSE)
