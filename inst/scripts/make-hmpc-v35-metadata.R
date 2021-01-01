
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
    SourceUrl = "https://doi.org/10.1038/nature11234,https://qiita.ucsd.edu/study/description/1927",
    DataProvider = "Qiita",
    Maintainer = "Felix GM Ernst <felix.gm.ernst@outlook.com>"
)

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Human Microbiome Project Consortium V3-5 dataset OTU count matrix",
                    Description = paste0("Count matrix for the 16S rRNA variable region 3-5 datasets of the Human Microbiome Project Consortium"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/hmpc-v35/counts.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Human Microbiome Project Consortium V3-5 dataset row data",
                    Description = paste0("Row data for the 16S rRNA variable region 3-5 datasets of the Human Microbiome Project Consortium"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/hmpc-v35/rowdata.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Human Microbiome Project Consortium V3-5 dataset sample data",
                    Description = paste0("Sample data for the 16S rRNA variable region 3-5 datasets of the Human Microbiome Project Consortium"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = "microbiomeDataSets/hmpc-v35/coldata.rds",
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Human Microbiome Project Consortium V3-5 dataset feature tree",
                    Description = paste0("Feature tree for the 16S rRNA variable region 3-5 datasets of the Human Microbiome Project Consortium"),
                    SourceType = "TXT",
                    RDataClass = NA,
                    DispatchClass = "Zip",
                    RDataPath = "microbiomeDataSets/hmpc-v35/rowtree.tre.gz",
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = "inst/extdata/metadata-hmpc-v35.csv", row.names = FALSE)
