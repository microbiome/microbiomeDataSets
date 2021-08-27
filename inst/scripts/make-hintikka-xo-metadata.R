# Base data for all data sets --------------------------------------------------
library(S4Vectors)

BiocVersion <- "3.14"
path <- paste0("microbiomeDataSets/",BiocVersion,"/")

df_Base <- DataFrame(
    BiocVersion = BiocVersion,
    SourceVersion = NA,
    Coordinate_1_based = TRUE,
    Species = "Rattus norvegicus",
    TaxonomyId = "10116",
    SourceVersion = Sys.time(),
    Genome = NA,
    SourceUrl = "https://ndownloader.figshare.com/files/",
    DataProvider = "University of Jyvaskyla",
    Maintainer = "Leo Lahti <leo.lahti@iki.fi>"
)


df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Hintikka XO sample data",
                    Description = paste0("Sample data for the HintikkaXO dataset"),
                    SourceType = "XLS/XLSX",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/coldata.rds"),
                    Tags = NA)),
		    
    cbind(df_Base,
          DataFrame(Title = "Hintikka XO microbiota row data",
                    Description = paste0("Row data for the HintikkaXO microbiota dataset"),
                    SourceType = "XLS/XLSX",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/microbiota_rowdata.rds"),
                    Tags = NA)),

    cbind(df_Base,
          DataFrame(Title = "Hintikka XO microbiota counts",
                    Description = paste0("Count matrix for the HintikkaXO microbiota dataset"),
                    SourceType = "XLS/XLSX",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/microbiota_counts.rds"),
                    Tags = NA)),		    

    cbind(df_Base,
          DataFrame(Title = "Hintikka XO metabolites",
                    Description = paste0("Count matrix for the HintikkaXO metabolites"),
                    SourceType = "XLS/XLSX",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/metabolites_nmr.rds"),
                    Tags = NA)),

    cbind(df_Base,
          DataFrame(Title = "Hintikka XO biomarkers",
                    Description = paste0("Data matrix for the HintikkaXO biomarkers"),
                    SourceType = "XLS/XLSX",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"hintikka-xo/biomarkers_signals.rds"),
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = paste0("../extdata/",BiocVersion,"/metadata-hintikka-xo.csv"),
          row.names = FALSE)

