
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
    SourceUrl = "http://msystems.asm.org/content/1/3/e00021-16,https://qiita.ucsd.edu/study/description/10394",
    DataProvider = "Qiita",
    Maintainer = "Felix GM Ernst <felix.gm.ernst@outlook.com>"
)

#

df <- rbind(
    cbind(df_Base,
          DataFrame(Title = "Song et al. 2012 dataset OTU count matrix",
                    Description = paste0("Count matrix for the Song et al. 2012 dataset"),
                    SourceType = "CSV",
                    RDataClass = "matrix",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"song-qa/counts.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Song et al. 2012 dataset row data",
                    Description = paste0("Row data for the Song et al. 2012 dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"song-qa/rowdata.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Song et al. 2012 dataset sample data",
                    Description = paste0("Sample data for the Song et al. 2012 dataset"),
                    SourceType = "CSV",
                    RDataClass = "DFrame",
                    DispatchClass = "Rds",
                    RDataPath = paste0(path,"song-qa/coldata.rds"),
                    Tags = NA)),
    cbind(df_Base,
          DataFrame(Title = "Song et al. 2012 dataset feature tree",
                    Description = paste0("Feature tree for the Song et al. 2012 dataset"),
                    SourceType = "TXT",
                    RDataClass = NA,
                    DispatchClass = "FilePath",
                    RDataPath = paste0(path,"song-qa/rowtree.tre.gz"),
                    Tags = NA))
)

df$Tags <- paste(df$Tags[!is.na(df$Tags)],"Microbiome",collapse = ":",sep="")

write.csv(df, file = paste0("../extdata/",BiocVersion,"/metadata-song-qa.csv"),
          row.names = FALSE)
