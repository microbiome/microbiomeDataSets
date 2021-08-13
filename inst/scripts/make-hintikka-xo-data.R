library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

temp_data <- "../extras/temp_data/"

# Read metabolite data
nmr <- read_excel(paste0(temp_data, "NMR_Quantification.xlsx")) %>%
         rename(Rat = "RAT ID") %>%
         rename(Group = "Group ID") %>%
	 select(-Group)

# Sample metadata
meta <- read_excel(paste0(temp_data, "Otut_abundanssit_metadata_ok.xlsx"),
           sheet = "Metadata") %>%
         rename(Rat = "Rotta ID") %>%
         rename(Group = "Ryhmä nro") %>%
         rename(Batch = "Sekvensointi ID") %>%
	 select(-Batch) %>% # Same as "Rat"
	 select(-Group) %>% # Not needed
         rename(Site = "Suolen osa") %>%
	 rename(Sample = Rat) %>% # These are sample IDs
         mutate(Rat = str_remove(Sample, "^[P|C]")) %>%	 # These are subject IDs 
	 mutate(Fat = factor(str_trim(str_remove(str_remove(Diet, "\\+ XOS"), "-fat")),
             levels = c("Low", "High"))) %>%
   mutate(XOS = as.numeric(str_replace(replace_na(str_match(Diet, "XOS"), 0), "XOS", "1")))



# Read microbiota data
ngs <- read_excel(paste0(temp_data, "Otut_abundanssit_metadata_ok.xlsx"),
  sheet = "OTU table siisti")

## Separate taxonomy table and abundances
tax <- ngs[, 1:7] %>%
         rename(OTU = "OTU ID") %>%
         rename(Phylum = "D1") %>%
         rename(Class = "D2") %>%
         rename(Order = "D3") %>%
         rename(Family = "D4") %>%
         rename(Genus = "D5") %>%
         rename(Species = "D6") 	 
tax <- tax[, c(2:7, 1)]
tax <- DataFrame(tax)
rownames(tax) <- tax$OTU

# Pick microbiota abundances
otu <- ngs[, 8:ncol(ngs)]

# Only Cecum data was used in the paper;
# separate metadata and blood+other measurements
inds <- which(meta$Site == "Cecum") # Check matching.txt
vars <- c("Sample", "Rat", "Site", "Diet", "Fat", "XOS")
meta_cecum <- as.data.frame(meta[inds, vars])
meta_cecum <- DataFrame(meta_cecum)
rownames(meta_cecum) <- meta_cecum$Sample
meta_cecum$Rat  <- as.factor(meta_cecum$Rat)
meta_cecum$Diet <- as.factor(meta_cecum$Diet)
meta_cecum$Fat  <- as.factor(meta_cecum$Fat)
#meta_cecum$XOS <- as.factor(meta_cecum$XOS)

# Manually checked that the sample order corresponds;
# let us rename the samples so they have same names in the different tables
otu_cecum <- as.matrix(otu[, inds])
rownames(otu_cecum) <- rownames(tax)
colnames(otu_cecum) <- meta_cecum$Sample

# Biomarkers
# Manually checked that the sample order corresponds;
# let us rename the samples so they have same names in the different tables
bm <- t(meta[inds, -match(vars, colnames(meta))])
bm <- as.matrix(bm)
colnames(bm) <- meta_cecum$Sample

# NMR
nmr$Rat <- NULL # Can be removed after matching
nmr <- t(nmr) # features x samples
nmr <- as.matrix(nmr)
colnames(nmr) <- meta_cecum$Sample

# Convert to MAE object
# build SummarizedExperiment objects
sem <- SummarizedExperiment(assays = list(counts=otu_cecum),
			    rowData = tax)
sen <- SummarizedExperiment(assays = list(nmr=nmr))
seb <- SummarizedExperiment(assays = list(signals=bm))

## Create a MultiAssayExperiment instance
ExpList <- ExperimentList(list(microbiota=sem,
                               metabolites=sen,
			       biomarkers=seb))
mae <- MultiAssayExperiment(experiments = ExpList,
                            colData = meta_cecum)
     
# Fetch the data from MAE object in order to ensure compabitility
# with ExperimentHub access
tax     <- rowData(mae[["microbiota"]])
counts  <- assay(mae[["microbiota"]], "counts")
nmr     <- assay(mae[["metabolites"]], "nmr")
bm      <- assay(mae[["biomarkers"]], "signals")
coldata <- colData(mae)

# Save the data components
path <- "../extras/microbiomeDataSets/3.14/hintikka-xo/"
saveRDS(tax,       file = paste0(path, "microbiota_rowdata.rds"))
saveRDS(counts,    file = paste0(path, "microbiota_counts.rds"))
saveRDS(nmr,       file = paste0(path, "metabolites_nmr.rds"))
saveRDS(bm,        file = paste0(path, "biomarkers_signals.rds"))
saveRDS(coldata,   file = paste0(path, "coldata.rds"))


