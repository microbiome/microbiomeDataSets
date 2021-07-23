library(readxl)
library(tidyverse)

# Read metabolite data
nmr <- read_excel("NMR_Quantification.xlsx") %>%
         rename(Rat = "RAT ID") %>%
         rename(Group = "Group ID") %>%
	 select(-Group)

# Sample metadata
meta <- read_excel("Otut_abundanssit_metadata_ok.xlsx", sheet = "Metadata") %>%
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
	 mutate(XOS = as.numeric(str_replace(replace_na(str_match(meta$Diet, "XOS"), 0), "XOS", "1")))
	 
# Read microbiota data
ngs <- read_excel("Otut_abundanssit_metadata_ok.xlsx", sheet = "OTU table siisti")
## Separate taxonomy table and abundances
tax <- ngs[, 1:7] %>%
         rename(OTU = "OTU ID") %>%
         rename(Phylum = "D1") %>%
         rename(Class = "D2") %>%
         rename(Order = "D3") %>%
         rename(Family = "D4") %>%
         rename(Genus = "D5") %>%
         rename(Species = "D6") 	 

otu <- ngs[, 8:ncol(ngs)]


# Only Cecum data was used in the paper; separate metadata and blood+other measurements
inds <- which(meta$Site == "Cecum") # Check matching.txt
vars <- c("Sample", "Rat", "Site", "Diet", "Fat", "XOS")
meta_cecum <- meta[inds, vars]

# Manually checked that the sample order corresponds;
# let us rename the samples so they have same names in the different tables
otu_cecum <- otu[, inds]
colnames(otu_cecum) <- meta_cecum$Sample

# Biomarkers
# Manually checked that the sample order corresponds;
# let us rename the samples so they have same names in the different tables
bm <- t(meta[inds, -match(vars, colnames(meta))])
colnames(bm) <- meta_cecum$Sample


# Convert to TSE object
nmr$Rat <- NULL # Can be removed after matching
nmr <- t(nmr) # features x samples

# Save the data components
saveRDS(otu_cecum, file = "microbiome_counts.rds")
saveRDS(tax, file = "microbiome_rowdata.rds")
saveRDS(meta_cecum, file = "coldata.rds")
saveRDS(nmr, file = "metabolites.rds")
saveRDS(bm, file = "biomarkers.rds")



