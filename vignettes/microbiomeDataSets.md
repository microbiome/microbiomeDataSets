---
title: "microbiomeDataSets"
output: 
  BiocStyle::html_document:
    fig_height: 7
    fig_width: 10
    toc: yes
    toc_depth: 2
    number_sections: true
vignette: >
  %\VignetteIndexEntry{microbiomeDataSets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(microbiomeDataSets)
```


# Microbiome example data sets

The data sets are primarily named by the first author of the
associated publication, together with a descriptive suffix. Aliases
are provided for some of the data sets.


### Intestinal microbiota profiling of 1006 Western adults

This data set from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) comes with 130 genus-like taxonomic groups across 1006 western adults with no reported health complications. Some subjects have also short time series.

This data set is based on the Human Intestinal Tract (HIT)Chip
phylogenetic 16S microarray [(Rajilić‐Stojanović _et al._
2009)](https://doi.org/10.1111/j.1462-2920.2009.01900.x). This
profiling technology differs from the more widely used 16S rRNA
amplicon sequencing.

The HITChip Atlas data set is available via ExperimentHub in
TreeSummarizedExperiment format. This data is also available in
phyloseq format through the microbiome R package, and via [Data
Dryad](http://doi.org/10.5061/dryad.pk75d) in tabular format.

The data can be loaded in R in TreeSummarizedExperiment format as
follows.


```r
# Data citation doi: 10.1038/ncomms5344
microbiomeDataSets::LahtiWAData()
#> class: TreeSummarizedExperiment 
#> dim: 130 1151 
#> metadata(0):
#> assays(1): counts
#> rownames(130): Actinomycetaceae Aerococcus ... Xanthomonadaceae
#>   Yersinia et rel.
#> rowData names(3): Phylum Family Genus
#> colnames(1151): Sample-1 Sample-2 ... Sample-1171 Sample-1172
#> colData names(10): age sex ... time sample
#> reducedDimNames(0):
#> altExpNames(0):
#> rowLinks: NULL
#> rowTree: NULL
#> colLinks: NULL
#> colTree: NULL

# Alias
# microbiomeDataSets::atlas1006()
```

### Diet swap between Rural and Western populations

A two-week diet swap study between western (USA) and traditional
(rural Africa) diets, reported in [O'Keefe et al. Nat. Comm. 6:6342,
2015](http://dx.doi.org/10.1038/ncomms7342). The data is also
available for download from [Data
Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). 

This data set is based on the Human Intestinal Tract (HIT)Chip
phylogenetic 16S microarray [(Rajilić‐Stojanović _et al._
2009)](https://doi.org/10.1111/j.1462-2920.2009.01900.x). This
profiling technology differs from the more widely used 16S rRNA
amplicon sequencing.


```r
# Data citation doi: 10.1038/ncomms7342
microbiomeDataSets::OKeefeDSData()
#> class: TreeSummarizedExperiment 
#> dim: 130 222 
#> metadata(0):
#> assays(1): counts
#> rownames(130): Actinomycetaceae Aerococcus ... Xanthomonadaceae
#>   Yersinia et rel.
#> rowData names(3): Phylum Family Genus
#> colnames(222): Sample-1 Sample-2 ... Sample-221 Sample-222
#> colData names(8): subject sex ... timepoint.within.group bmi_group
#> reducedDimNames(0):
#> altExpNames(0):
#> rowLinks: NULL
#> rowTree: NULL
#> colLinks: NULL
#> colTree: NULL

# Alias
# microbiomeDataSets::dietswap()
```

### Intestinal microbiota versus blood metabolites

Data set from [Lahti et al. PeerJ 1:e32,
2013](https://doi.org/10.7717/peerj.32) characterizes associations
between human intestinal microbiota and blood serum lipids. Note that
this data set contains an additional assay of lipid species, and is
therefore provided as MultiAssayExperiment object.

This data set is based on the Human Intestinal Tract (HIT)Chip
phylogenetic 16S microarray [(Rajilić‐Stojanović _et al._
2009)](https://doi.org/10.1111/j.1462-2920.2009.01900.x). This
profiling technology differs from the more widely used 16S rRNA
amplicon sequencing.

Load the data in R with:


```r
# Data citation doi: 10.7717/peerj.32
microbiomeDataSets::LahtiMLData()
#> snapshotDate(): 2020-04-27
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] microbiome: TreeSummarizedExperiment with 130 rows and 44 columns
#>  [2] lipids: SummarizedExperiment with 389 rows and 44 columns
#> Features:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DFrame
#>  sampleMap() - the sample availability DFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices

# Alias
# microbiomeDataSets::peerj32()
```
