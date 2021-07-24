# `microbiomeDataSets`

<!-- badges: start -->

<!-- badges: end -->

This R package is a collection of microbiome datasets published initially 
elsewhere. The data is available as 
[`TreeSummarizedExperiment`](https://doi.org/doi:10.18129/B9.bioc.TreeSummarizedExperiment)
or 
[`MultiAssayExperiment`](https://doi.org/doi:10.18129/B9.bioc.MultiAssayExperiment)
and a list of available dataset can be retrieved via the `availableDataSets()`
function.

The aim is to provide datasets for teaching, example workflows or comparative
efforts. If you have a dataset, which you like to see in this package, please
let us know and/or provide a PR for the datasets.

# Contribution

Feel free to contribute. Have a look at how existing datasets are organized and
prepared data accordingly. It is also good to get in touch at the earliest 
convenience to discuss any issues.

## Technical aspects

Let's use a git flow kind of approach. Development version should be done 
against the `master` branch and then merged to `master` for the next release. 
(https://guides.github.com/introduction/flow/)

Resources on how data is added to Bioconductor's ExperimentHub backend and
accessed are available from Bioconductor [ExperimentHub documentation](https://bioconductor.org/packages/release/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html) and in [Creating ExperimentHub Package](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/CreateAHubPackage.html).

Basic steps:

- Assemble a TreeSummarizedExperiment from the raw data (Install the latest version from GitHub)
- Save the individiual components as rds files in inst/extdata/hub/<dataset-name>/ (If you have questions, have a look at the other datasets or check with us through [online channels](microbiome.github.io))
- Prepare the metadata file, by creating a new metadata-<dataset-name>.R in inst/scripts and run the script to create inst/extdata/metadata-<dataset-name>.csv
- Make sure that the metadata files passes the check by running AnnotationHubData::makeAnnotationHubMetadata("../microbiomeDataSets","metadata-<dataset-name>.csv")


# Code of conduct

Please note that the microbiomeDataSets project is released with a 
[Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
