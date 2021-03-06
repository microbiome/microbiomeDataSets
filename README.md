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
accessed are available from Bioconductor 
[here](https://bioconductor.org/packages/release/bioc/vignettes/ExperimentHub/inst/doc/CreateAnExperimentHubPackage.html) 
and 
[here](https://bioconductor.org/packages/release/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html)

# Code of conduct

Please note that the microbiomeDataSets project is released with a 
[Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
