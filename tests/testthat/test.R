
context("Loading of example data")
test_that("example data", {
    tse <- LahtiMData()
    expect_s4_class(tse, "MultiAssayExperiment")
    mae <- LahtiMLData()
    expect_s4_class(mae, "TreeSummarizedExperiment")
})
