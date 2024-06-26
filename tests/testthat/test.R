context("Loading of GrieneisenTSData data")
test_that("Loading of GrieneisenTSData data", {
  tse <- GrieneisenTSData()
  expect_s4_class(tse, "TreeSummarizedExperiment")
})

context("Loading of LahtiMData data")
test_that("Loading of LahtiMData data", {
  mae <- LahtiMLData()
  expect_s4_class(mae, "MultiAssayExperiment")
  tse <- LahtiMData()
  expect_s4_class(tse, "TreeSummarizedExperiment")
})

context("Loading of OKeefeDSData data")
test_that("Loading of OKeefeDSData data", {
  tse <- OKeefeDSData()
  expect_s4_class(tse, "TreeSummarizedExperiment")
})

context("Loading of SilvermanAGutData data")
test_that("Loading of SilvermanAGutData data", {
  tse <- SilvermanAGutData()
  expect_s4_class(tse, "TreeSummarizedExperiment")
})

context("Loading of SongQAData data")
test_that("Loading of SongQAData data", {
  tse <- SongQAData()
  expect_s4_class(tse, "TreeSummarizedExperiment")
})

context("Loading of GrieneisenTS data")
test_that("Loading of a GrieneisenTS data", {
  tse <- GrieneisenTSData()
  expect_s4_class(tse, "TreeSummarizedExperiment")
})

context("Loading of SprockettTHData data")
test_that("Loading of SprockettTHData data", {
  tse <- SprockettTHData()
  expect_s4_class(tse, "TreeSummarizedExperiment")
})
