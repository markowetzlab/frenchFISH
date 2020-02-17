context("Testing conversion between FISHalyseR's outputted CSV to a frenchFISH probeCounts matrix")

test_that("convertFishalyserCsvToCountMatrix function works", {
  expect_error(convertFishalyserCsvToCountMatrix(5), "pathToFishalyserCsv must be a string")
  expect_error(convertFishalyserCsvToCountMatrix("/this/path/doesnt/exist/fishalyser_counts.csv"), "/this/path/doesnt/exist/fishalyser_counts.csv does not exist")
})

