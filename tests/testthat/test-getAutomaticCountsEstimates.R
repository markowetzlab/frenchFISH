context("Testing automatic spot counting correction")

cell_radius = 8
section_height = 4
automatic_counts = cbind(area = c(250, 300, 450),
                         red = c(0, 2, 4),
                         green = c(5, 3, 1),
                         blue = c(3, 0, 2))
automatic_counts_estimates <- getAutomaticCountsEstimates(automatic_counts, cell_radius, section_height)
# unknown failure: getAutomaticCountsEstimates(cbind(area=c(5,4,10), red=c(1,2,3), green=c(3,4,5), blue=c(3,0,1)), 8, 4)

test_that("getAutomaticCountsEstimates function works", {
  expect_true(is.data.frame(automatic_counts_estimates))
  expect_equal(dim(automatic_counts_estimates), c(3, 4))
  expect_equal(toString(automatic_counts_estimates$Probe), c("red, green, blue"))
  expect_equal(round(automatic_counts_estimates$lowCI, 3), c(2.810, 1.753, 0.815))
  expect_equal(round(automatic_counts_estimates$median, 3), c(5.895, 4.211, 2.526))
  expect_equal(round(automatic_counts_estimates$highCI, 3), c(12.365, 10.116, 7.833))

  #expect_equal(round(automatic_counts_estimates[1,1], 3), 2.810)
  #expect_equal(round(automatic_counts_estimates[1,2], 3), 5.895)
  #expect_equal(round(automatic_counts_estimates[1,3], 3), 12.365)
  #expect_equal(toString(automatic_counts_estimates[1,4]), "red")
  #expect_equal(round(automatic_counts_estimates[2,1], 3), 1.753)
  #expect_equal(round(automatic_counts_estimates[2,2], 3), 4.211)
  #expect_equal(round(automatic_counts_estimates[2,3], 3), 10.116)
  #expect_equal(toString(automatic_counts_estimates[2,4]), "green")
  #expect_equal(round(automatic_counts_estimates[3,1], 3), 0.815)
  #expect_equal(round(automatic_counts_estimates[3,2], 3), 2.526)
  #expect_equal(round(automatic_counts_estimates[3,3], 3), 7.833)
  #expect_equal(toString(automatic_counts_estimates[3,4]), "blue")
})

test_that("getAutomaticCountsEstimates function throws error messages when input arguments are invalid", {
  expect_error(getAutomaticCountsEstimates(automatic_counts, "radius", 4), "radius must be numeric")
  expect_error(getAutomaticCountsEstimates(automatic_counts, NA, 4), "radius must be numeric")
  expect_error(getAutomaticCountsEstimates(automatic_counts, 8, "height"), "height must be numeric")
  expect_error(getAutomaticCountsEstimates(automatic_counts, 8, NA), "height must be numeric")
  expect_error(getAutomaticCountsEstimates(automatic_counts, -8, 4), "radius must be greater than 0")
  expect_error(getAutomaticCountsEstimates(automatic_counts, 8, -4), "height must be greater than 0")
  expect_error(getAutomaticCountsEstimates("probeCounts", 8, 4), "probeCounts must be a matrix")
  expect_error(getAutomaticCountsEstimates(5, 8, 4), "probeCounts must be a matrix")
  expect_error(getAutomaticCountsEstimates(matrix(NA, nrow = 0, ncol = 0), 8, 4), "probeCounts must have at least one row")
  expect_error(getAutomaticCountsEstimates(cbind(area=c(500)), 8, 4), "probeCounts must have at least one column for nuclear blob area and one column for spot counts at a probe")
  expect_error(getAutomaticCountsEstimates(cbind(red=c(2), blue=c(0)), 8, 4), 'First column of probeCounts must be named "area" and contain nuclear blob areas')
  expect_error(getAutomaticCountsEstimates(cbind(red=c(2), area=c(500)), 8, 4), 'First column of probeCounts must be named "area" and contain nuclear blob areas')
  expect_error(getAutomaticCountsEstimates(cbind(Area=c(2), red=c(500)), 8, 4), 'First column of probeCounts must be named "area" and contain nuclear blob areas')
  #expect_error(getAutomaticCountsEstimates(cbind(area=c(500), red=c(NA)), 8, 4), "probeCounts cannot have any NA or NaN values")
  #expect_error(getAutomaticCountsEstimates(cbind(area=c(500), red=c(NaN)), 8, 4), "probeCounts cannot have any NA or NaN values")
  expect_error(getAutomaticCountsEstimates(cbind(area=c(500), red=c(Inf)), 8, 4), "probeCounts cannot have any Inf or -Inf values")
  expect_error(getAutomaticCountsEstimates(cbind(area=c(500), red=c(-Inf)), 8, 4), "probeCounts cannot have any Inf or -Inf values")
  expect_error(getAutomaticCountsEstimates(cbind(area=c(500), red=c(1.5)), 8, 4), "All non-NA/NaN counts in probeCounts must be non-negative integers")
  expect_error(getAutomaticCountsEstimates(cbind(area=c(500), red=c(-1)), 8, 4), "All non-NA/NaN counts in probeCounts must be non-negative integers")
  expect_error(getAutomaticCountsEstimates(cbind(area=c(500), red=c(-1.5)), 8, 4), "All non-NA/NaN counts in probeCounts must be non-negative integers")
  expect_error(getAutomaticCountsEstimates(cbind(area=c(-500), red=c(1)), 8, 4), "All values in area column must be greater than 0")
  expect_error(getAutomaticCountsEstimates(cbind(area=c(0), red=c(1)), 8, 4), "All values in area column must be greater than 0")
  expect_warning(getAutomaticCountsEstimates(cbind(area=c(500), red=c(0)), 8, 4), "All non-NA/NaN counts of red probe are 0. Estimated count for this probe may not be accurate")
  expect_warning(getAutomaticCountsEstimates(cbind(area=c(500), red=c(NA)), 8, 4), "All counts of red probe are NA or NaN. Estimated count will be NA")
})


