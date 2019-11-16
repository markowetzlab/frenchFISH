context("Testing automatic spot counting correction")

cr = 8
sh = 4
cm = data.frame(area = c(250, 300, 450), red = c(0, 2, 4), green = c(5, 3, 1), blue = c(3, 0, 2))
ccm <- getPPcountsEstimates(cm, cr, sh)

test_that("getPPcountsEstimates function works", {
  expect_true(is.data.frame(ccm))
  expect_equal(dim(ccm), c(3, 4))
  expect_equal(round(ccm[1,1], 3), 2.810)
  expect_equal(round(ccm[1,2], 3), 5.895)
  expect_equal(round(ccm[1,3], 3), 12.365)
  expect_equal(round(ccm[2,1], 3), 1.753)
  expect_equal(round(ccm[2,2], 3), 4.211)
  expect_equal(round(ccm[2,3], 3), 10.116)
  expect_equal(round(ccm[3,1], 3), 0.815)
  expect_equal(round(ccm[3,2], 3), 2.526)
  expect_equal(round(ccm[3,3], 3), 7.833)
})
