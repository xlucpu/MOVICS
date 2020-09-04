library(testthat)
library(MOVICS)
library(InterSIM)

test_that("get optimal cluster number from simulated dataset", {
  sim <- InterSIM(n.sample = 50)
  simdat <- lapply(sim[1:3],t)
  expect_equal(getClustNum(data = simdat)$N.clust, 3)
})
