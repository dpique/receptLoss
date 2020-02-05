library(tidyverse)

test_that("receptLoss returns a tibble", {
  exprMatrNml <- matrix(abs(rnorm(100, mean = 2)), nrow=10)
  exprMatrTum <- matrix(abs(rnorm(100)), nrow=10)
  geneNames <- paste0(letters[1:nrow(exprMatrNml)], 1:nrow(exprMatrNml))
  rownames(exprMatrNml) <- rownames(exprMatrTum) <- geneNames
  nSdBelow <- 2
  minPropPerGroup <- .2
  rl <- receptLoss(exprMatrNml, exprMatrTum, nSdBelow, minPropPerGroup)
  expect_true(object = is_tibble(rl))
})
