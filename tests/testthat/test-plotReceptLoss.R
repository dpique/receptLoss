
test_that("plotReceptLoss returns a ggplot", {
  exprMatrNml <- matrix(abs(rnorm(100, mean = 2)), nrow=10)
  exprMatrTum <- matrix(abs(rnorm(100)), nrow=10)
  geneNames <- paste0(letters[1:nrow(exprMatrNml)], 1:nrow(exprMatrNml))
  rownames(exprMatrNml) <- rownames(exprMatrTum) <- geneNames
  nSdBelow <- 2
  minPropPerGroup <- .2
  rl <- receptLoss(exprMatrNml, exprMatrTum, nSdBelow, minPropPerGroup)
  clrs <- c("#E78AC3", "#8DA0CB")
  p1 <- plotReceptLoss(exprMatrNml, exprMatrTum, rl,
                       geneName="g7", clrs=clrs)
  expect_true(object = identical(class(p1), c("gg", "ggplot")))
})
