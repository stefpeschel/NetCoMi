set.seed(123456)
data("amgut1.filt")

groups_diss <- sample(0:1, ncol(amgut1.filt), replace = TRUE)
groups_asso <- sample(0:1, nrow(amgut1.filt), replace = TRUE)

context("test netCompare")

assonet2 <- netConstruct(amgut1.filt, group = groups_asso,
                         filtTax = "highestVar",
                         filtTaxPar = list(highestVar = 50),
                         filtSamp = "totalReads",
                         filtSampPar = list(totalReads = 1000),
                         zeroMethod = "none", normMethod = "none",
                         measure = "pearson",
                         sparsMethod = "threshold", thresh = 0.3,
                         dissFunc = "signed",
                         seed = 20190101)

dissnet2 <- netConstruct(amgut1.filt, group = groups_diss,
                         filtTax = "totalReads",
                         filtTaxPar = list(highestVar = 1000),
                         filtSamp = "highestFreq",
                         filtSampPar = list(totalReads = 50),
                         zeroMethod = "none", normMethod = "none",
                         measure = "bray",
                         sparsMethod = "threshold", thresh = 0.3,
                         dissFunc = "signed",
                         seed = 20190101)

assoprops2 <- netAnalyze(assonet2, clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector")

dissprops2 <- netAnalyze(dissnet2, clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector")


context("association network; without permutation test")
netcomp_asso <- netCompare(assoprops2, permTest = FALSE)

context("association network; with permutation test")
netcomp_asso <- netCompare(assoprops2, permTest = TRUE, nPerm = 100, cores = 1L)

context("dissimilarity network; without permutation test")
netcomp_diss <- netCompare(dissprops2, permTest = FALSE)

context("dissimilarity network; with permutation test")
netcomp_diss <- netCompare(dissprops2, permTest = TRUE, nPerm = 100, cores = 1L)

context("test summary method")

summary(netcomp_asso)
summary(netcomp_diss)

#-------------------------------------------------------------------------------
context("differential network")
set.seed(123456)

context("permutation test")
diff_perm <- diffnet(assonet2, diffMethod = "permute", nPerm = 20, cores = 1L,
                     adjust = "none")

context("Fisher test")
diff_fisher <- diffnet(assonet2, diffMethod = "fisherTest", adjust = "none")

context("Discordant method")
diff_discord <- diffnet(assonet2, diffMethod = "discordant", adjust = "none")

context("test plot.diffnet")
plot(diff_perm)
plot(diff_fisher)
plot(diff_discord)


assonet3 <- netConstruct(assonet2$assoEst1,
                         data2 = assonet2$assoEst2,
                         dataType = "correlation",
                         filtTax = "highestVar",
                         filtTaxPar = list(highestVar = 50),
                         filtSamp = "totalReads",
                         filtSampPar = list(totalReads = 1000),
                         zeroMethod = "none", normMethod = "none",
                         measure = "pearson",
                         sparsMethod = "threshold", thresh = 0.3,
                         dissFunc = "signed",
                         seed = 20190101)

context("Fisher test")
diff_fisher <- diffnet(assonet3, diffMethod = "fisherTest", adjust = "none",
                       n1 = sum(groups_asso == 0), n2 = sum(groups_asso == 1))



#-------------------------------------------------------------------------------
context("Small sample size")

suppressWarnings(assonet2 <- netConstruct(amgut1.filt[c(1,2,3,6), ],
                                          group = groups_asso[c(1,2,3,6)],
                                          filtTax = "highestVar",
                                          filtTaxPar = list(highestVar = 50),
                                          zeroMethod = "none", normMethod = "none",
                                          measure = "pearson",
                                          sparsMethod = "threshold", thresh = 0.3,
                                          dissFunc = "signed",
                                          seed = 20190101))

assoprops2 <- netAnalyze(assonet2, clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector")

netcomp_asso <- netCompare(assoprops2, permTest = FALSE)

expect_error(netcomp_asso <- netCompare(assoprops2,
                                        permTest = TRUE,
                                        nPerm = 100, cores = 1L),
             "Possible number of permutations")

expect_error(diff_perm <- diffnet(assonet2,
                                  diffMethod = "permute",
                                  nPerm = 100, cores = 1L,
                                  adjust = "none"),
             "Possible number of permutations")


#-------------------------------------------------------------------------------
context("diffnet adjustment over non-zero associations")

make_diffnet_test_net <- function() {
  est1 <- est2 <- diag(4)
  dimnames(est1) <- dimnames(est2) <- list(paste0("t", 1:4), paste0("t", 1:4))

  est1[1, 2] <- est1[2, 1] <- 0.1
  est2[1, 2] <- est2[2, 1] <- 0.1

  est1[1, 3] <- est1[3, 1] <- -0.4
  est2[1, 3] <- est2[3, 1] <- 0.05

  est1[1, 4] <- est1[4, 1] <- 0.05
  est2[1, 4] <- est2[4, 1] <- 0.5

  est1[2, 3] <- est1[3, 2] <- 0.2
  est2[2, 3] <- est2[3, 2] <- 0.2

  est1[2, 4] <- est1[4, 2] <- 0.1
  est2[2, 4] <- est2[4, 2] <- 0.1

  est1[3, 4] <- est1[4, 3] <- 0.1
  est2[3, 4] <- est2[4, 3] <- 0.1

  sparse1 <- est1
  sparse2 <- est2
  sparse1[1, 2] <- sparse1[2, 1] <- 0
  sparse2[1, 2] <- sparse2[2, 1] <- 0
  sparse1[2, 4] <- sparse1[4, 2] <- 0
  sparse2[2, 4] <- sparse2[4, 2] <- 0
  sparse1[3, 4] <- sparse1[4, 3] <- 0
  sparse2[3, 4] <- sparse2[4, 3] <- 0

  out <- list(assoEst1 = est1,
              assoEst2 = est2,
              assoMat1 = sparse1,
              assoMat2 = sparse2,
              assoType = "association",
              parameters = list(dataType = "correlation"),
              normCounts1 = matrix(0, 20, 4),
              normCounts2 = matrix(0, 20, 4),
              countMat1 = matrix(0, 20, 4),
              countMat2 = matrix(0, 20, 4),
              groups = list(group1 = "1", group2 = "2"))
  class(out) <- "microNet"

  out
}

test_that("Fisher adjustment can exclude zero-zero association pairs", {
  diff_net <- diffnet(make_diffnet_test_net(),
                      diffMethod = "fisherTest",
                      adjust = "BH",
                      adjustNonzeroOnly = TRUE,
                      n1 = 20,
                      n2 = 20,
                      verbose = FALSE)

  expect_equal(diff_net$nTestsTotal, 6)
  expect_equal(diff_net$nTestsAdjusted, 3)
  expect_equal(sum(!diff_net$adjustIncludeVec), 3)
  expect_true(diff_net$adjustIncludeVec[["t3_t1"]])
  expect_false(diff_net$adjustIncludeVec[["t2_t1"]])
  expect_true(is.na(diff_net$pAdjustVec[["t2_t1"]]))
  expect_true(is.na(diff_net$pAdjustMat["t1", "t2"]))
  expect_equal(diff_net$diffAdjustMat["t1", "t2"], 0)
  expect_equal(diff_net$diffAdjustMat["t2", "t4"], 0)
  expect_equal(diff_net$diffAdjustMat["t3", "t4"], 0)
})

test_that("supplied permutation p-values use the non-zero adjustment mask", {
  net <- make_diffnet_test_net()
  net$parameters$dataType <- "counts"

  diff_net <- diffnet(net,
                      diffMethod = "permute",
                      adjust = "BH",
                      adjustNonzeroOnly = TRUE,
                      pvalsVec = c(0.01, 0.02, 0.03, 0.9, 0.04, 0.05),
                      verbose = FALSE)

  expect_equal(diff_net$nTestsTotal, 6)
  expect_equal(diff_net$nTestsAdjusted, 3)
  expect_true(diff_net$adjustIncludeVec[["t3_t1"]])
  expect_true(is.na(diff_net$pAdjustVec[["t2_t1"]]))
  expect_true(is.na(diff_net$pAdjustVec[["t4_t2"]]))
  expect_equal(diff_net$diffAdjustMat["t1", "t2"], 0)
})

test_that("diffnet summaries print for p-value, threshold, and discordant outputs", {
  net <- make_diffnet_test_net()

  diff_fisher <- diffnet(net,
                         diffMethod = "fisherTest",
                         adjust = "BH",
                         adjustNonzeroOnly = TRUE,
                         n1 = 20,
                         n2 = 20,
                         verbose = FALSE)
  expect_is(summary(diff_fisher), "summary.diffnet")
  expect_output(print(summary(diff_fisher)), "Differential Network Summary")

  diff_thresh <- diffnet(net,
                         diffMethod = "threshold",
                         diffThresh = 0.1,
                         verbose = FALSE)
  expect_is(summary(diff_thresh), "summary.diffnet")
  expect_output(print(summary(diff_thresh)), "Method-specific summary")

  diff_discord <- list(assoMat1 = net$assoEst1,
                       assoMat2 = net$assoEst2,
                       diffMat = abs(net$assoEst1 - net$assoEst2),
                       diffProbs = abs(net$assoEst1 - net$assoEst2),
                       diffMethod = "discordant",
                       adjusted = "none",
                       call = quote(diffnet(net, diffMethod = "discordant")))
  class(diff_discord) <- "diffnet"
  expect_is(summary(diff_discord), "summary.diffnet")
  expect_output(print(summary(diff_discord)), "Method-specific summary")
})

test_that("plot.diffnet returns with a message if there is nothing to plot", {
  net <- make_diffnet_test_net()
  diff_empty <- list(assoMat1 = net$assoEst1,
                     assoMat2 = net$assoEst2,
                     diffMat = matrix(0, 4, 4,
                                      dimnames = dimnames(net$assoEst1)),
                     diffAdjustMat = matrix(0, 4, 4,
                                            dimnames = dimnames(net$assoEst1)),
                     diffMethod = "fisherTest",
                     adjusted = "BH")
  class(diff_empty) <- "diffnet"

  expect_message(out <- plot(diff_empty),
                 "after multiple testing adjustment")
  expect_null(out)
})

test_that("plot.diffnet matches named layouts after node filtering", {
  net <- make_diffnet_test_net()
  diff_filtered <- list(assoMat1 = net$assoEst1,
                        assoMat2 = net$assoEst2,
                        diffMat = abs(net$assoEst1 - net$assoEst2),
                        diffAdjustMat = abs(net$assoEst1 - net$assoEst2),
                        diffMethod = "fisherTest",
                        adjusted = "BH",
                        groups = net$groups)
  class(diff_filtered) <- "diffnet"

  layout <- matrix(c(0, 0,
                     1, 0,
                     0, 1),
                   ncol = 2,
                   byrow = TRUE,
                   dimnames = list(c("t1", "t3", "t4"), NULL))
  labels <- setNames(c("one", "two", "three", "four"), paste0("t", 1:4))

  expect_silent(out <- plot(diff_filtered,
                            layout = layout,
                            labels = labels,
                            nodeFilter = "names",
                            nodeFilterPar = c("t1", "t3", "t4"),
                            rmSingles = FALSE,
                            legend = FALSE))
  expect_equal(rownames(out), c("t1", "t3", "t4"))
})

test_that("plot.diffnet lifts weak edge weights for width scaling", {
  net <- make_diffnet_test_net()
  weight <- .scaleDiffnetEdgesForPlot(abs(net$assoEst1 - net$assoEst2))

  expect_equal(max(weight), 1)
  expect_equal(weight["t1", "t3"], 1)
  expect_true(weight["t1", "t4"] >= 0.25)
  expect_true(weight["t1", "t4"] < 1)
})



