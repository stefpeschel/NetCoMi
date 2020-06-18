data("amgut1.filt")
data("amgut2.filt.phy")

groups_diss <- sample(0:1, ncol(amgut1.filt), replace = TRUE)
groups_asso <- sample(0:1, nrow(amgut1.filt), replace = TRUE)

context("netAnalyze with different clustering methods")

clusterMethods <- c("none", "hierarchical",
                    "cluster_fast_greedy",
                    "cluster_leading_eigen",
                    "cluster_optimal")


testnet <- netConstruct(amgut1.filt, group = groups_asso,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)


for(i in 1:length(clusterMethods)){

  context(clusterMethods[i])


  testprops<- netAnalyze(testnet,
                         clustMethod = clusterMethods[i],
                         hubPar = "eigenvector", hubQuant = 0.95)
}

summary(testprops, numbNodes = 10, showCentr = "degree", digits = 4)


context("netAnalyze with dissimilarity-based network")
testnet <- netConstruct(amgut1.filt, group = groups_diss,
                        filtTax = "totalReads",
                        filtTaxPar = list(totalReads = 1000),
                        filtSamp = "highestFreq",
                        filtSampPar = list(highestFreq = 30),
                        zeroMethod = "multRepl", normMethod = "none",
                        measure = "aitchison",
                        sparsMethod = "knn",
                        seed = 20190101)

testprops<- netAnalyze(testnet, clustMethod = "hierarchical",
                       hubPar = "eigenvector")

#===============================================================================
context("netAnalyze with and without lnormFit")

testnet <- netConstruct(amgut1.filt, group = groups_asso,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 30),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        dissFunc = "signed",
                        seed = 20190101)

testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector",
                       lnormFit = FALSE)

testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector",
                       lnormFit = TRUE)

#===============================================================================
context("netAnalyze with and without weighted degree")

testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector",
                       weightDeg = FALSE)

testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector",  hubQuant = 0.95,
                       weightDeg = TRUE)


#===============================================================================
context("netAnalyze with non-normalized centralities")

testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector", normDeg = FALSE,
                       normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
                       hubQuant = 0.95)





