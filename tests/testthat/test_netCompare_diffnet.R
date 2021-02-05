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

context("permutation test")

diff_perm <- diffnet(assonet2, diffMethod = "permute", nPerm = 100, cores = 1L,
                     adjust = "none")

context("Fisher test")
diff_fisher <- diffnet(assonet2, diffMethod = "fisherTest", adjust = "none")

context("test plot.diffnet")
plot(diff_perm)
plot(diff_fisher)



