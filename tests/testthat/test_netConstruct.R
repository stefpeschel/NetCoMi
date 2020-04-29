data("amgut1.filt")
data("amgut2.filt.phy")

context("netConstruct with different association measures")

measures <- c("pearson", "spearman", "bicor", "sparcc", "cclasso", "ccrepe",
              "propr","gcoda", "spieceasi", "spring" )

for(i in 1:length(measures)){

  context(measures[i])

  testnet <- netConstruct(amgut1.filt,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 20),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = measures[i],
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)

  testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector")
}

context("test plot.microNetProps")
plot(testprops)



context("netConstruct with phyloseq object as input")

library(phyloseq)

testnet <- netConstruct(amgut2.filt.phy,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 20),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

testprops<- netAnalyze(testnet,
                       clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector")


#===============================================================================
context("netConstruct with different dissimilarity measures")

measures <- c("euclidean", "bray", "kld", "jeffrey", "ckld",
              "jsd", "aitchison")

for(i in 1:length(measures)){
  context(measures[i])

  testnet1 <- netConstruct(amgut1.filt,
                          filtTax = "totalReads",
                          filtTaxPar = list(totalReads = 1000),
                          filtSamp = "highestFreq",
                          filtSampPar = list(highestFreq = 30),
                          zeroMethod = "none", normMethod = "none",
                          measure = measures[i],
                          sparsMethod = "knn", thresh = 0.3,
                          seed = 20190101)

  testprops<- netAnalyze(testnet1, clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector")
}

context("test plot.microNetProps")
plot(testprops)


#===============================================================================
context("netConstruct with different sparsification methods")

sparsMethod <- c("none", "t-test", "bootstrap", "threshold", "softThreshold")

for(i in 1:length(sparsMethod)){
  context(sparsMethod[i])

  testnet <- netConstruct(amgut1.filt,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 100),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson", adjust = "adaptBH",
                          sparsMethod = sparsMethod[i],
                          thresh = 0.3, softThreshPower = 10,
                          dissFunc = "signed", nboot = 1000, cores = 1,
                          seed = 20190101, logFile = NULL)

  testprops<- netAnalyze(testnet,
                            clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector", normDeg = FALSE,
                            hubQuant = 0.95, lnormFit = FALSE)
}


#===============================================================================
context("netConstruct with different methods for zero replacement")

zeroMethod <- c("none", "pseudo", "multRepl", "alrEM", "bayesMult")

for(i in 1:length(zeroMethod)){
  context(zeroMethod[i])

  testnet <- netConstruct(amgut1.filt,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 20),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = zeroMethod[i], normMethod = "clr",
                          measure = "pearson", adjust = "BY",
                          sparsMethod = "threshold", thresh = 0.3,
                          dissFunc = "signed", nboot = 20, cores = 4,
                          seed = 20190101)

  testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector")
}


#===============================================================================
context("netConstruct with different normalization methods")
normMethod <- c("none", "fractions", "TSS", "CSS", "COM", "rarefy", "clr")
zeroMethod <- c("none", "pseudo", "multRepl")

for(i in 1:length(normMethod)){
  for(z in 1:length(zeroMethod)){
   # context(paste0("normMethod: ", normMethod[i], "; zeroMethod: ", zeroMethod[z]))

    testnet <- netConstruct(amgut1.filt,
                            filtTax = "highestVar",
                            filtTaxPar = list(highestVar = 20),
                            filtSamp = "totalReads",
                            filtSampPar = list(totalReads = 1000),
                            zeroMethod = zeroMethod[z],
                            normMethod = normMethod[i],
                            measure = "pearson", adjust = "BY",
                            sparsMethod = "threshold", thresh = 0.3,
                            dissFunc = "signed", nboot = 20, cores = 1,
                            seed = 20190101)

    testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector")
  }

}


################################################################################
### two networks
################################################################################

groups_diss <- sample(0:1, ncol(amgut1.filt), replace = TRUE)
groups_asso <- sample(0:1, nrow(amgut1.filt), replace = TRUE)

context("netConstruct with different association measures")

measures <- c("pearson", "spearman", "bicor", "sparcc", "cclasso", "ccrepe",
              "propr","gcoda", "spieceasi", "spring" )

for(i in 1:length(measures)){
  context(measures[i])

  testnet <- netConstruct(amgut1.filt, group = groups_asso,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 20),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = measures[i],
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)

  testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector")
}

context("test plot.microNetProps")
plot(testprops)

#------------------------------------------------------------------------------
context("netConstruct with different dissimilarity measures")

measures <- c("euclidean", "bray", "kld", "jeffrey", "ckld",
              "jsd", "aitchison")

for(i in 1:length(measures)){
  context(measures[i])

  testnet1 <- netConstruct(amgut1.filt, group = groups_diss,
                           filtTax = "totalReads",
                           filtTaxPar = list(totalReads = 1000),
                           filtSamp = "highestFreq",
                           filtSampPar = list(highestFreq = 30),
                           zeroMethod = "none", normMethod = "none",
                           measure = measures[i],
                           sparsMethod = "knn", thresh = 0.3,
                           seed = 20190101)

  testprops<- netAnalyze(testnet1, clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector")
}

#===============================================================================
context("netConstruct with different sparsification methods")

sparsMethod <- c("none", "t-test", "bootstrap", "threshold", "softThreshold")

for(i in 1:length(sparsMethod)){
  context(sparsMethod[i])

  testnet <- netConstruct(amgut1.filt, group = groups_asso,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 100),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson", adjust = "adaptBH",
                          sparsMethod = sparsMethod[i], thresh = 0.3,
                          softThreshPower = c(8,10),
                          dissFunc = "signed", nboot = 1000, cores = 1,
                          seed = 20190101, logFile = NULL)

  testprops<- netAnalyze(testnet,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector", normDeg = FALSE,
                         hubQuant = 0.95, lnormFit = FALSE)
}


#===============================================================================
context("netConstruct with different methods for zero replacement")

zeroMethod <- c("none", "pseudo", "multRepl", "alrEM", "bayesMult")

for(i in 1:length(zeroMethod)){
  context(zeroMethod[i])

  testnet <- netConstruct(amgut1.filt,  group = groups_asso,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 20),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = zeroMethod[i], normMethod = "clr",
                          measure = "pearson", adjust = "BY",
                          sparsMethod = "threshold", thresh = 0.3,
                          dissFunc = "signed", nboot = 20, cores = 4,
                          seed = 20190101)

  testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector")
}


#===============================================================================
context("netConstruct with different normalization methods")
normMethod <- c("none", "fractions", "TSS", "CSS", "COM", "rarefy", "clr")
zeroMethod <- c("none", "pseudo", "multRepl")

for(i in 1:length(normMethod)){
  for(z in 1:length(zeroMethod)){

    context(paste0("normMethod: ", normMethod[i], "; zeroMethod: ", zeroMethod[z]))

    testnet <- netConstruct(amgut1.filt, group = groups_asso,
                            filtTax = "highestVar",
                            filtTaxPar = list(highestVar = 20),
                            filtSamp = "totalReads",
                            filtSampPar = list(totalReads = 1000),
                            zeroMethod = zeroMethod[z],
                            normMethod = normMethod[i],
                            measure = "pearson", adjust = "BY",
                            sparsMethod = "threshold", thresh = 0.3,
                            dissFunc = "signed", nboot = 20, cores = 4,
                            seed = 20190101)

    testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector")
  }
}

