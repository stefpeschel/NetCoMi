data("amgut1.filt")
data("amgut2.filt.phy")

context("netConstruct with different association measures")

measures <- c("pearson", "spearman", "bicor", "sparcc", "cclasso", "ccrepe",
              "propr","gcoda", "spieceasi_gl", "spieceasi_mb", "spring" )

for(i in 1:length(measures)){

  context(measures[i])
  
  measure.tmp <- measures[i]
  
  if(measure.tmp == "spieceasi_gl"){
    measurePar <- list(method = "glasso",
                      nlambda=5,
                      pulsar.params = list(rep.num=5))
    measure.tmp <- "spieceasi"
  } else if(measure.tmp == "spieceasi_mb"){
    measurePar <- list(method = "mb",
                       nlambda=5,
                       pulsar.params = list(rep.num=5))
    measure.tmp <- "spieceasi"
  } else if(measure.tmp == "spring"){
    measurePar <- list(nlambda = 5, rep.num = 5)
  } else if(measure.tmp == "gcoda"){
    measurePar <- list(nlambda = 5)
  } else{
    measurePar <- NULL
  }
  
  testnet <- netConstruct(amgut1.filt,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 20),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = measure.tmp,
                          measurePar = measurePar,
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)

  testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector")
  
  plot(testprops)
  mtext(measures[i], side = 3, cex = 1.5)
}

context("test plot.microNetProps")
plot(testprops)


context("netConstruct with association matrix as input")

testnet1 <- netConstruct(amgut1.filt,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 20),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

testthat::expect_that(netConstruct(testnet1$assoEst1, dataType = "correlation",
                         sparsMethod = "threshold", thresh = 0.3,
                         seed = 20190101)$assoMat1, 
                      testthat::equals(testnet1$assoMat1))

testthat::expect_that(netConstruct(testnet1$assoEst1, dataType = "correlation",
                         sparsMethod = "threshold", thresh = 0.3,
                         seed = 20190101)$adjaMat1, 
                      testthat::equals(testnet1$adjaMat1))


testnet1 <- netConstruct(amgut1.filt,
                         filtTax = "totalReads",
                         filtTaxPar = list(totalReads = 1000),
                         filtSamp = "highestFreq",
                         filtSampPar = list(highestFreq = 30),
                         zeroMethod = "none", normMethod = "none",
                         measure = "ckld",
                         sparsMethod = "knn", seed = 20190101)

testthat::expect_that(netConstruct(testnet1$dissEst1, dataType = "dissimilarity",
                         sparsMethod = "knn", seed = 20190101)$dissMat1, 
                      testthat::equals(testnet1$dissMat1))

testthat::expect_that(netConstruct(testnet1$dissEst1, dataType = "dissimilarity",
                         sparsMethod = "knn", seed = 20190101)$adjaMat1, 
                      testthat::equals(testnet1$adjaMat1))


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

plot(testprops)
mtext("phyloseq", side = 3, cex = 1.5)


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
  
  plot(testprops)
  mtext(measures[i], side = 3, cex = 1.5)
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
  
  plot(testprops)
  mtext(sparsMethod[i], side = 3, cex = 1.5)
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
  
  plot(testprops)
  mtext(zeroMethod[i], side = 3, cex = 1.5)
}


#===============================================================================
context("netConstruct with different normalization methods")
normMethod <- c("none", "fractions", "TSS", "CSS", "COM", "rarefy", "clr", "mclr")
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
    
    plot(testprops)
    mtext(paste0(normMethod[i], ", ", zeroMethod[z]), side = 3, cex = 1.5)
  }

}


################################################################################
### two networks
################################################################################

groups_diss <- sample(0:1, ncol(amgut1.filt), replace = TRUE)
groups_asso <- sample(0:1, nrow(amgut1.filt), replace = TRUE)

context("netConstruct with different association measures")

measures <- c("pearson", "spearman", "bicor", "sparcc", "cclasso", "ccrepe",
              "propr","gcoda", "spieceasi_gl", "spieceasi_mb", "spring" )

for(i in 1:length(measures)){
  
  context(measures[i])
  
  measure.tmp <- measures[i]
  
  if(measure.tmp == "spieceasi_gl"){
    measurePar <- list(method = "glasso",
                       nlambda=5,
                       pulsar.params = list(rep.num=5))
    measure.tmp <- "spieceasi"
  } else if(measure.tmp == "spieceasi_mb"){
    measurePar <- list(method = "mb",
                       nlambda=5,
                       pulsar.params = list(rep.num=5))
    measure.tmp <- "spieceasi"
  } else if(measure.tmp == "spring"){
    measurePar <- list(nlambda = 5, rep.num = 5)
  } else if(measure.tmp == "gcoda"){
    measurePar <- list(nlambda = 5)
  } else{
    measurePar <- NULL
  }
  
  testnet <- netConstruct(amgut1.filt, group = groups_asso,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 20),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = measure.tmp,
                          measurePar = measurePar,
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)

  testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector")
  
  plot(testprops)
  mtext(measures[i], side = 3, cex = 1.5)
  
  if(i < 7){
    netcomp_asso <- netCompare(testprops, permTest = TRUE, nPerm = 8, cores = 1L)
  }
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
                           filtSampPar = list(highestFreq = 100),
                           zeroMethod = "none", normMethod = "none",
                           measure = measures[i],
                           sparsMethod = "knn", thresh = 0.3,
                           seed = 20190101)

  testprops<- netAnalyze(testnet1, clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector")
  
  plot(testprops)
  mtext(measures[i], side = 3, cex = 1.5)
  
  netcomp_asso <- netCompare(testprops, permTest = TRUE, nPerm = 8, cores = 1L)
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
  
  plot(testprops)
  mtext(sparsMethod[i], side = 3, cex = 1.5)
  
  if(i != 3){
    netcomp_asso <- netCompare(testprops, permTest = TRUE, nPerm = 8, cores = 1L)
  }
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
  
  plot(testprops)
  mtext(zeroMethod[i], side = 3, cex = 1.5)
  
  netcomp_asso <- netCompare(testprops, permTest = TRUE, nPerm = 8, cores = 1L)
}


#===============================================================================
context("netConstruct with different normalization methods")
normMethod <- c("none", "fractions", "TSS", "CSS", "COM", "rarefy", "clr", "mclr")
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
    
    plot(testprops)
    mtext(paste0(normMethod[i], ", ", zeroMethod[z]), side = 3, cex = 1.5)
    
    netcomp_asso <- netCompare(testprops, permTest = TRUE, nPerm = 8, cores = 1L)
  }
}

