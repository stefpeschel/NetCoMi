set.seed(123456)

data("amgut1.filt")
data("amgut2.filt.phy")


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
# Filtering

context("netConstruct with different taxa filtering methods (single network)")

ftax <- c("none", "totalReads", "relFreq", "numbSamp", "highestVar", "highestFreq")

ftaxpar <- c(list(totalReads = 1000), list(relFreq = 0.05), 
             list(numbSamp = 10), list(highestVar = 50), list(highestFreq = 50))

dims <- c(127, 125, 3, 127, 50, 50)

for(i in 1:length(ftax)){
  testnet <- netConstruct(amgut1.filt,
                          filtTax = ftax[i],
                          filtTaxPar = ftaxpar[i-1],
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)
  expect_equal(dim(testnet$normCounts1)[2], dims[i])
}

#-------------------------------------------------------------------------------
context("netConstruct with different sample filtering methods (single network")

fsamp <- c("none", "totalReads", "numbTaxa", "highestFreq")

fsampar <- c(list(totalReads = 1000), list(numbTaxa = 50),list(highestFreq = 50))

dims <- c(289, 289, 282, 50)

for(i in 1:length(fsamp)){
  testnet <- netConstruct(amgut1.filt,
                          filtSamp = fsamp[i],
                          filtSampPar = fsampar[i-1],
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)
  expect_equal(dim(testnet$normCounts1)[1], dims[i])
}

#===============================================================================
# Association / dissimilarity measures

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


#-------------------------------------------------------------------------------
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
# Zero handling

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
# Normalization + zero handling

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
                            dissFunc = "signed", nboot = 20, 
                            seed = 20190101)
    
    testprops<- netAnalyze(testnet, clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector")
    
    plot(testprops)
    mtext(paste0(normMethod[i], ", ", zeroMethod[z]), side = 3, cex = 1.5)
  }
  
}

#===============================================================================
# Sparsification

context("netConstruct with different sparsification methods")

sparsMethod <- c("none", "t-test", #"bootstrap", 
                 "threshold", "softThreshold")

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
                          dissFunc = "signed", nboot = 1000, cores = 4,
                          seed = 20190101, logFile = NULL)

  testprops<- netAnalyze(testnet,
                            clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector", normDeg = FALSE,
                            hubQuant = 0.95, lnormFit = FALSE)

  plot(testprops)
  mtext(sparsMethod[i], side = 3, cex = 1.5)
}


#===============================================================================
# Multiple testing adjustment

context("netConstruct with different adjustment methods")

adjustMethod <- c("lfdr", "holm", "BH", "BY", rep("adaptBH", 5))
trueNullMethod <- c(rep("convest", 5), "lfdr", "mean", "hist", "farco")

for(i in 1:length(adjustMethod)){
  context(adjustMethod[i])
  
  testnet <- netConstruct(amgut1.filt,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 50),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson", 
                          sparsMethod = "t-test",
                          adjust = adjustMethod[i],
                          lfdrThresh = 0.2,
                          trueNullMethod = trueNullMethod[i],
                          dissFunc = "signed", nboot = 1000, 
                          seed = 20190101, logFile = NULL)
  
  testprops<- netAnalyze(testnet,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector", normDeg = FALSE,
                         hubQuant = 0.95, lnormFit = FALSE)
  
  plot(testprops)
  mtext(adjustMethod[i], side = 3, cex = 1.5)
}


#===============================================================================
# Soft-thresholding

context("netConstruct with different values for softThreshType")

softThreshType <- c("signed", "unsigned", "signed hybrid")

for(i in 1:length(softThreshType)){
  context(softThreshType[i])
  
  if(i == 1){
    softThreshPower <- 5
  } else{
    softThreshPower <- NULL
  }
  
  testnet <- netConstruct(amgut1.filt,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 50),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson", 
                          sparsMethod = "softThreshold",
                          softThreshType = softThreshType[i],
                          softThreshPower = softThreshPower,
                          dissFunc = "signed", 
                          seed = 20190101, logFile = NULL)
  
  testprops<- netAnalyze(testnet,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector", normDeg = FALSE,
                         hubQuant = 0.95, lnormFit = FALSE)
  
  plot(testprops)
  mtext(softThreshType[i], side = 3, cex = 1.5)
}

#===============================================================================
# K-nearest neighbor sparsification

context("netConstruct with knn sparsification")

knnMutual <- c(TRUE, FALSE)

for(i in 1:length(knnMutual)){
  context(knnMutual[i])

  testnet <- netConstruct(amgut1.filt,
                          filtSamp = "highestFreq",
                          filtSampPar = list(highestFreq = 50),
                          zeroMethod = "none", normMethod = "none",
                          measure = "aitchison", 
                          sparsMethod = "knn",
                          knnMutual = knnMutual[i],
                          dissFunc = "signed", 
                          seed = 20190101, logFile = NULL)
  
  testprops<- netAnalyze(testnet,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector", normDeg = FALSE,
                         hubQuant = 0.95, lnormFit = FALSE)
  
  plot(testprops)
  mtext(knnMutual[i], side = 3, cex = 1.5)
}

#===============================================================================
# Dissimilarity function

context("netConstruct with differenct values for dissFunc")

dissFunc <- c("signed", "unsigned", "signedPos", "TOMdiss")

for(i in 1:length(dissFunc)){
  context(dissFunc[i])
  
  testnet <- netConstruct(amgut1.filt,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 50),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson", 
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          dissFunc = dissFunc[i],
                          seed = 20190101)
  
  testprops<- netAnalyze(testnet,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector", normDeg = FALSE,
                         hubQuant = 0.95, lnormFit = FALSE)

  plot(testprops)
  mtext(dissFunc[i], side = 3, cex = 1.5)
}

#-------------------------------------------------------------------------------
# 'dissFunc' is a function

testfunc <- function(x){
  xvec <- x[lower.tri(x)]
  dissvec <- sqrt(0.5 * (1-xvec))
  
  dissMat <- x
  dissMat[lower.tri(dissMat)] <- dissvec
  dissMat[upper.tri(dissMat)] <- t(dissMat)[upper.tri(t(dissMat))]
  dissMat[x == 0] <- Inf
  diag(dissMat) <- 0
  
  return(dissMat)
}


context("'dissFunc' is a function")

testnet <- netConstruct(amgut1.filt,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson", 
                        sparsMethod = "threshold",
                        thresh = 0.3,
                        dissFunc = testfunc,
                        seed = 20190101)

testprops<- netAnalyze(testnet,
                       clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector", normDeg = FALSE,
                       hubQuant = 0.95, lnormFit = FALSE)

plot(testprops)
mtext("dissFunc as function", side = 3, cex = 1.5)

#===============================================================================
# Similarity function

testfunc <- function(x, power){
  1/(1 + x^power) 
}

context("simFunc")

testnet <- netConstruct(amgut1.filt,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson", 
                        sparsMethod = "threshold",
                        thresh = 0.3,
                        dissFunc = "signed",
                        simFunc = testfunc,
                        simFuncPar = list(2),
                        seed = 20190101)

testprops<- netAnalyze(testnet,
                       clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector", normDeg = FALSE,
                       hubQuant = 0.95, lnormFit = FALSE)

plot(testprops)
mtext("simFunc", side = 3, cex = 1.5)

#===============================================================================
# scaleDiss

context("scaleDiss")

scaleDiss <- c(TRUE, FALSE)

for(i in 1:length(scaleDiss)){
  context(scaleDiss[i])
  
  testnet <- netConstruct(amgut1.filt,
                          filtSamp = "highestFreq",
                          filtSampPar = list(highestFreq = 50),
                          zeroMethod = "none", normMethod = "none",
                          measure = "aitchison", 
                          sparsMethod = "knn",
                          dissFunc = "signed", 
                          scaleDiss = scaleDiss[i],
                          seed = 20190101, logFile = NULL)
  
  testprops<- netAnalyze(testnet,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector", normDeg = FALSE,
                         hubQuant = 0.95, lnormFit = FALSE)

  plot(testprops)
  mtext(scaleDiss[i], side = 3, cex = 1.5)
}

#===============================================================================
# Unweighted network

context("unweighted network")

testnet <- netConstruct(amgut1.filt,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson", 
                        sparsMethod = "threshold",
                        thresh = 0.3,
                        weighted = FALSE,
                        seed = 20190101)

testprops<- netAnalyze(testnet,
                       clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector", normDeg = FALSE,
                       hubQuant = 0.95, lnormFit = FALSE)

summary(testprops)

plot(testprops)
mtext("weighted", side = 3, cex = 1.5)



################################################################################
### two networks
################################################################################
set.seed(123456)
groups_diss <- sample(0:1, ncol(amgut1.filt), replace = TRUE)
groups_asso <- sample(0:1, nrow(amgut1.filt), replace = TRUE)

amgut_male <- phyloseq::subset_samples(amgut2.filt.phy, SEX == "male")
amgut_female <- phyloseq::subset_samples(amgut2.filt.phy, SEX == "female")

context("netConstruct with group vector (association)")
testnet <- netConstruct(amgut1.filt, group = groups_asso,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)


context("netConstruct with data and data2 (association)")
testnet <- netConstruct(amgut_male, amgut_female,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

context("netConstruct with group vector (dissimilarity)")
testnet <- netConstruct(amgut1.filt, group = groups_diss,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "aitchison",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

context("netConstruct with data and data2 (dissimilarity)")
amgut_diss1 <- amgut1.filt[, groups_diss == 0]
amgut_diss2 <- amgut1.filt[, groups_diss == 1]

testnet <- netConstruct(amgut_diss1, data2 = amgut_diss2,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "aitchison",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

#===============================================================================
# Matched-group designs

context("netConstruct with matched-group design")

testnet <- netConstruct(amgut1.filt[1:140, ], amgut1.filt[141:280, ],
                        matchDesign = c(1,1),
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

expect_error(
testnet <- netConstruct(amgut1.filt[1:140, ], amgut1.filt[141:281, ],
                        matchDesign = c(1,1),
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101))

testnet <- netConstruct(amgut1.filt[1:95, ], amgut1.filt[96:285, ],
                        matchDesign = c(1,2),
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

expect_error(
testnet <- netConstruct(amgut1.filt[1:140, ], amgut1.filt[141:280, ],
                        matchDesign = c(1,2),
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101))

#===============================================================================
# jointPrepro

context("netConstruct with jointPrepro")

testnet <- netConstruct(amgut1.filt, group = groups_asso,
                        jointPrepro = FALSE,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

testnet <- netConstruct(amgut_male, amgut_female,
                        jointPrepro = TRUE,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101)

expect_error(
testnet <- netConstruct(amgut1.filt, group = groups_diss,
                        jointPrepro = TRUE,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "aitchison",
                        sparsMethod = "threshold", thresh = 0.3,
                        seed = 20190101))


#===============================================================================
# Filtering

##context("netConstruct with different taxa filtering methods (two networks)")

ftax <- c("none", "totalReads", "relFreq", "numbSamp", "highestVar", "highestFreq")

ftaxpar <- c(list(totalReads = 1000), list(relFreq = 0.05), 
             list(numbSamp = 10), list(highestVar = 50), list(highestFreq = 50))

dims <- c(127, 125, 3, 127, 50, 50)

for(i in 1:length(ftax)){
  testnet <- netConstruct(amgut1.filt, group = groups_asso,
                          filtTax = ftax[i],
                          filtTaxPar = ftaxpar[i-1],
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)
  expect_equal(dim(testnet$normCounts1)[2], dims[i])
}

dims <- c(138, 123, 2, 138, 41, 42)

for(i in 1:length(ftax)){
  testnet <- netConstruct(amgut_male, amgut_female,
                          filtTax = ftax[i],
                          filtTaxPar = ftaxpar[i-1],
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)
  expect_equal(dim(testnet$normCounts1)[2], dims[i])
}


# dissimilarity networks
dims <- c(49, 49, 6, 49, 49, 49)
for(i in 1:length(ftax)){
  testnet <- netConstruct(amgut1.filt, group = groups_diss,
                          filtTax = ftax[i],
                          filtTaxPar = ftaxpar[i-1],
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "aitchison",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)
  
  expect_equal(dim(testnet$normCounts1)[2], dims[i])
}


#-------------------------------------------------------------------------------
##context("netConstruct with different sample filtering methods (two networks")

fsamp <- c("none", "totalReads", "numbTaxa", "highestFreq")

fsampar <- c(list(totalReads = 1000), list(numbTaxa = 30),list(highestFreq = 50))

dims <- c(289, 289, 288, 50)

for(i in 1:length(fsamp)){
  testnet <- netConstruct(amgut1.filt, group = groups_asso,
                          filtSamp = fsamp[i],
                          filtSampPar = fsampar[i-1],
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)
  
  expect_equal(dim(testnet$normCounts1)[1] + dim(testnet$normCounts2)[1], dims[i])
}

dims1 <- c(146, 127, 128, 50)
dims2 <- c(126, 114, 114, 50)

for(i in 1:length(fsamp)){
  testnet <- netConstruct(amgut_female, amgut_male,
                          filtSamp = fsamp[i],
                          filtSampPar = fsampar[i-1],
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)
  
  expect_equal(dim(testnet$normCounts1)[1], dims1[i])
  expect_equal(dim(testnet$normCounts2)[1], dims2[i])
}


dims <- c(289, 270, 228, 2)

for(i in 1:length(fsamp)){
  testnet <- netConstruct(amgut1.filt, group = groups_diss,
                          filtSamp = fsamp[i],
                          filtSampPar = fsampar[i-1],
                          zeroMethod = "none", normMethod = "none",
                          measure = "aitchison",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)
  
  expect_equal(dim(testnet$normCounts1)[1], dims[i])
}

#===============================================================================
# Association / dissimilarity measures

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
# Zero handling

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
# Normalization + zero handling

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


#===============================================================================
# Sparsification

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
                          dissFunc = "signed", nboot = 1000, 
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
# SoftThreshPower

context("netConstruct with given softThreshPower")

testnet <- netConstruct(amgut1.filt,
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        filtSamp = "totalReads",
                        filtSampPar = list(totalReads = 1000),
                        zeroMethod = "none", normMethod = "none",
                        measure = "pearson", 
                        sparsMethod = "softThreshold",
                        softThreshType = "signed",
                        softThreshPower = c(5,3),
                        dissFunc = "signed", 
                        seed = 20190101, logFile = NULL)

testprops<- netAnalyze(testnet,
                       clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector", normDeg = FALSE,
                       hubQuant = 0.95, lnormFit = FALSE)

plot(testprops)






