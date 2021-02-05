calc_diff_props <- function(adja1, adja2, dissMat1, dissMat2, assoMat1, assoMat2,
                            avDissIgnoreInf, sPathNorm, sPathAlgo, 
                            connectivity, normNatConnect,
                            weighted, clustMethod, clustPar, clustPar2, 
                            weightClustCoef, hubPar,
                            hubQuant, jaccQuant, lnormFit, weightDeg, 
                            normDeg, normBetw, normClose, normEigen, centrLCC,
                            testJacc = TRUE, jaccTestGreater = FALSE,
                            testRand = TRUE, nPermRand = 1000){

  isempty1 <- all(adja1[lower.tri(adja1)] == 0)
  isempty2 <- all(adja2[lower.tri(adja2)] == 0)

  if(isempty1 & isempty2) stop("There are no connected nodes in both networks.")

  # calculate network properties
  props1 <- calc_props(adjaMat = adja1, dissMat = dissMat1, assoMat = assoMat1, 
                       avDissIgnoreInf = avDissIgnoreInf,
                       sPathNorm = sPathNorm, sPathAlgo = sPathAlgo,
                       connectivity = connectivity,
                       normNatConnect = normNatConnect,
                       weighted = weighted, isempty = isempty1, 
                       clustMethod = clustMethod, clustPar = clustPar, 
                       weightClustCoef = weightClustCoef,
                       hubPar = hubPar, hubQuant = hubQuant, 
                       lnormFit = lnormFit, weightDeg = weightDeg, 
                       normDeg = normDeg, normBetw = normBetw, 
                       normClose = normClose, normEigen = normEigen, 
                       centrLCC = centrLCC,
                       jaccard = TRUE, jaccQuant = jaccQuant)
  
  props2 <- calc_props(adjaMat = adja2, dissMat = dissMat2, assoMat = assoMat2, 
                       avDissIgnoreInf = avDissIgnoreInf,
                       sPathNorm = sPathNorm, sPathAlgo = sPathAlgo,
                       connectivity = connectivity,
                       normNatConnect = normNatConnect,
                       weighted = weighted, isempty = isempty2, 
                       clustMethod = clustMethod, clustPar = clustPar2, 
                       weightClustCoef = weightClustCoef,
                       hubPar = hubPar, hubQuant = hubQuant, 
                       lnormFit = lnormFit, weightDeg = weightDeg, 
                       normDeg = normDeg, normBetw = normBetw, 
                       normClose = normClose, normEigen = normEigen, 
                       centrLCC = centrLCC,
                       jaccard = TRUE, jaccQuant = jaccQuant)
  
  #== differences in network properties ========================================
  
  ### centralities
  diffdeg <- sort(props1$deg - props2$deg[names(props1$deg)])
  diffbetw <- sort(props1$betw - props2$betw[names(props1$betw)])
  diffclose <- sort(props1$close - props2$close[names(props1$close)])
  diffeigen <- sort(props1$eigen - props2$eigen[names(props1$eigen)])

  # absolute differences
  absdiffdeg <- abs(props1$deg - props2$deg[names(props1$deg)])
  absdiffbetw <- abs(props1$betw - props2$betw[names(props1$betw)])
  absdiffclose <- abs(props1$close - props2$close[names(props1$close)])
  absdiffeigen <- abs(props1$eigen - props2$eigen[names(props1$eigen)])

  #--------------------------------------------------------------------------
  ### Average dissimilarity
  
  avDiss1 <- props1$avDiss
  avDiss2 <- props2$avDiss
  
  avDiss1_lcc <- props1$avDiss_lcc
  avDiss2_lcc <- props2$avDiss_lcc
  
  if(is.na(avDiss1)){
    avDiss1 <- 0
    avDiss1_lcc <- 0
  } 
  
  if(is.na(avDiss2)){
    avDiss2 <- 0
    avDiss2_lcc <- 0
  }
  
  diffdiss <- abs(avDiss1 - avDiss2)
  diffdiss_lcc <- abs(avDiss1_lcc - avDiss2_lcc)
  
  #--------------------------------------------------------------------------
  ### Average Path Length (or Mean Distance)

  avPath1 <- props1$avPath
  avPath2 <- props2$avPath
  
  avPath1_lcc <- props1$avPath_lcc
  avPath2_lcc <- props2$avPath_lcc
  
  if(is.na(avPath1)){
    avPath1 <- 0
    avPath1_lcc <- 0
  } 
  
  if(is.na(avPath2)){
    avPath2 <- 0
    avPath2_lcc <- 0
  }

  diffpath <- abs(avPath1 - avPath2)
  diffpath_lcc <- abs(avPath1_lcc - avPath2_lcc)

  #--------------------------------------------------------------------------
  ### Clustering coefficient

  clustCoef1 <- props1$clustCoef
  clustCoef2 <- props2$clustCoef
  
  clustCoef1_lcc <- props1$clustCoef_lcc
  clustCoef2_lcc <- props2$clustCoef_lcc

  if(is.na(clustCoef1)){
    clustCoef1 <- 0
    clustCoef1_lcc <- 0
  } 
  
  if(is.na(clustCoef2)){
    clustCoef2 <- 0
    clustCoef2_lcc <- 0
  } 
  
  diffclustcoef <- abs(clustCoef1 - clustCoef2)
  diffclustcoef_lcc <- abs(clustCoef1_lcc - clustCoef2_lcc)

  #--------------------------------------------------------------------------
  # differential Modularity
  if(clustMethod != "none"){
    modul1 <- props1$modul
    modul2 <- props2$modul
    
    modul1_lcc <- props1$modul_lcc
    modul2_lcc <- props2$modul_lcc
    
    diffmod <- abs(modul1 - modul2)
    diffmod_lcc <- abs(modul1_lcc - modul2_lcc)
    
  } else{
    modul1 <- modul2 <- modul1_lcc <- modul2_lcc <- NA
    diffmod <- diffmod_lcc <- NA
  }

  #--------------------------------------------------------------------------
  if(connectivity){
    # vertex connectivity
    diffvertconnect <- abs(props1$vertconnect - props2$vertconnect)
    diffvertconnect_lcc <- abs(props1$vertconnect_lcc - props2$vertconnect_lcc)
    
    # edge connectivity
    diffedgconnect <- abs(props1$edgeconnect - props2$edgeconnect)
    diffedgconnect_lcc <- abs(props1$edgeconnect_lcc - props2$edgeconnect_lcc)
  } else{
    diffvertconnect <- diffvertconnect_lcc <- NA
    diffedgconnect <- diffedgconnect_lcc <- NA
  }
  
  #--------------------------------------------------------------------------
  # natural connectivity
  diffnatconnect <- abs(props1$natConnect - props2$natConnect)
  diffnatconnect_lcc <- abs(props1$natConnect_lcc - props2$natConnect_lcc)
  
  #--------------------------------------------------------------------------
  # density (relative number of edges)
  diffdensity <- abs(props1$density - props2$density)
  diffdensity_lcc <- abs(props1$density_lcc - props2$density_lcc)
  
  #--------------------------------------------------------------------------
  # positive-to-negative ratio
  diffpep <- abs(props1$pep - props2$pep)
  diffpep_lcc <- abs(props1$pep_lcc - props2$pep_lcc)
  
  #--------------------------------------------------------------------------
  # number of connected components
  diffncomp <- abs(props1$nComp - props2$nComp)
  
  #--------------------------------------------------------------------------
  # size of the largest connected component
  difflccsize <- abs(props1$lccSize - props2$lccSize)
  
  #--------------------------------------------------------------------------
  # relative LCC size
  difflccsizerel <- abs(props1$lccSizeRel - props2$lccSizeRel)
  
  #--------------------------------------------------------------------------
  # adjusted Rand index for measuring similarity between two clusterings
  clust1 <- props1$clust
  clust2 <- props2$clust
  
  if(isempty1 || isempty2){
    randInd <- NA
  } else{
    randInd <- c(value = WGCNA::randIndex(table(clust1, clust2), adjust = TRUE), 
                 pval = NA)
  }

  # significance test for Rand index
  if(testRand){
    randPerm <- numeric(nPermRand)
    for(i in 1:nPermRand){
      clust1.tmp <- gtools::permute(clust1)
      clust2.tmp <- gtools::permute(clust2)
      randPerm[i] <- WGCNA::randIndex(table(clust1.tmp, clust2.tmp), adjust = TRUE)
    }
    randMean <- mean(randPerm)
    randSD <- sd(randPerm)
    normRandPerm <- (randPerm - randMean) / randSD
    normRand <- (randInd[1] - randMean) / randSD
    randInd["pval"] <- (sum(normRandPerm >= abs(normRand)) + 
                          sum(normRandPerm <= -abs(normRand))) / nPermRand
  }

  #--------------------------------------------------------------------------
  # Jaccard Index

  jaccDeg <-    calc_jaccard(props1$topdeg, props2$topdeg, sigTest = testJacc)
  jaccBetw <-   calc_jaccard(props1$topbetw, props2$topbetw, sigTest = testJacc)
  jaccClose <-  calc_jaccard(props1$topclose, props2$topclose, sigTest = testJacc)
  jaccEigen <-  calc_jaccard(props1$topeigen, props2$topeigen, sigTest = testJacc)
  jaccHub <-    calc_jaccard(props1$hubs, props2$hubs, sigTest = testJacc)
  
  output <- list(jaccDeg = jaccDeg,
                 jaccBetw =jaccBetw,
                 jaccClose = jaccClose,
                 jaccEigen = jaccEigen,
                 jaccHub = jaccHub, 
                 randInd = randInd,
                 diffsGlobal = list(diffnComp = diffncomp,
                                    diffavDiss = diffdiss,
                                    diffavPath = diffpath,
                                    diffDensity = diffdensity,
                                    diffVertConnect = diffvertconnect,
                                    diffEdgeConnect = diffedgconnect,
                                    diffNatConnect = diffnatconnect,
                                    diffPEP = diffpep,
                                    diffClustCoef = diffclustcoef,
                                    diffModul = diffmod),
                 diffsGlobalLCC = list(difflccSize = difflccsize,
                                       difflccSizeRel = difflccsizerel,
                                       diffavDiss = diffdiss_lcc,
                                       diffavPath = diffpath_lcc,
                                       diffDensity = diffdensity_lcc,
                                       diffVertConnect = diffvertconnect_lcc,
                                       diffEdgeConnect = diffedgconnect_lcc,
                                       diffNatConnect = diffnatconnect_lcc,
                                       diffPEP = diffpep_lcc,
                                       diffClustCoef = diffclustcoef_lcc,
                                       diffModul = diffmod_lcc),
                 diffsCentr = list(diffDeg = diffdeg, 
                                   diffBetw = diffbetw,
                                   diffClose = diffclose, 
                                   diffEigen = diffeigen),
                 absDiffsCentr = list(absDiffDeg = absdiffdeg, 
                                      absDiffBetw = absdiffbetw,
                                      absDiffClose = absdiffclose, 
                                      absDiffEigen = absdiffeigen),
                 props = list(deg1 = props1$deg, 
                              deg2 = props2$deg,
                              betw1 = props1$betw, 
                              betw2 = props2$betw,
                              close1 = props1$close, 
                              close2 = props2$close,
                              eigen1 = props1$eigen, 
                              eigen2 = props2$eigen,
                              hubs1 = props1$hubs, 
                              hubs2 = props2$hubs,
                              nComp1 = props1$nComp,
                              nComp2 = props2$nComp,
                              avDiss1 = avDiss1,
                              avDiss2 = avDiss2,
                              avPath1 = avPath1, 
                              avPath2 = avPath2,
                              density1 = props1$density, 
                              density2 = props2$density,
                              vertConnect1 = props1$vertconnect, 
                              vertConnect2 = props2$vertconnect,
                              edgeConnect1 = props1$edgeconnect, 
                              edgeConnect2 = props2$edgeconnect,
                              natConnect1 = props1$natConnect,
                              natConnect2 = props2$natConnect,
                              pep1 = props1$pep,
                              pep2 = props2$pep,
                              clustCoef1 = clustCoef1, 
                              clustCoef2 = clustCoef2,
                              modularity1 = modul1, 
                              modularity2 = modul2,
                              clust1 = clust1, 
                              clust2 = clust2),
                 propsLCC = list(lccSize1 = props1$lccSize,
                                 lccSize2 = props2$lccSize,
                                 lccSizeRel1 = props1$lccSizeRel,
                                 lccSizeRel2 = props2$lccSizeRel,
                                 avDiss1 = avDiss1_lcc,
                                 avDiss2 = avDiss2_lcc,
                                 avPath1 = avPath1_lcc, 
                                 avPath2 = avPath2_lcc,
                                 density1 = props1$density_lcc, 
                                 density2 = props2$density_lcc,
                                 vertConnect1 = props1$vertconnect_lcc, 
                                 vertConnect2 = props2$vertconnect_lcc,
                                 edgeConnect1 = props1$edgeconnect_lcc, 
                                 edgeConnect2 = props2$edgeconnect_lcc,
                                 natConnect1 = props1$natConnect_lcc,
                                 natConnect2 = props2$natConnect_lcc,
                                 pep1 = props1$pep_lcc,
                                 pep2 = props2$pep_lcc,
                                 clustCoef1 = clustCoef1_lcc, 
                                 clustCoef2 = clustCoef2_lcc,
                                 modularity1 = modul1_lcc, 
                                 modularity2 = modul2_lcc)
  )
  
  return(output)
}