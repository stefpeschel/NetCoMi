calc_diff_props <- function(adja1, adja2, dissMat1, dissMat2, weighted,
                          clustMethod, clustPar, hubPar,
                          hubQuant, jaccQuant, lnormFit, connect,
                          weightDeg, normDeg, normBetw, normClose, normEigen,
                          testJacc = TRUE, jaccTestGreater = FALSE,
                          testRand = TRUE, nPermRand = 1000){

  isempty1 <- all(adja1[lower.tri(adja1)] == 0)
  isempty2 <- all(adja2[lower.tri(adja2)] == 0)

  if(isempty1 & isempty2) stop("There are no connected nodes in both networks.")

  # calculate network properties
  props1 <- calc_props(adjaMat = adja1, dissMat = dissMat1, weighted = weighted,
                      isempty = isempty1, clustMethod = clustMethod,
                      clustPar = clustPar, hubPar = hubPar,
                      hubQuant = hubQuant, lnormFit = lnormFit,
                      connect = connect,
                      weightDeg = weightDeg, normDeg = normDeg,
                      normBetw = normBetw, normClose = normClose,
                      normEigen = normEigen, jaccard = TRUE, jaccQuant = jaccQuant)

  props2 <- calc_props(adjaMat = adja2, dissMat = dissMat2, weighted = weighted,
                      isempty = isempty2, clustMethod = clustMethod,
                      clustPar = clustPar, hubPar = hubPar,
                      hubQuant = hubQuant, lnormFit = lnormFit,
                      connect = connect,
                      weightDeg = weightDeg, normDeg = normDeg,
                      normBetw = normBetw, normClose = normClose,
                      normEigen = normEigen, jaccard = TRUE, jaccQuant = jaccQuant)



  # differences in network properties
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
  #### Average Path Length (or Mean Distance)

  avPath1 <- props1$avPath
  avPath2 <- props2$avPath

  if(is.na(avPath1)) avPath1 <- 0
  if(is.na(avPath2)) avPath2 <- 0

  diffpath <- abs(avPath1 - avPath2)

  #--------------------------------------------------------------------------
  ##### Clustering

  clustCoef1 <- props1$clustCoef
  clustCoef2 <- props2$clustCoef

  if(is.na(clustCoef1)) clustCoef1 <- 0
  if(is.na(clustCoef2)) clustCoef2 <- 0
  diffclust <- abs(clustCoef1 - clustCoef2)

  #--------------------------------------------------------------------------
  # differential Modularity
  if(clustMethod != "none"){
    modul1 <- props1$modul
    modul2 <- props2$modul
    diffmod <- abs(modul1 - modul2)
  } else{
    modul1 <- modul2 <- NULL
    diffmod <- NULL
  }

  #--------------------------------------------------------------------------
  # vertex connectivity
  diffvertconnect <- abs(props1$vertconnect - props2$vertconnect)

  # edge connectivity
  diffedgconnect <- abs(props1$edgeconnect - props2$edgeconnect)

  # relative number of edges(density)
  diffdensity <- abs(props1$density - props2$density)

  #--------------------------------------------------------------------------
  # adjusted Rand index for measuring similarity between two clusterings
  clust1 <- props1$clust
  clust2 <- props2$clust
  randInd <- c(value = randIndex(table(clust1, clust2), adjust = TRUE), pval = NA)

  # significance test for Rand index
  if(testRand){
    randPerm <- numeric(nPermRand)
    for(i in 1:nPermRand){
      clust1.tmp <- gtools::permute(clust1)
      clust2.tmp <- gtools::permute(clust2)
      randPerm[i] <- randIndex(table(clust1.tmp, clust2.tmp), adjust = TRUE)
    }
    randMean <- mean(randPerm)
    randSD <- sd(randPerm)
    normRandPerm <- (randPerm - randMean) / randSD
    normRand <- (randInd[1] - randMean) / randSD
    randInd["pval"] <- (sum(normRandPerm >= abs(normRand)) + sum(normRandPerm <= -abs(normRand))) / nPermRand
  }

  #--------------------------------------------------------------------------
  # Jaccard Index

  jaccDeg <- calc_jaccard(props1$topdeg, props2$topdeg, sigTest = testJacc)
  jaccBetw <- calc_jaccard(props1$topbetw, props2$topbetw, sigTest = testJacc)
  jaccClose <- calc_jaccard(props1$topclose, props2$topclose, sigTest = testJacc)
  jaccEigen <- calc_jaccard(props1$topeigen, props2$topeigen, sigTest = testJacc)
  jaccHub <- calc_jaccard(props1$hubs, props2$hubs, sigTest = testJacc)

  return(list(jaccDeg = jaccDeg, jaccBetw =jaccBetw,
              jaccClose = jaccClose,  jaccEigen = jaccEigen,
              jaccHub = jaccHub, diffPath = diffpath,
              diffClust = diffclust, diffModul = diffmod,
              diffVertConnect = diffvertconnect, diffEdgeConnect = diffedgconnect,
              diffDensity = diffdensity,  randInd = randInd,
              props = list(deg1 = props1$deg, deg2 = props2$deg,
                           betw1 = props1$betw, betw2 = props2$betw,
                           close1 = props1$close, close2 = props2$close,
                           eigen1 = props1$eigen, eigen2 = props2$eigen,
                           hubs1 = props1$hubs, hubs2 = props2$hubs2,
                           avPath1 = avPath1, avPath2 = avPath2,
                           clustCoef1 = clustCoef1, clustCoef2 = clustCoef2,
                           modul1 = modul1, modul2 = modul2,
                           clust1 = clust1, clust2 = clust2,
                           vertConnect1 = props1$vertconnect, vertConnect2 = props2$vertconnect,
                           edgeConnect1 = props1$edgeconnect, edgeConnect2 = props2$edgeconnect,
                           density1 = props1$density, density2 = props2$density),
              diffs = list(diffDeg = diffdeg, diffBetw = diffbetw,
                           diffClose = diffclose, diffEigen = diffeigen),
              absDiffs = list(absDiffDeg = absdiffdeg, absDiffBetw = absdiffbetw,
                              absDiffClose = absdiffclose, absDiffEigen = absdiffeigen)
  )
  )
}
