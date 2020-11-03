sparsify <- function(assoMat, countMat, sampleSize, measure, measurePar, assoType,
                     sparsMethod, thresh, alpha, adjust, lfdrThresh,
                     trueNullMethod, nboot, softThreshType, softThreshPower,
                     softThreshCut, kNeighbor, knnMutual, parallel, cores,
                     verbose, logFile = NULL, seed = NULL){

  if(sparsMethod == "threshold"){

    assoNew <- assoMat

    if(assoType == "dissimilarity"){
      if(thresh < 0) warning("'thresh' smaller than 0. No edges selected.")
      assoNew[assoNew > thresh] <- Inf

    } else{
      if(thresh > 1) warning("'thresh' not in [0,1]. No edges selected.")
      if(thresh < 0) warning("'thresh' not in [0,1]. All edges selected.")
      assoNew[abs(assoNew) < thresh] <- 0
    }

  } else if(sparsMethod %in% c("t-test", "bootstrap")){   #hypothesis testing

    assoVec <- assoMat[lower.tri(assoMat)]

    if(sparsMethod == "t-test"){

      if(is.null(countMat)){
        df <- sampleSize
      } else{
        df <- nrow(countMat) - 2
      }

      tstat <- (assoVec * sqrt(df)) / sqrt(1 - assoVec^2)
      pvals <- 2 * (1 - stats::pt(abs(tstat), df))
      pvals[pvals] <- 1

    } else{   #bootstrap procedure

      if(is.null(countMat)){
        stop("Count matrix needed for sparsification via bootstrapping.")
      }
      pvalsBoot <- boottest(countMat = countMat, assoMat = assoMat,
                            nboot = nboot, measure = measure,
                            measurePar = measurePar, parallel = parallel,
                            cores = cores, logFile = logFile, verbose = verbose,
                            seed = seed)

      pvals <- pvalsBoot[lower.tri(pvalsBoot)]
    }

    # adjust for multiple testing and identify links
    if(adjust == "none"){
      assoVec[pvals > alpha] <- 0

    } else{

      if(verbose %in% 2:3){
        message("\nAdjust for multiple testing via '", adjust, "' ... ",
                appendLF = FALSE)
      }

      verb.tmp <- ifelse(verbose == 3, TRUE, FALSE)
      pvals.adj <- multAdjust(pvals = pvals, adjust = adjust,
                              trueNullMethod = trueNullMethod,
                              verbose = verb.tmp)

      if(verbose %in% 2:3){
        message("Done.")
      }

      assoVec[pvals.adj > alpha] <- 0
    }

    assoNew <- assoMat
    assoNew[lower.tri(assoNew)] <- assoVec
    assoNew[upper.tri(assoNew)] <- t(assoNew)[upper.tri(t(assoNew))]

  } else if(sparsMethod == "softThreshold"){

    if(softThreshType == "signed"){
      simMat <- 0.5 + 0.5 * assoMat
    } else if(softThreshType == "unsigned"){
      simMat <- abs(assoMat)
    } else{
      simMat <- assoMat
      simMat[simMat < 0] <- 0
    }

    # power for soft thresholding
    if(is.null(softThreshPower)){

      diag(assoMat) <- 0

      if(verbose == 3){
        message("\nExecute pickSoftThreshold() ...")
        softthresh1 <- WGCNA::pickSoftThreshold(simMat, dataIsExpr = FALSE,
                                                networkType = softThreshType,
                                                RsquaredCut = softThreshCut,
                                                verbose = 5)
      } else{
        invisible(capture.output(softthresh1 <- 
                         WGCNA::pickSoftThreshold(simMat,
                                                  dataIsExpr = FALSE,
                                                  networkType = softThreshType,
                                                  RsquaredCut = softThreshCut,
                                                  verbose = 2)))
      }

      power <-  softthresh1$powerEstimate

      if(is.na(power)){
        warning(paste0("No power with R^2 above ", softThreshCut,
                       ". Power set to 1."))
        power <- 1
      }

      if(verbose >= 2){
        message("\n Estimated power: ", power)
      }


    } else{
      power <- softThreshPower
    }

    assoNew <- simMat^power

  } else if(sparsMethod == "knn"){

    knngraph <- cccd::nng(dx = assoMat, k = kNeighbor, mutual = knnMutual)
    knnMat <- as.matrix(igraph::as_adjacency_matrix(knngraph))

    if(!knnMutual){
      knnMat.tmp <- knnMat
      for (i in 1:nrow(knnMat)) {
        knnMat[i, ] <- knnMat.tmp[i, ] | knnMat.tmp[, i]
        knnMat[, i] <- knnMat[i, ]
      }
    }

    assoNew <- assoMat
    assoNew[knnMat == 0] <- Inf
    diag(assoNew) <- 0

  } else{
    assoNew <- assoMat
  }

  #-----------------------------------
  if(sparsMethod != "softThreshold"){
    power <- NULL
    simMat <- NULL
  }

  output <- list(assoNew = assoNew, power = power, simMat = simMat)
}

