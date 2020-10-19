#' @title Permutation Tests for Determining Differential Associations
#'
#' @description The function implements procedures to test whether pairs of taxa
#'   are differentially associated, whether a taxon is differentially associated
#'   to all other taxa, or whether two networks are differentially associated
#'   between two groups as proposed  by \cite{Gill et al.(2010)}.
#'
#' @param countMat1,countMat2 matrices containing microbiome data (read counts)
#'   of group 1 and group 2 (rows represent samples and columns taxonomic units,
#'   respectively).
#' @param assoMat1,assoMat2 association matrices corresponding to the two count
#'   matrices. The associations must have been estimated from the count matrices
#'   \code{countMat1} and \code{countMat2}.
#' @param measure a character indicating the association measure that has been
#'   used for determining the association matrices.
#' @param measurePar list with parameters passed to the function for association
#'   calculation.
#' @param method character vector indicating the tests to be performed. Possible
#'   values are \code{"connect.pairs"} (differentially correlated taxa pairs),
#'   \code{"connect.variables"} (one taxon to all other) and
#'   \code{"connect.network"} (differentially connected networks). By default,
#'   all three tests are conducted.
#' @param fisherTrans logical indicating whether the correlation values should
#'   be Fisher-transformed.
#' @param pvalsMethod currently only \code{"pseudo"} is available, where 1 is
#'   added to the number of permutations and the permutation test statistics
#'   being more extreme than the observed one in order to avoid zero p-values.
#' @param adjust multiple testing adjustment for the tests for differentially
#'   correlated pairs of taxa; possible values are "lfdr" (default) for local
#'   false discovery rate correction (via \code{\link[fdrtool]{fdrtool}}) or one
#'   of the methods provided by \code{\link[stats]{p.adjust}}
#' @param adjust2 multiple testing adjustment for the tests if a taxa pair is
#'   differentially correlated to all other taxa; possible methods are those
#'   provided by \code{\link[stats]{p.adjust}} (a few hundred tests are
#'   necessary for the local fdr correction)
#' @param trueNullMethod character indicating the method used for estimating the
#'   proportion of true null hypotheses from a vector of p-values. Used for the
#'   adaptive Benjamini-Hochberg method for multiple testing adjustment (chosen
#'   by \code{adjust = "adaptBH"}).
#' @param alpha significance level
#' @param lfdrThresh defines a threshold for the local fdr if "lfdr" is chosen
#'   as method for multiple testing correction; defaults to 0.2, which means
#'   that correlations with a corresponding local fdr less than or equal to 0.2
#'   are identified as significant
#' @param nPerm number of permutations
#' @param cores number of CPU cores (permutation tests are executed parallel)
#' @param verbose if \code{TRUE}, status messages and numbers of SparCC
#'   iterations are printed
#' @param logFile character string naming the log file within which the current
#'   iteration number is stored
#' @param seed an optional seed for reproducibility of the results
#' @param assoPerm list with two elements containing permutation association
#'   matrices for the two groups.
#' @references{\insertRef{gill2010statistical}{NetCoMi}}
#' @references{\insertRef{knijnenburg2009fewer}{NetCoMi}}
#'
#' @seealso \code{\link{diffnet}}

permtest_diff_asso <- function(countMat1, countMat2, assoMat1, assoMat2,
                           measure, measurePar,
                           method = c("connect.pairs", "connect.variables",
                                      "connect.network"), fisherTrans = TRUE,
                           pvalsMethod = "pseudo",
                           adjust = "lfdr", adjust2 = "holm",
                           trueNullMethod = "convest", alpha = 0.05,
                           lfdrThresh = 0.2, nPerm = 1000, cores = 4,
                           verbose = TRUE, logFile = "log.txt",
                           seed = NULL, assoPerm = NULL){

  if(!is.null(seed)){
    set.seed(seed)
    seeds <- sample.int(1e8, size = nPerm)
  }

  nVars <- ncol(assoMat1)
  n1 <- nrow(countMat1)
  n2 <- nrow(countMat2)
  n <- n1 + n2
  counts <- rbind(countMat1, countMat2)
  p <- ncol(assoMat1)

  #_____________________________________________________
  # generate teststatistics for permutated data
  if(verbose) message("Execute permutation tests ... ")

  if(cores > 1){
    cl <- makeCluster(cores, outfile = "")
    registerDoSNOW(cl)
    if(!is.null(logFile)) cat("", file=logFile, append=FALSE)
    '%do_or_dopar%' <- get('%dopar%')
  } else{
    '%do_or_dopar%' <- get('%do%')
  }

  if(verbose){
    progress <- function(p){
      progr <- round(p/nPerm*100)
      message(progr,"%\r",appendLF=FALSE)
    }
    message("0%\r",appendLF=FALSE)
    opts <- list(progress=progress)
  } else{
    opts <- list()
  }

  result <- foreach(p = 1:nPerm,
                    .packages = c("SPRING", "fdrtool",
                                               "SpiecEasi",
                                               "ccrepe",
                                               "vegan",
                                               "LaplacesDemon",
                                               "robCompositions",
                                               "propr",
                                               "zCompositions",
                                               "DESeq2",
                                               "compositions"),
                    .export = c("calc_association", "cclasso", "gcoda",
                                "diff_connect_pairs", "diff_connect_variables",
                                "diff_connect_network", "get_vec_names"),
                    .options.snow = opts) %do_or_dopar% {

                      if(!is.null(logFile)){
                        cat(paste("Iteration", p,"\n"),
                            file=logFile, append=TRUE)
                      }
                      if(verbose) progress(p)

                      if(!is.null(seed)) set.seed(seeds[p])

                      index <- sample(1:n, n)
                      countMat1.tmp <- counts[index[1:n1], ]
                      countMat2.tmp <- counts[index[(n1+1):n], ]

                      if(!is.null(assoPerm)){
                        assoMat1.tmp <- assoPerm[[1]][[p]]
                        assoMat2.tmp <- assoPerm[[2]][[p]]
                      } else{

                        assoMat1.tmp <- calc_association(countMat1.tmp,
                                                         measure = measure,
                                                         measurePar = measurePar,
                                                         verbose = FALSE)
                        assoMat2.tmp <- calc_association(countMat2.tmp,
                                                         measure = measure,
                                                         measurePar = measurePar,
                                                         verbose = FALSE)
                      }

                      returnlist <- list()

                      returnlist[["assoMat1"]] <- assoMat1.tmp
                      returnlist[["assoMat2"]] <- assoMat2.tmp

                      # teststatistics for simulated data
                      if("connect.pairs" %in% method){
                        connectPairs <- diff_connect_pairs(assoMat1.tmp,
                                                           assoMat2.tmp,
                                                           fisherTrans)
                        returnlist[["connectPairs"]] <- connectPairs
                      }
                      if("connect.variables" %in% method){
                        connectVariables <- diff_connect_variables(assoMat1.tmp,
                                                                   assoMat2.tmp,
                                                                   nVars,
                                                                   fisherTrans)
                        returnlist[["connectVariables"]] <- connectVariables
                      }
                      if("connect.network" %in% method){
                        connectNetwork <- diff_connect_network(assoMat1.tmp,
                                                               assoMat2.tmp,
                                                               nVars,
                                                               fisherTrans)
                        returnlist[["connectNetwork"]] <- connectNetwork
                      }

                      returnlist
                    }
  if(cores > 1) stopCluster(cl)
  #____________________________________________________

  output <- list()
  if(verbose) message("")

  if("connect.pairs" %in% method){
    # test statistic for original data
    connectPairsOrig <- diff_connect_pairs(assoMat1, assoMat2, fisherTrans)

    # test statistics for simulated data
    connectPairs <- matrix(NA, nPerm, length(connectPairsOrig),
                            dimnames = list(1:nPerm, names(connectPairsOrig)))
    for(i in 1:nPerm){
      connectPairs[i, ] <- result[[i]]$connectPairs
    }

    if(pvalsMethod == "pseudo"){

      pvalsVec <- sapply(1:ncol(connectPairs), function(i){
        (sum(connectPairs[, i] >= connectPairsOrig[i]) + 1) / (nPerm + 1)
      })

    } #else{
    #   pvalsVec <- calc_perm_pvals(tstatPerm = connectPairs,
    #                               tstatObs = connectPairsOrig,
    #                               nExceed = nExceed, ADalpha = 0.05)
    # }


    names(pvalsVec) <- names(connectPairsOrig)
    output[["pvalsVec"]] <- pvalsVec

    # adjust for multiple testing
    if(verbose & adjust != "none"){
      message("Adjust for multiple testing using '", adjust, "' ... ",
              appendLF = FALSE)
    }
    pAdjust <- multAdjust(pvals = pvalsVec, adjust = adjust,
                          trueNullMethod = trueNullMethod, verbose = verbose)
    if(verbose & adjust != "none") message("Done.")
    output[["pAdjustVec"]] <- pAdjust

    mat.tmp <- assoMat1
    mat.tmp[lower.tri(mat.tmp)] <- pvalsVec
    mat.tmp[upper.tri(mat.tmp)] <- t(mat.tmp)[upper.tri(t(mat.tmp))]
    pvalsMat <- mat.tmp
    output[["pvalsMat"]] <- pvalsMat

    if(adjust == "none"){
      output[["pAdjustMat"]] <- NULL
    } else{
      mat.tmp[lower.tri(mat.tmp)] <- output[["pAdjustVec"]]
      mat.tmp[upper.tri(mat.tmp)] <- t(mat.tmp)[upper.tri(t(mat.tmp))]
      output[["pAdjustMat"]] <- mat.tmp
    }

  }

  if("connect.variables" %in% method){

    connectVariablesOrig <- diff_connect_variables(assoMat1, assoMat2, nVars, fisherTrans)

    connectVariables <- matrix(NA, nPerm, length(connectVariablesOrig),
                                dimnames = list(1:nPerm,
                                                names(connectVariablesOrig)))
    for(i in 1:nPerm){
      connectVariables[i, ] <- result[[i]]$connectVariables
    }

    pvalsConnectVariables <- sapply(1:ncol(connectVariables), function(i){
      (sum(connectVariables[, i] >= connectVariablesOrig[i]) + 1) / (nPerm + 1)
    })

    names(pvalsConnectVariables) <- names(connectVariablesOrig)
    output[["pvalsConnectVariables"]] <- pvalsConnectVariables

    # adjust for multiple testing
    if(adjust2 == "none"){
      output[["pAdjustConnectVariables"]] <- NULL
    } else if(adjust == "lfdr"){
      lfdr <- fdrtool(pvalsConnectVariables, statistic = "pvalue", plot = FALSE)$lfdr
      output[["pAdjustConnectVariables"]] <- lfdr
    } else{
      p.adj <- p.adjust(pvalsConnectVariables, adjust2)
      output[["pAdjustConnectVariables"]] <- p.adj
    }
  }

  if("connect.network" %in% method){
    connectNetworkOrig <- diff_connect_network(assoMat1, assoMat2, nVars, fisherTrans)

    connectNetwork <- numeric(nPerm)

    for(i in 1:nPerm){
      connectNetwork[i] <- result[[i]]$connectNetwork
    }

    pvalConnectNetwork <- sum(connectNetwork >= connectNetworkOrig) / nPerm
    output[["pvalConnectNetwork"]] <- pvalConnectNetwork
  }


  assoPerm1 <- assoPerm2 <- list()
  for(i in 1:nPerm){
    assoPerm1[[i]] <- result[[i]]$assoMat1
    assoPerm2[[i]] <- result[[i]]$assoMat2
  }

  output[["assoPerm1"]] <- assoPerm1
  output[["assoPerm2"]] <- assoPerm2

  return(output)
}



# calculate test statistic for differential connectivity for all variable pairs
diff_connect_pairs <- function(assoMat1, assoMat2, fisherTrans = TRUE){

  # transform distance matrix to vector
  diag <- lower.tri(assoMat1, diag = FALSE)
  distvec1 <- assoMat1[diag]
  distvec2 <- assoMat2[diag]
  names(distvec1) <- names(distvec2) <- get_vec_names(assoMat1)

  if(fisherTrans){
    # Fisher transformation of correlation coefficients
    z1 <- atanh(distvec1)
    z2 <- atanh(distvec2)
    diff <- abs(z1 - z2)
  } else{
    diff <- abs(distvec1 - distvec2)
  }
  return(diff)
}


# calculate test statistic for difference in connectivity between a single
# variable and all other variables
diff_connect_variables <- function(assoMat1, assoMat2, nVars, fisherTrans = TRUE){

  if(fisherTrans){
    # Fisher transformation of correlation coefficients
    Z1 <- atanh(assoMat1)
    Z2 <- atanh(assoMat2)

    # build matrix with absolute differences of z-values
    diff <- abs(Z1 - Z2)
    diag(diff) <- 0
  } else{
    diff <- abs(assoMat1 - assoMat2)
  }

  ### test statistic
  d_vec <- numeric(nVars)

  # calculate test statistics for all variables
  for(i in 1:nVars){
    d_vec[i] <- sum(diff[i, -i]) / (nVars - 1)
  }
  names(d_vec) <- colnames(assoMat1)
  return(d_vec)
}


# calculate test statistic for difference in connectivity for the whole network
diff_connect_network <- function(assoMat1, assoMat2, nVars, fisherTrans = TRUE){

  if(fisherTrans){
    # Fisher transformation of correlation coefficients
    Z1 <- atanh(assoMat1)
    Z2 <- atanh(assoMat2)

    # build matrix with absolute differences of z-values
    diff <- abs(Z1 - Z2)
    diag(diff) <- 0
  } else{
    diff <- abs(assoMat1 - assoMat2)
  }

  # test statistic
  d <- (sum(diff) - sum(diag(diff))) / (nVars * (nVars-1))
  return(d)
}


