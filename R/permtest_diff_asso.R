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
#' @param countsJoint joint count matrices before preprocessing
#' @param normCounts1,normCounts2 normalized count matrices.
#' @param assoMat1,assoMat2 association matrices corresponding to the two count
#'   matrices. The associations must have been estimated from the count matrices
#'   \code{countMat1} and \code{countMat2}.
#' @param paramsNetConstruct parameters used for network construction.
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
#' @param matchDesign Numeric vector with two elements specifying an optional 
#'   matched-group (i.e. matched-pair) design, which is used for the permutation 
#'   tests in \code{\link{netCompare}} and \code{\link{diffnet}}. \code{c(1,1)} 
#'   corresponds to a matched-pair design. A 1:2 matching, for instance, is 
#'   defined by \code{c(1,2)}, which means that the first sample of group 1 is 
#'   matched to the first two samples of group 2 and so on. 
#'   The appropriate order of samples must be ensured. If 
#'   \code{NULL}, the group memberships are shuffled randomly while group sizes
#'   identical to the original data set are ensured.
#' @param callNetConstr call inherited from \code{netConstruct()}.
#' @param cores number of CPU cores (permutation tests are executed parallel)
#' @param verbose if \code{TRUE}, status messages and numbers of SparCC
#'   iterations are printed
#' @param logFile character string naming the log file within which the current
#'   iteration number is stored
#' @param seed an optional seed for reproducibility of the results
#' @param fileLoadAssoPerm  character giving the name (without extenstion) 
#'   or path of the file storing the "permuted" association/dissimilarity 
#'   matrices that have been exported by setting \code{storeAssoPerm} to 
#'   \code{TRUE}. Only used for permutation tests. Set to \code{NULL} if no 
#'   existing associations should be used.
#' @param fileLoadCountsPerm character giving the name (without extenstion) 
#'   or path of the file storing the "permuted" count matrices that have been 
#'   exported by setting \code{storeCountsPerm} to \code{TRUE}. 
#'   Only used for permutation tests, and if \code{fileLoadAssoPerm = NULL}. 
#'   Set to \code{NULL} if no existing count matrices should be used.
#' @param storeAssoPerm logical indicating whether the association (or
#'   dissimilarity) matrices for the permuted data should be stored in a file.
#'   The filename is given via \code{fileStoreAssoPerm}. If \code{TRUE}, 
#'   the computed "permutation" association/dissimilarity matrices can be reused
#'   via \code{fileLoadAssoPerm} to save runtime. Defaults to \code{FALSE}.
#' @param fileStoreAssoPerm character giving the file name to store a matrix
#'   containing a matrix with associations/dissimilarities for the permuted 
#'   data. Can also be a path.
#' @param storeCountsPerm logical indicating whether the permuted count matrices
#'   should be stored in an external file. Defaults to \code{FALSE}.
#' @param fileStoreCountsPerm character vector with two elements giving the 
#'   names of two files storing the permuted count matrices belonging to the 
#'   two groups.
#' @param assoPerm not used anymore.
#' @references{\insertRef{gill2010statistical}{NetCoMi}}
#' @references{\insertRef{knijnenburg2009fewer}{NetCoMi}}
#'
#' @seealso \code{\link{diffnet}}

permtest_diff_asso <- function(countMat1, countMat2, countsJoint,
                               normCounts1, normCounts2,
                               assoMat1, assoMat2, 
                               paramsNetConstruct,
                               method = c("connect.pairs", 
                                          "connect.variables",
                                          "connect.network"), 
                               fisherTrans = TRUE,
                               pvalsMethod = "pseudo",
                               adjust = "lfdr", 
                               adjust2 = "holm",
                               trueNullMethod = "convest", 
                               alpha = 0.05,
                               lfdrThresh = 0.2, 
                               nPerm = 1000,
                               matchDesign = NULL, 
                               callNetConstr = NULL,
                               cores = 4,
                               verbose = TRUE, 
                               logFile = "log.txt",
                               seed = NULL, 
                               fileLoadAssoPerm  = NULL,
                               fileLoadCountsPerm = NULL,
                               storeAssoPerm = FALSE,
                               fileStoreAssoPerm = "assoPerm",
                               storeCountsPerm = FALSE,
                               fileStoreCountsPerm = c("countsPerm1", 
                                                       "countsPerm2"),
                               assoPerm = NULL){

  measure <- paramsNetConstruct$measure
  measurePar <- paramsNetConstruct$measurePar
  zeroMethod <- paramsNetConstruct$zeroMethod
  zeroPar <- paramsNetConstruct$zeroPar
  needfrac <- paramsNetConstruct$needfrac
  needint <- paramsNetConstruct$needint
  normMethod <- paramsNetConstruct$normMethod
  normPar <- paramsNetConstruct$normPar
  jointPrepro <- paramsNetConstruct$jointPrepro
  
  
  if(!is.null(seed)){
    set.seed(seed)
  }

  if(jointPrepro){
    n1 <- nrow(normCounts1)
    n2 <- nrow(normCounts2)
    
    xbind <- rbind(normCounts1, normCounts2)
    
  } else{
    n1 <- nrow(countMat1)
    n2 <- nrow(countMat2)
    
    xbind <- rbind(countMat1, countMat2)
  }
  
  n <- n1 + n2
  nVar <- ncol(assoMat1)

  #_____________________________________________________

  if(!is.null(fileLoadAssoPerm)){
    stopifnot(is.character(fileLoadAssoPerm))
    stopifnot(length(fileLoadAssoPerm) == 1)
    
    fmat <- fm.open(filenamebase = fileLoadAssoPerm, 
                    readonly = TRUE)
    
    if(!(dim(fmat)[1] == nVar * nPerm && dim(fmat)[2] == nVar * 2)){
      stop("fileLoadAssoPerm has wrong dimensions. 
  Maybe 'nPerm' is not correct (must equal the number of permutation matrices in 'fileLoadAssoPerm').")
    }
    
    close(fmat)

  } else if(!is.null(fileLoadCountsPerm)){
    stopifnot(is.character(fileLoadCountsPerm))
    stopifnot(length(fileLoadCountsPerm) == 2)
    
    fmat_counts1 <- fm.open(filenamebase = fileLoadCountsPerm[1], 
                            readonly = TRUE)
    
    fmat_counts2 <- fm.open(filenamebase = fileLoadCountsPerm[2], 
                            readonly = TRUE)
    
    if(!(nrow(fmat_counts1) == nPerm * n1 && 
       ncol(fmat_counts1) == nVar)){
      stop("fileLoadCountsPerm has wrong dimensions. 
  Maybe 'nPerm' is not correct (must equal the number of permutation matrices in 'fileLoadCountsPerm').")
    }
    
    if(!(nrow(fmat_counts2) == nPerm * n2 && 
         ncol(fmat_counts2) == nVar)){
      stop("fileLoadCountsPerm has wrong dimensions. 
  Maybe 'nPerm' is not correct (must equal the number of permutation matrices in 'fileLoadCountsPerm').")
    }
    
    close(fmat_counts1)
    close(fmat_counts2)

  } else{
    perm_group_mat <- get_perm_group_mat(n1 = n1, n2 = n2, n = n, nPerm = nPerm, 
                                         matchDesign = matchDesign)
  }
  
  #-----------------
  if(storeCountsPerm){
    stopifnot(is.character(fileStoreCountsPerm))
    stopifnot(length(fileStoreCountsPerm) == 2)
    
    fmat_counts1 <- fm.create(filenamebase = fileStoreCountsPerm[1], 
                              nrow = (n1 * nPerm), ncol = nVar)
    fmat_counts2 <- fm.create(filenamebase = fileStoreCountsPerm[2], 
                              nrow = (n2 * nPerm), ncol = nVar)
    
    if(verbose){
      message("Files '", 
              paste0(fileStoreCountsPerm[1], ".bmat, "),
              paste0(fileStoreCountsPerm[1], ".desc.txt, "),
              paste0(fileStoreCountsPerm[2], ".bmat, and "),
              paste0(fileStoreCountsPerm[2], ".desc.txt created. "))
    }
    
    close(fmat_counts1)
    close(fmat_counts2)
  }
  
  #-----------------
  if(storeAssoPerm){
    stopifnot(is.character(fileStoreAssoPerm))
    stopifnot(length(fileStoreAssoPerm) == 1)
    
    fmat = fm.create(filenamebase = fileStoreAssoPerm, 
                     nrow = (nVar * nPerm), ncol = (2 * nVar))
    
    if(verbose){
      message("Files '", 
              paste0(fileStoreAssoPerm, ".bmat and "),
              paste0(fileStoreAssoPerm, ".desc.txt created. "))
    }
    
    close(fmat)
  }
  
  #_____________________________________________________
  # generate teststatistics for permutated data
  
  if(!is.null(seed)){
    seeds <- sample.int(1e8, size = nPerm)
  }
  
  if(verbose) message("Execute permutation tests ... ")

  if(cores > 1){
    cl <- makeCluster(cores, outfile = "")
    doSNOW::registerDoSNOW(cl)
    if(!is.null(logFile)) cat("", file=logFile, append=FALSE)
    '%do_or_dopar%' <- get('%dopar%')
  } else{
    '%do_or_dopar%' <- get('%do%')
  }

  if(verbose){
    pb <- utils::txtProgressBar(0, nPerm, style=3)
    
    progress <- function(n){
      utils::setTxtProgressBar(pb,n)
    }
    
    opts <- list(progress=progress)
  } else{
    opts <- list()
  }

  p <- NULL
  
  result <- foreach(p = 1:nPerm,
                    .packages = c("filematrix"),
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
                      
                      if(!is.null(assoPerm)){
                        # load permutation association matrices (old version)
                        assoMat1.tmp <- assoPerm[[1]][[p]]
                        assoMat2.tmp <- assoPerm[[2]][[p]]
                        count1.tmp <- count2.tmp <- NULL
                        dimnames(assoMat1.tmp) <- dimnames(assoMat1)
                        dimnames(assoMat2.tmp) <- dimnames(assoMat2)
                        
                      } else if(!is.null(fileLoadAssoPerm )){
                        
                        fmat <- fm.open(filenamebase = fileLoadAssoPerm, 
                                        readonly = TRUE)
                        
                        # load permutation asso/diss matrices
                        assoMat1.tmp <- fmat[(p-1) * nVar + (1:nVar), 
                                             1:nVar]
                        assoMat2.tmp <- fmat[(p-1) * nVar + (1:nVar), 
                                             nVar + (1:nVar)]
                        
                        dimnames(assoMat1.tmp) <- dimnames(assoMat1)
                        dimnames(assoMat2.tmp) <- dimnames(assoMat2)
                        
                        count1.tmp <- count2.tmp <- NULL
                        
                        close(fmat)
                        
                      } else{
                        
                        if(!is.null(fileLoadCountsPerm)){
                          
                          fmat_counts1 <- fm.open(filenamebase = 
                                                    fileLoadCountsPerm[1], 
                                                  readonly = TRUE)
                          
                          fmat_counts2 <- fm.open(filenamebase = 
                                                    fileLoadCountsPerm[2], 
                                                  readonly = TRUE)
                          
                          # load permutation count matrices
                          count1.tmp <- 
                            fmat_counts1[(p-1) * n1 + (1:n1), 1:nVar]
                          
                          count2.tmp <- 
                            fmat_counts2[(p-1) * n2 + (1:n2), 1:nVar]
                          
                          close(fmat_counts1)
                          close(fmat_counts2)
                          
                        } else{
                          # generate permutation count matrices
                          
                          count1.tmp <- 
                            xbind[which(perm_group_mat[p, ] == 1), ]
                          
                          count2.tmp <- 
                            xbind[which(perm_group_mat[p, ] == 2), ]
                          
                          if(!jointPrepro){
                            # zero treatment and normalization necessary if
                            # in network construction two count matrices 
                            # were given or dissimilarity network is created
                            
                            suppressMessages(
                              count1.tmp <- zero_treat(countMat = count1.tmp, 
                                                       zeroMethod = zeroMethod,
                                                       zeroParam = zeroPar, 
                                                       needfrac = needfrac,
                                                       needint = needint, 
                                                       verbose = FALSE))
                            
                            suppressMessages(
                              count2.tmp <- zero_treat(countMat = count2.tmp, 
                                                       zeroMethod = zeroMethod,
                                                       zeroParam = zeroPar, 
                                                       needfrac = needfrac,
                                                       needint = needint, 
                                                       verbose = FALSE))
                            
                            suppressMessages(
                              count1.tmp <- norm_counts(countMat = count1.tmp, 
                                                        normMethod = normMethod,
                                                        normParam = normPar, 
                                                        zeroMethod = zeroMethod,
                                                        needfrac = needfrac, 
                                                        verbose = FALSE))
                            
                            suppressMessages(
                              count2.tmp <- norm_counts(countMat = count2.tmp, 
                                                        normMethod = normMethod,
                                                        normParam = normPar, 
                                                        zeroMethod = zeroMethod,
                                                        needfrac = needfrac, 
                                                        verbose = FALSE))
                          }
                          
                          if(storeCountsPerm){
                            fmat_counts1 <- fm.open(filenamebase = 
                                                      fileStoreCountsPerm[1])
                            
                            fmat_counts2 <- fm.open(filenamebase = 
                                                      fileStoreCountsPerm[2])
                            
                            fmat_counts1[(p-1) * n1 + (1:n1), 
                                         1:nVar] <- count1.tmp
                            fmat_counts2[(p-1) * n2 + (1:n2), 
                                         1:nVar] <- count2.tmp
                            
                            close(fmat_counts1)
                            close(fmat_counts2)
                          }
                        }
                        
                        assoMat1.tmp <- calc_association(count1.tmp,
                                                         measure = measure,
                                                         measurePar = measurePar,
                                                         verbose = FALSE)
                        
                        assoMat2.tmp <- calc_association(count2.tmp,
                                                         measure = measure,
                                                         measurePar = measurePar,
                                                         verbose = FALSE)
                        
                        dimnames(assoMat1.tmp) <- dimnames(assoMat1)
                        dimnames(assoMat2.tmp) <- dimnames(assoMat2)
                        
                        if(storeAssoPerm){
                          fmat <- fm.open(filenamebase = fileStoreAssoPerm)
                          
                          fmat[(p-1) * nVar + (1:nVar), 
                               1:nVar] <- assoMat1.tmp
                          fmat[(p-1) * nVar + (1:nVar), 
                               nVar + (1:nVar)] <- assoMat2.tmp
                          
                          close(fmat)
                        }
                      }

                      returnlist <- list()

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
                                                                   nVar,
                                                                   fisherTrans)
                        returnlist[["connectVariables"]] <- connectVariables
                      }
                      if("connect.network" %in% method){
                        connectNetwork <- diff_connect_network(assoMat1.tmp,
                                                               assoMat2.tmp,
                                                               nVar,
                                                               fisherTrans)
                        returnlist[["connectNetwork"]] <- connectNetwork
                      }

                      returnlist
                    }
  
  if(cores > 1) stopCluster(cl)
  
  #____________________________________________________

  output <- list()
  
  if(verbose) close(pb) 

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
    
    output[["testStatData"]] <- connectPairsOrig
    
    output[["testStatPerm"]] <- connectPairs

  }

  if("connect.variables" %in% method){

    connectVariablesOrig <- diff_connect_variables(assoMat1, assoMat2, nVar, 
                                                   fisherTrans)

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
      lfdr <- fdrtool::fdrtool(pvalsConnectVariables, statistic = "pvalue", 
                               plot = FALSE)$lfdr
      output[["pAdjustConnectVariables"]] <- lfdr
      
    } else{
      p.adj <- p.adjust(pvalsConnectVariables, adjust2)
      output[["pAdjustConnectVariables"]] <- p.adj
    }

  }

  if("connect.network" %in% method){
    connectNetworkOrig <- diff_connect_network(assoMat1, assoMat2, nVar, 
                                               fisherTrans)

    connectNetwork <- numeric(nPerm)

    for(i in 1:nPerm){
      connectNetwork[i] <- result[[i]]$connectNetwork
    }

    pvalConnectNetwork <- sum(connectNetwork >= connectNetworkOrig) / nPerm
    output[["pvalConnectNetwork"]] <- pvalConnectNetwork
  }

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
diff_connect_variables <- function(assoMat1, assoMat2, nVar, fisherTrans = TRUE){

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
  d_vec <- numeric(nVar)

  # calculate test statistics for all variables
  for(i in 1:nVar){
    d_vec[i] <- sum(diff[i, -i]) / (nVar - 1)
  }

  names(d_vec) <- colnames(assoMat1)
  return(d_vec)
}


# calculate test statistic for difference in connectivity for the whole network
diff_connect_network <- function(assoMat1, assoMat2, nVar, fisherTrans = TRUE){

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
  d <- (sum(diff) - sum(diag(diff))) / (nVar * (nVar-1))
  return(d)
}


