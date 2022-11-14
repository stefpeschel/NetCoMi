#' @title Summary Method for Objects of Class microNetComp
#' @description  The main results returned by \code{\link{netCompare}} are
#'   printed in a well-arranged format.
#'
#' @param object object of class \code{microNetComp} returned by 
#'   \code{\link{netCompare}}.
#' @param groupNames character vector with two elements giving the group names
#'   for the two networks. If \code{NULL}, the names are adopted
#'   from \code{object}.
#' @param showCentr character vector indicating which centrality measures
#'   should be included in the summary. Possible values are "all", "degree",
#'   "betweenness", "closeness", "eigenvector" and "none".
#' @param numbNodes integer indicating for how many nodes the centrality
#'   values shall be printed. Defaults to 10 meaning that the first 10 taxa
#'   with highest absolute group difference of the specific centrality
#'   measure are shown.
#' @param showGlobal logical. If \code{TRUE}, global network properties for the 
#'   whole network are printed.
#' @param showGlobalLCC logical. If \code{TRUE}, global network properties for 
#'   the largest connected component are printed. If the network is connected 
#'   (number of components is 1) the global properties are printed only once 
#'   (if one of the arguments \code{showGlobal} and \code{showGlobalLCC}) is 
#'   \code{TRUE}.
#' @param showJacc logical. If \code{TRUE}, the Jaccard index is printed.
#' @param showRand logical. If \code{TRUE}, the adjusted Rand index (if 
#'   existent) is returned. 
#' @param showGCD logical. If \code{TRUE}, the Graphlet Correlation Distance 
#'   (if existent) is printed.
#' @param pAdjust logical. The permutation p-values (if existent) are adjusted 
#'   if \code{TRUE} (default) and not adjusted if \code{FALSE}.
#' @param digits integer giving the number of decimal places to which the 
#'   results are rounded. Defaults to 3L. 
#' @param digitsPval integer giving the number of decimal places to which the 
#'   p-values are rounded. Defaults to 6L.
#' @param ... not used.
#'
#' @seealso \code{\link{netCompare}}
#'
#' @method summary microNetComp
#' @rdname summarize.microNetComp
#' @export

summary.microNetComp <- function(object, 
                                 groupNames = NULL, 
                                 showCentr = "all", 
                                 numbNodes = 10L, 
                                 showGlobal = TRUE,
                                 showGlobalLCC = TRUE,
                                 showJacc = TRUE,
                                 showRand = TRUE,
                                 showGCD = TRUE,
                                 pAdjust = TRUE,
                                 digits = 3L, 
                                 digitsPval = 6L,
                                 ...) {
  
  #-----------------------------------------------------------------------------
  # Check input arguments
  
  if (!inherits(object, "microNetComp")) {
    stop('"object" must be of class "microNetComp".')
  }
  
  showCentr <- match.arg(showCentr, choices = c("all", "none", "degree", 
                                                "betweenness", "closeness", 
                                                "eigenvector"),
                         several.ok = TRUE)
  
  if ("none" %in% showCentr) stopifnot(length(showCentr) == 1)
  
  stopifnot(is.logical(showGlobal))
  stopifnot(is.logical(showGlobalLCC))
  stopifnot(is.logical(showJacc))
  stopifnot(is.logical(showRand))
  stopifnot(is.logical(showGCD))
  stopifnot(is.logical(pAdjust))
  
  stopifnot(is.numeric(numbNodes))
  numbNodes <- as.integer(numbNodes)
  stopifnot(numbNodes >= 1)
  numbNodes <- min(numbNodes, length(object$diffCentr$diffDeg))
  
  stopifnot(is.numeric(digits))
  digits <- as.integer(digits)
  
  stopifnot(is.numeric(digitsPval))
  digitsPval <- as.integer(digitsPval)
  
  if (is.null(groupNames)) {
    group1 <- paste0("group '" , object$groups$group1, "'")
    group2 <- paste0("group '" , object$groups$group2, "'")
    
  } else {
    if (!length(groupNames) == 2) {
      stop('Length of "groupNames" must be 2.')
    }
    group1 <- groupNames[1]
    group2 <- groupNames[2]
  }
  
  #=============================================================================
  # Jaccard index
  
  if (showJacc) {
    jacc <- numeric(5)
    p.greater <- p.less <- numeric(5)
    sig.greater <- sig.less <- character(5)
    
    tmp <- c("jaccDeg", "jaccBetw", "jaccClose", "jaccEigen", "jaccHub")
    
    for (i in 1:5) {
      jacc[i] <- round(object[[tmp[i]]]["jacc"], digits)
      
      p.greater[i] <- round(object[[tmp[i]]]["p.greater"], digitsPval)
      
      p.less[i] <- round(object[[tmp[i]]]["p.less"], digitsPval)
      
      sig.greater[i] <- .getSigCode(p.greater[i])
      
      sig.less[i] <- .getSigCode(p.less[i])
      
    }
    
    jaccmat <- data.frame(jacc, p.less, sig.less, p.greater, sig.greater)
    
    rownames(jaccmat) <- c("degree", "betweenness centr.", "closeness centr.", 
                           "eigenvec. centr.", "hub taxa")
    
    colnames(jaccmat) <- c("Jacc", "  P(<=Jacc)", "   ", "P(>=Jacc)", "  ")
    
  } else {
    jaccmat <- NULL
  }
  
  #=============================================================================
  # Rand index
  
  if (showRand) {
    rand <- matrix(0, nrow = 1, ncol = 2,
                   dimnames = list("", c("ARI", "      p-value")))
    rand <- as.data.frame(rand)
    
    rand[1,1] <- round(object$randInd[1], digits)
    rand[1,2] <- round(object$randInd[2], digitsPval)
    
    if (is.na(object$randInd[2])) {
      # no p-value computed
      rand <- matrix(0, nrow = 1, ncol = 2,
                     dimnames = list("ARI", 
                                     c("wholeNet", "      LCC")))
      rand <- as.data.frame(rand)
      
      rand[1, 1] <- round(object$randInd[1], digits)
      rand[1, 2] <- round(object$randIndLCC[1], digits)
      
    } else {
      rand <- matrix(0, nrow = 2, ncol = 2,
                     dimnames = list(c("ARI", "p-value"),
                                     c("wholeNet", "      LCC")))
      rand <- as.data.frame(rand)
      
      rand[1, 1] <- round(object$randInd[1], digits)
      rand[1, 2] <- round(object$randIndLCC[1], digits)
      rand[2, 1] <- round(object$randInd[2], digitsPval)
      rand[2, 2] <- round(object$randIndLCC[2], digitsPval)
    }
  } else {
    rand <- NULL
  }
  
  #=============================================================================
  # Graphlet Correlation Distance

  if (showGCD & !is.null(object$gcd)) {
    if (is.null(object$gcd$pval)) {
      gcd <- matrix(0, nrow = 1, ncol = 2,
                    dimnames = list("GCD", 
                                    c("wholeNet", "      LCC")))
      gcd <- as.data.frame(gcd)
      
      gcd[1, 1] <- round(object$gcd$gcd, digits)
      gcd[1, 2] <- round(object$gcdLCC$gcd, digits)
      
    } else {
      gcd <- matrix(0, nrow = 2, ncol = 4,
                    dimnames = list(c("GCD", "p-value"), 
                                    c("wholeNet", " ", 
                                      "      LCC", " ")))
      gcd <- as.data.frame(gcd)
      
      gcd[1, 1] <- round(object$gcd$gcd, digits)
      gcd[1, 2] <- " "
      gcd[1, 3] <- round(object$gcdLCC$gcd, digits)
      gcd[1, 4] <- " "
      
      
      gcd[2, 1] <- round(object$gcd$pval, digitsPval)
      gcd[2, 2] <- .getSigCode(object$gcd$pval)
      gcd[2, 3] <- round(object$gcdLCC$pval, digitsPval)
      gcd[2, 4] <- .getSigCode(object$gcdLCC$pval)
      
    }



    
  } else {
    gcd <- NULL
  }
  
  #=============================================================================
  # Global network properties
  
  if (object$properties$nComp1 == 1 && object$properties$nComp1 == 1) {
    is_disconnected <- FALSE
  } else {
    is_disconnected <- TRUE
  }
  
  if (!showGlobal & !showGlobalLCC) {
    showGlob <- showGlobLCC <- FALSE
    
  } else if (!is_disconnected) {
    showGlob <- showGlobLCC <- TRUE
    
  } else {
    showGlob <- showGlobal
    showGlobLCC <- showGlobalLCC
  }
  
 
  # Check if p-values for global properties are given
  if (is.null(object$pvalDiffGlobal)) {
    showGlobPvals <- FALSE
  } else {
    showGlobPvals <- TRUE
  }
  
  #-----------------------------------------------------------------------------
  # Whole network
  
  if (showGlob) {
    glob_rnames <- c("Number of components", 
                     "Clustering coefficient",
                     "Modularity",
                     "Positive edge percentage",
                     "Edge density",
                     "Natural connectivity")
    
    glob_names <- c("nComp", "clustCoef", "modularity", "pep", "density", 
                    "natConnect")
    
    glob_names2 <- c("nComp", "ClustCoef", "Modul", "PEP", "Density", 
                     "NatConnect")
    
    if (is.na(object$properties$modularity1)) {
      # exclude modularity
      glob_rnames <- glob_rnames[-3]
      glob_names <- glob_names[-3]
      glob_names2 <- glob_names2[-3]
    } 

    if (showGlobPvals) {
      propdiffs <- matrix(0, nrow = length(glob_rnames), ncol = 5,
                          dimnames = list(glob_rnames,
                                          c(group1, paste0("  ",
                                                           group2),
                                            "   abs.diff.",
                                            "    p-value", " ")))
      propdiffs <- as.data.frame(propdiffs)
      
    } else {
      propdiffs <- matrix(0, nrow = length(glob_rnames), ncol = 3,
                          dimnames = list(glob_rnames,
                                          c(group1, paste0("  ",
                                                           group2),
                                            "   difference")))
      propdiffs <- as.data.frame(propdiffs)
    }
    
    for (i in 1:length(glob_names)) {
      propdiffs[i,1] <- 
        round(as.numeric(object$properties[paste0(glob_names[i], 1)]), 
              digits)
      
      propdiffs[i,2] <- 
        round(as.numeric(object$properties[paste0(glob_names[i],2)]), digits)
      
      propdiffs[i,3] <- 
        round(object$diffGlobal[[paste0("diff", glob_names2[i])]], digits)
      
      if (showGlobPvals) {
        propdiffs[i,4] <- 
          round(object$pvalDiffGlobal[[paste0("pval", glob_names2[i])]], 
                digitsPval)
        
        propdiffs[i,5] <- 
          .getSigCode(object$pvalDiffGlobal[[paste0("pval", glob_names2[i])]])
      }
    }
    
  } else {
    propdiffs <- NULL
  }
  
  #-----------------------------------------------------------------------------
  # LCC
  
  if (showGlobLCC) {
    glob_rnames_lcc <- c("Relative LCC size",
                         "Clustering coefficient",
                         "Modularity",
                         "Positive edge percentage",
                         "Edge density",
                         "Natural connectivity",
                         "Vertex connectivity",
                         "Edge connectivity",
                         "Average dissimilarity*",
                         "Average path length**")
    
    glob_names_lcc <-c("lccSizeRel",
                       "clustCoef",
                       "modularity",
                       "pep",
                       "density",
                       "natConnect",
                       "vertConnect",
                       "edgeConnect",
                       "avDiss",
                       "avPath")
    
    glob_names_lcc2 <- c("lccSizeRel",
                         "ClustCoef",
                         "Modul",
                         "PEP",
                         "Density",
                         "NatConnect",
                         "VertConnect",
                         "EdgeConnect",
                         "avDiss",
                         "avPath")
    
    
    if (is.na(object$properties$modularity1)) {
      # exclude modularity
      glob_rnames_lcc <- glob_rnames_lcc[-3]
      glob_names_lcc <- glob_names_lcc[-3]
      glob_names_lcc2 <- glob_names_lcc2[-3]
      
      if (is.na(object$propertiesLCC$vertConnect1)) {
        # exclude connectivity measures
        glob_rnames_lcc <- glob_rnames_lcc[-c(6, 7)]
        glob_names_lcc <- glob_names_lcc[-c(6, 7)]
        glob_names_lcc2 <- glob_names_lcc2[-c(6, 7)]
      } 
      
    } else if (is.na(object$propertiesLCC$vertConnect1)) {
      # exclude connectivity measures
      glob_rnames_lcc <- glob_rnames_lcc[-c(7, 8)]
      glob_names_lcc <- glob_names_lcc[-c(7, 8)]
      glob_names_lcc2 <- glob_names_lcc2[-c(7, 8)]
    }

    if (showGlobPvals) {
      propdiffs_lcc <- matrix(0, nrow = length(glob_rnames_lcc), 
                              ncol = 5,
                              dimnames = list(glob_rnames_lcc,
                                              c(group1, paste0("  ",
                                                               group2),
                                                "   abs.diff.",
                                                "    p-value", " ")))
      propdiffs_lcc <- as.data.frame(propdiffs_lcc)
      
    } else {
      propdiffs_lcc <- matrix(0, nrow = length(glob_rnames_lcc), 
                              ncol = 3,
                              dimnames = list(glob_rnames_lcc,
                                              c(group1, paste0("  ",
                                                               group2),
                                                "   difference")))
      propdiffs_lcc <- as.data.frame(propdiffs_lcc)
    }
    
    for (i in 1:length(glob_names_lcc)) {
      propdiffs_lcc[i, 1] <-
        round(as.numeric(object$propertiesLCC[paste0(glob_names_lcc[i], 1)]),
              digits)
      
      propdiffs_lcc[i, 2] <-
        round(as.numeric(object$propertiesLCC[paste0(glob_names_lcc[i], 2)]),
              digits)
      
      propdiffs_lcc[i, 3] <-
        round(object$diffGlobalLCC[[paste0("diff", glob_names_lcc2[i])]],
              digits)
      
      if (showGlobPvals) {
        propdiffs_lcc[i, 4] <-
          round(object$pvalDiffGlobalLCC[[paste0("pval", glob_names_lcc2[i])]],
                digitsPval)
        propdiffs_lcc[i, 5] <-
          .getSigCode(object$pvalDiffGlobalLCC[[paste0("pval",
                                                   glob_names_lcc2[i])]])
      }
    }
  } else {
    propdiffs_lcc <- NULL
  }
  
  #-----------------------------------------------------------------------------
  # Combine for connected networks (with only one component)
  
  if (!is_disconnected & showGlob) {
    propdiffs <- rbind(propdiffs[1, , drop = FALSE], 
                       propdiffs_lcc[-1, , drop = FALSE])
  } 
  
  #============================================================
  # centrality measures
  
  topProps <- NULL
  if (showCentr[1] != "none") {
    topDegNames <- names(sort(abs(object$diffCentr$diffDeg), 
                              decreasing = TRUE))[1:numbNodes]
    topBetwNames <- names(sort(abs(object$diffCentr$diffBetw), 
                               decreasing = TRUE))[1:numbNodes]
    topCloseNames <- names(sort(abs(object$diffCentr$diffClose), 
                                decreasing = TRUE))[1:numbNodes]
    topEigenNames <- names(sort(abs(object$diffCentr$diffEigen), 
                                decreasing = TRUE))[1:numbNodes]
    
    cols <- ifelse(is.null(object$pvalDiffCentr), 3, 5)
    if (cols == 3) {
      cnames <- c(group1, group2, "abs.diff.")
    } else {
      if (pAdjust) {
        cnames <- c(group1, group2, "abs.diff.", "adj.p-value", " ")
      } else {
        cnames <- c(group1, group2, "abs.diff.", "p-value", " ")
      }
      
    }
    
    # degree
    if (any(c("all", "degree") %in% showCentr)) {
      topDeg <- matrix(0, nrow = numbNodes, ncol = cols,
                       dimnames = list(topDegNames, cnames))
      topDeg <- as.data.frame(topDeg)
      
      topDeg[, 1] <-
        round(object$properties$deg1[topDegNames], digits)
      
      topDeg[, 2] <-
        round(object$properties$deg2[topDegNames], digits)
      
      topDeg[, 3] <-
        round(abs(object$diffCentr$diffDeg[topDegNames]), digits)
      
      if (cols == 5) {
        if (pAdjust) {
          topDeg[, 4] <-
            round(object$pvalDiffCentrAdjust$pAdjustDiffDeg[topDegNames],
                  digitsPval)
        } else {
          topDeg[, 4] <-
            round(object$pvalDiffCentr$pvalDiffDeg[topDegNames],
                  digitsPval)
        }
        topDeg[, 5] <-
          sapply(1:numbNodes, function(i) {
            .getSigCode(topDeg[, 4][i])
          })
      }
    } else {
      topDeg <- NULL
    }
    
    
    # betweenness centrality
    if (any(c("all", "betweenness") %in% showCentr)) {
      topBetw <- matrix(0,
                        nrow = numbNodes,
                        ncol = cols,
                        dimnames = list(topBetwNames, cnames))
      topBetw <- as.data.frame(topBetw)
      
      topBetw[, 1] <-
        round(object$properties$betw1[topBetwNames], digits)
      
      topBetw[, 2] <-
        round(object$properties$betw2[topBetwNames], digits)
      
      topBetw[, 3] <-
        round(abs(object$diffCentr$diffBetw[topBetwNames]), digits)
      
      if (cols == 5) {
        if (pAdjust) {
          topBetw[, 4] <-
            round(object$pvalDiffCentrAdjust$pAdjustDiffBetw[topBetwNames],
                  digitsPval)
          
        } else {
          topBetw[, 4] <-
            round(object$pvalDiffCentr$pvalDiffBetw[topBetwNames],
                  digitsPval)
          
        }
        topBetw[, 5] <- sapply(1:numbNodes, function(i) {
          .getSigCode(topBetw[, 4][i])
        })
      }
    } else {
      topBetw <- NULL
    }
    
    # closeness centrality
    if (any(c("all", "closeness") %in% showCentr)) {
      topClose <- matrix(0, nrow = numbNodes, ncol = cols,
                         dimnames = list(topCloseNames, cnames))
      topClose <- as.data.frame(topClose)
      
      topClose[, 1] <-
        round(object$properties$close1[topCloseNames], digits)
      
      topClose[, 2] <-
        round(object$properties$close2[topCloseNames], digits)
      
      topClose[, 3] <-
        round(abs(object$diffCentr$diffClose[topCloseNames]), digits)
      
      if (cols == 5) {
        if (pAdjust) {
          topClose[, 4] <-
            round(object$pvalDiffCentrAdjust$pAdjustDiffClose[topCloseNames],
                  digitsPval)
        } else {
          topClose[, 4] <-
            round(object$pvalDiffCentr$pvalDiffClose[topCloseNames],
                  digitsPval)
        }
        topClose[, 5] <- sapply(1:numbNodes, function(i) {
          .getSigCode(topClose[, 4][i])})
      }
    } else {
      topClose <- NULL
    }
    
    # eigenvector centrality
    if (any(c("all", "eigenvector") %in% showCentr)) {
      topEigen <- matrix(0, nrow = numbNodes, ncol = cols,
                         dimnames = list(topEigenNames, cnames))
      topEigen <- as.data.frame(topEigen)
      
      topEigen[, 1] <-
        round(object$properties$eigen1[topEigenNames], digits)
      
      topEigen[, 2] <-
        round(object$properties$eigen2[topEigenNames], digits)
      
      topEigen[, 3] <-
        round(abs(object$diffCentr$diffEigen[topEigenNames]), digits)
      
      if (cols == 5) {
        if (pAdjust) {
          topEigen[, 4] <-
            round(object$pvalDiffCentrAdjust$pAdjustDiffEigen[topEigenNames],
                  digitsPval)
        } else {
          topEigen[, 4] <-
            round(object$pvalDiffCentr$pvalDiffEigen[topEigenNames],
                  digitsPval)
        }
        topEigen[, 5] <-
          sapply(1:numbNodes, function(i) {
            .getSigCode(topEigen[, 4][i])
          })
      }
    } else {
      topEigen <- NULL
    }
    
    topProps <- list(topDeg = topDeg,
                     topBetw = topBetw,
                     topClose = topClose,
                     topEigen = topEigen)
  }
  
  out <- list(jaccmat = jaccmat,
              propdiffs = propdiffs,
              propdiffs_lcc = propdiffs_lcc,
              showGlobPvals = showGlobPvals,
              rand = rand,
              gcd = gcd,
              properties = object$properties,
              topProps = topProps,
              pvalDiffCentr = object$pvalDiffCentr,
              is_disconnected = is_disconnected,
              paramsProperties = object$paramsProperties,
              call = object$call,
              callArgs = object$callArgs)
  
  class(out) <- "summary.microNetComp"
  
  return(out)
}



#' @title Print Summary for objects of class \code{microNetComp}
#'
#' @param x object of class \code{summary.microNetComp} (returned by 
#'   \code{\link{summary.microNetComp}}).
#' @param ... not used.
#'
#' @method print summary.microNetComp
#'
#' @rdname summarize.microNetComp
#' @export
print.summary.microNetComp <- function(x, ...) {
  cat("\nComparison of Network Properties\n")
  cat("----------------------------------")
  
  cat("\nCALL: \n")
  dput(x$call, control = NULL)
  
  showGlob <- !(is.null(x$propdiffs) & is.null(x$propdiffs_lcc))
  
  if (showGlob) {
    cat("\n______________________________")
    cat("\nGlobal network properties\n")
    cat("`````````````````````````\n")
  }

  if (!is.null(x$propdiffs_lcc) & x$is_disconnected) {
    cat("Largest connected component (LCC):\n")
    print(x$propdiffs_lcc)
    cat("\n")
  }
  
  if (!is.null(x$propdiffs)) {
    cat("Whole network:\n")
    print(x$propdiffs)
  }
  
  if (showGlob) {
    cat("-----\n")
    
    if (x$showGlobPvals) {
      cat("p-values: one-tailed test with null hypothesis diff=0\n")
    }
    
    cat(" *: Dissimilarity = 1 - edge weight\n")
    
    if (x$paramsProperties$sPathNorm) {
      cat("**: Path length = Units with average dissimilarity\n")
    } else {
      cat("**: Path length = Sum of dissimilarities along the path\n")
    }
  }

  if (!is.null(x$jaccmat)) {
    cat("\n______________________________")
    cat("\nJaccard index (similarity betw. sets of most central nodes)\n")
    cat("```````````````````````````````````````````````````````````\n")
    print(x$jaccmat)
    cat("-----\n")
    cat("Jaccard index in [0,1] (1 indicates perfect agreement)\n")
  }

  if (!is.null(x$rand)) {
    cat("\n______________________________")
    cat("\nAdjusted Rand index (similarity betw. clusterings)\n")
    cat("``````````````````````````````````````````````````\n")
    print(x$rand)
    cat("-----\n")
    cat("ARI in [-1,1] with ARI=1: perfect agreement betw. clusterings
                   ARI=0: expected for two random clusterings\n")
    
    if (nrow(x$rand) == 2) {
      nPerm <- x$callArgs$nPermRand
      cat("p-value: permutation test (n=", nPerm, 
          ") with null hypothesis ARI=0\n", sep = "")
    }
  }
  
  if (!is.null(x$gcd)) {
    cat("\n______________________________")
    cat("\nGraphlet Correlation Distance\n")
    cat("`````````````````````````````\n")
    print(x$gcd)
    cat("-----\n")
    cat("GCD >= 0 (GCD=0 indicates perfect agreement between GCMs)\n")
    if (nrow(x$gcd) == 2) {
      cat("p-value: permutation test with null hypothesis GCD=0\n")
    }
  } 
  
  if (!is.null(x$topProps)) {
    cat("\n______________________________")
    if (x$is_disconnected && x$paramsProperties$centrLCC) {
      cat("\nCentrality measures")
      cat("\n- In decreasing order")
      cat("\n- Centrality of disconnected components is zero\n")
      cat("````````````````````````````````````````````````")
    } else {
      cat("\nCentrality measures")
      cat("\n- In decreasing order")
      cat("\n- Computed for the whole network\n")
      cat("````````````````````````````````````")
    }
    
    
    if (!is.null(x$topProps$topDeg)) {
      if (x$paramsProperties$weightDeg) {
        cat("\nDegree (weighted):\n")
      } else if (x$paramsProperties$normDeg) {
        cat("\nDegree (normalized):\n")
      } else {
        cat("\nDegree (unnormalized):\n")
      }
      print(x$topProps$topDeg)
    }
    
    if (!is.null(x$topProps$topBetw)) {
      if (x$paramsProperties$normBetw) {
        cat("\nBetweenness centrality (normalized):\n")
      } else {
        cat("\nBetweenness centrality (unnormalized):\n")
      }
      print(x$topProps$topBetw)
    }
    
    if (!is.null(x$topProps$topClose)) {
      if (x$paramsProperties$normBetw) {
        cat("\nCloseness centrality (normalized):\n")
      } else {
        cat("\nCloseness centrality (unnormalized):\n")
      }
      print(x$topProps$topClose)
    }
    
    if (!is.null(x$topProps$topEigen)) {
      if (x$paramsProperties$normBetw) {
        cat("\nEigenvector centrality (normalized):\n")
      } else {
        cat("\nEigenvector centrality (unnormalized):\n")
      }
      print(x$topProps$topEigen)
    }
    
  }
  
  cat("\n_________________________________________________________\n")
  cat("Significance codes: ***: 0.001, **: 0.01, *: 0.05, .: 0.1\n")
}
