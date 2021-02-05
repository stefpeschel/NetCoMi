#' @title Summary Method for Objects of Class netProps
#' @description  The main results returned from \code{\link{netCompare}} are
#'   printed in a well-arranged format.
#'
#' @param object object of class \code{netProps} inheriting from a call to
#'   \code{\link{netCompare}}.
#' @param groupNames character vector with two elements giving the group names
#'   corresponding to the two networks. If \code{NULL}, the names are adopted
#'   from \code{object}
#' @param pAdjust logical. If \code{TRUE} the permutation p-values adjusted for
#'   multiple testing are shown.
#' @param showCentr character vector indicating for which centrality measures
#'   the results (values for both groups, difference and p-values resulting from
#'   permutation tests) shall be printed. Possible values are "all", "degree",
#'   "betweenness", "closeness", "eigenvector" and "none".
#' @param numbNodes integer indicating for how many nodes the centrality
#'   values shall be printed. Defaults to 10 which means that the first 10 taxa
#'   with highest absolute group difference of the specific centrality
#'   measure are shown.
#' @param digits integer giving the number of decimal places to which the 
#'   results are rounded. Defaults to 3L. 
#' @param digitsPval integer giving the number of decimal places to which the 
#'   p-values are rounded. Defaults to 6L.
#'
#' @seealso \code{\link{netCompare}}
#'
#' @method summary microNetComp
#' @rdname summarize.microNetComp
#' @export
summary.microNetComp <- function(object, groupNames = NULL, pAdjust = TRUE,
                                 showCentr = "all", 
                                 numbNodes = 10L, digits = 3L, 
                                 digitsPval = 6L, ...){

  showCentr <- match.arg(showCentr, choices = c("all", "none", "degree", 
                                                "betweenness", "closeness", 
                                                "eigenvector"),
                         several.ok = TRUE)
  if("none" %in% showCentr) stopifnot(length(showCentr) == 1)
  numbNodes <- as.integer(numbNodes)
  digits <- as.integer(digits)
  digitsPval <- as.integer(digitsPval)
  
  stopifnot(numbNodes >= 1)
  numbNodes <- min(as.integer(numbNodes), length(object$diffCentr$diffDeg))
  
  
  if(is.null(groupNames)){
    group1 <- paste0("group '" , object$groups$group1, "'")
    group2 <- paste0("group '" , object$groups$group2, "'")
    
  } else{
    group1 <- groupNames[1]
    group2 <- groupNames[2]
  }
  
  namesCentr <- c("degree", "betweenness centr.", "closeness centr.", 
                  "eigenvec. centr.", "hub taxa")
  
  # jaccard index
  jacc <- numeric(5)
  p.greater <- p.less <- numeric(5)
  sig.greater <- sig.less <- character(5)
  
  codesig <- function(x){
    if(x <= 0.001){
      return("***")
    } else if(x <= 0.01){
      return("** ")
    } else if(x <= 0.05){
      return("*  ")
    } else if(x <= 0.1){
      return(".  ")
    } else{
      return(" ")
    }
  }
  
  for(i in 1:5){
    jacc[i] <- round(object[[i]]["jacc"], digits)
    p.greater[i] <- round(object[[i]]["p.greater"], digitsPval)
    p.less[i] <- round(object[[i]]["p.less"], digitsPval)
    sig.greater[i] <- codesig(p.greater[i])
    sig.less[i] <- codesig(p.less[i])
  }
  jaccmat <- data.frame(jacc, p.less, sig.less, p.greater, sig.greater)
  rownames(jaccmat) <- namesCentr
  colnames(jaccmat) <- c("Jacc", "  P(<=Jacc)", "   ", "P(>=Jacc)", "  ")
  
  #============================================================
  # rand index
  rand <- as.data.frame(matrix(0, nrow = 1, ncol = 2,
                               dimnames = list("", c("ARI", "      p-value"))))
  rand[1,1] <- round(object$randInd[1], digits)
  rand[1,2] <- round(object$randInd[2], digitsPval)
  
  #============================================================
  # global network properties
  
  glob_rnames <- c("Number of components", 
                   "Clustering coefficient",
                   "Moduarity",
                   "Positive edge percentage",
                   "Edge density",
                   "Natural connectivity")
  
  glob_rnames_lcc <- c("Relative LCC size", 
                       "Clustering coefficient",
                       "Moduarity",
                       "Positive edge percentage",
                       "Edge density",
                       "Natural connectivity",
                       "Vertex connectivity",
                       "Edge connectivity",
                       "Average dissimilarity*",
                       "Average path length**")
  
  glob_names <- c("nComp", "clustCoef", "modularity", "pep", "density", 
                  "natConnect")
  
  glob_names2 <- c("nComp", "ClustCoef", "Modul", "PEP", "Density", 
                   "NatConnect")
  
  glob_names_lcc <- c("lccSizeRel", "clustCoef", "modularity", "pep", 
                      "density", "natConnect", "vertConnect", "edgeConnect",
                      "avDiss", "avPath")
  
  glob_names_lcc2 <- c("lccSizeRel", "ClustCoef", "Modul", "PEP", 
                       "Density", "NatConnect", "VertConnect", "EdgeConnect",
                       "avDiss", "avPath")
  
  
  if(is.na(object$properties$modularity1)){
    # exclude modularity
    glob_rnames <- glob_rnames[-3]
    glob_rnames_lcc <- glob_rnames_lcc[-3]
    glob_names <- glob_names[-3]
    glob_names2 <- glob_names2[-3]
    glob_names_lcc <- glob_names_lcc[-3]
    glob_names_lcc2 <- glob_names_lcc2[-3]
    
    if(is.na(object$propertiesLCC$vertConnect1)){
      # exclude connectivity measures
      glob_rnames_lcc <- glob_rnames_lcc[-c(6,7)]
      glob_names_lcc <- glob_names_lcc[-c(6,7)]
      glob_names_lcc2 <- glob_names_lcc2[-c(6,7)]
    } 
    
  } else if(is.na(object$propertiesLCC$vertConnect1)){
    # exclude connectivity measures
    glob_rnames_lcc <- glob_rnames_lcc[-c(7,8)]
    glob_names_lcc <- glob_names_lcc[-c(7,8)]
    glob_names_lcc2 <- glob_names_lcc2[-c(7,8)]
  }
  
  if(is.null(object$pvalDiffGlobal)){ 
    propdiffs <- as.data.frame(matrix(0, nrow = length(glob_rnames), ncol = 3,
                                      dimnames = list(glob_rnames,
                                                      c(group1, paste0("  ",
                                                                       group2),
                                                        "   difference"))))
    
    propdiffs_lcc <- as.data.frame(matrix(0, nrow = length(glob_rnames_lcc), ncol = 3,
                                          dimnames = list(glob_rnames_lcc,
                                                          c(group1, paste0("  ",
                                                                           group2),
                                                            "   difference"))))
    
  } else{
    propdiffs <- as.data.frame(matrix(0, nrow = length(glob_rnames), ncol = 5,
                                      dimnames = list(glob_rnames,
                                                      c(group1, paste0("  ",
                                                                       group2),
                                                        "   abs.diff.",
                                                        "    p-value", " "))))
    
    propdiffs_lcc <- as.data.frame(matrix(0, nrow = length(glob_rnames_lcc), ncol = 5,
                                          dimnames = list(glob_rnames_lcc,
                                                          c(group1, paste0("  ",
                                                                           group2),
                                                            "   abs.diff.",
                                                            "    p-value", " "))))
  }
  
  # whole network
  for(i in 1:length(glob_names)){
    propdiffs[i,1] <- round(as.numeric(object$properties[paste0(glob_names[i], 1)]), 
                            digits)
    propdiffs[i,2] <- round(as.numeric(object$properties[paste0(glob_names[i],2)]), 
                            digits)
    propdiffs[i,3] <- round(object$diffGlobal[[paste0("diff", glob_names2[i])]], digits)
    
    if(!is.null(object$pvalDiffGlobal)){
      propdiffs[i,4] <- round(object$pvalDiffGlobal[[paste0("pval", glob_names2[i])]], digitsPval)
      propdiffs[i,5] <- codesig(object$pvalDiffGlobal[[paste0("pval", glob_names2[i])]])
    }
  }
  
  # largest connected component
  for(i in 1:length(glob_names_lcc)){
    propdiffs_lcc[i,1] <- round(as.numeric(object$propertiesLCC[paste0(glob_names_lcc[i], 1)]), 
                                digits)
    propdiffs_lcc[i,2] <- round(as.numeric(object$propertiesLCC[paste0(glob_names_lcc[i],2)]), 
                                digits)
    propdiffs_lcc[i,3] <- round(object$diffGlobalLCC[[paste0("diff", glob_names_lcc2[i])]], digits)
    
    if(!is.null(object$pvalDiffGlobalLCC)){
      propdiffs_lcc[i,4] <- round(object$pvalDiffGlobalLCC[[paste0("pval", glob_names_lcc2[i])]], digitsPval)
      propdiffs_lcc[i,5] <- codesig(object$pvalDiffGlobalLCC[[paste0("pval", glob_names_lcc2[i])]])
    }
  }
  
  if(object$properties$nComp1 == 1 && object$properties$nComp1 == 1){
    propdiffs <- rbind(propdiffs[1, , drop = FALSE], 
                       propdiffs_lcc[-1, , drop = FALSE])
    is_disconnected <- FALSE
  } else{
    is_disconnected <- TRUE
  }
  
  #============================================================
  # centrality measures
  
  topProps <- NULL
  if(showCentr[1] != "none"){
    topDegNames <- names(sort(abs(object$diffCentr$diffDeg), 
                              decreasing = TRUE))[1:numbNodes]
    topBetwNames <- names(sort(abs(object$diffCentr$diffBetw), 
                               decreasing = TRUE))[1:numbNodes]
    topCloseNames <- names(sort(abs(object$diffCentr$diffClose), 
                                decreasing = TRUE))[1:numbNodes]
    topEigenNames <- names(sort(abs(object$diffCentr$diffEigen), 
                                decreasing = TRUE))[1:numbNodes]
    
    cols <- ifelse(is.null(object$pvalDiffCentr), 3, 5)
    if(cols == 3){
      cnames <- c(group1, group2, "abs.diff.")
    } else{
      if(pAdjust){
        cnames <- c(group1, group2, "abs.diff.", "adj.p-value", " ")
      } else{
        cnames <- c(group1, group2, "abs.diff.", "p-value", " ")
      }
      
    }
    
    # degree
    if(any(c("all", "degree") %in% showCentr)){
      topDeg <- as.data.frame(matrix(0, nrow = numbNodes, ncol = cols,
                                     dimnames = list(topDegNames, cnames)))
      topDeg[,1] <- round(object$properties$deg1[topDegNames], digits)
      topDeg[,2] <- round(object$properties$deg2[topDegNames], digits)
      topDeg[,3] <- round(abs(object$diffCentr$diffDeg[topDegNames]), digits)
      if(cols == 5){
        if(pAdjust){
          topDeg[,4] <- round(object$pvalDiffCentrAdjust$pAdjustDiffDeg[topDegNames],
                              digitsPval)
        } else{
          topDeg[,4] <- round(object$pvalDiffCentr$pvalDiffDeg[topDegNames],
                              digitsPval)
        }
        topDeg[,5] <- sapply(1:numbNodes, function(i){ codesig(topDeg[,4][i])})
      }
    } else{
      topDeg <- NULL
    }
    
    
    # betweenness centrality
    if(any(c("all", "betweenness") %in% showCentr)){
      topBetw <- as.data.frame(matrix(0, nrow = numbNodes, ncol = cols,
                                      dimnames = list(topBetwNames, cnames)))
      topBetw[,1] <- round(object$properties$betw1[topBetwNames], digits)
      topBetw[,2] <- round(object$properties$betw2[topBetwNames], digits)
      topBetw[,3] <- round(abs(object$diffCentr$diffBetw[topBetwNames]), digits)
      if(cols == 5){
        if(pAdjust){
          topBetw[,4] <- round(object$pvalDiffCentrAdjust$pAdjustDiffBetw[topBetwNames],
                               digitsPval)
        } else{
          topBetw[,4] <- round(object$pvalDiffCentr$pvalDiffBetw[topBetwNames],
                               digitsPval)
        }
        topBetw[,5] <- sapply(1:numbNodes, function(i){ codesig(topBetw[,4][i])})
      }
    } else{
      topBetw <- NULL
    }
    
    # closeness centrality
    if(any(c("all", "closeness") %in% showCentr)){
      topClose <- as.data.frame(matrix(0, nrow = numbNodes, ncol = cols,
                                       dimnames = list(topCloseNames, cnames)))
      topClose[,1] <- round(object$properties$close1[topCloseNames], digits)
      topClose[,2] <- round(object$properties$close2[topCloseNames], digits)
      topClose[,3] <- round(abs(object$diffCentr$diffClose[topCloseNames]), digits)
      if(cols == 5){
        if(pAdjust){
          topClose[,4] <- round(object$pvalDiffCentrAdjust$pAdjustDiffClose[topCloseNames],
                                digitsPval)
        } else{
          topClose[,4] <- round(object$pvalDiffCentr$pvalDiffClose[topCloseNames],
                                digitsPval)
        }
        topClose[,5] <- sapply(1:numbNodes, function(i){ codesig(topClose[,4][i])})
      }
    } else{
      topClose <- NULL
    }
    
    # eigenvector centrality
    if(any(c("all", "eigenvector") %in% showCentr)){
      topEigen <- as.data.frame(matrix(0, nrow = numbNodes, ncol = cols,
                                       dimnames = list(topEigenNames, cnames)))
      topEigen[,1] <- round(object$properties$eigen1[topEigenNames], digits)
      topEigen[,2] <- round(object$properties$eigen2[topEigenNames], digits)
      topEigen[,3] <- round(abs(object$diffCentr$diffEigen[topEigenNames]), digits)
      if(cols == 5){
        if(pAdjust){
          topEigen[,4] <- round(object$pvalDiffCentrAdjust$pAdjustDiffEigen[topEigenNames],
                                digitsPval)
        } else{
          topEigen[,4] <- round(object$pvalDiffCentr$pvalDiffEigen[topEigenNames],
                                digitsPval)
        }
        topEigen[,5] <- sapply(1:numbNodes, function(i){ codesig(topEigen[,4][i])})
      }
    } else{
      topEigen <- NULL
    }
    
    topProps <- list(topDeg = topDeg, topBetw = topBetw, topClose = topClose,
                     topEigen = topEigen)
  }

  structure(list(jaccmat = jaccmat, propdiffs = propdiffs, 
                 propdiffs_lcc = propdiffs_lcc, rand = rand, 
                 properties = object$properties, topProps = topProps,
                 pvalDiffCentr = object$pvalDiffCentr,
                 is_disconnected = is_disconnected,
                 paramsProperties = object$paramsProperties,
                 call = object$call),
            class = "summary.microNetComp")
}



#' @title Print Summary for objects of class \code{microNetComp}
#'
#' @param x object of class \code{summary.microNetComp} inheriting from a call
#'   to \code{\link{summary.microNetComp}}.
#' @param ... not used
#'
#' @method print summary.microNetComp
#'
#' @rdname summarize.microNetComp
#' @export
print.summary.microNetComp <- function(x, ...){
  
  cat("\nComparison of Network Properties\n")
  cat("----------------------------------")
  
  cat("\nCALL: \n")
  dput(x$call, control = NULL)
  
  cat("\n______________________________")
  cat("\nGlobal network properties\n")
  cat("`````````````````````````\n")
  
  if(x$is_disconnected){
    cat("Largest connected component (LCC):\n")
    print(x$propdiffs_lcc)
    
    cat("\nWhole network:\n")
  }
  
  print(x$propdiffs)
  
  cat("-----\n")
  if(ncol(x$propdiffs) == 5){
    cat("p-values: one-tailed test with null hypothesis diff=0\n")
  }
  
  cat(" *: Dissimilarity = 1 - edge weight\n")
  
  if(x$paramsProperties$sPathNorm){
    cat("**Path length: Units with average dissimilarity\n")
  } else{
    cat("**Path length: Sum of dissimilarities along the path\n")
  }
  
  cat("\n______________________________")
  
  cat("\nJaccard index (similarity betw. sets of most central nodes)\n")
  cat("``````````````````````````````````````````````````````````\n")
  print(x$jaccmat)
  cat("-----\n")
  cat("Jaccard index ranges from 0 (compl. different) to 1 (sets equal)\n")
  
  cat("\n______________________________")
  cat("\nAdjusted Rand index (similarity betw. clusterings)\n")
  cat("``````````````````````````````````````````````````\n")
  print(x$rand)
  cat("-----\n")
  cat("ARI in [-1,1] with ARI=1: perfect agreement betw. clusterings,
                   ARI=0: expected for two random clusterings\n")
  cat("p-value: two-tailed test with null hypothesis ARI=0\n")
  

  cat("\n______________________________")
  if(!is.null(x$topProps)){
    if(x$is_disconnected && x$paramsProperties$centrLCC){
      cat("\nCentrality measures")
      cat("\n- In decreasing order")
      cat("\n- Centrality of disconnected components is zero\n")
      cat("````````````````````````````````````````````````")
    } else{
      cat("\nCentrality measures")
      cat("\n- In decreasing order")
      cat("\n- Computed for the complete network\n")
      cat("````````````````````````````````````")
    }
    
    
    if(!is.null(x$topProps$topDeg)){
      if(x$paramsProperties$weightDeg){
        cat("\nDegree (weighted):\n")
      } else if(x$paramsProperties$normDeg){
        cat("\nDegree (normalized):\n")
      } else{
        cat("\nDegree (unnormalized):\n")
      }
      print(x$topProps$topDeg)
    }
    
    if(!is.null(x$topProps$topBetw)){
      if(x$paramsProperties$normBetw){
        cat("\nBetweenness centrality (normalized):\n")
      } else{
        cat("\nBetweenness centrality (unnormalized):\n")
      }
      print(x$topProps$topBetw)
    }
    
    if(!is.null(x$topProps$topClose)){
      if(x$paramsProperties$normBetw){
        cat("\nCloseness centrality (normalized):\n")
      } else{
        cat("\nCloseness centrality (unnormalized):\n")
      }
      print(x$topProps$topClose)
    }
    
    if(!is.null(x$topProps$topEigen)){
      if(x$paramsProperties$normBetw){
        cat("\nEigenvector centrality (normalized):\n")
      } else{
        cat("\nEigenvector centrality (unnormalized):\n")
      }
      print(x$topProps$topEigen)
    }
    
  }
  
  cat("\n_________________________________________________________\n")
  cat("Significance codes: ***: 0.001, **: 0.01, *: 0.05, .: 0.1\n")
}