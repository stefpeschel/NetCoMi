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
                                 showCentr = "all", numbNodes = 10L, digits = 3L, 
                                 digitsPval = 6L, ...){

  showCentr <- match.arg(showCentr, choices = c("all", "none", "degree", "betweenness",
                                                "closeness", "eigenvector"),
                         several.ok = TRUE)
  if("none" %in% showCentr) stopifnot(length(showCentr) == 1)
  numbNodes <- as.integer(numbNodes)
  digits <- as.integer(digits)
  digitsPval <- as.integer(digitsPval)

  stopifnot(numbNodes >= 1)
  numbNodes <- min(as.integer(numbNodes), length(object$diffs$diffDeg))


  if(is.null(groupNames)){
    group1 <- paste0("group '" , object$groups$group1, "'")
    group2 <- paste0("group '" , object$groups$group2, "'")

  } else{
    group1 <- groupNames[1]
    group2 <- groupNames[2]
  }


  names <- c("degree", "betweenness centr.", "closeness centr.", 
             "eigenvec. centr.", "hub taxa")

  if(is.na(object$properties$vertConnect1)){
    names2 <- c("avPath", "clustCoef", "modul", "density")
  } else{
    names2 <- c("avPath", "clustCoef", "modul", "vertConnect", 
                "edgeConnect", "density")
  }

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
  rownames(jaccmat) <- names
  colnames(jaccmat) <- c("Jacc", "  P(<=Jacc)", "   ", "P(>=Jacc)", "  ")

  #============================================================
  # rand index
  rand <- as.data.frame(matrix(0, nrow = 1, ncol = 2,
                               dimnames = list("", c("ARI", "      p-value"))))
  rand[1,1] <- round(object$randInd[1], digits)
  rand[1,2] <- round(object$randInd[2], digitsPval)

  #============================================================
  # global network properties

  if(is.na(object$properties$vertConnect1)){
    globProps <- c("average path length",
                   "clustering coeff.",
                   "modularity",
                   "edge density")
  } else{
    globProps <- c("average path length",
                   "clustering coeff.",
                   "modularity",
                   "vertex connectivity",
                   "edge connectivity",
                   "edge density")
  }

  if(length(object$avPath) == 2){
    propdiffs <- as.data.frame(matrix(0, nrow = length(globProps), ncol = 5,
                                      dimnames = list(globProps,
                                                      c(group1, paste0("  ",
                                                                       group2),
                                                        "   abs. diff.",
                                                        "    p-value", " "))))

    for(i in 1:length(globProps)){
      propdiffs[i,1] <- round(as.numeric(object$properties[paste0(names2[i],1)]), 
                              digits)
      propdiffs[i,2] <- round(as.numeric(object$properties[paste0(names2[i],2)]), 
                              digits)
      propdiffs[i,3] <- round(object[[names2[i]]][1], digits)
      propdiffs[i,4] <- round(object[[names2[i]]][2], digitsPval)
      propdiffs[i,5] <- codesig(object[[names2[i]]][2])
    }

  } else{
    propdiffs <- as.data.frame(matrix(0, nrow = length(globProps), ncol = 3,
                                      dimnames = list(globProps,
                                                      c(group1, paste0("  ",
                                                                       group2),
                                                        "   difference"))))
    for(i in 1:length(globProps)){
      propdiffs[i,1] <- round(as.numeric(object$properties[paste0(names2[i],1)]), 
                              digits)
      propdiffs[i,2] <- round(as.numeric(object$properties[paste0(names2[i],2)]), 
                              digits)
      propdiffs[i,3] <- round(object[[i+5]][1], digits)
    }
  }

  #============================================================
  # centrality measures

  topProps <- NULL
  if(showCentr[1] != "none"){
    l <- length(object$diffs$diffDeg)

    topDegNames <- names(sort(abs(object$diffs$diffDeg), 
                              decreasing = TRUE))[1:min(numbNodes, l)]
    topBetwNames <- names(sort(abs(object$diffs$diffBetw), 
                               decreasing = TRUE))[1:min(numbNodes, l)]
    topCloseNames <- names(sort(abs(object$diffs$diffClose), 
                                decreasing = TRUE))[1:min(numbNodes, l)]
    topEigenNames <- names(sort(abs(object$diffs$diffEigen), 
                                decreasing = TRUE))[1:min(numbNodes, l)]

    cols <- ifelse(is.null(object$pvalDiffCentr), 3, 5)
    if(cols == 3){
      cnames <- c(group1, group2, "abs. diff.")
    } else{
      if(pAdjust){
        cnames <- c(group1, group2, "abs. diff.", "adj. p-value", " ")
      } else{
        cnames <- c(group1, group2, "abs. diff.", "p-value", " ")
      }

    }

    # degree
    if(any(c("all", "degree") %in% showCentr)){
      topDeg <- as.data.frame(matrix(0, nrow = numbNodes, ncol = cols,
                                     dimnames = list(topDegNames, cnames)))
      topDeg[,1] <- round(object$properties$deg1[topDegNames], digits)
      topDeg[,2] <- round(object$properties$deg2[topDegNames], digits)
      topDeg[,3] <- round(abs(object$diffs$diffDeg[topDegNames]), digits)
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
      topBetw[,3] <- round(abs(object$diffs$diffBetw[topBetwNames]), digits)
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
      topClose[,3] <- round(abs(object$diffs$diffClose[topCloseNames]), digits)
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
      topEigen[,3] <- round(abs(object$diffs$diffEigen[topEigenNames]), digits)
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

  structure(list(jaccmat = jaccmat, propdiffs = propdiffs, rand = rand,
                 properties = object$properties, topProps = topProps,
                 pvalDiffCentr = object$pvalDiffCentr, call = object$call),
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

  cat("\n\nCALL: \n")
  dput(x$call, control = NULL)

  cat("\n\nJaccard index (similarity betw. sets of most central nodes):\n")
  cat("`````````````\n")
  print(x$jaccmat)
  cat("-----\n")
  cat("Jaccard index ranges from 0 (compl. different) to 1 (sets equal)\n")

  cat("\n\nGlobal network properties:\n")
  cat("``````````````````````````\n")
  print(x$propdiffs)

  if(ncol(x$propdiffs) == 5){
    cat("-----\n")
    cat("p-values: one-tailed test with null hypothesis diff=0\n")
  }

  cat("\n\nAdjusted Rand index (similarity betw. clusterings):\n")
  cat("```````````````````\n")
  print(x$rand)
  cat("-----\n")
  cat("ARI in [-1,1] with ARI=1: perfect agreement betw. clusterings,
                   ARI=0: expected for two random clusterings\n")
  cat("p-value: two-tailed test with null hypothesis ARI=0\n")


  if(!is.null(x$topProps)){
    cat("\n\nCentrality measures (sorted by decreasing difference):\n")
    cat("````````````````````")

    if(!is.null(x$topProps$topDeg)){
      cat("\nDegree:\n")
      print(x$topProps$topDeg)
    }

    if(!is.null(x$topProps$topBetw)){
      cat("\nBetweenness centrality:\n")
      print(x$topProps$topBetw)
    }

    if(!is.null(x$topProps$topClose)){
      cat("\nCloseness centrality:\n")
      print(x$topProps$topClose)
    }

    if(!is.null(x$topProps$topEigen)){
      cat("\nEigenvector centrality:\n")
      print(x$topProps$topEigen)
    }

  }

  cat("\n--------------------------------------------------------\n")
  cat("Significance codes: ***: 0.001, **: 0.01, *: 0.05, .: 0.1\n")
}
