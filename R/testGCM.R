#' @title Test GCM(s) for statistical significance
#' 
#' @description The function tests whether graphlet correlations (entries of 
#'   the GCM) are significantly different from zero.\cr\cr
#'   If two GCMs are given, the graphlet correlations of the two networks are 
#'   tested for being significantly different, i.e., Fishers z-test 
#'   is performed to test if the absolute differences between graphlet 
#'   correlations are significantly different from zero. 
#' @param obj1 object of class \code{GCM} or \code{GCD} returned by
#'   \code{\link{calcGCM}} or \code{\link{calcGCD}}. See details.
#' @param obj2 optional object of class \code{GCM} returned by 
#'   \code{\link{calcGCM}}. See details.
#' @param adjust character indicating the method used for multiple testing
#'   adjustment.
#'   Possible values are "lfdr" (default) for local
#'   false discovery rate correction (via \code{\link[fdrtool]{fdrtool}}),
#'   "adaptBH" for the adaptive Benjamini-Hochberg method \cite{(Benjamini and
#'   Hochberg, 2000)}, or one of the methods provided by
#'   \code{\link[stats]{p.adjust}}.
#' @param lfdrThresh defines a threshold for the local fdr if "lfdr" is chosen
#'   as method for multiple testing correction. Defaults to 0.2 meaning
#'   that differences with a corresponding local fdr less than or equal to 0.2
#'   are identified as significant.
#' @param trueNullMethod character indicating the method used for estimating the
#'   proportion of true null hypotheses from a vector of p-values. Used for the
#'   adaptive Benjamini-Hochberg method for multiple testing adjustment (chosen
#'   by \code{adjust = "adaptBH"}). Accepts the provided options of the
#'   \code{method} argument of \code{\link[limma]{propTrueNull}}: "convest"
#'   (default), "lfdr", "mean", and "hist". Can alternatively be "farco" for
#'   the "iterative plug-in method" proposed by \cite{Farcomeni (2007)}.
#' @param alpha numeric value between 0 and 1 giving the desired significance 
#'   level.
#' @param verbose logical. If \code{TRUE} (default), progress messages 
#'   are printed.
#'   
#' @details
#'   By applying Student's t-test to the Fisher-transformed correlations, 
#'   all entries of the GCM(s) are tested for being 
#'   significantly different from zero:\cr\cr
#'   H0:  gc_ij = 0  vs.  H1: gc_ij != 0,\cr\cr
#'   with gc_ij being the graphlet correlations.\cr\cr
#'   
#'   If both GCMs are given or \code{obj1} is of class \code{GCD}, the absolute
#'   differences between graphlet correlations are tested for being different 
#'   from zero using Fisher's z-test. The hypotheses are:\cr\cr
#'   H0: |d_ij| = 0 vs. H1: |d_ij| > 0,\cr\cr
#'   where d_ij = gc1_ij - gc2_ij
#' @return A list with the following elements:
#'   \tabular{ll}{
#'   \code{gcm1}\tab Graphlet Correlatoin Matrix GCM1\cr
#'   \code{pvals1}\tab Matrix with p-values (H0: gc1_ij = 0)\cr
#'   \code{padjust1}\tab Matrix with adjusted p-values\cr
#'   }
#'   \cr
#'   Additional elements if two GCMs are given:
#'   \tabular{ll}{
#'   \code{gcm2}\tab Graphlet Correlatoin Matrix GCM2\cr
#'   \code{pvals2}\tab Matrix with p-values (H0: gc2_ij = 0)\cr
#'   \code{padjust2}\tab Matrix with adjusted p-values\cr
#'   \code{diff}\tab Matrix with differences between graphlet correlations 
#'   (GCM1 - GCM2)\cr
#'   \code{absDiff}\tab Matrix with absolute differences between graphlet 
#'   correlations (|GCM1 - GCM2|)\cr
#'   \code{pvalsDiff}\tab Matrix with p-values (H0: |gc1_ij - gc2_ij| = 0)\cr
#'   \code{pAdjustDiff}\tab Matrix with adjusted p-values\cr
#'   \code{sigDiff}\tab Same as \code{diff}, but non-significant differences 
#'   are set to zero.\cr
#'   \code{sigAbsDiff}\tab Same as \code{absDiff}, but non-significant 
#'   values are set to zero.
#'   }
#' 
#' @examples 
#' # See help page of calcGCD()
#' ?calcGCD
#' 
#' @export

testGCM <- function(obj1, 
                    obj2 = NULL,
                    adjust = "adaptBH",
                    lfdrThresh = 0.2,
                    trueNullMethod = "convest",
                    alpha = 0.05,
                    verbose = TRUE) {
  
  if (!(inherits(obj1, "GCM") | inherits(obj1, "GCD"))) {
    stop("\"obj1\" must be of class \"GCM\" or \"GCD\" returned by ", 
         "calcGCM() or calcGCD().")
  }
  
  if (!is.null(obj2)) {
    if (!(inherits(obj2, "GCM"))) {
      stop("\"obj2\" must be of class \"GCM\" returned by calcGCM().")
    }
  }
  
  adjust <- match.arg(adjust, c(p.adjust.methods, "lfdr", "adaptBH"))
  
  if (!is.numeric(lfdrThresh) | length(lfdrThresh) != 1 | lfdrThresh > 1 |
      lfdrThresh < 0) {
    stop("\"lfdrThresh\" must be a numeric value in [0, 1].")
  }
  
  trueNullMethod <- match.arg(trueNullMethod, 
                              c("convest", "lfdr", "mean", "hist", "farco"))
  
  if (!is.numeric(alpha) | length(alpha) != 1 | alpha < 0 | alpha > 1) {
    stop("\"alpha\" must be a numeric value in [0, 1].")
  }
  
  stopifnot(is.logical(verbose))
  
  #-----------------------------------------------------------------------------
  # Initialize variables
  
  if (inherits(obj1, "GCD")) {
    gcm1 <- obj1$gcm1
    n1 <- nrow(obj1$ocount1)
    
    gcm2 <- obj1$gcm2
    n2 <- nrow(obj1$ocount2)
    
  } else {
    gcm1 <- obj1$gcm
    n1 <- nrow(obj1$ocount)
    
    if (!is.null(obj2)) {
      gcm2 <- obj2$gcm
      
      n2 <- nrow(obj2$ocount)
      
    } else {
      gcm2 <- n2 <- NULL
    }
  }
  
  # Matrix for p-values
  pMat <- gcm1
  pMat[,] <- 1
  
  #-----------------------------------------------------------------------------
  # Test GCM1
  
  if (verbose) {
    message("Perform Student's t-test for GCM1 ... ")
  }
  
  corvec1 <- gcm1[lower.tri(gcm1)]
  
  df <- n1 - 2
  
  tstat <- (corvec1 * sqrt(df)) / sqrt(1 - corvec1^2)
  
  pVec1 <- 2 * (1 - stats::pt(abs(tstat), df))
  
  pVec1[pVec1 > 1] <- 1
  
  pvals1 <- pMat
  pvals1[lower.tri(pvals1)] <- pVec1
  pvals1[upper.tri(pvals1)] <- t(pvals1)[upper.tri(t(pvals1))]
  
  out <- list()
  
  out$gcm1 <- gcm1
  out$pvals1 <- pvals1
  
  if (adjust != "none") {
    
    if (verbose) {
      message(" Adjust for multiple testing ... ")
    }
    
    paVec1 <- multAdjust(pvals = pVec1, 
                         adjust = adjust,
                         trueNullMethod = "lfdr",
                         verbose = verbose)
    
    pAdjust1 <- pMat
    pAdjust1[lower.tri(pAdjust1)] <- paVec1
    pAdjust1[upper.tri(pAdjust1)] <- t(pAdjust1)[upper.tri(t(pAdjust1))]
    
    out$pAdjust1 <- pAdjust1
  }
  
  if (verbose) message("Done.")
  
  #-----------------------------------------------------------------------------
  # Test GCM2
  
  if (!is.null(gcm2)) {
    
    if (verbose) {
      message("\nPerform Student's t-test for GCM2 ... ")
    }
    
    corvec2 <- gcm2[lower.tri(gcm2)]
    
    df <- n2 - 2
    
    tstat <- (corvec2 * sqrt(df)) / sqrt(1 - corvec2^2)
    
    pVec2 <- 2 * (1 - stats::pt(abs(tstat), df))
    
    pVec2[pVec2 > 1] <- 1
    
    pvals2 <- pMat
    pvals2[lower.tri(pvals2)] <- pVec2
    pvals2[upper.tri(pvals2)] <- t(pvals2)[upper.tri(t(pvals2))]
    
    out$gcm2 <- gcm2
    out$pvals2 <- pvals2
    
    if (adjust != "none") {
      if (verbose) {
        message(" Adjust for multiple testing ... ")
      }
      
      paVec2 <- multAdjust(pvals = pVec2, 
                           adjust = adjust,
                           trueNullMethod = trueNullMethod, 
                           verbose = verbose)
      
      pAdjust2 <- pMat
      pAdjust2[lower.tri(pAdjust2)] <- paVec2
      pAdjust2[upper.tri(pAdjust2)] <- t(pAdjust2)[upper.tri(t(pAdjust2))]
      
      out$pAdjust2 <- pAdjust2
    }
    
    if (verbose) message("Done.")
  } 
  
  #-----------------------------------------------------------------------------
  # Test GCM1 - GCM2 for significance
  
  if (!is.null(gcm2)) {
    
    if (verbose) {
      message("\nTest GCM1 and GCM2 for differences ... ")
    }
    
    # Avoid extreme z values
    corvec1[corvec1 > 0.9999] <- 0.9999
    corvec2[corvec2 > 0.9999] <- 0.9999
    corvec1[corvec1 < -0.9999] <- -0.9999
    corvec2[corvec2 < -0.9999] <- -0.9999
    
    # z-transform correlations
    z1 <- atanh(corvec1)
    z2 <- atanh(corvec2)
    
    diff_z <- (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
    
    pDiffVec <- 2 * (1 - pnorm(abs(diff_z)))
    
    pvalsDiff <- pMat
    pvalsDiff[lower.tri(pvalsDiff)] <- pDiffVec
    pvalsDiff[upper.tri(pvalsDiff)] <- t(pvalsDiff)[upper.tri(t(pvalsDiff))]
    
    # Matrix with differences in correlations
    diffs <- gcm1 - gcm2
    absDiffs <- abs(gcm1 - gcm2)
    diag(diffs) <- diag(absDiffs) <-0
    
    out$diff <- diffs
    out$absDiff <- absDiffs
    out$pvalsDiff <- pvalsDiff
    
    # Adjust for multiple testing
    if (adjust != "none") {
      
      if (verbose) {
        message(" Adjust for multiple testing ... ")
      }
      
      pAdjustVec <- multAdjust(pvals = pDiffVec, 
                               adjust = adjust,
                               trueNullMethod = trueNullMethod, 
                               verbose = verbose)
      
      # Matrix with adjusted p-values
      pAdjustDiff <- pMat
      pAdjustDiff[lower.tri(pAdjustDiff)] <- pAdjustVec
      pAdjustDiff[upper.tri(pAdjustDiff)] <- 
        t(pAdjustDiff)[upper.tri(t(pAdjustDiff))]
      
      out$pAdjustDiff <- pAdjustDiff
      
      pvalsDiff <- pAdjustDiff
    }
    
    sigDiffs <- diffs
    sigDiffs[pvalsDiff > alpha] <- 0
    
    sigAbsDiffs <- absDiffs
    sigAbsDiffs[pvalsDiff > alpha] <- 0
    
    out$sigDiff <- sigDiffs
    out$sigAbsDiff <- sigAbsDiffs
    
    if (verbose) message("Done.")
  }
  
  return(out)
}



