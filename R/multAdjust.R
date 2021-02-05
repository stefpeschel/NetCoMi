#' Multiple testing adjustment
#'
#' The functions adjusts a vector of p-values for multiple testing
#'
#' @param pvals numeric vector with p-values
#' @param adjust character specifying the method used for adjustment.
#'   Can be \code{"lfdr"}, \code{"adaptBH"}, or one of the methods provided by
#'   \code{\link[stats]{p.adjust}}.
#' @param trueNullMethod character indicating the method used for estimating the
#'   proportion of true null hypotheses from a vector of p-values. Used for the
#'   adaptive Benjamini-Hochberg method for multiple testing adjustment (chosen
#'   by \code{adjust = "adaptBH"}). Accepts the provided options of the
#'   \code{method} argument of \code{\link[limma]{propTrueNull}}:
#'   \code{"convest"}(default), \code{"lfdr"}, \code{"mean"}, and \code{"hist"}.
#'   Can alternatively be \code{"farco"} for
#'   the "iterative plug-in method" proposed by \cite{Farcomeni (2007)}.
#' @param verbose if \code{TRUE}, progress messages are returned.
#'
#' @references
#'   \insertRef{farcomeni2007some}{NetCoMi}
#'
#' @importFrom fdrtool fdrtool
#' @importFrom stats p.adjust
#' @export

multAdjust <- function(pvals, adjust, trueNullMethod, verbose){

  if(adjust == "lfdr"){
    
    if(verbose){
      message("")
      message("Execute fdrtool() ...")
    }

    pAdjust <- fdrtool::fdrtool(pvals, statistic = "pvalue", plot = FALSE,
                                verbose = verbose)$lfdr
    names(pAdjust) <- names(pvals)

  } else if(adjust == "adaptBH"){

    m <- length(pvals)
    ind <- m:1
    o <- order(pvals, decreasing = TRUE)
    ro <- order(o)

    if(trueNullMethod == "farco"){
      R <- 0
      iter <- TRUE

      while(iter){
        pTrueNull <- 1- (R/m)  # proportion of true null hypotheses
        pAdjust <- pmin(1, cummin(m * pTrueNull / ind * pvals[o]))[ro]
        R_new <- length(which(pAdjust < 0.05))
        iter <- R_new != R  # stop iteration if R_new==R
        R <- R_new
      }
      if(verbose){
        message("\n Proportion of true null hypotheses: ", round(pTrueNull, 2))
      }


    } else{
      # trueNullMethod must be one of "lfdr", "mean", "hist", or "convest"
      pTrueNull <- limma::propTrueNull(pvals, method = trueNullMethod)
      if(verbose){
        message("\n Proportion of true null hypotheses: ", round(pTrueNull, 2))
      }
      pAdjust <- pmin(1, cummin(m * pTrueNull / ind * pvals[o]))[ro]
    }

    names(pAdjust) <- names(pvals)

  } else{
    pAdjust <- stats::p.adjust(pvals, adjust)
  }

  return(pAdjust)
}
