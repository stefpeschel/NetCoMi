#' @title Estimate Correlations for Microbiome Data
#'
#' @description Function for estimating correlations between OTU pairs based on
#'   the method proposed by \cite{Friedman and Alm (2012)}.
#'
#' @details The function is written in C++ for a considerably increased runtime.
#'   The structure of the R code is adopted from r-sparcc package on GitHub.
#'
#' @param countData matrix containing microbiome data (read counts), where rows
#'   represent samples and columns the OTUs/taxa
#' @param exiter numeric; maximum number of iterations, where in each step the
#'   strongest correlation is excluded and the basis variances and correlations
#'   are re-estimated based on the reduced data
#' @param thresh numeric; threshold for excluding component pairs in the
#'   iteration steps (pairs are excluded if the magnitude of their correlation
#'   value is upon the threshold)
#' @param repEst numeric; estimation process is repeated 'repEst' times to
#'   account for sampling noise; new fraction samples are drawn from Dirichlet
#'   distribution in each iteration step
#' @param seed an optional seed for reproducibility of the results
#' @return \tabular{ll}{ \code{CORR}\tab matrix with estimated correlations\cr
#'   \code{COV}\tab estimated covariance matrix\cr \code{VBASIS}\tab vector with
#'   basis variances}
#'
#' @useDynLib NetCoMi, .registration = TRUE
#' @importFrom Rcpp sourceCpp evalCpp
#' @import RcppArmadillo
#' @importFrom gtools rdirichlet
#' @importFrom stats median
#'
#' @references
#'   \insertRef{friedman2012inferring}{NetCoMi}\cr
#'   \insertRef{sparccRimplement}{NetCoMi}
#'
#' @export
sparcc_cpp <- function(countData, exiter=20, thresh=0.1, repEst=100, seed = NULL){
  if(!is.null(seed)) set.seed(seed)

  x <- as.matrix(countData)
  xcol <- ncol(x)
  names <- colnames(x)


  sparcc_res <- C_sparcc(x, exiter = exiter, thresh = thresh, repEst = repEst)
  #sparcc_res <- .Call(`_NetCoMi_C_sparcc`, x, exiter, thresh, repEst)

  Vmat <- sparcc_res$Vmat
  Corlist <- array(unlist(sparcc_res$Corlist), dim = c(xcol, xcol, repEst))
  Covlist <- array(unlist(sparcc_res$Covlist), dim = c(xcol, xcol, repEst))
  fraclist <- array(unlist(sparcc_res$fraclist), dim = c(nrow(x), xcol, repEst))

  ## Compute variance basis and correlation
  vdef <- apply(Vmat,1,median)
  cor.def <- apply(Corlist, 1:2, median)
  fracs <- apply(fraclist, 1:2, median)
  dimnames(cor.def) <- list(names, names)
  dimnames(fracs) <- dimnames(x)

  ## Square root variances
  vdefsq <- vdef**0.5

  ## Compute covariance
  ttmp <- cor.def * vdefsq
  cov.def <- t(ttmp) * vdefsq
  dimnames(cov.def) <- list(names, names)

  ## Uncomment following lines for an alternative method
  ## x <- matrix(vdefsq,ncol=50,nrow=50, byrow=TRUE)
  ## y <- t(x)
  ## cov.def <- cor.def * x * y

  return(list(CORR = cor.def, COV = cov.def, VBASIS = vdef, fracs = fracs))
}
