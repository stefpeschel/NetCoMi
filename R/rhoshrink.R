#' @title Shrinkage estimator for proportionality rho
#'
#' @description Estimates the proportionality measure rho using shrinkage
#'   covariance estimation as proposed by Badri et al. (2020).
#'
#'   The implementation follows the accompanying R script from the
#'   NormCorr manuscript repository:
#'   \url{https://github.com/MichelleBadri/NormCorr-manuscript}. Since NetCoMi
#'   uses samples in rows and taxa/features in columns, the clr transformation
#'   is applied in the transpose-equivalent form.
#'
#' @param counts numeric matrix with samples in rows and taxa/features in columns.
#'   The matrix must contain positive values because a clr transformation is
#'   applied internally.
#' @param ... Additional arguments passed to \code{\link[corpcor]{cov.shrink}}.
#'
#' @return A list containing the following elements:
#' \tabular{ll}{
#'   \code{clrMat}\tab clr-transformed input matrix\cr
#'   \code{covMat}\tab shrinkage covariance matrix returned by
#'   \code{\link[corpcor]{cov.shrink}}\cr
#'   \code{assoMat}\tab matrix containing the estimated rho values}
#'
#' @references
#'   Badri M, Kurtz ZD, Bonneau R, Mueller CL (2020). "Shrinkage improves
#'   estimation of microbial associations under different normalization
#'   methods." \emph{NAR Genomics and Bioinformatics}, 2(4).
#'   \doi{10.1093/nargab/lqaa100}.\cr\cr
#'   Badri M, Kurtz ZD, Bonneau R, Mueller CL. "NormCorr-manuscript."
#'   \url{https://github.com/MichelleBadri/NormCorr-manuscript}.
#'
#' @export

rhoshrink <- function(counts, ...) {

  if (any(counts <= 0)) {
    stop("'counts' must contain positive values for clr transformation.")
  }

  counts_clr <- t(SpiecEasi::clr(t(counts), 2))
  shrunk_cov <- corpcor::cov.shrink(counts_clr, ...)
  p <- ncol(counts_clr)
  J <- matrix(rep(diag(shrunk_cov), p), p)
  rho <- 2 * shrunk_cov / (J + t(J))
  rho <- (rho + t(rho)) / 2
  rho <- matrix(as.numeric(rho), nrow = p, ncol = p)
  colnames(rho) <- rownames(rho) <- colnames(counts)

  list(clrMat = counts_clr, covMat = shrunk_cov, assoMat = rho)
}
