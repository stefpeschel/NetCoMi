#' @title Cumulative sum scaling factors
#'
#' @description Copy of the \code{cumNormMat} function from
#'   metagenomeSeq package. 
#'   
#' @importFrom matrixStats colQuantiles

cumNormMat <- function(obj, p = cumNormStatFast(obj), sl = 1000){
  x = obj
  xx = x
  xx[x == 0] <- NA
  qs = matrixStats::colQuantiles(xx, probs = p, na.rm = TRUE)
  newMat <- sapply(1:ncol(xx), function(i) {
    xx = (x[, i] - .Machine$double.eps)
    sum(xx[xx <= qs[i]])
  })
  nmat <- sweep(x, 2, newMat/sl, "/")
  return(nmat)
}