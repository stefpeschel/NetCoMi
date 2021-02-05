trans_to_diss <- function(x, dissFunc, dissFuncPar = NULL){

  if(is.function(dissFunc)){
    dissMat <- do.call(dissFunc, c(list(x), dissFuncPar))
    return(dissMat)

  } else if(dissFunc %in% c("signed", "unsigned", "signedPos")){

    xvec <- x[lower.tri(x)]

    if(dissFunc == "signed"){
      dissvec <- sqrt(0.5 * (1-xvec))
    } else if(dissFunc == "unsigned"){
      dissvec <- sqrt(1-xvec^2)
    } else {
      dissvec <- sqrt(0.5 * (1-xvec))
      dissvec[xvec < 0] <- 0
    }

    dissMat <- x
    dissMat[lower.tri(dissMat)] <- dissvec
    dissMat[upper.tri(dissMat)] <- t(dissMat)[upper.tri(t(dissMat))]
    dissMat[x == 0] <- Inf
    diag(dissMat) <- 0

  } else if(dissFunc == "TOMdiss"){
    if(is.null(dissFuncPar)){
      dissFuncPar$TOMType <- ifelse(any(x<0), "signed", "unsigned")
      dissFuncPar$verbose <- 0
    }
    diag(x) <- 0
    dissMat <- do.call(WGCNA::TOMdist, c(list(adjMat = x), dissFuncPar))
    dissMat[dissMat == 1] <- Inf
    dimnames(dissMat) <- dimnames(x)

  } else{
    stop("Argument 'dissFunc' not specified correctly.")
  }
  return(dissMat)
}


trans_to_sim <- function(x, simFunc, simFuncPar = NULL){
  if(is.function(simFunc)){
    simMat <- do.call(simFunc, c(list(x), simFuncPar))
    return(simMat)

  } else {
    x.tmp <- x[upper.tri(x)]
    x.tmp <- x[!is.infinite(x)]
    if(all(x.tmp <= 1)){
      simMat <- 1-x
      simMat[x == Inf] <- 0

    } else{
      simMat <- 1 / (1 + x)
    }

  }

  return(simMat)
}


trans_to_adja <- function(x, weighted){
  adjaMat <- x

  if(!weighted){
    adjaMat[adjaMat != 0] <- 1
  }
  return(adjaMat)
}


scale_diss <- function(x){
  xUpper <- x[upper.tri(x)]
  xScale <- (x - min(xUpper)) / (max(xUpper) - min(xUpper)  + 1)
  diag(xScale) <- 0
  xScale
}

