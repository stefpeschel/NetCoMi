#' @title Compute associations between taxa
#'
#' @description Computes associations between taxa or distances between subjects
#'   for a given read count matrix
#'
#' @param countMat numeric read count matrix, where rows are samples and columns
#'   are OTUs/taxa.
#' @param measure character giving the measure used for estimating associations
#'   or dissimilarities
#' @param measurePar optional list with parameters passed to the function for
#'   estimating associations/dissimilarities
#' @param verbose if \code{TRUE}, progress messages are returned.
#'
#' @importFrom stats cor
#' @importFrom WGCNA bicor
#' @importFrom SpiecEasi spiec.easi sparseiCov getOptCov getRefit symBeta
#'   getOptBeta sparcc
#' @importFrom SPRING SPRING
#' @importFrom vegan vegdist

calc_association <- function(countMat, measure, measurePar, verbose){

  if(measure == "pearson"){

    if(is.null(measurePar$use)) measurePar$use <- "complete.obs"

    assoMat <- stats::cor(countMat, use = measurePar$use, method = "pearson")

  } else if(measure == "spearman"){

    if(is.null(measurePar$use)) measurePar$use <- "complete.obs"

    assoMat <- stats::cor(countMat, use = measurePar$use, method = "spearman")

  } else if(measure == "bicor"){

    measurePar$x <- countMat
    measurePar$verbose <- ifelse(verbose == 3, 1, 0)
    if(verbose == 3){
      message("")
    }
    assoMat <- do.call(WGCNA::bicor, measurePar)

  } else if(measure == "sparcc"){

    if(is.null(measurePar$inner_iter)) measurePar$inner_iter <- 20
    if(is.null(measurePar$th)) measurePar$th <- 0.1
    if(is.null(measurePar$iter)) measurePar$iter <- 100
    
    assoMat <- SpiecEasi::sparcc(countMat,
                                  inner_iter = measurePar$inner_iter,
                                  th = measurePar$th,
                                  iter = measurePar$iter)$Cor
    
    colnames(assoMat) <- rownames(assoMat) <- colnames(countMat)

  } else if(measure == "cclasso"){

    measurePar$x <- countMat
    measurePar$counts <- FALSE
    measurePar$pseudo <- 0
    
    measurePar$verbose <- ifelse(verbose == 3, TRUE, FALSE)
    if(verbose == 3){
      message("")
    }

    assoMat <- do.call("cclasso", measurePar)$cor.w

  } else if(measure %in% c("reboot", "ccrepe")){

    measurePar$x <- countMat
    measurePar$verbose <- ifelse(verbose == 3, TRUE, FALSE)
    if(verbose == 3){
      message("")
    }

    if(is.null(measurePar$sim.score.args)){
      measurePar$sim.score.args <- list(method = "pearson")
    }

    assoMat <- do.call(ccrepe::ccrepe, measurePar)$sim.score
    diag(assoMat) <- 1

  } else if(measure == "spieceasi"){
    
    if(!is.null(measurePar$symBetaMode)){
      symBetaMode <- measurePar$symBetaMode
      measurePar$symBetaMode <- NULL
    } else{
      symBetaMode <- "ave"
    }

    measurePar$data <- countMat

    if(is.null(measurePar$method)) measurePar$method <- "mb"
    
    if(verbose == 3){
      measurePar$verbose <- TRUE
      message("")
    } else{
      measurePar$verbose <- FALSE
    }

    spiecres <- do.call(SpiecEasi::spiec.easi, measurePar)

    if(spiecres$est$method == "glasso"){
      secor <- stats::cov2cor(as.matrix(getOptCov(spiecres)))
      assoMat <- secor * SpiecEasi::getRefit(spiecres)
    } else{
      assoMat <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(spiecres), 
                                    mode = symBetaMode)
    }

    assoMat <- as.matrix(assoMat)

    colnames(assoMat) <- rownames(assoMat) <- colnames(countMat)
    diag(assoMat) <- 1

  } else if(measure == "spring"){
    
    if(!is.null(measurePar$symBetaMode)){
      symBetaMode <- measurePar$symBetaMode
      measurePar$symBetaMode <- NULL
    } else{
      symBetaMode <- "ave"
    }

    measurePar$data <- countMat

    if(is.null(measurePar$lambdaseq)) measurePar$lambdaseq <- "data-specific"
    if(is.null(measurePar$ncores)) measurePar$ncores <- 1
    
    if(verbose == 3){
      measurePar$verbose <- TRUE
      message("")
    } else{
      measurePar$verbose <- FALSE
    }

    springres <- do.call(SPRING::SPRING, measurePar)

    opt.K <- springres$output$stars$opt.index

    assoMat <- as.matrix(SpiecEasi::symBeta(springres$output$est$beta[[opt.K]],
                                            mode = symBetaMode))

    colnames(assoMat) <- rownames(assoMat) <- colnames(countMat)
    diag(assoMat) <- 1

  } else if(measure == "gcoda"){

    measurePar$x <- countMat
    measurePar$counts <- FALSE
    measurePar$pseudo <- 0
    
    measurePar$verbose <- ifelse(verbose == 3, TRUE, FALSE)
    if(verbose == 3){
      message("")
    }

    gcodares <- do.call("gcoda", measurePar)
    assoMat <- gcodares$opt.icov

    assoLow <- assoMat[lower.tri(assoMat, diag = FALSE)]
    assoUp <- t(assoMat)[lower.tri(assoMat, diag = FALSE)]
    assoMean <- rowMeans(cbind(assoLow, assoUp))
    assoMat <- matrix(0, ncol = ncol(assoMat), nrow = nrow(assoMat))
    assoMat[lower.tri(assoMat, diag = FALSE)] <- assoMean
    assoMat <- t(assoMat)
    assoMat[lower.tri(assoMat, diag = FALSE)] <- assoMean

    if(any(abs(assoMat) > 1)){
      warning("Estimated correlations outside [-1,1] detected (replaced by -1 or 1).")
      assoMat[assoMat < (-1)] <- -1
      assoMat[assoMat > 1] <- 1
    }

    colnames(assoMat) <- rownames(assoMat) <- colnames(countMat)

  } else if(measure == "propr"){

    measurePar$counts <- countMat
    if(is.null(measurePar$metric)) measurePar$metric <- "rho"

    assoMat <- suppressMessages(assoMat <- do.call(propr::propr, 
                                                   measurePar)@matrix)

  } else if(measure == "bray"){

    assoMat <- as.matrix(vegan::vegdist(countMat, method = "bray", 
                                        binary = FALSE, diag = FALSE, 
                                        upper = FALSE, na.rm = FALSE))

  } else if(measure == "euclidean"){

    assoMat <- as.matrix(vegan::vegdist(countMat, method = "euclidean",
                                        binary = FALSE, diag = FALSE, 
                                        upper = FALSE, na.rm = FALSE))

  } else if(measure %in% c("kld", "jeffrey")){

    assoMat <- matrix(0, ncol = nrow(countMat), nrow = nrow(countMat))

    if(is.null(measurePar$base)) measurePar$base <- exp(1)

    for(i in 1:nrow(countMat)){
      j = i + 1
      while(j <= nrow(countMat)){
        assoMat[i, j] <- LaplacesDemon::KLD(countMat[i, ], countMat[j, ],
                                            base = measurePar$base)$mean.sum.KLD
        j <- j + 1
      }
    }

    assoMat[lower.tri(assoMat)] <- t(assoMat)[lower.tri(t(assoMat))]

    if(measure == "jeffrey"){
      assoMat <- assoMat * 2
    }

    colnames(assoMat) <- rownames(assoMat) <- rownames(countMat)
    diag(assoMat) <- 1

  } else if(measure == "jsd"){

    assoMat <- matrix(0, ncol = nrow(countMat), nrow = nrow(countMat))

    if(is.null(measurePar$base)) measurePar$base <- exp(1)

    for(i in 1:nrow(countMat)){
      j = i + 1
      while(j <= nrow(countMat)){
        P <- countMat[i, ]
        Q <- countMat[j, ]
        M <- 0.5 * (P + Q)
        
        assoMat[i, j] <- 0.5 * 
          LaplacesDemon::KLD(P, M, base = measurePar$base)$sum.KLD.px.py +
          0.5 * LaplacesDemon::KLD(Q, M, base = measurePar$base)$sum.KLD.px.py
        
        j <- j + 1
      }
    }

    assoMat[lower.tri(assoMat)] <- t(assoMat)[lower.tri(t(assoMat))]

    colnames(assoMat) <- rownames(assoMat) <- rownames(countMat)
    diag(assoMat) <- 1

  } else if(measure == "ckld"){
    assoMat <- matrix(0, ncol = nrow(countMat), nrow = nrow(countMat))

    if(is.null(measurePar$base)) measurePar$base <- exp(1)

    for(i in 1:nrow(countMat)){
      j = i + 1
      while(j <= nrow(countMat)){
        A1 <- mean(countMat[i, ] / countMat[j, ])
        A2 <- mean(countMat[j, ] / countMat[i, ])
        assoMat[i, j] <- ncol(countMat) / 2 * log(A1 * A2,
                                                  base = measurePar$base)
        j <- j + 1
      }
    }

    assoMat[lower.tri(assoMat)] <- t(assoMat)[lower.tri(t(assoMat))]

    colnames(assoMat) <- rownames(assoMat) <- rownames(countMat)
    diag(assoMat) <- 1

  } else if(measure == "aitchison"){

    measurePar$x.f <- countMat
    measurePar$mar <- 1

    if(is.null(measurePar$base)) measurePar$base <- exp(1)
    
    countMat_clr<- t(do.call(SpiecEasi::clr, measurePar))

    assoMat <- as.matrix(vegan::vegdist(countMat_clr, method = "euclidean"))

  } else{
    stop("Possible association/dissimilarity measures are: 'pearson', 'spearman', 'bicor',
         'sparcc', 'cclasso', 'reboot', 'ccrepe', 'bray', 'kld', 'jeffrey', 'jsd', 'aitchison',
         'prop', 'spieceasi', 'spring', 'gcoda'.")
  }


  if(any(is.na(assoMat))){
    assoMat[is.na(assoMat)] <- 0
    warning("Association matrix contains NAs (replaced by zeros).")
  }

  return(assoMat)
}

