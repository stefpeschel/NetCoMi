#' @title Bootstrap Procedure for Testing Statistical Significance of
#'   Correlation Values
#'
#' @description Statistical significance of correlations between pairs of
#'   taxonomic units is tested using a bootstrap procedure as proposed by
#'   \cite{Friedman and Alm (2012)}.
#'
#' @param countMat matrix containing microbiome data (read counts) for which
#'   the correlations are calculated (rows represent samples, columns represent
#'   taxa)
#' @param assoMat matrix containing associations that have been estimated for
#'   \code{countMat}
#' @param nboot number of bootstrap samples
#' @param measure character specifying the method used for computing the
#'   associations between taxa.
#' @param measurePar list with parameters passed to the function for computing
#'   associations/dissimilarities. See details for the respective functions.
#' @param parallel if \code{TRUE}, the bootstrap iterations are executed
#'   parallel
#' @param cores number of CPU cores used for parallelization
#' @param logFile character defining a log file, where the number of iteration
#'   is stored. If NULL, no log file is created. wherein the current iteration
#'   numbers are stored
#' @param verbose logical; if \code{TRUE}, the iteration numbers are printed to
#'   the R console
#' @param seed an optional seed for reproducibility of the results
#' @return \tabular{ll}{ \code{pvals}\tab calculated p-values\cr
#'   \code{corrMat}\tab estimated correlation matrix}
#' @references{\insertRef{friedman2012inferring}{NetCoMi}}

boottest <- function(countMat, assoMat, nboot = 1000, measure, measurePar,
                     parallel = FALSE, cores = 4, logFile = NULL,
                     verbose = TRUE, seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  if(parallel){

    if(!is.null(seed)) seeds <- sample.int(1e8, size = nboot)

    cl <- makeCluster(cores, outfile = "")

    registerDoSNOW(cl)

    if(verbose %in% 2:3){
      message("")
      pb<-txtProgressBar(0, nboot, style=3)
      progress<-function(n){
        setTxtProgressBar(pb,n)
      }
      opts <- list(progress=progress)
    } else{
      opts <- list()
    }


    if(!is.null(logFile)) cat("", file= logFile, append=FALSE)

    reslist <- foreach(b = 1:nboot,
                       .packages = c("gtools", "vegan", "LaplacesDemon"),
                       .export = c("calc_association", "sparcc", "cclasso",
                                   "cclasso.sub", "gcoda"),
                       .options.snow = opts) %dopar% {

                         if(!is.null(seed)) set.seed(seeds[b])

                         if(verbose %in% 2:3) progress(b)

                         if(!is.null(logFile)){
                           cat(paste("Iteration", b,"\n"), file=logFile,
                               append=TRUE)
                         }

                         count.tmp <- sapply(1:ncol(countMat), function(i){
                           ind <- sample(1:nrow(countMat), nrow(countMat),
                                         replace = TRUE)
                           countMat[ind, i]
                         })
                         colnames(count.tmp) <- colnames(countMat)

                         assoMat.tmp <- calc_association(count.tmp,
                                                         measure = measure,
                                                         measurePar = measurePar)

                         assoMat.tmp
                       }

    if(verbose %in% 2:3) close(pb)

    stopCluster(cl)
    registerDoSEQ()

    invisible(gc)
    remove(cl)

  } else{   # no parallelization
    reslist <- list()

    if(verbose %in% 2:3){
      message("")
      pb<-txtProgressBar(0, nboot, style=3)
      progress<-function(n){
        setTxtProgressBar(pb,n)
      }
    }

    for(b in 1:nboot){

      if(verbose %in% 2:3){
        #progr <- round(b/nboot*100)
        #message(progr,"%\r",appendLF=FALSE)
        progress(b)
      }

      count.tmp <- sapply(1:ncol(countMat), function(i){
        ind <- sample(1:nrow(countMat), nrow(countMat),
                      replace = TRUE)
        countMat[ind, i]
      })
      colnames(count.tmp) <- colnames(countMat)

      assoMat.tmp <- calc_association(count.tmp, measure = measure,
                                      measurePar = measurePar)

      reslist[[b]] <- assoMat.tmp
    }
  }

  if(verbose %in% 2:3) close(pb)

  assoMat.orig <- assoMat

  # create matrices with logicals if bootstrap associations are at least as
  # extreme as associatins of original data
  sample.greater <- lapply(reslist, function(x) abs(x) >= abs(assoMat.orig))

  # p values are proportion of associations that are greater or equal
  pvals <- (Reduce('+', sample.greater) + 1) / (nboot + 1)
  dimnames(pvals) <- list(colnames(countMat), colnames(countMat))

  return(pvals)
}





