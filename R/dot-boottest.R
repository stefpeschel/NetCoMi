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
#'   \code{countMat}.
#' @param nboot number of bootstrap samples.
#' @param measure character specifying the method used for computing the
#'   associations between taxa.
#' @param measurePar list with parameters passed to the function for computing
#'   associations/dissimilarities. See details for the respective functions.
#' @param cores number of CPU cores used for parallelization.
#' @param logFile character defining a log file, where the number of iteration
#'   is stored. If NULL, no log file is created. wherein the current iteration
#'   numbers are stored.
#' @param verbose logical; if \code{TRUE}, the iteration numbers are printed to
#'   the R console
#' @param seed an optional seed for reproducibility of the results.
#' @param assoBoot list with bootstrap association matrices.
#' @return \tabular{ll}{ \code{pvals}\tab calculated p-values\cr
#'   \code{corrMat}\tab estimated correlation matrix}
#' @references{\insertRef{friedman2012inferring}{NetCoMi}}
#' @keywords internal

.boottest <- function(countMat, assoMat, nboot = 1000, measure, measurePar,
                      cores = 4, logFile = NULL, verbose = TRUE, seed = NULL, 
                      assoBoot = NULL) {
  
  if (!is.null(assoBoot) && !is.list(assoBoot) && assoBoot == TRUE) {
    returnAB <- TRUE
  } else {
    returnAB <- FALSE
  }
  
  if (is.null(assoBoot) || !is.list(assoBoot)) {

    if (!is.null(seed)) {
      seeds <- sample.int(1e8, size = nboot)
    } else {
      seeds <- NULL
    }
    
    if (cores > 1) {
      if (parallel::detectCores() < cores) cores <- parallel::detectCores()
      
      cl <- parallel::makeCluster(cores, outfile = "")
      
      doSNOW::registerDoSNOW(cl)
      
      '%do_or_dopar%' <- get('%dopar%')
      
    } else {
      '%do_or_dopar%' <- get('%do%')
    }
    
    if (verbose %in% 2:3) {
      message("")
      pb <- utils::txtProgressBar(0, nboot, style=3)
      progress <- function(n) {
        utils::setTxtProgressBar(pb,n)
      }
      opts <- list(progress=progress)
    } else {
      opts <- list()
    }
    
    if (!is.null(logFile)) cat("", file = logFile, append = FALSE)
    
    #---------------------------------------------------------------------------
    
    b <- NULL
    
    # Run foreach
    
    assoBoot <- foreach(
      b = 1:nboot,
      .packages = c("gtools", "vegan", "LaplacesDemon"),
      .export = c(
        ".calcAssociation",
        "sparcc",
        "cclasso",
        "cclasso.sub",
        "gcoda"
      ),
      .options.snow = opts
    ) %do_or_dopar% {
      
      if (!is.null(seed)) set.seed(seeds[b])
      
      if (verbose %in% 2:3) progress(b)
      
      if (!is.null(logFile)) {
        cat(paste("Iteration", b,"\n"), file=logFile,
            append=TRUE)
      }
      
      count.tmp <- sapply(1:ncol(countMat), function(i) {
        ind <- sample(1:nrow(countMat), nrow(countMat),
                      replace = TRUE)
        countMat[ind, i]
      })
      colnames(count.tmp) <- colnames(countMat)
      
      assoMat.tmp <- .calcAssociation(countMat = count.tmp,
                                      measure = measure,
                                      measurePar = measurePar,
                                      verbose = 0)
      
      assoMat.tmp
    }
    #---------------------------------------------------------------------------
    
    if (verbose %in% 2:3) close(pb)
    
    if (cores > 1) parallel::stopCluster(cl)
  }
  
  # create matrices with logicals if bootstrap associations are at least as
  # extreme as associatins of original data
  sample.greater <- lapply(assoBoot, function(x) abs(x) >= abs(assoMat))

  # p values are proportion of associations that are greater or equal
  pvals <- (Reduce('+', sample.greater) + 1) / (nboot + 1)
  dimnames(pvals) <- list(colnames(countMat), colnames(countMat))

  if (returnAB) {
    out <- list(pvals = pvals, assoBoot = assoBoot)
  } else {
    out <- list(pvals = pvals)
  }
  
  return(out)
}





