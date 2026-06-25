#' @keywords internal
.zeroTreat <- function(countMat, zeroMethod, zeroParam, 
                       needfrac, needint, verbose) {
  
  if (zeroMethod == "none") {
    countMat_repl <- countMat
    attributes(countMat_repl)$scale <- "counts"
    
  } else if (zeroMethod == "pseudoZO") {
    warning(
      "'zeroMethod = \"pseudoZO\"' is deprecated. Use ",
      "'zeroMethod = \"pseudo\", zeroPar = list(zerosOnly = TRUE)' instead.",
      call. = FALSE
    )
    if (is.null(zeroParam)) zeroParam <- list()
    zeroParam$zerosOnly <- TRUE
    zeroParam$random <- FALSE
    countMat_repl <- .pseudoReplace(countMat, zeroParam, needint, verbose)
    
  } else if (all(countMat != 0)) {
    message("Data contains no zeros.")
    countMat_repl <- countMat
    attributes(countMat_repl)$scale <- "counts"
    
  } else if (zeroMethod =="pseudo") {
    countMat_repl <- .pseudoReplace(countMat, zeroParam, needint, verbose)
    
  } else {
    
    rsums <- Matrix::rowSums(countMat)
    
    if (needfrac || needint) {
      countMat <- t(apply(countMat, 1, function(x) x/sum(x)))
    }
    
    if (zeroMethod == "multRepl") {
      # each zero is replaced by the same small amount (epsilon)
      
      zeroParam$X <- countMat
      if (is.null(zeroParam$label)) zeroParam$label <- 0
      
      if (is.null(zeroParam$dl)) {
        zeroParam$dl <- rep(0.001, ncol(countMat))
      } else if (is.numeric(zeroParam$dl)) {
        zeroParam$dl <- rep(zeroParam$dl, ncol(countMat))
      }
      
      if (is.null(zeroParam$z.warning)) {
        zeroParam$z.warning <- 1
      }
      
      if (verbose %in% 2:3) message("Execute multRepl() ... ", appendLF = FALSE)
      
      countMat_repl <- as.matrix(do.call(zCompositions::multRepl, zeroParam))
      
      if (verbose %in% 2:3) message("Done.")
      
      
    } else if (zeroMethod == "alrEM") {
      
      if (is.null(zeroParam$label)) zeroParam$label <- 0
      
      if (is.null(zeroParam$dl)) {
        dl <- 0.001
        zeroParam$dl <- rep(dl, ncol(countMat))
      } else if (is.numeric(zeroParam$dl)) {
        dl <- zeroParam$dl
        zeroParam$dl <- rep(dl, ncol(countMat))
      } else if (is.vector(zeroParam)) {
        dl <- zeroParam$dl[1]
      }
      
      if (is.null(zeroParam$ini.cov)) zeroParam$ini.cov <- "multRepl"
      
      zeroParam$suppress.print <- ifelse(verbose == 3, FALSE, TRUE)
      
      # number of zeros in each row
      nzeros <- apply(countMat, 1, function(x) sum(x == 0))
      
      if (all(nzeros > 0)) {
        minzeropos <- which(nzeros == min(nzeros))
        
        # multiplicative zero replacement in row with minimum number of zeros
        # to minimize covariance distortion
        rowminz <- countMat[minzeropos, ]
        nz <- sum(rowminz == 0)
        rowminz <- (1 - nz * dl) * rowminz
        rowminz[rowminz == 0] <- dl
        
        zeroParam$X <- rbind(countMat, rowminz)
        
        if (verbose %in% 2:3) message("Execute lrEM() ... ", appendLF = FALSE)
        if (verbose == 3) message("")
        
        countMat_repl <- as.matrix(do.call(zCompositions::lrEM, zeroParam))
        countMat_repl <- countMat_repl[-nrow(countMat_repl), ]
        
        if (verbose %in% 2:3) message("Done.")
        
      } else {
        zeroParam$X <- countMat
        if (verbose %in% 2:3) message("Execute lrEM() ... ", appendLF = FALSE)
        if (verbose == 3) message("")
        countMat_repl <- as.matrix(do.call(zCompositions::lrEM, zeroParam))
        if (verbose %in% 2:3) message("Done.")
      }
      
    } else if (zeroMethod == "bayesMult") {
      zeroParam$output <- "prop"
      if (is.null(zeroParam$method)) zeroParam$method <- "GBM"
      
      if (zeroParam$method == "GBM") {
        # GBM does not work if the data set contains taxa with a positive number 
        # of read counts in only one sample
        nonzerocols <- apply(countMat, 2, function(x) sum(x > 0))
        
        if (any(nonzerocols <= 1)) {
          torm <- which(nonzerocols <= 1)
          zeroParam$X <- countMat[, -torm]
          warning('Taxa with only one positive observation removed
                  (needed for bayesian zero replacement method "GBM".')
          
        } else {
          zeroParam$X <- countMat
        }
        
      } else {
        zeroParam$X <- countMat
      }
      
      zeroParam$suppress.print <- ifelse(verbose == 3, FALSE, TRUE)
      
      if (verbose %in% 2:3) message("Execute cmultRepl() ... ", 
                                    appendLF = FALSE)
      if (verbose == 3) message("")
      
      countMat_repl <- 
        tryCatch(do.call(zCompositions::cmultRepl, zeroParam),
                 error=function(e) {
                   print(paste("Function for zero replacement caused an error:", 
                               e))
                 })
      
      if (verbose %in% 2:3) message("Done.")
      
      countMat_repl <- as.matrix(countMat_repl)
      
    } 
    
    if (needfrac) {
      attributes(countMat_repl)$scale <- "fractions"
      
    } else if (needint) {
      countMat.tmp <- ceiling(countMat_repl * rsums)
      countMat_repl <- apply(countMat.tmp, 2, as.integer)
      rownames(countMat_repl) <- rownames(countMat.tmp)
      
      attributes(countMat_repl)$scale <- "integer"
      
      message("Counts coerced to integer mode (for normalization method). ",
              "Consider using 'zeroMethod = pseudo'.")
      
    } else {
      attributes(countMat_repl)$scale <- "pseudo-counts"
    }
    
    dimnames(countMat_repl) <- dimnames(countMat)
    attributes(countMat_repl)$rowsums <- rsums
  }
  
  
  return(countMat_repl)
}

.pseudoReplace <- function(countMat, zeroParam, needint, verbose) {
  if (is.null(zeroParam)) zeroParam <- list()
  
  if (is.null(zeroParam$pseudocount)) zeroParam$pseudocount <- 1
  if (is.null(zeroParam$zerosOnly)) zeroParam$zerosOnly <- FALSE
  if (is.null(zeroParam$random)) zeroParam$random <- FALSE
  
  if (!is.numeric(zeroParam$pseudocount) ||
      length(zeroParam$pseudocount) != 1 ||
      zeroParam$pseudocount <= 0) {
    stop("'zeroPar$pseudocount' must be a single positive numeric value.")
  }
  
  if (!is.logical(zeroParam$zerosOnly) ||
      length(zeroParam$zerosOnly) != 1 ||
      is.na(zeroParam$zerosOnly)) {
    stop("'zeroPar$zerosOnly' must be a single logical value.")
  }
  
  if (!is.logical(zeroParam$random) ||
      length(zeroParam$random) != 1 ||
      is.na(zeroParam$random)) {
    stop("'zeroPar$random' must be a single logical value.")
  }
  
  if (zeroParam$random && !zeroParam$zerosOnly) {
    stop("'zeroPar$random = TRUE' is only possible if 'zeroPar$zerosOnly = TRUE'.")
  }
  
  if (zeroParam$zerosOnly) {
    countMat_repl <- countMat
    nzeros <- sum(countMat_repl == 0)
    
    if (zeroParam$random) {
      countMat_repl[countMat_repl == 0] <- stats::runif(
        nzeros, min = .Machine$double.eps, max = zeroParam$pseudocount
      )
      
      if (verbose %in% 2:3) {
        message("Zero counts replaced by random values below ",
                zeroParam$pseudocount)
      }
    } else {
      countMat_repl[countMat_repl == 0] <- zeroParam$pseudocount
      
      if (verbose %in% 2:3) {
        message("Zero counts replaced by ", zeroParam$pseudocount)
      }
    }
  } else {
    countMat_repl <- countMat + zeroParam$pseudocount
    
    if (verbose %in% 2:3) {
      message("Pseudo count of ", zeroParam$pseudocount, " added.")
    }
  }
  
  if (needint) {
    if (zeroParam$pseudocount != 1) {
      message("Counts coerced to integer mode (needed for normalization ", 
              "method). \nConsider using unit pseudo counts for zero ", 
              "treatment: \n'zeroMethod = pseudo', ", 
              "zeroPar = list(pseudocount = 1)")
    }
    countMat.tmp <- ceiling(countMat_repl)
    countMat_repl <- apply(countMat.tmp, 2, as.integer)
    rownames(countMat_repl) <- rownames(countMat.tmp)
    
    attributes(countMat_repl)$scale <- "integer"
    
  } else {
    attributes(countMat_repl)$scale <- "pseudo-counts"
  }
  
  return(countMat_repl)
}



