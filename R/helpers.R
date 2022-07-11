#' @keywords internal
.firstUnequalElement <- function(x) {
  m <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  
  duplsum <- apply(m, 2, function(v) {
    sum(duplicated(v) | duplicated(v, fromLast = TRUE))
  })
  
  first_unequal <- which(duplsum < nrow(m))[1]
  
  all_unequal <- which(duplsum == 0)[1]
  
  return(list(first_unequal = first_unequal, all_unequal = all_unequal))
}

#-------------------------------------------------------------------------------
#' @keywords internal
.sigTestRand <- function(randInd, nPermRand, clust1, clust2) {
  randPerm <- numeric(nPermRand)
  
  for (i in 1:nPermRand) {
    clust1.tmp <- gtools::permute(clust1)
    clust2.tmp <- gtools::permute(clust2)
    randPerm[i] <- WGCNA::randIndex(table(clust1.tmp, clust2.tmp), 
                                    adjust = TRUE)
  }
  
  randMean <- mean(randPerm)
  randSD <- sd(randPerm)
  
  normRandPerm <- (randPerm - randMean) / randSD
  normRand <- (randInd[1] - randMean) / randSD
  
  pval <-  (sum(normRandPerm >= abs(normRand)) + 
              sum(normRandPerm <= -abs(normRand))) / nPermRand
  
  return(pval)
}


#-------------------------------------------------------------------------------
#' @keywords internal
# R code from match.arg()
# Used in functions for argument checking
.getMatchArgTxt <- function(arg, choices) {
  sprintf(
    ngettext(
      length(chs <- unique(choices[nzchar(choices)])),
      paste0("\"", arg, "\" should be %s."),
      paste0("\"", arg, "\" should be one of: %s.")
    ),
    paste(dQuote(chs), collapse = ", ")
  )
}

#-------------------------------------------------------------------------------
#' @keywords internal
# Check condition and add error to 'errs' if not fulfilled 
.checkArg <- function(cond, msg, errs) {
  if (!cond) {
    errs$nerr <- errs$nerr + 1
    errs$msg <- c(errs$msg, msg)
  }
  return(errs)
}

#-------------------------------------------------------------------------------
#' @keywords internal
# Compute empirical permutation p-value.
# Used in netCompare()
.calcPermPval <- function(tstat, tstatPerm, nPerm) {
  (sum(tstatPerm >= tstat) + 1) / (nPerm + 1)
}

#-------------------------------------------------------------------------------
#' @keywords internal
# Function for creating significance codes
# Used in summary.microNetComp
.getSigCode <- function(x) {
  if (x <= 0.001) {
    return("***")
  } else if (x <= 0.01) {
    return("** ")
  } else if (x <= 0.05) {
    return("*  ")
  } else if (x <= 0.1) {
    return(".  ")
  } else {
    return(" ")
  }
}
