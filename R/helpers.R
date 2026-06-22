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
# Turn two consecutive numbers in a vector into one element with 
# double-digit number. Used in editLabels().
.sing2doubDigit <- function(string) {
  
  trues <- suppressWarnings({!is.na(as.numeric(string))})
  
  ptrues <- which(trues)
  
  torm <- NULL
  
  if (sum(trues) > 1) {
    for (i in 1:(length(ptrues) - 1)) {
      if ((ptrues[i] + 1) == ptrues[i + 1]) {
        string[ptrues[i]] <- paste0(string[ptrues[i]], string[ptrues[i + 1]])
        torm <- c(torm, ptrues[i + 1])
      }
    }
    
    string <- string[-torm]
  }

  return(string)
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
# Used in netCompare() and .permTestDiffAsso()
.calcPermPval <- function(tstat, tstatPerm, nPerm, exact,
                          refinePermPvals = FALSE,
                          permApproxParams = list()) {
  pval <- if (exact) {
    sum(tstatPerm >= tstat) / nPerm
  } else {
    (sum(tstatPerm >= tstat) + 1) / (nPerm + 1)
  }

  if (refinePermPvals && !exact) {
    pval <- .permApproxPval(tstat = tstat,
                            tstatPerm = tstatPerm,
                            permApproxParams = permApproxParams)
  }

  pval
}

#-------------------------------------------------------------------------------
#' @keywords internal
# Refine a permutation p-value using permApprox::perm_approx().
.permApproxPval <- function(tstat, tstatPerm, permApproxParams = list()) {
  if (!requireNamespace("permApprox", quietly = TRUE)) {
    stop("Package \"permApprox\" is needed when refinePermPvals = TRUE. ",
         "Install it with installNetCoMiPacks() or install.packages().")
  }

  if (!is.list(permApproxParams)) {
    stop("\"permApproxParams\" must be a list.")
  }

  if (length(permApproxParams) > 0 &&
      (is.null(names(permApproxParams)) || any(names(permApproxParams) == ""))) {
    stop("\"permApproxParams\" must be a named list.")
  }
  
  call_args <- list(
    obs_stats = tstat,
    perm_stats = matrix(tstatPerm, ncol = 1),
    alternative = "greater",
    adjust_method = "none",
    verbose = FALSE
  )
  call_args[names(permApproxParams)] <- permApproxParams

  res <- tryCatch(
    do.call(permApprox::perm_approx, call_args),
    error = function(e) {
      stop("Calling permApprox::perm_approx() failed: ",
           conditionMessage(e), call. = FALSE)
    }
  )

  if (!is.list(res) || is.null(res$p_unadjusted)) {
    stop("permApprox::perm_approx() must return a list containing ",
         "\"p_unadjusted\".")
  }

  pval <- res$p_unadjusted

  if (!is.numeric(pval) || length(pval) != 1 || is.na(pval)) {
    stop("\"p_unadjusted\" returned by permApprox::perm_approx() must be a ",
         "single non-missing numeric p-value.")
  }

  pval
}

#-------------------------------------------------------------------------------
#' @keywords internal
# Function for creating significance codes
# Used in summary.microNetComp
.getSigCode <- function(x, code = c("***", "**", "*", ".", " "), 
                        addSpace = TRUE) {
  
  nc <- max(nchar(code))
  
  if (x > 0.1) {
    out <- code[5]
    
  } else {
    if (x <= 0.001) {
      out <- code[1]
    } else if (x <= 0.01) {
      out <- code[2]
    } else if (x <= 0.05) {
      out <- code[3]
    } else if (x <= 0.1) {
      out <- code[4]
    } 
    
    if (addSpace) {
      while(nchar(out) < nc) {
        out <- paste0(out, " ")
      }
    }
  }
  
  return(out)
}

#-------------------------------------------------------------------------------
#' @keywords internal
# Format p-values for summary output.
.formatPval <- function(x, digits = 6L) {
  out <- rep(NA_character_, length(x))
  is_na <- is.na(x)
  is_small <- !is_na & x <= 10^(-digits)

  out[is_small] <- formatC(x[is_small], format = "e", digits = 2)
  out[!is_na & !is_small] <- formatC(x[!is_na & !is_small],
                                     format = "f",
                                     digits = digits)
  out[is_na] <- NA_character_

  if (length(out) == 1L) out[[1L]] else out
}

#-------------------------------------------------------------------------------
#' @keywords internal
# Make permutation group matrix unique and remove original group vector
# Used in .getPermGroupMat
.cleanPermMat <- function(perm_group_mat, groupvec, exact = FALSE) {
  perm_group_mat <- unique.matrix(perm_group_mat)
  
  if (!exact) {
    # Remove original group vector from permutations (if existent)
    origInPerm <- sapply(1:nrow(perm_group_mat),
                         function(i) {
                           identical(perm_group_mat[i, ], groupvec)
                         })
    
    if (any(origInPerm)) {
      perm_group_mat <- perm_group_mat[!origInPerm, ]
    }
  }
  
  return(perm_group_mat)
}


#-------------------------------------------------------------------------------
#' @keywords internal
# Get the maximum number of permutations for matched set designs
# Used in .getPermGroupMat
.getMaxCombMatch <- function(matchDesign, n) {
  # size of matching sets
  set_size <- sum(matchDesign)
  
  # number of sets
  n_sets <- n / set_size
  
  # possible combinations within each set
  possib <- factorial(matchDesign[1] + matchDesign[2]) / 
    (factorial(matchDesign[1]) * factorial(matchDesign[2]))
  
  # maximum number of permutations
  maxcomb <- possib^n_sets
  
  return(maxcomb)
}


#-------------------------------------------------------------------------------
#' @keywords internal
#' @importFrom rlang as_function
# Based on a function provided by Antoine Fabri (alias moodymudskipper):
# https://stackoverflow.com/a/55182432/697473
# Example:
# sqrt(-1)
# .suppress_warnings(sqrt(-1), startsWith, "N")
# .suppress_warnings(sqrt(-1), startsWith, prefix = "N")

.suppress_warnings <- function(.expr, .f, ...) {
  eval.parent(
    substitute(
      withCallingHandlers( .expr, warning = function (w) {
        cm   <- conditionMessage(w)
        
        cond <- if(is.character(.f)) {
          grepl(.f, cm)
        } else {
          rlang::as_function(.f)(cm, ...)
        }
        
        if (cond) invokeRestart("muffleWarning")   
      })
    )
  )
}




#-------------------------------------------------------------------------------
#' @keywords internal
#' @importFrom graphics rect
# Add rectangles to a GCM heatmap with all orbs
# Used in netAnalyze()
.addOrbRect <- function() {
  graphics::rect(xleft   = c( 0.5,  1.5, 4.5,  7.5),
                 ybottom = c(11.5,  7.5, 4.5,  0.5),
                 xright  = c( 1.5,  4.5, 7.5, 11.5),
                 ytop    = c(10.5, 10.5, 7.5,  4.5),
                 lwd = 2, xpd = NA)
}
