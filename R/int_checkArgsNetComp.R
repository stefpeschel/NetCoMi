# Argument checking for function netCompare()

#' @keywords internal

.checkArgsNetComp <- function(args) {

  # Variable for collecting error messages
  errs <- list()
  errs$nerr <- 0
  #errs$msg <- NULL

  errs <- .checkArg(cond = "microNetProps" %in% class(args$x), 
                   msg = paste0("'x' must be of class 'microNetProps' ", 
                                "(returned by netAnalyze())"), 
                   errs = errs)

  errs <- .checkArg(cond = is.logical(args$permTest), 
                   msg = "'permTest' must be logical", 
                   errs = errs)
  
  errs <- .checkArg(cond = is.logical(args$gcd), 
                   msg = "'gcd' must be logical", 
                   errs = errs)

  errs <- .checkArg(cond = is.logical(args$verbose), 
                   msg = "'verbose' must be logical", 
                   errs = errs)

  errs <- .checkArg(cond = is.logical(args$storeAssoPerm), 
                    msg = "'storeAssoPerm' must be logical", 
                    errs = errs)

  errs <- .checkArg(cond = is.logical(args$storeCountsPerm), 
                   msg = "'storeCountsPerm' must be logical", 
                   errs = errs)

  errs <- .checkArg(cond = is.logical(args$returnPermProps), 
                   msg = "'returnPermProps' must be logical", 
                   errs = errs)

  errs <- .checkArg(cond = is.logical(args$returnPermCentr), 
                   msg = "'returnPermCentr' must be logical", 
                   errs = errs)

  if (!is.null(args$lnormFit)) {
    errs <- .checkArg(cond = is.logical(args$lnormFit), 
                     msg = "'lnormFit' must be logical", 
                     errs = errs)
  }

  errs <- .checkArg(cond = (args$jaccQuant >= 0 & args$jaccQuant <= 1), 
                   msg = "'jaccQuant' must be in [0,1]", 
                   errs = errs)

  errs <- .checkArg(cond = (is.numeric(args$nPerm) & args$nPerm >= 0), 
                   msg = "'nPerm' must be an integer >= 0", 
                   errs = errs)
  
  args$nPerm <- as.integer(args$nPerm)
  
  args$adjust <- match.arg(args$adjust, c(p.adjust.methods, "lfdr", "adaptBH"))

  args$trueNullMethod <- match.arg(args$trueNullMethod, c("convest", "lfdr", 
                                                          "mean", "hist", 
                                                          "farco"))

  errs <- .checkArg(cond = (is.numeric(args$nPermRand) & args$nPermRand >= 0), 
                   msg = "'nPermRand' must be an integer >= 0", 
                   errs = errs)
  
  args$nPermRand <- as.integer(args$nPermRand)

  errs <- .checkArg(cond = (is.numeric(args$cores) & args$cores >= 0), 
                   msg = "'cores' must be an integer >= 0", 
                   errs = errs)
  
  args$cores <- as.integer(args$cores)  

  if (!is.null(args$logFile)) {
    errs <- .checkArg(cond = is.character(args$logFile), 
                     msg = "'logFile' must be of type character", 
                     errs = errs)
  }

  
  if (args$permTest) {
    errs <- .checkArg(cond = args$x$paramsNetConstruct$dataType == "counts", 
                      msg = paste0("Permutation tests only possible if count ", 
                                   "tables were used for network construction."), 
                      errs = errs)
    
    if (!is.null(args$fileLoadAssoPerm)) {
      errs <- .checkArg(cond = is.character(args$fileLoadAssoPerm), 
                       msg = "'fileLoadAssoPerm' must be of type character", 
                       errs = errs)

      errs <- .checkArg(cond = length(args$fileLoadAssoPerm) == 1, 
                       msg = "'fileLoadAssoPerm' must be of length 1", 
                       errs = errs)
    }
    
    if (!is.null(args$fileLoadCountsPerm)) {
      errs <- .checkArg(cond = is.character(args$fileLoadCountsPerm), 
                       msg = "'fileLoadCountsPerm' must be of type character", 
                       errs = errs)
      
      if (!length(args$fileLoadCountsPerm) == 2) {
        nerr <- nerr + 1
        errmsg <- c(errmsg, "'fileLoadCountsPerm' must be of length 2")
      }
      
      errs <- .checkArg(cond = length(args$fileLoadCountsPerm) == 2, 
                       msg = "'fileLoadCountsPerm' must be of length 2", 
                       errs = errs)
    }
    
    if (args$storeCountsPerm) {
      errs <- .checkArg(cond = is.character(args$fileStoreCountsPerm), 
                       msg = "'fileStoreCountsPerm' must be of type character", 
                       errs = errs)

      errs <- .checkArg(cond = length(args$fileStoreCountsPerm) == 2, 
                       msg = "'fileStoreCountsPerm' must be of length 2", 
                       errs = errs)
    }
    
    if (args$storeAssoPerm) {
      errs <- .checkArg(cond = is.character(args$fileStoreAssoPerm), 
                       msg = "'fileStoreAssoPerm' must be of type character", 
                       errs = errs)
      
      errs <- .checkArg(cond = length(args$fileStoreAssoPerm) == 1, 
                       msg = "'fileStoreAssoPerm' must be of length 1", 
                       errs = errs)
    }
    
    if (!is.null(args$assoPerm)) {
      errs <- .checkArg(cond = (is.list(args$assoPerm) & 
                                 length(args$assoPerm) == 2), 
                       msg = "'assoPerm' must be a list of length 2", 
                       errs = errs)
    }
      
    
    if (!is.null(args$dissPerm)) {
      errs <- .checkArg(cond = (is.list(args$dissPerm) & 
                                 length(args$dissPerm) == 2), 
                       msg = "'dissPerm' must be a list of length 2", 
                       errs = errs)
    }
  }

  if (errs$nerr > 0) {
    # Get function call of netCompare()
    fn_call <- sys.call(-1)
    fn_call <- utils::capture.output(fn_call)
    
    # Remove white space
    fn_call <- gsub("\\s+", " ", fn_call)
    
    # Enumerate only for more than one errors
    if (errs$nerr == 1) {
      errvec <- paste0(errs$msg, collapse = "\n")
    } else{
      errvec <- paste0(1:errs$nerr, ": ", errs$msg, collapse = "\n")
    }
    
    stop("\nin ", fn_call , "\n", errvec, call.=FALSE)
  }

  return(args)
}




