# Argument checking for function createAssoPerm()

#' @keywords internal

.checkArgsCreateAP <- function(args) {

  # Variable for collecting error messages
  errs <- list()
  errs$nerr <- 0
  #errs$msg <- NULL
  
  #-------------------
  # class
  errs <- .checkArg(cond = inherits(args$x, "microNetProps") | 
                      inherits(args$x, "microNet"), 
                    msg = paste0("\"x\" must be of class \"microNet\" ", 
                                 "or \"microNetProps\" (returned by ", 
                                 "netConstruct() or netAnalyze())."), 
                    errs = errs)
                    
  #-------------------
  # computeAsso
  errs <- .checkArg(cond = is.logical(args$computeAsso), 
                    msg = "\"computeAsso\" must be logical", 
                    errs = errs)
  
  #-------------------
  # nPerm
  errs <- .checkArg(cond = (is.numeric(args$nPerm) & args$nPerm >= 0), 
                    msg = "\"nPerm\" must be an integer >= 0", 
                    errs = errs)
  
  args$nPerm <- as.integer(args$nPerm)
  
  #-------------------
  # cores
  errs <- .checkArg(cond = (is.numeric(args$cores) & args$cores >= 0), 
                    msg = "\"cores\" must be an integer >= 0", 
                    errs = errs)
  
  args$cores <- as.integer(args$cores) 
  
  #-------------------
  # seed
  if (!is.null(args$seed)) {
    errs <- .checkArg(cond = is.numeric(args$seed), 
                      msg = "\"seed\" must be numeric (is interpreted as integer)", 
                      errs = errs)
    
    args$seed <- as.integer(args$seed) 
  }

  #-------------------
  # fileStoreAssoPerm
  errs <- .checkArg(cond = is.character(args$fileStoreAssoPerm), 
                    msg = "\"fileStoreAssoPerm\" must be of type character", 
                    errs = errs)
  
  errs <- .checkArg(cond = length(args$fileStoreAssoPerm) == 1, 
                    msg = "\"fileStoreAssoPerm\" must be of length 1", 
                    errs = errs)
  
  #-------------------
  # append
  errs <- .checkArg(cond = is.logical(args$append), 
                    msg = "\"append\" must be logical", 
                    errs = errs)
  
  #-------------------
  # storeCountsPerm
  errs <- .checkArg(cond = is.logical(args$storeCountsPerm), 
                    msg = "\"storeCountsPerm\" must be logical", 
                    errs = errs)
  
  #-------------------
  # fileStoreCountsPerm
  if (args$storeCountsPerm) {
    errs <- .checkArg(cond = is.character(args$fileStoreCountsPerm), 
                      msg = "\"fileStoreCountsPerm\" must be of type character", 
                      errs = errs)
    
    errs <- .checkArg(cond = length(args$fileStoreCountsPerm) == 2, 
                      msg = "\"fileStoreCountsPerm\" must be of length 2", 
                      errs = errs)
  }
  
  #-------------------
  # logFile
  if (!is.null(args$logFile)) {
    errs <- .checkArg(cond = is.character(args$logFile), 
                      msg = "\"logFile\" must be of type character", 
                      errs = errs)
  }
  
  #-------------------
  # verbose
  errs <- .checkArg(cond = is.logical(args$verbose), 
                    msg = "\"verbose\" must be logical", 
                    errs = errs)
  
  #-----------------------------------------------------------------------------
  
  if (errs$nerr > 0) {
    # Get function call of netCompare()
    fn_call <- sys.call(-1)
    fn_call <- utils::capture.output(fn_call)
    
    # Remove white space
    fn_call <- gsub("\\s+", " ", fn_call)
    
    # Enumerate only for more than one errors
    if (errs$nerr == 1) {
      errvec <- paste0(errs$msg, collapse = "\n")
    } else {
      errvec <- paste0(1:errs$nerr, ": ", errs$msg, collapse = "\n")
    }
    
    # Temporarily change allowed length of messages
    tmp <- options("warning.length")
    options(warning.length = 5000L)
    stop("\nin ", fn_call , "\n", errvec, call.=FALSE)
    options(warning.length = tmp$warning.length)
  }
  
  return(args)
}




