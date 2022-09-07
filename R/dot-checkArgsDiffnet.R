# Argument checking for function diffnet()

#' @keywords internal

.checkArgsDiffnet <- function(args) {
  
  # Variable for collecting error messages
  errs <- list()
  errs$nerr <- 0
  #errs$msg <- NULL
  
  #-------------------
  # x
  errs <- 
    .checkArg(cond = inherits(args$x, "microNet"), 
              msg = paste0("\"x\" must be of class \"microNet\" returned by ", 
                           "netConstruct()."), 
              errs = errs)
  
  if (args$x$assoType == "dissimilarity") {
    stop("Differential network not implemented for dissimilarity-based ", 
         "networks.")
  }
  
  if (is.null(args$x$assoMat2)) {
    stop("\"x\" is a single network. Differential network cannot be computed.")
  }
  
  #-------------------
  # diffMethod
  choices <- c("discordant", "permute", "fisherTest")
  
  args$diffMethod <- try(match.arg(args$diffMethod, choices), silent = TRUE)
  
  if (inherits(args$diffMethod, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("diffMethod", choices), 
                errs = errs)
  }
  
  if (args$diffMethod == "permute" && args$x$parameters$dataType != "counts") {
    stop(paste0("Permutation test only possible if count matrices were used ", 
                "for network construction. Use Fisher test instead."))
  }
  
  if (args$diffMethod == "discordant" && 
      args$x$parameters$dataType != "counts") {
    stop(paste0("Discordant method only possible if count matrices were used ", 
                "for network construction. Use Fisher test instead."))
  }
  
  #-------------------
  # discordThresh
  errs <- 
    .checkArg(cond = is.numeric(args$discordThresh) & 
                length(args$discordThresh) == 1 & 
                args$discordThresh[1] >= 0 & args$discordThresh[1] <= 1, 
              msg = "\"discordThresh\" must be a numeric value in [0, 1].", 
              errs = errs)
  
  #-------------------
  # n1, n2
  if (args$diffMethod == "fisherTest" && 
      args$x$parameters$dataType != "counts") {
    errs <- 
      .checkArg(is.numeric(args$n1) & length(args$n1) == 1 &
                  is.numeric(args$n2) & length(args$n2) == 1,
                msg = paste0("\"n1\" and \"n2\" must be single numeric values ", 
                             "(needed because count matrices are missing)."), 
                errs = errs)
  }

  if (!is.null(args$n1)) {
    args$n1 <- as.integer(args$n1)
  }
  
  if (!is.null(args$n2)) {
    args$n2 <- as.integer(args$n2)
  }
  
  #-------------------
  # fisherTrans
  errs <- .checkArg(cond = is.logical(args$fisherTrans), 
                    msg = "\"fisherTrans\" must be logical.", 
                    errs = errs)
  #-------------------
  # nPerm
  errs <- .checkArg(cond = is.numeric(args$nPerm) & length(args$nPerm) == 1 &
                      args$nPerm[1] >= 0, 
                    msg = paste0("\"nPerm\" must be a single numeric value ", 
                                 "(is interpreted as integer)."), 
                    errs = errs)
  
  args$nPerm <- as.integer(args$nPerm) 
  
  #-------------------
  # permPvalsMethod
  if (args$permPvalsMethod != "pseudo") args$permPvalsMethod <- "pseudo"
  
  #-------------------
  # cores
  errs <- .checkArg(cond = (is.numeric(args$cores) & args$cores >= 0), 
                    msg = "\"cores\" must be an integer >= 0.", 
                    errs = errs)
  
  args$cores <- as.integer(args$cores)  
  
  #-------------------
  # verbose
  if (args$verbose %in% c(0,1)) {
    args$verbose <- as.logical(args$verbose)
    
  } else {
    errs <- 
      .checkArg(cond = is.logical(args$verbose) & length(args$verbose) == 1, 
                msg = "\"verbose\" must be logical.", 
                errs = errs)
  }
  
  #-------------------
  # logFile
  if (!is.null(args$logFile)) {
    errs <- 
      .checkArg(cond = is.character(args$logFile) & length(args$logFile) == 1, 
                msg = "\"logFile\" must be logical.", 
                errs = errs)
  }
  
  #-------------------
  # seed
  if (!is.null(args$seed)) {
    errs <- .checkArg(cond = is.numeric(args$seed), 
                      msg = "\"seed\" must be an integer >= 0.", 
                      errs = errs)
    
    args$seed <- as.integer(args$seed) 
  }
  
  #-------------------
  # alpha
  errs <- 
    .checkArg(cond = is.numeric(args$alpha) & length(args$alpha) == 1 & 
                args$alpha[1] >= 0, 
              msg = "\"alpha\" must be a numeric value in [0, 1].", 
              errs = errs)
  
  if (is.numeric(args$alpha) & length(args$alpha) == 1 & args$alpha > 1) {
    args$alpha <- args$alpha / 100
    message("\"alpha\" transformed to ", args$alpha, ".")
  }
  
  #-------------------
  # adjust
  choices <- c(p.adjust.methods, "lfdr", "adaptBH")
  
  args$adjust <- try(match.arg(args$adjust, choices), silent = TRUE)
  
  if (inherits(args$adjust, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("adjust", choices), 
                errs = errs)
  }
  
  #-------------------
  # lfdrThresh
  errs <- 
    .checkArg(cond = is.numeric(args$lfdrThresh) & 
                length(args$lfdrThresh) == 1 & 
                args$lfdrThresh[1] >= 0 & args$lfdrThresh[1] <= 1, 
              msg = "\"lfdrThresh\" must be a numeric value in [0, 1].", 
              errs = errs)
  
  #-------------------
  # trueNullMethod
  choices <- c("convest", "lfdr", "mean", "hist", "farco")
  
  args$trueNullMethod <- try(match.arg(args$trueNullMethod, choices), 
                             silent = TRUE)
  
  if (inherits(args$trueNullMethod, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("trueNullMethod", choices), 
                errs = errs)
  }
  
  #-------------------
  # pvalsVec
  
  #-------------------
  # fileLoadAssoPerm
  if (!is.null(args$fileLoadAssoPerm)) {
    errs <- .checkArg(cond = is.character(args$fileLoadAssoPerm), 
                      msg = "\"fileLoadAssoPerm\" must be of type character.", 
                      errs = errs)
    
    errs <- .checkArg(cond = length(args$fileLoadAssoPerm) == 1, 
                      msg = "\"fileLoadAssoPerm\" must be of length 1.", 
                      errs = errs)
  }
  
  #-------------------
  # fileLoadCountsPerm
  
  if (!is.null(args$fileLoadCountsPerm)) {
    errs <- .checkArg(cond = is.character(args$fileLoadCountsPerm), 
                      msg = "\"fileLoadCountsPerm\" must be of type character.", 
                      errs = errs)
    
    errs <- .checkArg(cond = length(args$fileLoadCountsPerm) == 2, 
                      msg = "\"fileLoadCountsPerm\" must be of length 2.", 
                      errs = errs)
  }
  
  #-------------------
  # storeAssoPerm
  errs <- .checkArg(cond = is.logical(args$storeAssoPerm), 
                    msg = "\"storeAssoPerm\" must be logical.", 
                    errs = errs)
  
  #-------------------
  # fileStoreAssoPerm
  if (args$storeAssoPerm) {
    errs <- .checkArg(cond = is.character(args$fileStoreAssoPerm), 
                      msg = "\"fileStoreAssoPerm\" must be of type character.", 
                      errs = errs)
    
    errs <- .checkArg(cond = length(args$fileStoreAssoPerm) == 1, 
                      msg = "\"fileStoreAssoPerm\" must be of length 1.", 
                      errs = errs)
  }
  
  #-------------------
  # storeCountsPerm
  errs <- .checkArg(cond = is.logical(args$storeCountsPerm), 
                    msg = "\"storeCountsPerm\" must be logical.", 
                    errs = errs)
  
  #-------------------
  # fileStoreCountsPerm
  if (args$storeCountsPerm) {
    errs <- .checkArg(cond = is.character(args$fileStoreCountsPerm), 
                      msg = "\"fileStoreCountsPerm\" must be of type character.", 
                      errs = errs)
    
    errs <- .checkArg(cond = length(args$fileStoreCountsPerm) == 2, 
                      msg = "\"fileStoreCountsPerm\" must be of length 2.", 
                      errs = errs)
  }
  
  #-----------------------------------------------------------------------------
  # Install limma package
  
  if (args$diffMethod != "discordant" && args$adjust == "adaptBH" && 
      !requireNamespace("limma", quietly = TRUE)) {
    
    message("Installing missing package \"limma\" ...")
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    
    BiocManager::install("limma", dependencies = TRUE)
    message("Done.")
    
    message("Check whether installed package can be loaded ...")
    requireNamespace("limma")
    message("Done.")
  }
  
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
