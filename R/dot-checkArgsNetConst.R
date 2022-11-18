# Argument checking for function netConstruct()

#' @keywords internal

.checkArgsNetConst <- function(args) {
  
  # Variable for collecting error messages
  errs <- list()
  errs$nerr <- 0
  #errs$msg <- NULL
  
  #-------------------
  # dataType
  
  choices <- c("counts", "phyloseq", "correlation", "partialCorr",
               "condDependence", "proportionality", "dissimilarity")
  
  args$dataType <- try(match.arg(args$dataType, choices), silent = TRUE)
  
  if (inherits(args$dataType, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("dataType", choices), 
                errs = errs)
  }
  
  if (args$dataType == "counts") {
    #-------------------
    # measure
    
    choices <- c("pearson",
                 "spearman",
                 "bicor",
                 "sparcc",
                 "cclasso",
                 "ccrepe",
                 "propr",
                 "spieceasi",
                 "spring",
                 "gcoda",
                 "euclidean",
                 "bray",
                 "kld",
                 "ckld",
                 "jeffrey",
                 "jsd",
                 "aitchison")
    
    args$measure <- try(match.arg(args$measure, choices), silent = TRUE)
    
    if (inherits(args$measure, "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("measure", choices), 
                  errs = errs)
    }
    
    #-------------------
    # Set assoType
    args$assoType <-
      if (args$measure %in% c("pearson", "spearman", "bicor", "sparcc",
                              "cclasso", "ccrepe")) {
        "correlation"
      } else if (args$measure %in% c("propr")) {
        "proportionality"
      } else if (args$measure %in% c("spieceasi", "spring", "gcoda")) {
        "condDependence"
      } else if (args$measure %in% c("euclidean", "kld", "jeffrey", "jsd",
                                     "ckld", "aitchison", "bray")) {
        "dissimilarity"
      }
  } else {
    args$measure <- "none"
    args$assoType <- args$dataType
  }
  
  #-------------------
  # Set distNet
  args$distNet <- ifelse(args$assoType == "dissimilarity", TRUE, FALSE)
  
  #-------------------
  # matchDesign
  if (!is.null(args$matchDesign)) {
    errs <- 
      .checkArg(cond = is.numeric(args$matchDesign) & 
                  length(args$matchDesign) == 2, 
                msg = paste0("\"matchDesign\" must a numeric vector ", 
                             "of length 2."), 
                errs = errs)
  }
  
  #-------------------
  # measurePar
  if (!is.null(args$measurePar)) {
    errs <- 
      .checkArg(cond = is.list(args$measurePar), 
                msg = paste0("\"measurePar\" must be a list."), 
                errs = errs)
  }
  
  #-------------------
  # jointPrepro
  if (!is.null(args$jointPrepro)) {
    errs <- 
      .checkArg(cond = is.logical(args$jointPrepro), 
                msg = paste0("\"jointPrepro\" must be logical."), 
                errs = errs)
  }
  
  #-------------------
  # filtTax
  
  choices <- c("totalReads", "relFreq", "numbSamp",
               "highestVar", "highestFreq", "none")
  
  args$filtTax <- try(match.arg(args$filtTax, choices, several.ok = TRUE), 
                      silent = TRUE)
  
  if (inherits(args$filtTax, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("filtTax", choices), 
                errs = errs)
  }
  
  if (args$filtTax[1] != "none") {
    choices <- c( "totalReads", "relFreq", "numbSamp", 
                  "highestVar", "highestFreq")
    
    names(args$filtTaxPar) <- try(match.arg(names(args$filtTaxPar), choices, 
                                            several.ok = TRUE), 
                                  silent = TRUE)
    
    if (inherits(names(args$filtTaxPar), "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("filtTaxPar", choices), 
                  errs = errs)
    }
  }
  
  #-------------------
  # filtSamp
  
  choices <- c("totalReads", "numbTaxa", "highestFreq",
               "none")
  
  args$filtSamp <- try(match.arg(args$filtSamp, choices, several.ok = TRUE), 
                       silent = TRUE)
  
  if (inherits(args$filtSamp, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("filtSamp", choices), 
                errs = errs)
  }
  
  if (args$filtSamp[1] != "none") {
    choices <- c( "totalReads", "numbTaxa", "highestFreq")
    
    names(args$filtSampPar) <- try(match.arg(names(args$filtSampPar), choices, 
                                             several.ok = TRUE), 
                                   silent = TRUE)
    
    if (inherits(names(args$filtSampPar), "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("filtSampPar", choices), 
                  errs = errs)
    }
  }
  
  
  #-------------------
  # Check filtSamp and matchDesign
  
  if ((!"none" %in% args$filtSamp) && !is.null(args$matchDesign)) {
    stop("Filtering samples is not possible for matched subjects.")
    
    errs <- 
      .checkArg(cond = FALSE, 
                msg = paste0("Filtering samples is not possible ", 
                             "for matched subjects."), 
                errs = errs)
  }
  
  #-------------------
  # zeroMethod
  
  choices <- c("none", "pseudo", "pseudoZO", "multRepl", "alrEM", "bayesMult")
  
  args$zeroMethod <- try(match.arg(args$zeroMethod, choices), silent = TRUE)
  
  if (inherits(args$zeroMethod, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("zeroMethod", choices), 
                errs = errs)
  }
  
  if (!is.null(args$zeroPar)) {
    errs <- 
      .checkArg(cond = is.list(args$zeroPar), 
                msg = paste0("\"zeroPar\" must be a list."), 
                errs = errs)
  }
  
  #-------------------
  # normMethod
  
  choices <- c("none", "fractions", "TSS", "CSS", "COM", 
               "rarefy", "VST", "clr", "mclr")
  
  args$normMethod <- try(match.arg(args$normMethod, choices), silent = TRUE)
  
  if (inherits(args$normMethod, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("normMethod", choices), 
                errs = errs)
  }
  
  if (!is.null(args$normPar)) {
    errs <- 
      .checkArg(cond = is.list(args$normPar), 
                msg = paste0("\"normPar\" must be a list."), 
                errs = errs)
  }
  
  #-------------------
  # sparsMethod
  
  choices <- c("none", "t-test", "bootstrap",
               "threshold", "softThreshold", "knn")
  
  args$sparsMethod <- try(match.arg(args$sparsMethod, choices), silent = TRUE)
  
  if (inherits(args$sparsMethod, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("sparsMethod", choices), 
                errs = errs)
  }
  
  #-------------------
  # threshold
  if (args$sparsMethod == "threshold") {
    errs <- 
      .checkArg(cond = is.numeric(args$thresh), 
                msg = paste0("\"thresh\" must be numeric."), 
                errs = errs)
  }
  
  #-------------------
  # alpha
  if (args$sparsMethod %in% c("t-test", "bootstrap")) {
    errs <- 
      .checkArg(cond = is.numeric(args$alpha), 
                msg = paste0("\"alpha\" must be numeric."), 
                errs = errs)
  }
  
  #-------------------
  # adjust
  if (args$sparsMethod %in% c("t-test", "bootstrap")) {
    choices <- c(p.adjust.methods, "lfdr", "adaptBH")
    
    args$adjust <- try(match.arg(args$adjust, choices), silent = TRUE)
    
    if (inherits(args$adjust, "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("adjust", choices), 
                  errs = errs)
    }
  }
  
  #-------------------
  # trueNullMethod
  if (args$adjust == "adaptBH") {
    choices <- c("farco", "lfdr", "mean", "hist", "convest")
    
    args$trueNullMethod <- try(match.arg(args$trueNullMethod, choices), 
                               silent = TRUE)
    
    if (inherits(args$trueNullMethod, "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("trueNullMethod", choices), 
                  errs = errs)
    }
  }
  
  #-------------------
  # lfdrThresh
  if (args$adjust == "lfdr") {
    errs <- 
      .checkArg(cond = is.numeric(args$lfdrThresh), 
                msg = paste0("\"lfdrThresh\" must be numeric."), 
                errs = errs)
  }
  
  if (args$sparsMethod == "bootstrap") {
    #-------------------
    # nboot
    errs <- 
      .checkArg(cond = is.numeric(args$nboot) & args$nboot >= 0, 
                msg = paste0("\"nboot\" must be a non-negative integer."), 
                errs = errs)
    
    args$nboot <- as.integer(args$nboot)
    
    #-------------------
    # assoBoot
    if (!is.null(args$assoBoot)) {
      errs <- 
        .checkArg(cond = is.logical(args$assoBoot) | is.list(args$assoBoot), 
                  msg = paste0("\"assoBoot\" must be either a numeric list or ", 
                               "logical value."), 
                  errs = errs)
    }

    #-------------------
    # cores
    errs <- .checkArg(cond = (is.numeric(args$cores) & args$cores >= 0), 
                      msg = "\"cores\" must be a non-negative integer.", 
                      errs = errs)
    
    args$cores <- as.integer(args$cores)  
    
    args$cores <- min(args$cores, parallel::detectCores())
    
    #-------------------
    # logFile
    if (!is.null(args$logFile)) {
      errs <- .checkArg(cond = is.character(args$logFile), 
                        msg = "\"logFile\" must be a character.", 
                        errs = errs)
    }
  }
  
  
  if (args$sparsMethod == "softThreshold") {
    #-------------------
    # softThreshType
    
    choices <- c("signed", "unsigned", "signed hybrid")
    
    args$softThreshType <- try(match.arg(args$softThreshType, choices), 
                               silent = TRUE)
    
    if (inherits(args$softThreshType, "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("softThreshType", choices), 
                  errs = errs)
    }
    
    #-------------------
    # softThreshPower
    if (!is.null(args$softThreshPower)) {
      errs <- .checkArg(cond = is.numeric(args$softThreshPower), 
                        msg = "\"softThreshPower\" must be numeric.", 
                        errs = errs)
      
      errs <- .checkArg(cond = length(args$softThreshPower) %in% 1:2, 
                        msg = "\"softThreshPower\" must have length 1 or 2.", 
                        errs = errs)
    }
    
    #-------------------
    # softThreshCut
    errs <- .checkArg(cond = is.numeric(args$softThreshCut), 
                      msg = "\"softThreshCut\" must be numeric.", 
                      errs = errs)
    errs <- .checkArg(cond = args$softThreshCut >= 0 & 
                        args$softThreshCut <= 1, 
                      msg = "\"softThreshCut\" must be in [0, 1].", 
                      errs = errs)
  }
  
  
  if (args$sparsMethod == "knn") {
    #-------------------
    # kNeighbor
    errs <- .checkArg(cond = is.numeric(args$kNeighbor), 
                      msg = "\"kNeighbor\" must be an integer.", 
                      errs = errs)
    
    args$kNeighbor <- as.integer(args$kNeighbor) 
    
    #-------------------
    # knnMutual
    errs <- .checkArg(cond = is.logical(args$knnMutual), 
                      msg = "\"knnMutual\" must be logical.", 
                      errs = errs)
  }
  
  #-------------------
  # dissFunc
  
  if (!is.function(args$dissFunc)) {
    choices <- c("signed", "unsigned", "signedPos", "TOMdiss")
    
    args$dissFunc <- try(match.arg(args$dissFunc, choices), 
                         silent = TRUE)
    
    if (inherits(args$dissFunc, "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("dissFunc", choices), 
                  errs = errs)
    }
  }
  
  #-------------------
  # dissFuncPar
  if (!is.null(args$dissFuncPar)) {
    errs <- .checkArg(cond = is.list(args$dissFuncPar), 
                      msg = "\"dissFuncPar\" must be a list.", 
                      errs = errs)
  }
  
  #-------------------
  # simFunc
  if (!is.null(args$simFunc)) {
    errs <- .checkArg(cond = is.function(args$simFunc), 
                      msg = "\"simFunc\" must be a function.", 
                      errs = errs)
  }
  
  #-------------------
  # simFuncPar
  if (!is.null(args$simFuncPar)) {
    errs <- .checkArg(cond = is.list(args$simFuncPar), 
                      msg = "\"simFuncPar\" must be a list.", 
                      errs = errs)
  }
  
  #-------------------
  # scaleDiss
  errs <- .checkArg(cond = is.logical(args$scaleDiss), 
                    msg = "\"scaleDiss\" must be logical.", 
                    errs = errs)
  
  #-------------------
  # weighted
  errs <- .checkArg(cond = is.logical(args$weighted), 
                    msg = "\"weighted\" must be logical.", 
                    errs = errs)
  
  #-------------------
  # verbose
  errs <- .checkArg(cond = is.logical(args$verbose) | 
                      (is.numeric(args$verbose) & args$verbose %in% c(0:3)), 
                    msg = "\"verbose\" must be logical or a value in 0:3.", 
                    errs = errs)
  
  args$verbose <- as.numeric(args$verbose)
  
  #-------------------
  # seed
  
  if (!is.null(args$seed)) {
    errs <- 
      .checkArg(cond = is.numeric(args$seed), 
                msg = "\"seed\" must be numeric (is interpreted as integer)", 
                errs = errs)
    
    args$seed <- as.integer(args$seed) 
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




