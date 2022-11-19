# Argument checking for function netAnalyze()

#' @keywords internal

.checkArgsNetAna <- function(args) {
  
  # Variable for collecting error messages
  errs <- list()
  errs$nerr <- 0
  #errs$msg <- NULL
  
  #-------------------
  # class
  errs <- 
    .checkArg(cond = inherits(args$net, "microNet"), 
              msg = paste0("\"net\" must be of class \"microNet\" returned by ", 
                           "netConstruct()."), 
              errs = errs)
  
  #-------------------
  # centrLCC
  errs <- 
    .checkArg(cond = is.logical(args$centrLCC), 
              msg = "\"centrLCC\" must be logical.", 
              errs = errs)
  
  #-------------------
  # avDissIgnoreInf
  errs <- 
    .checkArg(cond = is.logical(args$avDissIgnoreInf), 
              msg = "\"avDissIgnoreInf\" must be logical.", 
              errs = errs)
  
  #-------------------
  # sPathAlgo
  choices <- c("automatic", "unweighted", "dijkstra", "bellman-ford",  
               "johnson")
  
  args$sPathAlgo <- try(match.arg(args$sPathAlgo, choices), silent = TRUE)
  
  if (inherits(args$sPathAlgo, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("sPathAlgo", choices), 
                errs = errs)
  }
  
  #-------------------
  # sPathNorm
  errs <- 
    .checkArg(cond = is.logical(args$sPathNorm), 
              msg = "\"sPathNorm\" must be logical.", 
              errs = errs)
  
  #-------------------
  # normNatConnect
  errs <- 
    .checkArg(cond = is.logical(args$normNatConnect), 
              msg = "\"normNatConnect\" must be logical.", 
              errs = errs)
  
  #-------------------
  # clustMethod
  
  if (is.null(args$clustMethod)) {
    if (args$net$assoType == "dissimilarity") {
      args$clustMethod <- "hierarchical"
    } else {
      args$clustMethod <- "cluster_fast_greedy"
    }
    
  } else {
    choices <- c("none", 
                 "hierarchical",
                 "cluster_edge_betweenness",
                 "cluster_fast_greedy",
                 "cluster_leading_eigen",
                 "cluster_louvain",
                 "cluster_optimal",
                 "cluster_spinglass",
                 "cluster_walktrap")
    
    args$clustMethod <- try(match.arg(args$clustMethod, choices), silent = TRUE)
    
    if (inherits(args$clustMethod, "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("clustMethod", choices), 
                  errs = errs)
    }
  }
  
  #-------------------
  # clustPar
  if (!is.null(args$clustPar)) {
    errs <- 
      .checkArg(cond = is.list(args$clustPar), 
                msg = "\"clustPar\" must be a list", 
                errs = errs)
  }

  #-------------------
  # clustPar2
  if (is.null(args$clustPar2)) {
    args$clustPar2 <- args$clustPar
    
  } else {
    errs <- 
      .checkArg(cond = is.list(args$clustPar2), 
                msg = "\"clustPar2\" must be a list", 
                errs = errs)
  }
  
  #-------------------
  # weightClustCoef
  errs <- 
    .checkArg(cond = is.logical(args$weightClustCoef), 
              msg = "\"weightClustCoef\" must be logical.", 
              errs = errs)
  
  #-------------------
  # hubPar
  choices <- c("degree", "betweenness", "closeness", "eigenvector")
  
  args$hubPar <- try(match.arg(args$hubPar, choices, several.ok = TRUE), 
                     silent = TRUE)
  
  if (inherits(args$hubPar, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("hubPar", choices), 
                errs = errs)
  }

  #-------------------
  # hubQuant
  errs <- 
    .checkArg(cond = is.numeric(args$hubQuant) & length(args$hubQuant) == 1 & 
                args$hubQuant[1] >= 0 & args$hubQuant[1] <= 1, 
              msg = "\"hubQuant\" must be a single numeric value in [0,1]", 
              errs = errs)

  #-------------------
  # lnormFit
  errs <- 
    .checkArg(cond = is.logical(args$lnormFit), 
              msg = "\"lnormFit\" must be logical.", 
              errs = errs)
  
  #-------------------
  # weightDeg
  errs <- 
    .checkArg(cond = is.logical(args$weightDeg), 
              msg = "\"weightDeg\" must be logical.", 
              errs = errs)
  
  #-------------------
  # normDeg
  errs <- 
    .checkArg(cond = is.logical(args$normDeg), 
              msg = "\"normDeg\" must be logical.", 
              errs = errs)
  
  #-------------------
  # normBetw
  errs <- 
    .checkArg(cond = is.logical(args$normBetw), 
              msg = "\"normBetw\" must be logical.", 
              errs = errs)
  
  #-------------------
  # normClose
  errs <- 
    .checkArg(cond = is.logical(args$normClose), 
              msg = "\"normClose\" must be logical.", 
              errs = errs)
  
  #-------------------
  # normEigen
  errs <- 
    .checkArg(cond = is.logical(args$normEigen), 
              msg = "\"normEigen\" must be logical.", 
              errs = errs)
  
  #-------------------
  # connectivity
  errs <- 
    .checkArg(cond = is.logical(args$connectivity), 
              msg = "\"connectivity\" must be logical.", 
              errs = errs)
  
  #-------------------
  # graphlet
  errs <- 
    .checkArg(cond = is.logical(args$graphlet), 
              msg = "\"graphlet\" must be logical.", 
              errs = errs)
  
  if (args$graphlet) {
    #-------------------
    # orbits
    errs <- 
      .checkArg(cond = is.numeric(args$orbits), 
                msg = "\"orbits\" vector must be numeric.", 
                errs = errs)
    
    errs <- 
      .checkArg(length(args$orbits) >= 2 & 
                  length(args$orbits) <= 15 &
                  all(args$orbits %in% 0:14), 
                msg = "Only orbits 0 to 14 (from 4-node graphlets) are allowed.", 
                errs = errs)
    
    #-------------------
    # gcmHeat
    errs <- 
      .checkArg(cond = is.logical(args$gcmHeat), 
                msg = "\"gcmHeat\" must be logical.", 
                errs = errs)
    
    #-------------------
    # gcmHeatLCC
    errs <- 
      .checkArg(cond = is.logical(args$gcmHeatLCC), 
                msg = "\"gcmHeatLCC\" must be logical.", 
                errs = errs)
  }


  #-------------------
  # verbose
  errs <- .checkArg(cond = is.logical(args$verbose) | 
                      (is.numeric(args$verbose) & args$verbose %in% c(0:2)), 
                    msg = "\"verbose\" must be logical or a value in 0:2.", 
                    errs = errs)
  
  args$verbose <- as.numeric(args$verbose)
  
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




