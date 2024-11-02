# Argument checking for function plot.diffnet()
#' @keywords internal
.checkArgsPlotDiff<- function(args) {
  
  # Variable for collecting error messages
  errs <- list()
  errs$nerr <- 0
  #errs$msg <- NULL
  
  #-------------------
  # class
  errs <- 
    .checkArg(cond = inherits(args$x, "diffnet"), 
              msg = paste0("\"x\" must be of class \"diffnet\" returned by ", 
                           "diffnet()."), 
              errs = errs)
  
  #-------------------
  # layout
  if (is.null(args$layout)) {
    args$layout <- "spring"
    
  } else {
    layout.tmp <- try(match.fun(args$layout), silent = TRUE)
    
    if (inherits(layout.tmp, "try-error")) {
      if (!is.matrix(args$layout)) {
        errs <- 
          .checkArg(cond = args$layout %in% c("spring", "circle", "groups"), 
                    msg = paste0("\"layout\" must be one of: \"spring\", ", 
                                 "\"circle\", \"groups\"."), 
                    errs = errs)
      }
    } else {
      args$layout <- layout.tmp
    }
  }
  
  #-------------------
  # adjusted
  errs <- .checkArg(cond = is.logical(args$adjusted), 
                    msg = "\"adjusted\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # repulsion
  errs <- .checkArg(cond = is.numeric(args$repulsion) & 
                      length(args$repulsion) == 1, 
                    msg = "\"repulsion\" must be a single numeric value.", 
                    errs = errs)
  
  #-------------------
  # shortenLabels
  choices <- c("intelligent", "simple", "none")
  
  args$shortenLabels <- try(match.arg(args$shortenLabels, choices), 
                            silent = TRUE)
  
  if (inherits(args$shortenLabels, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("shortenLabels", choices), 
                errs = errs)
  }
  
  #-------------------
  # labelLength
  labelLength <- args$labelLength
  
  errs <- .checkArg(cond = is.numeric(labelLength) & length(labelLength) == 1, 
                    msg = "\"labelLength\" must be a single numeric value.", 
                    errs = errs)
  
  if (is.numeric(labelLength) & length(labelLength) == 1) {
    if (!is.integer(labelLength)) {
      if (labelLength %% 1 != 0) {
        warning("Argument \"labelLength\" coerced to integer.")
      }
      labelLength <-  as.integer(labelLength)
    }
  }
  
  args$labelLength <- labelLength
  
  #-------------------
  # labelPattern
  if (!is.null(args$labelPattern)) {
    errs <- 
      .checkArg(cond = is.vector(args$labelPattern) & 
                  length(args$labelPattern) %in% c(3, 5), 
                msg = "\"labelPattern\" must be a vector of length 3 or 5.", 
                errs = errs)
  }
  
  if (!is.null(args$charToRm)) {
    errs <- .checkArg(cond = is.vector(args$charToRm), 
                      msg = "\"charToRm\" must be a vector.", 
                      errs = errs)
  }
  
  #-------------------
  # labelScale
  errs <- .checkArg(cond = is.logical(args$labelScale) & 
                      length(args$labelScale) == 1, 
                    msg = "\"labelScale\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # labelFont
  errs <- .checkArg(cond = is.numeric(args$labelFont) & 
                      length(args$labelFont) == 1, 
                    msg = "\"labelFont\" must be a single numeric value.", 
                    errs = errs)
  
  #-------------------
  # rmSingles
  errs <- .checkArg(cond = is.logical(args$rmSingles) & 
                      length(args$rmSingles) == 1, 
                    msg = "\"rmSingles\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # nodeTransp
  errs <- .checkArg(cond = is.numeric(args$nodeTransp) & 
                      length(args$nodeTransp) == 1, 
                    msg = "\"nodeTransp\" must be a single numeric value.", 
                    errs = errs)
  
  errs <- .checkArg(cond = args$nodeTransp >= 0 & args$nodeTransp <= 100, 
                    msg = "\"nodeTransp\" must be in [0, 100].", 
                    errs = errs)
  
  #-------------------
  # edgeTransp
  errs <- .checkArg(cond = is.numeric(args$edgeTransp) & 
                      length(args$edgeTransp) == 1, 
                    msg = "\"edgeTransp\" must be a single numeric value.", 
                    errs = errs)
  
  errs <- .checkArg(cond = args$edgeTransp >= 0 & args$edgeTransp <= 100, 
                    msg = "\"edgeTransp\" must be in [0, 100].", 
                    errs = errs)
  
  #-------------------
  # borderWidth
  errs <- .checkArg(cond = is.numeric(args$borderWidth) & 
                      length(args$borderWidth) == 1 &
                      args$borderWidth[1] > 0, 
                    msg = "\"borderWidth\" must be a positive numeric value.", 
                    errs = errs)
  
  #-------------------
  # edgeWidth
  errs <- .checkArg(cond = is.numeric(args$edgeWidth) & 
                      length(args$edgeWidth) == 1 &
                      args$edgeWidth[1] >= 0, 
                    msg = paste0("\"edgeWidth\" must be a single ", 
                                 "numeric value >= 0."), 
                    errs = errs)
  
  #-------------------
  # edgeFilter
  choices <- c("none", "highestDiff")
  
  args$edgeFilter <- try(match.arg(args$edgeFilter, choices), 
                         silent = TRUE)
  
  if (inherits(args$edgeFilter, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("edgeFilter", choices), 
                errs = errs)
  }
  
  #-------------------
  # legend
  errs <- .checkArg(cond = is.logical(args$legend) & 
                      length(args$legend) == 1, 
                    msg = "\"legend\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  #legendGroupnames
  if (!is.null(args$legendGroupnames)) {
    errs <- .checkArg(cond = is.character(args$legendGroupnames) & 
                        length(args$legendGroupnames) == 2, 
                      msg = paste0("\"legendGroupnames\" must be a character ", 
                                   "vector of length 2."), 
                      errs = errs)
  }
  
  #-------------------
  #legendTitle
  if (!is.null(args$legendTitle)) {
    errs <- 
      .checkArg(cond = is.character(args$legendTitle) & 
                  length(args$legendTitle) == 1, 
                msg = paste0("\"legendTitle\" must be a single character ", 
                             "vector."), 
                errs = errs)
  }
  
  #-------------------
  #legendArgs
  if (!is.null(args$legendArgs)) {
    errs <- .checkArg(cond = is.list(args$legendArgs), 
                      msg = "\"legendArgs\" must be a list.", 
                      errs = errs)
  }
  #-------------------
  # cex
  errs <- .checkArg(cond = is.numeric(args$cexLegend) & 
                      length(args$cexLegend) == 1, 
                    msg = "\"cexLegend\" must be a single numeric value.", 
                    errs = errs)
  
  errs <- .checkArg(cond = is.numeric(args$cexNodes) & 
                      length(args$cexNodes) == 1, 
                    msg = "\"cexNodes\" must be a single numeric value.", 
                    errs = errs)
  
  errs <- .checkArg(cond = is.numeric(args$cexLabels) & 
                      length(args$cexLabels) == 1, 
                    msg = "\"cexLabels\" must be a single numeric value.", 
                    errs = errs)
  
  errs <- .checkArg(cond = is.numeric(args$cexTitle) & 
                      length(args$cexTitle) == 1, 
                    msg = "\"cexTitle\" must be a single numeric value.", 
                    errs = errs)
  
  
  # mar
  errs <- .checkArg(cond = is.numeric(args$mar) & length(args$mar) == 4, 
                    msg = "\"mar\" must be a numeric vector of length 4.", 
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
