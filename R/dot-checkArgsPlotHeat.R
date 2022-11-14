# Argument checking for function plotHeat()

#' @keywords internal

.checkArgsPlotHeat <- function(args) {
  if (!(is.matrix(args$mat) & is.numeric(args$mat))) {
    stop("\"mat\" must be a numeric matrix.")
  }
  
  if (!is.null(args$pmat)) {
    if (!(is.matrix(args$pmat) & is.numeric(args$pmat))) {
      stop("\"pmat\" must be a numeric matrix.")
    }
    
    if (is.null(dimnames(args$pmat)) || 
        !all.equal(rownames(args$pmat), rownames(args$mat)) ||
        !all.equal(colnames(args$pmat), colnames(args$mat))) {
      stop("Dimnames of \"mat\" and \"pmat\" must be equal.")
    }
  }
  
  args$type <- match.arg(args$type, c("mixed", "full", "lower", "upper"))
  
  args$textUpp <- match.arg(args$textUpp, 
                            c("mat", "sigmat", "pmat", "code", "none"))
  
  args$textLow <- match.arg(args$textLow, 
                            c("mat", "sigmat", "pmat", "code", "none"))
  
  if (args$type %in% c("mixed", "full", "upper")) {
    if (args$textUpp %in% c("pmat", "code", "sigmat") & is.null(args$pmat)) {
      stop("\"pmat\" is missing.")
    }
  }
  
  if (args$type %in% c("mixed", "lower")) {
    if (args$textLow %in% c("pmat", "code", "sigmat") & is.null(args$pmat)) {
      stop("\"pmat\" is missing.")
    }
  }
  
  choices <- c("circle", "square", "ellipse", "number", 
               "shade", "color", "pie")
  args$methUpp <- match.arg(args$methUpp, choices)
  args$methLow <- match.arg(args$methLow, choices)
  
  if (!is.logical(args$diag)) {
    stop("\"diag\" must be logical.")
  }
  
  if (!is.character(args$title)) {
    stop("\"title\" must be a character value.")
  }
  
  if (!(is.numeric(args$mar) & length(args$mar) == 4)) {
    stop("\"mar\" must be a numeric vector with 4 elements.")
  }
  
  args$labPos <- match.arg(args$labPos, c("lt", "ld", "td", "d", "n"))
  
  if (args$labPos == "d") {
    args$diag <- TRUE
  } else if (args$labPos == "ld" && args$type != "lower") {
    stop("type should be \"lower\" if labPos is \"ld\".")
  } else if (args$labPos == "td" && args$type != "upper") {
    stop("type should be \"upper\" if labPos is \"td\".")
  }
  
  if (!is.numeric(args$labCex)) {
    stop("\"labCex\" must be numeric.")
  }
  
  if (!is.numeric(args$textCex)) {
    stop("\"textCex\" must be numeric.")
  }
  
  if (!is.numeric(args$textFont)) {
    stop("\"textFont\" must be numeric.")
  }
  
  if (!is.numeric(args$digits)) {
    stop("\"digits\" must be an integer value.")
  }
  
  args$digits <- as.integer(args$digits)
  
  args$legendPos <- match.arg(args$legendPos, c("r", "b", "n"))
  
  # Sequential and diverging palettes from RColorBrewer function
  sequential <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", 
                  "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
                  "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")
  
  diverging <- c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", 
                 "RdYlGn", "Spectral")
  
  if (!is.null(args$colorPal)) {
    args$colorPal <- match.arg(args$colorPal, c(sequential, diverging))
  }
  
  if (!is.logical(args$addWhite)) {
    stop("\"addWhite\" must be logical.")
  }
  
  if (!is.numeric(args$nCol)) {
    stop("\"nCol\" must be an integer value.")
  }
  
  args$nCol <- as.integer(args$nCol)
  
  if (!is.null(args$colorLim) && 
      (!is.numeric(args$colorLim) | length(args$colorLim) != 2)) {
    stop("\"colorLim\" must be a numeric vector of length 2.")
  }
  
  if (!is.logical(args$revCol)) {
    stop("\"revCol\" must be logical.")
  }
  
  if (!is.null(args$argsUpp) && !is.list(args$argsUpp)) {
    stop("\"argsUpp\" must be a list.")
  }
  
  if (!is.null(args$argsLow) && !is.list(args$argsLow)) {
    stop("\"argsLow\" must be a list.")
  }
  
  return(args)
}