# Argument checking for function plot.microNetProps()

#' @keywords internal

.checkArgsPlotMNP <- function(args) {
  
  # Variable for collecting error messages
  errs <- list()
  errs$nerr <- 0
  #errs$msg <- NULL
  
  #-------------------
  # class
  if (!inherits(args$x, "microNetProps")) {
    stop("\"x\" must be of class \"microNetProps\" returned by netAnalyze().")
  }
  
  twoNets <- args$x$input$twoNets
  
  #-------------------
  # layout
  
  layout <- args$layout
  
  if (is.null(layout)) {
    layout <- "spring"
    
  } else {
    
    layout.tmp <- try(match.fun(layout), silent = TRUE)
    
    if (inherits(layout.tmp, "try-error")) {
      if (is.matrix(layout)) {
        if (!any(rownames(layout) %in% rownames(args$x$input$adjaMat1))) {
          warning("Rownames of layout matrix do not match node names.")
        }
      } else if (is.list(layout)) {
        errs <- .checkArg(cond = length(layout) == 2, 
                          msg = paste0("\"layout\" must have length 2."), 
                          errs = errs)
        
        if (!any(rownames(layout[[1]]) %in% rownames(args$x$input$adjaMat1))) {
          warning("Rownames of layout matrix do not match node names.")
        }
        
        if (!is.null(layout[[2]])) {
          if (!any(rownames(layout[[2]]) %in% rownames(args$x$input$adjaMat2))) {
            warning("Rownames of layout matrix do not match node names.")
          }
        }
        
      } else {
        errs <- .checkArg(cond = layout %in% c("spring", "circle", "groups"), 
                          msg = paste0("\"layout\" must be one of: \"spring\", ", 
                                       "\"circle\", \"groups\"."), 
                          errs = errs)
      }
    } else {
      layout <- layout.tmp
    }
  }
  args$layout <- layout
  
  #-------------------
  # sameLayout
  errs <- .checkArg(cond = is.logical(args$sameLayout), 
                    msg = "\"sameLayout\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # layoutGroup
  if (args$sameLayout & twoNets) {
    errs <- .checkArg(cond = args$layoutGroup %in% c(1,2, "union"), 
                      msg = paste0("\"layoutGroup\" must be one of: 1, 2, ",
                                   "\"union\"."), 
                      errs = errs)
  }
  
  #-------------------
  # repulsion
  errs <- .checkArg(cond = is.numeric(args$repulsion) & 
                      length(args$repulsion) == 1, 
                    msg = "\"repulsion\" must be a single numeric value.", 
                    errs = errs)
  
  #-------------------
  # groupNames
  if (!is.null(args$groupNames)) {
    errs <- 
      .checkArg(cond = is.vector(args$groupNames) & 
                  length(args$groupNames) == 2, 
                msg = "\"groupNames\" must be a vector with two elements.", 
                errs = errs)
  }
  
  #-------------------
  # groupsChanged
  errs <- .checkArg(cond = is.logical(args$groupsChanged) & 
                      length(args$groupsChanged) == 1, 
                    msg = "\"groupsChanged\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # labels
  
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
  
  #-------------------
  # charToRm
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
  # labelFile
  if (!is.null(args$labelFile)) {
    errs <- .checkArg(cond = is.character(args$labelFile) & 
                        length(args$labelFile) == 1, 
                      msg = "\"labelFile\" must be a single character value.", 
                      errs = errs)
  }
  
  #-------------------
  # nodeFilter
  choices <- c("none", "highestConnect", "highestDegree", "highestBetween",
               "highestClose", "highestEigen", "clustTaxon", "clustMin", 
               "names")
  
  args$nodeFilter <- try(match.arg(args$nodeFilter, choices), 
                         silent = TRUE)
  
  if (inherits(args$nodeFilter, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("nodeFilter", choices), 
                errs = errs)
  }
  
  #-------------------
  # nodeFilterPar
  
  if (args$nodeFilter %in% c("highestConnect", "clustMin",
                             "highestDegree", "highestBetween",
                             "highestClose", "highestEigen")) {
    
    errs <- .checkArg(cond = is.numeric(args$nodeFilterPar) & 
                        length(args$nodeFilterPar) == 1, 
                      msg = "\"nodeFilterPar\" must be a single numeric value.", 
                      errs = errs)
    
  } else if (args$nodeFilter == "clustTaxon") {
    errs <- .checkArg(cond = is.vector(args$nodeFilterPar) & 
                        is.character(args$nodeFilterPar), 
                      msg = "\"nodeFilterPar\" must be a character vector.", 
                      errs = errs)
  }
  
  #-------------------
  # rmSingles
  if (args$rmSingles[1] == TRUE) {
    args$rmSingles <- "all"
    
  } else if (args$rmSingles[1] == FALSE) {
    args$rmSingles <- "none"
    
  } else {
    choices <- c("all", "inboth", "none")
    
    args$rmSingles <- try(match.arg(args$rmSingles, choices), 
                          silent = TRUE)
    
    if (inherits(args$rmSingles, "try-error")) {
      errs <- 
        .checkArg(cond = FALSE, 
                  msg = .getMatchArgTxt("rmSingles", choices), 
                  errs = errs)
    }
  }
  
  #-------------------
  # nodeSize
  choices <- c("fix", "hubs", "degree", "betweenness", "closeness",
               "eigenvector", "counts", "normCounts", "TSS", "fractions", 
               "CSS", "COM", "rarefy", "VST", "clr", "mclr")
  
  args$nodeSize <- try(match.arg(args$nodeSize, choices), 
                       silent = TRUE)
  
  if (inherits(args$nodeSize, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("nodeSize", choices), 
                errs = errs)
  }
  
  #-------------------
  # normPar
  if (!is.null(args$normPar)) {
    errs <- .checkArg(cond = is.list(args$normPar), 
                      msg = "\"normPar\" must be a list.", 
                      errs = errs)
  }
  
  #-------------------
  # nodeSizeSpread
  errs <- .checkArg(cond = is.numeric(args$nodeSizeSpread) & 
                      length(args$nodeSizeSpread) == 1 &
                      args$nodeSizeSpread[1] >= 0, 
                    msg = "\"nodeSizeSpread\" must be a numeric value >= 0.", 
                    errs = errs)
  
  #-------------------
  # nodeColor
  if (!is.null(args$nodeColor)) {
    nodeColor.tmp <- try(match.arg(args$nodeColor,
                                   choices = c("cluster", "feature", 
                                               "colorVec")),
                         silent = TRUE)
    
    if (inherits(nodeColor.tmp, "try-error")) {
      if (!is.character(args$nodeColor)) {
        errs <- 
          .checkArg(cond = FALSE, 
                    msg = paste0("Possible values for \"nodeColor\" are: \n", 
                                 "\"cluster\", \"feature\", \"colorVec\", ", 
                                 "or a character defining a color."), 
                    errs = errs)
      }
    } else {
      args$nodeColor <- nodeColor.tmp
    }
  }
  
  #-------------------
  # colorVec
  
  #-------------------
  # featVecCol
  if (!is.null(args$featVecCol)) {
    errs <- .checkArg(cond = is.vector(args$featVecCol) | 
                        is.factor(args$featVecCol), 
                      msg = "\"featVecCol\" must be a vector.", 
                      errs = errs)
  }
  
  #-------------------
  # sameFeatCol
  errs <- .checkArg(cond = is.logical(args$sameFeatCol) & 
                      length(args$sameFeatCol) == 1, 
                    msg = "\"sameFeatCol\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # sameClustCol
  errs <- .checkArg(cond = is.logical(args$sameClustCol) & 
                      length(args$sameClustCol) == 1, 
                    msg = "\"sameClustCol\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # sameColThresh
  sameColThresh <- args$sameColThresh
  
  errs <- .checkArg(cond = is.numeric(sameColThresh) & 
                      length(sameColThresh) == 1, 
                    msg = "\"sameColThresh\" must be a single numeric value.", 
                    errs = errs)
  
  if (!is.integer(sameColThresh)) {
    if (sameColThresh %% 1 != 0) {
      warning("Argument \"sameColThresh\" coerced to integer.")
    }
    sameColThresh <-  as.integer(sameColThresh)
  }
  args$sameColThresh <- sameColThresh
  
  #-------------------
  # nodeShape
  errs <- 
    .checkArg(cond = all(args$nodeShape %in% c("circle", "square", 
                                               "triangle", "diamond")), 
              msg = paste0("\"nodeShape\" must have values in: ", 
                           "\"circle\", \"square\", \"triangle\", \"diamond\""), 
              errs = errs)
  
  #-------------------
  # featVecShape
  if (!is.null(args$featVecShape)) {
    errs <- 
      .checkArg(cond = is.vector(args$featVecShape) | 
                  is.factor(args$featVecShape), 
                msg = "\"featVecShape\" must be a vector.", 
                errs = errs)
  }
  
  
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
  # borderWidth
  errs <- 
    .checkArg(cond = is.numeric(args$borderWidth) & 
                length(args$borderWidth) == 1 &
                args$borderWidth[1] > 0, 
              msg = "\"borderWidth\" must be a positive numeric value.", 
              errs = errs)
  
  #-------------------
  # borderCol
  errs <- 
    .checkArg(cond = (is.numeric(args$borderCol) | 
                        is.character(args$borderCol)) & 
                length(args$borderCol) == 1, 
              msg = paste0("\"borderCol\" must be a single numeric or ", 
                           "character value."), 
              errs = errs)
  
  #-------------------
  # highlightHubs
  errs <- .checkArg(cond = is.logical(args$highlightHubs) & 
                      length(args$highlightHubs) == 1, 
                    msg = "\"highlightHubs\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # hubTransp
  if (!is.null(args$hubTransp)) {
    errs <- .checkArg(cond = is.numeric(args$hubTransp) & 
                        length(args$hubTransp) == 1, 
                      msg = "\"hubTransp\" must be a single numeric value.", 
                      errs = errs)
    
    errs <- .checkArg(cond = args$hubTransp >= 0 & args$hubTransp <= 100, 
                      msg = "\"hubTransp\" must be in [0, 100].", 
                      errs = errs)
  }
  
  #-------------------
  # hubLabelFont
  if (!is.null(args$hubLabelFont)) {
    errs <- .checkArg(cond = is.numeric(args$hubLabelFont) & 
                        length(args$hubLabelFont) == 1, 
                      msg = "\"hubLabelFont\" must be a single numeric value.", 
                      errs = errs)
  }
  
  #-------------------
  # hubBorderWidth
  if (!is.null(args$hubBorderWidth)) {
    errs <- .checkArg(cond = is.numeric(args$hubBorderWidth) & 
                        length(args$hubBorderWidth) == 1 &
                        args$hubBorderWidth[1] > 0, 
                      msg = paste0("\"hubBorderWidth\" must be a positive ", 
                                   "numeric value."), 
                      errs = errs)
  }
  
  #-------------------
  # hubBorderCol
  errs <- 
    .checkArg(cond = (is.numeric(args$hubBorderCol) | 
                        is.character(args$hubBorderCol)) & 
                length(args$hubBorderCol) == 1, 
              msg = paste0("\"hubBorderCol\" must be a single numeric or ", 
                           "character value."), 
              errs = errs)
  
  #-------------------
  # edgeFilter
  choices <- c("none", "threshold", "highestWeight")
  
  args$edgeFilter <- try(match.arg(args$edgeFilter, choices), 
                         silent = TRUE)
  
  if (inherits(args$edgeFilter, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("edgeFilter", choices), 
                errs = errs)
  }
  
  #-------------------
  # edgeFilterPar
  if (args$edgeFilter %in% c("threshold", "highestWeight")) {
    errs <- .checkArg(cond = is.numeric(args$edgeFilterPar) & 
                        length(args$edgeFilterPar) == 1 &
                        args$edgeFilterPar[1] >= 0, 
                      msg = paste0("\"edgeFilterPar\" must be a single ", 
                                   "numeric value >= 0."), 
                      errs = errs)
  }
  
  #-------------------
  # edgeInvisFilter
  choices <- c("none", "threshold", "highestWeight")
  
  args$edgeInvisFilter <- try(match.arg(args$edgeInvisFilter, choices), 
                              silent = TRUE)
  
  if (inherits(args$edgeInvisFilter, "try-error")) {
    errs <- 
      .checkArg(cond = FALSE, 
                msg = .getMatchArgTxt("edgeInvisFilter", choices), 
                errs = errs)
  }
  
  #-------------------
  # edgeInvisPar
  if (args$edgeInvisFilter %in% c("threshold", "highestWeight")) {
    errs <- .checkArg(cond = is.numeric(args$edgeInvisPar) & 
                        length(args$edgeInvisPar) == 1 &
                        args$edgeInvisPar[1] >= 0, 
                      msg = paste0("\"edgeInvisPar\" must be a single ", 
                                   "numeric value >= 0."), 
                      errs = errs)
  }
  #-------------------
  # edgeWidth
  errs <- .checkArg(cond = is.numeric(args$edgeWidth) & 
                      length(args$edgeWidth) == 1 &
                      args$edgeWidth[1] >= 0, 
                    msg = paste0("\"edgeWidth\" must be a single ", 
                                 "numeric value >= 0."), 
                    errs = errs)
  
  #-------------------
  # negDiffCol
  if (!is.null(args$colorNegAsso)) {
    warning("Name of \"colorNegAsso\" changed to \"negDiffCol\"")
    args$negDiffCol <- args$colorNegAsso
    args$colorNegAsso <- NULL
  }
  
  errs <- .checkArg(cond = is.logical(args$negDiffCol) & 
                      length(args$negDiffCol) == 1, 
                    msg = "\"negDiffCol\" must be a single logical value.", 
                    errs = errs)
  
  #-------------------
  # posCol
  if (!is.null(args$posCol)) {
    errs <- 
      .checkArg(cond = length(args$posCol) %in% 1:2, 
                msg = paste0("\"posCol\" must have one or two elements."), 
                errs = errs)
  }
  
  #-------------------
  # negCol
  if (!is.null(args$negCol)) {
    errs <- 
      .checkArg(cond = length(args$negCol) %in% 1:2, 
                msg = paste0("\"negCol\" must have one or two elements."), 
                errs = errs)
  }
  
  #-------------------
  # cut
  if (!is.null(args$cut)) {
    errs <- .checkArg(cond = is.numeric(args$cut) & 
                        length(args$cut) %in% 1:2, 
                      msg = paste0("\"cut\" must be numeric with one or two ", 
                                   "elements."), 
                      errs = errs)
  }
  
  #-------------------
  # edgeTranspLow
  errs <- .checkArg(cond = is.numeric(args$edgeTranspLow) & 
                      length(args$edgeTranspLow) == 1, 
                    msg = "\"edgeTranspLow\" must be a single numeric value.", 
                    errs = errs)
  
  errs <- .checkArg(cond = args$edgeTranspLow >= 0 & args$edgeTranspLow <= 100, 
                    msg = "\"edgeTranspLow\" must be in [0, 100].", 
                    errs = errs)
  
  #-------------------
  # edgeTranspHigh
  errs <- .checkArg(cond = is.numeric(args$edgeTranspHigh) & 
                      length(args$edgeTranspHigh) == 1, 
                    msg = "\"edgeTranspHigh\" must be a single numeric value.", 
                    errs = errs)
  
  errs <- .checkArg(cond = args$edgeTranspHigh >= 0 & 
                      args$edgeTranspHigh <= 100, 
                    msg = "\"edgeTranspHigh\" must be in [0, 100].", 
                    errs = errs)
  #-------------------
  # cexNodes
  errs <- .checkArg(cond = is.numeric(args$cexNodes) & 
                      length(args$cexHubs) == 1, 
                    msg = "\"cexNodes\" must be a single numeric value.", 
                    errs = errs)
  
  #-------------------
  # cexHubs
  errs <- .checkArg(cond = is.numeric(args$cexHubs) & 
                      length(args$cexHubs) == 1, 
                    msg = "\"cexHubs\" must be a single numeric value.", 
                    errs = errs)
  
  #-------------------
  # cexLabels
  errs <- .checkArg(cond = is.numeric(args$cexLabels) & 
                      length(args$cexLabels) == 1, 
                    msg = "\"cexLabels\" must be a single numeric value.", 
                    errs = errs)
  
  #-------------------
  # cexHubLabels
  if (!is.null(args$cexHubLabels)) {
    errs <- 
      .checkArg(cond = is.numeric(args$cexHubLabels) & 
                  length(args$cexHubLabels) == 1, 
                msg = "\"cexHubLabels\" must be a single numeric value.", 
                errs = errs)
  }
  
  #-------------------
  # cexTitle
  errs <- .checkArg(cond = is.numeric(args$cexTitle) & 
                      length(args$cexTitle) == 1, 
                    msg = "\"cexTitle\" must be a single numeric value.", 
                    errs = errs)
  
  #-------------------
  # showTitle
  if (is.null(args$showTitle)) {
    if (twoNets) {
      args$showTitle <- TRUE
    } else {
      args$showTitle <- FALSE
    }
  } else {
    errs <- .checkArg(cond = is.logical(args$showTitle) & 
                        length(args$showTitle) == 1, 
                      msg = "\"showTitle\" must be a single logical value.", 
                      errs = errs)
  }
  
  #-------------------
  # title1
  if (!is.null(args$title1)) {
    errs <- .checkArg(cond = is.character(args$title1) & 
                        length(args$title1) == 1, 
                      msg = "\"title1\" must be a single character value.", 
                      errs = errs)
  }
  
  #-------------------
  # title2
  if (!is.null(args$title2)) {
    errs <- .checkArg(cond = is.character(args$title2) & 
                        length(args$title2) == 1, 
                      msg = "\"title2\" must be a single character value.", 
                      errs = errs)
  }
  
  #-------------------
  # mar
  errs <- .checkArg(cond = is.numeric(args$mar) & length(args$mar) == 4, 
                    msg = "\"mar\" must be a numeric vector of length 4.", 
                    errs = errs)
  
  #-------------------
  # doPlot
  errs <- .checkArg(cond = is.logical(args$doPlot) & length(args$doPlot) == 1, 
                    msg = "\"doPlot\" must be a single logical value.", 
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

