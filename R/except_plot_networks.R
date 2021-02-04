except_plot_networks<- function(args){

  twoNets <- args$x$input$twoNets

  # x
  stopifnot(class(args$x) == "microNetProps")

  # layout
  layout <- args$layout
  if(is.null(layout)){
    layout <- "spring"
  } else{
    layout.tmp <- try(match.fun(layout), silent = TRUE)

    if(class(layout.tmp) == "try-error"){
      if(is.matrix(layout)){
        if(!any(rownames(layout) %in% rownames(args$x$input$adjaMat1))){
          warning("Rownames of layout matrix don't match node names.")
        }
      } else {
        stopifnot(layout %in% c("spring", "circle", "groups"))
      }
    } else{
      layout <- layout.tmp
    }
  }
  args$layout <- layout

  # sameLayout
  sameLayout <- args$sameLayout
  stopifnot(is.logical(sameLayout))

  # layoutGroup
  if(sameLayout & twoNets){
    stopifnot(args$layoutGroup %in% c(1,2, "union"))
  }

  # repulsion
  stopifnot(is.numeric(args$repulsion))

  # groupNames
  if(!is.null(args$groupNames)){
    if((!is.vector(args$groupNames)) || length(args$groupNames) != 2){
      stop("Argument 'groupNames' must be a vector with two elements.")
    }
  }

  # groupsChanged
  stopifnot(is.logical(args$groupsChanged))

  # shortenLabels
  args$shortenLabels <- match.arg(args$shortenLabels,
                                  choices = c("intelligent", "simple", "none"))

  # labelLength
  labelLength <- args$labelLength
  stopifnot(is.numeric(labelLength) || length(labelLength) == 1)
  if(!is.integer(labelLength)){
    if(labelLength %% 1 != 0){
      warning("Argument 'labelLength' coerced to integer.")
    }
    labelLength <-  as.integer(labelLength)
  }
  args$labelLength <- labelLength

  # labelPattern
  if(!is.null(args$labelPattern)){
    stopifnot(is.vector(args$labelPattern) & (length(args$labelPattern) == 3))
  }

  # labelScale
  stopifnot(is.logical(args$labelScale))

  # labelFont
  stopifnot(is.numeric(args$labelFont))
  if(!is.null(args$hubLabelFont)){
    stopifnot(is.numeric(args$hubLabelFont))
  }

  # nodeFilter
  args$nodeFilter <- match.arg(args$nodeFilter, choices = c("none", "highestConnect",
                                                  "highestDegree", "highestBetween",
                                                  "highestClose", "highestEigen",
                                                  "clustTaxon", "clustMin", "names"))

  # rmSingles
  rmSingles <- args$rmSingles
  if(rmSingles[1] == TRUE){
    rmSingles <- "all"
  } else if(rmSingles[1] == FALSE){
    rmSingles <- "none"
  } else{
    rmSingles <- match.arg(rmSingles, choices = c("all", "inboth", "none"))
  }
  args$rmSingles <- rmSingles

  # nodeSize
  args$nodeSize <- match.arg(args$nodeSize, c("fix", "hubs", "degree",
                                              "betweenness", "closeness",
                                              "eigenvector", "counts",
                                              "normCounts", "TSS", "fractions", 
                                              "CSS", "COM", "rarefy", "VST", 
                                              "clr", "mclr"))

  # nodeSizeSpread
  stopifnot(is.numeric(args$nodeSizeSpread) & args$nodeSizeSpread >= 0)

  # nodeColor
  if(!is.null(args$nodeColor)){
    nodeColor <- args$nodeColor
    nodeColor.tmp <- try(match.arg(nodeColor,
                                   choices = c("cluster", "feature", "colorVec")),
                         silent = TRUE)
    if(class(nodeColor.tmp) == "try-error"){
      if(!is.character(nodeColor)){
        stop("Possible values for 'nodeColor' are: 'cluster', 'feature', 'colorVec', or a character defining a color.")
      }
    } else{
      nodeColor <- nodeColor.tmp
    }
    args$nodeColor <- nodeColor
  }

  # sameClustCol
  stopifnot(is.logical(args$sameClustCol))

  # sameColThresh
  sameColThresh <- args$sameColThresh
  stopifnot(is.numeric(sameColThresh) || length(sameColThresh) == 1)
  if(!is.integer(sameColThresh)){
    if(sameColThresh %% 1 != 0){
      warning("Argument 'sameColThresh' coerced to integer.")
    }
    sameColThresh <-  as.integer(sameColThresh)
  }
  args$sameColThresh <- sameColThresh

  # transparencies
  stopifnot(is.numeric(args$nodeTransp))
  stopifnot(args$nodeTransp >= 0 & args$nodeTransp <= 100)
  if(!is.null(args$hubTransp)){
    stopifnot(is.numeric(args$hubTransp))
    stopifnot(args$hubTransp >= 0 & args$hubTransp <= 100)
  }
  stopifnot(is.numeric(args$edgeTranspLow))
  stopifnot(args$edgeTranspLow >= 0 & args$edgeTranspLow <= 100)
  stopifnot(is.numeric(args$edgeTranspHigh))
  stopifnot(args$edgeTranspHigh >= 0 & args$edgeTranspHigh <= 100)

  # Width
  stopifnot(is.numeric(args$borderWidth) & (args$borderWidth > 0))
  if(!is.null(args$hubBorderWidth)){
    stopifnot(is.numeric(args$hubBorderWidth) & (args$hubBorderWidth > 0))
  }
  stopifnot(is.numeric(args$edgeWidth) & (args$edgeWidth > 0))

  # borderCol
  stopifnot(length(args$borderCol) == 1)
  stopifnot(length(args$hubBorderCol) == 1)

  # highlightHubs
  stopifnot(is.logical(args$highlightHubs))

  # edgeFilter
  args$edgeFilter <- match.arg(args$edgeFilter, choices = c("none", "threshold", "highestWeight"))
  args$edgeInvisFilter <- match.arg(args$edgeInvisFilter,
                               choices = c("none", "threshold", "highestWeight"))

  # negDiffCol
  if(!is.null(args$colorNegAsso)){
    warning("Name of 'colorNegAsso' has changed to 'negDiffCol'")
    args$negDiffCol <- args$colorNegAsso
    args$colorNegAsso <- NULL
  }
  stopifnot(is.logical(args$negDiffCol))

  # posCol, negCol
  if(!is.null(args$posCol)){
    stopifnot(length(args$posCol) %in% c(1,2))
  }

  if(!is.null(args$negCol)){
    stopifnot(length(args$negCol) %in% c(1,2))
  }

  # cut
  if(!is.null(args$cut)){
    stopifnot(is.numeric(args$cut) & (length(args$cut) <= 2))
  }

  # cex
  stopifnot(is.numeric(args$cexHubs) & (length(args$cexHubs) == 1))
  stopifnot(is.numeric(args$cexNodes) & (length(args$cexNodes) == 1))
  stopifnot(is.numeric(args$cexLabels) & (length(args$cexLabels) == 1))
  stopifnot(is.numeric(args$cexTitle) & (length(args$cexTitle) == 1))

  # showTitle
  if(is.null(args$showTitle)){
    if(twoNets){
      args$showTitle <- TRUE
    } else{
      args$showTitle <- FALSE
    }
  } else{
    stopifnot(is.logical(args$showTitle))
  }

  # mar
  stopifnot(is.numeric(args$mar) & (length(args$mar) == 4))

  return(args)

}
