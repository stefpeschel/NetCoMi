argcheck_plot_diffnet<- function(args){

  # x
  stopifnot(class(args$x) == "diffnet")

  # layout
  layout <- args$layout
  if(is.null(layout)){
    layout <- "spring"
  } else{
    layout.tmp <- try(match.fun(layout), silent = TRUE)

    if(class(layout.tmp) == "try-error"){
      if(!is.matrix(layout)){
        stopifnot(layout %in% c("spring", "circle", "groups"))
      }
    } else{
      layout <- layout.tmp
    }
  }
  args$layout <- layout

  # repulsion
  stopifnot(is.numeric(args$repulsion))

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
    stopifnot(is.vector(args$labelPattern) & (length(args$labelPattern) %in% c(3, 5)))
  }


  # labelScale
  stopifnot(is.logical(args$labelScale))


  # labelFont
  stopifnot(is.numeric(args$labelFont))
  if(!is.null(args$hubLabelFont)){
    stopifnot(is.numeric(args$hubLabelFont))
  }


  # rmSingles
  stopifnot(is.logical(args$rmSingles))

  # nodeTransp
  stopifnot(is.numeric(args$nodeTransp))
  stopifnot(args$nodeTransp >= 0 & args$nodeTransp <= 100)

  # edgeTransp
  stopifnot(is.numeric(args$edgeTransp))
  stopifnot(args$edgeTransp >= 0 & args$edgeTransp <= 100)

  # Width
  stopifnot(is.numeric(args$borderWidth) & (args$borderWidth > 0))
  stopifnot(is.numeric(args$edgeWidth) & (args$edgeWidth > 0))

  # edgeFilter
  args$edgeFilter <- match.arg(args$edgeFilter,
                               choices = c("none", "highestDiff"))

  # legend
  stopifnot(is.logical(args$legend))

  # cex
  stopifnot(is.numeric(args$cexLegend) & (length(args$cexLegend) == 1))
  stopifnot(is.numeric(args$cexNodes) & (length(args$cexNodes) == 1))
  stopifnot(is.numeric(args$cexLabels) & (length(args$cexLabels) == 1))
  stopifnot(is.numeric(args$cexTitle) & (length(args$cexTitle) == 1))


  # mar
  stopifnot(is.numeric(args$mar) & (length(args$mar) == 4))


  return(args)

}
