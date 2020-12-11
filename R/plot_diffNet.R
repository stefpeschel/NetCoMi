#' Plot method for objects of class diffnet
#'
#' Plot method for objects of class \code{diffnet} inheriting from a call to
#' \code{\link{diffnet}}.
#'
#' @param x object of class \code{diffnet} containing the adjacency matrix 
#'   (absolute differences between associations)
#' @param adjusted logical indicating whether the adjacency matrix based on 
#'   adjusted p-values should be used. Defaults to \code{TRUE}. If \code{FALSE}, 
#'   the adjacency matrix is based on non-adjusted p-values. Ignored 
#'   for the discordant method.
#' @param layout indicates the layout used for defining node positions. Can be
#'   a character with one of the layouts provided by
#'   \code{\link[qgraph]{qgraph}}: \code{"spring"}(default), \code{"circle"},
#'   or \code{"groups"}. Alternatively, the layouts provided by igraph (see
#'   \code{\link[igraph:layout_]{layout\_}}) are accepted (must be given as
#'   character, e.g. \code{"layout_with_fr"}). Can also be a matrix with row
#'   number equal to the  number of nodes and two columns corresponding to the x
#'   and y coordinate.
#' @param repulsion positive numeric value indicating the strength of repulsive
#'   forces in the "spring" layout. Nodes are placed closer together for smaller
#'   values and further apart for higher values. See the \code{repulsion}
#'   argument of \code{\link[qgraph]{qgraph}}.
#' @param labels defines the node labels. Can be a character vector with an
#'   entry for each node. If \code{FALSE}, no labels are plotted. Defaults to
#'   the row/column names of the association matrices.
#' @param shortenLabels options to shorten node labels. Ignored if node labels
#'   are defined via \code{labels}. Possible options are:
#'   \describe{
#'   \item{\code{"intelligent"}}{Elements of \code{charToRm} are removed, labels
#'   are shortened to length \code{labelLength} and duplicates are removed by
#'   using \code{labelPattern}.}
#'   \item{\code{"simple"}}{Elements of \code{charToRm} are removed and labels
#'   are shortened to length \code{labelLength}.}
#'   \item{\code{"none"}}{Default. Original colnames of adjacency matrices are
#'   used} }
#' @param labelLength integer defining the length to which variable names shall
#'   be shortened if \code{shortenLabels} is used. Defaults to 6.
#' @param labelPattern vector of three elements, which is only used if argument
#'   \code{shortenLabels} is set to \code{"intelligent"}. If cutting a node label to
#'   length \code{labelLength} leads to duplicates, the label is shortened
#'   according to \code{labelPattern}, where the first entry gives the length of
#'   the first part, the second entry is used a separator, and the third entry
#'   is the length of the second part. Defaults to c(5, "'", 3). If the data
#'   contains, for example, two bacteria "Streptococcus" and "Streptomyces",
#'   they are by default shortened to "Strep'coc" and "Strep'myc".
#' @param charToRm vector with characters to remove from node names. Ignored if
#'   labels are given via \code{labels}.
#' @param labelScale logical. If \code{TRUE}, node labels are scaled according to
#'   node size
#' @param labelFont integer defining the font of node labels. Defaults to 1.
#' @param rmSingles logical. If \code{TRUE}, unconnected nodes are removed.
#' @param nodeColor character or numeric value specifying node colors. Can also
#'   be a vector with a color for each node.
#' @param nodeTransp an integer between 0 and 100 indicating the transparency of
#'   node colors. 0 means no transparency, 100 means full transparency. Defaults
#'   to 60.
#' @param borderWidth numeric specifying the width of node borders. Defaults to
#'   1.
#' @param borderCol character specifying the color of node borders. Defaults to
#'   "gray80"
#' @param edgeFilter character indicating whether and how edges should be
#'   filtered. Possible values are \code{"none"} (all edges are shown) and
#'   \code{"highestDiff"} (the first x edges with highest absolute difference
#'   are shown). x is defined by \code{edgeFilterPar}.
#' @param edgeFilterPar numeric value specifying the "x" in \code{edgeFilter}.
#' @param edgeWidth numeric specifying the edge width. See argument
#'   \code{"edge.width"} of \code{\link[qgraph]{qgraph}}.
#' @param edgeTransp an integer between 0 and 100 indicating the transparency of
#'   edge colors. 0 means no transparency (default), 100 means full transparency.
#' @param edgeCol character vector specifying the edge colors. Must be of length
#'   6 for the discordant method (default: c("hotpink", "aquamarine", "red",
#'   "orange", "green", "blue")) and of lengths 9 for permutation tests and
#'   Fisher's z-test (default: c("chartreuse2", "chartreuse4", "cyan",
#'   "magenta", "orange", "red", "blue", "black", "purple")).
#' @param title optional character string for the main title.
#' @param legend logical. If \code{TRUE}, a legend is plotted.
#' @param legendGroupnames a vector with two elements giving the group names
#'   shown in the legend.
#' @param legendTitle character specifying the legend title.
#' @param cexNodes numeric scaling node sizes. Defaults to 1.
#' @param cexLabels numeric scaling node labels. Defaults to 1.
#' @param cexTitle numeric scaling the title. Defaults to 1.2.
#' @param cexLegend numeric scaling the legend size. Defaults to 1.
#' @param mar a numeric vector of the form c(bottom, left, top, right) defining
#'   the plot margins. Works similar to the \code{mar} argument in
#'   \code{\link[graphics]{par}}. Defaults to c(2,2,4,6).
#' @param ... further arguments being passed to \code{\link[qgraph]{qgraph}}.
#' @param repulsion integer specifying repulse radius in the spring layout; for
#'   a value lower than 1, nodes are placed further apart
#' @seealso \code{\link{diffnet}}, \code{\link{netConstruct}}
#' @importFrom qgraph qgraph
#' @method plot diffnet
#' @export

plot.diffnet <- function(x,
                         adjusted = TRUE,
                         layout = NULL,
                         repulsion = 1,
                         labels = NULL,
                         shortenLabels = "none",
                         labelLength = 6,
                         labelPattern = NULL,
                         charToRm = NULL,
                         labelScale = TRUE,
                         labelFont = 1,
                         rmSingles = TRUE,
                         nodeColor = "gray90",
                         nodeTransp = 60,
                         borderWidth = 1,
                         borderCol = "gray80",
                         edgeFilter = "none",
                         edgeFilterPar = NULL,
                         edgeWidth = 1,
                         edgeTransp = 0,
                         edgeCol = NULL,
                         title = NULL,
                         legend = TRUE,
                         legendGroupnames = NULL,
                         legendTitle = NULL,
                         cexNodes = 1,
                         cexLabels = 1,
                         cexTitle = 1.2,
                         cexLegend = 1,
                         mar = c(2,2,4,6),
                         ...){

  inputArgs <- c(as.list(environment()), list(...))

  outputArgs <- except_plot_diffnet(inputArgs)
  
  for(i in 1:length(outputArgs)){
    assign(names(outputArgs)[i], outputArgs[[i]])
  }

  corrMat1 <- x$assoMat1
  corrMat2 <- x$assoMat2
  
  if(is.null(x$diffAdjustMat)){
    diffMat <- x$diffMat
    
  } else{
    if(adjusted){
      diffMat <- x$diffAdjustMat
    } else{
      diffMat <- x$diffMat
    }
  }
  
  if(all(diffMat == 0)){
    stop("Network is empty.")
  }

  if(edgeFilter != "none"){
    
    if(edgeFilter == "highestDiff"){
      diffabssort <- sort(abs(diffMat[lower.tri(diffMat)]), decreasing = TRUE)
      cutval <- diffabssort[edgeFilterPar]
      diffMat[diffMat < cutval] <- 0
    }
    
  }

  if(rmSingles){
    corrMat1.orig <- corrMat1
    corrMat2.orig <- corrMat2
    diffMat.orig <- diffMat

    zeros <- sapply(1:nrow(diffMat), function(i){ all(diffMat[i, ] == 0) })
    names(zeros) <- colnames(diffMat)


    if(any(zeros)){
      torm <- which(zeros == TRUE)
    } else{
      torm <- NULL
    }

    if(length(torm) != 0) corrMat1 <- corrMat1[-torm, -torm]
    if(length(torm) != 0) corrMat2 <- corrMat2[-torm, -torm]
    if(length(torm) != 0) diffMat <- diffMat[-torm, -torm]

    kept <- colnames(diffMat.orig)[which(colnames(diffMat.orig) %in% 
                                           colnames(diffMat))]

  }

  if(legend){
    if(is.null(legendTitle)){
      legendTitle = "Associations"
    }
    if(is.null(legendGroupnames)){

      legtitle1 <- paste0("group '" , x$groups[1], "'" )
      legtitle2 <- paste0("group '" , x$groups[2], "'" )

    } else{
      stopifnot(is.vector(legendGroupnames))
      stopifnot(length(legendGroupnames) == 2)
      legtitle1 <- legendGroupnames[1]
      legtitle2 <- legendGroupnames[2]
    }
  }

  #=============================================================================
  # define edge colors

  if(x$diffMethod == "discordant"){

    # create color matrix
    if(is.null(edgeCol)){
      edgeCol <- c("hotpink", "aquamarine", "red", "orange", "green", "blue")
    } else{
      stopifnot(length(col) == 6)
    }

    if(edgeTransp > 0){
      colVec <- col_to_transp(colVec, edgeTransp)
    }

    classMat <- x$classMat
    edgeColMat <- classMat
    
    for(i in 1:nrow(edgeColMat)){
      for(j in 1:ncol(edgeColMat)){
        if(classMat[i,j] %in% c(1,5,9)) edgeColMat[i,j] <- "black"
        if(classMat[i,j] == 2) edgeColMat[i,j] <- edgeCol[3]
        if(classMat[i,j] == 3) edgeColMat[i,j] <- edgeCol[5]
        if(classMat[i,j] == 4) edgeColMat[i,j] <- edgeCol[1]
        if(classMat[i,j] == 6) edgeColMat[i,j] <- edgeCol[6]
        if(classMat[i,j] == 7) edgeColMat[i,j] <- edgeCol[2]
        if(classMat[i,j] == 8) edgeColMat[i,j] <- edgeCol[4]
      }
    }


  } else{

    if(is.null(edgeCol)){
      edgeCol <- c("chartreuse2", "chartreuse4", "cyan", "magenta", "orange",
                   "red", "blue", "black", "purple")
      #plot(1:9, 1:9, col = edgeCol, cex=6, pch=16)
    } else{
      stopifnot(length(edgeCol) == 9)
    }


    # transform correlation matrices to vectors
    lowtri <- lower.tri(corrMat1, diag = FALSE)
    corrVec1 <- corrMat1[lowtri]
    corrVec2 <- corrMat2[lowtri]
    vector_names <- get_vec_names(t(corrMat1))
    names(corrVec1) <- vector_names
    names(corrVec2) <- vector_names

    colVec <- rep("black", length(corrVec1))
    if(any(corrVec1==0) || any(corrVec2 == 0)){
      colVec[corrVec1 > 0 & corrVec2 > 0] <- edgeCol[1]
      colVec[corrVec1 > 0 & corrVec2 == 0] <- edgeCol[2]
      colVec[corrVec1 > 0 & corrVec2 < 0] <- edgeCol[3]
      colVec[corrVec1 < 0 & corrVec2 > 0] <- edgeCol[4]
      colVec[corrVec1 < 0 & corrVec2 == 0] <- edgeCol[5]
      colVec[corrVec1 < 0 & corrVec2 < 0] <- edgeCol[6]
      colVec[corrVec1 == 0 & corrVec2 > 0] <- edgeCol[7]
      colVec[corrVec1 == 0 & corrVec2 == 0] <- edgeCol[8]
      colVec[corrVec1 == 0 & corrVec2 < 0] <- edgeCol[9]
    } else{
      edgeCol <- edgeCol[c(1,3,4,6)]
      colVec[corrVec1 > 0 & corrVec2 > 0] <- edgeCol[1]
      colVec[corrVec1 > 0 & corrVec2 < 0] <- edgeCol[2]
      colVec[corrVec1 < 0 & corrVec2 > 0] <- edgeCol[3]
      colVec[corrVec1 < 0 & corrVec2 < 0] <- edgeCol[4]
    }

    if(edgeTransp > 0){
      colVec <- col_to_transp(colVec, edgeTransp)
    }

    edgeColMat <- corrMat1
    edgeColMat[lower.tri(edgeColMat)] <- colVec
    edgeColMat <- t(edgeColMat)
    edgeColMat[lower.tri(edgeColMat)] <- colVec
  }



  if(rmSingles){
    edgeColMat <- edgeColMat[kept, kept]
  }

  #=============================================================================

  nodeSize <- (7*exp(-ncol(diffMat)/80)+1) * cexNodes

  if(is.null(title)){
    main <- "Differential network"
  } else if(title == FALSE){
    main <- ""
  } else{
    stopifnot(is.character(title))
    main <- title
  }

  #=============================================================================
  # rename taxa

  if(is.null(labels)){

    diffMat <- rename_taxa(diffMat, toRename = "both",
                           shortenLabels = shortenLabels,
                           labelLength = labelLength,
                           labelPattern = labelPattern,
                           charToRm = charToRm)
    labelsout <- rownames(diffMat)

  } else if(is.logical(labels)){
    labelsout <- labels
  } else{
    labelsout <- labels[kept]
  }

  #=============================================================================
  # node colors
  if(nodeTransp > 0){
    nodeColor <- col_to_transp(nodeColor, nodeTransp)
  }

  #=============================================================================


  q <- qgraph(diffMat, layout = layout, color = nodeColor,
              label.scale = labelScale, labels = labelsout,
              label.font = labelFont, label.cex = cexLabels, vsize = nodeSize,
              border.color = borderCol, border.width = borderWidth,
              edge.color = edgeColMat, edge.width = edgeWidth,
              repulsion = repulsion, mar = mar, ...)


  if(legend){

    if(x$diffMethod %in% c("discordant")){
      legend("topright", legend = c(paste0(legtitle1, "  ", legtitle2),
                                    "    0             -",
                                    "    0             +",
                                    "    -             0",
                                    "    -             +",
                                    "    +             0",
                                    "    +             -"),
             col = c("white", edgeCol),
             lty = rep(1,7), lwd = 2, cex = cexLegend, title = legendTitle)
    } else {
      if(length(edgeCol) == 9){
        legend("topright", legend = c(paste0(legtitle1, "  ", legtitle2),
                                      "    +             +",
                                      "    +             0",
                                      "    +             -",
                                      "    -             +",
                                      "    -             0",
                                      "    -             -",
                                      "    0             +",
                                      "    0             0",
                                      "    0             -"),
               col = c("white", edgeCol),
               lty = rep(1,9), lwd = 2, cex = cexLegend, title = legendTitle)
      } else{
        legend("topright", legend = c(paste0(legtitle1, "  ", legtitle2),
                                      "    +             +",
                                      "    +             -",
                                      "    -             +",
                                      "    -             -"),
               col = c("white", edgeCol),
               lty = rep(1,9), lwd = 2, cex = cexLegend, title = legendTitle)
      }
    }
  }

  if(main != ""){
    title(main = list(main, cex = cexTitle ))
  }

  lay <- q$layout
  rownames(lay) <- colnames(diffMat)
  invisible(lay)

}

