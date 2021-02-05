#' @title Plot Method for microNetProps Objects
#'
#' @description Plotting objects of class \code{microNetProps}.
#'
#' @aliases plot plot.microNetProps
#' @usage \method{plot}{microNetProps}(x,
#'   layout = "spring",
#'   sameLayout = FALSE,
#'   layoutGroup = "union",
#'   repulsion = 1,
#'   groupNames = NULL,
#'   groupsChanged = FALSE,
#'   labels = NULL,
#'   shortenLabels = "simple",
#'   labelLength = 6L,
#'   labelPattern = c(5,"'",3),
#'   charToRm = NULL,
#'   labelScale = TRUE,
#'   labelFont = 1,
#'   labelFile = NULL,
#'
#'   # Nodes:
#'   nodeFilter = "none",
#'   nodeFilterPar = NULL,
#'   rmSingles = "none",
#'   nodeSize = "fix",
#'   normPar = NULL,
#'   nodeSizeSpread = 4,
#'   nodeColor = "cluster",
#'   colorVec = NULL,
#'   featVecCol = NULL,
#'   sameFeatCol = TRUE,
#'   sameClustCol = TRUE,
#'   sameColThresh = 2L,
#'   nodeShape = NULL,
#'   featVecShape = NULL,
#'   nodeTransp = 60,
#'   borderWidth = 1,
#'   borderCol = "gray80",
#'
#'   # Hubs:
#'   highlightHubs = TRUE,
#'   hubTransp = NULL,
#'   hubLabelFont = NULL,
#'   hubBorderWidth = NULL,
#'   hubBorderCol = "black",
#'
#'   # Edges:
#'   edgeFilter = "none",
#'   edgeFilterPar = NULL,
#'   edgeInvisFilter = "none",
#'   edgeInvisPar = NULL,
#'   edgeWidth = 1,
#'   negDiffCol = TRUE,
#'   posCol = NULL,
#'   negCol = NULL,
#'   cut = NULL,
#'   edgeTranspLow = 0,
#'   edgeTranspHigh = 0,
#'
#'   # Additional arguments:
#'   cexNodes = 1,
#'   cexHubs = 1.2,
#'   cexLabels = 1,
#'   cexHubLabels = NULL,
#'   cexTitle = 1.2,
#'   showTitle = NULL,
#'   title1 = NULL,
#'   title2 = NULL,
#'   mar = c(1, 3, 3, 3),
#'   ...)
#'
#' @param x  object of class \code{microNetProps}
#' @param layout  indicates the layout used for defining node positions. Can be
#'   a character with one of the layouts provided by
#'   \code{\link[qgraph]{qgraph}}: \code{"spring"} (default), \code{"circle"},
#'   or \code{"groups"}. Alternatively, the layouts provided by igraph (see
#'   \code{\link[igraph:layout_]{layout\_}} are accepted (must be given as
#'   character, e.g. \code{"layout_with_fr"}). Can also be a matrix with row
#'   number equal to the number of nodes and two columns corresponding to the x
#'   and y coordinate.
#' @param sameLayout logical. Indicates whether the same layout should be used
#'   for both networks. Ignored if \code{x} contains only one network. See
#'   argument \code{layoutGroup}.
#' @param layoutGroup  numeric or character. Indicates the group, where the 
#'   layout is taken from if argument \code{sameLayout} is \code{TRUE}. 
#'   The layout is computed for group 1 (and adopted for group 2) if set to "1" 
#'   and it is computed for group 2 if set to "2". Can alternatively be set to
#'   "union" (default) to compute a union of both layouts, where the nodes are 
#'   placed as optimal as possible equally for both networks. 
#' @param repulsion positive numeric value indicating the strength of repulsive
#'   forces in the "spring" layout. Nodes are placed closer together for smaller
#'   values and further apart for higher values. See the \code{repulsion}
#'   argument of \code{\link[qgraph]{qgraph}}.
#' @param groupNames character vector with two entries naming the groups to
#'   which the networks belong. Defaults to the group names returned from
#'   \code{\link{netConstruct}}: If the data set is split according to a group
#'   variable, its factor levels (in increasing order) are used. Ignored if
#'   arguments \code{title1} and \code{title2} are set or if a single network is
#'   plotted.
#' @param groupsChanged logical. Indicates the order in which the networks are
#'   plotted. If \code{TRUE}, the order is exchanged. See details. Defaults to
#'   \code{FALSE}.
#' @param labels defines the node labels. Can be a named character vector, which 
#'   is used for both groups (then, the adjacency matrices in \code{x} must 
#'   contain the same variables). 
#'   Can also be list with two named vectors (names must match the row/column 
#'   names of the adjacency matrices). If \code{FALSE}, no labels are plotted. 
#'   Defaults to the row/column names of the adjacency matrices.
#' @param shortenLabels options to shorten node labels. Ignored if node labels
#'   are defined via \code{labels}. Possible options are:
#'   \describe{
#'   \item{\code{"intelligent"}}{Elements of \code{charToRm} are removed,
#'   labels are shortened to length \code{labelLength} and duplicates are
#'   removed by using \code{labelPattern}.}
#'   \item{\code{"simple"}}{Elements of \code{charToRm} are  removed and labels
#'   are shortened to length \code{labelLength}.}
#'   \item{\code{"none"}}{Default. Original colnames of adjacency matrices are
#'   used.} }
#' @param labelLength integer defining the length to which variable names shall
#'   be shortened if \code{shortenLabels} is used. Defaults to 6.
#' @param labelPattern vector of three elements, which is only used if argument
#'   \code{shortenLabels} is set to \code{"intelligent"}. If cutting a node label
#'   to length \code{labelLength} leads to duplicates, the label is shortened
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
#' @param labelFile optional character of the form "<file name>.txt" naming a 
#'   file where the original and renamed node labels are stored. The file is 
#'   stored into the current working directory.
#' @param nodeFilter character indicating whether and how nodes should be
#'   filtered. Possible values are:
#'   \describe{
#'   \item{\code{"none"}}{Default. All nodes are plotted.}
#'   \item{\code{"highestConnect"}}{x nodes with highest connectivity (sum of
#'   edge weights) are plotted.}
#'   \item{\code{"highestDegree"}, \code{"highestBetween"}, \code{"highestClose"}
#'   , \code{"highestEigen"}}{x nodes with highest
#'   degree/betweenness/closeness/eigenvector centrality are plotted.}
#'   \item{\code{"clustTaxon"}}{Only nodes belonging the same cluster as to
#'   variables that are given as character vector via  \code{nodeFilterPar}.}
#'   \item{\code{"clustMin"}}{Plotted are only nodes belonging to clusters with
#'   a minimum number of nodes of x.}
#'   \item{\code{"names"}}{Character vector with variable names to be plotted}}\cr
#'   Necessary parameters (e.g. "x") are given via the argument
#'   \code{nodeFilterPar}.
#' @param nodeFilterPar parameters needed for the filtering method defined by
#'   \code{nodeFilter}.
#' @param rmSingles character value indicating how to handle unconnected nodes.
#'   Possible values are \code{"all"} (all single nodes are deleted),
#'   \code{"inboth"} (only nodes that are unconnected in both networks are
#'   removed) or \code{"none"} (default; no nodes are removed). Cannot be set to
#'   \code{"all"}, if the same layout is used for both networks.
#' @param nodeSize  character indicating how node sizes should be determined.
#'   Possible values are: \describe{
#'   \item{\code{"fix"}}{Default. All nodes have same size (hub size can be
#'   defined separately via \code{cexHubs}).}
#'   \item{\code{"degree"}, \code{"betweenness"}, \code{"closeness"},
#'   \code{"eigenvector"}}{Size scaled according to node's centrality}
#'   \item{\code{"counts"}}{Size scaled according to the sum of counts (of
#'   microbes or samples, depending on what nodes express).}
#'   \item{\code{"normCounts"}}{Size scaled according to the sum of normalized
#'   counts (of microbes or samples), which are exported by
#'   \code{netConstruct}.}
#'   \item{\code{"TSS", "fractions", "CSS", "COM", "rarefy", "VST", "clr", 
#'   "mclr"}}{Size scaled according to the sum of normalized
#'   counts. Available are the same options as for \code{normMethod} in
#'   \code{\link{netConstruct}}. Parameters are set via \code{normPar}.}
#'   }
#' @param normPar list with parameters passed to the function for normalization
#'   if \code{nodeSize} is set to a normalization method. Used analogously to
#'   \code{normPar} of \code{\link{netConstruct}()}.
#' @param nodeSizeSpread positive numeric value indicating the spread of node
#'   sizes. The smaller the value, the more similar are the node sizes. Node 
#'   sizes are calculated by: (x - min(x)) / (max(x) - min(x)) * nodeSizeSpread 
#'   + cexNodes. For \code{nodeSizeSpread = 4} (default) and \code{cexNodes = 1}, 
#'   node sizes range from 1 to 5.
#' @param nodeColor a character specifying the node colors. Possible values are
#'   \code{"cluster"} (colors according to determined clusters),
#'   \code{"feature"} (colors according to node's features defined by
#'   \code{featVecCol}), \code{"colorVec"} (the vector \code{colorVec}). For the
#'   former two cases, the colors can be specified via \code{colorVec}. If
#'   \code{colorVec} is not defined, the \code{rainbow} function from
#'   \code{grDevices} package  is used. Also accepted is a character value
#'   defining a color, which is used for all nodes. If \code{NULL}, "grey40" is
#'   used for all nodes.
#' @param colorVec a vector or list with two vectors used to specify node colors. 
#'   Different usage depending on the "nodeColor" argument:
#'   \describe{
#'   \item{\code{nodeColor = "cluster"}}{\code{colorVec} must be a vector. 
#'   Depending on the \code{sameClustCol} argument, the colors are used only in 
#'   one or both networks. If the vector is not long enough, a warning is 
#'   returned and colors from \code{rainbow()} are used for the remaining 
#'   clusters.}
#'   \item{\code{nodeColor = "feature"}}{Defines a color for each level of 
#'   \code{featVecCol}. Can be a list with two vectors used for the two networks 
#'   (for a single network, only the first element is used) or a vector, which 
#'   is used for both groups if two networks are plotted.}
#'   \item{\code{nodeColor = "colorVec"}}{\code{colorVec} defines a color for 
#'   each node implying that its names must match the node's names (which is 
#'   also ensured if names match the colnames of the original count matrix). 
#'   Can be a list with two vectors used for the two networks (for a single 
#'   network, only the first element is used) or a vector, which is used for 
#'   both groups if two networks are plotted.}}
#' @param featVecCol  a vector with a feature for each node. Used for coloring
#'   nodes if \code{nodeColor} is set to \code{"feature"}. Is coerced to a
#'   factor. If \code{colorVec} is given, its length must be larger than or
#'   equal to the number of feature levels.
#' @param sameFeatCol logical indicating whether the same color should be used 
#'   for same features in both networks (only used if two networks are plotted, 
#'   \code{nodeColor = "feature"}, and no color vector/list is given (via 
#'   \code{featVecCol})).
#' @param sameClustCol if TRUE (default) and two networks are plotted, clusters
#'   having at least \code{sameColThresh} nodes in common have the same color.
#'   Only used if \code{nodeColor} is set to \code{"cluster"}.
#' @param sameColThresh indicates how many nodes a cluster must have in common
#'   in the two groups to have the same color. See argument \code{sameClustCol}.
#'   Defaults to 2.
#' @param nodeShape character vector specifying node shapes. Possible values
#'   are \code{"circle"} (default), \code{"square"}, \code{"triangle"}, and
#'   \code{"diamond"}. If \code{featVecShape} is not \code{NULL}, the length of
#'   \code{nodeShape} must equal the number of factor levels given by
#'   \code{featVecShape}. Then, each shape is assigned to one factor level (in
#'   increasing order). If \code{featVecShape} is \code{NULL}, the first shape
#'   is used for all nodes. See the example.
#' @param featVecShape a vector with a feature for each node. If not \code{NULL},
#'   a different node shape is used for each feature. Is coerced to factor mode.
#'   The maximum number of factor levels is 4 corresponding to the four possible
#'   shapes defined via \code{nodeShape}.
#' @param nodeTransp an integer between 0 and 100 indicating the transparency of
#'   node colors. 0 means no transparency, 100 means full transparency. Defaults
#'   to 60.
#' @param borderWidth numeric specifying the width of node borders. Defaults to
#'   1.
#' @param borderCol character specifying the color of node borders. Defaults to
#'   "gray80"
#' @param highlightHubs logical indicating if hubs should be highlighted. If
#'   \code{TRUE}, the following features can be defined separately for hubs:
#'   transparency (by \code{hubTransp}), label font (by \code{hubLabelFont}),
#'   border width (by \code{hubBorderWidth}), and border color (by
#'   \code{hubBorderCol}).
#' @param hubTransp numeric between 0 and 100 specifying the color transparency
#'   of hub nodes. See argument \code{nodeTransp}. Defaults to
#'   \code{0.5*nodeTransp}. Ignored if \code{highlightHubs} is \code{FALSE}.
#' @param hubLabelFont integer specifying the label font of hub nodes. Defaults
#'   to \code{2*labelFont}. Ignored if \code{highlightHubs} is \code{FALSE}.
#' @param hubBorderWidth numeric specifying the border width of hub nodes.
#'   Defaults to \code{2*borderWidth}. Ignored if \code{highlightHubs} is
#'   \code{FALSE}.
#' @param hubBorderCol character specifying the border color of hub nodes.
#'   Defaults to \code{"black"}. Ignored if \code{highlightHubs} is \code{FALSE}.
#' @param edgeFilter character indicating whether and how edges should be
#'   filtered. Possible values are \code{"none"} (default; all edges are shown),
#'   \code{"threshold"} (only edges with edge weight of at least x are shown),
#'   \code{"highestWeight"} (the first x edges with highest edge weight are
#'   shown). x is defined by \code{edgeFilterPar}, respectively.
#' @param edgeFilterPar numeric specifying the "x" in \code{edgeFilter}.
#' @param edgeInvisFilter similar to \code{edgeFilter} but the edges are removed
#'   only after computing the layout so that edge removal does not influence the
#'   layout. Defaults to \code{"none"}.
#' @param edgeInvisPar numeric specifying the "x" in \code{edgeInvisFilter}.
#' @param edgeWidth numeric specifying the edge width. See argument
#'   \code{"edge.width"} of \code{\link[qgraph]{qgraph}}.
#' @param negDiffCol logical indicating if edges with a negative corresponding
#'   association should be colored different. If \code{TRUE} (default), argument
#'   \code{posCol} is used for edges with positive association and \code{negCol}
#'   for those with negative association. If \code{FALSE} and for dissimilarity
#'   networks, only \code{posCol} is used.
#' @param posCol vector (character or numeric) with one or two elements
#'   specifying the color of edges with positive weight and also for edges with
#'   negative weight if \code{negDiffCol} is set to \code{FALSE}. The first
#'   element is used for edges with weight below \code{cut} and the second for
#'   edges with weight above \code{cut}. If a single value is given, it is used
#'   for both cases. Defaults to \code{c("#009900", "darkgreen")}.
#' @param negCol vector (character or numeric) with one or two elements
#'   specifying the color of edges with negative weight. The first
#'   element is used for edges with absolute weight below \code{cut} and the
#'   second for edges with absolute weight above \code{cut}. If a single value
#'   is given, it is used for both cases. Ignored if \code{negDiffCol} is
#'   \code{FALSE}. Defaults to \code{c("red", "#BF0000")}.
#' @param cut defines the \code{"cut"} parameter of 
#'   \code{\link[qgraph]{qgraph}}. Can
#'   be either a numeric value (is used for both groups if two networks are
#'   plotted) or vector of length two. The default is set analogous to that in
#'   \code{\link[qgraph]{qgraph}}: "0 for graphs with less then 20 nodes. For
#'   larger graphs the cut value is automatically chosen to be equal to the
#'   maximum of the 75th quantile of absolute edge strengths or the edge
#'   strength corresponding to 2n-th edge strength (n being the number of
#'   nodes.)"
#' @param edgeTranspLow numeric value between 0 and 100 specifying the
#'   transparency of edges with weight below \code{cut}. The higher this value,
#'   the higher the transparency.
#' @param edgeTranspHigh analogous to \code{edgeTranspLow}, but used for edges
#'   with weight ABOVE \code{cut}.
#' @param cexNodes numeric scaling node sizes. Defaults to 1.
#' @param cexHubs numeric scaling hub sizes. Only used if \code{nodeSize} is set
#'   to \code{"hubs"}.
#' @param cexLabels numeric scaling node labels. Defaults to 1.
#' @param cexHubLabels numeric scaling the node labels of hub nodes. Equals 
#'   \code{cexLabels} by default. Ignored, if \code{highlightHubs = FALSE}.
#' @param cexTitle numeric scaling title(s). Defaults to 1.2.
#' @param showTitle if \code{TRUE}, a title is shown for each network, which is
#'   either defined via \code{groupNames}, or \code{title1} and \code{title2}.
#'   Defaults to \code{TRUE} if two networks are plotted and \code{FALSE} for
#'   a single network.
#' @param title1 character giving a title for the first network.
#' @param title2 character giving a title for the second network (if existing).
#' @param mar a numeric vector of the form c(bottom, left, top, right) defining
#'   the plot margins. Works similar to the \code{mar} argument in
#'   \code{\link[graphics]{par}}. Defaults to c(1,3,3,3).
#' @param ... further arguments being passed to \code{\link[qgraph]{qgraph}}, 
#'   which is used for network plotting.
#'
#' @return Returns (invisibly) a list with the following elements: \tabular{ll}{
#'   \code{q1,q2}\tab the qgraph object(s)\cr
#'   \code{layout}\tab layout(s) specifying node positions\cr
#'   \code{nodecolor}\tab one or two vectors with node colors (one for each group) \cr
#'   \code{labels}\tab one or two vectors with node labels}
#'
#' @examples
#' # load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#'
#' # network construction (dissimilarity based network where nodes are subjects):
#' amgut_net <- netConstruct(amgut1.filt,
#'                           measure = "aitchison", filtSamp = "highestFreq",
#'                           filtSampPar = list(highestFreq = 50),
#'                           zeroMethod = "pseudo", sparsMethod = "knn")
#'
#' # network analysis:
#' amgut_props <- netAnalyze(amgut_net, clustMethod = "hierarchical",
#'                           clustPar = list(k=2))
#'
#' ### network plots
#' # clusters are used for node coloring: 
#' plot(amgut_props, nodeColor = "cluster")
#' 
#' # a higher repulsion places nodes with high edge weight closer together:
#' plot(amgut_props, nodeColor = "cluster", repulsion = 1.3)
#'
#' # a feature vector is used for node coloring:
#' set.seed(123456)
#' colVec <- sample(1:5, nrow(amgut1.filt), replace = TRUE)
#' names(colVec) <- rownames(amgut1.filt)
#'
#' plot(amgut_props, nodeColor = "feature", featVecCol = colVec,
#'      colorVec = heat.colors(5))
#'
#' # with a further feature vector for node shapes
#' shapeVec <- sample(1:3, nrow(amgut1.filt), replace = TRUE)
#' names(shapeVec) <- rownames(amgut1.filt)
#'
#' plot(amgut_props, nodeColor = "feature", featVecCol = colVec,
#'      colorVec = heat.colors(5), nodeShape = c("circle", "square", "diamond"),
#'      featVecShape = shapeVec, highlightHubs = FALSE)
#'
#' @importFrom  qgraph qgraph
#' @importFrom graphics par title
#' @importFrom grDevices rainbow
#' @importFrom stats quantile
#' @importFrom WGCNA labels2colors
#' @importFrom utils write.table
#' @method plot microNetProps
#' @export
plot.microNetProps <- function(x,
                               layout = "spring",
                               sameLayout = FALSE,
                               layoutGroup = "union",
                               repulsion = 1,
                               groupNames = NULL,
                               groupsChanged = FALSE,
                               labels = NULL,
                               shortenLabels = "simple",
                               labelLength = 6L,
                               labelPattern = c(5,"'",3),
                               charToRm = NULL,
                               labelScale = TRUE,
                               labelFont = 1,
                               labelFile = NULL,
                               #nodes
                               nodeFilter = "none",
                               nodeFilterPar = NULL,
                               rmSingles = "none",
                               nodeSize = "fix",
                               normPar = NULL,
                               nodeSizeSpread = 4,
                               nodeColor = "cluster",
                               colorVec = NULL,
                               featVecCol = NULL,
                               sameFeatCol = TRUE,
                               sameClustCol = TRUE,
                               sameColThresh = 2L,
                               nodeShape = NULL,
                               featVecShape = NULL,
                               nodeTransp = 60,
                               borderWidth = 1,
                               borderCol = "gray80",
                               #hubs
                               highlightHubs = TRUE,
                               hubTransp = NULL,
                               hubLabelFont = NULL,
                               hubBorderWidth = NULL,
                               hubBorderCol = "black",
                               #edges
                               edgeFilter = "none",
                               edgeFilterPar = NULL,
                               edgeInvisFilter = "none",
                               edgeInvisPar = NULL,
                               edgeWidth = 1,
                               negDiffCol = TRUE,
                               posCol = NULL,
                               negCol = NULL,
                               cut = NULL,
                               edgeTranspLow = 0,
                               edgeTranspHigh = 0,
                               cexNodes = 1,
                               cexHubs = 1.2,
                               cexLabels = 1,
                               cexHubLabels = NULL,
                               cexTitle = 1.2,
                               showTitle = NULL,
                               title1 = NULL,
                               title2 = NULL,
                               mar = c(1, 3, 3, 3),
                               ...){

  inputArgs <- c(as.list(environment()), list(...))

  outputArgs <- except_plot_networks(inputArgs)
  
  for(i in 1:length(outputArgs)){
    assign(names(outputArgs)[i], outputArgs[[i]])
  }

  xgroups <- x$input$groups
  isempty1 <- x$isempty$isempty1
  isempty2 <- x$isempty$isempty2

  # if twoNets==TRUE, two networks are plotted
  twoNets <- x$input$twoNets

  if(twoNets & !is.null(xgroups)){
    deletespace <- function(name){
      if(grepl(" ", strsplit(name, "")[[1]])[1]){
        gname <- strsplit(name, "")[[1]]
        name <- paste(gname[-1], collapse = "")
      }
      name
    }
    xgroups[1] <- deletespace(xgroups[1])
    xgroups[2] <- deletespace(xgroups[2])
  }

  opar <- par()
  
  adja1_orig <- x$input$adjaMat1
  
  #=============================================================================
  # Edge and node filtering

  adja1 <- filter_edges(adja = adja1_orig, edgeFilter = edgeFilter,
                            edgeFilterPar = edgeFilterPar)
  
  keep1 <- filter_nodes(adja = adja1, nodeFilter = nodeFilter,
                           nodeFilterPar = nodeFilterPar, layout = layout,
                           degree = x$centralities$degree1,
                           between = x$centralities$between1,
                           close = x$centralities$close1,
                           eigen = x$centralities$eigenv1,
                           cluster = x$clustering$clust1)

  if(twoNets){
    adja2_orig <- x$input$adjaMat2

    adja2 <- filter_edges(adja = adja2_orig, edgeFilter = edgeFilter,
                         edgeFilterPar = edgeFilterPar)

    keep2 <- filter_nodes(adja = adja2_orig, nodeFilter = nodeFilter,
                             nodeFilterPar = nodeFilterPar, layout = layout,
                             degree = x$centralities$degree2,
                             between = x$centralities$between2,
                             close = x$centralities$close2,
                             eigen = x$centralities$eigenv2,
                             cluster = x$clustering$clust2)

    if(length(keep1) == 0 & length(keep2) == 0){
      stop("No nodes remaining in both networks after node filtering.")
    }

    if(sameLayout){
      keep.tmp <- union(keep1, keep2)
      keep <- colnames(adja1)[colnames(adja1) %in% keep.tmp]
      adja1 <- adja1[keep, keep]
      adja2 <- adja2[keep, keep]
    } else{
      adja1 <- adja1[keep1, keep1]
      adja2 <- adja2[keep2, keep2]
    }

  } else{
    if(length(keep1) == 0){
      stop("No nodes remaining after node filtering.")
    }
    adja1 <- adja1[keep1, keep1]
    adja2 <- NULL
    adja2_orig <- NULL
  }


  if(rmSingles != "none"){
    adja1.tmp <- adja1
    diag(adja1.tmp) <- 0
    zeros1 <- sapply(1:nrow(adja1.tmp), function(i){ all(adja1.tmp[i, ] == 0) })
    names(zeros1) <- colnames(adja1.tmp)

    if(any(zeros1 == TRUE)){
      torm1 <- which(zeros1 == TRUE)
    } else{
      torm1 <- NULL
    }

    if(twoNets){
      adja2.tmp <- adja2
      diag(adja2.tmp) <- 0
      zeros2 <- sapply(1:nrow(adja2.tmp), function(i){ all(adja2.tmp[i, ] == 0) })
      names(zeros2) <- colnames(adja2.tmp)

      if(any(zeros1 == TRUE)){
        torm2 <- which(zeros2 == TRUE)
      } else{
        torm2 <- NULL
      }

      torm_both <- intersect(torm1, torm2)

      if(rmSingles == "inboth" & length(torm_both) != 0){
        adja1 <- adja1.tmp[-torm_both, -torm_both]
        adja2 <- adja2.tmp[-torm_both, -torm_both]

      } else if(rmSingles == "all"){
        if(sameLayout){
          if(length(torm_both) != 0){
            warning('Argument "rmSingles" has been set to "inboth" due to same layout.')
            adja1 <- adja1.tmp[-torm_both, -torm_both]
            adja2 <- adja2.tmp[-torm_both, -torm_both]
          }
        } else{
          if(length(torm1) != 0) adja1 <- adja1.tmp[-torm1, -torm1]
          if(length(torm2) != 0) adja2 <- adja2.tmp[-torm2, -torm2]

        }
      }

      isempty1 <- all(adja1[upper.tri(adja1)] == 0)
      isempty2 <- all(adja2[upper.tri(adja2)] == 0)

    } else{
      if(length(torm1) != 0) adja1 <- adja1.tmp[-torm1, -torm1]
      isempty1 <- all(adja1[upper.tri(adja1)] == 0)
      isempty2 <- NULL
      if(isempty1) stop("Network is empty and cannot be plotted.")
    }
  }

  kept1 <- which(colnames(adja1_orig) %in% colnames(adja1))
  kept2 <- which(colnames(adja2_orig) %in% colnames(adja2))

  #==========================================================================
  # Node colors

  if(nodeColor == "cluster"){
    
    if(!is.null(colorVec)){
      
      if(is.list(colorVec)){
        stop("'colorVec' must be a vector if clusters are used for node colors.
             Set 'sameColThresh' to a high value for different colors in the two networks.")
      }
      
      stopifnot(is.vector(colorVec))
    }

    clust1 <- x$clustering$clust1
    clust2 <- x$clustering$clust2

    if(!(is.null(clust1) & is.null(clust2))){
      clustcolors <- get_clust_cols(clust1 = clust1, clust2 = clust2,
                                    adja1 = adja1, adja2 = adja2,
                                    kept1 = kept1, kept2 = kept2,
                                    isempty1 = isempty1, isempty2 = isempty2,
                                    colorVec = colorVec,
                                    sameClustCol = sameClustCol,
                                    sameColThresh = sameColThresh,
                                    twoNets = twoNets)
      nodecol1 <- clustcolors$clustcol1
      nodecol2 <- clustcolors$clustcol2

    } else{
      message('No clusterings returned from "netAnalyze".')
      nodecol1 <- rep("grey40", ncol(adja1))
      if(twoNets){
        nodecol2 <- rep("grey40", ncol(adja2))
      } else{
        nodecol2 <- NULL
      }
    }

  } else if(nodeColor == "feature"){
    featVecCol <- as.factor(featVecCol)
    stopifnot(all(colnames(adja1) %in% names(featVecCol)))

    if(is.null(colorVec)){
      if(sameFeatCol){
        colorVec1 <- colorVec2 <- grDevices::rainbow(length(levels(featVecCol)))
      } else{
        cols <- grDevices::rainbow(2 * length(levels(featVecCol)))
        colorVec1 <- cols[seq.int(1L, length(cols), 2L)]
        colorVec2 <- cols[seq.int(2L, length(cols), 2L)]
      }
      
    } else{
      if(is.list(colorVec)){
        stopifnot(length(colorVec) == 2)
        colorVec1 <- colorVec[[1]]
        colorVec2 <- colorVec[[2]]
        
        if(length(colorVec1) < length(levels(featVecCol)) || 
           length(colorVec2) < length(levels(featVecCol))){
          stop("Length of color vector(s) (argument 'colorVec') shorter than
               number of levels of 'featVecCol'.")
        }
        
      } else{
        if(length(colorVec) < length(levels(featVecCol))){
          stop(paste("Argument 'featVecCol' has", length(levels(featVecCol)),
                     "levels but 'colorVec' has only length ", length(colorVec)))
        }
        
        colorVec1 <- colorVec2 <- colorVec
      }
    }

    feature1 <- featVecCol[colnames(adja1)]
    nodecol1 <- character(length(feature1))
    
    for(i in seq_along(levels(feature1))){
      nodecol1[feature1 == levels(feature1)[i]] <- colorVec1[i]
    }
    
    nodecol2 <- NULL

    if(twoNets){
      stopifnot(all(colnames(adja2) %in% names(featVecCol)))

      feature2 <- featVecCol[colnames(adja2)]
      nodecol2 <- character(length(feature2))
      
      for(i in seq_along(levels(feature2))){
        nodecol2[feature2 == levels(feature2)[i]] <- colorVec2[i]
      }
    }

  } else if(nodeColor == "colorVec"){
    
    if(is.list(colorVec)){
      stopifnot(length(colorVec) == 2)
      colorVec1 <- colorVec[[1]]
      colorVec2 <- colorVec[[2]]
      
    } else{
      colorVec1 <- colorVec2 <- colorVec
    }
    
    stopifnot(all(colnames(adja1) %in% names(colorVec1)))
    nodecol1 <- colorVec1[colnames(adja1)]
    
    nodecol2 <- NULL

    if(twoNets){
      stopifnot(all(colnames(adja2) %in% names(colorVec2)))
      nodecol2 <- colorVec2[colnames(adja2)]
    }

  } else if(is.character(nodeColor)){
    nodecol1 <- rep(nodeColor, ncol(adja1))
    nodecol2 <- NULL
    if(twoNets){
      nodecol2 <- rep(nodeColor, ncol(adja2))
    }

  } else{
    nodecol1 <- rep("grey40", ncol(adja1))
    if(twoNets){
      nodecol2 <- rep("grey40", ncol(adja2))
    }

  }

  if(nodeTransp > 0){
    nodecol1 <- colToTransp(nodecol1, nodeTransp)
    if(twoNets){
      nodecol2 <- colToTransp(nodecol2, nodeTransp)
    }
  }

  #==========================================================================
  # Node shapes

  if(is.null(featVecShape)){
    if(is.null(nodeShape)){
      nodeShape1 <- nodeShape2 <- "circle"
    } else{
      nodeShape2 <- nodeShape2 <- nodeShape[1]
    }

  } else{
    stopifnot(all(colnames(adja1) %in% names(featVecShape)))

    shape_fact <- as.factor(featVecShape)
    nlevel <- length(levels(shape_fact))

    if(nlevel > 4){
      stop("Argument 'featVecShape' may at most have four factor levels.")
    }

    if(is.null(nodeShape)){
      nodeShape <- c("circle", "square", "triangle", "diamond")

    } else{
      if(length(nodeShape) < nlevel){
        stop(paste("Argument 'nodeShape' must have length", nlevel,
                   "because 'featVecShape' has", nlevel, "levels."))
      }
    }

    shapeVec <- as.character(featVecShape)
    names(shapeVec) <- names(featVecShape)

    for(i in 1:nlevel){
      shapeVec[featVecShape == levels(shape_fact)[i]] <- nodeShape[i]
    }

    nodeShape1 <- shapeVec[colnames(adja1)]

    if(twoNets){
      stopifnot(all(colnames(adja2) %in% names(featVecShape)))
      nodeShape2 <- shapeVec[colnames(adja2)]
    }
  }

  #==========================================================================
  # Node sizes

  if(!is.numeric(nodeSize)){
    if(is.null(x$input$countMat1)){
      counts.tmp <- x$input$countsJoint
    } else{
      counts.tmp <- x$input$countMat1
    }
    nodeSize1 <- get_node_size(nodeSize = nodeSize, normPar = normPar,
                               nodeSizeSpread = nodeSizeSpread,
                               adja = adja1, countMat = counts.tmp,
                               normCounts = x$input$normCounts1,
                               assoType = x$input$assoType, kept = kept1,
                               cexNodes = cexNodes, cexHubs = cexHubs,
                               hubs = x$hubs$hubs1, 
                               highlightHubs = highlightHubs,
                               degree = x$centralities$degree1,
                               between = x$centralities$between1,
                               close = x$centralities$close1,
                               eigen = x$centralities$eigenv1)
    if(twoNets){
      if(is.null(x$input$countMat2)){
        counts.tmp <- x$input$countsJoint
      } else{
        counts.tmp <- x$input$countMat2
      }
      nodeSize2 <- get_node_size(nodeSize = nodeSize, normPar = normPar,
                                 nodeSizeSpread = nodeSizeSpread,
                                 adja = adja2, countMat = counts.tmp,
                                 normCounts = x$input$normCounts2, 
                                 assoType = x$input$assoType, kept = kept2,
                                 cexNodes = cexNodes, cexHubs = cexHubs,
                                 hubs = x$hubs$hubs2, 
                                 highlightHubs = highlightHubs,
                                 degree = x$centralities$degree2,
                                 between = x$centralities$between2,
                                 close = x$centralities$close2,
                                 eigen = x$centralities$eigenv2)
    }
    rm(counts.tmp)
  }

  #===============================================
  # Border colors and fonts

  border1 <- rep(borderCol, nrow(adja1))
  labelFont1 <- rep(labelFont, nrow(adja1))
  borderWidth1 <- rep(borderWidth, nrow(adja1))

  if(highlightHubs){
    if(is.null(hubLabelFont)){
      hubLabelFont <- labelFont * 2
    }
    if(is.null(hubBorderWidth)){
      hubBorderWidth <- borderWidth * 2
    }
    if(is.null(hubTransp)){
      hubTransp <- nodeTransp / 2
    }

    hubs1 <- x$hubs$hubs1
    hubs1 <- hubs1[hubs1 %in% colnames(adja1)]

    border1[match(hubs1, rownames(adja1))] <- hubBorderCol
    labelFont1[match(hubs1, rownames(adja1))] <- hubLabelFont
    borderWidth1[match(hubs1, rownames(adja1))] <- hubBorderWidth

    if(nodeTransp != hubTransp){
      hubcol1 <- colToTransp(nodecol1, hubTransp)
      nodecol1[rownames(adja1) %in% hubs1] <- hubcol1[rownames(adja1) %in% hubs1]
    }

  }

  if(twoNets){
    border2 <- rep(borderCol, nrow(adja2))
    labelFont2 <- rep(labelFont, nrow(adja2))
    borderWidth2 <- rep(borderWidth, nrow(adja2))

    if(highlightHubs){
      hubs2 <- x$hubs$hubs2
      hubs2 <- hubs2[hubs2 %in% colnames(adja2)]

      border2[match(hubs2, rownames(adja2))] <- hubBorderCol
      labelFont2[match(hubs2, rownames(adja2))] <- hubLabelFont
      borderWidth2[match(hubs2, rownames(adja2))] <- hubBorderWidth

      if(nodeTransp != hubTransp){
        hubcol2 <- colToTransp(nodecol2, hubTransp)
        nodecol2[rownames(adja2) %in% hubs2] <- hubcol2[rownames(adja2) %in% hubs2]
      }
    }
  }

  #==========================================================================
  # Title
  
  if(showTitle){
    if(twoNets){
      if(is.null(title1) || is.null(title2)){
        if(is.null(groupNames)){
          main1 = paste0("group '" , xgroups[1], "'" )
          main2 = paste0("group '" , xgroups[2], "'" )
        } else{
          stopifnot(length(groupNames) == 2)

          main1 = as.character(groupNames[1])
          main2 = as.character(groupNames[2])
        }
      } else{
        main1 <- title1
        main2 <- title2
      }
    } else{
      if(is.null(title1)){
        main1 <- ""
        warning('Argument "showTitle" is TRUE, but no title was defined.')
      } else{
        main1 <- title1
      }
    }

  }

  #==========================================================================
  # Define cut parameter for qgraph()

  if(!is.null(cut) & (length(cut) == 1)){
    stopifnot(is.numeric(cut))
    cut1 <- cut2 <- cut
  } else if(!is.null(cut) & (length(cut) == 2)){
    stopifnot(is.numeric(cut))
    cut1 <- cut[1]
    cut2 <- cut[2]
  } else{
    weights1 <- adja1[upper.tri(adja1)]
    weights1 <- weights1[abs(weights1) > 0]
    nNodes1 <- ncol(adja1)

    if(twoNets){
      weights2 <- adja2[upper.tri(adja2)]
      weights2 <- weights2[abs(weights2) > 0]
      nNodes2 <- ncol(adja2)
    } else{
      weights2 <- 0
      nNodes2 <- 1
    }

    if (length(weights1) > (3 * nNodes1) || length(weights2) > (3 * nNodes2)) {
      cut1 <- max(sort(abs(weights1), decreasing = TRUE)[2 * nNodes1],
                  quantile(abs(weights1), 0.75))
      cut2 <- ifelse(twoNets, 
                     max(sort(abs(weights2), decreasing = TRUE)[2 * nNodes2],
                         quantile(abs(weights2), 0.75)), 
                     0)


    } else if (length(weights1) > 1 || length(weights2) > 1){
      cut1 <- quantile(abs(weights1), 0.75)
      cut2 <- ifelse(twoNets, quantile(abs(weights2), 0.75), 0)

    } else {
      cut1 <- 0
      cut2 <- 0
    }
  }

  #==========================================================================
  # Layout

  if(twoNets & sameLayout & layoutGroup == "union"){
    graph1 <- igraph::graph_from_adjacency_matrix(adja1, weighted = TRUE)
    graph2 <- igraph::graph_from_adjacency_matrix(adja2, weighted = TRUE)
    
    graph_u <- igraph::union(graph1, graph2)
    
    # element-wise minimum of edge weights
    E(graph_u)$weight <- pmin(E(graph_u)$weight_1, 
                              E(graph_u)$weight_2, 
                              na.rm = TRUE) 
    
    graph_u <- igraph::delete_edge_attr(graph_u, "weight_1")
    graph_u <- igraph::delete_edge_attr(graph_u, "weight_2")
    
    if(is.function(layout)){
      lay1 <- layout(graph_u)

    } else{
      adja_u <- as.matrix(as_adjacency_matrix(graph_u, attr = "weight",
                                              sparse = T))
      
      lay1 <- qgraph(adja_u, color = nodecol1, layout = layout, 
                     vsize = nodeSize1, labels = colnames(adja1),
                     label.scale = labelScale, 
                     border.color = border1, border.width = borderWidth1,
                     label.font = labelFont1, label.cex = cexLabels,
                     edge.width = edgeWidth, repulsion = repulsion, cut = cut1,
                     shape = nodeShape1, DoNotPlot = TRUE, ...)$layout
    }
    
    rownames(lay1) <- rownames(adja1)
    lay2 <- lay1
    
  } else{
    # Group 1
    if(is.matrix(layout)){
      lay1 <- layout
      
    } else if(is.function(layout)){
      gr <- graph_from_adjacency_matrix(adja1, mode = "undirected",
                                        weighted = TRUE, diag = FALSE)
      lay1 <- layout(gr)
      rownames(lay1) <- colnames(adja1)
    } else{
      lay1 <- qgraph(adja1, color = nodecol1, layout = layout, vsize = nodeSize1,
                     labels = colnames(adja1),label.scale = labelScale,
                     border.color = border1, border.width = borderWidth1,
                     label.font = labelFont1, label.cex = cexLabels,
                     edge.width = edgeWidth, repulsion = repulsion, cut = cut1,
                     DoNotPlot = TRUE, shape = nodeShape1, ...)$layout
      rownames(lay1) <- rownames(adja1)
    }
    lay2 <- NULL
    
    #--------------------------------------------
    # Group 2

    if(twoNets){
      if(is.matrix(layout)){
        lay2 <- layout
        
      } else if(is.function(layout)){
        gr <- igraph::graph_from_adjacency_matrix(adja2, mode = "undirected",
                                                  weighted = TRUE, diag = FALSE)
        
        lay2 <- layout(gr)
        rownames(lay2) <- colnames(adja2)
        
      } else{
        lay2 <- qgraph(adja2, color = nodecol2, layout = layout, vsize = nodeSize2,
                       labels = colnames(adja2), label.scale = labelScale,
                       border.color = border2, border.width = borderWidth2,
                       label.font = labelFont2, label.cex = cexLabels,
                       edge.width = edgeWidth,  repulsion = repulsion, cut = cut2,
                       DoNotPlot = TRUE, shape = nodeShape2, ...)$layout
        
        rownames(lay2) <- rownames(adja2)
        
      }
      
      if(sameLayout){
        if(layoutGroup == 1){
          lay2 <- lay1
        } else{
          lay1 <- lay2
        }
      }
    }
  }

  #==========================================================================
  # Edge colors

  if(is.null(posCol)){
    poscol1 <- "#009900"
    poscol2 <- "darkgreen"
  } else{
    if(length(posCol) == 2){
      poscol1 <- posCol[1]
      poscol2 <- posCol[2]
    } else{
      poscol1 <- poscol2 <- posCol
    }

  }

  if(is.null(negCol)){
    negcol1 <- "red"
    negcol2 <- "#BF0000"
  } else if(negCol == FALSE){
    negcol1 <- poscol1
    negcol2 <- poscol2
  } else{
    if(length(negCol) == 2){
      negcol1 <- negCol[1]
      negcol2 <- negCol[2]
    } else{
      negcol1 <- negcol2 <- negCol
    }
  }

  if(edgeTranspLow > 0){  # transparency for values below cut
    poscol1 <- colToTransp(poscol1, edgeTranspLow)
    negcol1 <- colToTransp(negcol1, edgeTranspLow)
  }

  if(edgeTranspHigh > 0){  # transparency for values above cut
    poscol2 <- colToTransp(poscol2, edgeTranspHigh)
    negcol2 <- colToTransp(negcol2, edgeTranspHigh)
  }

  if(!is.null(x$input$assoEst1)){
    assoMat1 <- x$input$assoEst1
  } else{
    assoMat1 <- x$input$dissScale1
  }

  #--------------------------------------------
  # Edge colors (group 1)

  if(x$input$assoType == "dissimilarity"){
    assoMat1 <- x$input$dissEst1
  } else{
    assoMat1 <- x$input$assoEst1
  }

  if(negDiffCol){
    colmat1 <- matrix(NA, ncol(assoMat1), ncol(assoMat1), 
                      dimnames = dimnames(assoMat1))
    colmat1[assoMat1 < 0] <- negcol1
    colmat1[assoMat1 >= 0] <- poscol1

    colmat1 <- colmat1[kept1, kept1]

    colmat1[abs(adja1) >= cut1 & colmat1 == poscol1] <- poscol2
    colmat1[abs(adja1) >= cut1 & colmat1 == negcol1] <- negcol2
    
  } else{
    colmat1 <- matrix(poscol1, ncol(assoMat1), ncol(assoMat1), 
                      dimnames = dimnames(assoMat1))
    colmat1 <- colmat1[kept1, kept1]
    colmat1[abs(adja1) >= cut1 & colmat1 == poscol1] <- poscol2
  }

  #--------------------------------------------
  
  # # plot dissimilarity values instead of similarities
  # if(plotdiss){
  #   adja1 <- 1-adja1
  #   adja1[adja1 == 1] <- 0
  #   
  #   adja2 <- 1-adja2
  #   adja2[adja2 == 1] <- 0
  # }

  #--------------------------------------------
  # Labels (group 1)

  labels1.orig <- rownames(adja1)
  
  if(is.null(labels)){

    adja.tmp <- rename_taxa(adja1, toRename = "both", 
                            shortenLabels = shortenLabels,
                            labelLength = labelLength, 
                            labelPattern = labelPattern,
                            charToRm = charToRm)
    
    labels1 <- rownames(adja.tmp)
    
  } else if(is.logical(labels)){
    labels1 <- labels
    
  } else if(is.list(labels)){
    stopifnot(all(colnames(adja1) %in% names(labels[[1]])))
    labels1 <- labels[[1]][colnames(adja1)]
    
  } else if(is.vector(labels)){
    stopifnot(all(colnames(adja1) %in% names(labels)))
    labels1 <- labels[colnames(adja1)]
    
  } else{
    stop("Argument 'labels' must be either a named vector with a label for each 
         node, or a list of length 2 naming the nodes in each network.")
  }
  
  #--------------------------------------------
  # Label size (group 1)

  if(highlightHubs){
    if(is.null(cexHubLabels)){
      cexHubLabels <- cexLabels
    }
    
    cexLabels1 <- rep(cexLabels, length(labels1))
    cexLabels1[match(hubs1, rownames(adja1))] <- cexHubLabels
    
  } else{
    cexLabels1 <- cexLabels
  }
  
  #--------------------------------------------
  # Filter edges without influencing the layout (group 1)

  if(edgeInvisFilter != "none"){
    adja1 <- filter_edges(adja1, edgeFilter = edgeInvisFilter,
                         edgeFilterPar = edgeInvisPar)
  }

  #==========================================================================
  if(twoNets){
    # Edge colors (group 2)

    if(x$input$assoType == "dissimilarity"){
      assoMat2 <- x$input$dissEst2
    } else{
      assoMat2 <- x$input$assoEst2
    }


    if(negDiffCol){
      colmat2 <- matrix(NA, ncol(assoMat2), ncol(assoMat2),
                        dimnames = dimnames(assoMat2))
      
      colmat2[assoMat2 < 0] <- negcol1
      colmat2[assoMat2 >= 0] <- poscol1

      colmat2 <- colmat2[kept2, kept2]

      colmat2[abs(adja2) >= cut2 & colmat2 == poscol1] <- poscol2
      colmat2[abs(adja2) >= cut2 & colmat2 == negcol1] <- negcol2
      
    } else{
      colmat2 <- matrix(poscol1, ncol(assoMat2), ncol(assoMat2),
                        dimnames = dimnames(assoMat2))
      
      colmat2 <- colmat2[kept2, kept2]
      colmat2[abs(adja2) >= cut2 & colmat2 == poscol1] <- poscol2
    }

    #--------------------------------------------
    # Labels (group 2)
    labels2.orig <- rownames(adja2)
    
    if(is.null(labels)){
      adja.tmp <- rename_taxa(adja2, toRename = "both", 
                              shortenLabels = shortenLabels,
                              labelLength = labelLength, 
                              labelPattern = labelPattern,
                              charToRm = charToRm)
      labels2 <- rownames(adja.tmp)

    } else if(is.logical(labels)){
      labels2 <- labels
      
    } else if(is.list(labels)){
      stopifnot(all(colnames(adja2) %in% names(labels[[2]])))
      labels2 <- labels[[2]][colnames(adja2)]
      
    } else if(is.vector(labels)){
      stopifnot(all(colnames(adja2) %in% names(labels)))
      labels2 <- labels[colnames(adja2)]
    }

    # store label names to file
    if(!is.null(labelFile) && !is.logical(labels)){
      stopifnot(is.character(labelFile))
      
      labelMat <- cbind(labels1, labels1.orig)
      labelMat <- rbind(c("Renamed", "Original"), labelMat)
      labelMat <- rbind(c(" ", " "), labelMat)
      labelMat <- rbind(c("Network 1:", " "), labelMat)
      
      dframe <- data.frame(labelMat, stringsAsFactors=FALSE)
      
      # apply format over each column for alignment
      if(shortenLabels != "none"){
        formWidth <- max(suppressWarnings(sum(as.numeric(labelPattern), 
                                              na.rm = TRUE)),
                         labelLength) + 4
      } else{
        formWidth <- 12
      }
      
      dframe <- apply(dframe, 2, format, width = formWidth)
      
      utils::write.table(dframe, labelFile, quote=FALSE, row.names=FALSE, 
                         col.names = FALSE, append = FALSE)
      
      labelMat <- cbind(labels2, labels2.orig)
      labelMat <- rbind(c("Renamed", "Original"), labelMat)
      labelMat <- rbind(c(" ", " "), labelMat)
      labelMat <- rbind(c("Network 2:", " "), labelMat)
      labelMat <- rbind(c(" ", " "), labelMat)
      
      dframe <- data.frame(labelMat, stringsAsFactors=FALSE)
      
      # apply format over each column for alignment
      dframe <- apply(dframe, 2, format, width = formWidth)
      
      utils::write.table(dframe, labelFile, quote=FALSE, row.names=FALSE,
                         col.names = FALSE, append = TRUE)
    }
    
    #--------------------------------------------
    # Label size (group 1)
    
    if(highlightHubs){
      cexLabels2 <- rep(cexLabels, length(labels2))
      cexLabels2[match(hubs2, rownames(adja2))] <- cexHubLabels
      
    } else{
      cexLabels2 <- cexLabels
    }

    #--------------------------------------------
    # Filter edges without influencing the layout (group 2)

    if(edgeInvisFilter != "none"){
      adja2 <- filter_edges(adja2, edgeFilter = edgeInvisFilter,
                           edgeFilterPar = edgeInvisPar)
    }

  } else{ #single network

    # store label names to file
    if(!is.null(labelFile) && !is.logical(labels)){
      stopifnot(is.character(labelFile))
      
      labelMat <- cbind(labels1, labels1.orig)
      labelMat <- rbind(c("Renamed", "Original"), labelMat)

      dframe <- data.frame(labelMat, stringsAsFactors=FALSE)
      
      # apply format over each column for alignment
      if(shortenLabels != "none"){
        formWidth <- max(suppressWarnings(sum(as.numeric(labelPattern), 
                                              na.rm = TRUE)),
                         labelLength) + 4
      } else{
        formWidth <- 12
      }
      dframe <- apply(dframe, 2, format, width = formWidth)
      
      utils::write.table(dframe, labelFile, quote=FALSE, row.names=FALSE, 
                         col.names = FALSE, append = FALSE)
    }
  }

  #=============================================================================
  # Plot network(s)

  if(twoNets){
      par(mfrow = c(1,2))

    if(groupsChanged){
      q2 <- qgraph(adja2, color = nodecol2, layout = lay2, vsize = nodeSize2,
                   labels = labels2, label.scale = labelScale,
                   border.color = border2, border.width = borderWidth2,
                   label.font = labelFont2, label.cex = cexLabels2,
                   edge.width = edgeWidth, edge.color = colmat2,
                   repulsion = repulsion, cut = cut2,
                   mar = mar, shape = nodeShape2, ...)
      if(showTitle) title(main = list(main2, cex = cexTitle))
    }

    q1 <- qgraph(adja1, color = nodecol1, layout = lay1, vsize = nodeSize1,
                 labels = labels1,label.scale = labelScale,
                 border.color = border1, border.width = borderWidth1,
                 label.font = labelFont1, label.cex = cexLabels1,
                 edge.width = edgeWidth, edge.color = colmat1,
                 repulsion = repulsion, cut = cut1,
                 mar = mar,  shape = nodeShape1, ...)
    if(showTitle) title(main = list(main1, cex = cexTitle))

    if(!groupsChanged){
      q2 <-  qgraph(adja2, color = nodecol2, layout = lay2, vsize = nodeSize2,
                    labels = labels2, label.scale = labelScale,
                    border.color = border2, border.width = borderWidth2,
                    label.font = labelFont2, label.cex = cexLabels2,
                    edge.width = edgeWidth, edge.color = colmat2,
                    repulsion = repulsion, cut = cut2,
                    mar = mar,  shape = nodeShape2, ...)
      if(showTitle) title(main = list(main2, cex = cexTitle))
    }

    #---------------------------------------------
  } else{

    q1 <- qgraph(adja1, color = nodecol1, layout = lay1, vsize = nodeSize1,
                 labels = labels1, label.scale = labelScale,
                 border.color = border1, border.width = borderWidth1,
                 label.font = labelFont1, label.cex = cexLabels1,
                 edge.width = edgeWidth, edge.color = colmat1,
                 repulsion = repulsion, cut = cut1,
                 mar = mar,  shape = nodeShape1, ...)
    if(showTitle) title(main = list(main1, cex = cexTitle))
    q2 <- NULL
    labels2 <- NULL
  }

  layout.out <- list(layout1 = lay1, layout2 = lay2)
  nodecol.out <- list(nodecol1 = nodecol1, nodecol2 = nodecol2)
  labels.out <- list(labels1 = labels1, labels2 = labels2)


  output <- list(q1 = q1, q2 = q2, layout = layout.out,
                 nodecolor = nodecol.out, labels = labels.out)

  #------------------------------------------------------------------
  #------------------------------------------------------------------

  par(mfrow = opar$mfrow, mar = opar$mar)

  invisible(output)


}

