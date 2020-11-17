#' @title Microbiome Network Analysis
#'
#' @description Determine network properties for objects of class
#'   \code{microNet}.
#'
#' @param net object of class \code{microNet} inheriting from a call to
#'   \code{\link{netConstruct}}
#' @param clustMethod character indicating the clustering algorithm. Possible
#'   values are "hierarchical" for a hierarchical algorithm based on
#'   dissimilarity values, or the clustering methods provided by the igraph
#'   package (see \code{\link[igraph]{communities}} for possible methods).
#'   Defaults to \code{"cluster_fast_greedy"}.
#' @param clustPar list with parameters passed to the clustering functions.
#'   If hierarchical clustering is used, the parameters are passed to
#'   \code{\link[stats]{hclust}} as well as \code{\link[stats]{cutree}}.
#' @param clustPar2 optional list with clustering parameters for the second
#'   network. If \code{NULL} and \code{net} contains two networks,
#'   \code{clustPar} is used for the second network as well.
#' @param hubPar character vector with one or more elements (centrality 
#'   measures) used for identifying hub nodes. Possible values are \code{degree},
#'   \code{betweenness}, \code{closeness}, and \code{eigenvector}. If multiple
#'   measures are given, hubs are nodes with highest centrality for all selected
#'   measures. See details.
#' @param hubQuant quantile used for determining hub nodes. Defaults to 0.95.
#' @param lnormFit hubs are nodes with a centrality value above the 95\%
#'   quantile of the fitted log-normal distribution (if \code{lnormFit = TRUE})
#'   or of the empirical distribution of centrality values (if
#'   \code{lnormFit = FALSE}, which is default).
#' @param connect logical indicating whether edge and vertex connectivity should
#'   be calculated. Might be disabled to reduce execution time.
#' @param weightDeg if \code{TRUE}, the weighted degree is used (see
#'   \code{\link[igraph]{strength}}). Is automatically set to TRUE for a fully
#'   connected network.
#' @param normDeg,normBetw,normClose,normEigen if \code{TRUE}, a normalized
#'   version of the respective centrality values is returned. By default, all
#'   centralities are normalized.
#' @param verbose logical. If \code{TRUE} (default), status messages are shown.
#' @details \strong{Identification of hub nodes:}\cr\cr
#'   Hubs are the nodes with highest centrality values for one or more
#'   centrality measure.\cr\cr The "highest values" regarding a centrality
#'   measure are defined as values lying above a certain quantile (defined by
#'   \code{hubQuant}) either of the empirical distribution of the centralities
#'   or of the fitted log-normal distribution (using
#'   \code{\link[MASS]{fitdistr}}).\cr\cr If \code{clustPar} contains multiple
#'   measures, the centrality values of a hub node must be above the given
#'   quantile for all these measures.
#'
#' @return An object of class \code{microNetProps} containing the following
#'   elements for both groups, respectively: \tabular{ll}{
#'   \code{clustering}\tab determined clusters\cr
#'   \code{centralities}\tab determined degree, betweenness, closeness, and
#'   eigenvector centrality values\cr
#'   \code{hubs}\tab names of the hub nodes (taxa or subjects)\cr
#'   \code{globalProps}\tab values for global network properties: average path
#'   length, global clustering coefficient, modularity, vertex connectivity,
#'   edge connectivity, density (ratio of the number of edges and the number of
#'   possible edges)\cr
#'   \code{paramsProperties}\tab Given parameters used for computing network
#'   properties \cr
#'   \code{input}\tab input adopted from \code{net} \cr
#'   \code{paramsNetConstruct }\tab parameters used for network construction by
#'   \code{netConstruct}\cr
#'   \code{isempty}\tab indicates if the networks are empty (do not contain any
#'   edges)\cr}
#'
#' @examples
#' # load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#'
#' # network construction:
#' amgut_net1 <- netConstruct(amgut1.filt, measure = "pearson",
#'                            filtTax = "highestVar",
#'                            filtTaxPar = list(highestVar = 50),
#'                            zeroMethod = "pseudo", normMethod = "clr")
#'
#' # network analysis:
#' amgut_props1 <- netAnalyze(amgut_net1, clustMethod = "cluster_fast_greedy",
#'                            hubPar = "eigenvector")
#' amgut_props2 <- netAnalyze(amgut_net1, clustMethod = "cluster_fast_greedy",
#'                            hubPar = c("degree", "betweenness", "closeness"))
#'
#' summary(amgut_props1, showCentr = "eigenvector", numbNodes = 15L, digits = 3L)
#' summary(amgut_props2)
#' 
#' # network plot:
#' plot(amgut_props1)
#' plot(amgut_props2)
#'
#'
#' # dissimilarity-based network (where nodes are subjects):
#' amgut_net3 <- netConstruct(amgut1.filt, measure = "aitchison",
#'                            filtSamp = "highestFreq",
#'                            filtSampPar = list(highestFreq = 30),
#'                            zeroMethod = "multRepl", sparsMethod = "knn")
#'
#' amgut_props3 <- netAnalyze(amgut_net3, clustMethod = "hierarchical",
#'                            clustPar = list(k = 3))
#'
#' plot(amgut_props3)
#'
#' @seealso \code{\link{netConstruct}} for network construction,
#'   \code{\link{netCompare}} for network comparison,
#'   \code{\link{diffnet}} for constructing differential networks.
#' @import igraph
#' @importFrom MASS fitdistr
#' @export

netAnalyze <- function(net,
                       clustMethod = NULL,
                       clustPar = NULL,
                       clustPar2 = NULL,
                       hubPar = "eigenvector",
                       hubQuant = 0.95,
                       lnormFit = FALSE,
                       connect = TRUE,
                       weightDeg = FALSE,
                       normDeg = TRUE,
                       normBetw = TRUE,
                       normClose = TRUE,
                       normEigen = TRUE,
                       verbose = TRUE){
  x <- net
  stopifnot(class(x) == "microNet")

  clustMethod <- match.arg(clustMethod, c("none", "hierarchical",
                                          "cluster_edge_betweenness",
                                          "cluster_fast_greedy",
                                          "cluster_leading_eigen",
                                          "cluster_louvain",
                                          "cluster_optimal",
                                          "cluster_spinglass",
                                          "cluster_walktrap"))

  if(is.null(clustPar2)) clustPar2 <- clustPar

  hubPar <- match.arg(hubPar, c("degree", "betweenness", "closeness",
                                "eigenvector"), several.ok = TRUE)

  if(!is.null(clustPar)) stopifnot(is.list(clustPar))
  stopifnot(is.logical(lnormFit))
  stopifnot(is.logical(weightDeg))
  stopifnot(is.logical(normDeg))
  stopifnot(is.logical(normBetw))
  stopifnot(is.logical(normClose))
  stopifnot(is.logical(normEigen))

  twoNets <- x$twoNets
  groups <- x$groups
  adja1 <- x$adjaMat1

  #=============================================================================

  isempty1 <- all(adja1[lower.tri(adja1)] == 0)
  isempty2 <- NULL
  if(!twoNets & isempty1) stop("Network is empty.")

  if(all(adja1[lower.tri(adja1)] > 0) & !weightDeg){
    if(verbose){
      message(paste0('Weighted degree used (unweighted degree not meaningful ',
                     'for a fully connected network).'))
    }

    weightDeg <- TRUE
    if(normDeg){
      normDeg <- FALSE
      if(verbose) message("Argument 'normDeg' ignored for weighted degree.")
    }
  }


  if(twoNets){
    adja2 <- x$adjaMat2

    isempty2 <- all(adja2[lower.tri(adja2)] == 0)
    if(isempty1 & isempty2) stop("There are no connected nodes in both networks.")

    if(all(adja2[lower.tri(adja2)] > 0) & !weightDeg){
      if(verbose){
        message(paste0('Weighted degree used (unweighted degree not meaningful ',
                       'for a fully connected network).'))
      }
      weightDeg <- TRUE
    }
  }


  props1 <- calc_props(adjaMat = adja1, dissMat = x$dissMat1,
                      weighted = x$parameters$weighted, isempty = isempty1,
                      clustMethod = clustMethod, clustPar = clustPar,
                      hubPar = hubPar, hubQuant = hubQuant, connect = connect,
                      lnormFit = lnormFit, weightDeg = weightDeg,
                      normDeg = normDeg, normBetw = normBetw,
                      normClose = normClose, normEigen = normEigen)
  props2 <- NULL

  if(twoNets){
    props2 <- calc_props(adjaMat = adja2, dissMat = x$dissMat2,
                        weighted = x$parameters$weighted, isempty = isempty2,
                        clustMethod = clustMethod, clustPar = clustPar2,
                        hubPar = hubPar, hubQuant = hubQuant, connect = connect,
                        lnormFit = lnormFit, weightDeg = weightDeg,
                        normDeg = normDeg, normBetw = normBetw,
                        normClose = normClose, normEigen = normEigen)
  }

  output <- list(clustering = list(clust1 = props1$clust,
                                   clust2 = props2$clust,
                                   tree1 = props1$tree,
                                   tree2 = props2$tree),
                 centralities = list(degree1 = props1$deg,
                                     degree2 = props2$deg,
                                     between1 = props1$betw,
                                     between2 = props2$betw,
                                     close1 = props1$close,
                                     close2 = props2$close,
                                     eigenv1 = props1$eigen,
                                     eigenv2 = props2$eigen),
                 hubs = list(hubs1 = props1$hubs, hubs2 = props2$hubs),
                 globalProps = list(avPath1 = props1$avPath,
                                    avPath2 = props2$avPath,
                                    clustCoef1 = props1$clustCoef,
                                    clustCoef2 = props2$clustCoef,
                                    modularity1 = props1$modul,
                                    modularity2 = props2$modul,
                                    vertConnect1 = props1$vertconnect,
                                    vertConnect2 = props2$vertconnect,
                                    edgeConnect1 = props1$edgeconnect,
                                    edgeConnect2 = props2$edgeconnect,
                                    density1 = props1$density,
                                    density2 = props2$density),
                 paramsProperties = list(clustMethod = clustMethod,
                                         clustPar = clustPar,
                                         clustPar2 = clustPar2,
                                         hubPar = hubPar,
                                         hubQuant = hubQuant,
                                         lnormFit = lnormFit,
                                         connect = connect,
                                         weightDeg = weightDeg,
                                         normDeg = normDeg,
                                         normBetw = normBetw,
                                         normClose = normClose,
                                         normEigen = normEigen),
                 input = list(assoMat1 = x$assoMat1,
                              assoMat2 = x$assoMat2,
                              dissMat1 = x$dissMat1,
                              dissMat2 = x$dissMat2,
                              simMat1 = x$simMat1,
                              simMat2 = x$simMat2,
                              adjaMat1 = x$adjaMat1,
                              adjaMat2 = x$adjaMat2,
                              assoEst1 = x$assoEst1,
                              assoEst2 = x$assoEst2,
                              dissEst1 = x$dissEst1,
                              dissEst2 = x$dissEst2,
                              dissScale1 = x$dissScale1,
                              dissScale2 = x$dissScale2,
                              countMat1 = x$countMat1,
                              countMat2 = x$countMat2,
                              normCounts1 = x$normCounts1,
                              normCounts2 = x$normCounts2,
                              twoNets = twoNets,
                              groups = groups,
                              matchDesign = x$matchDesign,
                              assoType = x$assoType,
                              softThreshPower = x$softThreshPower,
                              sampleSize = x$sampleSize),
                 paramsNetConstruct = x$parameters,
                 isempty = list(isempty1 = isempty1, isempty2 = isempty2))

  class(output) <- "microNetProps"
  return(output)
}



