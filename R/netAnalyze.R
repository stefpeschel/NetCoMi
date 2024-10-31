#' @title Microbiome Network Analysis
#'
#' @description Determine network properties for objects of class
#'   \code{microNet}.
#'   
#' @usage netAnalyze(net,
#'            # Centrality related:
#'            centrLCC = TRUE,
#'            weightDeg = FALSE,
#'            normDeg = TRUE,
#'            normBetw = TRUE,
#'            normClose = TRUE,
#'            normEigen = TRUE,
#'            
#'            # Cluster related:
#'            clustMethod = NULL,
#'            clustPar = NULL,
#'            clustPar2 = NULL,
#'            weightClustCoef = TRUE,
#'            
#'            # Hub related:
#'            hubPar = "eigenvector",
#'            hubQuant = 0.95,
#'            lnormFit = FALSE,
#'            
#'            # Graphlet related:
#'            graphlet = TRUE,
#'            orbits = c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1),
#'            gcmHeat = TRUE,
#'            gcmHeatLCC = TRUE,
#'            
#'            # Further arguments:
#'            avDissIgnoreInf = FALSE,
#'            sPathAlgo = "dijkstra",
#'            sPathNorm = TRUE,
#'            normNatConnect = TRUE,
#'            connectivity = TRUE,
#'            verbose = 1
#'            )
#'   
#' @details 
#' \strong{Definitions:}\cr
#'   \describe{
#'   \item{(Connected) Component}{Subnetwork where any two nodes are connected 
#'   by a path.}
#'   \item{Number of components}{Number of connected components. Since a single 
#'   node is connected to itself by the trivial path, each single node is a 
#'   component.}
#'   \item{Largest connected component (LCC)}{The connected component with 
#'   highest number of nodes.}
#'   \item{Shortest paths}{Computed using \code{\link[igraph]{distances}}. The 
#'   algorithm is defined via \code{sPathAlgo}. Normalized shortest paths (if 
#'   \code{sPathNorm} is \code{TRUE}) are calculated by dividing the shortest 
#'   paths by the average dissimilarity (see below).}
#'   }
#'   \strong{Global network properties:}\cr
#'   \describe{
#'   \item{Relative LCC size}{= (# nodes in the LCC) / (# nodes in the complete 
#'   network)}
#'   \item{Clustering Coefficient}{The weighted (global) clustering coefficient 
#'   is the arithmetic mean of the local clustering coefficient defined by 
#'   Barrat et al. (computed by \code{\link[igraph]{transitivity}} with 
#'   type = "barrat"), where NAs are ignored. \cr 
#'   The unweighted (global) clustering coefficient is computed using 
#'   \code{\link[igraph]{transitivity}} with type = "global".}
#'   \item{Modularity}{The modularity score for the determined clustering is 
#'   computed using \code{\link[igraph]{modularity.igraph}}.}
#'   \item{Positive edge percentage}{Percentage of edges with positive estimated
#'   association of the total number of edges.}
#'   \item{Edge density}{Computed using \code{\link[igraph]{edge_density}}.}
#'   \item{Natural connectivity}{Computed using 
#'   \code{\link[pulsar]{natural.connectivity}}. The "norm" parameter is 
#'   defined by \code{normNatConnect}.}
#'   \item{Vertex / Edge connectivity}{Computed using 
#'   \code{\link[igraph]{vertex_connectivity}} and 
#'   \code{\link[igraph]{edge_connectivity}}. Both equal zero for a 
#'   disconnected network.}
#'   \item{Average dissimilarity}{Computed as the mean of dissimilarity values 
#'   (lower triangle of \code{dissMat}). By \code{avDissIgnoreInf} is specified 
#'   whether to ignore infinite dissimilarities. The average dissimilarity of an
#'   empty network is 1.}
#'   \item{Average path length}{Computed as the mean of shortest paths 
#'   (normalized or unnormalized). The av. path length of an empty network is 1.}
#'   }
#'   \strong{Clustering algorithms:}\cr
#'   \describe{
#'   \item{Hierarchical clustering}{Based on dissimilarity values. Computed 
#'   using \code{\link[stats]{hclust}} and \code{\link[stats]{cutree}}.}
#'   \item{cluster_optimal}{Modularity optimization. See 
#'   \code{\link[igraph]{cluster_optimal}}.}
#'   \item{cluster_fast_greedy}{Fast greedy modularity optimization. See 
#'   \code{\link[igraph]{cluster_fast_greedy}}.}
#'   \item{cluster_louvain}{Multilevel optimization of modularity. See 
#'   \code{\link[igraph]{cluster_louvain}}.}
#'   \item{cluster_edge_betweenness}{Based on edge betweenness. Dissimilarity
#'   values are used. See \code{\link[igraph]{cluster_edge_betweenness}}.}
#'   \item{cluster_leading_eigen}{Based on leading eigenvector of the community
#'   matrix. See \code{\link[igraph]{cluster_leading_eigen}}.}
#'   \item{cluster_spinglass}{Find communities via spin-glass model and 
#'   simulated annealing. See \code{\link[igraph]{cluster_spinglass}}.}
#'   \item{cluster_walktrap}{Find communities via short random walks. See
#'   \code{\link[igraph]{cluster_walktrap}}.}
#'   }
#'   \strong{Hubs:}\cr
#'   Hubs are nodes with highest centrality values for one or more
#'   centrality measures. The "highest values" regarding a centrality
#'   measure are defined as values lying above a certain quantile (defined by
#'   \code{hubQuant}) either of the empirical distribution of the centralities 
#'   (if \code{lnormFit = FALSE}) or of the fitted log-normal distribution 
#'   (if \code{lnormFit = TRUE}; \code{\link[MASS]{fitdistr}} is used for 
#'   fitting). The quantile is set using \code{hubQuant}.\cr
#'   If \code{clustPar} contains multiple measures, the centrality values of a 
#'   hub node must be above the given quantile for all measures at the same time.
#'   \cr\cr
#'   \strong{Centrality measures:}\cr
#'   Via \code{centrLCC} is decided whether centralities should be calculated 
#'   for the whole network or only for the largest connected component. In the 
#'   latter case (\code{centrLCC = FALSE}), nodes outside the LCC have a 
#'   centrality value of zero.
#'   \describe{
#'   \item{Degree}{The unweighted degree (normalized and unnormalized) is 
#'   computed using \code{\link[igraph]{degree}}, and the weighted degree 
#'   using \code{\link[igraph]{strength}}.}
#'   \item{Betweenness centrality}{The unnormalized and normalized betweenness 
#'   centrality is computed using \code{\link[igraph]{betweenness}}.}
#'   \item{Closeness centrality}{Unnormalized: closeness = sum(1/shortest paths)
#'   \cr Normalized: closeness_unnorm = closeness / (# nodes – 1)}
#'   \item{Eigenvector centrality}{If \code{centrLCC == FALSE} and the network 
#'   consists of more than one components: The eigenvector centrality (EVC) is 
#'   computed for each component separately (using 
#'   \code{\link[igraph]{eigen_centrality}}) and scaled according to component 
#'   size to overcome the fact that nodes in smaller components have a higher 
#'   EVC. If \code{normEigen == TRUE}, the EVC values are divided by the maximum 
#'   EVC value. EVC of single nodes is zero.\cr\cr
#'   Otherwise, the EVC is computed for the LCC using 
#'   \code{\link[igraph]{eigen_centrality}} (scale argument is set according to 
#'   \code{normEigen}).}
#'   }
#'   \strong{Graphlet-based properties:}\cr
#'   \describe{
#'   \item{Orbit counts}{Count of node orbits in graphlets with 2 to 4 nodes. 
#'   See Hocevar and Demsar (2016) for details. The \code{\link[orca]{count4}} 
#'   function from \code{orca} package is used for orbit counting.
#'   }
#'   \item{Graphlet Correlation Matrix (GCM)}{Matrix with Spearman's 
#'   correlations between the network's (non-redundant) node orbits 
#'   (Yaveroglu et al., 2014).}
#'   }
#'   By default, only the 11 non-redundant orbits are used. These are grouped 
#'   according to their role: orbit 0 represents the degree, orbits (2, 5, 7) 
#'   represent nodes within a chain, orbits (8, 10, 11) represent nodes in a 
#'   cycle, and orbits (6, 9, 4, 1) represent a terminal node.
#'
#' @param net object of class \code{microNet} (returned by 
#'   \code{\link{netConstruct}}).
#' @param centrLCC logical indicating whether to compute centralities only 
#'   for the largest connected component (LCC). If \code{TRUE} 
#'   (default), centrality values of disconnected components are zero. 
#' @param avDissIgnoreInf logical indicating whether to ignore infinities when 
#'   calculating the average dissimilarity. If \code{FALSE} (default), infinity 
#'   values are set to 1.
#' @param sPathAlgo character indicating the algorithm used for computing
#'   the shortest paths between all node pairs. \code{\link[igraph]{distances}} 
#'   (igraph) is used for shortest path calculation. 
#'   Possible values are: "unweighted", "dijkstra" (default), "bellman-ford", 
#'   "johnson", or "automatic" (the fastest suitable algorithm is used). The 
#'   shortest paths are needed for the average (shortest) path length and 
#'   closeness centrality.
#' @param sPathNorm logical. If \code{TRUE} (default), shortest paths are 
#'   normalized by average dissimilarity (only connected nodes are considered),  
#'   i.e., a path is interpreted as steps with average dissimilarity. 
#'   If \code{FALSE}, the shortest path is the minimum sum of dissimilarities 
#'   between two nodes.
#' @param normNatConnect logical. If \code{TRUE} (default), the normalized 
#'   natural connectivity is returned.
#' @param clustMethod character indicating the clustering algorithm. Possible
#'   values are \code{"hierarchical"} for a hierarchical algorithm based on
#'   dissimilarity values, or the clustering methods provided by the igraph
#'   package (see \code{\link[igraph]{communities}} for possible methods).
#'   Defaults to \code{"cluster_fast_greedy"} for association-based networks and
#'   to \code{"hierarchical"} for sample similarity networks.
#' @param clustPar list with parameters passed to the clustering functions.
#'   If hierarchical clustering is used, the parameters are passed to
#'   \code{\link[stats]{hclust}} and \code{\link[stats]{cutree}} (default is 
#'   \code{list(method = "average", k = 3)}.
#' @param clustPar2 same as \code{clustPar} but for the second network. 
#'   If \code{NULL} and \code{net} contains two networks,
#'   \code{clustPar} is used for the second network as well.
#' @param weightClustCoef logical indicating whether (global) clustering 
#'   coefficient should be weighted (\code{TRUE}, default) or unweighted 
#'   (\code{FALSE}).
#' @param hubPar character vector with one or more elements (centrality 
#'   measures) used for identifying hub nodes. Possible values are \code{degree},
#'   \code{betweenness}, \code{closeness}, and \code{eigenvector}. If multiple
#'   measures are given, hubs are nodes with highest centrality for all selected
#'   measures. See details.
#' @param hubQuant quantile used for determining hub nodes. Defaults to 0.95.
#' @param lnormFit hubs are nodes with a centrality value above the 95\%
#'   quantile of the fitted log-normal distribution (if \code{lnormFit = TRUE})
#'   or of the empirical distribution of centrality values 
#'   (\code{lnormFit = FALSE}; default).
#' @param weightDeg logical. If \code{TRUE}, the weighted degree is used (see
#'   \code{\link[igraph]{strength}}). Default is \code{FALSE}. 
#'   Is automatically set to \code{TRUE} for a fully connected (dense) network.
#' @param normDeg,normBetw,normClose,normEigen logical. If \code{TRUE} 
#'   (default for all measures), a normalized version of the respective 
#'   centrality values is returned.
#' @param connectivity logical. If \code{TRUE} (default), edge and vertex 
#'   connectivity are calculated. Might be disabled to reduce execution time.
#' @param graphlet logical. If \code{TRUE} (default), graphlet-based network 
#'   properties are computed: orbit counts as defined by \code{orbits} 
#'   and the corresponding Graphlet Correlation Matrix (\code{gcm}).
#' @param orbits numeric vector with integers from 0 to 14 defining the orbits 
#'   used for calculating the GCM. Minimum length is 2. 
#'   Defaults to c(0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11), 
#'   thus excluding redundant orbits such as the orbit o3.
#' @param gcmHeat logical indicating if a heatmap of the GCM(s) should be 
#'   plotted. Default is \code{TRUE}.
#' @param gcmHeatLCC logical. The GCM heatmap is plotted for the LCC if 
#'   \code{TRUE} (default) and for the whole network if \code{FALSE}. 
#' @param verbose integer indicating the level of verbosity. Possible values:
#'   \code{"0"}: no messages, \code{"1"}: only important messages,
#'   \code{"2"}(default): all progress messages are shown. Can also be logical.
#'
#' @return An object of class \code{microNetProps} containing the following
#'   elements: \tabular{ll}{
#'   \code{lccNames1, lccNames2}\tab Names of nodes in the largest connected 
#'   component(s).\cr
#'   \code{compSize1, compSize2}\tab Matrix/matrices with component sizes (1st 
#'   row: sizes; 2nd row: number of components with the respective size)\cr
#'   \code{clustering}\tab Determined clusters in the whole network (and 
#'   corresponding trees if hierarchical clustering is used)\cr
#'   \code{clusteringLCC}\tab Clusters (and optional trees) of the
#'   largest connected component.\cr
#'   \code{centralities}\tab Centrality values\cr
#'   \code{hubs}\tab Names of hub nodes\cr
#'   \code{globalProps}\tab Global network properties of the whole network.\cr
#'   \code{globalPropsLCC}\tab Global network properties of the largest 
#'   component.\cr
#'   \code{graphlet}\tab Graphlet-based properties (orbit counts and GCM).\cr
#'   \code{graphletLCC}\tab Graphlet-based properties of the largest connected 
#'   component.\cr
#'   \code{paramsProperties}\tab Given parameters used for network analysis\cr
#'   \code{paramsNetConstruct}\tab Parameters used for network construction 
#'   (inherited from \code{\link{netConstruct}}).\cr
#'   \code{input}\tab Input inherited from \code{\link{netConstruct}}.\cr
#'   \code{isempty}\tab Indicates whether network(s) is/are empty.
#'   }
#'
#' @examples
#' # Load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#'
#' # Network construction
#' amgut_net1 <- netConstruct(amgut1.filt, measure = "pearson",
#'                            filtTax = "highestVar",
#'                            filtTaxPar = list(highestVar = 50),
#'                            zeroMethod = "pseudoZO", normMethod = "clr",
#'                            sparsMethod = "threshold", thresh = 0.4)
#'
#' # Network analysis
#' 
#' # Using eigenvector centrality as hub score
#' amgut_props1 <- netAnalyze(amgut_net1, clustMethod = "cluster_fast_greedy",
#'                            hubPar = "eigenvector")
#'                            
#' summary(amgut_props1, showCentr = "eigenvector", numbNodes = 15L, digits = 3L)
#' 
#' # Using degree, betweenness and closeness centrality as hub scores
#' amgut_props2 <- netAnalyze(amgut_net1, clustMethod = "cluster_fast_greedy",
#'                            hubPar = c("degree", "betweenness", "closeness"))
#'
#' summary(amgut_props2, showCentr = "all",  numbNodes = 5L, digits = 5L)
#' 
#' # Calculate centralities only for the largest connected component
#' amgut_props3 <- netAnalyze(amgut_net1, centrLCC = TRUE, 
#'                            clustMethod = "cluster_fast_greedy",
#'                            hubPar = "eigenvector")
#'
#' summary(amgut_props3, showCentr = "none", clusterLCC = TRUE)
#' 
#' # Network plot
#' plot(amgut_props1)
#' plot(amgut_props2)
#' plot(amgut_props3)
#' 
#' #----------------------------------------------------------------------------
#' # Plot the GCM heatmap
#' plotHeat(mat = amgut_props1$graphletLCC$gcm1,
#'          pmat = amgut_props1$graphletLCC$pAdjust1,
#'          type = "mixed",
#'          title = "GCM", 
#'          colorLim = c(-1, 1),
#'          mar = c(2, 0, 2, 0))

#' # Add rectangles
#' graphics::rect(xleft   = c( 0.5,  1.5, 4.5,  7.5),
#'                ybottom = c(11.5,  7.5, 4.5,  0.5),
#'                xright  = c( 1.5,  4.5, 7.5, 11.5),
#'                ytop    = c(10.5, 10.5, 7.5,  4.5),
#'                lwd = 2, xpd = NA)
#' 
#' text(6, -0.2, xpd = NA, 
#'      "Significance codes:  ***: 0.001;  **: 0.01;  *: 0.05")
#' 
#' #----------------------------------------------------------------------------
#' # Dissimilarity-based network (where nodes are subjects)
#' amgut_net4 <- netConstruct(amgut1.filt, measure = "aitchison",
#'                            filtSamp = "highestFreq",
#'                            filtSampPar = list(highestFreq = 30),
#'                            zeroMethod = "multRepl", sparsMethod = "knn")
#'
#' amgut_props4 <- netAnalyze(amgut_net4, clustMethod = "hierarchical",
#'                            clustPar = list(k = 3))
#'
#' plot(amgut_props4)
#'
#' @seealso \code{\link{netConstruct}} for network construction,
#'   \code{\link{netCompare}} for network comparison,
#'   \code{\link{diffnet}} for constructing differential networks,
#'   \code{\link{plot.microNetProps}} for the plot method, and
#'   \code{\link{summary.microNetProps}} for the summary method.
#'   
#' @references
#'   \insertRef{hocevar2016computation}{NetCoMi}\cr\cr
#'   \insertRef{yaveroglu2014revealing}{NetCoMi}
#'   
#' @import igraph
#' @importFrom MASS fitdistr
#' @importFrom graphics layout text rect
#' @export

netAnalyze <- function(net,
                       centrLCC = TRUE,
                       weightDeg = FALSE,
                       normDeg = TRUE,
                       normBetw = TRUE,
                       normClose = TRUE,
                       normEigen = TRUE,
                       clustMethod = NULL,
                       clustPar = NULL,
                       clustPar2 = NULL,
                       weightClustCoef = TRUE,
                       hubPar = "eigenvector",
                       hubQuant = 0.95,
                       lnormFit = FALSE,
                       graphlet = TRUE,
                       orbits = c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1),
                       gcmHeat = TRUE,
                       gcmHeatLCC = TRUE,
                       avDissIgnoreInf = FALSE,
                       sPathAlgo = "dijkstra",
                       sPathNorm = TRUE,
                       normNatConnect = TRUE,
                       connectivity = TRUE,
                       verbose = 1) {
  
  # Check input arguments
  argsIn <- as.list(environment())

  if (verbose == 2) message("Checking input arguments ... ", appendLF = FALSE)
  
  argsOut <- .checkArgsNetAna(argsIn)
  
  for (i in 1:length(argsOut)) {
    assign(names(argsOut)[i], argsOut[[i]])
  }
  
  if (verbose == 2) message("Done.")
  
  x <- net

  twoNets <- x$twoNets
  groups <- x$groups
  adja1 <- x$adjaMat1
  
  #=============================================================================
  
  isempty1 <- all(adja1[lower.tri(adja1)] == 0)
  isempty2 <- NULL
  if (!twoNets & isempty1) stop("Network is empty.")
  
  if (all(adja1[lower.tri(adja1)] > 0) & !weightDeg) {
    if (verbose > 0) {
      message(paste0('Weighted degree used (unweighted degree not meaningful ',
                     'for a fully connected network).'))
    }
    
    weightDeg <- TRUE
  }
  
  
  if (twoNets) {
    adja2 <- x$adjaMat2
    
    isempty2 <- all(adja2[lower.tri(adja2)] == 0)
    if (isempty1 & isempty2) 
      stop("There are no connected nodes in both networks.")
    
    if (all(adja2[lower.tri(adja2)] > 0) & !weightDeg) {
      if (verbose > 0) {
        message(paste0('Weighted degree used (unweighted degree not meaningful ',
                       'for a fully connected network).'))
      }
      weightDeg <- TRUE
    }
  }
  
  if (weightDeg && normDeg) {
    normDeg <- FALSE
    if (verbose > 0) {
      message("Argument 'normDeg' set to FALSE")
      message("(no normalization implemented for weighted degree).")
    } 
  }
  
  if (verbose == 2 && twoNets) {
    message("Compute network properties for group 1 ...")
  }
  
  props1 <- .calcProps(adjaMat = adja1,
                       dissMat = x$dissMat1,
                       assoMat = x$assoMat1,
                       centrLCC = centrLCC,
                       avDissIgnoreInf = avDissIgnoreInf,
                       sPathNorm = sPathNorm,
                       sPathAlgo = sPathAlgo,
                       normNatConnect = normNatConnect,
                       weighted = x$parameters$weighted,
                       isempty = isempty1,
                       clustMethod = clustMethod,
                       clustPar = clustPar,
                       weightClustCoef = weightClustCoef,
                       hubPar = hubPar,
                       hubQuant = hubQuant,
                       lnormFit = lnormFit,
                       connectivity = connectivity,
                       graphlet = graphlet,
                       orbits = orbits,
                       weightDeg = weightDeg,
                       normDeg = normDeg,
                       normBetw = normBetw,
                       normClose = normClose,
                       normEigen = normEigen,
                       verbose = verbose)
  
  props2 <- NULL
  
  if (twoNets) {
    if (verbose == 2) {
      message("__________________________________________")
      message("Compute network properties for group 2 ...")
    }
    
    props2 <- .calcProps(adjaMat = adja2,
                         dissMat = x$dissMat2,
                         assoMat = x$assoMat2,
                         centrLCC = centrLCC,
                         avDissIgnoreInf = avDissIgnoreInf,
                         sPathNorm = sPathNorm,
                         sPathAlgo = sPathAlgo,
                         normNatConnect = normNatConnect,
                         weighted = x$parameters$weighted,
                         isempty = isempty2,
                         clustMethod = clustMethod,
                         clustPar = clustPar2,
                         weightClustCoef = weightClustCoef,
                         hubPar = hubPar,
                         hubQuant = hubQuant,
                         lnormFit = lnormFit,
                         connectivity = connectivity,
                         graphlet = graphlet,
                         orbits = orbits,
                         weightDeg = weightDeg,
                         normDeg = normDeg,
                         normBetw = normBetw,
                         normClose = normClose,
                         normEigen = normEigen,
                         verbose = verbose)
  }
  
  #=============================================================================
  # Plot GCM(s)
  
  if (graphlet && gcmHeat) {
    opar <- par() 
    
    if (twoNets) {
      mat <- matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE)

      graphics::layout(mat = mat)
      
      gobj1 <- if (gcmHeatLCC) props1$graphlet_lcc else props1$graphlet
      gobj2 <- if (gcmHeatLCC) props2$graphlet_lcc else props2$graphlet
      
      plotHeat(mat = gobj1$gcm,
               pmat = gobj1$pAdjust,
               type = "mixed",
               title = "GCM1", 
               colorLim = c(-1, 1),
               mar = c(0, 0, 2, 0))
      .addOrbRect()
      
      plotHeat(mat = gobj2$gcm,
               pmat = gobj2$pAdjust,
               type = "mixed",
               title = "GCM2", 
               colorLim = c(-1, 1),
               mar = c(0, 0, 2, 0))
      .addOrbRect()
      
      gcmtest <- testGCM(gobj1, gobj2, verbose = FALSE)
      
      plotHeat(mat = gcmtest$diff, 
               pmat = gcmtest$pAdjustDiff,
               type = "mixed", 
               title = "Differences between GCs (GCM1-GCM2)",
               mar = c(0, 0, 2, 0))
      
      .addOrbRect()

      graphics::text(14, 1, "Significance codes:\n***: 0.001;  **: 0.01;  *: 0.05", 
           xpd = NA, pos = 4)

      par(mfrow = opar$mfrow)
      
    } else {
      gobj1 <- if (gcmHeatLCC) props1$graphlet_lcc else props1$graphlet
      
      plotHeat(mat = gobj1$gcm,
               pmat = gobj1$pAdjust,
               type = "mixed",
               title = "GCM", 
               colorLim = c(-1, 1),
               mar = c(2, 0, 2, 0))
      
      .addOrbRect()
      
      graphics::text(6, -0.2, 
           "Significance codes:  ***: 0.001;  **: 0.01;  *: 0.05", xpd = NA)
    }
  }
  
  #=============================================================================
  output <- list()
  
  output$lccNames1 <- props1$lccNames
  
  output$compSize1 <- props1$compSize
  
  if (twoNets) {
    output$lccNames2 <- props2$lccNames
    output$compSize2 <- props2$compSize
  }
  
  output$clustering <- list(clust1 = props1$clust,
                            clust2 = props2$clust,
                            tree1 = props1$tree,
                            tree2 = props2$tree)
  
  output$clusteringLCC <- list(clust1 = props1$clust_lcc,
                               clust2 = props2$clust_lcc,
                               tree1 = props1$tree_lcc,
                               tree2 = props2$tree_lcc)
  
  output$centralities <- list(degree1 = props1$deg,
                              degree2 = props2$deg,
                              between1 = props1$betw,
                              between2 = props2$betw,
                              close1 = props1$close,
                              close2 = props2$close,
                              eigenv1 = props1$eigen,
                              eigenv2 = props2$eigen)
  
  output$hubs <- list(hubs1 = props1$hubs, hubs2 = props2$hubs)
  
  output$globalProps <- list(nComp1 = props1$nComp,
                             nComp2 = props2$nComp,
                             avDiss1 = props1$avDiss,
                             avDiss2 = props2$avDiss,
                             avPath1 = props1$avPath,
                             avPath2 = props2$avPath,
                             clustCoef1 = props1$clustCoef,
                             clustCoef2 = props2$clustCoef,
                             modularity1 = props1$modul,
                             modularity2 = props2$modul,
                             vertConnect1 = props1$vertconnect,
                             vertConnect2 = props2$vertconnect,
                             edgeConnect1 = props1$edgeconnect,
                             edgeConnect2 = props2$edgeconnect,
                             natConnect1 = props1$natConnect,
                             natConnect2 = props2$natConnect,
                             density1 = props1$density,
                             density2 = props2$density,
                             pep1 = props1$pep,
                             pep2 = props2$pep)
  
  output$globalPropsLCC <- list(lccSize1 = props1$lccSize, 
                                lccSize2 = props2$lccSize,
                                lccSizeRel1 = props1$lccSizeRel,
                                lccSizeRel2 = props2$lccSizeRel,
                                avDiss1 = props1$avDiss_lcc,
                                avDiss2 = props2$avDiss_lcc,
                                avPath1 = props1$avPath_lcc,
                                avPath2 = props2$avPath_lcc,
                                clustCoef1 = props1$clustCoef_lcc,
                                clustCoef2 = props2$clustCoef_lcc,
                                modularity1 = props1$modul_lcc,
                                modularity2 = props2$modul_lcc,
                                vertConnect1 = props1$vertconnect_lcc,
                                vertConnect2 = props2$vertconnect_lcc,
                                edgeConnect1 = props1$edgeconnect_lcc,
                                edgeConnect2 = props2$edgeconnect_lcc,
                                natConnect1 = props1$natConnect_lcc,
                                natConnect2 = props2$natConnect_lcc,
                                density1 = props1$density_lcc,
                                density2 = props2$density_lcc,
                                pep1 = props1$pep_lcc,
                                pep2 = props2$pep_lcc)
  
  if (graphlet) {
    output$graphlet <- list(ocount1 = props1$graphlet$ocount,
                            ocount2 = props2$graphlet$ocount,
                            gcm1 = props1$graphlet$gcm,
                            gcm2 = props2$graphlet$gcm,
                            pvals1 = props1$graphlet$pvals,
                            pvals2 = props2$graphlet$pvals,
                            pAdjust1 = props1$graphlet$pAdjust,
                            pAdjust2 = props2$graphlet$pAdjust)
    
    output$graphletLCC <- list(ocount1 = props1$graphlet_lcc$ocount,
                            ocount2 = props2$graphlet_lcc$ocount,
                            gcm1 = props1$graphlet_lcc$gcm,
                            gcm2 = props2$graphlet_lcc$gcm,
                            pvals1 = props1$graphlet_lcc$pvals,
                            pvals2 = props2$graphlet_lcc$pvals,
                            pAdjust1 = props1$graphlet_lcc$pAdjust,
                            pAdjust2 = props2$graphlet_lcc$pAdjust)
  }
  
  output$paramsProperties <- list(centrLCC = centrLCC,
                                  avDissIgnoreInf = avDissIgnoreInf,
                                  sPathNorm = sPathNorm,
                                  sPathAlgo = sPathAlgo,
                                  normNatConnect = normNatConnect,
                                  connectivity = connectivity,
                                  clustMethod = clustMethod,
                                  clustPar = clustPar,
                                  clustPar2 = clustPar2,
                                  weightClustCoef = weightClustCoef,
                                  hubPar = hubPar,
                                  hubQuant = hubQuant,
                                  lnormFit = lnormFit,
                                  weightDeg = weightDeg,
                                  normDeg = normDeg,
                                  normBetw = normBetw,
                                  normClose = normClose,
                                  normEigen = normEigen)
  
  output$paramsNetConstruct <- x$parameters
  
  output$input <- list(assoMat1 = x$assoMat1,
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
                       countsJoint = x$countsJoint,
                       normCounts1 = x$normCounts1,
                       normCounts2 = x$normCounts2,
                       twoNets = twoNets,
                       groups = groups,
                       matchDesign = x$matchDesign,
                       assoType = x$assoType,
                       softThreshPower = x$softThreshPower,
                       sampleSize = x$sampleSize,
                       weighted = x$parameters$weighted, 
                       call = x$call)
  
  output$isempty <- list(isempty1 = isempty1, isempty2 = isempty2)
  
  class(output) <- "microNetProps"
  return(output)
}


