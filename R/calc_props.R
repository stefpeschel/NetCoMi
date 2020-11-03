#' @title Calculate network properties
#'
#' @description Calculates network properterties for a given adjacency matrix
#'
#' @param adjaMat adjacency matrix
#' @param dissMat dissimilarity matrix
#' @param weighted indicated whether the network is weighted
#' @param isempty indicator whether the network contains any edges.
#' @param clustMethod character indicating the clustering algorithm.
#' @param clustPar list with parameters passed to the clustering functions.
#' @param hubPar character vector with one or more centrality measures used
#'   for identifying hub nodes. Possible values are \code{degree},
#'   \code{betweenness}, \code{closeness}, and \code{eigenvector}.
#' @param hubQuant quantile used for determining hub nodes.
#' @param lnormFit hubs are nodes with a centrality value above the 95\%
#'   quantile of the fitted log-normal distribution (if \code{lnormFit = TRUE})
#'   or of the empirical distribution of centrality values.
#' @param connect logical indicating whether edge and vertex connectivity should
#'   be calculated. Might be disabled to reduce execution time.
#' @param weightDeg if \code{TRUE}, the weighted degree is used.
#' @param normDeg,normBetw,normClose,normEigen if \code{TRUE}, a normalized
#'   version of the respective centrality values is returned. By default, all
#'   centralities are normalized.
#' @param jaccard shall the Jaccard index be calculated?
#' @param jaccQuant quantile for the Jaccard index
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom stats hclust as.dist cutree qlnorm quantile


calc_props <- function(adjaMat, dissMat, weighted, isempty, clustMethod, clustPar,
                      hubPar, hubQuant, lnormFit, connect,
                      weightDeg, normDeg, normBetw, normClose, normEigen,
                      jaccard = FALSE, jaccQuant = NULL){
  if(isempty){
    output <- list(clust = NULL, hubs = NULL, deg = 0, betw = 0, close = 0,
                   eigen = 0, avPath = 0, clustCoef = 0, modul = 0,
                   vertconnect = 0, edgeconnect = 0, density = 0)
    return(output)
  }

  #== create igraph objects =================================================
  # create network from adjacency matrix
  net <- igraph::graph_from_adjacency_matrix(adjaMat, weighted=T,
                                             mode="undirected", diag=F)

  dissMat_noinf <- dissMat
  dissMat_noinf[is.infinite(dissMat_noinf)] <- 0
  distnet <- igraph::graph_from_adjacency_matrix(dissMat_noinf, weighted=T,
                                                 mode="undirected", diag=F)

  #== clustering ============================================================

  clust <- NULL
  tree <- NULL
  if(clustMethod != "none"){
    if(clustMethod == "hierarchical"){

      if(is.null(clustPar$method)){
        clustPar$method <- "average"
      }

      dissMat.tmp <- dissMat
      dissMat.tmp[is.infinite(dissMat.tmp)] <- 1
      tree <- stats::hclust(stats::as.dist(dissMat.tmp), method = clustPar$method)

      if(is.null(clustPar$k) & is.null(clustPar$h)){
        clustPar$k <- 3
      }
      clust <- do.call(stats::cutree, list(tree = tree, 
                                           k = clustPar$k, 
                                           h = clustPar$h))
      names(clust) <- rownames(adjaMat)

    } else{

      if(clustMethod == "cluster_edge_betweenness"){
        clustres <- do.call(clustMethod,  c(list(graph = net,
                                                 weights = E(distnet)$weight),
                                            clustPar))
      } else{
        clustres <- do.call(clustMethod, c(list(net), clustPar))
      }

      clust <- clustres$membership
      names(clust) <- clustres$names

      cltab <- table(clust)

      # cluster 0 assigned to elements in a cluster with only one element
      clust[clust %in% which(cltab == 1)] <- 0
    }

  }

  #== centrality measures ===================================================

  ### degree
  if(weightDeg){
    deg <- deg_unnorm <- igraph::strength(net)

  } else{
    deg <- igraph::degree(net, normalized = normDeg)
    deg_unnorm <- igraph::degree(net)
  }

  #-------------------------------
  ### betweenness centrality (based on distances)

  betw <- igraph::betweenness(distnet, normalized = normBetw)
  betw_unnorm <- igraph::betweenness(distnet)


  #-------------------------------
  ### closeness centrality (based on distances)

  # compute shortest paths
  sPath <- igraph::distances(distnet, algorithm = "dijkstra")
  sPath_rev <- 1/sPath
  diag(sPath_rev) <- NA

  close_unnorm <- sapply(1:nrow(sPath_rev), function(i) sum(sPath_rev[i,], na.rm = TRUE))
  names(close_unnorm) <- colnames(adjaMat)

  if(normClose){
    # normalize closeness by n-1 (unconnected nodes are ignored)
    vnumb_connected <- sum(close_unnorm != 0)
    close <- close_unnorm / (vnumb_connected - 1)

  } else{
    close <- close_unnorm
  }

  #-------------------------------
  ### Eigenvector centrality

  dg <- decompose.graph(net)
  if(length(dg) > 1){
    dgcount <- unlist(lapply(dg, vcount))
    dgcount[dgcount == 1] <- 0
    vnumb <- sum(unlist(dgcount))

    ev <- numeric(0)
    for(i in seq_along(dg)){
      ev <- c(ev, igraph::eigen_centrality(dg[[i]], scale = FALSE)$vector * 
                (dgcount[i] / vnumb))
    }

    eigen_unnorm <- ev[colnames(adjaMat)]

    if(normEigen){
      eigen <- eigen_unnorm / max(eigen_unnorm)
    } else{
      eigen <- eigen_unnorm
    }

  } else{
    eigen <- igraph::eigen_centrality(net, scale = normEigen)$vector
    eigen_unnorm <- igraph::eigen_centrality(net, scale = FALSE)$vector
  }

  #-------------------------------
  ### global network properties

  # average path length
  sPath[is.infinite(sPath)] <- NA

  avPath <- mean(sPath, na.rm = TRUE)
  if(is.na(avPath)) avPath <- 0

  # clustering coefficient
  clustCoef <- igraph::transitivity(net, type = "global")
  if(is.na(clustCoef)) clustCoef <- 0

  # modularity
  if(clustMethod != "none"){
    modul <- igraph::modularity(net, (clust+1))
  } else{
    modul <- NA
  }


  if(connect){
    # vertex connectivity
    vertconnect <- igraph::vertex_connectivity(net)

    # edge connectivity
    edgeconnect <- igraph::edge_connectivity(net)
  } else{
    vertconnect <- NA
    edgeconnect <- NA
  }


  # relative number of edges(density)
  density <- igraph::edge_density(net)

  #== hubs and Jaccard index ================================================


  if(lnormFit){
    # identify nodes with highest centrality value
    pdeg <- try(MASS::fitdistr(deg_unnorm[deg_unnorm>0], "lognormal")$estimate, silent = TRUE)
    
    if(class(pdeg) == "try-error"){
      topdeg <- hubdeg <- NULL
      
    } else{
      hubdeg <- names(deg_unnorm[deg_unnorm > stats::qlnorm(hubQuant, pdeg[1], pdeg[2])])
      
      if(jaccard){
        topdeg <- names(deg_unnorm[deg_unnorm > stats::qlnorm(jaccQuant, pdeg[1], pdeg[2])])
      } else{
        topdeg <- NULL
      }
    }

    pbetw <- try(MASS::fitdistr(betw_unnorm[betw_unnorm>0], "lognormal")$estimate, silent = TRUE)
    
    if(class(pbetw) == "try-error"){
      topbetw <- hubbetw <- NULL
    } else{
      hubbetw <- names(betw_unnorm[betw_unnorm > stats::qlnorm(hubQuant, pbetw[1], pbetw[2])])
      if(jaccard){
        topbetw <- names(betw_unnorm[betw_unnorm > stats::qlnorm(jaccQuant, pbetw[1], pbetw[2])])
      } else{
        topbetw <- NULL
      }
    }

    pclose <- try(MASS::fitdistr(close_unnorm[close_unnorm>0], "lognormal")$estimate, silent = TRUE)
    
    if(class(pclose) == "try-error"){
      topclose <- hubclose <- NULL
    } else{
      hubclose <- names(close_unnorm[close_unnorm > stats::qlnorm(hubQuant, pclose[1], pclose[2])])
      if(jaccard){
        topclose <- names(close_unnorm[close_unnorm > stats::qlnorm(jaccQuant, pclose[1], pclose[2])])
      } else{
        topclose <- NULL
      }
    }

    peigen <- try(MASS::fitdistr(eigen_unnorm[eigen_unnorm>0], "lognormal")$estimate, silent = TRUE)
    
    if(class(peigen) == "try-error"){
      topeigen <- hubeigen <- NULL
    } else{
      hubeigen <- names(eigen_unnorm[eigen_unnorm > stats::qlnorm(hubQuant, peigen[1], peigen[2])])
      if(jaccard){
        topeigen <- names(eigen_unnorm[eigen_unnorm > stats::qlnorm(jaccQuant, peigen[1], peigen[2])])
      } else{
        topeigen <- NULL
      }
    }

  } else{
    hubdeg <- names(deg[deg > quantile(deg, hubQuant)])
    if(jaccard){
      topdeg <- names(deg[deg > quantile(deg, jaccQuant)])
    } else{
      topdeg <- NULL
    }

    hubbetw <- names(betw[betw > quantile(betw, hubQuant)])
    if(jaccard){
      topbetw <- names(betw[betw > quantile(betw, jaccQuant)])
    } else{
      topbetw <- NULL
    }

    hubclose <- names(close[close > quantile(close, hubQuant)])
    if(jaccard){
      topclose <- names(close[close > quantile(close, jaccQuant)])
    } else{
      topclose <- NULL
    }

    hubeigen <- names(eigen[eigen > quantile(eigen, hubQuant)])
    if(jaccard){
      topeigen <- names(eigen[eigen > quantile(eigen, jaccQuant)])
    } else{
      topeigen <- NULL
    }
  }

  ###############
  # identify hub nodes
  hubparams <-  list()
  if("degree" %in% hubPar){
    hubparams <- c(hubparams, list(hubdeg))
  }
  if("betweenness" %in% hubPar){
    hubparams <- c(hubparams, list(hubbetw))
  }
  if("closeness" %in% hubPar){
    hubparams <- c(hubparams, list(hubclose))
  }
  if("eigenvector" %in% hubPar){
    hubparams <- c(hubparams, list(hubeigen))
  }

  hubs <- Reduce(intersect, hubparams)

  #========================================================================

  output <- list(clust = clust, tree = tree, deg = deg, deg_unnorm = deg_unnorm,
                 betw = betw, betw_unnorm = betw_unnorm,
                 close = close, close_unnorm = close_unnorm,
                 eigen = eigen, eigen_unnorm = eigen_unnorm,
                 avPath = avPath, clustCoef = clustCoef, modul = modul,
                 vertconnect = vertconnect, edgeconnect = edgeconnect,
                 density = density, hubs = hubs, topdeg = topdeg,
                 topbetw = topbetw, topclose = topclose, topeigen = topeigen)

}
