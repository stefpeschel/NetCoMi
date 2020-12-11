#' @title Calculate network properties
#'
#' @description Calculates network properterties for a given adjacency matrix
#'
#' @param adjaMat adjacency matrix
#' @param dissMat dissimilarity matrix
#' @param assoMat association matrix
#' @param avDissIgnoreInf logical indicating whether to ignore infinities when 
#'   calculating the average dissimilarity. If \code{FALSE}, infinity values 
#'   are set to 1.
#' @param sPathNorm logical. If TRUE, shortest paths are normalized by 
#'   average dissimilarity.
#' @param sPathAlgo character indicating the algorithm used for shortest path
#'   calculation.
#' @param connectivity logical indicating whether edge and vertex connectivity 
#'   should be calculated. Might be disabled to reduce execution time.
#' @param normNatConnect logical indicating whether natural connectivity should
#'   be normalized.
#' @param weighted indicated whether the network is weighted
#' @param isempty indicator whether the network contains any edges.
#' @param clustMethod character indicating the clustering algorithm.
#' @param clustPar list with parameters passed to the clustering functions.
#' @param weightClustCoef logical. If \code{TRUE}, the weighted clustering
#'   coefficient is used.
#' @param hubPar character vector with one or more centrality measures used
#'   for identifying hub nodes. Possible values are \code{degree},
#'   \code{betweenness}, \code{closeness}, and \code{eigenvector}.
#' @param hubQuant quantile used for determining hub nodes.
#' @param lnormFit hubs are nodes with a centrality value above the 95\%
#'   quantile of the fitted log-normal distribution (if \code{lnormFit = TRUE})
#'   or of the empirical distribution of centrality values.
#' @param weightDeg if \code{TRUE}, the weighted degree is used.
#' @param normDeg,normBetw,normClose,normEigen if \code{TRUE}, a normalized
#'   version of the respective centrality values is returned. By default, all
#'   centralities are normalized.
#' @param centrLCC logical indicating whether centralities should only be 
#'   computed for the largest connected component (LCC). If \code{TRUE} 
#'   (default), centrality values of disconnected components are zero. 
#' @param jaccard shall the Jaccard index be calculated?
#' @param jaccQuant quantile for the Jaccard index
#' @param verbose integer indicating the level of verbosity. Possible values:
#'   \code{"0"}: no messages, \code{"1"}: only important messages,
#'   \code{"2"}(default): all progress messages are shown. Can also be logical.
#'
#' @importFrom igraph graph_from_adjacency_matrix decompose.graph
#' @importFrom stats hclust as.dist cutree qlnorm quantile
#' @importFrom pulsar natural.connectivity


calc_props <- function(adjaMat, dissMat, assoMat, avDissIgnoreInf,
                       sPathNorm, sPathAlgo, connectivity, normNatConnect, 
                       weighted, isempty, clustMethod, clustPar, 
                       weightClustCoef, hubPar, hubQuant, lnormFit, 
                       weightDeg, normDeg, normBetw, normClose, normEigen,
                       centrLCC, jaccard = FALSE, jaccQuant = NULL, 
                       verbose = 0){
  
  if(isempty){
    output <- list(clust = NULL, tree = NULL, deg = 0, deg_unnorm = 0,
                   betw = 0, betw_unnorm = 0, close = 0, close_unnorm = 0,
                   eigen = 0, eigen_unnorm = 0, lccSize = 0,
                   avDiss = 0, avDiss_lcc = 0, avPath = 0, avPath_lcc = 0,
                   vertconnect = 0, vertconnect_lcc = 0, edgeconnect = 0, 
                   edgeconnect_lcc = 0, clustCoef = 0, clustCoef_lcc = 0,
                   density = 0, density_lcc = 0, modul = 0, modul_lcc = 0,
                   pep = 0, pep_lcc = 0,  hubs = NULL, topdeg = NULL, 
                   topbetw = NULL, topclose = NULL, topeigen = NULL)
    return(output)
  }

  #== Create igraph objects and decompose graph ================================
  # Create graph from adjacency matrix

  if(verbose == 2){
    message("Create igraph objects ... ", appendLF = FALSE)
  }
  
  net <- igraph::graph_from_adjacency_matrix(adjaMat, weighted=T,
                                             mode="undirected", diag=F)
  
  # decomposed graph
  dg_net <- igraph::decompose.graph(net)
  
  nNodes <- ncol(adjaMat)

  # size and number of the connected components
  dgcount <- unlist(lapply(dg_net, igraph::vcount))
  nComp <- length(dgcount)
  compSize <- t(as.matrix(table(dgcount)))
  compSize <- rbind(as.numeric(colnames(compSize)), compSize[1, ])
  dimnames(compSize) <- list(c("size:", "   #:"), rep("", ncol(compSize)))
  compSize <- compSize[, order(compSize[1,], decreasing = TRUE), drop = FALSE]
  
  # Largest connected component (LCC)
  indLCC <- which.max(unlist(lapply(dg_net, function(x) length(igraph::V(x)))))
  
  net_lcc <- dg_net[[indLCC]]
  lccSize <- length(igraph::V(net_lcc))
  lccSizeRel <- lccSize / nNodes
  
  # Adjacency of the LCC
  adjaMat_lcc <- as.matrix(igraph::as_adjacency_matrix(net_lcc, attr="weight"))
  diag(adjaMat_lcc) <- 1
  
  # Names of nodes in the LCC
  lccNames <- colnames(adjaMat_lcc)
  
  if(!weighted){
    dissMat[!is.infinite(dissMat)] <- 1
    diag(dissMat) <- 0
  }
  
  # Dissimilarity of the LCC
  dissMat_lcc <- dissMat[lccNames, lccNames]
  
  #== Graph objects of dissimilarity matrices ==================================
  
  dissMatnoInf <- dissMat
  dissMatnoInf[is.infinite(dissMatnoInf)] <- 0
  
  dissMatnoInf_lcc <- dissMat_lcc
  dissMatnoInf_lcc[is.infinite(dissMatnoInf_lcc)] <- 0
  
  # Whole network
  dissnet <- igraph::graph_from_adjacency_matrix(dissMatnoInf, 
                                                 weighted=T,
                                                 mode="undirected", 
                                                 diag=F)
  
  # LCC
  dissnet_lcc <- igraph::graph_from_adjacency_matrix(dissMatnoInf_lcc, 
                                                    weighted=T,
                                                    mode="undirected", 
                                                    diag=F)
  
  if(verbose == 2){
    message("Done.")
  }
  
  #== Clustering ============================================================

  clust <- clust_lcc <- NULL
  tree <- tree_lcc <- NULL
  
  if(clustMethod != "none"){
    
    if(verbose == 2){
      message("Compute clustering (", clustMethod, ") ... ", appendLF = FALSE)
    }
    
    if(clustMethod == "hierarchical"){

      if(is.null(clustPar$method)){
        clustPar$method <- "average"
      }

      dissMat.tmp <- dissMat
      dissMat.tmp[is.infinite(dissMat.tmp)] <- 1
      
      dissMat_lcc.tmp <- dissMat_lcc
      dissMat_lcc.tmp[is.infinite(dissMat_lcc.tmp)] <- 1
      
      tree <- stats::hclust(stats::as.dist(dissMat.tmp), 
                            method = clustPar$method)

      tree_lcc <- stats::hclust(as.dist(dissMat_lcc.tmp), 
                                method = clustPar$method)
      rm(dissMat.tmp, dissMat_lcc.tmp)
      
      if(is.null(clustPar$k) & is.null(clustPar$h)){
        clustPar$k <- 3
      }
      
      clust <- do.call(stats::cutree, list(tree = tree, 
                                           k = clustPar$k, 
                                           h = clustPar$h))
      
      clust_lcc <- do.call(stats::cutree, list(tree = tree_lcc, 
                                               k = clustPar$k, 
                                               h = clustPar$h))
      
      names(clust) <- rownames(adjaMat)
      names(clust_lcc) <- rownames(adjaMat_lcc)

    } else{

      if(clustMethod == "cluster_edge_betweenness"){
        clustres <- do.call(getExportedValue("igraph", clustMethod),  
                            c(list(graph = net, weights = E(dissnet)$weight),
                              clustPar))
        
        clustres_lcc <- do.call(getExportedValue("igraph", clustMethod),  
                            c(list(graph = net_lcc, 
                                   weights = E(dissnet_lcc)$weight),
                              clustPar))
      } else{
        clustres <- do.call(getExportedValue("igraph", clustMethod), 
                            c(list(net), clustPar))
        clustres_lcc <- do.call(getExportedValue("igraph", clustMethod), 
                            c(list(net_lcc), clustPar))
      }

      clust <- clustres$membership
      names(clust) <- clustres$names

      # LCC
      clust_lcc <- clustres_lcc$membership
      names(clust_lcc) <- clustres_lcc$names
    }
    
    cltab <- table(clust)
    cltab_lcc <- table(clust_lcc)
    
    # cluster 0 assigned to elements in a cluster with only one element
    clust[clust %in% which(cltab == 1)] <- 0
    clust_lcc[clust_lcc %in% which(cltab_lcc == 1)] <- 0
    
    if(verbose == 2){
      message("Done.")
    }
    
  }
  

  #== Shortest paths ===========================================================
  
  if(verbose == 2){
    message("Compute shortest paths ... ", appendLF = FALSE)
  }
  
  # Whole network
  sPath <- igraph::distances(dissnet, algorithm = sPathAlgo)
  
  # LCC
  sPath_lcc <- igraph::distances(dissnet_lcc, algorithm = sPathAlgo)
  
  if(verbose == 2){
    message("Done.")
  }
  
  #== Average dissimilarity/distance ===========================================

  if(weighted){
    dissVec <- dissMat[lower.tri(dissMat)]
    dissVec_lcc <- dissMat_lcc[lower.tri(dissMat_lcc)]
    
    if(avDissIgnoreInf){
      dissVec[is.infinite(dissVec)] <- NA
      dissVec_lcc[is.infinite(dissVec_lcc)] <- NA
    } else{
      dissVec[is.infinite(dissVec)] <- 1
      dissVec_lcc[is.infinite(dissVec_lcc)] <- 1
    }
    
    ### average dissimilarity
    avDiss <- mean(dissVec, na.rm = TRUE)
    avDiss_lcc <- mean(dissVec_lcc, na.rm = TRUE)
    
    # "normalized" shortest paths (units with average dissimilarity)
    if(sPathNorm){
      sPath <- sPath / avDiss
      sPath_lcc <- sPath_lcc / avDiss_lcc
    }
  } else{
    avDiss <- avDiss_lcc <- 1
  }
  
  #== global network properties ================================================
  
  if(verbose == 2){
    message("\nCompute global properties:")
    
    message("Average dissimilarity ... ", appendLF = FALSE)
  }
  
  ### Average shortest path length 
  
  sPathVec <- sPath[lower.tri(sPath)]
  sPathVec[is.infinite(sPathVec)] <- NA
  avPath <- mean(sPathVec, na.rm = TRUE)
  if(is.na(avPath)) avPath <- 0
  
  sPathVec_lcc <- sPath_lcc[lower.tri(sPath_lcc)]
  avPath_lcc <- mean(sPathVec_lcc, na.rm = TRUE)
  if(is.na(avPath_lcc)) avPath_lcc <- 0
  
  if(verbose == 2){
    message("Done.")
  }
  
  #-------------------------------
  ### Connectivity 
  
  if(verbose == 2){
    message("Edge / vertex connectivity ... ", appendLF = FALSE)
  }
  
  if(connectivity){
    # vertex connectivity
    vertconnect <- igraph::vertex_connectivity(net)
    vertconnect_lcc <- igraph::vertex_connectivity(net_lcc)
    
    # edge connectivity
    edgeconnect <- igraph::edge_connectivity(net)
    edgeconnect_lcc <- igraph::edge_connectivity(net_lcc)
  } else{
    vertconnect <- vertconnect_lcc <- NA
    edgeconnect <- edgeconnect_lcc <- NA
  }
  
  if(verbose == 2){
    message("Done.")
  }
  
  #-------------------------------
  ### Natural connectivity
  
  if(verbose == 2){
    message("Natural connectivity ... ", appendLF = FALSE)
  }
  
  natConnect <- pulsar::natural.connectivity(adjaMat, norm = normNatConnect)
  natConnect_lcc <- pulsar::natural.connectivity(adjaMat_lcc,
                                                 norm = normNatConnect)
  
  if(verbose == 2){
    message("Done.")
  }
  
  #-------------------------------
  ### Clustering coefficient

  if(verbose == 2){
    message("Clustering coefficient ... ", appendLF = FALSE)
  }
  
  if(weightClustCoef){
    # Complete network
    clustCoef.tmp <- igraph::transitivity(net, type = "barrat")
    clustCoef <- mean(clustCoef.tmp, na.rm = TRUE)

    # LCC
    clustCoef.tmp <- igraph::transitivity(net_lcc, type = "barrat")
    clustCoef_lcc <- mean(clustCoef.tmp, na.rm = TRUE)
    
    rm(clustCoef.tmp)
  } else{
    # Complete network
    clustCoef <- igraph::transitivity(net, type = "global")
    if(is.na(clustCoef)) clustCoef <- 0
    
    # LCC
    clustCoef_lcc <- igraph::transitivity(net_lcc, type = "global")
    if(is.na(clustCoef_lcc)) clustCoef_lcc <- 0
  }

  
  if(verbose == 2){
    message("Done.")
  }
  
  #-------------------------------
  ### Modularity
  
  if(verbose == 2){
    message("Modularity ... ", appendLF = FALSE)
  }
  
  # Complete network
  if(clustMethod != "none"){
    modul <- igraph::modularity(net, (clust+1))
  } else{
    modul <- NA
  }
  
  # LCC
  if(clustMethod != "none"){
    modul_lcc <- igraph::modularity(net_lcc, (clust_lcc+1))
  } else{
    modul_lcc <- NA
  }
  
  if(verbose == 2){
    message("Done.")
  }
  
  #-------------------------------
  ### Density (relative number of edges)
  
  if(verbose == 2){
    message("Density ... ", appendLF = FALSE)
  }
  
  density <- igraph::edge_density(net)
  density_lcc <- igraph::edge_density(net_lcc)
  
  if(verbose == 2){
    message("Done.")
  }
  #-------------------------------
  ### Ratio of positive to negative associations

  if(verbose == 2){
    message("P-N-Ratio ... ", appendLF = FALSE)
  }
  
  # Complete network
  if(is.null(assoMat)){ # no negative dissimilarities possible
    pep <- 100
    pep_lcc <- 100
    
  } else{
    edge_all <- sum(assoMat[lower.tri(assoMat)] != 0)
    edge_pos <- sum(assoMat[lower.tri(assoMat)] > 0)
    pep <- edge_pos / edge_all * 100

    # LCC
    assoMat_lcc <- assoMat[rownames(dissMat_lcc), colnames(dissMat_lcc)]
    
    edge_all <- sum(assoMat_lcc[lower.tri(assoMat_lcc)] != 0)
    edge_pos <- sum(assoMat_lcc[lower.tri(assoMat_lcc)] > 0)
    pep_lcc <- edge_pos / edge_all * 100
  }

  if(verbose == 2){
    message("Done.")
  }
  
  #== Centrality measures ======================================================

  if(verbose == 2){
    message("\nCompute centralities:")
    
    message("Degree ... ", appendLF = FALSE)
  }
  
  ### degree
  if(weightDeg){
    deg <- deg_unnorm <- igraph::strength(net)
    deg_lcc <- strength(net_lcc)

  } else{
    deg <- igraph::degree(net, normalized = normDeg)
    deg_unnorm <- igraph::degree(net)
  }
  
  if(centrLCC){
    deg[!names(deg) %in% lccNames] <- 0
    deg_unnorm[!names(deg_unnorm) %in% lccNames] <- 0
  }

  if(verbose == 2){
    message("Done.")
  }
  
  #-------------------------------
  ### betweenness centrality (based on distances)

  if(verbose == 2){
    message("Betweenness centrality ... ", appendLF = FALSE)
  }
  
  if(centrLCC){
    betw <- betw_unnorm <- rep(0, ncol(adjaMat))
    names(betw) <- names(betw_unnorm) <- colnames(adjaMat)
    
    betw.tmp <- igraph::betweenness(dissnet_lcc, normalized = normBetw)
    betw_unnorm.tmp <- igraph::betweenness(dissnet_lcc)
    
    betw[names(betw.tmp)] <- betw.tmp
    betw_unnorm[names(betw_unnorm.tmp)] <- betw_unnorm.tmp
    
  } else{
    betw <- igraph::betweenness(dissnet, normalized = normBetw)
    betw_unnorm <- igraph::betweenness(dissnet)
  }
  
  betw[is.nan(betw)] <- 0
  betw_unnorm[is.nan(betw_unnorm)] <- 0

  if(verbose == 2){
    message("Done.")
  }
  
  #-------------------------------
  ### closeness centrality (based on distances)

  if(verbose == 2){
    message("Closeness centrality ... ", appendLF = FALSE)
  }
  
  if(centrLCC){
    close_unnorm <- rep(0, ncol(adjaMat))
    names(close_unnorm) <- colnames(adjaMat)
    
    sPath_rev <- 1/sPath_lcc
    diag(sPath_rev) <- NA
    
    close_unnorm.tmp <- sapply(1:nrow(sPath_rev), function(i){
      sum(sPath_rev[i,], na.rm = TRUE)
    })
    names(close_unnorm.tmp) <- lccNames
    close_unnorm[lccNames] <- close_unnorm.tmp
    
    if(normClose){
      # normalize closeness by n-1 
      close <- close_unnorm / (lccSize - 1)
    } else{
      close <- close_unnorm
    }

  } else{
    #if(nComp > 1 && sPathNorm){
    #warning("Normalized shortest paths depend on average dissimilarity, which
    #may be biased because infinite paths between unconnected nodes are ignored. 
    #Hence, closeness centrality may also be biased and should either be based 
    #on non-normalized shortest paths or be computed only for the LCC.")
    #}
    
    sPath_rev <- 1/sPath
    diag(sPath_rev) <- NA
    
    close_unnorm <- sapply(1:nrow(sPath_rev), function(i){
      sum(sPath_rev[i,], na.rm = TRUE)
    })
    names(close_unnorm) <- colnames(adjaMat)
    
    if(normClose){
      # normalize closeness by n-1 
      close <- close_unnorm / (nNodes - 1)
    } else{
      close <- close_unnorm
    }
  }

  if(verbose == 2){
    message("Done.")
  }
  
  #-------------------------------
  ### Eigenvector centrality

  if(verbose == 2){
    message("Eigenvector centrality ... ", appendLF = FALSE)
  }
  
  if(!centrLCC && nComp > 1){
    
    dgcount[dgcount == 1] <- 0
    
    vnumb <- sum(unlist(dgcount))

    ev <- numeric(0)
    
    for(i in seq_along(dg_net)){
      ev <- c(ev, 
              igraph::eigen_centrality(dg_net[[i]], 
                                       scale = FALSE)$vector * (dgcount[i] / vnumb))
    }

    eigen_unnorm <- ev[colnames(adjaMat)]

    if(normEigen){
      eigen <- eigen_unnorm / max(eigen_unnorm)
    } else{
      eigen <- eigen_unnorm
    }

  } else{
    eigen.tmp <- igraph::eigen_centrality(net_lcc, scale = normEigen)$vector
    eigen_unnorm.tmp <- igraph::eigen_centrality(net_lcc, scale = FALSE)$vector
    
    eigen <- eigen_unnorm <- numeric(nNodes)
    names(eigen) <- names(eigen_unnorm) <- colnames(adjaMat)
    eigen[names(eigen.tmp)] <- eigen.tmp
    eigen_unnorm[names(eigen_unnorm.tmp)] <- eigen_unnorm.tmp
  }
  
  if(verbose == 2){
    message("Done.")
  }
  
  #== hubs and Jaccard index ================================================

  if(verbose == 2){
    if(lnormFit){
      message("\nCompute hubs (based on log-normal quantiles... ", appendLF = FALSE)
    } else{
      message("\nCompute hubs (based on empirical quantiles) ... ", appendLF = FALSE)
    }
    
  }
  
  if(lnormFit){
    # identify nodes with highest centrality value
    pdeg <- try(MASS::fitdistr(deg_unnorm[deg_unnorm>0], "lognormal")$estimate, 
                silent = TRUE)
    
    if(class(pdeg) == "try-error"){
      topdeg <- hubdeg <- NULL
      
    } else{
      hubdeg <- names(deg_unnorm[deg_unnorm > stats::qlnorm(hubQuant, 
                                                            pdeg[1], 
                                                            pdeg[2])])
      
      if(jaccard){
        topdeg <- names(deg_unnorm[deg_unnorm > stats::qlnorm(jaccQuant, 
                                                              pdeg[1], 
                                                              pdeg[2])])
      } else{
        topdeg <- NULL
      }
    }

    pbetw <- try(MASS::fitdistr(betw_unnorm[betw_unnorm>0], "lognormal")$estimate, 
                 silent = TRUE)
    
    if(class(pbetw) == "try-error"){
      topbetw <- hubbetw <- NULL
    } else{
      hubbetw <- names(betw_unnorm[betw_unnorm > stats::qlnorm(hubQuant, 
                                                               pbetw[1], 
                                                               pbetw[2])])
      if(jaccard){
        topbetw <- names(betw_unnorm[betw_unnorm > stats::qlnorm(jaccQuant, 
                                                                 pbetw[1], 
                                                                 pbetw[2])])
      } else{
        topbetw <- NULL
      }
    }

    pclose <- try(MASS::fitdistr(close_unnorm[close_unnorm>0], "lognormal")$estimate, 
                  silent = TRUE)
    
    if(class(pclose) == "try-error"){
      topclose <- hubclose <- NULL
    } else{
      hubclose <- names(close_unnorm[close_unnorm > stats::qlnorm(hubQuant, 
                                                                  pclose[1], 
                                                                  pclose[2])])
      if(jaccard){
        topclose <- names(close_unnorm[close_unnorm > stats::qlnorm(jaccQuant, 
                                                                    pclose[1], 
                                                                    pclose[2])])
      } else{
        topclose <- NULL
      }
    }

    peigen <- try(MASS::fitdistr(eigen_unnorm[eigen_unnorm>0], "lognormal")$estimate, 
                  silent = TRUE)
    
    if(class(peigen) == "try-error"){
      topeigen <- hubeigen <- NULL
    } else{
      hubeigen <- names(eigen_unnorm[eigen_unnorm > stats::qlnorm(hubQuant, 
                                                                  peigen[1], 
                                                                  peigen[2])])
      if(jaccard){
        topeigen <- names(eigen_unnorm[eigen_unnorm > stats::qlnorm(jaccQuant, 
                                                                    peigen[1], 
                                                                    peigen[2])])
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
  
  if(verbose == 2){
    message("Done.")
  }

  #========================================================================

  output <- list(nComp = nComp, compSize= compSize, lccNames = lccNames,
                 clust = clust, tree = tree, 
                 clust_lcc = clust_lcc, tree_lcc = tree_lcc,
                 deg = deg, deg_unnorm = deg_unnorm,
                 betw = betw, betw_unnorm = betw_unnorm,
                 close = close, close_unnorm = close_unnorm,
                 eigen = eigen, eigen_unnorm = eigen_unnorm,
                 lccSize = lccSize, lccSizeRel = lccSizeRel,
                 avDiss = avDiss, avDiss_lcc = avDiss_lcc,
                 avPath = avPath, avPath_lcc = avPath_lcc,
                 vertconnect = vertconnect, vertconnect_lcc = vertconnect_lcc,
                 edgeconnect = edgeconnect, edgeconnect_lcc = edgeconnect_lcc,
                 natConnect = natConnect, natConnect_lcc = natConnect_lcc,
                 clustCoef = clustCoef, clustCoef_lcc = clustCoef_lcc,
                 density = density, density_lcc = density_lcc,
                 modul = modul, modul_lcc = modul_lcc,
                 pep = pep, pep_lcc = pep_lcc,
                 hubs = hubs, 
                 topdeg = topdeg, topbetw = topbetw, 
                 topclose = topclose, topeigen = topeigen)
  
  return(output)
}

