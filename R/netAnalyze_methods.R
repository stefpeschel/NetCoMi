#' @title Summary Method for Objects of Class microNetProps
#' @description  The main results returned by \code{\link{netAnalyze}} are
#'   printed in a well-arranged format.
#'
#' @param object object of class \code{microNetProps} (returned by 
#'   \code{\link{netAnalyze}}).
#' @param groupNames character vector with two elements giving the group names
#'   corresponding to the two networks. If \code{NULL}, the names are adopted
#'   from \code{object}. Ignored if \code{object} contains a single network.
#' @param showCompSize logical. If \code{TRUE}, the component sizes are printed.
#' @param showGlobal logical. If \code{TRUE}, global network properties for the 
#'   whole network are printed.
#' @param showGlobalLCC logical. If \code{TRUE}, global network properties for 
#'   the largest connected component are printed. If the network is connected 
#'   (number of components is 1) the global properties are printed only once 
#'   (if one of the arguments \code{showGlobal} and \code{showGlobalLCC}) is 
#'   \code{TRUE}.
#' @param showCluster logical. If \code{TRUE}, the cluster(s) are printed.
#' @param clusterLCC logical. If \code{TRUE}, clusters are printed only for the
#'   largest connected component. Defaults to \code{FALSE} (whole network).
#' @param showHubs logical. If \code{TRUE}, the detected hubs are printed.
#' @param showCentr character vector indicating for which centrality measures
#'   the results shall be printed. Possible values are "all", "degree",
#'   "betweenness", "closeness", "eigenvector" and "none".
#' @param numbNodes integer indicating for how many nodes the centrality
#'   values shall be printed. Defaults to 10L for a single network and 5L for 
#'   two networks. Thus, in case of a single network, the first 10 
#'   nodes with highest centrality value of the specific centrality measure
#'   are shown. If \code{object} contains two networks, for each centrality 
#'   measure a splitted matrix is shown where the upper part contains the 
#'   highest values of the first group, and the lower part the highest values of 
#'   the second group.
#' @param digits integer giving the number of decimal places to which the 
#'   results are rounded. Defaults to 5L.
#' @param ... not used.
#'
#' @seealso \code{\link{netConstruct}}, \code{\link{netAnalyze}}
#'
#' @method summary microNetProps
#' @rdname summarize.microNetProps
#' @export

summary.microNetProps <- function(object, 
                                  groupNames = NULL,
                                  showCompSize = TRUE,
                                  showGlobal = TRUE,
                                  showGlobalLCC = TRUE,
                                  showCluster = TRUE,
                                  clusterLCC = FALSE,
                                  showHubs = TRUE,
                                  showCentr = "all", 
                                  numbNodes = NULL,
                                  digits = 5L,
                                  ...) {
  
  #-----------------------------------------------------------------------------
  # Check input arguments
  
  if (!inherits(object, "microNetProps")) {
    stop('"object" must be of class "microNetProps".')
  }
  
  twonets <- object$input$twoNets
  
  if (twonets) {
    if (is.null(groupNames)) {
      group1 <- paste0("group '" , object$input$groups[1], "'")
      group2 <- paste0("group '" , object$input$groups[2], "'")
      
    } else {
      
      if (!length(groupNames) == 2) {
        stop('Length of "groupNames" must be 2.')
      }
      
      group1 <- groupNames[1]
      group2 <- groupNames[2]
    }
    
  } else {
    group1 <- group2 <- NULL
  }
  
  showCentr <- match.arg(showCentr, choices = c("all", "none", "degree",
                                                "betweenness", "closeness",
                                                "eigenvector"),
                         several.ok = TRUE)
  
  if ("none" %in% showCentr) stopifnot(length(showCentr) == 1)
  
  
  if (is.null(numbNodes)) {
    numbNodes <- ifelse(twonets, 5L, 10L)
    
  } else {
    stopifnot(is.numeric(numbNodes))
    stopifnot(numbNodes >= 1)
    
    numbNodes <- min(as.integer(numbNodes), length(object$centralities$degree1))
  }
  
  stopifnot(is.logical(showCompSize))
  stopifnot(is.logical(showGlobal))
  stopifnot(is.logical(showGlobalLCC))
  stopifnot(is.logical(showCluster))
  stopifnot(is.logical(clusterLCC))
  
  digits <- as.integer(digits)
  stopifnot(digits >= 0)
  
  #=============================================================================
  # Global network properties
  
  nComp1 <- object$globalProps$nComp1
  nComp2 <- ifelse(twonets, object$globalProps$nComp2, 1)
  
  if (nComp1 == 1 && nComp2 == 1) {
    is_disconnected <- FALSE
  } else {
    is_disconnected <- TRUE
  }
  
  if (!showGlobal & !showGlobalLCC) {
    showGlob <- showGlobLCC <- FALSE
    
  } else if (!is_disconnected) {
    showGlob <- showGlobLCC <- TRUE
    
  } else {
    showGlob <- showGlobal
    showGlobLCC <- showGlobalLCC
  }
  
  #-----------------------------------------------------------------------------
  # Whole network
  
  if (showGlob) {
    glob_rnames <- c("Number of components", 
                     "Clustering coefficient",
                     "Modularity",
                     "Positive edge percentage",
                     "Edge density",
                     "Natural connectivity")
    
    glob_names <- c("nComp", 
                    "clustCoef", 
                    "modularity", 
                    "pep", 
                    "density", 
                    "natConnect")
    
    if (is.na(object$globalProps$modularity1)) {
      # exclude modularity
      glob_rnames <- glob_rnames[-3]
      glob_names <- glob_names[-3]
    } 
    
    if (twonets) {
      glob_probs <-
        as.data.frame(matrix(0, nrow = length(glob_rnames),
                             ncol = 2, dimnames = list(glob_rnames,
                                                       c(group1, group2))
        ))
      
    } else {
      glob_probs <- 
        as.data.frame(matrix(0,  nrow = length(glob_rnames), 
                             ncol = 1, dimnames = list(glob_rnames, " ")))
    }
    
    for (i in 1:length(glob_names)) {
      glob_probs[i, 1] <- 
        round(as.numeric(object$globalProps[paste0(glob_names[i], 1)]), 
              digits = digits)
    }
    
    if (twonets) {
      for (i in 1:length(glob_names)) {
        glob_probs[i, 2] <- 
          round(as.numeric(object$globalProps[paste0(glob_names[i], 2)]), 
                digits = digits)
      }
    }
    
  } else {
    glob_probs <- NULL
  }
  
  #-----------------------------------------------------------------------------
  # LCC
  if (showGlobLCC) {
    glob_rnames_lcc <- c("Relative LCC size",
                         "Clustering coefficient",
                         "Modularity",
                         "Positive edge percentage",
                         "Edge density",
                         "Natural connectivity",
                         "Vertex connectivity",
                         "Edge connectivity",
                         "Average dissimilarity*",
                         "Average path length**")
    
    glob_names_lcc <-c("lccSizeRel",
                       "clustCoef",
                       "modularity",
                       "pep",
                       "density",
                       "natConnect",
                       "vertConnect",
                       "edgeConnect",
                       "avDiss",
                       "avPath")
    
    if (is.na(object$globalProps$modularity1)) {
      # exclude modularity
      glob_rnames_lcc <- glob_rnames_lcc[-3]
      glob_names_lcc <- glob_names_lcc[-3]
      
      if (is.na(object$globalPropsLCC$vertConnect1)) {
        # exclude connectivity measures
        glob_rnames_lcc <- glob_rnames_lcc[-c(6,7)]
        glob_names_lcc <- glob_names_lcc[-c(6,7)]
      } 
    } else if (is.na(object$globalPropsLCC$vertConnect1)) {
      # exclude connectivity measures
      glob_rnames_lcc <- glob_rnames_lcc[-c(7,8)]
      glob_names_lcc <- glob_names_lcc[-c(7,8)]
    }
    
    if (twonets) {
      glob_probs_lcc <- 
        as.data.frame(matrix(0, nrow = length(glob_rnames_lcc), 
                             ncol = 2, 
                             dimnames = list(glob_rnames_lcc,
                                             c(group1, group2))))
    } else {
      glob_probs_lcc <- 
        as.data.frame(matrix(0, nrow = length(glob_rnames_lcc), 
                             ncol = 1, 
                             dimnames = list(glob_rnames_lcc, " ")))
    }
    
    for (i in 1:length(glob_names_lcc)) {
      glob_probs_lcc[i, 1] <- round(as.numeric(
        object$globalPropsLCC[paste0(glob_names_lcc[i], 1)]), digits = digits)
    }
    
    if (twonets) {
      for (i in 1:length(glob_names_lcc)) {
        glob_probs_lcc[i, 2] <- round(as.numeric(
          object$globalPropsLCC[paste0(glob_names_lcc[i], 2)]), digits = digits)
      }
      
    }
    
  } else {
    glob_probs_lcc <- NULL
  }
  
  #-----------------------------------------------------------------------------
  # Combine for connected networks (with only one component)
  
  if (!is_disconnected & showGlob) {
    glob_probs <- rbind(glob_probs[1, , drop = FALSE], 
                        glob_probs_lcc[-1, , drop = FALSE])
  } 
  
  #=============================================================================
  # clustering
  
  if (showCluster) {
    clust <- list()
    if (clusterLCC) {
      clusttab1 <- table(object$clusteringLCC$clust1)
    } else {
      clusttab1 <- table(object$clustering$clust1)
    }
    
    clust1 <- matrix(0, nrow = 2, ncol = length(clusttab1), 
                     dimnames = list(c("name:", "   #:"),
                                     rep("", length(clusttab1))))
    clust1[1, ] <- as.numeric(names(clusttab1))
    clust1[2, ] <- clusttab1
    
    clust[[1]] <- clust1
    if (twonets) {
      if (clusterLCC) {
        clusttab2 <- table(object$clusteringLCC$clust2)
      } else {
        clusttab2 <- table(object$clustering$clust2)
      }
      
      clust2 <- matrix(0, nrow = 2, ncol = length(clusttab2), 
                       dimnames = list(c("name:", "   #:"),
                                       rep("", length(clusttab2))))
      clust2[1, ] <- as.numeric(names(clusttab2))
      clust2[2, ] <- clusttab2
      clust[[2]] <- clust2
    }
    
  } else {
    clust <- NULL
  }
  
  #=============================================================================
  # hubs
  
  if (showHubs) {
    if (twonets) {
      hubs1 <- sort(object$hubs$hubs1)
      hubs2 <- sort(object$hubs$hubs2)
      
      if (length(hubs1) != length(hubs2)) {
        diff <- length(hubs1) - length(hubs2)
        if (diff > 0) {
          hubs2 <- c(hubs2, rep("", diff))
        } else {
          hubs1 <- c(hubs1, rep("", abs(diff)))
        }
      } 
      
      hubmat <- cbind(hubs1, hubs2)
      
      dimnames(hubmat) <- list(rep("", nrow(hubmat)),
                               c(group1, group2))
      hubs <- as.data.frame(hubmat)
    } else {
      hubmat <- as.matrix(sort(object$hubs$hubs1))
      dimnames(hubmat) <- list(rep("", nrow(hubmat)), "")
      hubs <- hubmat
    }
    
  } else {
    hubs <- NULL
  }
  
  #=============================================================================
  # centrality measures
  
  top_centr <- NULL
  
  if (showCentr[1] != "none") {
    
    l <- length(object$centralities$degree1)
    
    centr_names1 <- c("degree1", "between1", "close1", "eigenv1")
    centr_names2 <- c("degree2", "between2", "close2", "eigenv2")
    centr_names <- c("degree", "betweenness", "closeness", "eigenvector")
    
    top_centr_names <- list()
    
    for (i in 1:length(centr_names1)) {
      top_centr_names[[centr_names1[i]]] <- 
        names(sort(object$centralities[[centr_names1[i]]], 
                   decreasing = TRUE))[1:min(numbNodes, l)]
    }
    
    if (twonets) {
      for (i in 1:length(centr_names2)) {
        top_centr_names[[centr_names2[i]]] <- 
          names(sort(object$centralities[[centr_names2[i]]], 
                     decreasing = TRUE))[1:min(numbNodes, l)]
      }
    }
    
    top_centr <- list()
    
    for (i in 1:length(centr_names)) {
      if (any(c("all", centr_names[i]) %in% showCentr)) {
        
        cn1 <- centr_names1[i]
        cn2 <- centr_names2[i]
        
        if (twonets) {
          mat1 <-
            as.data.frame(matrix(
              0,
              nrow = numbNodes ,
              ncol = 2,
              dimnames = list(top_centr_names[[cn1]],
                              c(group1, group2))
            ))
          
        } else {
          mat1 <-
            as.data.frame(matrix(
              0,
              nrow = numbNodes ,
              ncol = 1,
              dimnames = list(top_centr_names[[cn1]],
                              " ")
            ))
        }
        
        mat1[, 1] <-
          round(object$centralities[[cn1]][top_centr_names[[cn1]]],
                digits = digits)
        
        
        if (twonets) {
          mat1[, 2] <-
            round(object$centralities[[cn2]][top_centr_names[[cn1]]],
                  digits = digits)
          
          mat2 <-
            as.data.frame(matrix(
              0,
              nrow = numbNodes ,
              ncol = 2,
              dimnames = list(top_centr_names[[cn2]],
                              c(group1, group2))
            ))
          
          mat2[, 1] <-
            round(object$centralities[[cn1]][top_centr_names[[cn2]]],
                  digits = digits)
          
          mat2[, 2] <-
            round(object$centralities[[cn2]][top_centr_names[[cn2]]],
                  digits = digits)
          
          rnames <- c(rownames(mat1), "", rownames(mat2))
          
          mat.tmp <- rbind(mat1, c("______","______"), 
                           mat2, make.row.names = FALSE)
          
          mat <- cbind(rnames, mat.tmp)
          
          colnames(mat) <- c(" ", colnames(mat)[2:3])
        } else {
          mat <- mat1
        }
        
        top_centr[[cn1]] <- mat
      }
    }
  }
  
  #============================================================
  # Component sizes
  
  if (showCompSize) {
    if (twonets) {
      compSize <- list(object$compSize1, object$compSize2)
      names(compSize) <- c(group1, group2)
    } else {
      compSize <- object$compSize1
    }
  } else {
    compSize <- NULL
  }
  
  out <- list(glob_probs = glob_probs,
              glob_probs_lcc = glob_probs_lcc,
              clust = clust,
              hubs = hubs,
              central = top_centr,
              group1 = group1,
              group2 = group2,
              is_disconnected = is_disconnected,
              compSize = compSize,
              clusterLCC = clusterLCC,
              paramsProperties = object$paramsProperties,
              twoNets = object$input$twoNets,
              call = object$call)
  
  class(out) <-  "summary.microNetProps"
  
  return(out)
}


#' @title Print Summary for objects of class \code{microNetProps}
#'
#' @param x object of class \code{summary.microNetProps} (returned by 
#'   to \code{\link{summary.microNetProps}}).
#' @param ... not used.
#'
#' @method print summary.microNetProps
#'
#' @rdname summarize.microNetProps
#' @export
print.summary.microNetProps <- function(x, ...) {
  
  if (!is.null(x$compSize)) {
    cat("\nComponent sizes\n")
    
    if (is.list(x$compSize)) {
      cat("```````````````\n")
      cat(names(x$compSize)[1], ":", sep = "")
      print(x$compSize[[1]])
      cat(names(x$compSize)[2], ":", sep = "")
      print(x$compSize[[2]])
    } else {
      cat("```````````````")
      print(x$compSize)
    }
  }
  
  #-----------------------------------------------------------------------------
  # Global properties
  
  showGlob <- !(is.null(x$glob_probs) & is.null(x$glob_probs_lcc))
  
  if (showGlob) {
    cat("______________________________")
    cat("\nGlobal network properties\n")
    cat("`````````````````````````\n")
  }
  
  if (!is.null(x$glob_probs_lcc) & x$is_disconnected) {
    cat("Largest connected component (LCC):\n")
    print(x$glob_probs_lcc)
    
    if (!is.null(x$glob_probs)) {
      cat("\n")
    }
  }
  
  if (!is.null(x$glob_probs)) {
    cat("Whole network:\n")
    print(x$glob_probs)
  }
  
  if (showGlob) {
    cat("-----\n")
    cat("*: Dissimilarity = 1 - edge weight")
    
    if (x$paramsProperties$sPathNorm) {
      cat("\n**: Path length = Units with average dissimilarity\n")
    } else {
      cat("\n**: Path length = Sum of dissimilarities along the path\n")
    }
  }
  
  #-----------------------------------------------------------------------------
  # Cluster
  
  if (!is.null(x$clust)) {
    if (ncol(x$clust[[1]]) != 0) {
      cat("\n______________________________")
      if (x$is_disconnected & x$clusterLCC) {
        cat("\nClusters")
        cat("\n- In the LCC")
        cat(paste0("\n- Algorithm: ", x$paramsProperties$clustMethod, "\n"))
        cat(paste(replicate(13+nchar(x$paramsProperties$clustMethod), "`"), 
                  collapse = ""), "\n")
        
      } else {
        cat("\nClusters")
        cat("\n- In the whole network")
        cat(paste0("\n- Algorithm: ", x$paramsProperties$clustMethod, "\n"))
        cat(paste(replicate(13+nchar(x$paramsProperties$clustMethod), "`"), 
                  collapse = ""), "\n")
      } 
      
      if (length(x$clust) == 2) {
        cat(x$group1, ":", sep = "")
        print(x$clust[[1]])
        cat("\n", x$group2, ":", sep = "")
        print(x$clust[[2]])
      } else {
        print(x$clust[[1]])
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  # Hubs
  if (!is.null(x$hubs)) {
    cat("\n______________________________")
    cat("\nHubs\n")
    cat("- In alphabetical/numerical order")
    
    if (x$paramsProperties$lnormFit) {
      cat("\n- Based on log-normal quantiles of centralities\n")
    } else {
      cat("\n- Based on empirical quantiles of centralities\n")
    }
    
    if (nrow(x$hubs) == 0) {
      cat("```````````````````````````````````````````````\n")
      cat("No hubs detected.")
    } else if (ncol(x$hubs) == 2) {
      cat("```````````````````````````````````````````````\n")
      print(x$hubs, row.names = FALSE, quote = FALSE)
    } else {
      cat("```````````````````````````````````````````````")
      print(x$hubs, quote = FALSE)
    }
  }
  
  #-----------------------------------------------------------------------------
  # Centralities
  
  if (!is.null(x$central)) {
    show_rownames <- !x$twoNets
    
    cat("\n______________________________")
    
    if (x$is_disconnected && x$paramsProperties$centrLCC) {
      cat("\nCentrality measures")
      cat("\n- In decreasing order")
      cat("\n- Centrality of disconnected components is zero\n")
      cat("````````````````````````````````````````````````")
    } else {
      cat("\nCentrality measures")
      cat("\n- In decreasing order")
      cat("\n- Computed for the complete network\n")
      cat("````````````````````````````````````")
    }
    
    if (!is.null(x$central$degree1)) {
      
      if (x$paramsProperties$weightDeg) {
        cat("\nDegree (weighted):\n")
      } else if (x$paramsProperties$normDeg) {
        cat("\nDegree (normalized):\n")
      } else {
        cat("\nDegree (unnormalized):\n")
      }
      
      print(x$central$degree1, row.names = show_rownames)
    }
    
    if (!is.null(x$central$between1)) {
      
      if (x$paramsProperties$normBetw) {
        cat("\nBetweenness centrality (normalized):\n")
      } else {
        cat("\nBetweenness centrality (unnormalized):\n")
      }
      
      print(x$central$between1, row.names = show_rownames)
    }
    
    if (!is.null(x$central$close1)) {
      
      if (x$paramsProperties$normBetw) {
        cat("\nCloseness centrality (normalized):\n")
      } else {
        cat("\nCloseness centrality (unnormalized):\n")
      }
      
      print(x$central$close1, row.names = show_rownames)
    }
    
    if (!is.null(x$central$eigenv1)) {
      
      if (x$paramsProperties$normBetw) {
        cat("\nEigenvector centrality (normalized):\n")
      } else {
        cat("\nEigenvector centrality (unnormalized):\n")
      }
      
      print(x$central$eigenv1, row.names = show_rownames)
    }
    
  }
}
