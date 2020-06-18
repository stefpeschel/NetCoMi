#' @title Summary Method for Objects of Class microNetProps
#' @description  The main results returned from \code{\link{netAnalyze}} are
#'   printed in a well-arranged format.
#'
#' @param object object of class \code{microNetProps} inheriting from a call to
#'   \code{\link{netAnalyze}}.
#' @param groupNames character vector with two elements giving the group names
#'   corresponding to the two networks. If \code{NULL}, the names are adopted
#'   from \code{object}. Ignored if \code{object} contains a single network.
#' @param showCentr character vector indicating for which centrality measures
#'   the results shall be printed. Possible values are "all", "degree",
#'   "betweenness", "closeness", "eigenvector" and "none".
#' @param numbNodes integer indicating for how much taxa the centrality
#'   values shall be printed. Defaults to 10L for a single network and 5L for 
#'   two networks. That means that in case of a single network the first 10 
#'   nodes with highest centrality value of the specific centrality measure
#'   are shown. If \code{object} contains two networks, for each centrality 
#'   measure a splitted matrix is shown where the upper part contains the 
#'   highest values of the first group, and the lower part the highest values of 
#'   the second group.
#' @param digits integer giving the number of decimal places to which the 
#'   results are rounded. Defaults to 5.
#' @param ... not used.
#'
#' @seealso \code{\link{netConstruct}}, \code{\link{netAnalyze}}
#'
#' @method summary microNetProps
#' @rdname summarize.microNetProps
#' @export
summary.microNetProps <- function(object, groupNames = NULL, showCentr = "all", 
                                  numbNodes = NULL, digits = 5, ...){

  showCentr <- match.arg(showCentr, choices = c("all", "none", "degree",
                                                "betweenness", "closeness",
                                                "eigenvector"),
                         several.ok = TRUE)
  
  if("none" %in% showCentr) stopifnot(length(showCentr) == 1)
  numbNodes <- as.integer(numbNodes)

  if(is.numeric(numbNodes)){
    stopifnot(numbNodes >= 1 & numbNodes <= length(object$centralities$degree1))
  }
  
  twonets <- object$input$twoNets
  
  if(is.null(numbNodes)){
    numbNodes <- ifelse(twonets, 5L, 10L)
  }
  
  if(twonets){
    if(is.null(groupNames)){
      group1 <- paste0("group '" , object$input$groups[1], "'")
      group2 <- paste0("group '" , object$input$groups[2], "'")
      
    } else{
      group1 <- groupNames[1]
      group2 <- groupNames[2]
    }
  } else{
    group1 <- group2 <- NULL
  }


  #============================================================
  # global network properties

  if(is.na(object$globalProps$vertConnect1)){
    glob_rownames <- c("average path length",
                       "clustering coeff.",
                       "modularity",
                       "edge density")
    glob_names1 <- c("avPath1", "clustCoef1", "modularity1", "density1")
    glob_names2 <- c("avPath2", "clustCoef2", "modularity2", "density2")
  } else{
    glob_rownames <- c("average path length  ",
                       "clustering coeff.",
                       "modularity",
                       "edge density",
                       "vertex connectivity",
                       "edge connectivity")
    glob_names1 <- c("avPath1", "clustCoef1", "modularity1", "density1", 
                     "vertConnect1", "edgeConnect1")
    glob_names2 <- c("avPath2", "clustCoef2", "modularity2", "density2", 
                     "vertConnect2", "edgeConnect2")
  }
  
  if(is.na(object$globalProps$modularity1)){
    # exclude modularity
    glob_rownames <- glob_rownames[-3]
    glob_names1 <- glob_names1[-3]
    glob_names2 <- glob_names2[-3]
  }
  

  
  if(twonets){
    glob_probs <- as.data.frame(matrix(0, nrow = length(glob_rownames), ncol = 2, 
                                       dimnames = list(glob_rownames,
                                                       c(group1, group2))))
  } else{
    glob_probs <- as.data.frame(matrix(0, nrow = length(glob_rownames), ncol = 1, 
                                       dimnames = list(glob_rownames,
                                                       " ")))
  }
  
  for(i in 1:length(glob_rownames)){
    glob_probs[i, 1] <- round(as.numeric(object$globalProps[glob_names1[i]]), 
                              digits = digits)
  }
  
  if(twonets){
    for(i in 1:length(glob_rownames)){
      glob_probs[i, 2] <- round(as.numeric(object$globalProps[glob_names2[i]]), 
                                digits = digits)
    }
    
  }
  
  #============================================================
  # clustering
  
  clust <- list()
  clusttab1 <- table(object$clustering$clust1)
  clust1 <- matrix(0, nrow = 2, ncol = length(clusttab1), 
                   dimnames = list(c("name:", "freq:"),
                                   rep("", length(clusttab1))))
  clust1[1, ] <- as.numeric(names(clusttab1))
  clust1[2, ] <- clusttab1
  
  clust[[1]] <- clust1
  if(twonets){
    clusttab2 <- table(object$clustering$clust2)
    clust2 <- matrix(0, nrow = 2, ncol = length(clusttab2), 
                     dimnames = list(c("name:", "freq:"),
                                     rep("", length(clusttab2))))
    clust2[1, ] <- as.numeric(names(clusttab2))
    clust2[2, ] <- clusttab2
    clust[[2]] <- clust2
  }
  
  #============================================================
  # hubs
  
  if(twonets){
    hubs1 <- object$hubs$hubs1
    hubs2 <- object$hubs$hubs2
    
    if(length(hubs1) != length(hubs2)){
      diff <- length(hubs1) - length(hubs2)
      if(diff > 0){
        hubs2 <- c(hubs2, rep("", diff))
      } else{
        hubs1 <- c(hubs1, rep("", abs(diff)))
      }
    } 
    
    hubmat <- cbind(object$hubs$hubs1, object$hubs$hubs2)
    
    dimnames(hubmat) <- list(rep("", nrow(hubmat)),
                             c(group1, group2))
    hubs <- as.data.frame(hubmat)
  } else{
    hubmat <- as.matrix(object$hubs$hubs1)
    dimnames(hubmat) <- list(rep("", nrow(hubmat)), "")
    hubs <- hubmat
  }

  
  #============================================================
  # centrality measures
  
  top_centr <- NULL
  
  if(showCentr[1] != "none"){
    
    l <- length(object$centralities$degree1)
    
    centr_names1 <- c("degree1", "between1", "close1", "eigenv1")
    centr_names2 <- c("degree2", "between2", "close2", "eigenv2")
    centr_names <- c("degree", "betweenness", "closeness", "eigenvector")
    
    top_centr_names <- list()
    
    for(i in 1:length(centr_names1)){
      top_centr_names[[centr_names1[i]]] <- 
        names(sort(object$centralities[[centr_names1[i]]], 
                   decreasing = TRUE))[1:min(numbNodes, l)]
    }
    
    if(twonets){
      for(i in 1:length(centr_names2)){
        top_centr_names[[centr_names2[i]]] <- 
          names(sort(object$centralities[[centr_names2[i]]], 
                     decreasing = TRUE))[1:min(numbNodes, l)]
      }
    }
    
    top_centr <- list()

    for(i in 1:length(centr_names)){
      if(any(c("all", centr_names[i]) %in% showCentr)){
        
        if(twonets){
          mat1 <- 
            as.data.frame(matrix(0, nrow = numbNodes , ncol = 2,
                                 dimnames = list(top_centr_names[[centr_names1[i]]], 
                                                 c(group1, group2))))
        } else{
          mat1 <- 
            as.data.frame(matrix(0, nrow = numbNodes , ncol = 1,
                                 dimnames = list(top_centr_names[[centr_names1[i]]], 
                                                 " ")))
        }

        mat1[, 1] <- 
          round(object$centralities[[centr_names1[i]]][top_centr_names[[centr_names1[i]]]], 
                digits = digits)
        
        
        if(twonets){
          mat1[, 2] <- 
            round(object$centralities[[centr_names2[i]]][top_centr_names[[centr_names1[i]]]], 
                  digits = digits)
          
          mat2 <- 
            as.data.frame(matrix(0, nrow = numbNodes , ncol = 2,
                                 dimnames = list(top_centr_names[[centr_names2[i]]], 
                                                 c(group1, group2))))
          
          mat2[, 1] <- 
            round(object$centralities[[centr_names1[i]]][top_centr_names[[centr_names2[i]]]], 
                  digits = digits)
          mat2[, 2] <- 
            round(object$centralities[[centr_names2[i]]][top_centr_names[[centr_names2[i]]]], 
                  digits = digits)
          
          rnames <- c(rownames(mat1), "", rownames(mat2))
          
          data.frame(mat1, mat2)
          mat.tmp <- rbind(mat1, c("______","______"), mat2, make.row.names = FALSE)
          mat <- cbind(rnames, mat.tmp)
          colnames(mat) <- c(" ", colnames(mat)[2:3])
        } else{
          mat <- mat1
        }
        
        top_centr[[centr_names1[i]]] <- mat
      }
    }
  }
  


  structure(list(glob_probs = glob_probs, clust = clust, hubs = hubs,
                 central = top_centr, group1 = group1, group2 = group2, 
                 call = object$call),
            class = "summary.microNetProps")
}


#' @title Print Summary for objects of class \code{microNetProps}
#'
#' @param x object of class \code{summary.microNetProps} inheriting from a call
#'   to \code{\link{summary.microNetProps}}.
#' @param ... not used
#'
#' @method print summary.microNetProps
#'
#' @rdname summarize.microNetProps
#' @export
print.summary.microNetProps <- function(x, ...){
  cat("\nGlobal network properties:\n")
  cat("``````````````````````````\n")
  print(x$glob_probs)

  if(ncol(x$clust[[1]]) != 0){
    cat("\n\nClusters:\n")
    cat("`````````\n")
    if(length(x$clust) == 2){
      cat(x$group1, ":", sep = "")
      print(x$clust[[1]])
      cat("\n", x$group2, ":", sep = "")
      print(x$clust[[2]])
    } else{
      print(x$clust[[1]])
    }
  }

  cat("\n\nHubs:\n")
  if(ncol(x$hubs) == 2){
    cat("`````\n")
    print(x$hubs, row.names = FALSE, quote = TRUE)
  } else{
    cat("`````")
    print(x$hubs)
  }
  
  if(!is.null(x$central)){
    show_rownames <- ifelse(ncol(x$central$degree1) == 3, FALSE, TRUE)
    
    cat("\n\nCentrality measures (in decreasing order):\n")
    cat("```````````````````")
    
    if(!is.null(x$central$degree1)){
      cat("\nDegree:\n")
      print(x$central$degree1, row.names = show_rownames)
    }
    
    if(!is.null(x$central$between1)){
      cat("\nBetweenness centrality:\n")
      print(x$central$between1, row.names = show_rownames)
    }
    
    if(!is.null(x$central$close1)){
      cat("\nCloseness centrality:\n")
      print(x$central$close1, row.names = show_rownames)
    }
    
    if(!is.null(x$central$eigenv1)){
      cat("\nEigenvector centrality:\n")
      print(x$central$eigenv1, row.names = show_rownames)
    }
    
  }
}




