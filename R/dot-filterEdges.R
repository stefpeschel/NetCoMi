.filterEdges <- function(adja, assoEst, dissEst, edgeFilter, edgeFilterPar) {
  if (edgeFilter == "threshold") {
    
    if (is.null(dissEst)) { # if association network
      assoEst <- assoEst[rownames(adja), colnames(adja)]
      
      adja[abs(assoEst) < edgeFilterPar] <- 0
      
    } else { # if dissimilarity network
      dissEst <- dissEst[rownames(adja), colnames(adja)]
      
      adja[dissEst > edgeFilterPar] <- 0
    }

  } else if (edgeFilter == "highestWeight") {
    
    adjaSort <- sort(adja[lower.tri(adja)], decreasing = TRUE)
    
    if (edgeFilterPar > length(adjaSort)) {
      stop(paste0("'edgeFilterPar' must be smaller than maximum number ", 
                  "of edges (", length(adjaSort), ")."))
    }
    
    adja[adja < adjaSort[edgeFilterPar]] <- 0
  }
  
  return(adja)
}

