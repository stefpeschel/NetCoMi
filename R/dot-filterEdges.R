.filterEdges <- function(adja, edgeFilter, edgeFilterPar) {
  
  if (edgeFilter == "threshold") {
    
    adja[abs(adja) < edgeFilterPar] <- 0
    
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

