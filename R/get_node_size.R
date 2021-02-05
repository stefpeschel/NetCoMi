get_node_size <- function(nodeSize, normPar, nodeSizeSpread, adja, countMat, 
                          normCounts, assoType, kept, cexNodes, cexHubs, hubs, 
                          highlightHubs, degree, between,  close, eigen){
  
  getsize <- function(x, nodeSizeSpread, cexNodes){
    (x - min(x)) / (max(x) - min(x)) * nodeSizeSpread + cexNodes
  }
  
  if(nodeSize == "fix"){
    nodeSize <- (7*exp(-ncol(adja)/80)+1) * cexNodes
    
    if(cexNodes != cexHubs & highlightHubs){
      nodeSize.tmp <- nodeSize
      hubSize <- (7*exp(-ncol(adja)/80)+1) * cexHubs
      
      hubs <- hubs[hubs %in% colnames(adja)]
      
      nodeSize <- rep(nodeSize.tmp, ncol(adja))
      nodeSize[match(hubs, rownames(adja))] <- hubSize
    }
    
  } else if(nodeSize == "degree"){
    degree <- abs(degree[kept])
    nodeSize <- getsize(degree, nodeSizeSpread, cexNodes)
    
  } else if(nodeSize == "betweenness"){
    between <- abs(between[kept])
    nodeSize <- getsize(between, nodeSizeSpread, cexNodes)
    
  } else if(nodeSize == "closeness"){
    close <- abs(close[kept])
    nodeSize <- getsize(close, nodeSizeSpread, cexNodes)
    
  } else if(nodeSize == "eigenvector"){
    
    eigen <- abs(eigen[kept])
    nodeSize <- getsize(eigen, nodeSizeSpread, cexNodes)
    
  } else if(nodeSize == "counts"){
    
    if(is.null(countMat)){
      stop("Count matrix must be available for node sizes based on counts.")
    }
    
    if(assoType != "dissimilarity"){
      absfreq <- apply(countMat, 2, stats::median)[kept]
    } else{
      absfreq <- apply(countMat, 1, stats::median)[kept]
    }
    
    names(absfreq) <- colnames(adja)
    absfreq <- absfreq[colnames(adja)]
    
    nodeSize <- absfreq/max(absfreq) * nodeSizeSpread * cexNodes + 1
    
  } else if(nodeSize == "normCounts"){
    
    if(assoType != "dissimilarity"){
      normCounts_sum <- colSums(normCounts)[kept]
      
    } else{ # dissimilarity network
      normCounts_sum <- rowSums(normCounts)[kept]
      warning("Node sizes based on normalized counts not meaningful for dissimilarity networks.")
    }
    
    nodeSize <- getsize(normCounts_sum, nodeSizeSpread, cexNodes)
  } else if(nodeSize %in% c("TSS", "fractions", "CSS", "COM", "rarefy", "VST", 
                            "clr", "mclr")){
    if(is.null(countMat)){
      stop("Count matrix must be available for node sizes based on normalized counts.")
    }
    
    if(nodeSize == "clr"){
      if(any(countMat == 0)) countMat <- countMat + 1
    }
    
    normCounts <- norm_counts(countMat, normMethod = nodeSize, 
                              normParam = normPar, zeroMethod = "none", 
                              needfrac = FALSE, verbose = FALSE)
    
    if(assoType != "dissimilarity"){ 
      normCounts_sum <- colSums(normCounts)[kept]
      
    } else{ # dissimilarity network
      normCounts_sum <- rowSums(normCounts)[kept]
      warning("Node sizes based on normalized counts not meaningful for dissimilarity networks.")
    }
    
    nodeSize <- getsize(normCounts_sum, nodeSizeSpread, cexNodes)
    
  }
  
  return(nodeSize)
  
}
