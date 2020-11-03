get_node_size <- function(nodeSize, nodeSizeSpread, adja, countMat, normCounts,
                        kept, cexNodes, cexHubs, hubs, highlightHubs,
                        degree, between,  close, eigen){
  
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

  } else if(nodeSize %in% c("counts")){
    if(ncol(adja) == ncol(countMat)){
      absfreq <- apply(countMat, 2, stats::median)[kept]
    } else{
      absfreq <- apply(countMat, 1, stats::median)[kept]
    }

    names(absfreq) <- colnames(adja)
    absfreq <- absfreq[colnames(adja)]

    nodeSize <- absfreq/max(absfreq) * nodeSizeSpread * cexNodes + 1

  } else if(nodeSize %in% c("normCounts")){

    if(ncol(adja) == ncol(countMat)){
      normfreq <- apply(normCounts, 2, stats::median) + 0.01
    } else{
      normfreq <- apply(normCounts, 1, stats::median) + 0.01
    }

    normfreq <- normfreq[kept]
    names(normfreq) <- colnames(adja)

    greater <- FALSE
    while(!greater){
      normfreq <- normfreq * 10
      greater <- all(normfreq > 1)
    }

    normfreq <- normfreq[colnames(adja)]
    nodeSize <- normfreq/max(normfreq) * nodeSizeSpread * cexNodes + 1
  }

  return(nodeSize)

}
