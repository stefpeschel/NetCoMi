#' @keywords internal
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

#' @keywords internal
.filterNodes <- function(adja, nodeFilter, nodeFilterPar, layout,
                         degree, between, close, eigen, cluster) {
  
  adja.alltax <- adja
  
  if (!is.null(layout) & is.matrix(layout)) {
    keep <- colnames(adja)[which(colnames(adja) %in% rownames(layout))]
    
  } else if (nodeFilter != "none") {
    
    if (nodeFilter == "highestConnect") {
      adja.tmp <- adja
      
      diag(adja.tmp) <- 0
      conct <- Matrix::rowSums(abs(adja.tmp))
      
      conct <- names(sort(conct))[1:nodeFilterPar]
      keep <- colnames(adja)[which(colnames(adja) %in% conct)]
      
    } else if (nodeFilter == "highestDegree") {
      sel <- names(sort(degree, decreasing = TRUE)[1:nodeFilterPar])
      keep <- colnames(adja)[which(colnames(adja) %in% sel)]
      
    } else if (nodeFilter == "highestBetween") {
      sel <- names(sort(between, decreasing = TRUE)[1:nodeFilterPar])
      
      keep <- colnames(adja)[which(colnames(adja) %in% sel)]
      
    } else if (nodeFilter == "highestClose") {
      sel <- names(sort(close, decreasing = TRUE)[1:nodeFilterPar])
      
      keep <- colnames(adja)[which(colnames(adja) %in% sel)]
      
    } else if (nodeFilter == "highestEigen") {
      sel <- names(sort(eigen, decreasing = TRUE)[1:nodeFilterPar])
      
      keep <- colnames(adja)[which(colnames(adja) %in% sel)]
      
    } else if (nodeFilter == "clustTaxon") {
      stopifnot(all(nodeFilterPar %in% colnames(adja)))
      
      selClust <- cluster[nodeFilterPar]
      keep <- names(cluster[cluster %in% selClust])
      #keep <- names(cluster) %in% selnodes
      
    } else if (nodeFilter == "clustMin") {
      clusttab <- table(cluster)
      selclust <- names(clusttab[clusttab >= nodeFilterPar & 
                                   names(clusttab) != 0])
      
      keep <- names(cluster[cluster %in% selclust])
      
    } else if (nodeFilter == "names") {
      stopifnot(all(nodeFilterPar %in% colnames(adja)))
      keep <- colnames(adja)[which(colnames(adja) %in% nodeFilterPar)]
    }
    
  } else {
    keep <- colnames(adja)
  }
  
  # names_alltaxa <- matrix(NA, nrow = ncol(adja.alltax1), ncol = 2)
  # names_alltaxa[, 1] <- colnames(adja.alltax1)
  # names_alltaxa[, 2] <- colnames(x$input$adja1)
  # rownames(names_alltaxa) <- names_alltaxa[,1]
  #
  # unitnames <- union(colnames(adja), colnames(adja2))
  # names_selected_taxa <- matrix(NA, nrow = length(unitnames), ncol = 2)
  # names_selected_taxa[, 1] <- unitnames
  # names_selected_taxa[, 2] <- names_alltaxa[unitnames, 2]
  # rownames(names_alltaxa) <- NULL
  #
  # if (labelsToFile == "all") {
  #   write.matrix(names_alltaxa, file="taxalabels.txt")
  # } else if (labelsToFile == "selected") {
  #   write.matrix(names_selected_taxa, file="taxalabels.txt")
  # }
  #
  #
  # colnames(names_alltaxa) <- colnames(names_selected_taxa) <- c("shortened", 
  #                                                               "original")
  #
  # taxalabels <- list(all_taxa = names_alltaxa, 
  #                    selected_taxa = names_selected_taxa)
  
  rm(adja.alltax)
  
  return(keep = keep)
  
}


#' @keywords internal
.filterSamples <- function(countMat, filter, filterParam) {
  
  countMat_orig <- countMat
  
  if ("highestFreq" %in% filter) {
    highestFreq <-
      ifelse(is.null(filterParam$highestFreq),
             100,
             filterParam$highestFreq)
    highestFreq <- min(nrow(countMat), highestFreq)
    
    if (length(filter) > 1) {
      stop('Filter method "highestFreq" not comparable with other ', 
           'filter methods.')
    }
    names_highFreq <- names(sort(Matrix::rowSums(countMat), 
                                 decreasing = TRUE)[1:highestFreq])
    #rmRows <- which(!rownames(countMat) %in% names_highFreq)
    keepRows <- names_highFreq
    
  } else {
    
    if ("totalReads" %in% filter) {
      totalReads <-
        ifelse(is.null(filterParam$totalReads),
               1000,
               filterParam$totalReads)
      idx_totalreads <- which(Matrix::rowSums(countMat) >= totalReads)
    } else {
      idx_totalreads <- 1:nrow(countMat)
    }
    
    if ("numbTaxa" %in% filter) {
      numbTaxa <- ifelse(is.null(filterParam$numbTaxa), 
                         0.1, 
                         filterParam$numbTaxa)
      stopifnot(numbTaxa > 0)
      if (numbTaxa < 1) {
        numbTaxa <- round(numbTaxa * nrow(countMat))
      }
      
      idx_numbTaxa = numeric(0)
      
      for (i in 1:nrow(countMat)) {
        if (length(which(!countMat[i,] == 0)) >= numbTaxa) {
          idx_numbTaxa <- append(idx_numbTaxa, i)
        }
      }
    } else {
      idx_numbTaxa <- 1:nrow(countMat)
    }
    
    keepRows <- rownames(countMat[intersect(idx_totalreads, idx_numbTaxa), ])
  }
  
  return(keepRows)
}


#' @keywords internal
.filterTaxa <- function(countMat, filter, filterParam) {
  
  countMat_orig <- countMat
  #---------------------------------------------------------------------------
  # filter taxa of interest
  
  if ("highestVar" %in% filter) {
    highestVar <-
      ifelse(is.null(filterParam$highestVar),
             100,
             filterParam$highestVar)
    
    highestVar <- min(ncol(countMat), highestVar)
    
    if (length(filter) > 1) {
      stop('Filter method "highestVar" not comparable with other ', 
           'filter methods.')
    }
    var_taxa <- apply(countMat, 2, var)
    keepCols <- names(sort(var_taxa, decreasing = TRUE)[1:highestVar])
    
    #countMat <- countMat[, names_highvar]
    
  } else if ("highestFreq" %in% filter) {
    highestFreq <-
      ifelse(is.null(filterParam$highestFreq),
             100,
             filterParam$highestFreq)
    
    highestFreq <- min(ncol(countMat), highestFreq)
    
    if (length(filter) > 1) {
      stop('Filter method "highestFreq" not comparable with other ', 
           'filter methods.')
    }
    keepCols <- names(sort(Matrix::colSums(countMat), 
                           decreasing = TRUE)[1:highestFreq])
    #countMat <- countMat[, names_highFreq]
    
  } else {
    if ("relFreq" %in% filter) {
      relFreq <-
        ifelse(is.null(filterParam$relFreq),
               0.01,
               filterParam$relFreq)
      idx_relfreq <- which(Matrix::colSums(countMat) >= 
                             sum(Matrix::colSums(countMat)) * relFreq)
    } else {
      idx_relfreq <- 1:ncol(countMat)
    }
    
    
    if ("totalReads" %in% filter) {
      totalReads <-
        ifelse(is.null(filterParam$totalReads),
               1000,
               filterParam$totalReads)
      idx_totalreads <- which(Matrix::colSums(countMat) >= totalReads)
    } else {
      idx_totalreads <- 1:ncol(countMat)
    }
    
    if ("numbSamp" %in% filter) {
      numbSamp <- ifelse(is.null(filterParam$numbSamp), 0.1, 
                         filterParam$numbSamp)
      
      stopifnot(numbSamp > 0)
      
      if (numbSamp < 1) {
        numbSamp <- round(numbSamp * nrow(countMat))
      }
      
      idx_numbSamp = numeric(0)
      
      for (i in 1:ncol(countMat)) {
        if (length(which(!countMat[, i] == 0)) >= numbSamp) {
          idx_numbSamp <- append(idx_numbSamp, i)
        }
      }
    } else {
      idx_numbSamp <- 1:ncol(countMat)
    }
    
    keepCols <- colnames(countMat[, intersect(intersect(idx_relfreq,
                                                        idx_totalreads), 
                                              idx_numbSamp)])
  }
  
  #rmCols <- which(!colnames(countMat_orig) %in% colnames(countMat))
  
  return(keepCols)
}
