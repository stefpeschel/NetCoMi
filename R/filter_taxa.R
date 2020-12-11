filter_taxa <- function(countMat, filter, filterParam){

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
      stop('Filter method "highestVar" not comparable with other filter methods.')
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
      stop('Filter method "highestFreq" not comparable with other filter methods.')
    }
    keepCols <- names(sort(Matrix::colSums(countMat), 
                           decreasing = TRUE)[1:highestFreq])
    #countMat <- countMat[, names_highFreq]

  } else{
    if ("relFreq" %in% filter) {
      relFreq <-
        ifelse(is.null(filterParam$relFreq),
               0.01,
               filterParam$relFreq)
      idx_relfreq <-
        which(Matrix::colSums(countMat) >= sum(Matrix::colSums(countMat)) * relFreq)
    } else{
      idx_relfreq <- 1:ncol(countMat)
    }


    if ("totalReads" %in% filter) {
      totalReads <-
        ifelse(is.null(filterParam$totalReads),
               1000,
               filterParam$totalReads)
      idx_totalreads <- which(Matrix::colSums(countMat) >= totalReads)
    } else{
      idx_totalreads <- 1:ncol(countMat)
    }

    if ("numbSamp" %in% filter) {
      numbSamp <- ifelse(is.null(filterParam$numbSamp), 0.1, filterParam$numbSamp)
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
    } else{
      idx_numbSamp <- 1:ncol(countMat)
    }

    keepCols <- colnames(countMat[, intersect(intersect(idx_relfreq,
                                                        idx_totalreads), idx_numbSamp)])
  }

  #rmCols <- which(!colnames(countMat_orig) %in% colnames(countMat))

  return(keepCols)
}
