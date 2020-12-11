filter_samples <- function(countMat, filter, filterParam){

  countMat_orig <- countMat

  if ("highestFreq" %in% filter) {
    highestFreq <-
      ifelse(is.null(filterParam$highestFreq),
             100,
             filterParam$highestFreq)
    highestFreq <- min(nrow(countMat), highestFreq)

    if (length(filter) > 1) {
      stop('Filter method "highestFreq" not comparable with other filter methods.')
    }
    names_highFreq <- names(sort(Matrix::rowSums(countMat), 
                                 decreasing = TRUE)[1:highestFreq])
    #rmRows <- which(!rownames(countMat) %in% names_highFreq)
    keepRows <- names_highFreq

  } else{

    if ("totalReads" %in% filter) {
      totalReads <-
        ifelse(is.null(filterParam$totalReads),
               1000,
               filterParam$totalReads)
      idx_totalreads <- which(Matrix::rowSums(countMat) >= totalReads)
    } else{
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
    } else{
      idx_numbTaxa <- 1:nrow(countMat)
    }

    keepRows <- rownames(countMat[intersect(idx_totalreads, idx_numbTaxa), ])
  }

  return(keepRows)
}
