filter_nodes <- function(adja, nodeFilter, nodeFilterPar, layout,
                          degree, between, close, eigen, cluster){

  adja.alltax <- adja

  if(!is.null(layout) & is.matrix(layout)){
    keep <- colnames(adja)[which(colnames(adja) %in% rownames(layout))]

  } else if(nodeFilter != "none"){

    if(nodeFilter == "highestConnect"){
      adja.tmp <- adja

      diag(adja.tmp) <- 0
      conct <- Matrix::rowSums(abs(adja.tmp))

      conct <- names(sort(conct))[1:nodeFilterPar]
      keep <- colnames(adja)[which(colnames(adja) %in% conct)]

    } else if(nodeFilter == "highestDegree"){
      sel <- names(sort(degree, decreasing = TRUE)[1:nodeFilterPar])
      keep <- colnames(adja)[which(colnames(adja) %in% sel)]

    } else if(nodeFilter == "highestBetween"){
      sel <- names(sort(between, decreasing = TRUE)[1:nodeFilterPar])

      keep <- colnames(adja)[which(colnames(adja) %in% sel)]

    } else if(nodeFilter == "highestClose"){
      sel <- names(sort(close, decreasing = TRUE)[1:nodeFilterPar])

      keep <- colnames(adja)[which(colnames(adja) %in% sel)]

    } else if(nodeFilter == "highestEigen"){
      sel <- names(sort(eigen, decreasing = TRUE)[1:nodeFilterPar])

      keep <- colnames(adja)[which(colnames(adja) %in% sel)]

    } else if(nodeFilter == "clustTaxon"){
      stopifnot(all(nodeFilterPar %in% colnames(adja)))

      selClust <- cluster[nodeFilterPar]
      keep <- names(cluster[cluster %in% selClust])
      #keep <- names(cluster) %in% selnodes

    } else if(nodeFilter == "clustMin"){
      clusttab <- table(cluster)
      selclust <- names(clusttab[clusttab >= nodeFilterPar & names(clusttab) != 0])

      keep <- names(cluster[cluster %in% selclust])

    } else if(nodeFilter == "names"){
      stopifnot(all(nodeFilterPar %in% colnames(adja)))
      keep <- colnames(adja)[which(colnames(adja) %in% nodeFilterPar)]
    }

  } else{
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
  # if(labelsToFile == "all"){
  #   write.matrix(names_alltaxa, file="taxalabels.txt")
  # } else if(labelsToFile == "selected"){
  #   write.matrix(names_selected_taxa, file="taxalabels.txt")
  # }
  #
  #
  # colnames(names_alltaxa) <- colnames(names_selected_taxa) <- c("shortened", "original")
  #
  # taxalabels <- list(all_taxa = names_alltaxa, selected_taxa = names_selected_taxa)

  rm(adja.alltax)
  
  return(keep = keep)

}

