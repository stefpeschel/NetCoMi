rename_taxa <- function(x, toRename = c("both", "cols", "rows"),
                       shortenLabels = c("intelligent", "simple", "none"),
                       labelLength = 6, labelPattern = NULL,
                       charToRm = NULL){

  toRename <- match.arg(toRename)
  shortenLabels <- match.arg(shortenLabels)

  if(is.null(labelPattern)){
    labelPattern <- c(4, "'", 3)
  } else{
    stopifnot(length(labelPattern) == 3)
  }

  if(toRename %in% c("both", "cols")){
    labels <- colnames(x)
  } else{
    labels <- rownames(x)
  }


  if(shortenLabels == "intelligent"){

    for(char in charToRm){
      labels <- gsub(char, "", labels)
    }

    labels[labels == ""] <- "-"
    shortlabels <- substring(labels, 1,labelLength)

    dupli <- which(duplicated(shortlabels))

    # find duplicates in variable names
    while(length(dupli) > 0){
      ind <- which(shortlabels == shortlabels[dupli[1]])
      dupnames <- strsplit(labels[ind], "")

      lvec <- unlist(lapply(dupnames, length))
      l <- max(lvec)
      if(min(lvec) != max(lvec)){
        for(i in 1:length(dupnames)){
          if(lvec[i] < l){
            dupnames[[i]] <- c(dupnames[[i]], rep(" ", l-lvec[i]))
          }
        }
      }

      pos <- first_unequal_element(dupnames[[1]], dupnames[[2]])
      cut1 <- as.numeric(labelPattern[1])
      cut2 <- as.numeric(labelPattern[3])
      cut2 <- min(cut2, l-pos+1)


      for(k in 1:length(dupnames)){
        if(dupnames[[k]][pos] == "["){
          shortlabels[ind[k]] <- paste0(paste(dupnames[[k]][1:cut1], collapse = ""),
                                        paste(dupnames[[k]][pos:(pos+cut2+2)], collapse = ""))
        } else{
          str1 <- paste(dupnames[[k]][1:cut1], collapse = "")
          str2 <- labelPattern[2]
          str3 <- paste(dupnames[[k]][pos:(pos+cut2-1)], collapse = "")
          if(str3 == " "){
            shortlabels[ind[k]] <- str1
          } else{
            shortlabels[ind[k]] <- paste0(str1, str2, str3)
          }

        }
      }
      dupli <- dupli[!dupli %in% ind]
    }

    for(i in 1:length(shortlabels)){
      splitname <- strsplit(shortlabels[i], "")
      if("[" %in% splitname[[1]]) shortlabels[i] <- paste0(shortlabels[i], "]")
    }

    labels <- shortlabels

  } else if(shortenLabels == "simple"){

    for(char in charToRm){
      labels <- gsub(char, "", labels)
    }

    labels[labels == ""] <- "-"
    labels <- substr(labels, 1, labelLength)
    
  } else if(!is.null(charToRm)){

    for(char in charToRm){
      labels <- gsub(char, "", labels)
    }

  }

  if(toRename == "both"){
    dimnames(x) <- list(labels, labels)
  } else if(toRename == "cols"){
    colnames(x) <- labels
  } else{
    rownames(x) <- labels
  }

  return(x)
}


