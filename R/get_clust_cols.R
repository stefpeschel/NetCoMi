get_clust_cols <- function(clust1, clust2, adja1, adja2, kept1, kept2, isempty1,
                           isempty2, colorVec, sameClustCol, sameColThresh,
                           twoNets){

  if(twoNets){
    clust1 <- clust1[kept1]
    clust2 <- clust2[kept2]

    #-----------------------------------------
    # Ensure that cluster names are numbers
    if(!isempty1){
      noclust1 <- which(clust1 %in% c(0, NA))
      clusttab <- table(clust1)
      cnames <- names(clusttab)
      cnames <- cnames[cnames != "0"]
      for(i in 1:length(cnames)){
        clust1[clust1 == cnames[i]] <- i
      }
    }

    if(!isempty2){
      noclust2 <- which(clust2 %in% c(0, NA))

      clusttab <- table(clust2)
      cnames <- names(clusttab)
      cnames <- cnames[cnames != "0"]
      for(i in 1:length(cnames)){
        clust2[clust2 == cnames[i]] <- i
      }
    }
    #-----------------------------------------

    if(!(isempty1 | isempty2)){

      clustcols <- grDevices::rainbow(max(clust1) + max(clust2))

      if(!is.null(colorVec)){
        clustcols <- clustcols[!clustcols %in% colorVec]
        clustcols <- c(colorVec, clustcols)
      }

      if(sameClustCol){
        nc1 <- max(clust1)
        nc2 <- max(clust2)
        clustcol1 <- clust1
        clustcol2 <- clust2

        cl1 <- 1:nc1
        cl2 <- 1:nc2

        matchmat <- matrix(0, nrow = nc1, ncol = nc2)

        for(i in 1:nc1){
          matchmat[i, ] <- sapply(1:nc2, function(x){
            length( intersect(names(clust1[clust1 == i]),
                              names(clust2[clust2 == x])))
          })
        }

        maxgreater <- max(matchmat) >= sameColThresh

        i = 1
        while(maxgreater){
          ismax <- which(matchmat == max(matchmat), arr.ind = TRUE)[1,]
          clustcol1[clust1 == ismax[1]] <- clustcols[i]
          clustcol2[clust2 == ismax[2]] <- clustcols[i]
          matchmat[ismax[1], ] <- 0
          matchmat[, ismax[2]] <- 0
          cl1 <- cl1[cl1 != ismax[1]]
          cl2 <- cl2[cl2 != ismax[2]]

          maxgreater <- max(matchmat) >= sameColThresh
          i = i+1
        }

        for(c in cl1){
          clustcol1[clust1 == c] <- clustcols[i]
          i = i+1
        }
        for(c in cl2){
          clustcol2[clust2 == c] <- clustcols[i]
          i = i+1
        }
        
        if(!is.null(colorVec) && (i-1) > length(colorVec)){
          
          warning(i-1, " colors needed but 'colorVec' has only length ", 
                  length(colorVec), ". Missing colors are filled up with colors from rainbow().")
          
        }

        clustcol1[noclust1] <- "grey80"
        clustcol2[noclust2] <- "grey80"

      } else{

        clustcols <- c("grey80", clustcols)
        
        nc1 <- max(clust1)
        clust1 <- clust1 + 1
        
        clust2 <- clust2 + nc1 + 1
        clust2[noclust2] <- 1

        clustcol1 <- clustcols[clust1]
        clustcol2 <- clustcols[clust2]
      }

    } else{

      if(isempty1){
        clustcol1 <- rep("grey80", ncol(adja1))

      } else{
        clust1 <- clust1
        nc1 <- max(clust1)

        clustcols <- grDevices::rainbow(nc1)

        if(length(noclust1) != 0){
          clustcols <- c("grey80", clustcols)
          clust1 <- clust1 + 1
        }
        clustcol1 <- clustcols[clust1]
        names(clustcol1) <- names(clust1)
      }

      if(isempty2){
        clustcol2 <- rep("grey80", ncol(adja2))
      } else{
        clust2 <- clust2
        nc2 <- max(clust2)
        
        clustcols <- grDevices::rainbow(nc2)

        if(length(noclust2) != 0){
          clustcols <- c("grey80", clustcols)
          clust2 <- clust2 + 1
        }
        clustcol2 <- clustcols[clust2]
        names(clustcol2) <- names(clust2)
      }
    }
    names(clustcol1) <- names(clust1)
    names(clustcol2) <- names(clust2)

  } else{
    clust1 <- clust1[kept1]

    if(!is.null(clust1)){
      noclust <- which(clust1 == 0)

      clusttab <- table(clust1)
      cnames <- names(clusttab)
      cnames <- cnames[cnames != "0"]
      for(i in 1:length(cnames)){
        clust1[clust1 == cnames[i]] <- i
      }

      clustcols <- grDevices::rainbow(max(clust1))

      if(!is.null(colorVec)){
        if(length(colorVec) < max(clust1)){
          warning("'colorVec' too short (number of clusters = ", max(clust1), ").")
        }
        clustcols <- clustcols[!clustcols %in% colorVec]
        clustcols <- c(colorVec, clustcols)
      }
      if(length(noclust) != 0){
        clustcols <- c("grey80", clustcols)
        clust1 <- clust1 + 1
      }
      clustcol1 <- clustcols[clust1]
      names(clustcol1) <- names(clust1)
    } else{
      warning('No clusterings returned from "netProperties".')
      clust1 <- NULL
      clustcol1 <- rep("grey40", ncol(adja1))
    }

    clustcol2 <- NULL

  }


    return(list(clustcol1 = clustcol1, clustcol2 = clustcol2))
}
