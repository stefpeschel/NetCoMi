norm_counts <- function(countMat, normMethod, normParam, zeroMethod, needfrac,
                        verbose){


  if(normMethod == "none"){
    countMat_norm <- countMat
    if(needfrac){
      attributes(countMat_norm)$scale <- "fractions"
      if(verbose %in% 2:3) message("Counts normalized by total sum scaling.")
    } else{
      attributes(countMat_norm)$scale <- "counts"
    }

    return(countMat_norm)
  }


  if(normMethod %in% c("rarefy", "VST") & attributes(countMat)$scale != "integer"){
    countMat.tmp <- countMat
    countMat <- apply(countMat.tmp, 2, as.integer)
    rownames(countMat) <- rownames(countMat.tmp)
  }

  #-----------------------------------------------------------------------------
  if(normMethod %in% c("fractions", "TSS")){

    if(attributes(countMat)$scale != "fractions"){
      countMat_norm <- t(apply(countMat, 1, function(x) x/sum(x)))
      attributes(countMat_norm)$scale <- "fractions"
    }else{
      countMat_norm <- countMat
    }

    if(verbose %in% 2:3) message("Counts normalized by total sum scaling.")

  } else if(normMethod == "CSS"){

    normParam$obj <- t(countMat)
    if(is.null(normParam$p)) normParam$p <- 0.75
    if(is.null(normParam$sl)) normParam$sl <- 100

    if(verbose %in% 2:3){
      message("Execute cumNormMat() for cumulative sum scaling ... ",
              appendLF = FALSE)
    }
    countMat_norm <- t(do.call("cumNormMat", normParam))
    if(verbose %in% 2:3) message("Done.")

    attributes(countMat_norm)$scale <- "CSS normalized"

  } else if(normMethod == "COM"){

    countMat_norm <- t(apply(countMat, 1, function(x){
      (x * min(rowSums(countMat)) / sum(x))
    }))

    if(verbose %in% 2:3) message("Counts normalized by common sum scaling.")

    attributes(countMat_norm)$scale <- "COM normalized"

  } else if(normMethod == "rarefy"){

    normParam$x <- countMat
    if(is.null(normParam$sample)) normParam$sample <- min(rowSums(countMat))

    if(verbose %in% 2:3){
      message("Execute rrarefy() for rarefaction ... ", appendLF = FALSE)
    }
    countMat_norm <- do.call("rrarefy", normParam)
    if(verbose %in% 2:3) message("Done.")

    attributes(countMat_norm)$scale <- "rarefied"

  } else if(normMethod == "VST"){
    
    normParam$object <- t(countMat)
    
    if(verbose %in% 2:3){
      message("Execute varianceStabilizingTransformation() for VST normalization ... ",
              appendLF = FALSE)
    }
    if(zeroMethod == "none"){
      countMat_norm_t <- try(do.call("varianceStabilizingTransformation",
                                     normParam), silent = TRUE)
      if(class(countMat_norm) == "try-error"){
        stop("Every variable contains at least one zero. ",
             "VST normalization not possible without zero replacement.")
      }
      
    } else{
      countMat_norm_t <- do.call("varianceStabilizingTransformation", normParam)
    }
    
    countMat_norm <- t(countMat_norm_t)
    rm(countMat_norm_t)
    
    if(verbose %in% 2:3) message("Done.")
    
    attributes(countMat_norm)$scale <- "VST normalized"
    
  } else if(normMethod == "clr"){

    if(attributes(countMat)$scale != "fractions"){
      countMat <- t(apply(countMat, 1, function(x) x/sum(x)))
    }

    if(verbose %in% 2:3){
      message("Counts transformed to fractions.")
    }

    normParam$x <- countMat
    if(verbose %in% 2:3){
      message("Execute cenLR() for clr transformation ... ",
              appendLF = FALSE)
    }

    if(is.null(normParam$base)) normParam$base <- exp(1)
    countMat_norm <- do.call("cenLR", normParam)$x.clr

    if(verbose %in% 2:3) message("Done.")

    attributes(countMat_norm)$scale <- "clr transformed"
  } else{
    warning("No normalization conducted. ",
    "'normMethod' must be one of 'none', 'pseudo', 'multRepl', 'alrEM', 'bayesMult'.")
  }

  return(countMat_norm)
}
