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
      message("Execute cumNormMat(){metagenomeSeq} ... ",
              appendLF = FALSE)
    }
    countMat_norm <- t(do.call(metagenomeSeq::cumNormMat, normParam))
    if(verbose %in% 2:3) message("Done.")

    attributes(countMat_norm)$scale <- "CSS normalized"

  } else if(normMethod == "COM"){
    
    if(verbose %in% 2:3){
      message("Normalize counts by common sum scaling ... ",
              appendLF = FALSE)
    }

    countMat_norm <- t(apply(countMat, 1, function(x){
      (x * min(rowSums(countMat)) / sum(x))
    }))

    if(verbose %in% 2:3) message("Done.")

    attributes(countMat_norm)$scale <- "COM normalized"

  } else if(normMethod == "rarefy"){

    normParam$x <- countMat
    if(is.null(normParam$sample)){
      normParam$sample <- min(Matrix::rowSums(countMat))
    }

    if(verbose %in% 2:3){
      message("Execute rrarefy(){vegan} ... ", appendLF = FALSE)
    }
    
    countMat_norm <- do.call(vegan::rrarefy, normParam)
    
    if(verbose %in% 2:3) message("Done.")

    attributes(countMat_norm)$scale <- "rarefied"

  } else if(normMethod == "VST"){
    
    normParam$object <- t(countMat)
    
    if(verbose %in% 2:3){
      message("Execute varianceStabilizingTransformation(){DESeq2} ... ",
              appendLF = FALSE)
    }
    if(zeroMethod == "none"){
      countMat_norm_t <- try(do.call(DESeq2::varianceStabilizingTransformation,
                                     normParam), silent = TRUE)
      
      if(class(countMat_norm_t) == "try-error"){
        stop("Every variable contains at least one zero. ",
             "VST normalization not possible without zero replacement.")
      }
      
    } else{
      countMat_norm_t <- do.call(DESeq2::varianceStabilizingTransformation, 
                                 normParam)
    }
    
    countMat_norm <- t(countMat_norm_t)
    rm(countMat_norm_t)
    
    if(verbose %in% 2:3) message("Done.")
    
    attributes(countMat_norm)$scale <- "VST normalized"
    
  } else if(normMethod == "clr"){

    normParam$x.f <- countMat
    normParam$mar <- 1
    if(is.null(normParam$base)) normParam$base <- exp(1)
    
    if(verbose %in% 2:3){
      message("Execute clr(){SpiecEasi} ... ", appendLF = FALSE)
    }
    
    countMat_norm <- t(do.call(SpiecEasi::clr, normParam))

    if(verbose %in% 2:3) message("Done.")

    attributes(countMat_norm)$scale <- "clr transformed"
    
  } else if(normMethod == "mclr"){
    
    normParam$dat <- countMat
    
    if(verbose %in% 2:3){
      message("Execute mclr(){SPRING} ... ", appendLF = FALSE)
    }
    
    countMat_norm <- do.call(SPRING::mclr, normParam)
    
    if(verbose %in% 2:3) message("Done.")
    
    attributes(countMat_norm)$scale <- "mclr transformed"
    
  } else{
    warning("No normalization conducted. ",
    "'normMethod' must be one of 'none', 'fractions', 'TSS', ", 
    "'CSS', 'COM', 'rarefy', 'VST', 'clr', 'mclr'.")
  }

  return(countMat_norm)
}
