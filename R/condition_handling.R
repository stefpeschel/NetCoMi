condition_handling <- function(dataType, assoType, data2, measure, normMethod,
         zeroMethod, sparsMethod, dissFunc, sampleSize, verbose){


  # set to TRUE, if a measure or normMethod needs fractions or integer as input:
  needfrac <- FALSE
  needint <- FALSE

  #-----------------------------------------------------------------------------
  # warnings regarding compositionality
  if(measure %in% c("kld", "jsd", "jeffrey")){
    if(verbose > 0){
      message("Attention! The chosen measure is not robust to compositional effects.\n")
    }
  }

  if(measure %in% c("pearson", "spearman", "bicor")){
    if(!normMethod %in% c("VST", "clr", "mclr")){
      if(verbose > 0){
        message("Attention! The chosen combination of association measure
and normalization is not robust to compositional effects.\n")
      }
    }
  }

  if(measure %in% c("bray", "euclidean")){
    if(!normMethod %in% c("VST", "clr", "mclr")){
      if(verbose > 0){
        message("Attention! The chosen combination of dissimilarity measure
and normalization is not robust to compositional effects.\n")
      }
    }
  }

  #-----------------------------------------------------------------------------
  # exception handling - normalization and zero replacement

  msg <- character(0)

  if (dataType == "counts") {

    if (measure %in% c("sparcc", "spring", "spieceasi")){
      if (zeroMethod != "none") {
        zeroMethod <- "none"
        msg <- c(msg, paste0("Zero handling included in '", measure, "'."))
      }

      if (normMethod != "none") {
        normMethod <- "none"
        msg <- c(msg, paste0("Normalization ignored for measure '",
                             measure, "'."))
      }

      #-------------------------------------------------------------------------
    } else if (measure %in% c("cclasso", "gcoda")){
      if (zeroMethod == "none"){
        zeroMethod <- "multRepl"
        msg <- c(msg, paste0("Zero replacement needed for measure '",
                             measure, "'. 'multRepl' used."))
      }

      if(!normMethod %in% c("none", "fractions")){
        msg <- c(msg, paste0("Normalization ignored for measure '",
                             measure, "'."))
      }

      if(zeroMethod == "pseudo"){
        normMethod <- "fractions"
      } else{
        normMethod <- "none"
        needfrac <- TRUE
      }

      #-------------------------------------------------------------------------
    } else if (measure %in% c("aitchison", "ckld")){
      if (zeroMethod == "none"){
        zeroMethod <- "multRepl"
        msg <- c(msg, paste0("Zero replacement needed for measure '",
                             measure, "'. 'multRepl' used."))
      }

      if(!normMethod %in% c("fractions")){
        msg <- c(msg, paste0("Counts normalized to fractions for measure '",
                             measure, "'."))
      }

      if(zeroMethod == "pseudo"){
        normMethod <- "fractions"
      } else{
        normMethod <- "none"
        needfrac <- TRUE
      }

      #-------------------------------------------------------------------------
    } else if(measure == "ccrepe") {
      if (normMethod != "fractions"){
        msg <- c(msg, paste0("Measure '", measure, "' needs fractions as input. ",
                             "'normMethod' changed to 'fractions'."))
      }
      if(zeroMethod %in% c("none", "pseudo")){
        normMethod <- "fractions"
      } else{
        normMethod <- "none"
        needfrac <- TRUE
      }

      #-------------------------------------------------------------------------
    } else if(measure %in% c("propr")){

      if (zeroMethod == "none"){
        zeroMethod <- "multRepl"
        msg <- c(msg, paste0("Zero replacement needed for measure '",
                             measure, "'. '", zeroMethod, "' used."))
      }

      if (normMethod != "none") {
        normMethod <- "none"
        msg <- c(msg, paste0("Normalization ignored for measure '",
                             measure, "'."))
      }

      #-------------------------------------------------------------------------
    } else if(measure %in% c("kld", "jeffrey", "jsd")){

      if (zeroMethod == "none") {
        if(normMethod %in% c("VST", "rarefy")){
          zeroMethod <- "pseudo"
        } else{
          zeroMethod <- "multRepl"
        }
        msg <- c(msg, paste0("Zero replacement needed for measure '",
                             measure, "'. '", zeroMethod, "' used."))
      }

    }

    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    if(normMethod %in% c("TSS", "fractions") & !zeroMethod %in% c("none", "pseudo")){
      normMethod <- "none"
      needfrac <- TRUE
    }

    if(normMethod %in% c("VST", "rarefy") & !zeroMethod %in% c("none")){
      needint <- TRUE
    }

    if(normMethod == "clr"){
      if(zeroMethod == "none"){
        zeroMethod <- "multRepl"
        msg <- c(msg, paste0("Zero replacement needed for clr transformation. '",
                             zeroMethod, "' used."))
      }

      needfrac <- TRUE
    }

    if (measure %in% c("spring", "spieceasi", "gcoda")) {
      if (sparsMethod != "none") {
        sparsMethod <- "none"
        msg <- c(msg, paste0("Sparsification included in '", measure, "'."))
      }
      
    }

  }

  #-----------------------------------------------------------------------------
  # exception handling: sparsification and transformation

  if (assoType == "dissimilarity") {
    if (!sparsMethod %in% c("none", "threshold", "knn")) {
      stop(
        "Sparsification method '", sparsMethod,
        "' not implemented for dissimilarities.
        Possible options are: 'none', 'threshold', and 'knn'."
      )
    }

  }

  if (assoType == "proportionality") {
    if (!sparsMethod %in% c("none", "threshold")) {
      stop(
        "Edge selection method '", sparsMethod,
        "' not implemented for proportionality.
        Possible options are: 'none' and 'threshold'."
      )
    }
  }

  if (dataType != "counts") {
    if (sparsMethod == "t-test") {
      if (is.null(sampleSize)) {
        stop("Sample size necessary for Student's t-test.")
      }

      stopifnot(is.vector(sampleSize))

      if (is.null(data2)) {
        if (length(sampleSize) > 1) {
          warning("Length of 'sampleSize' > 1. Only the first element is taken.")
          sampleSize <- sampleSize[1]
        }
      } else{
        stopifnot(length(sampleSize) == 2)
      }
    }

    if (sparsMethod == "bootstrap" & dataType != "counts") {
      stop("Count matrix needed for bootstrapping.")
    }

  }

  if(sparsMethod == "softThreshold"){
    if(assoType != "correlation"){
      stop("Sparsification method 'softThreshold' only possible for correlations.")
    }

    if(dissFunc != "TOMdiss"){
      dissFunc <- "TOMdiss"
      msg <- c(msg, paste0("Only TOM dissimilarity (dissFunc = 'TOMdiss') ",
                           "possible for soft thresholding."))
    }
  }

  if(length(msg) != 0 & verbose > 0){
    msg_new <- as.vector(rbind(msg, rep("\n", length(msg))))
    message("Infos about changed arguments:")
    message(msg_new)

  }

  #-----------------------------------------------------------------------------

  output <- list(zeroMethod = zeroMethod, normMethod = normMethod,
                 sparsMethod = sparsMethod, dissFunc = dissFunc,
                 sampleSize = sampleSize, needfrac = needfrac, needint = needint)
}
