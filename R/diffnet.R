#' @title Constructing Differential Networks for Microbiome Data
#'
#' @description Constructs a differential network for objects of class
#'   \code{microNet}. Three methods for identifying differentially associated
#'   taxa are provided: Fisher's z-test, permutation test, discordant method.
#'
#' @details \strong{Permutation procedure:}\cr
#'   The null hypothesis of these tests is defined as
#'   \deqn{H_0: a1_ij - a2_ij = 0,} where \eqn{a1_ij} and \eqn{a2_ij} denote the
#'   association between taxon i and j in group 1 and 2, respectively.\cr
#'   To generate a sampling distribution of the differences under \eqn{H_0},
#'   the group labels are randomly reassigned to the samples while the group
#'   sizes are kept. The associations are then re-estimated for each permuted
#'   data set. The p-values are calculated as the proportion of
#'   "permutation-differences" being larger than the observed difference. A
#'   pseudo-count is added to the numerator and denominator in order to avoid
#'   zero p-values. The p-values should be adjusted for multiple testing.
#'
#' @param x an object of class \code{microNet} inheriting from a call to
#'   \code{\link{netConstruct}}
#' @param diffMethod character string indicating the method used for determining
#'   differential associations. Possible values are \code{"permute"} (default)
#'   for performing permutation tests according to \cite{Gill et al. (2010)},
#'   \code{"discordant"}, which calls \code{\link[discordant]{discordantRun}}
#'   \cite{(Siska and Kechris, 2016)}, and \code{"fisherTest"} for Fisher's
#'   z-test \cite{(Fisher , 1992)}.
#' @param discordThresh numeric value in [0,1].
#'   Only used for the discordant method. Specifies
#'   a threshold for the posterior probability that a pair of
#'   taxa is differentially correlated between the groups. Taxa pairs with a
#'   posterior above this threshold are connected in the network. Defaults to
#'   0.8.
#' @param fisherTrans logical. If \code{TRUE} (default), Fisher-transformed
#'   correlations are used for permutation tests.
#' @param nPerm integer giving the number of permutations for the permutation
#'   tests. Defaults to 1000L.
#' @param permPvalsMethod character indicating the method used for determining
#'   p-values for permutation tests. Currently, "pseudo" is the only available
#'   option (see details).
#' @param cores integer indicating the number of CPU cores used for
#'   permutation tests. If cores > 1, the tests are performed parallel.
#'   Is limited to the number of available CPU cores determined by
#'   \code{\link[parallel]{detectCores}}. Defaults to 1L (no parallelization).
#' @param verbose logical. If \code{TRUE} (default), progress messages are shown.
#' @param logFile a character string defining the name of a log file, which is
#'   created when permutation tests are conducted (therein the current iteration
#'   numbers are stored). Defaults to \code{NULL} so that no file is created.
#' @param seed integer giving a seed for reproducibility of the results.
#' @param alpha numeric value between 0 and 1 giving the significance level.
#'   Significantly different correlations are connected in the network. Defaults
#'   to 0.05.
#' @param adjust character indicating the method used for multiple testing
#'   adjustment for the tests for differentially correlated pairs of taxa.
#'   Possible values are "lfdr" (default) for local
#'   false discovery rate correction (via \code{\link[fdrtool]{fdrtool}}),
#'   "adaptBH" for the adaptive Benjamini-Hochberg method \cite{(Benjamini and
#'   Hochberg, 2000)}, or one of the methods provided by
#'   \code{\link[stats]{p.adjust}}.
#' @param lfdrThresh defines a threshold for the local fdr if "lfdr" is chosen
#'   as method for multiple testing correction. Defaults to 0.2 meaning
#'   that correlations with a corresponding local fdr less than or equal to 0.2
#'   are identified as significant.
#' @param trueNullMethod character indicating the method used for estimating the
#'   proportion of true null hypotheses from a vector of p-values. Used for the
#'   adaptive Benjamini-Hochberg method for multiple testing adjustment (chosen
#'   by \code{adjust = "adaptBH"}). Accepts the provided options of the
#'   \code{method} argument of \code{\link[limma]{propTrueNull}}: "convest"
#'   (default), "lfdr", "mean", and "hist". Can alternatively be "farco" for
#'   the "iterative plug-in method" proposed by \cite{Farcomeni (2007)}.
#' @param pvalsVec vector with p-values used for permutation tests. Can be used
#'   for performing another method for multiple testing adjustment without
#'   executing the complete permutation process again. See the example.
#' @param fileLoadAssoPerm character giving the name (without extension) 
#'   or path of the file storing the "permuted" association/dissimilarity 
#'   matrices that have been exported by setting \code{storeAssoPerm} to 
#'   \code{TRUE}. Only used for permutation tests. Set to \code{NULL} if no 
#'   existing associations should be used.
#' @param fileLoadCountsPerm character giving the name (without extension) 
#'   or path of the file storing the "permuted" count matrices that have been 
#'   exported by setting \code{storeCountsPerm} to \code{TRUE}. 
#'   Only used for permutation tests, and if \code{fileLoadAssoPerm = NULL}. 
#'   Set to \code{NULL} if no existing count matrices should be used.
#' @param storeAssoPerm logical indicating whether the association (or
#'   dissimilarity) matrices for the permuted data should be stored in a file.
#'   The filename is given via \code{fileStoreAssoPerm}. If \code{TRUE}, 
#'   the computed "permutation" association/dissimilarity matrices can be reused
#'   via \code{fileLoadAssoPerm} to save runtime. Defaults to \code{FALSE}.
#' @param fileStoreAssoPerm character giving the file name to store a matrix
#'   containing a matrix with associations/dissimilarities for the permuted 
#'   data. Can also be a path.
#' @param storeCountsPerm logical indicating whether the permuted count matrices
#'   should be stored in an external file. Defaults to \code{FALSE}.
#' @param fileStoreCountsPerm character vector with two elements giving the 
#'   names of two files storing the permuted count matrices belonging to the 
#'   two groups.
#' @param assoPerm only needed for output generated with NetCoMi v1.0.1! A 
#'   list with two elements used for the permutation procedure.
#'   Each entry must contain association matrices for \code{"nPerm"}
#'   permutations. This can be either the \code{"assoPerm"} value as part of the
#'   output returned from \code{diffnet} or from \code{\link{netCompare}}. See
#'   the example.
#' @return The function returns an object of class \code{diffnet}. Depending on
#'   the performed test method, the output contains the following
#'   elements:\cr\cr
#'   \strong{Permutation tests:}
#'   \tabular{ll}{
#'   \code{diffMat}\tab matrix with absolute differences of associations that are
#'   significantly different from zero; optional adjacency matrix\cr
#'   \code{diffAdjustMat}\tab matrix with absolute differences of associations 
#'   that are significantly different from zero (after multiple testing 
#'   correction); optional adjacency matrix\cr
#'   \code{pvalsVec}\tab vector with p-values\cr
#'   \code{pAdjustVec}\tab vector with adjusted p-values\cr
#'   \code{pvalsMat}\tab matrix with p-values\cr
#'   \code{pAdjustMat}\tab matrix with adjusted p-values\cr
#'   \code{testStatData}\tab vector with test statistics (absolute differences 
#'   of associations) for the original data\cr
#'   \code{testStatPerm}\tab matrix with test statistics (absolute differences 
#'   of associations) for the permuted data\cr
#'   \code{assoMat1,assoMat2}\tab matrices with estimated associations (of the
#'   original data)}
#'   \strong{Discordant:}
#'   \tabular{ll}{
#'   \code{assoMat1,assoMat2}\tab matrices with estimated correlations\cr
#'   \code{diffMat}\tab adjacency matrix (absolute difference of correlations)\cr
#'   \code{classMat}\tab matrix with classes assigned to a taxa pair\cr
#'   \code{diffProbs}\tab matrix with posterior probabilities that a taxa pair
#'   is differentially correlated between the groups}
#'   \strong{Fisher's z-test:}
#'   \tabular{ll}{
#'   \code{diffMat}\tab matrix with absolute differences of associations that are
#'   significantly different from zero; optional adjacency matrix\cr
#'   \code{diffAdjustMat}\tab matrix with absolute differences of associations 
#'   that are significantly different from zero (after multiple testing 
#'   correction); optional adjacency matrix\cr
#'   \code{pvalsVec}\tab vector with p-values\cr
#'   \code{pAdjustVec}\tab vector with adjusted p-values\cr
#'   \code{pvalsMat}\tab matrix with p-values\cr
#'   \code{pAdjustMat}\tab matrix with adjusted p-values\cr
#'   \code{assoMat1,assoMat2}\tab matrices with estimated associations (of the
#'   original data)}
#'
#' @examples
#' # Load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#' 
#' # Generate a random group vector
#' set.seed(123456)
#' group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)
#' 
#' # Network construction:
#' amgut_net <- netConstruct(amgut1.filt, group = group,
#'                           measure = "pearson",
#'                           filtTax = "highestVar",
#'                           filtTaxPar = list(highestVar = 30),
#'                           zeroMethod = "pseudo", normMethod = "clr")
#' 
#' #---------------------
#' # Differential network
#' 
#' # Fisher's z-test
#' amgut_diff1 <- diffnet(amgut_net, diffMethod = "fisherTest")
#' 
#' # Network contains no differentially correlated taxa:
#' \dontrun{
#'   plot(amgut_diff1)
#' }
#' 
#' # Without multiple testing correction (statistically not correct!)
#' amgut_diff2 <- diffnet(amgut_net, diffMethod = "fisherTest", adjust = "none")
#' plot(amgut_diff2)
#' 
#' \dontrun{
#'   # Permutation test (permutation matrices are stored)
#'   amgut_diff3 <- diffnet(amgut_net, diffMethod = "permute", nPerm = 1000L,
#'                          cores = 4L, adjust = "lfdr",
#'                          storeCountsPerm = TRUE, 
#'                          fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
#'                          storeAssoPerm = TRUE,
#'                          fileStoreAssoPerm = "assoPerm",
#'                          seed = 123456)
#'   
#'   # Use the p-values again (different adjustment method possible), but without
#'   # re-estimating the associations
#'   amgut_diff4 <- diffnet(amgut_net, diffMethod = "permute", nPerm = 1000L,
#'                          adjust = "none", pvalsVec = amgut_diff3$pvalsVec)
#'   x11()
#'   plot(amgut_diff4)
#'   
#'   # Use the permutation associations again (same result as amgut_diff4)
#'   amgut_diff5 <- diffnet(amgut_net, diffMethod = "permute", nPerm = 1000L,
#'                          adjust = "none", 
#'                          fileLoadAssoPerm = "assoPerm")
#'   x11()
#'   plot(amgut_diff5)
#'   
#'   # Use the permuted count matrices again (same result as amgut_diff4)
#'   amgut_diff6 <- diffnet(amgut_net, diffMethod = "permute", nPerm = 1000L,
#'                          adjust = "none", 
#'                          fileLoadCountsPerm = c("countsPerm1", "countsPerm2"),
#'                          seed = 123456)
#'   x11()
#'   plot(amgut_diff6)
#' }
#'
#' @seealso \code{\link{plot.diffnet}}
#' @references \insertRef{benjamini2000adaptive}{NetCoMi} \cr
#'   \insertRef{discordant2016}{NetCoMi} \cr
#'   \insertRef{farcomeni2007some}{NetCoMi} \cr
#'   \insertRef{fisher1992statistical}{NetCoMi} \cr
#'   \insertRef{gill2010statistical}{NetCoMi}
#' @importFrom Biobase ExpressionSet
#' @importFrom Rdpack reprompt
#' @importFrom stats p.adjust.methods
#' @importFrom stats pnorm
#' @export

diffnet <- function(x, diffMethod = "permute", discordThresh = 0.8,
                    fisherTrans = TRUE, nPerm = 1000L,
                    permPvalsMethod = "pseudo",
                    cores = 1L, verbose = TRUE, logFile = NULL,
                    seed = NULL, alpha = 0.05, adjust = "lfdr",
                    lfdrThresh = 0.2, trueNullMethod = "convest",
                    pvalsVec = NULL,
                    fileLoadAssoPerm  = NULL,
                    fileLoadCountsPerm = NULL,
                    storeAssoPerm = FALSE,
                    fileStoreAssoPerm = "assoPerm",
                    storeCountsPerm = FALSE,
                    fileStoreCountsPerm = c("countsPerm1", "countsPerm2"), 
                    assoPerm = NULL){

  stopifnot(class(x) == "microNet")

  if(x$assoType == "dissimilarity"){
    stop("Differential network not implemented for dissimilarity-based networks.")
  }

  if(is.null(x$groups)){
    stop("'net' is a single network. A group vector must be passed to 'NetConstruct()'
         for network comparison.")
  }

  diffMethod <- match.arg(diffMethod, choices = c("discordant", "permute",
                                                  "fisherTest"))
  if(discordThresh < 0 || discordThresh > 1){
    stop("'discordThresh' must be in [0,1].")
  }

  nPerm <- as.integer(nPerm)

  if(permPvalsMethod != "pseudo") permPvalsMethod <- "pseudo"

  if(verbose %in% c(0,1)){
    verbose <- as.logical(verbose)
  } else{
    stopifnot(is.logical(verbose))
  }

  if(!is.null(logFile)) stopifnot(is.character(logFile))

  stopifnot(is.numeric(alpha))
  if(alpha < 0){
    stop("Significance level 'alpha' must be in [0,1].")
  }
  if(alpha > 1){
    alpha <- alpha / 100
    message("'alpha' transformed to ", alpha, ".")
  }

  adjust <- match.arg(adjust, c(p.adjust.methods, "lfdr", "adaptBH"))

  if(lfdrThresh < 0 || lfdrThresh > 1){
    stop("'lfdrThresh' must be in [0,1].")
  }

  trueNullMethod <- match.arg(trueNullMethod, c("convest", "lfdr", "mean",
                                                "hist", "farco"))
  
  if(diffMethod != "discordant" && adjust == "adaptBH" && 
     !requireNamespace("limma", quietly = TRUE)){
    
    message("Installing missing package 'limma' ...")
    
    if(!requireNamespace("BiocManager", quietly = TRUE)){
      utils::install.packages("BiocManager")
    }
    
    BiocManager::install("limma", dependencies = TRUE)
    message("Done.")
    
    message("Check whether installed package can be loaded ...")
    requireNamespace("limma")
    message("Done.")
  }

  #-----------------------------------------------------------------------------

  countMat1 <- x$countMat1
  countMat2 <- x$countMat2
  countsJoint <- x$countsJoint
  
  normCounts1 <- x$normCounts1
  normCounts2 <- x$normCounts2

  assoMat1 <- x$assoEst1
  assoMat2 <- x$assoEst2

  if(diffMethod == "discordant"){

    if(!requireNamespace("discordant", quietly = TRUE)){
      
      message("Installing missing package 'discordant' ...")
      
      if(!requireNamespace("BiocManager", quietly = TRUE)){
        utils::install.packages("BiocManager")
      }
      
      BiocManager::install("discordant", dependencies = TRUE)
      message("Done.")
      
      message("Check whether installed package can be loaded ...")
      requireNamespace("discordant")
      message("Done.")
    }
    
    if(!is.null(seed)) set.seed(seed)

    # create object of class ExpressionSet
    x_expr <- Biobase::ExpressionSet(assayData = t(rbind(countMat1, countMat2)))

    groups <- c(rep(1, nrow(countMat1)), rep(2, nrow(countMat2)))

    # transform correlation matrices to vectors
    lowtri <- lower.tri(assoMat1, diag = FALSE)
    corrVector1 <- assoMat1[lowtri]
    corrVector2 <- assoMat2[lowtri]
    vector_names <- get_vec_names(t(countMat1))
    names(corrVector1) <- vector_names
    names(corrVector2) <- vector_names

    # Erzeugen der Klassen mit Wahrscheinlichkeiten mittels 'discordant'-Methode
    discord <- discordant::discordantRun(corrVector1, corrVector2, x_expr)

    # Matrix mit zugeordneten Klassen (Klassen mit höchster Wsk.)
    classMat <- discord$classMatrix
    classMat[upper.tri(classMat)] <- t(classMat)[upper.tri(t(classMat))]
    diag(classMat) <- 1

    # Matrix mit Wsk. dass sich die Korrelationen zwischen den Gruppen unterscheiden
    # (entspricht aufsummierten Wahrscheinlichkeiten für Klassen 2,3,4,6,7,8)
    diffProbs <- discord$discordPPMatrix
    diffProbs[upper.tri(diffProbs)] <- t(diffProbs)[upper.tri(t(diffProbs))]
    diag(diffProbs) <- 0

    diffMat <- abs(assoMat1 - assoMat2)
    diffMat[diffProbs < discordThresh] <- 0
    output <- list(assoMat1 = assoMat1, assoMat2 = assoMat2,
                   diffMat = diffMat, classMat = classMat, diffProbs = diffProbs)

  } else if(diffMethod == "permute"){
    
    matchDesign <- x$matchDesign
    callNetConstr <- x$call

    pvalsVecInput <- pvalsVec

    if(is.null(pvalsVec)){

      cores <- as.integer(cores)

      if(cores > 1){
        if(parallel::detectCores() < cores) cores <- parallel::detectCores()
      }

      permResult <- permtest_diff_asso(countMat1 = countMat1,
                                       countMat2 = countMat2,
                                       countsJoint = countsJoint,
                                       normCounts1 = normCounts1, 
                                       normCounts2 = normCounts2,
                                       assoMat1 = assoMat1,
                                       assoMat2 = assoMat2,
                                       paramsNetConstruct = x$parameters,
                                       method = "connect.pairs",
                                       fisherTrans = fisherTrans,
                                       pvalsMethod = permPvalsMethod,
                                       adjust = adjust, adjust2 = "none",
                                       alpha = alpha, lfdrThresh = lfdrThresh,
                                       verbose = verbose, nPerm = nPerm,
                                       matchDesign = matchDesign,
                                       callNetConstr = callNetConstr,
                                       cores = cores, logFile = logFile,
                                       seed = seed, 
                                       fileLoadAssoPerm = fileLoadAssoPerm,
                                       fileLoadCountsPerm = fileLoadCountsPerm,
                                       storeAssoPerm = storeAssoPerm,
                                       fileStoreAssoPerm = fileStoreAssoPerm,
                                       storeCountsPerm = storeCountsPerm,
                                       fileStoreCountsPerm = fileStoreCountsPerm,
                                       assoPerm = assoPerm)

      pvalsVec <- permResult$pvalsVec
      pAdjustVec <- permResult$pAdjustVec
      
      testStatData <- permResult$testStatData
      testStatPerm <- permResult$testStatPerm

      nExceedsVec <- permResult$nExceedsVec

    } else{

      # adjust for multiple testing
      if(verbose & adjust != "none"){
        message("Adjust for multiple testing using '", adjust, "' ... ",
                appendLF = FALSE)
      }
      pAdjustVec <- multAdjust(pvals = pvalsVec, adjust = adjust,
                            trueNullMethod = trueNullMethod, verbose = verbose)
      if(verbose & adjust != "none") message("Done.")

      testStatData <- NULL
      testStatPerm <- NULL
      
      nExceedsVec <- NULL
    }

    diffMat <- diffAdjustMat <- abs(assoMat1 - assoMat2)
    diffVec <- diffVecAdjust <- diffMat[lower.tri(diffMat)]
    
    # identify links
    diffVec[pvalsVec > alpha] <- 0
    
    if(adjust == "none"){
      diffVecAdjust <- diffVec
    } else if(adjust == "lfdr"){
      diffVecAdjust[pAdjustVec > lfdrThresh] <- 0
    } else{
      diffVecAdjust[pAdjustVec > alpha] <- 0
    }

    diffMat[lower.tri(diffMat)] <- diffVec
    diffMat[upper.tri(diffMat)] <- t(diffMat)[upper.tri(t(diffMat))]
    
    diffAdjustMat[lower.tri(diffAdjustMat)] <- diffVecAdjust
    diffAdjustMat[upper.tri(diffAdjustMat)] <- 
      t(diffAdjustMat)[upper.tri(t(diffAdjustMat))]
    
    pvalsMat <- diffMat
    pvalsMat[lower.tri(pvalsMat)] <- pvalsVec
    pvalsMat[upper.tri(pvalsMat)] <- t(pvalsMat)[upper.tri(t(pvalsMat))]
    
    pAdjustMat <- diffMat
    pAdjustMat[lower.tri(pAdjustMat)] <- pAdjustVec
    pAdjustMat[upper.tri(pAdjustMat)] <- t(pAdjustMat)[upper.tri(t(pAdjustMat))]
    
    output = list()

    output[["diffMat"]] <- diffMat
    output[["diffAdjustMat"]] <- diffAdjustMat
    
    output[["pvalsVec"]] <- pvalsVec
    output[["pAdjustVec"]] <- pAdjustVec
    #output[["nExceedsVec"]] <- nExceedsVec
    
    output[["pvalsMat"]] <- pvalsMat
    output[["pAdjustMat"]] <- pAdjustMat
    
    output[["testStatData"]] <- testStatData
    output[["testStatPerm"]] <- testStatPerm
    
    output[["assoMat1"]] <- assoMat1
    output[["assoMat2"]] <- assoMat2

  }else{ #Fisher's z-test

    assoVec1 <- as.vector(assoMat1[lower.tri(assoMat1)])
    assoVec2 <- as.vector(assoMat2[lower.tri(assoMat2)])

    assoVec1[assoVec1 == 1] <- 0.9999
    assoVec2[assoVec2 == 1] <- 0.9999
    assoVec1[assoVec1 == -1] <- -0.9999
    assoVec2[assoVec2 == -1] <- -0.9999

    n1 <- nrow(normCounts1)
    n2 <- nrow(normCounts2)
    z1 <- atanh(assoVec1)
    z2 <- atanh(assoVec2)
    diff_z <- (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
    pvalsVec <- 2 * (1 - pnorm(abs(diff_z)))

    # adjust for multiple testing
    if(verbose & adjust != "none"){
      message("Adjust for multiple testing using '", adjust, "' ... ",
              appendLF = FALSE)
    }
    
    pAdjustVec <- multAdjust(pvals = pvalsVec, adjust = adjust,
                             trueNullMethod = trueNullMethod, verbose = verbose)
    
    if(verbose & adjust != "none") message("Done.")

    diffMat <- diffAdjustMat <- abs(assoMat1 - assoMat2)
    diag(diffMat) <- diag(diffAdjustMat) <- 0
    
    diffVec <- diffVecAdjust <- diffMat[lower.tri(diffMat)]
    
    diffVec[pvalsVec > alpha] <- 0

    # identify links
    if(adjust == "none"){
      diffVecAdjust <- diffVec
      
    } else if(adjust == "lfdr"){
      diffVecAdjust[pAdjustVec > lfdrThresh] <- 0
    } else{
      diffVecAdjust[pAdjustVec > alpha] <- 0
    }

    diffMat[lower.tri(diffMat)] <- diffVec
    diffMat[upper.tri(diffMat)] <- t(diffMat)[upper.tri(t(diffMat))]
    
    diffAdjustMat[lower.tri(diffAdjustMat)] <- diffVecAdjust
    diffAdjustMat[upper.tri(diffAdjustMat)] <- 
      t(diffAdjustMat)[upper.tri(t(diffAdjustMat))]
    
    
    pvalsMat <- diffMat
    pvalsMat[lower.tri(pvalsMat)] <- pvalsVec
    pvalsMat[upper.tri(pvalsMat)] <- t(pvalsMat)[upper.tri(t(pvalsMat))]
    
    pAdjustMat <- diffMat
    pAdjustMat[lower.tri(pAdjustMat)] <- pAdjustVec
    pAdjustMat[upper.tri(pAdjustMat)] <- t(pAdjustMat)[upper.tri(t(pAdjustMat))]
    
    output <- list()
    
    output[["diffMat"]] <- diffMat
    
    output[["diffAdjustMat"]] <- diffAdjustMat
    
    output[["pvalsVec"]] <- pvalsVec
    output[["pAdjustVec"]] <- pAdjustVec
    #output[["nExceedsVec"]] <- nExceedsVec
    
    output[["pvalsMat"]] <- pvalsMat
    output[["pAdjustMat"]] <- pAdjustMat
    
    output[["assoMat1"]] <- assoMat1
    output[["assoMat2"]] <- assoMat2

  }

  if(verbose & all(diffMat == 0)){
    message("No differentially associated taxa detected.")
  }

  output[["groups"]] <- x$groups
  output[["diffMethod"]] <- diffMethod
  output[["call"]] <- match.call()
  class(output) <- "diffnet"
  return(output)
}





