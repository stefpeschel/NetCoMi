#' @title Group Comparison of Network Properties
#'
#' @description Calculate and compare network properties for microbial networks
#'   using Jaccard's index, the Rand index, and permutation tests.
#'
#' @param x object of class \code{"microNetProps"} (inheriting from a call to
#'   \code{\link[NetCoMi]{netAnalyze}}).
#' @param permTest logical. If \code{TRUE}, a permutation test is conducted in
#'   order to test centrality measures and global network properties for group
#'   differences. Defaults to \code{FALSE}. May lead to a considerably increased
#'   execution time!
#' @param lnormFit logical indicating whether a log-normal distribution should
#'   be fitted to the calculated centrality values for determining Jaccard's
#'   index (see details). If \code{NULL} (default), the value is adopted
#'   from the input, so that the same method is used which has already been used
#'   for determining the hub nodes.
#' @param jaccQuant numeric value between 0 and 1 specifying the quantile which
#'   is used as threshold to identify the most central nodes for each centrality
#'   measure. The resulting sets of nodes are used to calculate Jaccard's index
#'   (see details). Defaults to 0.75.
#' @param nPerm number of permutations.
#' @param adjust character indicating the method used for multiple testing
#'   adjustment of the permutation p-values. Possible values are \code{"lfdr"}
#'   (default) for local false discovery rate correction (via
#'   \code{\link[fdrtool]{fdrtool}}), \code{"adaptBH"} for the adaptive
#'   Benjamini-Hochberg method \cite{(Benjamini and Hochberg, 2000)},  or one of
#'   the methods provided by \code{\link[stats]{p.adjust}} (see
#'   \code{p.adjust.methods()}).
#' @param trueNullMethod character indicating the method used for estimating the
#'   proportion of true null hypotheses from a vector of p-values. Used for the
#'   adaptive Benjamini-Hochberg method for multiple testing adjustment (chosen
#'   by \code{adjust = "adaptBH"}). Accepts the provided options of the
#'   \code{method} argument of \code{\link[limma]{propTrueNull}}:
#'   \code{"convest"}(default), \code{"lfdr"}, \code{"mean"}, and \code{"hist"}.
#'   Can alternatively be \code{"farco"} for the "iterative plug-in method"
#'   proposed by \cite{Farcomeni (2007)}.
#' @param nPermRand number of permutations used for testing the Rand index for
#'   being significantly different from a random assignment of nodes to the
#'   clusters. Execution time is not significantly increased, even for a high
#'   number of permutations. Defaults to 1000L.
#' @param cores integer indicating the number of CPU cores used for
#'   permutation tests. If cores > 1, the tests are performed parallel.
#'   Is limited to the number of available CPU cores determined by
#'   \code{\link[parallel]{detectCores}}. Defaults to 1L (no parallelization).
#' @param logFile character string naming the log file within which the current
#'   iteration number is stored (if permutation tests are performed). Defaults
#'   to \code{NULL} so that no log file is generated.
#' @param seed integer giving a seed for reproducibility of the results.
#' @param verbose logical. If \code{TRUE} (default), status messages are shown.
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
#' @param storeAssoPerm logical indicating whether the association/dissimilarity 
#'   matrices for the permuted data should be stored in a file.
#'   The filename is given via \code{fileStoreAssoPerm}. If \code{TRUE}, 
#'   the computed "permutation" association/dissimilarity matrices can be reused
#'   via \code{fileLoadAssoPerm} to save runtime. Defaults to \code{FALSE}.
#'   Ignored if \code{fileLoadAssoPerm} is not \code{NULL}.
#' @param fileStoreAssoPerm character giving the file name to store a matrix
#'   containing a matrix with associations/dissimilarities for the permuted 
#'   data. Can also be a path.
#' @param storeCountsPerm logical indicating whether the permuted count matrices
#'   should be stored in an external file. Defaults to \code{FALSE}.
#'   Ignored if \code{fileLoadCountsPerm} is not \code{NULL}.
#' @param fileStoreCountsPerm character vector with two elements giving the 
#'   names of two files storing the permuted count matrices belonging to the 
#'   two groups.
#' @param returnPermProps logical. If \code{TRUE}, the global properties and 
#'   their absolute differences for the permuted data are returned.
#' @param returnPermCentr logical. If \code{TRUE}, the centralities and 
#'   their absolute differences for the permuted data are returned.
#' @param assoPerm only needed for output generated with NetCoMi v1.0.1! A list 
#'   with two elements used for the permutation procedure.
#'   Each entry must contain association matrices for \code{"nPerm"}
#'   permutations. This can be the \code{"assoPerm"} value as part of the
#'   output either returned from \code{diffnet} or from
#'   \code{\link{netCompare}}. 
#' @param dissPerm only needed for output generated with NetCoMi v1.0.1! 
#'   Usage analog to \code{assoPerm} if a dissimilarity measure has been used 
#'   for network construction.
#' @details
#'   \strong{Permutation procedure:}\cr
#'   Used for testing centrality measures and global network properties for
#'   group differences.\cr
#'   The null hypothesis of these tests is defined as
#'   \deqn{H_0: c1_i - c2_i = 0,} where \eqn{c1_i} and
#'   \eqn{c2_i} denote the centrality values of taxon i in group 1 and 2,
#'   respectively.\cr
#'   To generate a sampling distribution of the differences under \eqn{H_0},
#'   the group labels are randomly reassigned to the samples while the group
#'   sizes are kept. The associations are then re-estimated for each permuted
#'   data set. The p-values are calculated as the proportion of
#'   "permutation-differences" being larger than the observed difference. A
#'   pseudo-count is added to the numerator and denominator in order to avoid
#'   zero p-values. The p-values should be adjusted for multiple testing.
#'
#'   \strong{Jaccard's index:}\cr
#'   Jaccard's index expresses for each centrality measure how equal the sets of
#'   most central nodes are among the two networks.\cr
#'   These sets are defined as nodes with a centrality value above a defined
#'   quantile (via \code{jaccQuant}) either of the empirical distribution of the
#'   centrality values (\code{lnormFit = FALSE}) or of a fitted log-normal
#'   distribution (\code{lnormFit = TRUE}).\cr
#'   The index ranges from 0 to 1, where 1 means the sets of most central nodes
#'   are exactly equal in both networks and a value of 0 implicates that the
#'   most central nodes are completely different.\cr
#'   The index is calculated as suggested by \cite{Real and Vargas (1996)}.
#'   \cr\cr
#'   \strong{Rand index:}
#'   The Rand index is used to express whether the determined clusterings are
#'   equal in both groups. The adjusted Rand index (ARI) ranges from -1 to 1,
#'   where 1 indicates that the two clusterings are exactly equal. The expected
#'   index value for two random clusterings is 0. The implemented test procedure
#'   is in accordance with the explanations in \cite{Qannari et al. (2014)},
#'   where a p-value below the alpha levels means that ARI is significantly
#'   higher than expected for two random clusterings.
#' @return Returned is an object of class \code{microNetComp} with the following
#'   elements:\cr
#'   \tabular{ll}{
#'   \code{jaccDeg,jaccBetw,jaccClose,jaccEigen}\tab Values of Jaccard's index
#'   for the centrality measures\cr
#'   \code{jaccHub}\tab Jaccard index for the sets of hub nodes\cr
#'   \code{randInd}\tab Calculated Rand index\cr
#'   \code{properties}\tab List with calculated network properties\cr
#'   \code{propertiesLCC}\tab List with calculated network properties of the 
#'   largest connected component (LCC)\cr
#'   \code{diffGlobal}\tab Vectors with differences of global properties\cr
#'   \code{diffGlobalLCC}\tab Vectors with differences of global properties for
#'   the LCC\cr
#'   \code{diffCent}\tab Vectors with differences of the centrality values\cr
#'   \code{countMatrices}\tab The two count matrices returned
#'   from \code{netConstruct}\cr
#'   \code{assoMatrices}\tab The two association matrices returned
#'   from \code{netConstruct}\cr
#'   \code{dissMatrices}\tab The two dissimilarity matrices returned
#'   from \code{netConstruct}\cr
#'   \code{adjaMatrices}\tab The two adjacency matrices returned
#'   from \code{netConstruct}\cr
#'   \code{groups}\tab Group names returned from \code{netConstruct}\cr
#'   \code{paramsProperties}\tab Parameters used for network analysis
#'   }
#'   \strong{Additional output if permutation tests are conducted:}
#'  \tabular{ll}{
#'  \code{pvalDiffGlobal}\tab P-values of the tests for differential global
#'  properties\cr
#'  \code{pvalDiffGlobalLCC}\tab P-values of the tests for differential
#'  global properties in the LCC\cr
#'  \code{pvalDiffCentr}\tab P-values of the tests for differential centrality
#'  values\cr
#'  \code{pvalDiffCentrAdjust}\tab Adjusted p-values of the tests for
#'  differential centrality values\cr
#'  \code{permDiffGlobal}\tab \code{nPerm} x 10 matrix containing the absolute 
#'  differences of the ten global network properties (computed for the whole 
#'  network) for all \code{nPerm} permutations\cr
#'  \code{permDiffGlobalLCC}\tab \code{nPerm} x 11 matrix containing the 
#'  absolute differences of the eleven global network properties (computed for 
#'  the LCC) for all \code{nPerm} permutations\cr
#'  \code{permDiffCentr}\tab List with absolute differences of the four 
#'  centrality measures for all \code{nPerm} permutations. Each list contains 
#'  a \code{nPerm} x \code{nNodes} matrix.
#'  }
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
#' # Network analysis:
#' amgut_props <- netAnalyze(amgut_net, clustMethod = "cluster_fast_greedy")
#' 
#' # Network plot:
#' plot(amgut_props, sameLayout = TRUE)
#' 
#' #--------------------------
#' # Network comparison:
#' 
#' # Without permutation tests:
#' amgut_comp1 <- netCompare(amgut_props, permTest = FALSE)
#' summary(amgut_comp1)
#' 
#' \donttest{
#'   # With permutation tests (with only 100 permutations to decrease runtime):
#'   amgut_comp2 <- netCompare(amgut_props, 
#'                             permTest = TRUE, nPerm = 100L, cores = 1L, 
#'                             storeCountsPerm = TRUE, 
#'                             fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
#'                             storeAssoPerm = TRUE, 
#'                             fileStoreAssoPerm = "assoPerm",
#'                             seed = 123456)
#'   
#'   # Rerun with a different adjustment method ...
#'   # ... using the stored permutation count matrices
#'   amgut_comp3 <- netCompare(amgut_props, adjust = "BH",
#'                             permTest = TRUE, nPerm = 100L, 
#'                             fileLoadCountsPerm = c("countsPerm1", "countsPerm2"),
#'                             seed = 123456)
#'   
#'   # ... using the stored permutation association matrices
#'   amgut_comp4 <- netCompare(amgut_props, adjust = "BH",
#'                             permTest = TRUE, nPerm = 100L, 
#'                             fileLoadAssoPerm = "assoPerm",
#'                             seed = 123456)
#'   
#'   # amgut_comp3 and amgut_comp4 should be equal
#'   all.equal(amgut_comp3$adjaMatrices, amgut_comp4$adjaMatrices)
#'   all.equal(amgut_comp3$properties, amgut_comp4$properties)
#'   
#'   summary(amgut_comp2)
#'   summary(amgut_comp3)
#'   summary(amgut_comp4)
#'   
#'   #--------------------------
#'   # Use 'createAssoPerm' to create "permuted" count and association matrices
#'   createAssoPerm(amgut_props, nPerm = 100, 
#'                  computeAsso = TRUE,
#'                  fileStoreAssoPerm = "assoPerm",
#'                  storeCountsPerm = TRUE, 
#'                  fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
#'                  append = FALSE, seed = 123456)
#'   
#'   amgut_comp5 <- netCompare(amgut_props, permTest = TRUE, nPerm = 100L, 
#'                             fileLoadAssoPerm = "assoPerm")
#'   
#'   all.equal(amgut_comp3$properties, amgut_comp5$properties)
#'   
#'   summary(amgut_comp5)
#' }
#'
#' @seealso \code{\link{summary.microNetComp}}, \code{\link{netConstruct}},
#'   \code{\link{netAnalyze}}
#' @references \insertRef{benjamini2000adaptive}{NetCoMi} \cr
#'   \insertRef{farcomeni2007some}{NetCoMi} \cr
#'   \insertRef{gill2010statistical}{NetCoMi} \cr
#'   \insertRef{qannari2014significance}{NetCoMi} \cr
#'   \insertRef{real1996probabilistic}{NetCoMi}
#' @import foreach doSNOW filematrix
#' @importFrom snow makeCluster stopCluster
#' @importFrom stats sd
#' @importFrom WGCNA randIndex
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

netCompare <- function(x,
                       permTest = FALSE,
                       lnormFit = NULL,
                       jaccQuant = 0.75,
                       nPerm = 1000L,
                       adjust = "adaptBH",
                       trueNullMethod = "convest",
                       nPermRand = 1000L,
                       cores = 1L,
                       logFile = NULL,
                       seed = NULL, 
                       verbose = TRUE,
                       fileLoadAssoPerm  = NULL,
                       fileLoadCountsPerm = NULL,
                       storeAssoPerm = FALSE,
                       fileStoreAssoPerm = "assoPerm",
                       storeCountsPerm = FALSE,
                       fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
                       returnPermProps = FALSE,
                       returnPermCentr = FALSE,
                       assoPerm = NULL, dissPerm = NULL){

  stopifnot(class(x) =="microNetProps")
  stopifnot(is.logical(permTest))
  if(!is.null(lnormFit)) stopifnot(is.logical(lnormFit))
  stopifnot(jaccQuant >= 0 & jaccQuant <= 1)
  nPerm <- as.integer(nPerm)

  adjust <- match.arg(adjust, c(p.adjust.methods, "lfdr", "adaptBH"))
  adjustPerm <- adjust

  trueNullMethod <- match.arg(trueNullMethod, c("convest", "lfdr", "mean",
                                                "hist", "farco"))
  
  nPermRand <- as.integer(nPermRand)
  cores <- as.integer(cores)
  if(!is.null(logFile)) stopifnot(is.character(logFile))
  if(!is.null(assoPerm)) stopifnot(is.list(assoPerm) & length(assoPerm) == 2)

  if(!is.null(dissPerm)) stopifnot(is.list(dissPerm) & length(dissPerm) == 2)

  #-----------------------------------------------------------------------------
  avDissIgnoreInf <- sPathNorm <- sPathAlgo <- connectivity <- normNatConnect <- 
  weighted <- clustMethod <- clustPar <- clustPar2 <- weightClustCoef <- 
  hubPar <- hubQuant <- weightDeg <- normDeg <- normBetw <- normClose <- 
  normEigen <- centrLCC <- softThreshPower <- softThreshCut <- measure <- 
  measurePar <- scaleDiss <- sparsMethod <- thresh <- alpha <- lfdrThresh <- 
  nboot <- softThreshType <- kNeighbor <- knnMutual <- dissFunc <- 
  dissFuncPar <- simFunc <- simFuncPar <- pvalavDiss <- pvalavPath <- 
  pvalDensity <- pvalEdgeConnect <- pvalVertConnect <- pvalNatConnect <- 
  pvalPEP <- pvalClustCoef <- pvalModul <- pvalavDiss_lcc <- 
  pvalavPath_lcc <- pvalDensity_lcc <- pvalEdgeConnect_lcc <- 
  pvalVertConnect_lcc <- pvalNatConnect_lcc <- pvalPEP_lcc <- 
  pvalClustCoef_lcc <- pvalModul_lcc <- p <- zeroMethod <- zeroPar <- 
  needfrac <- needint <- jointPrepro <- normMethod <- normPar <- NULL
  
  for(i in 1:length(x$paramsProperties)){
    assign(names(x$paramsProperties)[i], x$paramsProperties[[i]])
  }

  parnames <- c("zeroMethod", "zeroPar", "needfrac", "needint", "jointPrepro",
                "normMethod", "normPar", "measure", "measurePar", 
                "sparsMethod", "thresh", "adjust", "alpha", "lfdrThresh", 
                "nboot", "softThreshType", "softThreshCut", "softThreshPower", 
                "kNeighbor", "knnMutual", "dissFunc", "dissFuncPar",
                "simFunc", "simFuncPar", "scaleDiss", "weighted")

  for(i in 1:length(parnames)){
    assign(parnames[i], x$paramsNetConstruct[[parnames[i]]])
  }

  assoType <- x$input$assoType
  distNet <- ifelse(assoType == "dissimilarity", TRUE, FALSE)
  
  xgroups <- x$input$groups
  matchDesign <- x$input$matchDesign
  sampleSize <- x$input$sampleSize
  twoNets <- x$input$twoNets
  callNetConstr <- x$input$call

  count1 <- x$input$countMat1
  count2 <- x$input$countMat2
  countsJoint <- x$input$countsJoint
  
  count_norm1 <- x$input$normCounts1
  count_norm2 <- x$input$normCounts2

  assoMat1 <- x$input$assoMat1
  assoMat2 <- x$input$assoMat2

  dissMat1 <- x$input$dissMat1
  dissMat2 <- x$input$dissMat2

  adja1 <- x$input$adjaMat1
  adja2 <- x$input$adjaMat2

  #---------------------------------------------------------------------------

  if(!twoNets){
    stop("Network comparison not possible because a single network
                    has been constructed.")
  }

  if(all(adja1[upper.tri(adja1)] == 0)){
    stop(paste0("Network properties cannot be compared because network of group'",
                xgroups[1], "' is empty."))
  }
  if(all(adja2[upper.tri(adja2)] == 0)){
    stop(paste0("Network properties cannot be compared because network of group'",
                xgroups[2], "' is empty."))
  }

  #---------------------------------------------------------------------------
  # calculate and compare network properties
  
  if(!is.null(seed)) set.seed(seed)
  
  if(permTest & verbose) message("Calculate network properties ... ",
                                 appendLF = FALSE)

  props <- calc_diff_props(adja1 = adja1, adja2 = adja2,
                           dissMat1 = dissMat1, dissMat2 = dissMat2,
                           assoMat1 = assoMat1, assoMat2 = assoMat2,
                           avDissIgnoreInf = avDissIgnoreInf,
                           sPathNorm = sPathNorm, sPathAlgo = sPathAlgo,
                           connectivity = connectivity, 
                           normNatConnect = normNatConnect,
                           weighted = weighted, clustMethod = clustMethod, 
                           clustPar = clustPar, clustPar2 = clustPar2,
                           weightClustCoef = weightClustCoef,
                           hubPar = hubPar,  hubQuant = hubQuant,
                           jaccQuant = jaccQuant, lnormFit = lnormFit,
                           weightDeg = weightDeg,
                           normDeg = normDeg,  normBetw = normBetw,
                           normClose = normClose, normEigen = normEigen,
                           centrLCC = centrLCC, nPermRand = nPermRand)
  
  if(permTest & verbose) message("Done.")
  
  output <- list(jaccDeg = props$jaccDeg,
                 jaccBetw = props$jaccBetw,
                 jaccClose = props$jaccClose,
                 jaccEigen = props$jaccEigen,
                 jaccHub = props$jaccHub,
                 randInd = props$randInd,
                 diffGlobal = props$diffsGlobal,
                 diffGlobalLCC = props$diffsGlobalLCC,
                 diffCentr = props$diffsCentr,
                 properties = props$props,
                 propertiesLCC = props$propsLCC,
                 countMatrices = list(count1 = count1,
                                      count2 = count2),
                 assoMatrices = list(assoMat1 = assoMat1,
                                     assoMat2 = assoMat2),
                 dissMatrices = list(dissMat1 = dissMat1,
                                     dissMat2 = dissMat2),
                 adjaMatrices = list(adja1 = adja1,
                                     adja2 = adja2),
                 groups = list(group1 = xgroups[1],
                               group2 = xgroups[2]),
                 paramsProperties = x$paramsProperties,
                 call = match.call())

  #-----------------------------------------------------------------------------
  # generate teststatistics for permutated data

  if(permTest){

    if(!is.null(seed)) set.seed(seed)
    
    if(assoType == "dissimilarity"){
      n1 <- ncol(count1)
      n2 <- ncol(count2)
      n <- n1 + n2
      xbind <- cbind(count1, count2)
      
    } else{
      if(jointPrepro){
        n1 <- nrow(count_norm1)
        n2 <- nrow(count_norm2)
        
        xbind <- rbind(count_norm1, count_norm2)
        
      } else{
        n1 <- nrow(count1)
        n2 <- nrow(count2)
        
        xbind <- rbind(count1, count2)
      }
      
      n <- n1 + n2
    }
    
    
    nVar = ncol(adja1)
    
    #---------------------------------------------------------------------------

    if(!is.null(fileLoadAssoPerm)){
      stopifnot(is.character(fileLoadAssoPerm))
      stopifnot(length(fileLoadAssoPerm) == 1)
      
      if(storeAssoPerm && identical(fileLoadAssoPerm, fileStoreAssoPerm)){
        storeAssoPerm <- FALSE
      }

    } else if(!is.null(fileLoadCountsPerm)){
      stopifnot(is.character(fileLoadCountsPerm))
      stopifnot(length(fileLoadCountsPerm) == 2)
      
      if(storeCountsPerm && identical(fileLoadCountsPerm, fileStoreCountsPerm)){
        storeCountsPerm <- FALSE
      }

    } else{
      perm_group_mat <- get_perm_group_mat(n1 = n1, n2 = n2, n = n, 
                                           nPerm = nPerm, 
                                           matchDesign = matchDesign)
    }
    
    #-----------------
    if(storeCountsPerm){
      stopifnot(is.character(fileStoreCountsPerm))
      stopifnot(length(fileStoreCountsPerm) == 2)
      
      fmat_counts1 <- fm.create(filenamebase = fileStoreCountsPerm[1], 
                                nrow = (n1 * nPerm), ncol = nVar)
      fmat_counts2 <- fm.create(filenamebase = fileStoreCountsPerm[2], 
                                nrow = (n2 * nPerm), ncol = nVar)
      
      if(verbose){
        message("Files '", 
                paste0(fileStoreCountsPerm[1], ".bmat, "),
                paste0(fileStoreCountsPerm[1], ".desc.txt, \n  "),
                paste0(fileStoreCountsPerm[2], ".bmat, and "),
                paste0(fileStoreCountsPerm[2], ".desc.txt created. "))
      }
      
      close(fmat_counts1)
      close(fmat_counts2)
    }
    
    #-----------------
    if(storeAssoPerm){
      stopifnot(is.character(fileStoreAssoPerm))
      stopifnot(length(fileStoreAssoPerm) == 1)
      
      fmat = fm.create(filenamebase = fileStoreAssoPerm, 
                       nrow = (nVar * nPerm), ncol = (2 * nVar))
      
      if(verbose){
        message("Files '", 
                paste0(fileStoreAssoPerm, ".bmat and "),
                paste0(fileStoreAssoPerm, ".desc.txt created. "))
      }
      
      close(fmat)
    }
    
    #---------------------------------------------------------------------------
    if(verbose) message("Execute permutation tests ... ")
    if(!is.null(seed)){
      seeds <- sample.int(1e8, size = nPerm)
    } else{
      seeds <- NULL
    }

    if(cores > 1){
      if(parallel::detectCores() < cores) cores <- parallel::detectCores()

      cl <- snow::makeCluster(cores, outfile = "")
      doSNOW::registerDoSNOW(cl)
      '%do_or_dopar%' <- get('%dopar%')

    } else{
      '%do_or_dopar%' <- get('%do%')
    }

    if(verbose){
      pb <- utils::txtProgressBar(0, nPerm, style=3)
      
      progress <- function(n){
        utils::setTxtProgressBar(pb,n)
      }
      
      opts <- list(progress=progress)
    } else{
      opts <- list()
    }

    if(!is.null(logFile)) cat("", file=logFile, append=FALSE)

    p <- NULL
    propsPerm <- foreach(p = 1:nPerm,
                         .packages = c("filematrix"),
                         .export = c("calc_association", "cclasso", "gcoda",
                                     "sparsify", "multAdjust", "trans_to_diss",
                                     "trans_to_sim", "trans_to_adja",
                                     "scale_diss", "norm_counts", "zero_treat",
                                     "calc_props", "calc_diff_props",
                                     "calc_jaccard"),
                         .options.snow = opts) %do_or_dopar% {

                           if(!is.null(logFile)){
                             cat(paste("iteration", p,"\n"),
                                 file=logFile, append=TRUE)
                           }
                           
                           if(verbose) progress(p)

                           if(!is.null(seed)) set.seed(seeds[p])

                           if(!is.null(assoPerm)){
                             # load permutation association matrices (old version)
                             assoMat1.tmp <- assoPerm[[1]][[p]]
                             assoMat2.tmp <- assoPerm[[2]][[p]]
                             count1.tmp <- count2.tmp <- NULL
                             
                           } else if(!is.null(dissPerm)){
                             # load permutation dissimilarity matrices (old version)
                             assoMat1.tmp <- dissPerm[[1]][[p]]
                             assoMat2.tmp <- dissPerm[[2]][[p]]
                             count2.tmp <- count2.tmp <- NULL
                             
                           } else if(!is.null(fileLoadAssoPerm)){
                             
                             fmat <- fm.open(filenamebase = fileLoadAssoPerm, 
                                             readonly = TRUE)
                             
                             # load permutation asso/diss matrices
                             assoMat1.tmp <- fmat[(p-1) * nVar + (1:nVar), 
                                                  1:nVar]
                             assoMat2.tmp <- fmat[(p-1) * nVar + (1:nVar), 
                                                  nVar + (1:nVar)]
                             
                             dimnames(assoMat1.tmp) <- dimnames(assoMat1)
                             dimnames(assoMat2.tmp) <- dimnames(assoMat2)
                             
                             count1.tmp <- count2.tmp <- NULL
                             
                             close(fmat)
                             
                           } else{
                             
                             if(!is.null(fileLoadCountsPerm)){
                               
                               fmat_counts1 <- fm.open(filenamebase = 
                                                         fileLoadCountsPerm[1], 
                                                       readonly = TRUE)
                               
                               fmat_counts2 <- fm.open(filenamebase = 
                                                         fileLoadCountsPerm[2], 
                                                       readonly = TRUE)
                               
                               # load permutation count matrices
                               count1.tmp <- 
                                 fmat_counts1[(p-1) * n1 + (1:n1), 1:nVar]
                               
                               count2.tmp <- 
                                 fmat_counts2[(p-1) * n2 + (1:n2), 1:nVar]
                               
                               close(fmat_counts1)
                               close(fmat_counts2)
                               
                             } else{
                               # generate permutation count matrices
                               if(assoType == "dissimilarity"){
                                 count1.tmp <- xbind[ ,which(perm_group_mat[p, ] == 1)]
                                 count2.tmp <- xbind[ ,which(perm_group_mat[p, ] == 2)]
                                 
                               } else{
                                 count1.tmp <- xbind[which(perm_group_mat[p, ] == 1), ]
                                 count2.tmp <- xbind[which(perm_group_mat[p, ] == 2), ]
                               }
                               
                               if(!jointPrepro){
                                 # zero treatment and normalization necessary if
                                 # in network construction two count matrices 
                                 # were given or dissimilarity network is created
                                 
                                 suppressMessages(
                                   count1.tmp <- zero_treat(countMat = count1.tmp, 
                                                            zeroMethod = zeroMethod,
                                                            zeroParam = zeroPar, 
                                                            needfrac = needfrac,
                                                            needint = needint, 
                                                            verbose = FALSE))
                                 
                                 suppressMessages(
                                   count2.tmp <- zero_treat(countMat = count2.tmp, 
                                                            zeroMethod = zeroMethod,
                                                            zeroParam = zeroPar, 
                                                            needfrac = needfrac,
                                                            needint = needint, 
                                                            verbose = FALSE))
                                 
                                 suppressMessages(
                                   count1.tmp <- norm_counts(countMat = count1.tmp, 
                                                             normMethod = normMethod,
                                                             normParam = normPar, 
                                                             zeroMethod = zeroMethod,
                                                             needfrac = needfrac, 
                                                             verbose = FALSE))
                                 
                                 suppressMessages(
                                   count2.tmp <- norm_counts(countMat = count2.tmp, 
                                                             normMethod = normMethod,
                                                             normParam = normPar, 
                                                             zeroMethod = zeroMethod,
                                                             needfrac = needfrac, 
                                                             verbose = FALSE))
                               }
                               
                               if(storeCountsPerm){
                                 fmat_counts1 <- fm.open(filenamebase = 
                                                           fileStoreCountsPerm[1])
                                 
                                 fmat_counts2 <- fm.open(filenamebase = 
                                                           fileStoreCountsPerm[2])
                                 
                                 fmat_counts1[(p-1) * n1 + (1:n1), 
                                              1:nVar] <- count1.tmp
                                 fmat_counts2[(p-1) * n2 + (1:n2), 
                                              1:nVar] <- count2.tmp
                                 
                                 close(fmat_counts1)
                                 close(fmat_counts2)
                               }
                             }
                             
                             assoMat1.tmp <- calc_association(count1.tmp,
                                                              measure = measure,
                                                              measurePar = measurePar,
                                                              verbose = FALSE)
                             assoMat2.tmp <- calc_association(count2.tmp,
                                                              measure = measure,
                                                              measurePar = measurePar,
                                                              verbose = FALSE)
                             
                             if(storeAssoPerm){
                               fmat <- fm.open(filenamebase = fileStoreAssoPerm)
                               
                               fmat[(p-1) * nVar + (1:nVar), 
                                    1:nVar] <- assoMat1.tmp
                               fmat[(p-1) * nVar + (1:nVar), 
                                    nVar + (1:nVar)] <- assoMat2.tmp
                               
                               close(fmat)
                             }
                           }
                           
                           #----------------------------------------------------
                           
                           if(distNet && scaleDiss){
                             assoMat1.tmp <- scale_diss(assoMat1.tmp)
                             assoMat2.tmp <- scale_diss(assoMat1.tmp)
                           }
                           
                           #----------------------------------------------------
                           # sparsification and transformation
                           
                           if(x$paramsNetConstruct$sparsMethod == "softThreshold"){
                             if(length(softThreshPower) < 2){
                               power1 <- power2 <- softThreshPower
                             } else{
                               power1 <- softThreshPower[1]
                               power2 <- softThreshPower[2]
                             }
                             
                             if(length(softThreshCut) < 2){
                               stcut1 <- stcut2 <- softThreshCut
                             } else{
                               stcut1 <- softThreshCut[1]
                               stcut2 <- softThreshCut[2]
                             }
                           }

                           sparsReslt <- sparsify(assoMat = assoMat1.tmp,
                                                  countMat = count1.tmp,
                                                  sampleSize = n1,
                                                  measure = measure,
                                                  measurePar = measurePar,
                                                  assoType = assoType,
                                                  sparsMethod = sparsMethod,
                                                  thresh = thresh[1],
                                                  alpha = alpha[1],
                                                  adjust = adjust,
                                                  lfdrThresh = lfdrThresh,
                                                  trueNullMethod = trueNullMethod,
                                                  nboot = nboot,
                                                  softThreshType = softThreshType,
                                                  softThreshPower = power1,
                                                  softThreshCut = stcut1,
                                                  kNeighbor = kNeighbor,
                                                  knnMutual = knnMutual,
                                                  parallel = FALSE,
                                                  cores = 1L,
                                                  verbose = FALSE)

                           if(distNet){
                             dissMat1.tmp <- sparsReslt$assoNew
                             assoMat1.tmp <- NULL
                           } else{
                             assoMat1.tmp <- sparsReslt$assoNew
                             dissMat1.tmp <- trans_to_diss(x = assoMat1.tmp,
                                                           dissFunc = dissFunc,
                                                           dissFuncPar = dissFuncPar)
                           }

                           simMat1.tmp <- trans_to_sim(x = dissMat1.tmp,
                                                       simFunc = simFunc,
                                                       simFuncPar = simFuncPar)
                           adja1.tmp <- trans_to_adja(x = simMat1.tmp,
                                                      weighted = weighted)
                           
                           # network 2
                           sparsReslt <- sparsify(assoMat = assoMat2.tmp,
                                                  countMat = count2.tmp,
                                                  sampleSize = n2,
                                                  measure = measure,
                                                  assoType = assoType,
                                                  sparsMethod = sparsMethod,
                                                  thresh = thresh[2],
                                                  alpha = alpha[2],
                                                  adjust = adjust,
                                                  lfdrThresh = lfdrThresh,
                                                  trueNullMethod = trueNullMethod,
                                                  nboot = nboot,
                                                  softThreshType = softThreshType,
                                                  softThreshPower = power2,
                                                  softThreshCut = stcut2,
                                                  kNeighbor = kNeighbor,
                                                  knnMutual = knnMutual,
                                                  parallel = FALSE,
                                                  cores = 1L,
                                                  verbose = FALSE)

                           if(distNet){
                             dissMat2.tmp <- sparsReslt$assoNew
                             assoMat2.tmp <- NULL
                           } else{
                             assoMat2.tmp <- sparsReslt$assoNew
                             dissMat2.tmp <- trans_to_diss(x = assoMat2.tmp,
                                                           dissFunc = dissFunc,
                                                           dissFuncPar = dissFuncPar)
                           }

                           simMat2.tmp <- trans_to_sim(x = dissMat2.tmp,
                                                       simFunc = simFunc,
                                                       simFuncPar = simFuncPar)
                           adja2.tmp <- trans_to_adja(x = simMat2.tmp,
                                                      weighted = weighted)
                           
                           dimnames(adja1.tmp) <- dimnames(adja1)
                           dimnames(adja2.tmp) <- dimnames(adja2)
                           dimnames(assoMat1.tmp) <- dimnames(assoMat1)
                           dimnames(assoMat2.tmp) <- dimnames(assoMat2)
                           dimnames(dissMat1.tmp) <- dimnames(dissMat1)
                           dimnames(dissMat2.tmp) <- dimnames(dissMat2)
                           
                           #----------------------------------------------------
                           # compute network properties

                           prop.tmp <- calc_diff_props(adja1 = adja1.tmp,
                                                       adja2 = adja2.tmp,
                                                       dissMat1 = dissMat1.tmp,
                                                       dissMat2 = dissMat2.tmp,
                                                       assoMat1 = assoMat1.tmp,
                                                       assoMat2 = assoMat2.tmp,
                                                       avDissIgnoreInf =
                                                         avDissIgnoreInf,
                                                       sPathNorm = sPathNorm,
                                                       sPathAlgo = sPathAlgo,
                                                       connectivity = connectivity,
                                                       normNatConnect = normNatConnect,
                                                       weighted = weighted,
                                                       clustMethod = clustMethod,
                                                       clustPar = clustPar,
                                                       clustPar2 = clustPar2,
                                                       weightClustCoef =
                                                         weightClustCoef,
                                                       hubPar = hubPar,
                                                       hubQuant = hubQuant,
                                                       jaccQuant = jaccQuant,
                                                       lnormFit = lnormFit,
                                                       weightDeg = weightDeg,
                                                       normDeg = normDeg,
                                                       normBetw = normBetw,
                                                       normClose = normClose,
                                                       normEigen = normEigen,
                                                       centrLCC = centrLCC,
                                                       nPermRand = nPermRand,
                                                       testJacc = FALSE,
                                                       testRand = FALSE)

                           out <- list(diffsGlobal = NULL,
                                       diffsGlobalLCC = NULL,
                                       absDiffsCentr = NULL)

                           for(i in c("diffavDiss", "diffavPath", "diffDensity",
                                      "diffVertConnect", "diffEdgeConnect",
                                      "diffNatConnect", "diffPEP",
                                      "diffClustCoef",  "diffModul")){
                             out$diffsGlobal[[i]] <- prop.tmp$diffsGlobal[[i]]
                             out$diffsGlobalLCC[[i]] <- prop.tmp$diffsGlobalLCC[[i]]
                           }

                           out$diffsGlobal[["diffnComp"]] <-
                             prop.tmp$diffsGlobal$diffnComp

                           out$diffsGlobalLCC[["difflccSize"]] <-
                             prop.tmp$diffsGlobalLC$difflccSize

                           out$diffsGlobalLCC[["difflccSizeRel"]] <-
                             prop.tmp$diffsGlobalLC$difflccSizeRel

                           out$absDiffsCentr <- prop.tmp$absDiffsCentr
                           
                           out$props <- prop.tmp$props
                           out$propsLCC <- prop.tmp$propsLCC

                           out
                         }

    if(verbose){
      # close progress bar
      close(pb)  
      if(cores > 1){
        # stop cluster
        message("Stopping socket cluster ... ", appendLF = FALSE)
      }
    }

    if(cores > 1) snow::stopCluster(cl)
    
    if(verbose){
      message("Done.")
    }

    
    if(verbose){
      message("Calculating p-values ... ", appendLF = FALSE)
    }

    #---------------------------------------------------------------------------
    # create matrices for differences of global properties

    permDiffsGlobal <- matrix(0, nrow = nPerm, ncol = 10)
    permDiffsGlobal_lcc <- matrix(0, nrow = nPerm, ncol = 11)
    
    names <- c("avDiss", "avPath", "Density", "VertConnect", "EdgeConnect", 
                  "NatConnect", "PEP", "ClustCoef", "Modul")
    names_diff <- paste0("diff", names)
    
    colnames(permDiffsGlobal) <- c("diffnComp", names_diff)
    
    colnames(permDiffsGlobal_lcc) <- c("difflccSize", "difflccSizeRel", 
                                       names_diff)
    
    for(i in 1:9){
      for(b in 1:nPerm){
        permDiffsGlobal[b,names_diff[i]] <- as.numeric(propsPerm[[b]]$diffsGlobal[names_diff[i]])
        permDiffsGlobal_lcc[b,names_diff[i]] <- as.numeric(propsPerm[[b]]$diffsGlobalLC[names_diff[i]])
      }
    }
    
    for(b in 1:nPerm){
      permDiffsGlobal[b, "diffnComp"] <- 
        as.numeric(propsPerm[[b]]$diffsGlobal["diffnComp"])
      
      permDiffsGlobal_lcc[b, "difflccSize"] <- 
        as.numeric(propsPerm[[b]]$diffsGlobalLCC["difflccSize"])
      
      permDiffsGlobal_lcc[b, "difflccSizeRel"] <- 
        as.numeric(propsPerm[[b]]$diffsGlobalLCC["difflccSizeRel"])
    }
    
    #---------------------------------------------------------------------------
    # create matrices for differences of centralities

    absDiffsPermDeg <- absDiffsPermBetw <- 
      absDiffsPermClose <- absDiffsPermEigen <- 
      permDeg1 <- permDeg2 <- permBetw1 <- permBetw2 <- 
      permClose1 <- permClose2 <- permEigen1 <- permEigen2 <-
      matrix(0, nrow = nPerm, ncol = ncol(adja1))

    for(b in 1:nPerm){
      absDiffsPermDeg[b,]   <- propsPerm[[b]]$absDiffsCentr$absDiffDeg
      absDiffsPermBetw[b,]  <- propsPerm[[b]]$absDiffsCentr$absDiffBetw
      absDiffsPermClose[b,] <- propsPerm[[b]]$absDiffsCentr$absDiffClose
      absDiffsPermEigen[b,] <- propsPerm[[b]]$absDiffsCentr$absDiffEigen
      
      permDeg1[b, ]  <- propsPerm[[b]]$props$deg1
      permDeg2[b, ]  <- propsPerm[[b]]$props$deg2
      permBetw1[b, ]  <- propsPerm[[b]]$props$betw1
      permBetw2[b, ]  <- propsPerm[[b]]$props$betw2
      permClose1[b, ] <- propsPerm[[b]]$props$close1
      permClose2[b, ] <- propsPerm[[b]]$props$close2
      permEigen1[b, ] <- propsPerm[[b]]$props$eigen1
      permEigen2[b, ] <- propsPerm[[b]]$props$eigen2
    }
    
    permDiffsCentr <- list(degree = absDiffsPermDeg,
                           between = absDiffsPermBetw,
                           close = absDiffsPermClose,
                           eigen = absDiffsPermEigen)
    
    permPropsCentr1 <- list(degree = permDeg1, 
                           between = permBetw1,
                           close = permClose1,
                           eigen = permEigen1)
    
    permPropsCentr2 <- list(degree = permDeg2, 
                            between = permBetw2,
                            close = permClose2,
                            eigen = permEigen2)
    
    #---------------------------------------------------------------------------
    # create matrices for global properties of the permuted data

    permPropsGlobal1 <- matrix(0, nrow = nPerm, ncol = 10)
    permPropsGlobal1_lcc <- matrix(0, nrow = nPerm, ncol = 11)
    
    permPropsGlobal2 <- matrix(0, nrow = nPerm, ncol = 10)
    permPropsGlobal2_lcc <- matrix(0, nrow = nPerm, ncol = 11)
    
    propnames <- c("avDiss", "avPath", "Density", "VertConnect", 
                   "EdgeConnect", "NatConnect", "PEP", "ClustCoef", 
                   "Modul")
    propnames2 <- c("avDiss", "avPath", "density", "vertConnect", 
                    "edgeConnect", "natConnect", "pep", "clustCoef", 
                    "modularity")
    
    colnames(permPropsGlobal1) <- c("nComp", propnames)
    colnames(permPropsGlobal2) <- c("nComp", propnames)
    
    colnames(permPropsGlobal1_lcc) <- c("lccSize", "lccSizeRel", propnames)
    colnames(permPropsGlobal2_lcc) <- c("lccSize", "lccSizeRel", propnames)
 
    for(i in 1:9){
      for(b in 1:nPerm){
        permPropsGlobal1[b,i] <- 
          as.numeric(propsPerm[[b]]$props[paste0(propnames2[i], 1)])
        
        permPropsGlobal2[b,i] <- 
          as.numeric(propsPerm[[b]]$props[paste0(propnames2[i], 2)])
        
        permPropsGlobal1_lcc[b,i] <- 
          as.numeric(propsPerm[[b]]$propsLCC[paste0(propnames2[i], 1)])
        
        permPropsGlobal2_lcc[b,i] <- 
          as.numeric(propsPerm[[b]]$propsLCC[paste0(propnames2[i], 2)])
      }
    }
    

    for(b in 1:nPerm){
      permPropsGlobal1[b, "nComp"] <- 
        as.numeric(propsPerm[[b]]$props["nComp1"])
      
      permPropsGlobal2[b, "nComp"] <- 
        as.numeric(propsPerm[[b]]$props["nComp2"])
      
      permPropsGlobal1_lcc[b, "lccSize"] <- 
        as.numeric(propsPerm[[b]]$propsLCC["lccSize1"])
      
      permPropsGlobal2_lcc[b, "lccSize"] <- 
        as.numeric(propsPerm[[b]]$propsLCC["lccSize2"])
      
      permPropsGlobal1_lcc[b, "lccSizeRel"] <- 
        as.numeric(propsPerm[[b]]$propsLCC["lccSizeRel1"])
      
      permPropsGlobal2_lcc[b, "lccSizeRel"] <- 
        as.numeric(propsPerm[[b]]$propsLCC["lccSizeRel2"])
    }
    
    
    #---------------------------------------------------------------------------
    # compute p-values
    
    pvalnames <- paste0("pval", names)
    
    for(i in seq_along(pvalnames)){
      assign(pvalnames[i], 
             (sum(permDiffsGlobal[, names_diff[i]] >= 
                    props$diffsGlobal[[names_diff[i]]]) + 1) / (nPerm + 1))
      
      assign(paste0(pvalnames[i], "_lcc"), 
             (sum(permDiffsGlobal_lcc[, names_diff[i]] >= 
                    props$diffsGlobalLCC[[names_diff[i]]]) + 1) / (nPerm + 1))
    }
    
    pvalnComp <- (sum(permDiffsGlobal[, "diffnComp"] >= 
                        props$diffsGlobal[["diffnComp"]]) + 1) / (nPerm + 1)
    
    pvallccSize <- (sum(permDiffsGlobal_lcc[, "difflccSize"] >= 
                          props$diffsGlobalLCC[["difflccSize"]]) + 1) / (nPerm + 1)
    
    pvallccSizeRel <- (sum(permDiffsGlobal_lcc[, "difflccSizeRel"] >= 
                             props$diffsGlobalLCC[["difflccSizeRel"]]) + 1) / (nPerm + 1)

    pvalDiffDeg <- sapply(1:ncol(adja1), function(i){
      (sum(absDiffsPermDeg[, i] >= props$absDiffs$absDiffDeg[i]) + 1) / (nPerm + 1)
    })
    pvalDiffBetw <- sapply(1:ncol(adja1), function(i){
      (sum(absDiffsPermBetw[, i] >= props$absDiffs$absDiffBetw[i]) + 1) / (nPerm + 1)
    })
    pvalDiffClose <- sapply(1:ncol(adja1), function(i){
      (sum(absDiffsPermClose[, i] >= props$absDiffs$absDiffClose[i]) + 1) / (nPerm + 1)
    })
    pvalDiffEigen <- sapply(1:ncol(adja1), function(i){
      (sum(absDiffsPermEigen[, i] >= props$absDiffs$absDiffEigen[i]) + 1) / (nPerm + 1)
    })

    names(pvalDiffDeg) <- colnames(adja1)
    names(pvalDiffBetw) <- colnames(adja1)
    names(pvalDiffClose) <- colnames(adja1)
    names(pvalDiffEigen) <- colnames(adja1)

    if(verbose) message("Done.")

    if(verbose & adjustPerm != "none"){
      message("Adjust for multiple testing using '", adjustPerm, "' ... ",
              appendLF = FALSE)
    }

    # adjust for multiple testing
    
    if(adjust == "adaptBH" && !requireNamespace("limma", quietly = TRUE)){
      
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
    
    pAdjustDiffDeg <- multAdjust(pvals = pvalDiffDeg, adjust = adjustPerm,
                                 trueNullMethod = trueNullMethod,
                                 verbose = FALSE)
    
    pAdjustDiffBetw <- multAdjust(pvals = pvalDiffBetw, adjust = adjustPerm,
                                  trueNullMethod = trueNullMethod,
                                  verbose = FALSE)
    
    pAdjustDiffClose <- multAdjust(pvals = pvalDiffClose, adjust = adjustPerm,
                                   trueNullMethod = trueNullMethod,
                                   verbose = FALSE)
    
    pAdjustDiffEigen <- multAdjust(pvals = pvalDiffEigen, adjust = adjustPerm,
                                   trueNullMethod = trueNullMethod,
                                   verbose = FALSE)
    
    if(verbose & adjustPerm != "none") message("Done.")
    
    output$pvalDiffGlobal <- list(pvalnComp = pvalnComp,
                                  pvalavDiss = pvalavDiss, 
                                  pvalavPath = pvalavPath,
                                  pvalDensity = pvalDensity,
                                  pvalEdgeConnect = pvalEdgeConnect,
                                  pvalVertConnect = pvalVertConnect,
                                  pvalNatConnect = pvalNatConnect,
                                  pvalPEP = pvalPEP,
                                  pvalClustCoef = pvalClustCoef,
                                  pvalModul = pvalModul)
    
    output$pvalDiffGlobalLCC <- list(pvallccSize = pvallccSize,
                                     pvallccSizeRel = pvallccSizeRel,
                                     pvalavDiss = pvalavDiss_lcc, 
                                     pvalavPath = pvalavPath_lcc,
                                     pvalDensity = pvalDensity_lcc,
                                     pvalEdgeConnect = pvalEdgeConnect_lcc,
                                     pvalVertConnect = pvalVertConnect_lcc,
                                     pvalNatConnect = pvalNatConnect_lcc,
                                     pvalPEP = pvalPEP_lcc,
                                     pvalClustCoef = pvalClustCoef_lcc,
                                     pvalModul = pvalModul_lcc)
    
    output$pvalDiffCentr <- list(pvalDiffDeg = pvalDiffDeg,
                                 pvalDiffBetw = pvalDiffBetw,
                                 pvalDiffClose = pvalDiffClose,
                                 pvalDiffEigen = pvalDiffEigen)
    
    output$pvalDiffCentrAdjust <- list(pAdjustDiffDeg = pAdjustDiffDeg,
                                       pAdjustDiffBetw = pAdjustDiffBetw,
                                       pAdjustDiffClose = pAdjustDiffClose,
                                       pAdjustDiffEigen = pAdjustDiffEigen)
    
    if(returnPermProps){
      output$permPropsGlobal1 <- permPropsGlobal1
      output$permPropsGlobal2 <- permPropsGlobal2
      output$permPropsGlobalLCC1 <- permPropsGlobal1_lcc
      output$permPropsGlobalLCC2 <- permPropsGlobal2_lcc
      output$permDiffGlobal <- permDiffsGlobal
      output$permDiffGlobalLCC <- permDiffsGlobal_lcc
    }
    
    if(returnPermCentr){
      output$permPropsCentr1 <- permPropsCentr1
      output$permPropsCentr2 <- permPropsCentr2
      output$permDiffCentr <- permDiffsCentr
    }
  } 
  
  class(output) <- "microNetComp"
  return(output)
}




