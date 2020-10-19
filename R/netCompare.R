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
#' @param storeAssoPerm logical indicating whether the association (or
#'   dissimilarity) matrices for the permuted data should be returned. If
#'   \code{TRUE} (default), \code{assoPerm} (or \code{dissPerm}) as part of the
#'   output can be passed to \code{netCompare} again to save runtime.
#' @param assoPerm a list with two elements used for the permutation procedure.
#'   Each entry must contain association matrices for \code{"nPerm"}
#'   permutations. This can be the \code{"assoPerm"} value as part of the
#'   output either returned from \code{diffnet} or from
#'   \code{\link{netCompare}}. See the example.
#' @param dissPerm usage analog to \code{assoPerm} if a dissimilarity measure
#'   has been used for network construction.
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
#'   data set. The p-values are calcutated as the proportion of
#'   "permutation-differences" being larger than the observed difference. A
#'   pseudocount is added to the numerator and denominator in order to avoid
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
#' @return Returned is an object of class \code{microNetComp}. Important
#'   elements are the following:\cr
#'   \strong{With and without permutation tests:}
#'   \tabular{ll}{
#'   \code{jaccDeg,jaccBetw,jaccClose,jaccEigen}\tab values of Jaccard's index
#'   for the centrality measures\cr
#'   \code{jaccHub}\tab Jaccard index for the sets of hub nodes\cr
#'   \code{randInd}\tab calculated Rand index\cr
#'   \code{properties}\tab a list with calculated network properties for both
#'   groups\cr
#'   \code{diffs}\tab vectors with differences between the centrality values\cr
#'   \code{countMatrices}\tab the two count matrices returned
#'   from \code{netConstruct}\cr
#'   \code{assoMatrices}\tab the two association matrices returned
#'   from \code{netConstruct}\cr
#'   \code{dissMatrices}\tab the two dissimilarity matrices returned
#'   from \code{netConstruct}\cr
#'   \code{adjaMatrices}\tab the two adjacency matrices returned
#'   from \code{netConstruct}\cr
#'   \code{groups}\tab the group names returned from \code{netConstruct}
#'   }
#'   \strong{Additional output if permutation tests are conducted:}
#'   \tabular{ll}{
#'   \code{avPath}\tab difference between average path lengths with
#'   corresponding p-value\cr
#'   \code{clustCoef}\tab difference between the global clustering coefficients
#'   with corresponding p-value\cr
#'   \code{modul}\tab difference between modularity values for the
#'   determined clusterings with corresponding p-value\cr
#'   \code{vertConnect}\tab difference between the vertex connectivity values
#'   with corresponding p-value\cr
#'   \code{edgeConnect}\tab difference between the edge connectivity values
#'   with corresponding p-value\cr
#'   \code{density}\tab difference between density values with corresponding
#'   p-value\cr
#'   \code{pvalDiffCentr}\tab p-values of the tests for differential centrality
#'   values\cr
#'   \code{pvalDiffCentrAdjust}\tab adjusted p-values of the tests for
#'   differential centrality values\
#'   }
#'   \strong{Additional output if NO permutation tests are conducted:}
#'   \tabular{ll}{
#'   \code{diffPath}\tab difference between average path lengths\cr
#'   \code{diffClust}\tab difference between the global clustering
#'   coefficients\cr
#'   \code{diffModul}\tab difference between modularity values for the
#'   determined clusterings\cr
#'   \code{diffVertConnect}\tab difference between the vertex connectivity
#'   values\cr
#'   \code{diffEdgeConnect}\tab difference between the edge connectivity
#'   values\cr
#'   \code{diffDensity}\tab difference between the density values
#'   }
#' @examples
#' # load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#'
#' # generate a random group vector
#' set.seed(123456)
#' group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)
#'
#' # network construction:
#' amgut_net <- netConstruct(amgut1.filt, group = group,
#'                           measure = "pearson",
#'                           filtTax = "highestVar",
#'                           filtTaxPar = list(highestVar = 30),
#'                           zeroMethod = "pseudo", normMethod = "clr")
#'
#' # network analysis:
#' amgut_props <- netAnalyze(amgut_net, clustMethod = "cluster_fast_greedy",
#'                           hubPar = "eigenvector")
#'
#' # network plot:
#' plot(amgut_props, sameLayout = TRUE)
#'
#' # network comparison:
#' # without permutation tests:
#' amgut_comp1 <- netCompare(amgut_props, permTest = FALSE)
#' summary(amgut_comp1)
#' summary(amgut_comp1, showCentr = "degree", numbNodes = 20)
#'
#' \donttest{
#' # with permutation tests (with only 100 permutations to decrease runtime):
#' amgut_comp2 <- netCompare(amgut_props, permTest = TRUE, nPerm = 100L, cores = 4L)
#'
#' # the estimated association matrices for the permuted data can be passed to
#' # netCompare in order to reduce execution time:
#' amgut_comp3 <- netCompare(amgut_props, permTest = TRUE, nPerm = 100L,
#'                           assoPerm = amgut_comp2$assoPerm)
#'
#' summary(amgut_comp2)
#' summary(amgut_comp3)
#' }
#'
#' @seealso \code{\link{summary.microNetComp}}, \code{\link{netConstruct}},
#'   \code{\link{netAnalyze}}
#' @references \insertRef{benjamini2000adaptive}{NetCoMi} \cr
#'   \insertRef{farcomeni2007some}{NetCoMi} \cr
#'   \insertRef{gill2010statistical}{NetCoMi} \cr
#'   \insertRef{qannari2014significance}{NetCoMi} \cr
#'   \insertRef{real1996probabilistic}{NetCoMi}
#' @import foreach doSNOW
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
                       seed = NULL, verbose = TRUE,
                       storeAssoPerm = TRUE,
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

  weighted <- clustMethod <- clustPar <- hubPar <- hubQuant <- weightDeg <- NULL
  normDeg <- normBetw <- normClose <- normEigen <- NULL
  softThreshPower <- softThreshCut <- softThreshType <- NULL
  kNeighbor <- knnMutual <- NULL
  measure <- measurePar <- scaleDiss <- sparsMethod <- thresh <- alpha <- NULL
  lfdrThresh <- nboot <-  dissFunc <- dissFuncPar <- simFunc <- simFuncPar <- NULL

  for(i in 1:length(x$paramsProperties)){
    assign(names(x$paramsProperties)[i], x$paramsProperties[[i]])
  }

  parnames <- c("zeroMethod", "zeroParam", "normMethod", "normParam",
                "measure", "measurePar", "sparsMethod", "thresh", "adjust",
                "alpha", "lfdrThresh", "nboot", "softThreshType",
                "softThreshCut", "softThreshPower", "kNeighbor", "knnMutual",
                "dissFunc", "dissFuncPar",
                "simFunc", "simFuncPar", "scaleDiss", "weighted")

  for(i in 1:length(parnames)){
    assign(parnames[i], x$paramsNetConstruct[[parnames[i]]])
  }

  assoType <- x$input$assoType
  distNet <- ifelse(assoType == "dissimilarity", TRUE, FALSE)

  xgroups <- x$input$groups
  sampleSize <- x$input$sampleSize
  twoNets <- x$input$twoNets

  count1 <- x$input$countMat1
  count2 <- x$input$countMat2

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


  if(!is.null(seed)) set.seed(seed)


  #---------------------------------------------------------------------------
  # calculate and compare network properties
  if(permTest & verbose) message("Calculate network properties ... ",
                                 appendLF = FALSE)

  props <- calc_diff_props(adja1 = adja1, adja2 = adja2,
                         dissMat1 = dissMat1, dissMat2 = dissMat2,
                         weighted = weighted,
                         clustMethod = clustMethod, clustPar = clustPar,
                         hubPar = hubPar,  hubQuant = hubQuant,
                         jaccQuant = jaccQuant, lnormFit = lnormFit,
                         connect = connect, weightDeg = weightDeg,
                         normDeg = normDeg,  normBetw = normBetw,
                         normClose = normClose, normEigen = normEigen,
                         nPermRand = nPermRand)
  if(permTest & verbose) message("Done.")


  #---------------------------------------------------------------------------
  # generate teststatistics for permutated data

  if(permTest){

    if(verbose) message("Execute permutation tests ... ")

    if(!is.null(seed)){
      seeds <- sample.int(1e8, size = nPerm)
    } else{
      seeds <- NULL
    }

    if(assoType == "dissimilarity"){
      n1 <- ncol(count1)
      n2 <- ncol(count2)
      n <- n1 + n2
      xbind <- cbind(count1, count2)
    } else{
      n1 <- nrow(count1)
      n2 <- nrow(count2)
      n <- n1 + n2
      xbind <- rbind(count1, count2)
    }

    if(cores > 1){
      if(parallel::detectCores() < cores) cores <- parallel::detectCores()

      cl <- makeCluster(cores, outfile = "")
      registerDoSNOW(cl)
      '%do_or_dopar%' <- get('%dopar%')

    } else{
      '%do_or_dopar%' <- get('%do%')
    }

    if(verbose){
      pb<-txtProgressBar(0, nPerm, style=3)
      
      progress<-function(n){
        setTxtProgressBar(pb,n)
      }
      
      opts <- list(progress=progress)
    } else{
      opts <- list()
    }

    if(!is.null(logFile)) cat("", file=logFile, append=FALSE)

    p <- NULL
    propsPerm <- foreach(p = 1:nPerm,
                         .packages = c("igraph", "WGCNA", "SPRING",
                                                    "fdrtool",
                                                    "SpiecEasi",
                                                    "ccrepe",
                                                    "vegan",
                                                    "LaplacesDemon",
                                                    "robCompositions",
                                                    "propr",
                                                    "zCompositions",
                                                    "DESeq2",
                                                    "compositions"),
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
                             assoMat1.tmp <- assoPerm[[1]][[p]]
                             assoMat2.tmp <- assoPerm[[2]][[p]]
                             count1.tmp <- count2.tmp <- NULL
                           } else if(!is.null(dissPerm)){
                             assoMat1.tmp <- dissPerm[[1]][[p]]
                             assoMat2.tmp <- dissPerm[[2]][[p]]
                             count2.tmp <- count2.tmp <- NULL
                           } else{

                             index <- sample(1:n, n)
                             if(assoType == "dissimilarity"){
                               count1.tmp <- xbind[ ,index[1:n1]]
                               count2.tmp <- xbind[ ,index[(n1+1):n]]
                             } else{
                               count1.tmp <- xbind[index[1:n1], ]
                               count2.tmp <- xbind[index[(n1+1):n], ]
                             }

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

                             assoMat1.tmp <- calc_association(count1.tmp,
                                                              measure = measure,
                                                              measurePar = measurePar,
                                                              verbose = FALSE)
                             assoMat2.tmp <- calc_association(count2.tmp,
                                                              measure = measure,
                                                              measurePar = measurePar,
                                                              verbose = FALSE)
                           }

                           if(distNet){
                             dissEst1.tmp <- assoMat1.tmp
                             dissEst2.tmp <- assoMat2.tmp
                             if(scaleDiss){
                               assoMat1.tmp <- scale_diss(dissEst1.tmp)
                               assoMat2.tmp <- scale_diss(dissEst2.tmp)
                             }
                             assoEst1.tmp <- assoEst2.tmp <- NULL
                           } else{
                             assoEst1.tmp <- assoMat1.tmp
                             assoEst2.tmp <- assoMat2.tmp
                             dissEst1.tmp <- dissEst2.tmp <- NULL

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

                           prop.tmp <- calc_diff_props(adja1 = adja1.tmp,
                                                       adja2 = adja2.tmp,
                                                       dissMat1 = dissMat1.tmp,
                                                       dissMat2 = dissMat2.tmp,
                                                       weighted = weighted,
                                                       clustMethod = clustMethod,
                                                       clustPar = clustPar,
                                                       hubPar = hubPar,
                                                       hubQuant = hubQuant,
                                                       jaccQuant = jaccQuant,
                                                       lnormFit = lnormFit,
                                                       connect = connect,
                                                       weightDeg = weightDeg,
                                                       normDeg = normDeg,
                                                       normBetw = normBetw,
                                                       normClose = normClose,
                                                       normEigen = normEigen,
                                                       nPermRand = nPermRand,
                                                       testJacc = FALSE,
                                                       testRand = FALSE)

                           out <- list()
                           for(i in c("diffPath", "diffClust",
                                      "diffModul", "diffVertConnect",
                                      "diffEdgeConnect",
                                      "diffDensity", "absDiffs")){
                             out[[i]] <- prop.tmp[[i]]
                           }

                           if(storeAssoPerm){
                             if(distNet){
                               out$dissEst1 <- dissEst1.tmp
                               out$dissEst2 <- dissEst2.tmp
                             } else{
                               out$assoEst1 <- assoEst1.tmp
                               out$assoEst2 <- assoEst2.tmp
                             }
                           }
                           out
                         }

    if(verbose){
      close(pb)
      message("Stopping socket cluster ... ", appendLF = FALSE)
    }

    if(cores > 1) snow::stopCluster(cl)
    if(verbose){
      message("Done.")
      message("Calculating p-values ... ", appendLF = FALSE)
    }

    results <- matrix(0, nrow = nPerm, ncol = 6)
    selnames <- c("diffPath", "diffClust", "diffModul", "diffVertConnect",
                  "diffEdgeConnect", "diffDensity")

    for(i in 1:6){
      for(b in 1:nPerm){
        results[b,i] <- as.numeric(propsPerm[[b]][selnames[i]])
      }
    }

    absDiffsPermDeg <- matrix(0, nrow = nPerm, ncol = ncol(adja1))
    absDiffsPermBetw <- matrix(0, nrow = nPerm, ncol = ncol(adja1))
    absDiffsPermClose <- matrix(0, nrow = nPerm, ncol = ncol(adja1))
    absDiffsPermEigen <- matrix(0, nrow = nPerm, ncol = ncol(adja1))

    for(i in 1:nPerm){
      absDiffsPermDeg[i,] <- propsPerm[[i]]$absDiffs$absDiffDeg
      absDiffsPermBetw[i,] <- propsPerm[[i]]$absDiffs$absDiffBetw
      absDiffsPermClose[i,] <- propsPerm[[i]]$absDiffs$absDiffClose
      absDiffsPermEigen[i,] <- propsPerm[[i]]$absDiffs$absDiffEigen
    }


    if(storeAssoPerm){
      if(distNet){
        dissEstPerm1 <- dissEstPerm2 <- list()
        for(i in 1:nPerm){
          dissEstPerm1[[i]] <- propsPerm[[i]]$dissEst1
          dissEstPerm2[[i]] <- propsPerm[[i]]$dissEst2
        }
        assoEstPerm1 <- assoEstPerm2 <- NULL
      } else{
        assoEstPerm1 <- assoEstPerm2 <- list()
        for(i in 1:nPerm){
          assoEstPerm1[[i]] <- propsPerm[[i]]$assoEst1
          assoEstPerm2[[i]] <- propsPerm[[i]]$assoEst2
        }
        dissEstPerm1 <- dissEstPerm2 <- NULL
      }
    } else{
      assoEstPerm1 <- assoEstPerm2 <- NULL
      dissEstPerm1 <- dissEstPerm2 <- NULL
    }

    pvalPath <- (sum(results[,1] >= props$diffPath) + 1) / (nPerm + 1)
    pvalClust <- (sum(results[,2] >= props$diffClust) + 1) / (nPerm + 1)
    pvalModul <- (sum(results[,3] >= props$diffModul) + 1) / (nPerm + 1)
    pvalVertConnect <- (sum(results[,4] >= props$diffVertConnect) + 1) / (nPerm + 1)
    pvalEdgeConnect <- (sum(results[,5] >= props$diffEdgeConnect) + 1) / (nPerm + 1)
    pvalDensity <- (sum(results[,6] >= props$diffDensity) + 1) / (nPerm + 1)


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

    output <- list(jaccDeg = props$jaccDeg,
                   jaccBetw = props$jaccBetw,
                   jaccClose = props$jaccClose,
                   jaccEigen = props$jaccEigen,
                   jaccHub = props$jaccHub,
                   avPath = c(diff = props$diffPath, pval = pvalPath),
                   clustCoef = c(diff = props$diffClust, pval = pvalClust),
                   modul = c(diff = props$diffModul, pval = pvalModul),
                   vertConnect = c(diff = props$diffVertConnect,
                                   pval = pvalVertConnect),
                   edgeConnect = c(diff = props$diffEdgeConnect,
                                   pval = pvalEdgeConnect),
                   density = c(diff = props$diffDensity,
                               pval = pvalDensity),
                   randInd = props$randInd,
                   properties = props$props,
                   diffs = props$diffs,
                   pvalDiffCentr = list(pvalDiffDeg = pvalDiffDeg,
                                    pvalDiffBetw = pvalDiffBetw,
                                    pvalDiffClose = pvalDiffClose,
                                    pvalDiffEigen = pvalDiffEigen),
                   pvalDiffCentrAdjust = list(pAdjustDiffDeg = pAdjustDiffDeg,
                                          pAdjustDiffBetw = pAdjustDiffBetw,
                                          pAdjustDiffClose = pAdjustDiffClose,
                                          pAdjustDiffEigen = pAdjustDiffEigen),
                   countMatrices = list(count1 = count1,
                                        count2 = count2),
                   assoMatrices = list(assoMat1 = assoMat1,
                                       assoMat2 = assoMat2),
                   dissMatrices = list(dissMat1 = dissMat1,
                                       dissMat2 = dissMat2),
                   adjaMatrices = list(adja1 = adja1,
                                       adja2 = adja2),
                   assoPerm = list(assoEstPerm1 = assoEstPerm1,
                                   assoEstPerm2 = assoEstPerm2),
                   dissPerm = list(dissEstPerm1 = dissEstPerm1,
                                   dissEstPerm2 = dissEstPerm2),
                   groups = list(group1 = xgroups[1], group2 = xgroups[2]),
                   call = match.call())
  } else{
    output <- list(jaccDeg = props$jaccDeg,
                   jaccBetw = props$jaccBetw,
                   jaccClose = props$jaccClose,
                   jaccEigen = props$jaccEigen,
                   jaccHub = props$jaccHub,
                   diffPath = props$diffPath,
                   diffClust = props$diffClust,
                   diffModul = props$diffModul,
                   diffVertConnect = props$diffVertConnect,
                   diffEdgeConnect = props$diffEdgeConnect,
                   diffDensity = props$diffDensity,
                   randInd = props$randInd,
                   properties = props$props,
                   diffs = props$diffs,
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
                   call = match.call())
  }
  class(output) <- "microNetComp"
  return(output)
}




