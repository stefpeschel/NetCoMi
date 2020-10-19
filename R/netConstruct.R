#' @title Constructing Networks for Microbiome Data
#'
#' @description Constructing microbial association networks and dissimilarity
#'   based networks (where nodes are subjects) from compositional count data.
#'
#' @details The function enables the construction of either a single network or
#'   two networks that can be compared using the function
#'   \code{\link{netCompare}}. Four types of association  measures are
#'   available: correlation, conditional dependence, proportionality, and
#'   dissimilarity.  Depending on the measure, nodes are either taxa or subjects
#'   in the resulting network: In association-based networks (correlation,
#'   partial correlation, conditional dependence, proportionality) nodes are
#'   taxa, whereas in dissimilarity-based networks nodes are subjects.\cr\cr
#'   In order to conduct a network comparison, two networks can be  constructed
#'   by: \tabular{ll}{
#'   (1)\tab passing a group vector to \code{group} (of length
#'   \code{nrow(data)} for association networks and of length \code{ncol(data)}
#'   for dissimilarity-based networks), or\cr
#'   (2)\tab passing a second count matrix to \code{data2} (column names must
#'   match for constructing association networks and row names must match for
#'   dissimilarity-based networks), or \cr
#'   (3)\tab passing a second association/dissimilarity matrix to \code{data2}.}
#'   \cr
#'   The object returned from \code{netConstruct} can either be passed to
#'   \code{\link{netAnalyze}} for network analysis, or to
#'   \code{\link{diffnet}} to construct a differential network from the
#'   estimated associations.\cr\cr
#'   \strong{Association measures:}
#'   \tabular{ll}{
#'   Argument \tab Function\cr
#'   \code{"pearson"}\tab \code{\link[stats]{cor}} \cr
#'   \code{"spearman"}\tab \code{\link[stats]{cor}} \cr
#'   \code{"bicor"}\tab \code{\link[WGCNA]{bicor}} \cr
#'   \code{"sparcc"}\tab \code{\link[SpiecEasi]{sparcc}} \cr
#'   \code{"cclasso"}\tab \code{\link[NetCoMi]{cclasso}} \cr
#'   \code{"ccrepe"}\tab \code{\link[ccrepe]{ccrepe}} \cr
#'   \code{"spieceasi"}\tab \code{\link[SpiecEasi]{spiec.easi}} \cr
#'   \code{"spring"}\tab \code{\link[SPRING]{SPRING}} \cr
#'   \code{"gcoda"}\tab \code{\link[NetCoMi]{gcoda}} \cr
#'   \code{"propr"}\tab \code{\link[propr]{propr}}}
#'   \cr
#'   \strong{Dissimilarity measures:}
#'   \tabular{lll}{
#'   Argument \tab Function \tab Measure\cr
#'   \code{"euclidean"} \tab \code{\link[vegan]{vegdist}} \tab
#'   Euclidean distance \cr
#'   \code{"bray"}\tab \code{\link[vegan]{vegdist}} \tab
#'   Bray-Curtis  dissimilarity \cr
#'   \code{"kld"}\tab \code{\link[LaplacesDemon]{KLD}} \tab
#'   Kullback-Leibler divergence \cr
#'   \code{"jeffrey"}\tab \code{\link[LaplacesDemon]{KLD}} \tab
#'   Jeffrey divergence\cr
#'   \code{"jsd"}\tab \code{\link[LaplacesDemon]{KLD}} \tab
#'   Jensen-Shannon divergence \cr
#'   \code{"ckld"}\tab \code{\link[base]{log}} \tab
#'   Compositional Kullback-Leibler divergence \cr
#'   \code{"aitchison"}\tab \code{\link[vegan]{vegdist}},
#'   \code{\link[robCompositions]{cenLR}}\tab
#'   Aitchison distance}
#'   Definitions:\cr
#'   \describe{
#'   \item{Kullback-Leibler divergence:}{Since KLD is not symmetric,
#'   0.5 * (KLD(p(x)||p(y)) + KLD(p(y)||p(x))) is returned.}
#'   \item{Jeffrey divergence:}{Jeff = KLD(p(x)||p(y)) + KLD(p(y)||p(x))}
#'   \item{Jensen-Shannon divergence:}{JSD = 0.5 KLD(P||M) + 0.5 KLD(Q||M),
#'   where P=p(x), Q=p(y), and M=0.5(P+Q).}
#'   \item{Compositional Kullback-Leibler divergence:}{cKLD(x,y) =
#'   p/2 * log(A(x/y) * A(y/x)), where A(x/y) is the arithmetic mean of the
#'   vector of ratios x/y.}
#'   \item{Aitchison distance:}{Euclidean distance of the clr-transformed data.}
#'   }\cr
#'   \strong{Normalization methods:}
#'   \tabular{lll}{
#'   Argument \tab Method \tab Function\cr
#'   \code{"TSS"} \tab Total sum scaling \tab t(apply(countMat, 1, function(x)
#'   x/sum(x)))\cr
#'   \code{"CSS"} \tab Cumulative sum scaling \tab
#'   \code{\link[metagenomeSeq]{cumNormMat}}\cr
#'   \code{"COM"} \tab Common sum scaling \tab
#'   t(apply(countMat, 1, function(x) x * min(rowSums(countMat)) / sum(x)))\cr
#'   \code{"rarefy"} \tab Rarefying \tab \code{\link[vegan:rarefy]{rrarefy}}\cr
#'   \code{"VST"} \tab Variance stabilizing transformation\tab
#'   \code{\link[DESeq2]{varianceStabilizingTransformation}}\cr
#'   \code{"clr"} \tab Centered log-ratio transformation\tab
#'   \code{\link[robCompositions]{cenLR}}
#'   }
#'   These methods (except rarefying) are described in
#'   \cite{Badri et al.(2018)}.
#'   \cr\cr
#'   \strong{Transforming associations to dissimilarities:}
#'   \tabular{ll}{
#'   Argument \tab Function\cr
#'   \code{"signed"} \tab sqrt(0.5 * (1 - x))\cr
#'   \code{"unsigned"} \tab sqrt(1 - x^2)\cr
#'   \code{"signedPos"} \tab  diss <- sqrt(0.5 * (1-x))\cr
#'     \tab diss[x < 0] <- 0\cr
#'   \code{"TOMdiss"} \tab \code{\link[WGCNA:TOMsimilarity]{TOMdist}}
#'   }
#'
#' @param data numeric matrix. Can be a count matrix (rows are samples, columns
#'   are OTUs/taxa), a phyloseq object, or an association/dissimilarity matrix
#'   (\code{dataType} must be set).
#' @param data2 optional numeric matrix corresponding to group 2. Can be either
#'   a second count matrix/phyloseq object or a second association/dissimilarity
#'   matrix.
#' @param dataType character indicating the data type. Defaults to "counts",
#'   which means that \code{data} (and data2) is a count matrix or object of
#'   class \code{\link[phyloseq:phyloseq-class]{phyloseq}}. Further options
#'   are "correlation", "partialCorr" (partial correlation), "condDependence"
#'   (conditional dependence), "proportionality" and "dissimilarity".
#' @param group optional binary vector used or splitting the data into two
#'   groups. If \code{group} is \code{NULL} (default), a single network is
#'   constructed. See details.
#' @param measure character specifying the method used for either computing the
#'   associations between taxa or dissimilarities between subjects.
#'   Ignored if \code{data} is not a count matrix (if \code{dataType} is not set
#'   to \code{"counts"}). Available measures are:
#'   \code{"pearson"}, \code{"spearman"}, \code{"bicor"}, \code{"sparcc"},
#'   \code{"cclasso"}, \code{"ccrepe"}, \code{"spieceasi"} (default),
#'   \code{"spring"}, \code{"gcoda"} and \code{"propr"} as association measures,
#'   and \code{"euclidean"}, \code{"bray"}, \code{"kld"}, \code{"jeffrey"},
#'   \code{"jsd"}, \code{"ckld"}, and \code{"aitchison"} as dissimilarity
#'   measures. Parameters are set via \code{measurePar}.
#' @param measurePar list with parameters passed to the function for computing
#'   associations/dissimilarities. See details for the respective functions.
#' @param filtTax character indicating how taxa shall be filtered. Possible
#'   options are:
#'   \describe{
#'   \item{\code{"none"}}{Default. All taxa are kept.}
#'   \item{\code{"totalReads"}}{Keep taxa with a total number
#'   of reads of at least x.}
#'   \item{\code{"relFreq"}}{Keep taxa whose number of reads is at
#'   least x\% of the total number of reads.}
#'   \item{\code{"numbSamp"}}{Keep taxa observed in at least x samples.}
#'   \item{\code{"highestVar"}}{Keep the x taxa with highest variance.}
#'   \item{\code{"highestFreq"}}{Keep the x taxa with highest frequency.}}
#'   Except for "highestVar" and "highestFreq", different filter methods can be
#'   combined. The values x are set via \code{filtTaxPar}.
#' @param filtTaxPar list with parameters for the filter methods given by
#'   \code{filtTax}. Possible list entries are: \code{"totalReads"} (int),
#'   \code{"relFreq"} (value in [0,1]), \code{"numbSamp"} (int),
#'   \code{"highestVar"} (int), \code{"highestFreq"} (int).
#' @param filtSamp character indicating how samples shall be filtered. Possible
#'   options are: \describe{
#'   \item{\code{"none"}}{Default. All samples are kept.}
#'   \item{\code{"totalReads"}}{Keep samples with a total number of reads of at
#'   least x.}
#'   \item{\code{"numbTaxa"}}{Keep samples for which at least
#'   x taxa are observed.}
#'   \item{\code{"highestFreq"}}{Keep the x samples with highest frequency.}}
#'   Except for "highestFreq", different filter methods can be
#'   combined. The values x are set via \code{filtSampPar}.
#' @param filtSampPar list with parameters for the filter methods given by
#'   \code{filtSamp}. Possible list entries are: \code{"totalReads"} (int),
#'   \code{"numbTaxa"} (int), \code{"highestFreq"} (int).
#' @param zeroMethod character indicating the method used for zero replacement.
#'   Possible values are: \code{"none"} (default), \code{"pseudo"}
#'   (a unit pseudo count is added to the original count data),
#'   \code{"multRepl"} (multiplicative replacement using
#'   \code{\link[zCompositions]{multRepl}}), \code{"alrEM"} (modified EM
#'   alr-algorithm using \code{\link[zCompositions]{lrEM}}), \code{"bayesMult"}
#'   (Bayesian-multiplicative replacement using
#'   \code{\link[zCompositions]{cmultRepl}}). The corresponding parameters are
#'   set via \code{zeroPar}. \code{zeroMethod} is ignored if the approach for
#'   calculating the associations/dissimilarity includes zero handling.
#'   Defaults to \code{"multRepl"} or \code{"pseudo"} (depending on the expected
#'   input of the normalization function and measure) if zero replacement is
#'   required.
#' @param zeroPar list with parameters passed to the function for zero
#'   replacement (\code{zeroMethod}). See the help page of the respective
#'   function for details.
#' @param normMethod character indicating the normalization method (in order to
#'   make counts of different samples comparable). Possible options are:
#'   \code{"none"} (default), \code{"TSS"} (or \code{"fractions"}), \code{"CSS"},
#'   \code{"COM"}, \code{"rarefy"}, \code{"VST"}, \code{"clr"}.
#'   The corresponding parameters are set via \code{normPar}.
#' @param normPar list with parameters passed to the function for normalization
#'   (defined by \code{normMethod}).
#' @param sparsMethod character indicating the method used for sparsification
#'   (selected edges that are connected in the network). Available methods are:
#'   \describe{
#'   \item{\code{"none"}}{Leads to a fully connected network}
#'   \item{\code{"t-test"}}{Student's t-test. Significance level and multiple
#'   testing adjustment is specified via \code{alpha} and \code{adjust}.}
#'   \item{\code{"bootstrap"}}{Bootstrap procedure as described in
#'   \cite{Friedman and Alm (2012)}. Corresponding arguments are
#'   \code{nboot}, \code{cores}, and \code{logFile}.}
#'   \item{\code{"threshold"}}{Default. Selected are taxa
#'   pairs with an absolute association/dissimilarity greater than or equal to
#'   the threshold defined via \code{thresh}.}
#'   \item{\code{"softThreshold"}}{Soft thresholding method according to
#'   \cite{Zhang and Horvath (2005)} available in the
#'   \code{\link[WGCNA:pickSoftThreshold]{WGCNA}} package. Corresponding
#'   arguments are \code{softThreshType}, \code{softThreshPower}, and
#'   \code{softThreshCut}.}
#'   \item{\code{"knn"}}{Construct a k-nearest neighbor or mutual k-nearest
#'   neighbor graph using \code{\link[cccd]{nng}}. Corresponding
#'   arguments are \code{kNeighbor}, and \code{knnMutual}.}}
#' @param thresh numeric value defining the threshold if \code{sparsMethod} is
#'   set to \code{"threshold"}. Defaults to 0.3.
#' @param alpha significance niveau. Only used if Student's t-test or bootstrap
#'   procedure is used as sparsification method. Defaults to 0.05.
#' @param adjust character indicating the method used for multiple testing
#'   adjustment (if tudent's t-test or bootstrap procedure is used for edge
#'   selection). Possible values are \code{"lfdr"} (default) for local
#'   false discovery rate correction (via \code{\link[fdrtool]{fdrtool}}),
#'   \code{"adaptBH"} for the adaptive Benjamini-Hochberg method
#'   \cite{(Benjamini and Hochberg, 2000)}, or one of the methods provided by
#'   \code{\link[stats]{p.adjust}} (see \code{p.adjust.methods()}.
#' @param trueNullMethod character indicating the method used for estimating the
#'   proportion of true null hypotheses from a vector of p-values. Used for the
#'   adaptive Benjamini-Hochberg method for multiple testing adjustment (chosen
#'   by \code{adjust = "adaptBH"}). Accepts the provided options of the
#'   \code{method} argument of \code{\link[limma]{propTrueNull}}:
#'   \code{"convest"}(default), \code{"lfdr"}, \code{"mean"}, and \code{"hist"}.
#'   Can alternatively be \code{"farco"} for
#'   the "iterative plug-in method" proposed by \cite{Farcomeni (2007)}.
#' @param lfdrThresh threshold for local FDR (if \code{adjust} is set to
#'   \code{"locfdr"}). Defaults to 0.2 meaning that associations with a
#'   corresponding local FDR less than or equal to 0.2 are identified as
#'   significant.
#' @param nboot integer indicating the number of bootstrap samples, if
#'   bootstrappng is used as sparsification method.
#' @param cores integer indicating the number of CPU cores used for
#'   bootstrapping. If cores > 1, bootstrapping is performed parallel.
#'   \code{cores} is limited to the number of available CPU cores determined by
#'   \code{\link[parallel]{detectCores}}. Then, core arguments of the function 
#'   used for association estimation (if provided) should be set to 1.
#' @param logFile if bootstrapping is used as sparsification method, a log file
#'   containing the iteration numbers is stored into the current working
#'   directory. Defaults to \code{"log.txt"}. If\code{ NULL}, no log file is
#'   created.
#' @param softThreshType character indicating the method used for transforming
#'   correlations to similarities if soft thresholding is used as sparsification
#'   method (\code{sparsMethod = "softThreshold"}). Possible values are
#'   \code{"signed"}, \code{"unsigned"}, and \code{"signed hybrid"} (according
#'   to the available options for the argument \code{type} of
#'   \code{\link[WGCNA]{adjacency}} from \code{WGCNA} package).
#' @param softThreshPower power for soft thresholding. Only used if
#'   \code{edgeSelect} is set to \code{"softThreshold"}. Expects either a single
#'   numeric value (used for a single or both networks) or a vector with two
#'   values (one for each network). If no power is set, it is computed using
#'   \code{\link[WGCNA]{pickSoftThreshold}}, where the argument
#'   \code{softThreshCut} is needed in addition.
#' @param softThreshCut numeric value between 0 and 1 indicating the desired
#'   minimum scale free topology fitting index (corresponds to the argument
#'   "RsquaredCut" in \code{\link[WGCNA]{pickSoftThreshold}}). Defaults to 0.8.
#' @param kNeighbor integer specifying the number of neighbors if the k-nearest
#'   neighbor method is used for sparsification.
#' @param knnMutual logical used for k-nearest neigbor sparsification. If
#'   \code{TRUE}, the neighbors must be mutual.
#' @param dissFunc method used for transforming associations to dissimilarities.
#'   Can be a character with one of the following values: \code{"signed"}
#'   (default), \code{"unsigned"}, \code{"signedPos"}, \code{"TOMdiss"}.
#'   Alternatively, a function is accepted with the association matrix as first
#'   argument and optional further arguments, which can be set viel
#'   \code{dissFuncPar}. Ignored for dissimilarity measures. See details.
#' @param dissFuncPar optional list with parameters if a function is passed to
#'   \code{dissFunc}.
#' @param simFunc alternative function for transforming dissimilarities to
#'   similarities. Defaults to f(x)=1-x for dissimilarities in [0,1], and
#'   f(x)=1/(1 + x) otherwise.
#' @param simFuncPar optional list with parameters for the function passed to
#'   \code{simFunc}.
#' @param scaleDiss logical. Indicates whether dissimilarity values should be
#'   scaled to [0,1] by (x - min(dissEst)) / (max(dissEst) - min(dissEst)),
#'   where dissEst is the matrix with estimated dissimilarities.
#'   Defaults to \code{TRUE}.
#' @param weighted logical. If \code{TRUE}, similarity values are used as
#'   adjacencies. \code{FALSE} leads to a binary adjacency matrix whose entries
#'   equal 1 for (sparsified) similarity values > 0, and 0 otherwise.
#' @param sampleSize integer giving the number of samples that have been used
#'   for computing the association matrix. Only needed if an association matrix
#'   is given instead of a count matrix and if, in addition, Student's t-test is
#'   used for edge selection.
#' @param verbose integer indicating the level of verbosity. Possible values:
#'   \code{"0"}: no messages, \code{"1"}: only important messages shown,
#'   \code{"2"}(default): all progress messages shown, \code{"3"} messages
#'   returned from external functions are shown in addition. Can also be logical.
#' @param seed integer giving a seed for reproducibility of the results.
#' @return An object of class \code{microNet} containing the following elements:
#'   \tabular{ll}{
#'   \code{assoMat1, assoMat2}\tab Sparsified associations (\code{NULL} for
#'   dissimlarity based networks)\cr
#'   \code{dissMat1, dissMat2}\tab Sparsified dissimilarities (for association
#'   networks, these are the sparsified associations transformed to
#'   dissimilarities)\cr
#'   \code{simMat1, simMat2}\tab Sparsified similarities\cr
#'   \code{adjamat1, adjamat2}\tab Adjacency matrices\cr
#'   \code{assoEst1, assoEst2}\tab Estimated associations (\code{NULL} for
#'   dissimlarity based networks)\cr
#'   \code{dissEst1, dissEst2}\tab Estimated dissimilarities (\code{NULL} for
#'   association networks)\cr
#'   \code{dissScale1, dissScale2}\tab Scaled dissimilarities (\code{NULL} for
#'   association networks)\cr
#'   \code{countMat1, countMat2}\tab Count matrices, where taxa and samples are
#'   filtered according to \code{filtTax} and \code{filtSamp}\cr
#'   \code{normCounts1, normCounts2}\tab Counts that are normalized according to
#'   \code{normMethod}\cr
#'   \code{groups}\tab Names of the factor levels according to which the groups
#'   have been built\cr
#'   \code{softThreshPower}\tab Determined (or given) power for
#'   soft-thresholding.\cr
#'   \code{assoType}\tab Data type (either given by \code{dataType} or
#'   determined from \code{measure})\cr
#'   \code{twoNets}\tab Indicates whether two networks have been constructed\cr
#'   \code{parameters}\tab Parameters used for network construction}
#'
#' @examples
#' # load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#' data("amgut2.filt.phy")
#'
#' # network construction:
#' amgut_net1 <- netConstruct(amgut2.filt.phy, measure = "spieceasi",
#'                            measurePar = list(method = "mb",
#'                            pulsar.params = list(rep.num = 10)),
#'                            filtTax = "highestVar",
#'                            filtTaxPar = list(highestVar = 50),
#'                            filtSamp = "totalReads",
#'                            filtSampPar = list(totalReads = 1000))
#'
#' amgut_props1 <- netAnalyze(amgut_net1, clustMethod = "cluster_fast_greedy")
#'
#' plot(amgut_props1)
#'
#' # constructing two networks according to a random group variable
#' set.seed(123456)
#' group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)
#' amgut_net2 <- netConstruct(amgut1.filt, group = group,
#'                            measure = "pearson",
#'                            filtTax = "highestVar",
#'                            filtTaxPar = list(highestVar = 50),
#'                            filtSamp = "totalReads",
#'                            filtSampPar = list(totalReads = 1000),
#'                            zeroMethod = "multRepl", normMethod = "clr")
#' amgut_props2 <- netAnalyze(amgut_net2, clustMethod = "cluster_fast_greedy")
#' plot(amgut_props2)
#'
#' @seealso \code{\link{netAnalyze}} for analyzing the constructed
#'   network(s), \code{\link{netCompare}} for network comparison,
#'   \code{\link{diffnet}} for constructing differential networks.
#' @references
#'   \insertRef{badri2018normalization}{NetCoMi}\cr
#'   \insertRef{benjamini2000adaptive}{NetCoMi}\cr
#'   \insertRef{farcomeni2007some}{NetCoMi}\cr
#'   \insertRef{friedman2012inferring}{NetCoMi} \cr
#'   \insertRef{WGCNApackage}{NetCoMi}\cr
#'   \insertRef{zhang2005general}{NetCoMi}
#' @importFrom Rdpack reprompt
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom compositions clr
#' @importFrom robCompositions cenLR
#' @importFrom vegan vegdist rrarefy
#' @importFrom Matrix nearPD
#' @importFrom zCompositions multRepl lrEM cmultRepl
#' @importFrom cccd nng
#' @importFrom stats var complete.cases pt
#' @importFrom utils capture.output
#' @importFrom WGCNA pickSoftThreshold TOMdist
#' @import phyloseq
#' @import MASS
#' @export

netConstruct <- function(data,
                         data2 = NULL,
                         dataType = "counts",
                         group = NULL,
                         measure = "spieceasi",
                         measurePar = NULL,
                         filtTax = "none",
                         filtTaxPar = NULL,
                         filtSamp = "none",
                         filtSampPar = NULL,
                         zeroMethod = "none",
                         zeroPar = NULL,
                         normMethod = "none",
                         normPar = NULL,
                         sparsMethod = "t-test",
                         thresh = 0.3,
                         alpha = 0.05,
                         adjust = "adaptBH",
                         trueNullMethod = "convest",
                         lfdrThresh = 0.2,
                         nboot = 1000L,
                         cores = 1L,
                         logFile = "log.txt",
                         softThreshType = "signed",
                         softThreshPower = NULL,
                         softThreshCut = 0.8,
                         kNeighbor = 3L,
                         knnMutual = FALSE,
                         dissFunc = "signed",
                         dissFuncPar = NULL,
                         simFunc = NULL,
                         simFuncPar = NULL,
                         scaleDiss = TRUE,
                         weighted = TRUE,
                         sampleSize = NULL,
                         verbose = 2,
                         seed = NULL) {

  dataType <-  match.arg(dataType, choices = c("counts",
                                               "phyloseq",
                                               "correlation",
                                               "partialCorr",
                                               "condDependence",
                                               "proportionality",
                                               "dissimilarity"))

  if(dataType == "phyloseq"){
    dataType <- "counts"
  }

  if (dataType =="counts") {
    
    measure <- match.arg(measure,
                         choices = c("pearson", "spearman", "bicor",
                                     "sparcc", "cclasso", "ccrepe",
                                     "propr",
                                     "spieceasi", "spring", "gcoda",
                                     "euclidean", "bray", "kld", "ckld",
                                     "jeffrey", "jsd", "aitchison"))

    if(class(data) == "phyloseq"){
      data_phy <- data
      otutab <- data@otu_table
      if(attributes(otutab)$taxa_are_rows){
        data <- t(otutab@.Data)
      } else{
        data <- otutab@.Data
      }

    } else{
      if (all(colnames(data) %in% rownames(data))) {
        warning(paste0("Row names and column names of 'data' are equal. ",
                       "Ensure 'data' is a count matrix."))
      }
    }

    assoType <-
      if (measure %in% c("pearson", "spearman", "bicor", "sparcc",
                         "cclasso", "ccrepe")) {
        "correlation"
      } else if (measure %in% c("propr")) {
        "proportionality"
      } else if (measure %in% c("spieceasi", "spring", "gcoda")) {
        "condDependence"
      } else if (measure %in% c("euclidean", "kld", "jeffrey", "jsd",
                                "ckld", "aitchison", "bray")) {
        "dissimilarity"
      }

  } else{
    measure <- "none"
    assoType <- dataType
  }

  distNet <- ifelse(assoType == "dissimilarity", TRUE, FALSE)


  filtTax <- match.arg(filtTax,
                      choices = c("totalReads", "relFreq", "numbSamp",
                                  "highestVar", "highestFreq", "none"),
                      several.ok = TRUE)

  filtSamp <- match.arg(filtSamp,
                        choices = c("totalReads", "numbTaxa", "highestFreq",
                                    "none"), several.ok = TRUE)

  sparsMethod <- match.arg(sparsMethod,
                           choices = c("none", "t-test", "bootstrap",
                                       "threshold", "softThreshold", "knn"))

  dissFunc <- match.arg(dissFunc, choices = c("signed", "unsigned", "signedPos",
                                              "TOMdiss"))

  softThreshType <- match.arg(softThreshType,
                              choices = c("signed", "unsigned", "signed hybrid"))

  zeroMethod <- match.arg(zeroMethod,
                          choices = c("none", "pseudo",
                                      "multRepl", "alrEM", "bayesMult"))

  normMethod <- match.arg(normMethod,
                          choices = c("none", "fractions", "TSS",
                                      "CSS", "COM", "rarefy", "VST", "clr"))

  adjust <- match.arg(adjust, c(p.adjust.methods, "lfdr", "adaptBH"))

  trueNullMethod <- match.arg(trueNullMethod, c("farco", "lfdr", "mean",
                                                "hist", "convest"))

  if (filtTax[1] != "none") {
    names(filtTaxPar) <- match.arg(names(filtTaxPar),
                                    choices = c( "totalReads",
                                                 "relFreq",
                                                 "numbSamp",
                                                 "highestVar",
                                                 "highestFreq"),
                                    several.ok = TRUE)
  }

  if (filtSamp[1] != "none") {
    names(filtSampPar) <- match.arg(names(filtSampPar),
                                   choices = c( "totalReads",
                                                "numbTaxa",
                                                "highestFreq"),
                                   several.ok = TRUE)
  }

  verbose <- as.numeric(verbose)


  cond <- condition_handling(dataType = dataType, assoType = assoType,
                            data2 = data2, measure = measure,
                            normMethod = normMethod, zeroMethod = zeroMethod,
                            sparsMethod = sparsMethod, dissFunc = dissFunc,
                            sampleSize = sampleSize, verbose = verbose)

  sparsMethod <- cond$sparsMethod
  dissFunc <- cond$dissFunc
  zeroMethod <- cond$zeroMethod
  normMethod <- cond$normMethod
  sampleSize <- cond$sampleSize
  needfrac <- cond$needfrac
  needint <- cond$needint

  if(cores > 1){
    parallel <- TRUE
    if(parallel::detectCores() < cores) cores <- parallel::detectCores()
  } else{
    parallel <- FALSE
  }
  
  
  if(!is.null(seed)) set.seed(seed)

 #==============================================================================

  # shall two networks be constructed?
  twoNets <- ifelse(is.null(data2) & is.null(group), FALSE, TRUE)

  if(twoNets){
    if(length(thresh) == 1) thresh <- c(thresh, thresh)
    if(length(alpha) == 1) alpha <- c(alpha, alpha)
    if(length(lfdrThresh) == 1) lfdrThresh <- c(lfdrThresh, lfdrThresh)
    if(length(softThreshPower) == 1) softThreshPower <- c(softThreshPower, softThreshPower)
    if(length(softThreshCut) == 1) softThreshCut <- c(softThreshCut, softThreshCut)
    if(length(sampleSize) == 1) sampleSize <- c(sampleSize, sampleSize)
  }

  # if data contains counts:
  if (dataType == "counts") {

    if(!(filtTax[1] == "none" & filtSamp[1] == "none") & verbose %in% 1:3){
      message("Data filtering ...")
    }

    # coerce data to numeric
    data_orig <- data
    data <-  t(apply(data_orig, 1, function(x) as.numeric(x)))
    colnames(data) <- colnames(data_orig)

    if(twoNets){
      if(!is.null(group)){
        if(distNet){
          stopifnot((is.vector(group) || is.factor(group)) &
                      length(group) == ncol(data))
          group <- as.numeric(group)
          names(group) <- colnames(data)
          if (any(is.na(data))){
            if(verbose %in% 1:3) message("Samples with NAs removed.")
            data <- data[complete.cases(data), ]
          }
        } else{
          stopifnot((is.vector(group) || is.factor(group)) &
                      length(group) == nrow(data))
          group <- as.numeric(group)
          names(group) <- rownames(data)

          if (any(is.na(data))){
            if(verbose %in% 1:3) message("Samples with NAs removed.")
            data_tmp <- cbind(data, group)
            data_tmp <- data_tmp[complete.cases(data_tmp), ]
            data <- data_tmp[, 1:(ncol(data))]
            group <- data_tmp[, ncol(data_tmp)]
          }
        }
      } else{

        if(class(data2) == "phyloseq"){
          otutab <- data2@otu_table
          if(attributes(otutab)$taxa_are_rows){
            data2 <- t(otutab@.Data)
          } else{
            data2 <- otutab@.Data
          }
        }

        if(distNet){
          if(!all(rownames(data) %in% rownames(data2))){
            if(verbose > 0) message("Intersection of samples selected.")
          }
          sel <- intersect(rownames(data), rownames(data2))

          if(length(sel) == 0) stop("Data sets contain different samples")

          data <- data[sel, ]
          data2 <- data2[sel, ]

          if (any(is.na(data)) || any(is.na(data2))){
            if(verbose %in% 1:3) message("Samples with NAs removed.")
            keep <- intersect(which(complete.cases(data)),
                              which(complete.cases(data2)))
            data <- data[keep, ]
            data2 <- data2[keep,]
          }
        } else{
          if(!identical(colnames(data), colnames(data2))){

            if(!all(colnames(data) %in% colnames(data2))){
              if(verbose > 0) message("Intersection of taxa selected.")
            }

            sel <- intersect(colnames(data), colnames(data2))

            if(length(sel) == 0) stop("Data sets contain different taxa.")

            data <- data[ , sel]
            data2 <- data2[, sel]
          }
          if (any(is.na(data)) || any(is.na(data2))){
            if(verbose %in% 1:3) message("Samples with NAs removed.")
            data <- data[complete.cases(data), ]
            data2 <- data2[complete.cases(data2), ]
          }
        }

        # coerce data2 to numeric
        data2_orig <- data2
        data2 <-  t(apply(data2_orig, 1, function(x) as.numeric(x)))
        colnames(data2) <- colnames(data2_orig)
      }
    } else{
      if (any(is.na(data))){
        if(verbose %in% 1:3) message("Samples with NAs removed.")
        data <- data[complete.cases(data), ]
      }
    }

    #============================

    keepRows1 <- filter_samples(countMat = data, filter = filtSamp,
                                filterParam = filtSampPar)

    if(twoNets){
      if(!is.null(group)){
        if(length(keepRows1) != nrow(data)){
          countMat1 <- data[keepRows1, ]
          if(!distNet) group <- group[keepRows1]
          if(verbose %in% 2:3) message(dim(data)[1] - dim(countMat1)[1],
                                              " samples removed.")
        } else{
          countMat1 <- data
        }
      } else{
        keepRows2 <- filter_samples(countMat = data2,  filter = filtSamp,
                                    filterParam = filtSampPar)
        if(distNet){
          keepRows <- intersect(keepRows1, keepRows2)
          countMat1 <- data[keepRows, ]
          countMat2 <- data2[keepRows, ]
          if(verbose %in% 2:3) message(dim(data)[1] - dim(countMat1)[1],
                                       " samples removed in each data sets.")
        } else{
          countMat1 <- data[keepRows1, ]
          countMat2 <- data2[keepRows2, ]
          if(verbose %in% 2:3){
            message(dim(data)[1] - dim(countMat1)[1],
                    " samples removed in data set 1.")
            message(dim(data2)[1] - dim(countMat2)[1],
                    " samples removed in data set 2.")
          }
        }
      }

    } else{# single network
      countMat1 <- data[keepRows1, ]
      if(verbose %in% 2:3) message(dim(data)[1] - dim(countMat1)[1],
                                   " samples removed.")
    }

    keepCols1 <- filter_taxa(countMat = countMat1, filter = filtTax,
                             filterParam = filtTaxPar)

    if(twoNets){
      if(!is.null(group)){
        if(length(keepCols1) != ncol(countMat1)){
          countMat1 <- countMat1[, keepCols1]
          if(verbose %in% 2:3) message(dim(data)[2] - dim(countMat1)[2],
                                       " taxa removed.")
          if(distNet) group <- group[keepCols1]
        }

      } else{
        keepCols2 <- filter_taxa(countMat = countMat2,  filter = filtTax,
                                 filterParam = filtTaxPar)
        if(!distNet){
          keepCols <- intersect(keepCols1, keepCols2)
          countMat1 <- countMat1[, keepCols]
          countMat2 <- countMat2[, keepCols]
          if(verbose %in% 2:3) message(dim(data)[2] - dim(countMat1)[2],
                                       " taxa removed in each data set.")
        } else{
          countMat1 <- countMat1[, keepCols1]
          countMat2 <- countMat2[, keepCols2]
          if(verbose %in% 2:3){
            message(dim(data)[2] - dim(countMat1)[2],
                    " samples removed in data set 1.")
            message(dim(data2)[2] - dim(countMat2)[2],
                    " samples removed in data set 2.")
          }
        }
      }

    } else{# single network
      countMat1 <- countMat1[, keepCols1]
      if(verbose %in% 2:3) message(dim(data)[2] - dim(countMat1)[2],
                                   " taxa removed.")
    }


    # remove samples with zero overall sum
    rs <- rowSums(countMat1)
    if (any(rs == 0)) {
      rmRows <- which(rs == 0)
      if(verbose %in% 2:3){
        message(paste(length(rmRows), "rows with zero sum removed."))
      }
      countMat1 <- countMat1[-rmRows,]
      if(!is.null(group) & !distNet & length(rmRows)!=0){
        group <- group[-rmRows]
      }
    }


    if(twoNets & !is.null(data2)){
      rs <- rowSums(countMat2)
      if (any(rs == 0)) {
        rmRows <- which(rs == 0)
        countMat2 <- countMat2[-rmRows,]

        if(distNet & length(rmRows)!=0){
          countMat1 <- countMat1[-rmRows,]
          if(verbose %in% 2:3){
            message(paste(length(rmRows),
                          "rows with zero sum removed in both data sets."))
          }
        } else{
          if(verbose %in% 2:3){
            message(paste(length(rmRows),
                          "rows with zero sum removed in data set 2."))
          }
        }
      }
    }

    if(twoNets & !is.null(data2)){
      if(verbose %in% 1:3){
        message(ncol(countMat1), " taxa and ", nrow(countMat1),
                " samples remaining in data set 1.")
        message(ncol(countMat2), " taxa and ", nrow(countMat2),
                " samples remaining in data set 2.")
      }

    } else{
      if(verbose %in% 1:3) message(ncol(countMat1), " taxa and ",
                                   nrow(countMat1), " samples remaining.")
    }



    #---------------------------------------------------------------------------
    # zero treatment
    if (zeroMethod != "none") {
      if(verbose %in% 2:3) message("\nZero treatment:")
      countMat1 <- zero_treat(countMat = countMat1, zeroMethod = zeroMethod,
                              zeroParam = zeroPar, needfrac = needfrac,
                              needint = needint, verbose = verbose)

    } else{
      attributes(countMat1)$scale <- "counts"
    }

    #---------------------------------------------------------------------------
    # normalization

    if(verbose %in% 2:3 & (normMethod != "none" || needfrac)){
      message("\nNormalization:")
    }
    count1_norm <- norm_counts(countMat = countMat1, normMethod = normMethod,
                               normParam = normPar, zeroMethod = zeroMethod,
                               needfrac = needfrac, verbose = verbose)
    countMat1 <- count1_norm
    count2_norm <- NULL
    #---------------------------------------------------------------------------

    if(!is.null(group)){

      # split count data into two groups
      if(distNet){
        splitcount <- split(as.data.frame(t(countMat1)), as.factor(group))
        if (length(splitcount) != 2)
          stop("Argument 'group' has to be binary.")
        groups <- names(splitcount)
        count1 <- t(as.matrix(splitcount[[1]]))
        count2 <- t(as.matrix(splitcount[[2]]))

      } else{
        splitcount <- split(as.data.frame(countMat1), as.factor(group))
        if (length(splitcount) != 2)
          stop("Argument 'group' has to be binary.")
        groups <- names(splitcount)
        count1 <- as.matrix(splitcount[[1]])
        count2 <- as.matrix(splitcount[[2]])
      }

      sampleSize <- c(nrow(count1), nrow(count2))

    } else if(!is.null(data2)){
      count1 <- countMat1
      count2 <- countMat2

      # zero treatment
      if (zeroMethod != "none") {
        if(verbose %in% 2:3) message("\nZero treatment for data2:")
        count2 <- zero_treat(countMat = count2, zeroMethod = zeroMethod,
                             zeroParam = zeroPar, needfrac = needfrac,
                             needint = needint, verbose = verbose)

      } else{
        attributes(count2)$scale <- "counts"
      }

      # normalization
      if(verbose %in% 2:3 & (normMethod != "none" || needfrac)){
        message("\nNormalization for data2:")
      }
      count2_norm <- norm_counts(countMat = count2, normMethod = normMethod,
                                 normParam = normPar, zeroMethod = zeroMethod,
                                 needfrac = needfrac, verbose = verbose)
      count2 <- count2_norm

      if(distNet){
        keep <- intersect(rownames(count1), rownames(count2))
        count1 <- count1[keep, ]
        count2 <- count2[keep, ]
      }else{
        keep <- intersect(colnames(count1), colnames(count2))
        count1 <- count1[, keep]
        count2 <- count2[, keep]
      }

      groups <- c("1", "2")
      sampleSize <- c(nrow(count1), nrow(count2))

    } else{
      count1 <- countMat1
      count2 <- NULL
      groups <- NULL
      sampleSize <- nrow(count1)
    }



    if(verbose %in% 2:3){
      if(distNet){
        txt.tmp <- "dissimilarities"
      } else{
        txt.tmp <- "associations"
      }
      message("\nCalculate '", measure, "' ", txt.tmp, " ... ", appendLF = FALSE)
    }

    assoMat1 <- calc_association(countMat = count1, measure = measure,
                                 measurePar = measurePar, verbose = verbose)

    if(verbose %in% 2:3) message("Done.")


    if(twoNets){

      if(verbose %in% 2:3){
        message("\nCalculate ", txt.tmp, " for group 2 ... ",
                appendLF = FALSE)
      }

      assoMat2 <- calc_association(countMat = count2, measure = measure,
                                   measurePar = measurePar, verbose = verbose)
      if(verbose %in% 2:3) message("Done.")

    } else{
      assoMat2 <- NULL
    }

  } else{
    assoMat1 <- data
    assoMat2 <- data2
    count1 <- NULL
    count2 <- NULL
    count1_norm <- NULL
    count2_norm <- NULL
    groups <- NULL
  }

  if(distNet){

    dissEst1 <- assoMat1
    
    if(any(is.infinite(dissEst1)) & scaleDiss){
      scaleDiss <- FALSE
      warning("Dissimilarity matrix contains infinite values and cannot be scaled to [0,1].")
    }

    if(scaleDiss){
      assoUpper <- assoMat1[upper.tri(assoMat1)]
      assoMat1 <- (assoMat1 - min(assoUpper)) / (max(assoUpper) - min(assoUpper))
      diag(assoMat1) <- 0
    }

    dissScale1 <- assoMat1

    if(verbose %in% 2:3){
      if(sparsMethod != "none"){
        message("\nSparsify dissimilarities via '", sparsMethod, "' ... ",
                appendLF = FALSE)
      }
    }

    sparsReslt <- sparsify(assoMat = assoMat1, countMat = count1,
                           sampleSize = sampleSize[1], measure = measure,
                           measurePar = measurePar,
                           assoType = assoType, sparsMethod = sparsMethod,
                           thresh = thresh[1], alpha = alpha[1], adjust = adjust,
                           lfdrThresh = lfdrThresh[1],
                           trueNullMethod = trueNullMethod, nboot = nboot,
                           softThreshType = softThreshType,
                           softThreshPower = softThreshPower[1],
                           softThreshCut = softThreshCut[1], parallel = parallel,
                           cores = cores, logFile = logFile, kNeighbor = kNeighbor,
                           knnMutual = knnMutual, verbose = verbose, seed = seed)
    if(verbose %in% 2:3 & sparsMethod != "none") message("Done.")

    assoMat1 <- NULL
    assoEst1 <- NULL
    dissMat1 <- sparsReslt$assoNew
    power1 <- sparsReslt$power
    simMat1 <- trans_to_sim(x = dissMat1, simFunc = simFunc, simFuncPar = simFuncPar)
    adjaMat1 <- trans_to_adja(x = simMat1, weighted = weighted)

    if(twoNets){
      dissEst2 <- assoMat2

      if(scaleDiss){
        assoUpper <- assoMat2[upper.tri(assoMat2)]
        assoMat2 <- (assoMat2 - min(assoUpper)) / (max(assoUpper) - min(assoUpper))
      }
      dissScale2 <- assoMat2

      if(verbose %in% 2:3){
        if(sparsMethod != "none"){
          message("\nSparsify dissimilarities for group 2 ... ",
                  appendLF = FALSE)
        }
      }
      sparsReslt <- sparsify(assoMat = assoMat2, countMat = count2,
                             sampleSize = sampleSize[2], measure = measure,
                             measurePar = measurePar,
                             assoType = assoType, sparsMethod = sparsMethod,
                             thresh = thresh[2], alpha = alpha[2], adjust = adjust,
                             lfdrThresh = lfdrThresh[2],
                             trueNullMethod = trueNullMethod, nboot = nboot,
                             softThreshType = softThreshType,
                             softThreshPower = softThreshPower[2],
                             softThreshCut = softThreshCut[2], parallel = parallel,
                             cores = cores, logFile = logFile, kNeighbor = kNeighbor,
                             knnMutual = knnMutual, verbose = verbose, seed = seed)
      if(verbose %in% 2:3 & sparsMethod != "none") message("Done.")

      assoMat2 <- NULL
      assoEst2 <- NULL
      dissMat2 <- sparsReslt$assoNew
      power2 <- sparsReslt$power
      simMat2 <- trans_to_sim(x = dissMat2, simFunc = simFunc,
                              simFuncPar = simFuncPar)
      adjaMat2 <- trans_to_adja(x = simMat2, weighted = weighted)

    } else{
      dissEst2 <- dissScale2 <- assoMat2 <- assoEst2 <- dissMat2 <-
        power2 <- simMat2 <- adjaMat2 <- NULL
    }


  } else{ # association network
    if(verbose %in% 2:3){
      if(sparsMethod != "none"){
        message("\nSparsify associations via '", sparsMethod, "' ... ",
                appendLF = FALSE)
      }
    }
    sparsReslt <- sparsify(assoMat = assoMat1, countMat = count1,
                           sampleSize = sampleSize[1], measure = measure,
                           measurePar = measurePar,
                           assoType = assoType, sparsMethod = sparsMethod,
                           thresh = thresh[1], alpha = alpha[1], adjust = adjust,
                           lfdrThresh = lfdrThresh[1],
                           trueNullMethod = trueNullMethod, nboot = nboot,
                           softThreshType = softThreshType,
                           softThreshPower = softThreshPower[1],
                           softThreshCut = softThreshCut[1], parallel = parallel,
                           cores = cores, logFile = logFile, kNeighbor = kNeighbor,
                           knnMutual = knnMutual, verbose = verbose, seed = seed)
    if(verbose %in% 2:3 & sparsMethod != "none") message("Done.")
    assoEst1 <- assoMat1
    assoMat1 <- sparsReslt$assoNew
    power1 <- sparsReslt$power
    dissEst1 <- dissScale1 <- NULL

    dissMat1 <- trans_to_diss(x = assoMat1, dissFunc = dissFunc,
                            dissFuncPar = dissFuncPar)

    if(sparsMethod == "softThreshold"){
      simMat1 <- sparsReslt$simMat
      adjaMat1 <- assoMat1
    } else{
      simMat1 <- trans_to_sim(x = dissMat1, simFunc = simFunc,
                              simFuncPar = simFuncPar)

      adjaMat1 <- trans_to_adja(x = simMat1, weighted = weighted)
    }

    if(twoNets){
      if(verbose %in% 2:3){
        if(sparsMethod != "none"){
          message("\nSparsify associations for group 2 ... ",
                  appendLF = FALSE)
        }
      }
      sparsReslt <- sparsify(assoMat = assoMat2, countMat = count2,
                             sampleSize = sampleSize[2], measure = measure,
                             measurePar = measurePar,
                             assoType = assoType, sparsMethod = sparsMethod,
                             thresh = thresh[2], alpha = alpha[2],
                             adjust = adjust, lfdrThresh = lfdrThresh[2],
                             trueNullMethod = trueNullMethod, nboot = nboot,
                             softThreshType = softThreshType,
                             softThreshPower = softThreshPower[2],
                             softThreshCut = softThreshCut[2],
                             parallel = parallel, cores = cores,
                             logFile = logFile, kNeighbor = kNeighbor,
                             knnMutual = knnMutual, verbose = verbose,
                             seed = seed)
      if(verbose %in% 2:3 & sparsMethod != "none") message("Done.")
      assoEst2 <- assoMat2
      assoMat2 <- sparsReslt$assoNew
      power2 <- sparsReslt$power
      dissEst2 <- dissScale2 <- NULL

      dissMat2 <- trans_to_diss(x = assoMat2, dissFunc = dissFunc,
                              dissFuncPar = dissFuncPar)

      if(sparsMethod == "softThreshold"){
        simMat2 <- sparsReslt$simMat
        adjaMat2 <- assoMat2
      } else{
        simMat2 <- trans_to_sim(x = dissMat2, simFunc = simFunc,
                                simFuncPar = simFuncPar)

        adjaMat2 <- trans_to_adja(x = simMat2, weighted = weighted)
      }

    } else{
      dissEst2 <- dissScale2 <- assoMat2 <- assoEst2 <- dissMat2 <-
        power2 <- simMat2 <- adjaMat2 <- NULL
    }

  }


  #=============================================================================
  output <- list()
  output$assoMat1 <- assoMat1
  output$assoMat2 <- assoMat2
  output$dissMat1 <- dissMat1
  output$dissMat2 <- dissMat2
  output$simMat1 <- simMat1
  output$simMat2 <- simMat2
  output$adjaMat1 <- adjaMat1
  output$adjaMat2 <- adjaMat2

  output$assoEst1 <- assoEst1
  output$assoEst2 <- assoEst2
  output$dissEst1 <- dissEst1
  output$dissEst2 <- dissEst2
  output$dissScale1 <- dissScale1
  output$dissScale2 <- dissScale2

  output$countMat1 <- count1
  output$countMat2 <- count2
  output$normCounts1 <- count1_norm
  output$normCounts2 <- count2_norm
  output$groups <- groups
  output$sampleSize <- sampleSize
  output$softThreshPower <- list( power1 = power1, power2 = power2) # calculated power
  output$assoType <- assoType
  output$twoNets <- twoNets

  output$parameters <- list(
    dataType = dataType,
    group = group,
    filtTax = filtTax,
    filtTaxPar = filtTaxPar,
    filtSamp = filtSamp,
    filtSampPar = filtSampPar,
    zeroMethod = zeroMethod,
    zeroPar = zeroPar,
    normMethod = normMethod,
    normPar = normPar,
    measure = measure,
    measurePar = measurePar,
    sparsMethod = sparsMethod,
    thresh = thresh,
    adjust = adjust,
    alpha = alpha,
    lfdrThresh = lfdrThresh,
    nboot = nboot,
    softThreshType = softThreshType,
    softThreshPower = softThreshPower,
    softThreshCut = softThreshCut,
    kNeighbor = kNeighbor,
    knnMutual = knnMutual,
    dissFunc = dissFunc,
    dissFuncPar = dissFuncPar,
    simFunc = simFunc,
    simFuncPar = simFuncPar,
    scaleDiss = scaleDiss,
    weighted = weighted,
    sampleSize = sampleSize
  )

  class(output) <- "microNet"
  return(output)

}
