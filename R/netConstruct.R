#' @title Constructing Networks for Microbiome Data
#'
#' @description Constructing microbial association networks and dissimilarity
#'   based networks (where nodes are subjects) from compositional count data.
#'   
#' @usage netConstruct(data,
#'              data2 = NULL,
#'              dataType = "counts",
#'              group = NULL,
#'              matchDesign = NULL,
#'              taxRank = NULL,
#'              
#'              # Association/dissimilarity measure:
#'              measure = "spieceasi",
#'              measurePar = NULL,
#'              
#'              # Preprocessing:
#'              jointPrepro = NULL,
#'              filtTax = "none",
#'              filtTaxPar = NULL,
#'              filtSamp = "none",
#'              filtSampPar = NULL,
#'              zeroMethod = "none",
#'              zeroPar = NULL,
#'              normMethod = "none",
#'              normPar = NULL,
#'              
#'              # Sparsification:
#'              sparsMethod = "t-test",
#'              thresh = 0.3,
#'              alpha = 0.05,
#'              adjust = "adaptBH",
#'              trueNullMethod = "convest",
#'              lfdrThresh = 0.2,
#'              nboot = 1000L,
#'              assoBoot = NULL,
#'              cores = 1L,
#'              logFile = "log.txt",
#'              softThreshType = "signed",
#'              softThreshPower = NULL,
#'              softThreshCut = 0.8,
#'              kNeighbor = 3L,
#'              knnMutual = FALSE,
#'              
#'              # Transformation:
#'              dissFunc = "signed",
#'              dissFuncPar = NULL,
#'              simFunc = NULL,
#'              simFuncPar = NULL,
#'              scaleDiss = TRUE,
#'              weighted = TRUE,
#'              
#'              # Further arguments:
#'              sampleSize = NULL,
#'              verbose = 2,
#'              seed = NULL
#'              )
#'
#' @details The object returned by \code{netConstruct} can either be passed to
#'   \code{\link{netAnalyze}} for network analysis, or to
#'   \code{\link{diffnet}} to construct a differential network from the
#'   estimated associations.
#'   \cr
#'   The function enables the construction of either a \strong{single network} 
#'   or \strong{two networks}. The latter can be compared using the function
#'   \code{\link{netCompare}}. 
#'   \cr\cr
#'   The network(s) can either be based on \strong{associations} (correlation,
#'   partial correlation / conditional dependence, proportionality) or 
#'   \strong{dissimilarities}. Several measures are available, respectively, 
#'   to estimate associations or dissimilarities using \code{netConstruct}. 
#'   Alternatively, a pre-generated association or dissimilarity matrix is 
#'   accepted as input to start the workflow (argument \code{dataType} must be 
#'   set appropriately). 
#'   Depending on the measure, network nodes are either taxa or subjects: 
#'   In association-based networks nodes are taxa, whereas in 
#'   dissimilarity-based networks nodes are subjects.
#'   \cr
#'   \cr
#'   In order to perform a \strong{network comparison}, the following options 
#'   for constructing two networks are available:
#'   \enumerate{
#'     \item Passing the combined count matrix to \code{data} and a group 
#'     vector to \code{group} (of length \code{nrow(data)} for association 
#'     networks and of length \code{ncol(data)} for dissimilarity-based 
#'     networks).
#'     \item Passing the count data for group 1 to \code{data} (matrix or
#'     phyloseq object) and the count data for group 2 to \code{data2} (matrix 
#'     or phyloseq object). For association networks, the column names must 
#'     match, and for dissimilarity networks the row names.
#'     \item Passing an association/dissimilarity matrix for group 1 to 
#'     \code{data} and an association/dissimilarity matrix for group 2 to 
#'     \code{data2}.
#'   }
#'   \cr
#'   \strong{Group labeling:}\cr 
#'   If two networks are generated, the network belonging to \code{data} 
#'   is always denoted by "group 1" and the network belonging to \code{data2} 
#'   by "group 2".\cr
#'   If a group vector is used for splitting the data into two groups, the group
#'   names are assigned according to the order of group levels. If \code{group}
#'   contains the levels 0 and 1, for instance, "group 1" is assigned to level 0
#'   and "group 2" is assigned to level 1. \cr
#'   In the network plot, group 1 is shown on the left and group 2 on the 
#'   right if not defined otherwise (see \code{\link{plot.microNetProps}}).
#'   \cr\cr
#'   \strong{Association measures}
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
#'   \strong{Dissimilarity measures}
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
#'   \strong{Methods for zero replacement}
#'   \tabular{lll}{
#'   Argument \tab Method \tab Function\cr
#'   \code{"none"} \tab No zero replacement (only available if no zero
#'   replacement is needed for the chosen normalization method and 
#'   association/dissimilarity measure).\tab -\cr
#'   \code{"pseudo"} \tab A pseudo count (defined by \code{pseudocount} as 
#'   optional element of \code{zeroPar}) is added to all counts. A unit zero 
#'   count is used by default.\tab -\cr
#'   \code{"pseudoZO"} \tab A pseudo count (defined by \code{pseudocount} as 
#'   optional element of \code{zeroPar}) is added to zero counts only. 
#'   A unit zero count is used by default.\tab -\cr
#'   \code{"multRepl"} \tab Multiplicative simple replacement \tab
#'  \code{\link[zCompositions:multRepl]{multRepl}}\cr
#'   \code{"alrEM"} \tab Modified EM alr-algorithm 
#'   \tab \code{\link[zCompositions:lrEM]{lrEM}}\cr
#'   \code{"bayesMult"} \tab Bayesian-multiplicative replacement\tab
#'   \code{\link[zCompositions:cmultRepl]{cmultRepl}}\cr
#'   }
#'   \strong{Normalization methods}
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
#'   \code{\link[SpiecEasi]{clr}}\cr
#'   \code{"mclr"} \tab Modified central log ratio transformation\tab
#'   \code{\link[SPRING]{mclr}}
#'   }
#'   These methods (except for rarefying) are described in
#'   \cite{Badri et al.(2020)}.
#'   \cr\cr
#'   \strong{Transformation methods}\cr
#'   Functions used for transforming associations into dissimilarities:
#'   \tabular{ll}{
#'   Argument \tab Function\cr
#'   \code{"signed"} \tab sqrt(0.5 * (1 - x))\cr
#'   \code{"unsigned"} \tab sqrt(1 - x^2)\cr
#'   \code{"signedPos"} \tab  diss <- sqrt(0.5 * (1-x))\cr
#'     \tab diss[x < 0] <- 0\cr
#'   \code{"TOMdiss"} \tab \code{\link[WGCNA:TOMsimilarity]{TOMdist}}
#'   }
#' @param data numeric matrix. Can be a count matrix (rows are samples, columns
#'   are OTUs/taxa), a phyloseq object, or an association/dissimilarity matrix
#'   (\code{dataType} must be set).
#'   a second count matrix/phyloseq object or a second association/dissimilarity
#'   matrix.
#' @param data2 optional numeric matrix used for constructing a second 
#'   network (belonging to group 2). Can be either a second count 
#'   matrix/phyloseq object or a second association/dissimilarity matrix.
#' @param dataType character indicating the data type. Defaults to "counts",
#'   which means that \code{data} (and data2) is a count matrix or object of
#'   class \code{\link[phyloseq:phyloseq-class]{phyloseq}}. Further options
#'   are "correlation", "partialCorr" (partial correlation), "condDependence"
#'   (conditional dependence), "proportionality" and "dissimilarity".
#' @param group optional binary vector used for splitting the data into two
#'   groups. If \code{group} is \code{NULL} (default) and \code{data2} is not 
#'   set, a single network is constructed. See 'Details.'
#' @param matchDesign Numeric vector with two elements specifying an optional 
#'   matched-group (i.e. matched-pair) design, which is used for the permutation 
#'   tests in \code{\link{netCompare}} and \code{\link{diffnet}}. \code{c(1,1)} 
#'   corresponds to a matched-pair design. A 1:2 matching, for instance, is 
#'   defined by \code{c(1,2)}, which means that the first sample of group 1 is 
#'   matched to the first two samples of group 2 and so on. 
#'   The appropriate order of samples must be ensured. If 
#'   \code{NULL}, the group memberships are shuffled randomly while group sizes
#'   identical to the original data set are ensured. 
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
#'   associations/dissimilarities. See 'Details' for the respective functions. 
#'   For SpiecEasi or SPRING as association measure, an additional list element
#'   "symBetaMode" is accepted to define the "mode" argument of 
#'   \code{\link[SpiecEasi]{symBeta}}.
#' @param jointPrepro logical indicating whether data preprocessing (filtering,
#'   zero treatment, normalization) should be done for the combined data sets,
#'   or each data set separately. Ignored if a single network is constructed.
#'   Defaults to \code{TRUE} if \code{group} is given, and to \code{FALSE} if
#'   \code{data2} is given. Joint preprocessing is not possible for 
#'   dissimilarity networks. 
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
#'   \item{\code{"highestFreq"}}{Keep the x taxa with highest frequency.}
#'   }
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
#' @param taxRank character indicating the taxonomic rank at which the network 
#'   should be constructed. Only used if data (and data 2) is a phyloseq object. 
#'   The given rank must match one of the column names of the taxonomic table 
#'   (the \code{@tax_table} slot of the phyloseq object). Taxa names of the 
#'   chosen taxonomic rank must be unique (consider using the function 
#'   \code{\link{renameTaxa}} to make them unique). If a phyloseq object is 
#'   given and \code{taxRank = NULL}, the row names of the OTU table are used 
#'   as node labels.
#' @param zeroMethod character indicating the method used for zero replacement.
#'   Possible values are: \code{"none"} (default), \code{"pseudo"}, 
#'   \code{"pseudoZO"}, \code{"multRepl"}, \code{"alrEM"}, \code{"bayesMult"}. 
#'   See 'Details'. The corresponding parameters are set via \code{zeroPar}. 
#'   \code{zeroMethod} is ignored if the approach for
#'   calculating the associations/dissimilarity includes zero handling.
#'   Defaults to \code{"multRepl"} or \code{"pseudo"} (depending on the expected
#'   input of the normalization function and measure) if zero replacement is
#'   required.
#' @param zeroPar list with parameters passed to the function for zero
#'   replacement (\code{zeroMethod}). See the help page of the respective
#'   function for details. If \code{zeroMethod} is \code{"pseudo"} or 
#'   \code{"pseudoZO"}, the pseudo count can be specified via 
#'   \code{zeroPar = list(pseudocount = x)} (where x is numeric).
#' @param normMethod character indicating the normalization method (to
#'   make counts of different samples comparable). Possible options are:
#'   \code{"none"} (default), \code{"TSS"} (or \code{"fractions"}), \code{"CSS"},
#'   \code{"COM"}, \code{"rarefy"}, \code{"VST"}, \code{"clr"}, and 
#'   \code{"mclr"}. See 'Details'. The corresponding parameters are set via 
#'   \code{normPar}.
#' @param normPar list with parameters passed to the function for normalization
#'   (defined by \code{normMethod}).
#' @param sparsMethod character indicating the method used for sparsification
#'   (selected edges that are connected in the network). Available methods are:
#'   \describe{
#'   \item{\code{"none"}}{Leads to a fully connected network}
#'   \item{\code{"t-test"}}{Default. Associations being significantly different 
#'   from zero are selected using Student's t-test. Significance level and 
#'   multiple testing adjustment is specified via \code{alpha} and 
#'   \code{adjust}. \code{sampleSize} must be set if \code{dataType} is not 
#'   "counts".}
#'   \item{\code{"bootstrap"}}{Bootstrap procedure as described in
#'   \cite{Friedman and Alm (2012)}. Corresponding arguments are
#'   \code{nboot}, \code{cores}, and \code{logFile}. Data type must be 
#'   "counts".}
#'   \item{\code{"threshold"}}{Selected are taxa
#'   pairs with an absolute association/dissimilarity greater than or equal to
#'   the threshold defined via \code{thresh}.}
#'   \item{\code{"softThreshold"}}{Soft thresholding method according to
#'   \cite{Zhang and Horvath (2005)} available in the
#'   \code{\link[WGCNA:pickSoftThreshold]{WGCNA}} package. Corresponding
#'   arguments are \code{softThreshType}, \code{softThreshPower}, and
#'   \code{softThreshCut}.}
#'   \item{\code{"knn"}}{Construct a k-nearest neighbor or mutual k-nearest
#'   neighbor graph using \code{\link[cccd]{nng}}. Corresponding
#'   arguments are \code{kNeighbor}, and \code{knnMutual}. Available for 
#'   dissimilarity networks only.}}
#' @param thresh numeric vector with one or two elements defining the threshold 
#'   used for sparsification if \code{sparsMethod} is set to \code{"threshold"}. 
#'   If two networks are constructed and one value is given, it is used 
#'   for both groups. Defaults to 0.3.
#' @param alpha numeric vector with one or two elements indicating the 
#'   significance level. Only used if Student's t-test or bootstrap
#'   procedure is used as sparsification method. If two networks are constructed 
#'   and one value is given, it is used for both groups. Defaults to 0.05.
#' @param adjust character indicating the method used for multiple testing
#'   adjustment (if Student's t-test or bootstrap procedure is used for edge
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
#' @param lfdrThresh numeric vector with one or two elements defining the 
#'   threshold(s) for local FDR correction (if \code{adjust = "locfdr"}). 
#'   Defaults to 0.2 meaning that associations with a
#'   corresponding local FDR less than or equal to 0.2 are identified as
#'   significant. If two networks are constructed and one value is given, it is 
#'   used for both groups.
#' @param nboot integer indicating the number of bootstrap samples, if
#'   bootstrapping is used as sparsification method.
#' @param assoBoot logical or list. Only relevant for bootstrapping. 
#'   Set to \code{TRUE} if a list (\code{assoBoot}) with bootstrap association 
#'   matrices should be returned. Can also be a list with bootstrap 
#'   association matrices, which are used for sparsification. See the example.
#' @param cores integer indicating the number of CPU cores used for
#'   bootstrapping. If cores > 1, bootstrapping is performed parallel.
#'   \code{cores} is limited to the number of available CPU cores determined by
#'   \code{\link[parallel]{detectCores}}. Then, core arguments of the function 
#'   used for association estimation (if provided) should be set to 1.
#' @param logFile character defining a log file to which the iteration numbers 
#'   are stored if bootstrapping is used for sparsification. The file is written 
#'   to the current working directory. Defaults to \code{"log.txt"}. 
#'   If\code{ NULL}, no log file is created.
#' @param softThreshType character indicating the method used for transforming
#'   correlations into similarities if soft thresholding is used as sparsification
#'   method (\code{sparsMethod = "softThreshold"}). Possible values are
#'   \code{"signed"}, \code{"unsigned"}, and \code{"signed hybrid"} (according
#'   to the available options for the argument \code{type} of
#'   \code{\link[WGCNA]{adjacency}} from \code{WGCNA} package).
#' @param softThreshPower numeric vector with one or two elements defining the 
#'   power for soft thresholding. Only used if
#'   \code{edgeSelect = "softThreshold"}. If two networks are constructed and 
#'   one value is given, it is used for both groups. If no power is set, it is 
#'   computed using \code{\link[WGCNA]{pickSoftThreshold}}, where the argument
#'   \code{softThreshCut} is needed in addition.
#' @param softThreshCut numeric vector with one or two elements (each between 0 
#'   and 1) indicating the desired minimum scale free topology fitting index 
#'   (corresponds to the argument "RsquaredCut" in 
#'   \code{\link[WGCNA]{pickSoftThreshold}}). Defaults to 0.8.
#'   If two networks are constructed and one value is given, it is 
#'   used for both groups.
#' @param kNeighbor integer specifying the number of neighbors if the k-nearest
#'   neighbor method is used for sparsification. Defaults to 3L.
#' @param knnMutual logical used for k-nearest neighbor sparsification. If
#'   \code{TRUE}, the neighbors must be mutual. Defaults to \code{FALSE}.
#' @param dissFunc method used for transforming associations into 
#'   dissimilarities. Can be a character with one of the following values: 
#'   \code{"signed"}(default), \code{"unsigned"}, \code{"signedPos"}, 
#'   \code{"TOMdiss"}.
#'   Alternatively, a function is accepted with the association matrix as first
#'   argument and optional further arguments, which can be set via
#'   \code{dissFuncPar}. Ignored for dissimilarity measures. See 'Details.'
#' @param dissFuncPar optional list with parameters if a function is passed to
#'   \code{dissFunc}.
#' @param simFunc function for transforming dissimilarities into
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
#' @param sampleSize numeric vector with one or two elements giving the number 
#'   of samples that have been used for computing the association matrix. 
#'   Only needed if an association matrix
#'   is given instead of a count matrix and if, in addition, Student's t-test is
#'   used for edge selection. If two networks are constructed and one value is 
#'   given, it is used for both groups.
#' @param verbose integer indicating the level of verbosity. Possible values:
#'   \code{"0"}: no messages, \code{"1"}: only important messages,
#'   \code{"2"}(default): all progress messages, \code{"3"} messages returned 
#'   by external functions are shown in addition. Can also be logical.
#' @param seed integer giving a seed for reproducibility of the results.
#' @return An object of class \code{microNet} containing the following elements:
#'   \tabular{ll}{
#'   \code{edgelist1, edgelist2}\tab Edge list with the following columns:
#'   \itemize{
##'  \item \code{v1}, \code{v2}: names of adjacent nodes/vertices
##'  \item \code{asso}: estimated association (only for association networks)
##'  \item \code{diss}: dissimilarity
##'  \item \code{sim}: similarity (only for unweighted networks)
##'  \item \code{adja}: adjacency (equals similarity for 
#'   weighted networks)
##'  }\cr
#'   \code{assoMat1, assoMat2}\tab Sparsified associations (\code{NULL} for
#'   dissimilarity based networks)\cr
#'   \code{dissMat1, dissMat2}\tab Sparsified dissimilarities (for association
#'   networks, these are the sparsified associations transformed into
#'   dissimilarities)\cr
#'   \code{simMat1, simMat2}\tab Sparsified similarities\cr
#'   \code{adjaMat1, adjaMat2}\tab Adjacency matrices\cr
#'   \code{assoEst1, assoEst2}\tab Estimated associations (\code{NULL} for
#'   dissimilarity based networks)\cr
#'   \code{dissEst1, dissEst2}\tab Estimated dissimilarities (\code{NULL} for
#'   association networks)\cr
#'   \code{dissScale1, dissScale2}\tab Scaled dissimilarities (\code{NULL} for
#'   association networks)\cr
#'   \code{assoBoot1, assoBoot2}\tab List with association matrices for the 
#'   bootstrap samples. Returned if bootstrapping is used for 
#'   sparsification and \code{assoBoot} is \code{TRUE}.\cr
#'   \code{countMat1, countMat2}\tab Count matrices after filtering but before
#'   zero replacement and normalization. Only returned if \code{jointPrepro}
#'   is \code{FALSE} or for a single network.\cr
#'   \code{countsJoint}\tab Joint count matrix after filtering but before
#'   zero replacement and normalization. Only returned if \code{jointPrepro}
#'   is \code{TRUE}.\cr
#'   \code{normCounts1, normCounts2}\tab Count matrices after zero handling and
#'   normalization\cr
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
#' # Load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#' data("amgut2.filt.phy")
#'
#' # Single network with the following specifications:
#' # - Association measure: SpiecEasi
#' # - SpiecEasi parameters are defined via 'measurePar' 
#' #   (check ?SpiecEasi::spiec.easi for available options)
#' # - Note: 'rep.num' should be higher for real data sets
#' # - Taxa filtering: Keep the 50 taxa with highest variance
#' # - Sample filtering: Keep samples with a total number of reads 
#' #   of at least 1000
#' 
#' net1 <- netConstruct(amgut2.filt.phy, 
#'                      measure = "spieceasi",
#'                      measurePar = list(method = "mb",
#'                                        pulsar.params = list(rep.num = 10),
#'                                        symBetaMode = "ave"),
#'                      filtTax = "highestVar",
#'                      filtTaxPar = list(highestVar = 50),
#'                      filtSamp = "totalReads",
#'                      filtSampPar = list(totalReads = 1000),
#'                      sparsMethod = "none",
#'                      normMethod = "none",
#'                      verbose = 3)
#'
#' # Network analysis (see ?netAnalyze for details)
#' props1 <- netAnalyze(net1, clustMethod = "cluster_fast_greedy")
#'
#' # Network plot (see ?plot.microNetProps for details)
#' plot(props1)
#' 
#' #----------------------------------------------------------------------------
#' # Same network as before but on genus level and without taxa filtering
#' 
#' amgut.genus.phy <- phyloseq::tax_glom(amgut2.filt.phy, taxrank = "Rank6")
#' 
#' dim(phyloseq::otu_table(amgut.genus.phy))
#' 
#' # Rename taxonomic table and make Rank6 (genus) unique
#' amgut.genus.renamed <- renameTaxa(amgut.genus.phy, pat = "<name>", 
#'                                   substPat = "<name>_<subst_name>(<subst_R>)",
#'                                   numDupli = "Rank6")
#'                                   
#' net_genus <- netConstruct(amgut.genus.renamed, 
#'                           taxRank = "Rank6",
#'                           measure = "spieceasi",
#'                           measurePar = list(method = "mb",
#'                                             pulsar.params = list(rep.num = 10),
#'                                             symBetaMode = "ave"),
#'                           filtSamp = "totalReads",
#'                           filtSampPar = list(totalReads = 1000),
#'                           sparsMethod = "none",
#'                           normMethod = "none",
#'                           verbose = 3)
#' 
#' # Network analysis
#' props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")
#' 
#' # Network plot (with some modifications)
#' plot(props_genus, 
#'      shortenLabels = "none",
#'      labelScale = FALSE,
#'      cexLabels = 0.8)
#' 
#' #----------------------------------------------------------------------------
#' # Single network with the following specifications:
#' # - Association measure: Pearson correlation
#' # - Taxa filtering: Keep the 50 taxa with highest frequency
#' # - Sample filtering: Keep samples with a total number of reads of at least 
#' #  1000 and with at least 10 taxa with a non-zero count
#' # - Zero replacement: A pseudo count of 0.5 is added to all counts
#' # - Normalization: clr transformation
#' # - Sparsification: Threshold = 0.3 
#' #  (an edge exists between taxa with an estimated association >= 0.3)
#' 
#' net2 <- netConstruct(amgut2.filt.phy, 
#'                      measure = "pearson",
#'                      filtTax = "highestFreq",
#'                      filtTaxPar = list(highestFreq = 50),
#'                      filtSamp = c("numbTaxa", "totalReads"),
#'                      filtSampPar = list(totalReads = 1000, numbTaxa = 10),
#'                      zeroMethod = "pseudo", 
#'                      zeroPar = list(pseudocount = 0.5),
#'                      normMethod = "clr",
#'                      sparsMethod = "threshold",
#'                      thresh = 0.3,
#'                      verbose = 3)
#'
#' # Network analysis
#' props2 <- netAnalyze(net2, clustMethod = "cluster_fast_greedy")
#'
#' plot(props2)
#' 
#' #----------------------------------------------------------------------------
#' # Constructing and analyzing two networks 
#' # - A random group variable is used for splitting the data into two groups
#' 
#' set.seed(123456)
#' group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)
#' 
#' # Option 1: Use the count matrix and group vector as input:
#' net3 <- netConstruct(amgut1.filt,  
#'                      group = group,
#'                      measure = "pearson",
#'                      filtTax = "highestVar",
#'                      filtTaxPar = list(highestVar = 50),
#'                      filtSamp = "totalReads",
#'                      filtSampPar = list(totalReads = 1000),
#'                      zeroMethod = "multRepl", 
#'                      normMethod = "clr",
#'                      sparsMethod = "t-test")
#'                      
#' # Option 2: Pass the count matrix of group 1 to 'data' 
#' #           and that of group 2 to 'data2'
#' # Note: Argument 'jointPrepro' is set to FALSE by default (the data sets 
#' # are filtered separately and the intersect of filtered taxa is kept, 
#' # which leads to less than 50 taxa in this example).
#' 
#' amgut1 <- amgut1.filt[group == 1, ]
#' amgut2 <- amgut1.filt[group == 2, ]
#' 
#' net3 <- netConstruct(data = amgut1, 
#'                      data2 = amgut2,
#'                      measure = "pearson",
#'                      filtTax = "highestVar",
#'                      filtTaxPar = list(highestVar = 50),
#'                      filtSamp = "totalReads",
#'                      filtSampPar = list(totalReads = 1000),
#'                      zeroMethod = "multRepl", 
#'                      normMethod = "clr",
#'                      sparsMethod = "t-test")
#'                            
#' # Network analysis
#' # Note: Please zoom into the GCM plot or open a new window using:
#' # x11(width = 10, height = 10)
#' props3 <- netAnalyze(net3, clustMethod = "cluster_fast_greedy")
#' 
#' # Network plot (same layout is used in both groups)
#' plot(props3, sameLayout = TRUE)
#' 
#' # The two networks can be compared with NetCoMi's function netCompare().
#' 
#' #----------------------------------------------------------------------------
#' # Example of using the argument "assoBoot"
#' 
#' # This functionality is useful for splitting up a large number of bootstrap 
#' # replicates and run the bootstrapping procedure iteratively.
#' 
#' niter <- 5
#' nboot <- 1000
#' # Overall number of bootstrap replicates: 5000
#' 
#' # Use a different seed for each iteration
#' seeds <- sample.int(1e8, size = niter)
#' 
#' # List where all bootstrap association matrices are stored
#' assoList <- list()
#' 
#' for (i in 1:niter) {
#'   # assoBoot is set to TRUE to return the bootstrap association matrices
#'   net <- netConstruct(amgut1.filt,
#'                       filtTax = "highestFreq",
#'                       filtTaxPar = list(highestFreq = 50),
#'                       filtSamp = "totalReads",
#'                       filtSampPar = list(totalReads = 0),
#'                       measure = "pearson",
#'                       normMethod = "clr",    
#'                       zeroMethod = "pseudoZO",
#'                       sparsMethod = "bootstrap",
#'                       cores = 1,
#'                       nboot = nboot,
#'                       assoBoot = TRUE,
#'                       verbose = 3,
#'                       seed = seeds[i])
#'   
#'   assoList[(1:nboot) + (i - 1) * nboot] <- net$assoBoot1
#' }
#' 
#' # Construct the actual network with all 5000 bootstrap association matrices
#' net_final <- netConstruct(amgut1.filt,
#'                           filtTax = "highestFreq",
#'                           filtTaxPar = list(highestFreq = 50),
#'                           filtSamp = "totalReads",
#'                           filtSampPar = list(totalReads = 0),
#'                           measure = "pearson",
#'                           normMethod = "clr",    
#'                           zeroMethod = "pseudoZO",
#'                           sparsMethod = "bootstrap",
#'                           cores = 1,
#'                           nboot = nboot * niter,
#'                           assoBoot = assoList,
#'                           verbose = 3)
#' 
#' # Network analysis
#' props <- netAnalyze(net_final, clustMethod = "cluster_fast_greedy")
#' 
#' # Network plot
#' plot(props)
#'
#' @seealso \code{\link{netAnalyze}} for analyzing the constructed
#'   network(s), \code{\link{netCompare}} for network comparison,
#'   \code{\link{diffnet}} for constructing differential networks.
#' @references
#'   \insertRef{badri2020normalization}{NetCoMi}\cr\cr
#'   \insertRef{benjamini2000adaptive}{NetCoMi}\cr\cr
#'   \insertRef{farcomeni2007some}{NetCoMi}\cr\cr
#'   \insertRef{friedman2012inferring}{NetCoMi} \cr\cr
#'   \insertRef{WGCNApackage}{NetCoMi}\cr\cr
#'   \insertRef{zhang2005general}{NetCoMi}
#' @importFrom Rdpack reprompt
#' @importFrom vegan vegdist rrarefy
#' @importFrom Matrix nearPD colSums rowSums
#' @importFrom stats var complete.cases pt
#' @importFrom SPRING SPRING mclr
#' @importFrom utils capture.output install.packages installed.packages
#' @importFrom WGCNA pickSoftThreshold TOMdist
#' @import phyloseq
#' @export

netConstruct <- function(data,
                         data2 = NULL,
                         dataType = "counts",
                         group = NULL,
                         matchDesign = NULL,
                         taxRank = NULL,
                         measure = "spieceasi",
                         measurePar = NULL,
                         jointPrepro = NULL,
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
                         assoBoot = NULL,
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
  
  # Check input arguments
  argsIn <- as.list(environment())
  
  # Initialize variables (to pass devtools check)
  assoType <- distNet <- needfrac <- needint <- NULL
  
  if (verbose %in% 2:3) {
    message("Checking input arguments ... ", 
            appendLF = FALSE)
  }
  
  argsOut <- .checkArgsNetConst(argsIn)
  
  for (i in 1:length(argsOut)) {
    assign(names(argsOut)[i], argsOut[[i]])
  }
  
  if (verbose %in% 2:3) message("Done.")
  
  #-----------------------------------------------------------------------------

  if (dataType == "phyloseq") {
    dataType <- "counts"
  }
  
  if (dataType =="counts") {

    # handle phyloseq objects
    
    if (inherits(data, "phyloseq")) {
      
      otutab <- phyloseq::otu_table(data)
      
      if (!is.null(taxRank)) {
        taxtab <- as(tax_table(data), "matrix")
      }
      
      if (attributes(otutab)$taxa_are_rows) {
        data <- t(as(otutab, "matrix"))
      } else {
        data <- as(otutab, "matrix")
      }
      
      if (!is.null(taxRank)) {
        
        if (!taxRank %in% colnames(taxtab)) {
          stop('Argument "taxRank" must match column names of the taxonomic ', 
               'table:\n', paste(colnames(taxtab), collapse = ", "))
        }
        
        if (any(duplicated(taxtab[, taxRank]))) {
          stop(paste0("Taxa names of chosen taxonomic rank must be unique. ",
                      "Consider using NetCoMi's function renameTaxa()."))
        }
        colnames(data) <- taxtab[, taxRank]
      }
      
    } else {
      if (is.null(rownames(data))) {
        message("Row names are numbered because sample names were missing.\n")
        rownames(data) <- 1:nrow(data)
      }
      
      if (is.null(colnames(data))) {
        message("Column names are numbered because taxa names were missing.\n")
        colnames(data) <- 1:ncol(data)
      }
      
      if (identical(colnames(data), rownames(data))) {
        warning(paste0("Row names and column names of 'data' are equal. ",
                       "Ensure 'data' is a count matrix."))
      }
    }
    
    if (!is.null(data2)) {
      if (inherits(data2, "phyloseq")) {
        
        otutab <- phyloseq::otu_table(data2)
        
        if (!is.null(taxRank)) {
          taxtab <- as(tax_table(data2), "matrix")
        }
        
        if (attributes(otutab)$taxa_are_rows) {
          data2 <- t(as(otutab, "matrix"))
        } else {
          data2 <- as(otutab, "matrix")
        }
        
        if (!is.null(taxRank)) {
          if (!taxRank %in% colnames(taxtab)) {
            stop('Argument "taxRank" must match column names of the taxonomic ', 
                 'table:\n', paste(colnames(taxtab), collapse = ", "))
          }
          
          if (any(duplicated(taxtab[, taxRank]))) {
            stop(paste0("Taxa names of chosen taxonomic rank must be unique. ",
                        "Consider using NetCoMi's function renameTaxa()."))
          }
          colnames(data2) <- taxtab[, taxRank]
        }
      }
      
      if (is.null(rownames(data2))) {
        message("Row names of 'data2' are numbered because sample names ", 
                "were missing.\n")
        rownames(data2) <- 1:nrow(data2)
      }
      
      if (is.null(colnames(data2))) {
        message("Column names of 'data2' are numbered because taxa names ", 
                "were missing.\n")
        colnames(data2) <- 1:ncol(data2)
      }
      
      if (identical(colnames(data2), rownames(data2))) {
        warning(paste0("Row names and column names of 'data2' are equal. ",
                       "Ensure 'data2' is a count matrix."))
      }
    }
    
  } else {
    # Catch the case that row and column names are missing
    
    if (xor(is.null(rownames(data)), is.null(colnames(data))) || 
        !identical(rownames(data), colnames(data))) {
      stop("Row and column names must match.")
    }
    
    if (is.null(rownames(data))) {
      message("Row and columns names of 'data' missing. ", 
              "Numbers are used instead.")
      rownames(data) <- colnames(data) <- 1:nrow(data)
    }
    
    
    if (!is.null(data2)) {
      if (xor(is.null(rownames(data2)), is.null(colnames(data2))) || 
          !identical(rownames(data2), colnames(data2))) {
        stop("Row and column names must match.")
      }
      
      if (is.null(rownames(data2))) {
        message("Row and columns names of 'data2' missing. ", 
                "Numbers are used instead.")
        rownames(data2) <- colnames(data2) <- 1:nrow(data2)
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  # check whether arguments are compatible and change if necessary
  
  plausCheck <- .checkPlausNetConst(dataType = dataType, assoType = assoType,
                             data2 = data2, measure = measure,
                             normMethod = normMethod, zeroMethod = zeroMethod,
                             sparsMethod = sparsMethod, dissFunc = dissFunc,
                             sampleSize = sampleSize, verbose = verbose)

  for (i in 1:length(plausCheck)) {
    assign(names(plausCheck)[i], plausCheck[[i]])
  }

  #-----------------------------------------------------------------------------
  # install missing packages
  
  .checkPackNetConst(measure = measure,
                     zeroMethod = zeroMethod,
                     normMethod = normMethod,
                     sparsMethod = sparsMethod,
                     adjust = adjust)
  
  #-----------------------------------------------------------------------------
  # set seed
  
  if (!is.null(seed)) set.seed(seed)
  
  #-----------------------------------------------------------------------------
  # set 'jointPrepro'
  
  # shall two networks be constructed?
  twoNets <- ifelse(is.null(data2) & is.null(group), FALSE, TRUE)
  
  if (twoNets) {
    if (!is.null(group) && !is.null(data2)) {
      stop("Only one of the arguments 'group' and 'data2' may be defined.")
    }
    
    if (!is.null(group)) {
      if (is.null(jointPrepro)) {
        if (distNet) {
          jointPrepro <- FALSE
        } else {
          jointPrepro <- TRUE
        }
      }
    } else if (!is.null(data2) && is.null(jointPrepro)) {
      jointPrepro <- FALSE
    }
    
    if (jointPrepro && distNet) {
      stop("'jointPrepro' is TRUE but data must be normalized separately ", 
           "for dissimilarity measures.")
    }
    
    #--------------------------------
    # set parameters needed for sparsification
    
    if (length(thresh) == 1) thresh <- c(thresh, thresh)
    if (length(alpha) == 1) alpha <- c(alpha, alpha)
    if (length(lfdrThresh) == 1) lfdrThresh <- c(lfdrThresh, lfdrThresh)
    if (length(softThreshPower) == 1) softThreshPower <- c(softThreshPower, 
                                                           softThreshPower)
    if (length(softThreshCut) == 1) softThreshCut <- c(softThreshCut, 
                                                       softThreshCut)
    if (length(sampleSize) == 1) sampleSize <- c(sampleSize, sampleSize)
    
  }
  
  #=============================================================================
  #=============================================================================
  # data preprocessing
  
  # if data contains counts:
  if (dataType == "counts") {
    
    if (!(filtTax[1] == "none" & filtSamp[1] == "none") & verbose %in% 1:3) {
      message("Data filtering ...")
    }
    
    # coerce data to numeric
    countMat1 <-  t(apply(data, 1, function(x) as.numeric(x)))
    colnames(countMat1) <- colnames(data)
    
    if (!is.null(data2)) {
      countMat2 <-  t(apply(data2, 1, function(x) as.numeric(x)))
      colnames(countMat2) <- colnames(data2)
    }
    #---------------------------------------------------------------------------
    
    countMatJoint <- countsJointOrig <- NULL
    
    if (twoNets) {
      if (!is.null(group)) {
        
        # remove NAs
        if (distNet) {
          if (!(is.vector(group) || is.factor(group))) {
            stop("'group' must be of type vector or factor.")
          }
          
          if (length(group) != ncol(countMat1)) {
            stop("Length of 'group' must match the number of columns of 'data'.")
          }
          
          group <- as.numeric(group)
          
          if (is.null(names(group))) {
            names(group) <- colnames(countMat1)
            
          } else {
            if (!all(colnames(countMat1) %in% names(group))) {
              stop("Names of 'group' must match column names of 'data'.")
            }
          }
          
          # remove samples NAs
          if (any(is.na(countMat1))) {
            countMat1 <- countMat1[complete.cases(countMat1), , drop = FALSE]
            if (verbose %in% 1:3) message("Samples with NAs removed.")
          }
          
        } else { # assoNet
          
          if (!(is.vector(group) || is.factor(group))) {
            stop("'group' must be of type vector or factor.")
          }
          
          if (length(group) != nrow(countMat1)) {
            stop("Length of 'group' must match the number of rows of 'data'.")
          }
          
          if (is.character(group)) {
            group <- as.factor(group)
          }
          
          group <- as.numeric(group)
          
          if (is.null(names(group))) {
            names(group) <- rownames(countMat1)
            
          } else {
            if (!all(rownames(countMat1) %in% names(group))) {
              stop("Names of 'group' must match column names of 'data'.")
            }
          }
          
          # remove samples with NAs (group has to be adapted)
          if (any(is.na(countMat1))) {
            if (!is.null(matchDesign)) {
              stop("Data set contains NAs. ", 
                   "Cannot be removed if a matched-group design is used.")
            }
            
            data_tmp <- cbind(countMat1, group)
            data_tmp <- data_tmp[complete.cases(data_tmp), , drop = FALSE]
            countMat1 <- data_tmp[, 1:(ncol(countMat1))]
            group <- data_tmp[, ncol(data_tmp)]
            
            if (verbose %in% 1:3) message("Samples with NAs removed.")
          }
          
          if (!is.null(matchDesign)) {
            if (! (matchDesign[2] / matchDesign[1]) == 
                (table(group)[2] / table(group)[1])) {
              stop("Group vector not consistent with matched-group design.")
            }
          }
          
        }
        
        #-----------------------------------------------------------------------
        # split or join data sets according to 'jointPrepro'
        
        if (jointPrepro) {
          countMatJoint  <- countMat1
          countMat2 <- NULL
          
        } else {
          if (distNet) {
            splitcount <- split(as.data.frame(t(countMat1)), as.factor(group))
            
            if (length(splitcount) != 2)
              stop("Argument 'group' has to be binary.")
            
            groups <- names(splitcount)
            
            countMat1 <- t(as.matrix(splitcount[[1]]))
            countMat2 <- t(as.matrix(splitcount[[2]]))
            
          } else { # assonet
            splitcount <- split(as.data.frame(countMat1), as.factor(group))
            
            if (length(splitcount) != 2)
              stop("Argument 'group' has to be binary.")
            
            groups <- names(splitcount)
            
            countMat1 <- as.matrix(splitcount[[1]])
            countMat2 <- as.matrix(splitcount[[2]])
          }
        }
        
      } else { # data2 given
        groups <- c("1", "2")
        
        if (distNet) {
          if (!identical(colnames(countMat1), colnames(countMat2))) {
            if (!all(rownames(countMat1) %in% rownames(countMat2))) {
              if (verbose > 0) message("Intersection of samples selected.")
            }
            sel <- intersect(rownames(countMat1), rownames(countMat2))
            
            if (length(sel) == 0) stop("Data sets contain different samples")
            
            countMat1 <- countMat1[sel, , drop = FALSE]
            countMat2 <- countMat2[sel, , drop = FALSE]
          }
          
          # remove samples with NAs
          if (any(is.na(countMat1)) || any(is.na(countMat2))) {
            if (verbose %in% 1:3) message("Samples with NAs removed.")
            keep <- intersect(which(complete.cases(countMat1)),
                              which(complete.cases(countMat2)))
            countMat1 <- countMat1[keep, , drop = FALSE]
            countMat2 <- countMat2[keep, , drop = FALSE]
          }
          
        } else { #assoNet
          if (!identical(colnames(countMat1), colnames(countMat2))) {
            if (!all(colnames(countMat1) %in% colnames(countMat2))) {
              if (verbose > 0) message("Intersection of taxa selected.")
            }
            
            sel <- intersect(colnames(countMat1), colnames(countMat2))
            
            if (length(sel) == 0) stop("Data sets contain different taxa.")
            
            countMat1 <- countMat1[ , sel, drop = FALSE]
            countMat2 <- countMat2[ , sel, drop = FALSE]
          }
          
          # remove samples with NAs
          if (any(is.na(countMat1)) || any(is.na(countMat2))) {
            if (!is.null(matchDesign)) {
              stop("Data contain NAs. ", 
                   "Cannot be removed if a matched-group design is used.")
            }
            
            if (verbose %in% 1:3) message("Samples with NAs removed.")
            countMat1 <- countMat1[complete.cases(countMat1), , drop = FALSE]
            countMat2 <- countMat2[complete.cases(countMat2), , drop = FALSE]
          }
          
          if (!is.null(matchDesign)) {
            if (! (matchDesign[2] / matchDesign[1]) == 
                (nrow(countMat2) / nrow(countMat1))) {
              stop("Sample sizes not consistent with matched-group design.")
            }
          }
        }
        
        # join data sets if 'jountPrepro' is TRUE
        if (jointPrepro) {
          
          # Remove duplicates in sample names
          if (any(rownames(countMat2) %in% rownames(countMat1))) {
            if (verbose %in% 1:3) {
              message("\"*\" Added to duplicated sample names in group 2.")
            }
            dupidx <- which(rownames(countMat2) %in% rownames(countMat1))
            rownames(countMat2)[dupidx] <- 
              paste0(rownames(countMat2)[dupidx], "*")
          }
          
          countMatJoint <- rbind(countMat1, countMat2)
          
          n1 <- nrow(countMat1)
          n2 <- nrow(countMat2)
        }
      }
      
      
    } else { # single network
      
      if (any(is.na(data))) {
        if (verbose %in% 1:3) message("Samples with NAs removed.")
        data <- data[complete.cases(data), , drop = FALSE]
      }
      
      countMatJoint <- countMat1
      countMat2 <- NULL
    }
    
    #---------------------------------------------------------------------------
    # sample filtering
    
    if (!twoNets || jointPrepro) {
      keepRows <- .filterSamples(countMat = countMatJoint, filter = filtSamp,
                                 filterParam = filtSampPar)
      
      if (length(keepRows) == 0) {
        stop("No samples remaining after filtering.")
      }
      
      if (length(keepRows) != nrow(countMatJoint)) {
        n_old <- nrow(countMatJoint)
        countMatJoint <- countMatJoint[keepRows, , drop = FALSE]
        if (!distNet) group <- group[keepRows]
        
        if (verbose %in% 2:3) {
          message(n_old - nrow(countMatJoint), " samples removed.")
        }
      }
      
    } else {
      n_old1 <- nrow(countMat1)
      n_old2 <- nrow(countMat2)
      
      keepRows1 <- .filterSamples(countMat = countMat1, filter = filtSamp,
                                  filterParam = filtSampPar)
      
      keepRows2 <- .filterSamples(countMat = countMat2, filter = filtSamp,
                                  filterParam = filtSampPar)
      
      if (distNet) {
        keepRows <- intersect(keepRows1, keepRows2)
        countMat1 <- countMat1[keepRows, , drop = FALSE]
        countMat2 <- countMat2[keepRows, , drop = FALSE]
        
        if (n_old1 - nrow(countMat1) != 0 && verbose %in% 2:3) {
          message(n_old1 - nrow(countMat1),
                  " samples removed in each data set.")
        } 
        
        if (length(keepRows) == 0) {
          stop("No samples remaining after filtering and building the ",
               "intercept of the remaining samples.")
        }
        
      } else {
        countMat1 <- countMat1[keepRows1, , drop = FALSE]
        countMat2 <- countMat2[keepRows2, , drop = FALSE]
        
        if (n_old1 - nrow(countMat1) != 0 || n_old2 - nrow(countMat2)) {
          if (verbose %in% 2:3) {
            message(n_old1 - nrow(countMat1),
                    " samples removed in data set 1.")
            message(n_old2 - nrow(countMat2),
                    " samples removed in data set 2.")
          }
        }
      }
      
      if (length(keepRows1) == 0) {
        stop("No samples remaining in group 1 after filtering.")
      }
      
      if (length(keepRows2) == 0) {
        stop("No samples remaining in group 2 after filtering.")
      }
      
    }
    
    
    #---------------------------------------------------------------------------
    # taxa filtering
    
    if (!twoNets || jointPrepro) {
      keepCols <- .filterTaxa(countMat = countMatJoint, filter = filtTax,
                              filterParam = filtTaxPar)
      
      if (length(keepCols) != ncol(countMatJoint)) {
        p_old <- ncol(countMatJoint)
        countMatJoint <- countMatJoint[, keepCols, drop = FALSE]
        if (verbose %in% 2:3) message(p_old - ncol(countMatJoint),
                                      " taxa removed.")
        if (distNet) group <- group[keepCols]
      }
      
      if (length(keepCols) == 0) {
        stop("No taxa remaining after filtering.")
      }
      
    } else {
      p_old1 <- ncol(countMat1)
      p_old2 <- ncol(countMat2)
      keepCols1 <- .filterTaxa(countMat = countMat1, filter = filtTax,
                               filterParam = filtTaxPar)
      keepCols2 <- .filterTaxa(countMat = countMat2,  filter = filtTax,
                               filterParam = filtTaxPar)
      
      if (!distNet) {
        keepCols <- intersect(keepCols1, keepCols2)
        countMat1 <- countMat1[, keepCols, drop = FALSE]
        countMat2 <- countMat2[, keepCols, drop = FALSE]
        
        if (length(keepCols) == 0) {
          stop("No taxa remaining after filtering and building the intercept ",
               "of the remaining taxa.")
        }
        
        if (p_old1 - dim(countMat1)[2] != 0 && verbose %in% 2:3) {
          message(p_old1 - dim(countMat1)[2],
                  " taxa removed in each data set.")
        }
        
      } else {
        countMat1 <- countMat1[, keepCols1, drop = FALSE]
        countMat2 <- countMat2[, keepCols2, drop = FALSE]
        
        if (p_old1 - dim(countMat1)[2] != 0 || p_old2 - dim(countMat2)[2]!= 0) {
          if (verbose %in% 2:3) {
            message(p_old1 - dim(countMat1)[2],
                    " taxa removed in data set 1.")
            message(p_old2 - dim(countMat2)[2],
                    " taxa removed in data set 2.")
          }
        }
      }
      
      if (length(keepCols1) == 0) {
        stop("No samples remaining in group 1 after filtering.")
      }
      
      if (length(keepCols2) == 0) {
        stop("No samples remaining in group 2 after filtering.")
      }
    }
    
    #---------------------------------------------------------------------------
    # remove samples with zero overall sum
    
    rmZeroSum <- function(countMat, group, matchDesign) {
      rs <- Matrix::rowSums(countMat)
      if (any(rs == 0)) {
        rmRows <- which(rs == 0)
        
        if (!is.null(matchDesign)) {
          stop(paste0("The following samples have an overall sum of zero ",
                      "but cannot be removed if a matched-group design is used: "),
               paste0(rmRows, sep = ","))
        }
        
        countMat <- countMat[-rmRows, , drop = FALSE]
        if (!is.null(group) & !distNet & length(rmRows)!=0) {
          group <- group[-rmRows]
        }
      }
      return(list(countMat = countMat, group = group))
    }
    
    if (!twoNets || jointPrepro) {
      n_old <- nrow(countMatJoint)
      
      rmZeroSum_res <- rmZeroSum(countMatJoint, group = group, 
                                 matchDesign = matchDesign)
      countMatJoint <- rmZeroSum_res$countMat
      group <- rmZeroSum_res$group
      
      if (verbose %in% 2:3 && n_old != nrow(countMatJoint)) {
        message(paste(n_old - nrow(countMatJoint), 
                      "rows with zero sum removed."))
      }
      
    } else {
      n_old1 <- nrow(countMat1)
      n_old2 <- nrow(countMat2)
      
      countMat1 <- rmZeroSum(countMat1, group = group, 
                             matchDesign = matchDesign)$countMat
      
      countMat2 <- rmZeroSum(countMat2, group = group, 
                             matchDesign = matchDesign)$countMat
      
      
      if (distNet) {
        keep <- intersect(rownames(countMat1), rownames(countMat2))
        countMat1 <- countMat1[keep, , drop = FALSE]
        countMat2 <- countMat2[keep, , drop = FALSE]
        
        if (verbose %in% 2:3 && n_old1 != nrow(countMat1)) {
          message(paste(n_old1 - nrow(countMat1), 
                        "rows with zero sum removed in both groups."))
        }
        
      } else {
        if (verbose %in% 2:3) {
          if (n_old1 != nrow(countMat1)) {
            message(paste(n_old1 - nrow(countMat1), 
                          "rows with zero sum removed in group 1."))
          }
          
          if (n_old2 != nrow(countMat2)) {
            message(paste(n_old2 - nrow(countMat2), 
                          "rows with zero sum removed in group 2."))
          }
        }
      }
    }
    
    #---------------------------------------------------------------------------
    # message with remaining samples and taxa
    if (!twoNets || jointPrepro) {
      if (verbose %in% 1:3) message(ncol(countMatJoint), " taxa and ",
                                    nrow(countMatJoint), " samples remaining.")
      
    } else {
      if (verbose %in% 1:3) {
        message(ncol(countMat1), " taxa and ", nrow(countMat1),
                " samples remaining in group 1.")
        message(ncol(countMat2), " taxa and ", nrow(countMat2),
                " samples remaining in group 2.")
      }
    }
    
    #---------------------------------------------------------------------------
    # store original counts and sample sizes
    
    if (twoNets) {
      if (jointPrepro) {
        countsJointOrig <- countMatJoint
        attributes(countsJointOrig)$scale <- "counts"
        
        if (!is.null(data2)) {
          n1 <- sum(rownames(countMatJoint) %in% rownames(countMat1))
          n2 <- sum(rownames(countMatJoint) %in% rownames(countMat2))
        }
        
        countsOrig1 <- countsOrig2 <- NULL
      } else {
        countsOrig1 <- countMat1
        countsOrig2 <- countMat2
        attributes(countsOrig1)$scale <- "counts"
        attributes(countsOrig2)$scale <- "counts"
      }
      
    } else {
      countsOrig1 <- countMatJoint
      attributes(countsOrig1)$scale <- "counts"
      countsOrig2 <- NULL
    }
    
    #---------------------------------------------------------------------------
    # zero treatment
    if (zeroMethod != "none") {
      if (!twoNets || jointPrepro) {
        if (verbose %in% 2:3) message("\nZero treatment:")
        countMatJoint <- .zeroTreat(countMat = countMatJoint, 
                                    zeroMethod = zeroMethod, 
                                    zeroParam = zeroPar, 
                                    needfrac = needfrac, 
                                    needint = needint, 
                                    verbose = verbose)
      } else {
        if (verbose %in% 2:3) message("\nZero treatment in group 1:")
        countMat1 <- .zeroTreat(countMat = countMat1, 
                                zeroMethod = zeroMethod, zeroParam = zeroPar, 
                                needfrac = needfrac, needint = needint, 
                                verbose = verbose)
        
        if (verbose %in% 2:3) message("\nZero treatment in group 2:")
        countMat2 <- .zeroTreat(countMat = countMat2, 
                                zeroMethod = zeroMethod, zeroParam = zeroPar, 
                                needfrac = needfrac, needint = needint, 
                                verbose = verbose)
      }
    } else {
      attributes(countMat1)$scale <- "counts"
      if (!is.null(countMat2)) {
        attributes(countMat2)$scale <- "counts"
      }
      if (!is.null(countMatJoint)) {
        attributes(countMatJoint)$scale <- "counts"
      }
    }
    
    #---------------------------------------------------------------------------
    # normalization
    
    
    if (!twoNets || jointPrepro) {
      if (verbose %in% 2:3 & (normMethod != "none" || needfrac)) {
        message("\nNormalization:")
      }
      countMatJoint <- .normCounts(countMat = countMatJoint, 
                                   normMethod = normMethod,
                                   normParam = normPar, 
                                   zeroMethod = zeroMethod,
                                   needfrac = needfrac, 
                                   verbose = verbose)
      
      if (!twoNets) {
        counts1 <- countMatJoint
        counts2 <- NULL
        
        groups <- NULL
        sampleSize <- nrow(counts1)
      }
      
    } else {
      if (verbose %in% 2:3 & (normMethod != "none" || needfrac)) {
        message("\nNormalization in group 1:")
      }
      countMat1 <- .normCounts(countMat = countMat1, normMethod = normMethod,
                               normParam = normPar, zeroMethod = zeroMethod,
                               needfrac = needfrac, verbose = verbose)
      
      if (verbose %in% 2:3 & (normMethod != "none" || needfrac)) {
        message("\nNormalization in group 2:")
      }
      countMat2 <- .normCounts(countMat = countMat2, normMethod = normMethod,
                               normParam = normPar, zeroMethod = zeroMethod,
                               needfrac = needfrac, verbose = verbose)
    }
    
    #---------------------------------------------------------------------------
    if (twoNets) {
      
      if (jointPrepro) {
        
        if (!is.null(group)) {
          splitcount <- split(as.data.frame(countMatJoint), as.factor(group))
          
          if (length(splitcount) != 2)
            stop("Argument 'group' has to be binary.")
          
          groups <- names(splitcount)
          
          counts1 <- as.matrix(splitcount[[1]])
          counts2 <- as.matrix(splitcount[[2]])
          
          sampleSize <- c(nrow(counts1), nrow(counts2))
          
          rm(countMatJoint, countMat1)
        } else {
          counts1 <- countMatJoint[1:n1, , drop = FALSE]
          counts2 <- countMatJoint[(n1+1):(n1+n2), , drop = FALSE]
          
          sampleSize <- c(nrow(counts1), nrow(counts2))
        }
      } else {
        counts1 <- countMat1
        counts2 <- countMat2
        
        if (distNet) {
          keep <- intersect(rownames(counts1), rownames(counts2))
          counts1 <- counts1[keep, , drop = FALSE]
          counts2 <- counts2[keep, , drop = FALSE]
        }else {
          keep <- intersect(colnames(counts1), colnames(counts2))
          counts1 <- counts1[ , keep, drop = FALSE]
          counts2 <- counts2[ , keep, drop = FALSE]
        }
      }
      
    } else { # single network
      rm(countMat1, countMatJoint)
    }
    
    #===========================================================================
    #===========================================================================
    # association / dissimilarity estimation
    
    if (verbose %in% 2:3) {
      if (distNet) {
        txt.tmp <- "dissimilarities"
      } else {
        txt.tmp <- "associations"
      }
      message("\nCalculate '", measure, "' ", txt.tmp, " ... ", 
              appendLF = FALSE)
    }
    
    assoMat1 <- .calcAssociation(countMat = counts1, measure = measure,
                                 measurePar = measurePar, verbose = verbose)
    
    if (verbose %in% 2:3) message("Done.")
    
    
    if (twoNets) {
      
      if (verbose %in% 2:3) {
        message("\nCalculate ", txt.tmp, " in group 2 ... ",
                appendLF = FALSE)
      }
      
      assoMat2 <- .calcAssociation(countMat = counts2, measure = measure,
                                   measurePar = measurePar, verbose = verbose)
      if (verbose %in% 2:3) message("Done.")
      
    } else {
      assoMat2 <- NULL
    }
    
  } else {
    assoMat1 <- data
    assoMat2 <- data2
    counts1 <- NULL
    counts2 <- NULL
    countsJointOrig <- NULL
    countsOrig1 <- NULL
    countsOrig2 <- NULL
    groups <- NULL
  }
  
  if (distNet) {
    
    dissEst1 <- assoMat1
    
    if (any(is.infinite(dissEst1)) & scaleDiss) {
      scaleDiss <- FALSE
      warning("Dissimilarity matrix contains infinite values and cannot be ", 
              "scaled to [0,1].")
    }
    
    if (scaleDiss) {
      assoUpper <- assoMat1[upper.tri(assoMat1)]
      
      if (length(assoUpper) == 1) {
        warning("Network consists of only two nodes.")
        assoMat1[1,2] <- assoMat1[2,1] <- 1
        
      } else {
        assoMat1 <- (assoMat1 - min(assoUpper)) / 
          (max(assoUpper) - min(assoUpper))
        
        diag(assoMat1) <- 0
      }
    }
    
    dissScale1 <- assoMat1
    
    if (verbose %in% 2:3) {
      if (sparsMethod != "none") {
        message("\nSparsify dissimilarities via '", sparsMethod, "' ... ",
                appendLF = FALSE)
      }
    }
    
    sparsReslt <- .sparsify(assoMat = assoMat1,
                            countMat = counts1,
                            sampleSize = sampleSize[1],
                            measure = measure,
                            measurePar = measurePar,
                            assoType = assoType,
                            sparsMethod = sparsMethod,
                            thresh = thresh[1],
                            alpha = alpha[1],
                            adjust = adjust,
                            lfdrThresh = lfdrThresh[1],
                            trueNullMethod = trueNullMethod,
                            nboot = nboot,
                            assoBoot = assoBoot,
                            softThreshType = softThreshType,
                            softThreshPower = softThreshPower[1],
                            softThreshCut = softThreshCut[1],
                            cores = cores,
                            logFile = logFile,
                            kNeighbor = kNeighbor,
                            knnMutual = knnMutual,
                            verbose = verbose,
                            seed = seed)
    
    if (verbose %in% 2:3 & sparsMethod != "none") message("Done.")
    
    assoMat1 <- NULL
    assoEst1 <- NULL
    dissMat1 <- sparsReslt$assoNew
    power1 <- sparsReslt$power
    simMat1 <- .transToSim(x = dissMat1, simFunc = simFunc, 
                            simFuncPar = simFuncPar)
    adjaMat1 <- .transToAdja(x = simMat1, weighted = weighted)
    
    if (twoNets) {
      dissEst2 <- assoMat2
      
      if (scaleDiss) {
        assoUpper <- assoMat2[upper.tri(assoMat2)]
        
        if (length(assoUpper) == 1) {
          warning("Network consists of only two nodes.")
          assoMat2[1,2] <- assoMat2[2,1] <- 1
          
        } else {
          assoMat2 <- (assoMat2 - min(assoUpper)) / 
            (max(assoUpper) - min(assoUpper))
          
          diag(assoMat2) <- 0
        }
        
      }
      dissScale2 <- assoMat2
      
      if (verbose %in% 2:3) {
        if (sparsMethod != "none") {
          message("\nSparsify dissimilarities in group 2 ... ",
                  appendLF = FALSE)
        }
      }
      
      sparsReslt <- .sparsify(assoMat = assoMat2,
                              countMat = counts2,
                              sampleSize = sampleSize[2],
                              measure = measure,
                              measurePar = measurePar,
                              assoType = assoType,
                              sparsMethod = sparsMethod,
                              thresh = thresh[2],
                              alpha = alpha[2],
                              adjust = adjust,
                              lfdrThresh = lfdrThresh[2],
                              trueNullMethod = trueNullMethod,
                              nboot = nboot,
                              assoBoot = assoBoot,
                              softThreshType = softThreshType,
                              softThreshPower = softThreshPower[2],
                              softThreshCut = softThreshCut[2],
                              cores = cores,
                              logFile = logFile,
                              kNeighbor = kNeighbor,
                              knnMutual = knnMutual,
                              verbose = verbose,
                              seed = seed)
      
      if (verbose %in% 2:3 & sparsMethod != "none") message("Done.")
      
      assoMat2 <- NULL
      assoEst2 <- NULL
      dissMat2 <- sparsReslt$assoNew
      power2 <- sparsReslt$power
      simMat2 <- .transToSim(x = dissMat2, simFunc = simFunc,
                              simFuncPar = simFuncPar)
      adjaMat2 <- .transToAdja(x = simMat2, weighted = weighted)
      
    } else {
      dissEst2 <- dissScale2 <- assoMat2 <- assoEst2 <- dissMat2 <-
        power2 <- simMat2 <- adjaMat2 <- NULL
    }
    
    assoBoot1 <- assoBoot2 <- NULL
    
  } else { # association network
    if (verbose %in% 2:3) {
      if (sparsMethod != "none") {
        message("\nSparsify associations via '", sparsMethod, "' ... ",
                appendLF = FALSE)
      }
    }
    
    sparsReslt <- .sparsify(assoMat = assoMat1,
                            countMat = counts1,
                            sampleSize = sampleSize[1],
                            measure = measure,
                            measurePar = measurePar,
                            assoType = assoType,
                            sparsMethod = sparsMethod,
                            thresh = thresh[1],
                            alpha = alpha[1],
                            adjust = adjust,
                            lfdrThresh = lfdrThresh[1],
                            trueNullMethod = trueNullMethod,
                            nboot = nboot,
                            assoBoot = assoBoot,
                            softThreshType = softThreshType,
                            softThreshPower = softThreshPower[1],
                            softThreshCut = softThreshCut[1],
                            cores = cores,
                            logFile = logFile,
                            kNeighbor = kNeighbor,
                            knnMutual = knnMutual,
                            verbose = verbose,
                            seed = seed)
    
    if (verbose %in% 2:3 & sparsMethod != "none") message("Done.")
    assoEst1 <- assoMat1
    assoMat1 <- sparsReslt$assoNew
    power1 <- sparsReslt$power
    dissEst1 <- dissScale1 <- NULL
    
    dissMat1 <- .transToDiss(x = assoMat1, dissFunc = dissFunc,
                              dissFuncPar = dissFuncPar)
    
    if (sparsMethod == "softThreshold") {
      simMat1 <- sparsReslt$simMat
      adjaMat1 <- assoMat1
      assoBoot1 <- NULL
      
    } else {
      simMat1 <- .transToSim(x = dissMat1, simFunc = simFunc,
                              simFuncPar = simFuncPar)
      
      adjaMat1 <- .transToAdja(x = simMat1, weighted = weighted)
      
      if (sparsMethod == "bootstrap" && !is.null(assoBoot) && 
          !is.list(assoBoot) && assoBoot == TRUE) {
        assoBoot1 <- sparsReslt$assoBoot
      } else {
        assoBoot1 <- NULL
      }
    }
    
    if (twoNets) {
      if (verbose %in% 2:3) {
        if (sparsMethod != "none") {
          message("\nSparsify associations in group 2 ... ",
                  appendLF = FALSE)
        }
      }
      
      sparsReslt <- .sparsify(assoMat = assoMat2,
                              countMat = counts2,
                              sampleSize = sampleSize[2],
                              measure = measure,
                              measurePar = measurePar,
                              assoType = assoType,
                              sparsMethod = sparsMethod,
                              thresh = thresh[2],
                              alpha = alpha[2],
                              adjust = adjust,
                              lfdrThresh = lfdrThresh[2],
                              trueNullMethod = trueNullMethod,
                              nboot = nboot,
                              assoBoot = assoBoot,
                              softThreshType = softThreshType,
                              softThreshPower = softThreshPower[2],
                              softThreshCut = softThreshCut[2],
                              cores = cores,
                              logFile = logFile,
                              kNeighbor = kNeighbor,
                              knnMutual = knnMutual,
                              verbose = verbose,
                              seed = seed)
      
      if (verbose %in% 2:3 & sparsMethod != "none") message("Done.")
      assoEst2 <- assoMat2
      assoMat2 <- sparsReslt$assoNew
      power2 <- sparsReslt$power
      dissEst2 <- dissScale2 <- NULL
      
      dissMat2 <- .transToDiss(x = assoMat2, dissFunc = dissFunc,
                                dissFuncPar = dissFuncPar)
      
      if (sparsMethod == "softThreshold") {
        simMat2 <- sparsReslt$simMat
        adjaMat2 <- assoMat2
        assoBoot2 <- NULL
        
      } else {
        simMat2 <- .transToSim(x = dissMat2, simFunc = simFunc,
                                simFuncPar = simFuncPar)
        
        adjaMat2 <- .transToAdja(x = simMat2, weighted = weighted)
        
        if (sparsMethod == "bootstrap" && !is.null(assoBoot) && 
            !is.list(assoBoot) && assoBoot == TRUE) {
          assoBoot2 <- sparsReslt$assoBoot
        } else {
          assoBoot2 <- NULL
        }
      }
      
    } else {
      dissEst2 <- dissScale2 <- assoMat2 <- assoEst2 <- dissMat2 <-
        power2 <- simMat2 <- adjaMat2 <- assoBoot2 <- NULL
    }
    
  }
  
  #=============================================================================
  # Create edge list
  g <- graph_from_adjacency_matrix(adjaMat1, weighted = TRUE,
                                   mode = "undirected", diag = FALSE)
  
  if (is.null(E(g)$weight)) {
    isempty1 <- TRUE
    edgelist1 <- NULL
    
  } else {
    isempty1 <- FALSE
    
    edgelist1 <- data.frame(get.edgelist(g))
    colnames(edgelist1) <- c("v1", "v2")
    
    if (!is.null(assoMat1)) {
      edgelist1$asso <- sapply(1:nrow(edgelist1), function(i) {
        assoMat1[edgelist1[i, 1], edgelist1[i, 2]]
      })
    }
    
    edgelist1$diss <- sapply(1:nrow(edgelist1), function(i) {
      dissMat1[edgelist1[i, 1], edgelist1[i, 2]]
    })
    
    if (all(adjaMat1 %in% c(0,1))) {
      edgelist1$sim <-sapply(1:nrow(edgelist1), function(i) {
        simMat1[edgelist1[i, 1], edgelist1[i, 2]]
      })
    }
    
    edgelist1$adja <- sapply(1:nrow(edgelist1), function(i) {
      adjaMat1[edgelist1[i, 1], edgelist1[i, 2]]
    })
  }
  
  if (twoNets) {
    # Create edge list
    g <- graph_from_adjacency_matrix(adjaMat2, weighted = TRUE, 
                                     mode = "undirected", diag = FALSE)
    
    if (is.null(E(g)$weight)) {
      isempty2 <- TRUE
      edgelist2 <- NULL
      
    } else {
      isempty2 <- FALSE
      
      edgelist2 <- data.frame(get.edgelist(g))
      colnames(edgelist2) <- c("v1", "v2")
      
      if (!is.null(assoMat2)) {
        edgelist2$asso <- sapply(1:nrow(edgelist2), function(i) {
          assoMat2[edgelist2[i, 1], edgelist2[i, 2]]
        })
      }
      
      edgelist2$diss <- sapply(1:nrow(edgelist2), function(i) {
        dissMat2[edgelist2[i, 1], edgelist2[i, 2]]
      })
      
      if (all(adjaMat2 %in% c(0, 1))) {
        edgelist2$sim <- sapply(1:nrow(edgelist2), function(i) {
          simMat2[edgelist2[i, 1], edgelist2[i, 2]]
        })
      }
      
      edgelist2$adja <- sapply(1:nrow(edgelist2), function(i) {
        adjaMat2[edgelist2[i, 1], edgelist2[i, 2]]
      })
    }
    
    if (isempty1 && verbose > 0) {
      message("\nNetwork 1 has no edges.")
    }
    if (isempty2 && verbose > 0) {
      message("Network 2 has no edges.")
    }
  } else {
    edgelist2 <- NULL
    if (isempty1 && verbose > 0) {
      message("\nNetwork has no edges.")
    }
  }
  
  #=============================================================================
  output <- list()
  
  output$edgelist1 <- edgelist1
  output$edgelist2 <- edgelist2
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
  
  output$assoBoot1 <- assoBoot1
  output$assoBoot2 <- assoBoot2
  
  output$countMat1 <- countsOrig1
  output$countMat2 <- countsOrig2
  if (!is.null(countsJointOrig)) output$countsJoint <- countsJointOrig
  output$normCounts1 <- counts1
  output$normCounts2 <- counts2
  output$groups <- groups
  output$matchDesign <- matchDesign
  output$sampleSize <- sampleSize
  output$softThreshPower <- list(power1 = power1, 
                                 power2 = power2) # calculated power
  output$assoType <- assoType
  output$twoNets <- twoNets
  
  output$parameters <- list(
    dataType = dataType,
    group = group,
    filtTax = filtTax,
    filtTaxPar = filtTaxPar,
    filtSamp = filtSamp,
    filtSampPar = filtSampPar,
    jointPrepro = jointPrepro,
    zeroMethod = zeroMethod,
    zeroPar = zeroPar,
    needfrac = needfrac,
    needint = needint,
    normMethod = normMethod,
    normPar = normPar,
    measure = measure,
    measurePar = measurePar,
    sparsMethod = sparsMethod,
    thresh = thresh,
    adjust = adjust,
    trueNullMethod = trueNullMethod,
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
  
  output$call = match.call()
  
  class(output) <- "microNet"
  return(output)
  
}
