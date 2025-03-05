#' @title Graphlet Correlation Distance (GCD)
#' 
#' @description Computes the Graphlet Correlation Distance (GCD) - a 
#'   graphlet-based distance measure - between two networks.
#'   
#'   Following Yaveroglu et al. (2014), the GCD is defined as
#'   the Euclidean distance of the upper triangle values of the Graphlet 
#'   Correlation Matrices (GCM) of two networks, which are defined by their 
#'   adjacency matrices. 
#'   The GCM of a network is a matrix with Spearman's correlations between 
#'   the network's node orbits (Hocevar and Demsar, 2016). 
#'   
#'   The function considers only orbits for graphlets with up to four nodes. 
#'   Orbit counts are determined using the function \code{\link[orca]{count4}} 
#'   from \code{orca} package.
#'   
#'   Unobserved orbits would lead to NAs in the correlation matrix, which is 
#'   why a row with pseudo counts of 1 is added to the orbit count matrices 
#'   (\code{ocount1} and \code{ocount2}).
#'   
#'   The function is based on R code provided by Theresa Ullmann 
#'   (\url{https://orcid.org/0000-0003-1215-8561}).
#'   
#' @param adja1,adja2 adjacency matrices (numeric) defining the two networks 
#'   between which the GCD shall be calculated.
#' @param orbits numeric vector with integers from 0 to 14 defining the graphlet 
#'   orbits to use for GCD calculation. Minimum length is 2. Defaults to 
#'   c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1), thus excluding redundant orbits such
#'   as the orbit o3. See details.
#'   
#' @details 
#'   By default, only the 11 non-redundant orbits are used. These are grouped 
#'   according to their role: orbit 0 represents the degree, orbits (2, 5, 7) 
#'   represent nodes within a chain, orbits (8, 10, 11) represent nodes in a 
#'   cycle, and orbits (6, 9, 4, 1) represent a terminal node.
#' 
#' @return An object of class \code{gcd} containing the following elements:
#'   \tabular{ll}{
#'   \code{gcd}\tab Graphlet Correlation Distance between the two networks\cr
#'   \code{ocount1, ocount2}\tab Orbit counts \cr
#'   \code{gcm1, gcm2}\tab Graphlet Correlation Matrices
#'   }
#' 
#' @examples 
#' library(phyloseq)
#' 
#' # Load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut2.filt.phy")
#' 
#' # Split data into two groups: with and without seasonal allergies
#' amgut_season_yes <- phyloseq::subset_samples(amgut2.filt.phy, 
#'                                       SEASONAL_ALLERGIES == "yes")
#' amgut_season_no <- phyloseq::subset_samples(amgut2.filt.phy, 
#'                                      SEASONAL_ALLERGIES == "no")
#' 
#' # Make sample sizes equal to ensure comparability
#' n_yes <- phyloseq::nsamples(amgut_season_yes)
#' ids_yes <- phyloseq::get_variable(amgut_season_no, "X.SampleID")[1:n_yes]
#' 
#' amgut_season_no <- phyloseq::subset_samples(amgut_season_no, X.SampleID %in% ids_yes)
#' 
#' # Network construction
#' net <- netConstruct(amgut_season_yes,
#'                     amgut_season_no, 
#'                     filtTax = "highestFreq",
#'                     filtTaxPar = list(highestFreq = 50),
#'                     measure = "pearson",
#'                     normMethod = "clr",
#'                     zeroMethod = "pseudoZO",
#'                     sparsMethod = "thresh",
#'                     thresh = 0.5)
#' 
#' # Get adjacency matrices
#' adja1 <- net$adjaMat1
#' adja2 <- net$adjaMat2
#' 
#' # Network visualization
#' props <- netAnalyze(net)
#' plot(props, rmSingles = TRUE, cexLabels = 1.7)
#' 
#' # Calculate the GCD
#' gcd <- calcGCD(adja1, adja2)
#' 
#' gcd
#' 
#' # Orbit counts
#' head(gcd$ocount1)
#' head(gcd$ocount2)
#' 
#' # GCMs
#' gcd$gcm1
#' gcd$gcm2
#' 
#' # Test Graphlet Correlations for significant differences
#' gcmtest <- testGCM(gcd)
#' 
#' ### Plot heatmaps
#' # GCM 1 (with significance code in the lower triangle)
#' plotHeat(gcmtest$gcm1, 
#'          pmat = gcmtest$pAdjust1,
#'          type = "mixed")
#' 
#' # GCM 2 (with significance code in the lower triangle)
#' plotHeat(gcmtest$gcm2, 
#'          pmat = gcmtest$pAdjust2,
#'          type = "mixed")
#' 
#' # Difference GCM1-GCM2 (with p-values in the lower triangle)
#' plotHeat(gcmtest$diff, 
#'          pmat = gcmtest$pAdjustDiff,
#'          type = "mixed",
#'          textLow = "pmat")
#'   
#' @references
#'   \insertRef{hocevar2016computation}{NetCoMi}\cr\cr
#'   \insertRef{yaveroglu2014revealing}{NetCoMi}
#'   
#' @seealso \code{\link{calcGCM}}, \code{\link{testGCM}}
#' 
#' @export

calcGCD <- function(adja1, adja2,
                    orbits = c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)) {
  
  if (!is.numeric(orbits)) {
    stop("\"orbits\" vector must be numeric.")
  }
  
  if (length(orbits) < 2) {
    stop("At least two orbits must be selected to compute the GCD.")
  }
  
  if (any(orbits < 0) | any(orbits > 14) | length(orbits) > 15) {
    stop("Only orbits 0 to 14 (from 4-node graphlets) are allowed.")
  }
  
  if (!all(orbits %% 1 == 0)) {
    stop("Elements of \"orbits\" must be whole numbers.")
  }
  
  if (!nrow(adja1) == ncol(adja1)) {
    stop("Numbers of rows and columns of \"adja1\" differ.", 
         "\"adja1\" must be an adjacency matrix.")
  }
  
  if (!nrow(adja2) == ncol(adja2)) {
    stop("Numbers of rows and columns of \"adja2\" differ.", 
         "\"adja2\" must be an adjacency matrix.")
  }
  
  #=============================================================================

  # Compute Graphlet Correlation Matrices
  gcm1List <- calcGCM(adja1)
  gcm2List <- calcGCM(adja2)
  
  gcm1 <- gcm1List$gcm
  gcm2 <- gcm2List$gcm
  
  ocount1 <- gcm1List$ocount
  ocount2 <- gcm2List$ocount
  
  # GC vectors
  corvec1 <- gcm1[upper.tri(gcm1)]
  corvec2 <- gcm2[upper.tri(gcm2)]
  
  # Compute Graphlet Correlation Distance
  gcd <- sqrt(sum((corvec1 - corvec2)^2))

  out <- list(gcd = gcd, 
              ocount1 = ocount1, 
              ocount2 = ocount2, 
              gcm1 = gcm1, 
              gcm2 = gcm2)

  class(out) <- "GCD"
  
  return(out)
}


#' @title Print method for GCD objects
#' 
#' @param x object of class \code{GCD} (returned by to \code{\link{calcGCD}}).
#' @param ... not used.
#' 
#' @method print GCD
#' @rdname print.GCD
#' @export
print.GCD <- function(x, ...) {
  cat("GCD: ", round(x$gcd, 5), "\n")
}

