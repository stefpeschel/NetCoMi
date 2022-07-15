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
#' @param orbits numeric vector with integers from 0 to 14 defining the orbits 
#'   of interest. Minimum length is 2. Defaults to 
#'   c(0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11), thus excluding redundant orbits such
#'   as the orbit o3.  
#' 
#' @return An object of class \code{gcd} containing the following elements:
#'   \tabular{ll}{
#'   \code{gcd}\tab Graphlet Correlation Distance between the two networks\cr
#'   \code{ocount1, ocount2}\tab Orbit counts \cr
#'   \code{gcm1, gcm2}\tab Graphlet 
#'   Correlation Matrices
#'   }
#' 
#' @examples 
#' # Load data set from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#' 
#' # Generate a random group vector
#' set.seed(123456)
#' group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)
#' 
#' # Network construction
#' net <- netConstruct(amgut1.filt, 
#'                     group = group,
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
#' # Calculate GCD
#' gcdlist <- calcGCD(adja1, adja2)
#' 
#' gcdlist
#' 
#' # Orbit counts
#' gcdlist$ocount1
#' gcdlist$ocount2
#' 
#' # GCMs
#' gcdlist$gcm1
#' gcdlist$gcm2
#'   
#' @references
#'   \insertRef{hocevar2016computation}{NetCoMi}\cr\cr
#'   \insertRef{yaveroglu2014revealing}{NetCoMi}
#'   
#' @import orca
#' 
#' @export

calcGCD <- function(adja1, adja2, orbits = c(0:2, 4:11)) {
  
  if (!is.numeric(orbits)) {
    stop("'orbits' vector must be numeric.")
  }
  
  if (length(orbits) < 2) {
    stop("At least two orbits must be selected to compute the GCD.")
  }
  
  if (any(orbits < 0) | any(orbits > 14) | length(orbits) > 15) {
    stop("Only orbits 0 to 14 (from 4-node graphlets) are allowed.")
  }
  
  if (!all(orbits %% 1 == 0)) {
    stop("Elements of 'orbits' must be whole numbers.")
  }
  
  if (!nrow(adja1) == ncol(adja1)) {
    stop("Numbers of rows and columns of 'adja1' differ.", 
         "'adja1' must be an adjacency matrix.")
  }
  
  if (!nrow(adja2) == ncol(adja2)) {
    stop("Numbers of rows and columns of 'adja2' differ.", 
         "'adja2' must be an adjacency matrix.")
  }
  
  # Compute Graphlet Correlation Matrices
  gcm1.list <- calcGCM(adja1)
  gcm2.list <- calcGCM(adja2)
  
  gcm1 <- gcm1.list$gcm
  gcm2 <- gcm2.list$gcm
  
  ocount1 <- gcm1.list$ocount
  ocount2 <- gcm2.list$ocount
  
  # GC vectors
  corvec1 <- gcm1[upper.tri(gcm1)]
  corvec2 <- gcm2[upper.tri(gcm2)]
  
  # Compute Graphlet Correlation Distance
  gcd <- sqrt(sum((corvec1 - corvec2)^2))
  
  out <- list(gcd = gcd, 
              ocount1 = ocount1, 
              ocount2 = ocount2, 
              gcm1 = gcm1, gcm2 = gcm2)
  
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

