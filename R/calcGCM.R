#' @title Graphlet Correlation Matrix (GCM)
#' 
#' @description Computes the Graphlet Correlation Matrix (GCM) of a network, 
#'   given as adjacency matrix.
#'   
#'   The GCM of a network is a matrix with Spearman's correlations between 
#'   the network's node orbits (Hocevar and Demsar, 
#'   2016; Yaveroglu et al., 2014). 
#'   
#'   The function considers only orbits for graphlets with up to four nodes. 
#'   Orbit counts are determined using the function \code{\link[orca]{count4}} 
#'   from \code{orca} package.
#'   
#'   Unobserved orbits would lead to NAs in the correlation matrix, which is 
#'   why a row with pseudo counts of 1 is added to the orbit count matrix 
#'   (\code{ocount}).
#'   
#'   The function is based on R code provided by Theresa Ullmann 
#'   (\url{https://orcid.org/0000-0003-1215-8561}).
#'   
#' @param adja adjacency matrix (numeric) defining the network
#'   for which the GCM should be calculated.
#' @param orbits numeric vector with integers from 0 to 14 defining the graphlet 
#'   orbits to use for GCM calculation. Minimum length is 2. Defaults to 
#'   c(0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11), thus excluding redundant orbits such
#'   as the orbit o3.  
#'   
#' @details 
#'   By default, only the 11 non-redundant orbits are used. These are grouped 
#'   according to their role: orbit 0 represents the degree, orbits (2, 5, 7) 
#'   represent nodes within a chain, orbits (8, 10, 11) represent nodes in a 
#'   cycle, and orbits (6, 9, 4, 1) represent a terminal node.
#'   
#' @return A list with the following elements:
#'   \tabular{ll}{
#'   \code{gcm}\tab Graphlet Correlation Matrix\cr
#'   \code{ocount}\tab Orbit counts
#'   }
#' 
#' @examples 
#' # Load data set from American Gut Project (from SpiecEasi package)
#' data("amgut1.filt")
#' 
#' # Network construction
#' net <- netConstruct(amgut1.filt, 
#'                     filtTax = "highestFreq",
#'                     filtTaxPar = list(highestFreq = 50),
#'                     measure = "pearson",
#'                     normMethod = "clr",
#'                     zeroMethod = "pseudoZO",
#'                     sparsMethod = "thresh",
#'                     thresh = 0.5)
#' 
#' # Get adjacency matrices
#' adja <- net$adjaMat1
#' 
#' # Network visualization
#' props <- netAnalyze(net)
#' plot(props, rmSingles = TRUE, cexLabels = 1.7)
#' 
#' # Calculate Graphlet Correlation Matrix (GCM)
#' gcm <- calcGCM(adja)
#' 
#' gcm
#' 
#' # Plot heatmap of the GCM
#' plotHeat(gcm$gcm)
#'   
#' @references
#'   \insertRef{hocevar2016computation}{NetCoMi}\cr\cr
#'   \insertRef{yaveroglu2014revealing}{NetCoMi}
#'   
#' @seealso \code{\link{calcGCD}}, \code{\link{testGCM}}
#'   
#' @importFrom orca count4
#' 
#' @export

calcGCM <- function(adja, 
                    orbits = c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1)) {
  
  # Install missing packages
  if (!"orca" %in% utils::installed.packages()[,"Package"]) {
    message("Installing missing package: orca\n")
    install.packages("orca", dependencies = TRUE)
  }
  
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
  
  if (!nrow(adja) == ncol(adja)) {
    stop("Numbers of rows and columns of 'adja' differ.", 
         "'adja' must be an adjacency matrix.")
  }
  
  # Get orbit counts
  ocount <- .getOrbcounts(adja = adja, orbits = orbits)
  
  # A pseudo count of 1 is added to all columns
  # (to avoid NAs in the correlation matrix due to unobserved orbits)
  pseudo <- rep(1, ncol(ocount))
  
  # Pseudo count of 2 is used if all entries of an orb are 1 (no variation)
  for (i in 1:ncol(ocount)) {
    if (all(ocount[, i] == 1)) {
      pseudo[i] <- 2
    }
  }
  
  ocount <- rbind(ocount, pseudo)
  
  # Compute Graphlet Correlation Matrix
  gcm <- suppressWarnings(cor(ocount, method = "spearman"))
  
  out <- list(gcm = gcm, ocount = ocount)
  
  class(out) <- "GCM"
  
  return(out)
}


#' @keywords internal
.getOrbcounts <- function(adja, orbits = 0:14) {
  
  if (all(adja[lower.tri(adja)] == 0)) {
    # Handle adjacency with only zero entries
    ocounts <- matrix(0, nrow = 2, ncol = 15)
    colnames(ocounts) <- paste0("O", 0:14)
    
  } else {
    net <- igraph::graph_from_adjacency_matrix(adja, weighted=T,
                                               mode="undirected", diag=F)
    
    edgelist <- igraph::get.edgelist(net, names = FALSE)
    edgelist <- apply(edgelist, 2, as, Class = "integer")
    
    # Transform into matrix if network consists of one edge only
    
    if (is.vector(edgelist)) {
      edgelist <- t(as.matrix(edgelist))
    }
    
    # Get orbit counts
    ocounts <- orca::count4(edgelist)
  }

  # Select orbits of interest
  orbits <- paste0("O", orbits)
  ocounts <- ocounts[, orbits]
  
  # Fill up count matrix so that nrow equals number of nodes
  nnodes <- ncol(adja)
  
  if (nrow(ocounts) < nnodes) {
    ocounts <- rbind(ocounts, matrix(0, nrow = nnodes-nrow(ocounts), 
                                     ncol = ncol(ocounts)))
  }
  
  rownames(ocounts) <- rownames(adja)
  
  return(ocounts)
}

