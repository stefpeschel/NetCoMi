% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcGCM.R
\name{calcGCM}
\alias{calcGCM}
\title{Graphlet Correlation Matrix (GCM)}
\usage{
calcGCM(adja, orbits = c(0, 2, 5, 7, 8, 10, 11, 6, 9, 4, 1))
}
\arguments{
\item{adja}{adjacency matrix (numeric) defining the network
for which the GCM should be calculated.}

\item{orbits}{numeric vector with integers from 0 to 14 defining the graphlet 
orbits to use for GCM calculation. Minimum length is 2. Defaults to 
c(0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11), thus excluding redundant orbits such
as the orbit o3.}
}
\value{
A list with the following elements:
  \tabular{ll}{
  \code{gcm}\tab Graphlet Correlation Matrix\cr
  \code{ocount}\tab Orbit counts
  }
}
\description{
Computes the Graphlet Correlation Matrix (GCM) of a network, 
  given as adjacency matrix.
  
  The GCM of a network is a matrix with Spearman's correlations between 
  the network's node orbits (Hocevar and Demsar, 
  2016; Yaveroglu et al., 2014). 
  
  The function considers only orbits for graphlets with up to four nodes. 
  Orbit counts are determined using the function \code{\link[orca]{count4}} 
  from \code{orca} package.
  
  Unobserved orbits would lead to NAs in the correlation matrix, which is 
  why a row with pseudo counts of 1 is added to the orbit count matrix 
  (\code{ocount}).
  
  The function is based on R code provided by Theresa Ullmann 
  (\url{https://orcid.org/0000-0003-1215-8561}).
}
\details{
By default, only the 11 non-redundant orbits are used. These are grouped 
  according to their role: orbit 0 represents the degree, orbits (2, 5, 7) 
  represent nodes within a chain, orbits (8, 10, 11) represent nodes in a 
  cycle, and orbits (6, 9, 4, 1) represent a terminal node.
}
\examples{
# Load data set from American Gut Project (from SpiecEasi package)
data("amgut1.filt")

# Network construction
net <- netConstruct(amgut1.filt, 
                    filtTax = "highestFreq",
                    filtTaxPar = list(highestFreq = 50),
                    measure = "pearson",
                    normMethod = "clr",
                    zeroMethod = "pseudoZO",
                    sparsMethod = "thresh",
                    thresh = 0.5)

# Get adjacency matrices
adja <- net$adjaMat1

# Network visualization
props <- netAnalyze(net)
plot(props, rmSingles = TRUE, cexLabels = 1.7)

# Calculate Graphlet Correlation Matrix (GCM)
gcm <- calcGCM(adja)

gcm

# Plot heatmap of the GCM
plotHeat(gcm$gcm)
  
}
\references{
\insertRef{hocevar2016computation}{NetCoMi}\cr\cr
  \insertRef{yaveroglu2014revealing}{NetCoMi}
}
\seealso{
\code{\link{calcGCD}}, \code{\link{testGCM}}
}
