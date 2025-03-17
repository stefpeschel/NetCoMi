knitr::opts_chunk$set(fig.width = 8, fig.height = 8)

# Load data sets from American Gut Project (from SpiecEasi package)
data("amgut1.filt")
data("amgut2.filt.phy")

# Single network with the following specifications:
# - Association measure: SpiecEasi
# - SpiecEasi parameters are defined via 'measurePar'
#   (check ?SpiecEasi::spiec.easi for available options)
# - Note: 'rep.num' should be higher for real data sets
# - Taxa filtering: Keep the 50 taxa with highest variance
# - Sample filtering: Keep samples with a total number of reads
#   of at least 1000

net1 <- netConstruct(amgut2.filt.phy,
                     measure = "spieceasi",
                     measurePar = list(method = "mb",
                                       pulsar.params = list(rep.num = 10),
                                       symBetaMode = "ave"),
                     filtTax = "highestVar",
                     filtTaxPar = list(highestVar = 50),
                     filtSamp = "totalReads",
                     filtSampPar = list(totalReads = 1000),
                     sparsMethod = "none",
                     normMethod = "none",
                     verbose = 3)
# Output returned by spiec.easi()
spiec_output <- net1$measureOut1

# Network analysis (see ?netAnalyze for details)
props1 <- netAnalyze(net1, clustMethod = "cluster_fast_greedy")

# Network plot (see ?plot.microNetProps for details)
plot(props1)

#----------------------------------------------------------------------------
# Same network as before but on genus level and without taxa filtering

amgut.genus.phy <- phyloseq::tax_glom(amgut2.filt.phy, taxrank = "Rank6")

dim(phyloseq::otu_table(amgut.genus.phy))

# Rename taxonomy table and make Rank6 (genus) unique
amgut.genus.renamed <- renameTaxa(amgut.genus.phy, pat = "<name>",
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Rank6")

net_genus <- netConstruct(amgut.genus.renamed,
                          taxRank = "Rank6",
                          measure = "spieceasi",
                          measurePar = list(method = "mb",
                                            pulsar.params = list(rep.num = 10),
                                            symBetaMode = "ave"),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          sparsMethod = "none",
                          normMethod = "none",
                          verbose = 3)

# Network analysis
props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")

# Network plot (with some modifications)
plot(props_genus,
     shortenLabels = "none",
     labelScale = FALSE,
     cexLabels = 0.8)

#----------------------------------------------------------------------------
# Single network with the following specifications:
# - Association measure: Pearson correlation
# - Taxa filtering: Keep the 50 taxa with highest frequency
# - Sample filtering: Keep samples with a total number of reads of at least
#  1000 and with at least 10 taxa with a non-zero count
# - Zero replacement: A pseudo count of 0.5 is added to all counts
# - Normalization: clr transformation
# - Sparsification: Threshold = 0.3
#  (an edge exists between taxa with an estimated association >= 0.3)

net2 <- netConstruct(amgut2.filt.phy,
                     measure = "pearson",
                     filtTax = "highestFreq",
                     filtTaxPar = list(highestFreq = 50),
                     filtSamp = c("numbTaxa", "totalReads"),
                     filtSampPar = list(totalReads = 1000, numbTaxa = 10),
                     zeroMethod = "pseudo",
                     zeroPar = list(pseudocount = 0.5),
                     normMethod = "clr",
                     sparsMethod = "threshold",
                     thresh = 0.3,
                     verbose = 3)

# Network analysis
props2 <- netAnalyze(net2, clustMethod = "cluster_fast_greedy")

plot(props2)

#----------------------------------------------------------------------------
# Example of using the argument "assoBoot"

# This functionality is useful for splitting up a large number of bootstrap
# replicates and run the bootstrapping procedure iteratively.

niter <- 5
nboot <- 1000
# Overall number of bootstrap replicates: 5000

# Use a different seed for each iteration
seeds <- sample.int(1e8, size = niter)

# List where all bootstrap association matrices are stored
assoList <- list()

for (i in 1:niter) {
  # assoBoot is set to TRUE to return the bootstrap association matrices
  net <- netConstruct(amgut1.filt,
                      filtTax = "highestFreq",
                      filtTaxPar = list(highestFreq = 50),
                      filtSamp = "totalReads",
                      filtSampPar = list(totalReads = 0),
                      measure = "pearson",
                      normMethod = "clr",
                      zeroMethod = "pseudoZO",
                      sparsMethod = "bootstrap",
                      cores = 1,
                      nboot = nboot,
                      assoBoot = TRUE,
                      verbose = 1, # Set to 2 for progress bar
                      seed = seeds[i])
  
  assoList[(1:nboot) + (i - 1) * nboot] <- net$assoBoot1
}

# Construct the actual network with all 5000 bootstrap association matrices
net_final <- netConstruct(amgut1.filt,
                          filtTax = "highestFreq",
                          filtTaxPar = list(highestFreq = 50),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 0),
                          measure = "pearson",
                          normMethod = "clr",
                          zeroMethod = "pseudoZO",
                          sparsMethod = "bootstrap",
                          cores = 1,
                          nboot = nboot * niter,
                          assoBoot = assoList,
                          verbose = 1) # Set to 2 for progress bar

# Network analysis
props <- netAnalyze(net_final, clustMethod = "cluster_fast_greedy")

# Network plot
plot(props)

#----------------------------------------------------------------------------
knitr::opts_chunk$set(fig.width = 16, fig.height = 8)
#----------------------------------------------------------------------------
# Constructing and analyzing two networks
# - A random group variable is used for splitting the data into two groups

set.seed(123456)
group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)

# Option 1: Use the count matrix and group vector as input:
net3 <- netConstruct(amgut1.filt,
                     group = group,
                     measure = "pearson",
                     filtTax = "highestVar",
                     filtTaxPar = list(highestVar = 50),
                     filtSamp = "totalReads",
                     filtSampPar = list(totalReads = 1000),
                     zeroMethod = "multRepl",
                     normMethod = "clr",
                     sparsMethod = "t-test")

# Option 2: Pass the count matrix of group 1 to 'data'
#           and that of group 2 to 'data2'
# Note: Argument 'jointPrepro' is set to FALSE by default (the data sets
# are filtered separately and the intersect of filtered taxa is kept,
# which leads to less than 50 taxa in this example).

amgut1 <- amgut1.filt[group == 1, ]
amgut2 <- amgut1.filt[group == 2, ]

net3 <- netConstruct(data = amgut1,
                     data2 = amgut2,
                     measure = "pearson",
                     filtTax = "highestVar",
                     filtTaxPar = list(highestVar = 50),
                     filtSamp = "totalReads",
                     filtSampPar = list(totalReads = 1000),
                     zeroMethod = "multRepl",
                     normMethod = "clr",
                     sparsMethod = "t-test")

# Network analysis
# Note: Please zoom into the GCM plot or open a new window using:
# x11(width = 10, height = 10)
props3 <- netAnalyze(net3, clustMethod = "cluster_fast_greedy")

# Network plot (same layout is used in both groups)
plot(props3, sameLayout = TRUE)

# The two networks can be compared with NetCoMi's function netCompare().