knitr::opts_chunk$set(fig.width = 8, fig.height = 8)

# Load data sets from American Gut Project (from SpiecEasi package)
data("amgut1.filt")

# Network construction
amgut_net1 <- netConstruct(amgut1.filt, measure = "pearson",
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 50),
                           zeroMethod = "pseudoZO", normMethod = "clr",
                           sparsMethod = "threshold", thresh = 0.4)

# Network analysis

# Using eigenvector centrality as hub score
amgut_props1 <- netAnalyze(amgut_net1, clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector")

summary(amgut_props1, showCentr = "eigenvector", numbNodes = 15L, digits = 3L)

# Using degree, betweenness and closeness centrality as hub scores
amgut_props2 <- netAnalyze(amgut_net1, clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "betweenness", "closeness"))

summary(amgut_props2, showCentr = "all",  numbNodes = 5L, digits = 5L)

# Calculate centralities only for the largest connected component
amgut_props3 <- netAnalyze(amgut_net1, centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector")

summary(amgut_props3, showCentr = "none", clusterLCC = TRUE)

# Network plot
plot(amgut_props1)
plot(amgut_props2)
plot(amgut_props3)

#----------------------------------------------------------------------------
# Plot the GCM heatmap
plotGCM(amgut_props1)

#----------------------------------------------------------------------------
# Dissimilarity-based network (where nodes are subjects)
amgut_net4 <- netConstruct(amgut1.filt, measure = "aitchison",
                           filtSamp = "highestFreq",
                           filtSampPar = list(highestFreq = 30),
                           zeroMethod = "multRepl", sparsMethod = "knn")

amgut_props4 <- netAnalyze(amgut_net4, clustMethod = "hierarchical",
                           clustPar = list(k = 3))

plot(amgut_props4)
