knitr::opts_chunk$set(fig.width = 8, fig.height = 8)

# Load data sets from American Gut Project (from SpiecEasi package)
data("amgut1.filt")

# Network construction
amgut_net <- netConstruct(amgut1.filt, measure = "pearson",
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 50),
                          zeroMethod = "pseudoZO", normMethod = "clr",
                          sparsMethod = "threshold", thresh = 0.3)

# Network analysis
amgut_props <- netAnalyze(amgut_net)

### Network plots ###
# Clusters are used for node coloring:
plot(amgut_props,
     nodeColor = "cluster")

# Remove singletons
plot(amgut_props,
     nodeColor = "cluster",
     rmSingles = TRUE)

# A higher repulsion places nodes with high edge weight closer together
plot(amgut_props,
     nodeColor = "cluster",
     rmSingles = TRUE,
     repulsion = 1.2)

# A feature vector is used for node coloring
# (this could be a vector with phylum names of the ASVs)
set.seed(123456)
featVec <- sample(1:5, nrow(amgut1.filt), replace = TRUE)

# Names must be equal to ASV names
names(featVec) <- colnames(amgut1.filt)

plot(amgut_props,
     rmSingles = TRUE,
     nodeColor = "feature",
     featVecCol = featVec,
     colorVec = heat.colors(5))

# Use a further feature vector for node shapes
shapeVec <- sample(1:3, ncol(amgut1.filt), replace = TRUE)
names(shapeVec) <- colnames(amgut1.filt)

plot(amgut_props,
     rmSingles = TRUE,
     nodeColor = "feature",
     featVecCol = featVec,
     colorVec = heat.colors(5),
     nodeShape = c("circle", "square", "diamond"),
     featVecShape = shapeVec,
     highlightHubs = FALSE)