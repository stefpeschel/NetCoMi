knitr::opts_chunk$set(fig.width = 10, fig.height = 6)

# Load data sets from American Gut Project (from SpiecEasi package)
data("amgut1.filt")

# Generate a random group vector
set.seed(123456)
group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)

# Network construction:
amgut_net <- netConstruct(amgut1.filt, group = group,
                          measure = "pearson",
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 30),
                          zeroMethod = "pseudoZO", normMethod = "clr")

#---------------------
# Differential network

# Fisher's z-test
amgut_diff1 <- diffnet(amgut_net, diffMethod = "fisherTest")

# Network contains no differentially correlated taxa:
\dontrun{
  plot(amgut_diff1)
}

# Without multiple testing correction (statistically not correct!)
amgut_diff2 <- diffnet(amgut_net, diffMethod = "fisherTest", adjust = "none")
plot(amgut_diff2)

\dontrun{
  # Permutation test (permutation matrices are stored)
  amgut_diff3 <- diffnet(amgut_net, 
                         diffMethod = "permute", 
                         nPerm = 1000L,
                         cores = 4L, 
                         adjust = "lfdr",
                         storeCountsPerm = TRUE,
                         fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
                         storeAssoPerm = TRUE,
                         fileStoreAssoPerm = "assoPerm",
                         seed = 123456)
  
  # Use the p-values again (different adjustment method possible), but without
  # re-estimating the associations
  amgut_diff4 <- diffnet(amgut_net, 
                         diffMethod = "permute", 
                         nPerm = 1000L,
                         adjust = "none", 
                         pvalsVec = amgut_diff3$pvalsVec)
  x11()
  plot(amgut_diff4)
  
  # Use the permutation associations again (same result as amgut_diff4)
  amgut_diff5 <- diffnet(amgut_net, 
                         diffMethod = "permute", 
                         nPerm = 1000L,
                         adjust = "none",
                         fileLoadAssoPerm = "assoPerm")
  x11()
  plot(amgut_diff5)
  
  # Use the permuted count matrices again (same result as amgut_diff4)
  amgut_diff6 <- diffnet(amgut_net, 
                         diffMethod = "permute", 
                         nPerm = 1000L,
                         adjust = "none",
                         fileLoadCountsPerm = c("countsPerm1", "countsPerm2"),
                         seed = 123456)
  x11()
  plot(amgut_diff6)
}