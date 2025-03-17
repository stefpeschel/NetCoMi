knitr::opts_chunk$set(fig.width = 16, fig.height = 8)

# Load data sets from American Gut Project (from SpiecEasi package)
data("amgut2.filt.phy")

# Split data into two groups: with and without seasonal allergies
amgut_season_yes <- phyloseq::subset_samples(amgut2.filt.phy,
                                      SEASONAL_ALLERGIES == "yes")
amgut_season_no <- phyloseq::subset_samples(amgut2.filt.phy,
                                     SEASONAL_ALLERGIES == "no")

amgut_season_yes
amgut_season_no

# Filter the 121 samples (sample size of the smaller group) with highest
# frequency to make the sample sizes equal and thus ensure comparability.
n_yes <- phyloseq::nsamples(amgut_season_yes)

# Network construction
amgut_net <- netConstruct(data = amgut_season_yes,
                          data2 = amgut_season_no,
                          measure = "pearson",
                          filtSamp = "highestFreq",
                          filtSampPar = list(highestFreq = n_yes),
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 30),
                          zeroMethod = "pseudoZO", normMethod = "clr")

# Network analysis
# Note: Please zoom into the GCM plot or open a new window using:
# x11(width = 10, height = 10)
amgut_props <- netAnalyze(amgut_net, clustMethod = "cluster_fast_greedy")

# Network plot
plot(amgut_props,
     sameLayout = TRUE,
     title1 = "Seasonal allergies",
     title2 = "No seasonal allergies")

#--------------------------
# Network comparison

# Without permutation tests
amgut_comp1 <- netCompare(amgut_props, permTest = FALSE)
summary(amgut_comp1)

\donttest{
  # With permutation tests (with only 100 permutations to decrease runtime)
  amgut_comp2 <- netCompare(amgut_props,
                            permTest = TRUE,
                            nPerm = 100L,
                            cores = 1L,
                            storeCountsPerm = TRUE,
                            fileStoreCountsPerm = c("countsPerm1",
                                                    "countsPerm2"),
                            storeAssoPerm = TRUE,
                            fileStoreAssoPerm = "assoPerm",
                            seed = 123456)

# Rerun with a different adjustment method ...
# ... using the stored permutation count matrices
amgut_comp3 <- netCompare(amgut_props, adjust = "BH",
                          permTest = TRUE, nPerm = 100L,
                          fileLoadCountsPerm = c("countsPerm1",
                                                 "countsPerm2"),
                          seed = 123456)

# ... using the stored permutation association matrices
amgut_comp4 <- netCompare(amgut_props, adjust = "BH",
                          permTest = TRUE, nPerm = 100L,
                          fileLoadAssoPerm = "assoPerm",
                          seed = 123456)

# amgut_comp3 and amgut_comp4 should be equal
all.equal(amgut_comp3$adjaMatrices, amgut_comp4$adjaMatrices)
all.equal(amgut_comp3$properties, amgut_comp4$properties)

summary(amgut_comp2)
summary(amgut_comp3)
summary(amgut_comp4)

#--------------------------
# Use 'createAssoPerm' to create "permuted" count and association matrices
createAssoPerm(amgut_props, nPerm = 100,
               computeAsso = TRUE,
               fileStoreAssoPerm = "assoPerm",
               storeCountsPerm = TRUE,
               fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
               append = FALSE, seed = 123456)

amgut_comp5 <- netCompare(amgut_props, permTest = TRUE, nPerm = 100L,
                          fileLoadAssoPerm = "assoPerm")

all.equal(amgut_comp3$properties, amgut_comp5$properties)

summary(amgut_comp5)
}