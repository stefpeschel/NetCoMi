# Load data sets from American Gut Project (from SpiecEasi package)
data("amgut2.filt.phy")

# Split data into two groups: with and without seasonal allergies
amgut_season_yes <- phyloseq::subset_samples(amgut2.filt.phy, 
                                      SEASONAL_ALLERGIES == "yes")
amgut_season_no <- phyloseq::subset_samples(amgut2.filt.phy, 
                                     SEASONAL_ALLERGIES == "no")

# Network construction
net <- netConstruct(data = amgut_season_yes,
                    data2 = amgut_season_no,
                    measure = "pearson",
                    filtTax = "highestVar",
                    filtTaxPar = list(highestVar = 50),
                    zeroMethod = "pseudoZO", 
                    normMethod = "clr",
                    sparsMethod = "thresh",
                    dissFunc = "signed",
                    thresh = 0.3,
                    seed = 123456)

# Get adjacency matrices
adja1 <- net$adjaMat1
adja2 <- net$adjaMat2

# Network visualization
props <- netAnalyze(net)

plot(props, rmSingles = TRUE, cexLabels = 1.7, edgeWidth = 1.2)

gcm1 <- calcGCM(adja1)
gcm2 <- calcGCM(adja2)
gcd <- calcGCD(adja1, adja2)


gcddiff <- testGCM(obj1 = gcm1, obj2 = gcm2, adjust = "adaptBH")

plotHeat(gcm1$ocount, textUpp = "none", labPos = "lt")

plotHeat(mat = gcm1$gcm, diag = FALSE)

plotHeat(mat = gcm2$gcm)

plotHeat(mat = gcddiff$diff)

plotHeat(mat = gcddiff$absDiff)


# error
expect_error(
  plotHeat(mat = gcddiff$diff,
           type = "mixed",
           textUpp = "mat",
           textLow = "pmat"))

#===============================================================================
# Test the plotHeat function
testfunc <- function(type = c("mixed", "full", "lower", "upper"), ...) {

  upper <- lower <- c("mat", "sigmat", "pmat", "code", "none")

  if ("mixed" %in% type) {
    t <- "mixed"
    
    for (u in upper) {
      for (l in lower) {
        plotHeat(mat = gcddiff$diff,
                pmat = gcddiff$pvalsDiff,
                type = t,
                textUpp = u,
                textLow = l,
                title = paste0(t, ", u=", u, ", l=", l), 
                mar = c(0,0,1,0),
                ...)
      }
    }
  } 
  
  if ("full" %in% type) {
    t <- "full"
    l <- ""
    for (u in upper) {
      plotHeat(mat = gcddiff$diff,
              pmat = gcddiff$pvalsDiff,
              type = t,
              textUpp = u,
              title = paste0(t, ", u=", u, ", l=", l), 
              mar = c(0,0,1,0),
              ...)
    }
  }
  
  if ("upper" %in% type) {
    t <- "upper"
    l <- "--"
    for (u in upper) {
      plotHeat(mat = gcddiff$diff,
              pmat = gcddiff$pvalsDiff,
              type = t,
              textUpp = u,
              title = paste0(t, ", u=", u, ", l=", l), 
              mar = c(0,0,1,0),
              ...)
    }
  }
  
  if ("lower" %in% type) {
    t <- "lower"
    u <- "--"
    for (l in lower) {
      plotHeat(mat = gcddiff$diff,
              pmat = gcddiff$pvalsDiff,
              type = t,
              textLow = l,
              title = paste0(t, ", u=", u, ", l=", l), 
              mar = c(0,0,1,0),
              ...)
    }
  }
}

testfunc()
testfunc(diag = FALSE)
testfunc(labPos = "d")
testfunc("mixed", methUpp = "circle")
testfunc("mixed", methLow = "circle")
testfunc(type = "lower", labPos = "ld", diag = FALSE)
testfunc(type = "upper", labPos = "td", diag = FALSE)
testfunc("mixed", labPos = "n")

plotHeat(mat = gcm1$gcm, colorPal = "Blues")
plotHeat(mat = gcm1$gcm, colorPal = "PiYG")

plotHeat(mat = gcm1$gcm, colorPal = "PiYG", colorLim = c(-0.5, 1))

plotHeat(mat = gcm1$gcm, 
        pmat = gcddiff$pvalsDiff,
        type = "mixed",
        colorPal = "PiYG", 
        revCol = TRUE)

plotHeat(mat = gcm1$gcm, 
        pmat = gcddiff$pvalsDiff,
        type = "mixed",
        nCol = 10)

plotHeat(mat = gcm1$gcm, 
        pmat = gcddiff$pvalsDiff,
        type = "mixed",
        colorPal = "PiYG",
        nCol = 11,
        addWhite = FALSE)

plotHeat(mat = gcm1$gcm, 
        pmat = gcddiff$pvalsDiff,
        type = "mixed",
        textUpp = "code",
        textLow = "mat",
        diag = FALSE,
        color = rev(hcl.colors(50, palette = "viridis")))


plotHeat(mat = gcm1$gcm, 
        pmat = gcddiff$pvalsDiff,
        type = "mixed",
        argsUpp = list(method = "circle"),
        argsLow = list(method = "pie"))
