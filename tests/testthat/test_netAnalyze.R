set.seed(123456)
data("amgut1.filt")
data("amgut2.filt.phy")

groups_diss <- sample(0:1, ncol(amgut1.filt), replace = TRUE)
groups_asso <- sample(0:1, nrow(amgut1.filt), replace = TRUE)

net_asso_single <- netConstruct(amgut1.filt,
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 20),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 1000),
                          zeroMethod = "none", normMethod = "none",
                          measure = "pearson",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)

net_unweight_single <- netConstruct(amgut1.filt,
                                filtTax = "highestVar",
                                filtTaxPar = list(highestVar = 20),
                                filtSamp = "totalReads",
                                filtSampPar = list(totalReads = 1000),
                                zeroMethod = "none", normMethod = "none",
                                measure = "pearson",
                                sparsMethod = "threshold", thresh = 0.3,
                                seed = 20190101,
                                weighted = FALSE)

net_diss_single <- netConstruct(amgut1.filt,
                          filtTax = "totalReads",
                          filtTaxPar = list(totalReads = 1000),
                          filtSamp = "highestFreq",
                          filtSampPar = list(highestFreq = 20),
                          zeroMethod = "none", normMethod = "none",
                          measure = "aitchison",
                          sparsMethod = "threshold", thresh = 0.3,
                          seed = 20190101)


net_asso_two <- netConstruct(amgut1.filt, group = groups_asso,
                                filtTax = "highestVar",
                                filtTaxPar = list(highestVar = 20),
                                filtSamp = "totalReads",
                                filtSampPar = list(totalReads = 1000),
                                zeroMethod = "none", normMethod = "none",
                                measure = "pearson",
                                sparsMethod = "threshold", thresh = 0.3,
                                seed = 20190101)

net_unweight_two <- netConstruct(amgut1.filt, group = groups_asso,
                                    filtTax = "highestVar",
                                    filtTaxPar = list(highestVar = 20),
                                    filtSamp = "totalReads",
                                    filtSampPar = list(totalReads = 1000),
                                    zeroMethod = "none", normMethod = "none",
                                    measure = "pearson",
                                    sparsMethod = "threshold", thresh = 0.3,
                                    seed = 20190101,
                                    weighted = FALSE)

net_diss_two <- netConstruct(amgut1.filt, group = groups_diss,
                                filtTax = "totalReads",
                                filtTaxPar = list(totalReads = 1000),
                                filtSamp = "highestFreq",
                                filtSampPar = list(highestFreq = 150),
                                zeroMethod = "none", normMethod = "none",
                                measure = "aitchison",
                                sparsMethod = "threshold", thresh = 0.3,
                                seed = 20190101)

networks <- c("net_asso_single", "net_diss_single", "net_unweight_single",
              "net_asso_two", "net_diss_two", "net_unweight_two")

#===============================================================================
# centrLCC

context("netAnalyze: Test centrLCC")

centrLCC <- c(TRUE, FALSE)

for(i in 1:length(centrLCC)){
  context(centrLCC[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           centrLCC = centrLCC[i])
    summary(testprops)
  }
}

#===============================================================================
# avDissIgnoreInf

context("netAnalyze: Test avDissIgnoreInf ")

avDissIgnoreInf <- c(TRUE, FALSE)

for(i in 1:length(avDissIgnoreInf)){
  context(avDissIgnoreInf[i])

  for(net in networks){
    #cat("avDissIgnoreInf: ", avDissIgnoreInf[i], "\n")
    #print(net)
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           avDissIgnoreInf = avDissIgnoreInf[i])
    
    summary(testprops)
    #print(testprops$globalProps$avDiss1)
    #print(testprops$globalPropsLCC$avDiss1)
  }
}

#===============================================================================
# sPathAlgo

context("netAnalyze: Test sPathAlgo")

sPathAlgo <- c("unweighted", "dijkstra", "bellman-ford", "johnson", "automatic")

for(i in 1:length(sPathAlgo)){
  context(sPathAlgo[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           sPathAlgo = sPathAlgo[i])
    
    summary(testprops)
  }
}

#===============================================================================
# sPathNorm

context("netAnalyze: Test sPathNorm")

sPathNorm <- c(TRUE, FALSE)

for(i in 1:length(sPathNorm)){
  context(sPathNorm[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           sPathNorm = sPathNorm[i])
    
    summary(testprops)
  }
}

#===============================================================================
# normNatConnect

context("netAnalyze: Test normNatConnect")

normNatConnect <- c(TRUE, FALSE)

for(i in 1:length(normNatConnect)){
  context(normNatConnect[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           normNatConnect = normNatConnect[i])
    
    summary(testprops)
  }
}

#===============================================================================
# connectivity

context("netAnalyze: Test connectivity")

connectivity <- c(TRUE, FALSE)

for(i in 1:length(connectivity)){
  context(connectivity[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           connectivity = connectivity[i])
    
    summary(testprops)
  }
}


#===============================================================================
# Clustering

context("netAnalyze: Test clustering methods")

clustMethod <- c("none", 
                 "hierarchical",
                 #"cluster_optimal",
                 "cluster_fast_greedy",
                 "cluster_louvain",
                 #"cluster_edge_betweenness",
                 "cluster_leading_eigen",
                 #"cluster_spinglass",
                 "cluster_walktrap")

for(i in 1:length(clustMethod)){
  context(clustMethod[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = clustMethod[i],
                           hubPar = "eigenvector", hubQuant = 0.95)
    
    summary(testprops)
  }
}

context("netAnalyze: Test hierarchical clustering")

for(i in 1:length(clustMethod)){
  context("hierarchical")
  
  for(net in networks[4:6]){
    testprops<- netAnalyze(get(net),
                           clustMethod = "hierarchical",
                           clustPar = list(method = "complete", k = 3),
                           clustPar2 = list(method = "complete", k = 2),
                           hubPar = "eigenvector", hubQuant = 0.95)
    
    summary(testprops)
  }
}

#===============================================================================
# weightClustCoef

context("netAnalyze: Test weightClustCoef")

weightClustCoef <- c(TRUE, FALSE)

for(i in 1:length(weightClustCoef)){
  context(weightClustCoef[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           weightClustCoef  = weightClustCoef[i])
    
    summary(testprops)
  }
}

#===============================================================================
# hubPar

context("netAnalyze: Test hubPar with multiple hub parameters")

hubPar <- c("degree", "betweenness", "closeness")

for(net in networks){
  testprops<- netAnalyze(get(net),
                         clustMethod = "cluster_fast_greedy",
                         hubPar = hubPar, hubQuant = 0.95)
  
  summary(testprops)
}

context("netAnalyze: Test hubPar with eigencentr")

hubPar <- c("eigenvector")

for(net in networks){
  testprops<- netAnalyze(get(net),
                         clustMethod = "cluster_fast_greedy",
                         hubPar = hubPar, hubQuant = 0.95)
  
  summary(testprops)
}

#===============================================================================
# lnormFit

context("netAnalyze: Test lnormFit ")

lnormFit <- c(TRUE, FALSE)

for(i in 1:length(lnormFit)){
  context(lnormFit[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           lnormFit = lnormFit[i])
    
    summary(testprops)
  }
}

#===============================================================================
# weightDeg

context("netAnalyze: Test weightDeg ")

weightDeg <- c(TRUE, FALSE)

for(i in 1:length(weightDeg)){
  context(weightDeg[i])
  
  for(net in networks){
    testprops<- netAnalyze(get(net),
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector", hubQuant = 0.95,
                           weightDeg = weightDeg[i])
    
    summary(testprops, showCentr = "degree")
  }
}


#===============================================================================
context("netAnalyze: with non-normalized centralities")

for(net in networks){
testprops<- netAnalyze(get(net), clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector", centrLCC = TRUE,
                       normDeg = FALSE,
                       normBetw = FALSE, normClose = FALSE, normEigen = FALSE,
                       hubQuant = 0.95)

summary(testprops)
}


#===============================================================================
context("plot.microNetProps: test different layouts")

testprops<- netAnalyze(net_asso_two, clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector")

plot(testprops, sameLayout = FALSE)
plot(testprops, sameLayout = TRUE)
plot(testprops, sameLayout = TRUE, layoutGroup = 1)
plot(testprops, sameLayout = TRUE, layout = "layout_with_fr")

graph.tmp <- graph_from_adjacency_matrix(testprops$input$adjaMat1, weighted = TRUE)
lay.tmp <- layout_with_fr(graph.tmp)
rownames(lay.tmp) <- rownames(testprops$input$adjaMat1)
plot(testprops, layout = lay.tmp)

# unweighted network
testprops<- netAnalyze(net_unweight_two, clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector")

plot(testprops, sameLayout = FALSE)
plot(testprops, sameLayout = TRUE)
plot(testprops, sameLayout = TRUE, layoutGroup = 1)
plot(testprops, sameLayout = TRUE, layout = "layout_with_fr")

graph.tmp <- graph_from_adjacency_matrix(testprops$input$adjaMat1, weighted = TRUE)
lay.tmp <- layout_with_fr(graph.tmp)
rownames(lay.tmp) <- rownames(testprops$input$adjaMat1)
plot(testprops, layout = lay.tmp)


# dissimilarity network
testprops<- netAnalyze(net_diss_two, clustMethod = "cluster_fast_greedy",
                       hubPar = "eigenvector")

plot(testprops, sameLayout = FALSE)
plot(testprops, sameLayout = TRUE, repulsion = 0.7)
plot(testprops, sameLayout = TRUE, layoutGroup = 1)
plot(testprops, sameLayout = TRUE, layout = "layout_with_fr")

graph.tmp <- graph_from_adjacency_matrix(testprops$input$adjaMat1, weighted = TRUE)
lay.tmp <- layout_with_fr(graph.tmp)
rownames(lay.tmp) <- rownames(testprops$input$adjaMat1)
plot(testprops, layout = lay.tmp)










