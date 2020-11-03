checkPack <- function(measure, zeroMethod, normMethod, sparsMethod, adjust){
  
  needpack <- character()

  if(measure == "ccrepe"){
    needpack <- c(needpack, "ccrepe")
  } else if(measure == "propr"){
    needpack <- c(needpack, "propr")
  } else if(measure %in% c("kld", "jeffrey", "jsd")){
    needpack <- c(needpack, "LaplacesDemon")
  }

  if(zeroMethod %in% c("multRepl", "alrEM", "bayesMult")){
    needpack <- c(needpack, "zCompositions")
  }
  
  if(normMethod == "CSS"){
    needpack <- c(needpack, "metagenomeSeq")
  } else if(normMethod == "VST"){
    needpack <- c(needpack, "DESeq2")
  } 

  if(sparsMethod %in% c("t-test", "bootstrap") && adjust == "adaptBH"){
    needpack <- c(needpack, "limma")
    
  } else if(sparsMethod == "knn"){
    needpack <- c(needpack, "cccd")
  }
  
  instpack <- needpack[!needpack %in% utils::installed.packages()[,"Package"]]

  if(length(instpack) != 0){
    
    if(length(instpack) > 1){
      message("Installing missing packages: ", 
              paste0(instpack[1:(length(instpack)-1)], ", ", 
                     instpack[length(instpack)]))
    } else{
      message("Installing missing package: ", instpack)
    }
    
    if(!requireNamespace("BiocManager", quietly = TRUE)){
      utils::install.packages("BiocManager")
    }
    
    BiocManager::install(pkgs = instpack, dependencies = TRUE)
    
    message("Done.")
    
    message("Check whether installed package(s) can be loaded ...")
    
    for(pack in instpack){
      requireNamespace(pack)
    }
    
    message("Done.")
  }

}
