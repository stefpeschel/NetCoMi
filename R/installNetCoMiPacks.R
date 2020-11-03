#' @title Install all packages used within NetCoMi
#'
#' @description This function installs the R packages used in NetCoMi not listed
#'   as dependencies or imports in NetCoMi's description file. These are 
#'   optional packages only needed in certain network construction settings.\cr
#'   
#'   \code{BiocManager::\link[BiocManager]{install}} is used for installation 
#'   since it installs or updates Bioconductor as well as CRAN packages.\cr
#'   
#'   Installed CRAN packages:
#'   \itemize{
#'   \item cccd
#'   \item LaplacesDemon
#'   \item propr
#'   \item zCompositions
#'   }
#'   
#'   Installed Bioconductor packages:
#'   \itemize{
#'   \item ccrepe
#'   \item DESeq2
#'   \item discordant
#'   \item limma
#'   \item metagenomeSeq
#'   }
#'   
#'   If not installed via this function, the packages are installed by the 
#'   respective NetCoMi functions when needed.
#'
#' @param onlyMissing logical. If \code{TRUE} (default), 
#'   \code{\link[utils]{installed.packages}} is used to read out the packages
#'   installed in the given library and only missing packages are installed. If
#'   \code{FALSE}, all packages are installed or updated (if already installed).
#' @param lib character vector giving the library directories where to install 
#'   missing packages. If \code{NULL}, the first element of 
#'   \code{\link{.libPaths}} is used.
#' @param ... Additional arguments used by \code{\link[BiocManager]{install}} or
#'   \code{\link[utils]{install.packages}}.
#' 
#' @export

installNetCoMiPacks <- function(onlyMissing = TRUE, lib = NULL, ...){
  
  if(is.null(lib)) lib <- .libPaths()[1]
  
  needpack <- c("cccd", "ccrepe", "DESeq2", "discordant", "LaplacesDemon", 
                "limma", "metagenomeSeq", "propr", "zCompositions")
  
  if(onlyMissing){
    instpack <- needpack[!needpack %in% utils::installed.packages()[,"Package"]]
  } else{
    instpack <- needpack
  }

  if(length(instpack) != 0){
    
    if(length(instpack) > 1){
      message("Installing packages using BiocManager: ", 
              paste0(instpack[1:(length(instpack)-1)], ", "), 
              instpack[length(instpack)])
    } else{
      message("Installing package using BiocManager: ", instpack)
    }
    
    if(!requireNamespace("BiocManager", quietly = TRUE)){
      utils::install.packages("BiocManager")
    }
    
    BiocManager::install(pkgs = instpack, lib = lib, ...)  
    
    message("Done.")
    message("Check whether installed package(s) can be loaded ...")
    
    for(pack in instpack){
      requireNamespace(pack)
    }
    
    message("Done.")
    
  } else{
    message("All packages already installed.")
  }

}