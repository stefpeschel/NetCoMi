#' @keywords internal
.getPermWorkerPackages <- function(measure,
                                   zeroMethod = "none",
                                   normMethod = "none",
                                   includePrepro = TRUE,
                                   useFilematrix = TRUE) {
  packages <- character()
  
  if (isTRUE(useFilematrix)) {
    packages <- c(packages, "filematrix")
  }
  
  packages <- c(
    packages,
    switch(measure,
           bicor = "WGCNA",
           corshrink = "corpcor",
           sparcc = "SpiecEasi",
           ccrepe = "ccrepe",
           reboot = "ccrepe",
           spieceasi = "SpiecEasi",
           spring = c("SPRING", "SpiecEasi"),
           propr = "propr",
           rhoshrink = c("corpcor", "SpiecEasi"),
           bray = "vegan",
           euclidean = "vegan",
           kld = "LaplacesDemon",
           jeffrey = "LaplacesDemon",
           jsd = "LaplacesDemon",
           aitchison = c("SpiecEasi", "vegan"),
           character())
  )
  
  if (isTRUE(includePrepro)) {
    if (zeroMethod %in% c("multRepl", "alrEM", "bayesMult")) {
      packages <- c(packages, "zCompositions")
    }
    
    packages <- c(
      packages,
      switch(normMethod,
             CSS = "metagenomeSeq",
             rarefy = "vegan",
             VST = "DESeq2",
             clr = "SpiecEasi",
             mclr = "SPRING",
             character())
    )
  }
  
  unique(packages)
}
