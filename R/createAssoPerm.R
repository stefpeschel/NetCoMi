#' @title Create and store association matrices for permuted data
#'
#' @description The function creates and returns a matrix with permuted group 
#'   labels and saves association matrices computed for the permuted data to 
#'   an external file.
#'
#' @param x object of class \code{"microNet"} or \code{"microNetProps"} 
#'   (returned by \code{\link[NetCoMi]{netConstruct}} or 
#'   \code{\link[NetCoMi]{netAnalyze}}).
#' @param computeAsso logical indicating whether the association matrices should
#'   be computed. If \code{FALSE}, only the permuted group labels are computed 
#'   and returned.
#' @param nPerm integer indicating the number of permutations.
#' @param cores integer indicating the number of CPU cores used for
#'   permutation tests. If cores > 1, the tests are performed in parallel.
#'   Is limited to the number of available CPU cores determined by
#'   \code{\link[parallel]{detectCores}}. Defaults to 1L (no parallelization).
#' @param seed integer giving a seed for reproducibility of the results.
#' @param permGroupMat an optional matrix with permuted group labels 
#' (with nPerm rows and n1+n2 columns). 
#' @param fileStoreAssoPerm character giving the name of a file to which the 
#'   matrix with associations/dissimilarities of the permuted data is saved. 
#'   Can also be a path.
#' @param append logical indicating whether existing files (given by 
#'   fileStoreAssoPerm and fileStoreCountsPerm) should be extended. 
#'   If \code{TRUE}, a new file is created only if the file is not existing. 
#'   If \code{FALSE}, a new file is created in any case.
#' @param storeCountsPerm logical indicating whether the permuted count matrices
#'   should be saved to an external file. Defaults to \code{FALSE}.
#'   Ignored if \code{fileLoadCountsPerm} is not \code{NULL}.
#' @param fileStoreCountsPerm character vector with two elements giving the 
#'   names of two files storing the permuted count matrices belonging to the 
#'   two groups.
#' @param logFile character string naming the log file to which the current 
#'   iteration number is written. Defaults
#'   to \code{NULL} so that no log file is generated.
#' @param verbose logical. If \code{TRUE} (default), status messages are shown.
#' 
#' @return Invisible object: Matrix with permuted group labels.
#' 
#' @examples 
#' \donttest{
#'   # Load data sets from American Gut Project (from SpiecEasi package)
#'   data("amgut1.filt")
#' 
#'   # Generate a random group vector
#'   set.seed(123456)
#'   group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)
#' 
#'   # Network construction:
#'   amgut_net <- netConstruct(amgut1.filt, group = group,
#'                             measure = "pearson",
#'                             filtTax = "highestVar",
#'                             filtTaxPar = list(highestVar = 30),
#'                             zeroMethod = "pseudoZO", normMethod = "clr")
#' 
#'   # Network analysis:
#'   amgut_props <- netAnalyze(amgut_net, clustMethod = "cluster_fast_greedy")
#' 
#'   # Use 'createAssoPerm' to create "permuted" count and association matrices,
#'   # which can be reused by netCompare() and diffNet()
#'   # Note: 
#'   # createAssoPerm() accepts objects 'amgut_net' and 'amgut_props' as input
#'   
#'   createAssoPerm(amgut_props, nPerm = 100L, 
#'                  computeAsso = TRUE,
#'                  fileStoreAssoPerm = "assoPerm",
#'                  storeCountsPerm = TRUE, 
#'                  fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
#'                  append = FALSE, seed = 123456)
#'   
#'   # Run netcompare using the stored permutation count matrices 
#'   # (association matrices are still computed within netCompare):
#'   amgut_comp1 <- netCompare(amgut_props, permTest = TRUE, nPerm = 100L, 
#'                             fileLoadCountsPerm = c("countsPerm1", 
#'                                                    "countsPerm2"),
#'                             seed = 123456)
#'                             
#'   # Run netcompare using the stored permutation association matrices:
#'   amgut_comp2 <- netCompare(amgut_props, permTest = TRUE, nPerm = 100L, 
#'                             fileLoadAssoPerm = "assoPerm")
#'   
#'   summary(amgut_comp1)
#'   summary(amgut_comp2)
#'   all.equal(amgut_comp1$properties, amgut_comp2$properties)
#'   
#'   # Run diffnet using the stored permutation count matrices in diffnet()
#'   diff1 <- diffnet(amgut_net, diffMethod = "permute", nPerm = 100L, 
#'                   fileLoadCountsPerm = c("countsPerm1", "countsPerm2"))
#'                   
#'   # Run diffnet using the stored permutation association matrices 
#'   diff2 <- diffnet(amgut_net, diffMethod = "permute", nPerm = 100L, 
#'                   fileLoadAssoPerm = "assoPerm")
#'                  
#'  #plot(diff1)
#'  #plot(diff2)
#'  # Note: Networks are empty (no significantly different associations) 
#'  # for only 100 permutations
#' }
#' @import foreach doSNOW filematrix
#' @importFrom stats sd
#' @importFrom WGCNA randIndex
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export

createAssoPerm <- function(x,
                           computeAsso = TRUE,
                           nPerm = 1000L,
                           cores = 1L,
                           seed = NULL, 
                           permGroupMat = NULL,
                           fileStoreAssoPerm = "assoPerm",
                           append = TRUE,
                           storeCountsPerm = FALSE,
                           fileStoreCountsPerm = c("countsPerm1", 
                                                   "countsPerm2"),
                           logFile = NULL, 
                           verbose = TRUE) {
  
  # Check input arguments
  args_in <- as.list(environment())
  
  args_out <- .checkArgsCreateAP(args_in)
  
  for (i in 1:length(args_out)) {
    assign(names(args_out)[i], args_out[[i]])
  }
  
  #-----------------------------------------------------------------------------
  if (inherits(x, "microNet")) {
    
    # Parameters used for network construction
    parNC <- x$parameters
    
    assoType <- x$assoType
    xgroups <- x$groups
    matchDesign <- x$matchDesign
    twoNets <- x$twoNets
    callNetConstr <- x$call
    
    count1 <- x$countMat1
    count2 <- x$countMat2
    countsJoint <- x$countsJoint
    
    count_norm1 <- x$normCounts1
    count_norm2 <- x$normCounts2
    
  } else { # class=microNetProps
    
    parNC <- x$paramsNetConstruct
    
    assoType <- x$input$assoType
    xgroups <- x$input$groups 
    matchDesign <- x$input$matchDesign
    twoNets <- x$input$twoNets
    callNetConstr <- x$input$call
    
    count1 <- x$input$countMat1
    count2 <- x$input$countMat2
    countsJoint <- x$input$countsJoint
    
    count_norm1 <- x$input$normCounts1
    count_norm2 <- x$input$normCounts2
  }

  distNet <- ifelse(assoType == "dissimilarity", TRUE, FALSE)

  if (!twoNets) {
    stop("Input contains only a single network.")
  }
  
  #---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  # create combined data matrix
  if (assoType == "dissimilarity") {
    n1 <- ncol(count1)
    n2 <- ncol(count2)
    n <- n1 + n2
    nVar <- nrow(count1)
    
    xbind <- cbind(count1, count2)
    
  } else {
    if (parNC$jointPrepro) {
      n1 <- nrow(count_norm1)
      n2 <- nrow(count_norm2)
      nVar <- ncol(count_norm1)
      
      xbind <- rbind(count_norm1, count_norm2)
      
    } else {
      n1 <- nrow(count1)
      n2 <- nrow(count2)
      nVar <- ncol(count1)
      
      xbind <- rbind(count1, count2)
    }
    
    n <- n1 + n2
  }
  
  
  #---------------------------------------
  # create matrix with permuted group labels
  
  if (is.null(permGroupMat)) {
    
    if (verbose) {
      message("Create matrix with permuted group labels ... ", appendLF = FALSE)
    }
    
    perm_group_mat <- .getPermGroupMat(n1 = n1, n2 = n2, n = n, nPerm = nPerm, 
                                         matchDesign = matchDesign)
    
    if (verbose) {
      message("Done.")
    }
    
  } else {
    stopifnot(ncol(permGroupMat) == n)
    
    perm_group_mat <- permGroupMat
  }
  
  #---------------------------------------
  # compute and store association matrices
  
  if (computeAsso) {
    
    stopifnot(is.character(fileStoreAssoPerm))
    stopifnot(length(fileStoreAssoPerm) == 1)
    
    if (append) {
      fmat <- suppressWarnings(try(fm.open(filenamebase = fileStoreAssoPerm), 
                                   silent = TRUE))
      
      if (inherits(fmat, "try-error")) {
        fmat <- fm.create(filenamebase = fileStoreAssoPerm, 
                          nrow = (nVar * nPerm), ncol = (2 * nVar))
        
        if (verbose) {
          message("New files '", 
                  paste0(fileStoreAssoPerm, ".bmat and "),
                  paste0(fileStoreAssoPerm, ".desc.txt created. "))
        }
        
        startIndex <- 0
        
      } else {
        fmat.tmp <- as.matrix(fmat)
        close(fmat)
        
        fmat <- fm.create(filenamebase = fileStoreAssoPerm, 
                          nrow = nrow(fmat.tmp) + (nVar * nPerm), 
                          ncol = (2 * nVar))
        
        fmat[1:nrow(fmat.tmp), 1:(2 * nVar)] <- fmat.tmp
        
        startIndex <- nrow(fmat.tmp)
      }
      
    } else {
      fmat <- fm.create(filenamebase = fileStoreAssoPerm, 
                        nrow = (nVar * nPerm), ncol = (2 * nVar))
      startIndex <- 0
      
      if (verbose) {
        message("Files '", 
                paste0(fileStoreAssoPerm, ".bmat and "),
                paste0(fileStoreAssoPerm, ".desc.txt created. "))
      }
    }
    
    close(fmat)
    
    if (storeCountsPerm) {
      stopifnot(is.character(fileStoreCountsPerm))
      stopifnot(length(fileStoreCountsPerm) == 2)
      
      if (append) {
        fmat_counts1 <- 
          suppressWarnings(try(fm.open(filenamebase = fileStoreCountsPerm[1]), 
                               silent = TRUE))
        
        fmat_counts2 <- 
          suppressWarnings(try(fm.open(filenamebase = fileStoreCountsPerm[2]), 
                               silent = TRUE))
        
        if (inherits(fmat_counts1, "try-error") ||  
            inherits(fmat_counts2, "try-error")) {
          fmat_counts1 <- fm.create(filenamebase = fileStoreCountsPerm[1], 
                                    nrow = (n1 * nPerm), ncol = nVar)
          
          fmat_counts2 <- fm.create(filenamebase = fileStoreCountsPerm[2], 
                                    nrow = (n2 * nPerm), ncol = nVar)
          
          if (verbose) {
            message("New files '", 
                    paste0(fileStoreCountsPerm[1], ".bmat, "),
                    paste0(fileStoreCountsPerm[1], ".desc.txt, "),
                    paste0(fileStoreCountsPerm[2], ".bmat, and "),
                    paste0(fileStoreCountsPerm[2], ".desc.txt created. "))
          }
          
          startIndex_counts1 <- startIndex_counts2 <- 0
          
        } else {
          
          fmat_counts1.tmp <- as.matrix(fmat_counts1)
          fmat_counts2.tmp <- as.matrix(fmat_counts2)
          close(fmat_counts1)
          close(fmat_counts2)
          
          fmat_counts1 <- fm.create(filenamebase = fileStoreCountsPerm[1], 
                                    nrow = nrow(fmat_counts1.tmp) + (n1 * nPerm), 
                                    ncol = nVar)
          
          fmat_counts2 <- fm.create(filenamebase = fileStoreCountsPerm[2], 
                                    nrow = nrow(fmat_counts2.tmp) + (n2 * nPerm), 
                                    ncol = nVar)
          
          fmat_counts1[1:nrow(fmat_counts1.tmp), 1:nVar] <- fmat_counts1.tmp
          fmat_counts2[1:nrow(fmat_counts2.tmp), 1:nVar] <- fmat_counts2.tmp
          
          startIndex_counts1 <- nrow(fmat_counts1.tmp)
          startIndex_counts2 <- nrow(fmat_counts2.tmp)
        }
        
      } else {
        fmat_counts1 <- fm.create(filenamebase = fileStoreCountsPerm[1], 
                                  nrow = (n1 * nPerm), ncol = nVar)
        
        fmat_counts2 <- fm.create(filenamebase = fileStoreCountsPerm[2], 
                                  nrow = (n2 * nPerm), ncol = nVar)
        
        startIndex_counts1 <- startIndex_counts2 <- 0
        
        if (verbose) {
          message("Files '", 
                  paste0(fileStoreCountsPerm[1], ".bmat, "),
                  paste0(fileStoreCountsPerm[1], ".desc.txt, "),
                  paste0(fileStoreCountsPerm[2], ".bmat, and "),
                  paste0(fileStoreCountsPerm[2], ".desc.txt created. "))
        }
      }
      
      close(fmat_counts1)
      close(fmat_counts2)
    }
    
    if (!is.null(seed)) {
      seeds <- sample.int(1e8, size = nPerm)
    } else {
      seeds <- NULL
    }
    
    if (cores > 1) {
      if (parallel::detectCores() < cores) cores <- parallel::detectCores()
      
      if (verbose) {
        message("Starting socket cluster to compute permutation ", 
                "associations in parallel ... ", 
                appendLF = FALSE)
      }
      
      cl <- parallel::makeCluster(cores, outfile = "")
      registerDoSNOW(cl)
      '%do_or_dopar%' <- get('%dopar%')
      
      if (verbose) {
        message("Done.")
      }
      
    } else {
      if (verbose) {
        message("Compute permutation associations ... ")
      }
      '%do_or_dopar%' <- get('%do%')
    }
    
    if (verbose) {
      pb<-txtProgressBar(0, nPerm, style=3)
      
      progress<-function(n) {
        setTxtProgressBar(pb,n)
      }
      
      opts <- list(progress=progress)
    } else {
      opts <- list()
    }
    
    if (!is.null(logFile)) cat("", file=logFile, append=FALSE)
    
    p <- NULL
    propsPerm <- foreach(
      p = 1:nPerm,
      .packages = c(
        "WGCNA",
        "SPRING",
        "SpiecEasi",
        "ccrepe",
        "vegan",
        "LaplacesDemon",
        "propr",
        "DESeq2",
        "filematrix"
      ),
      .export = c(".calcAssociation", "cclasso", "gcoda"),
      .options.snow = opts
    ) %do_or_dopar% {
      
      if (!is.null(logFile)) {
        cat(paste("iteration", p,"\n"),
            file=logFile, append=TRUE)
      }
      
      if (verbose) progress(p)
      
      if (!is.null(seed)) set.seed(seeds[p])
      
      if (assoType == "dissimilarity") {
        count1.tmp <- xbind[ ,which(perm_group_mat[p, ] == 1)]
        count2.tmp <- xbind[ ,which(perm_group_mat[p, ] == 2)]
      } else {
        count1.tmp <- xbind[which(perm_group_mat[p, ] == 1), ]
        count2.tmp <- xbind[which(perm_group_mat[p, ] == 2), ]
      }
      
      if (!parNC$jointPrepro) {
        # zero treatment and normalization necessary if
        # in network construction two count matrices 
        # were given or dissimilarity network is created
        
        suppressMessages(
          count1.tmp <- .zeroTreat(countMat = count1.tmp, 
                                   zeroMethod = parNC$zeroMethod,
                                   zeroParam = parNC$zeroPar, 
                                   needfrac = parNC$needfrac,
                                   needint = parNC$needint, 
                                   verbose = FALSE))
        
        suppressMessages(
          count2.tmp <- .zeroTreat(countMat = count2.tmp, 
                                   zeroMethod = parNC$zeroMethod,
                                   zeroParam = parNC$zeroPar, 
                                   needfrac = parNC$needfrac,
                                   needint = parNC$needint, 
                                   verbose = FALSE))
        
        suppressMessages(
          count1.tmp <- .normCounts(countMat = count1.tmp, 
                                    normMethod = parNC$normMethod,
                                    normParam = parNC$normPar, 
                                    zeroMethod = parNC$zeroMethod,
                                    needfrac = parNC$needfrac, 
                                    verbose = FALSE))
        
        suppressMessages(
          count2.tmp <- .normCounts(countMat = count2.tmp, 
                                    normMethod = parNC$normMethod,
                                    normParam = parNC$normPar, 
                                    zeroMethod = parNC$zeroMethod,
                                    needfrac = parNC$needfrac, 
                                    verbose = FALSE))
      }
      
      if (storeCountsPerm) {
        fmat_counts1 <- fm.open(filenamebase = fileStoreCountsPerm[1])
        fmat_counts2 <- fm.open(filenamebase = fileStoreCountsPerm[2])
        
        fmat_counts1[startIndex_counts1 + (p-1) * n1 + (1:n1), 
                     1:nVar] <- count1.tmp
        
        fmat_counts2[startIndex_counts2 + (p-1) * n2 + (1:n2), 
                     1:nVar] <- count2.tmp
        
        close(fmat_counts1)
        close(fmat_counts2)
      }
      
      assoMat1.tmp <- .calcAssociation(count1.tmp,
                                       measure = parNC$measure,
                                       measurePar = parNC$measurePar,
                                       verbose = FALSE)
      assoMat2.tmp <- .calcAssociation(count2.tmp,
                                       measure = parNC$measure,
                                       measurePar = parNC$measurePar,
                                       verbose = FALSE)
      
      fmat <- fm.open(filenamebase = fileStoreAssoPerm)
      
      fmat[startIndex + (p-1) * nVar + (1:nVar), 
           1:nVar] <- assoMat1.tmp
      fmat[startIndex + (p-1) * nVar + (1:nVar), 
           nVar + (1:nVar)] <- assoMat2.tmp
      
      close(fmat)
    }
    
    if (verbose) {
      close(pb)
      message("Done.")
    }
    
    
    if (cores > 1) {
      if (verbose) {
        message("Stopping socket cluster ... ", appendLF = FALSE)
      }
      parallel::stopCluster(cl)
      if (verbose) {
        message("Done.")
      }
    } 
    
  }
  
  invisible(perm_group_mat)
}
