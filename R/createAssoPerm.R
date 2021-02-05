#' @title Create and store association matrices for permuted data
#'
#' @description The function creates and returns a matrix with permuted group 
#'   labels and stores association matrices computed for the permuted data in 
#'   an external file.
#'
#' @param x object of class \code{"microNetProps"} (inheriting from a call to
#'   \code{\link[NetCoMi]{netAnalyze}}).
#' @param computeAsso logical indicating whether the association matrices should
#'   be computed. If \code{FALSE}, only the permuted group labels are computed 
#'   and returned.
#' @param nPerm number of permutations.
#' @param cores integer indicating the number of CPU cores used for
#'   permutation tests. If cores > 1, the tests are performed parallel.
#'   Is limited to the number of available CPU cores determined by
#'   \code{\link[parallel]{detectCores}}. Defaults to 1L (no parallelization).
#' @param seed integer giving a seed for reproducibility of the results.
#' @param permGroupMat an optional matrix with permuted group labels 
#' (with nPerm rows and n1+n2 columns). 
#' @param fileStoreAssoPerm character giving the file name to store a matrix
#'   containing a matrix with associations/dissimilarities for the permuted 
#'   data. Can also be a path.
#' @param append logical indicating whether existing files (given by 
#'   fileStoreAssoPerm and fileStoreCountsPerm) should be extended. If \code{TRUE}, a new
#'   file is only created if the file is not existing. If \code{FALSE}, a new 
#'   file is created in any case.
#' @param storeCountsPerm logical indicating whether the permuted count matrices
#'   should be stored in an external file. Defaults to \code{FALSE}.
#' @param fileStoreCountsPerm character vector with two elements giving the 
#'   names of two files storing the permuted count matrices belonging to the 
#'   two groups.
#' @param logFile character string naming the log file within which the current
#'   iteration number is stored (if permutation tests are performed). Defaults
#'   to \code{NULL} so that no log file is generated.
#' @param verbose logical. If \code{TRUE} (default), status messages are shown.
#' @return Invisible object: Matrix with permuted group labels.
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
#'                             zeroMethod = "pseudo", normMethod = "clr")
#' 
#'   # Network analysis:
#'   amgut_props <- netAnalyze(amgut_net, clustMethod = "cluster_fast_greedy")
#' 
#'   # Use 'createAssoPerm' to create "permuted" count and association matrices,
#'   # which can be reused by netCompare() and diffNet()
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
#'                             fileLoadCountsPerm = c("countsPerm1", "countsPerm2"),
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
#' @importFrom snow makeCluster stopCluster
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
                           fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
                           logFile = NULL, 
                           verbose = TRUE){

  stopifnot(class(x) =="microNetProps")

  nPerm <- as.integer(nPerm)
  
  cores <- as.integer(cores)
  
  if(!is.null(logFile)) stopifnot(is.character(logFile))
  
  measure <- measurePar <- zeroMethod <- zeroPar <- needfrac <- needint <- 
    normMethod <- normPar <- jointPrepro <- NULL
  
  parnames <- c("measure", "measurePar", "zeroMethod", "zeroPar", "needfrac", 
                "needint", "normMethod", "normPar", "jointPrepro")
  
  for(i in 1:length(parnames)){
    assign(parnames[i], x$paramsNetConstruct[[parnames[i]]])
  }

  assoType <- x$input$assoType
  distNet <- ifelse(assoType == "dissimilarity", TRUE, FALSE)
  
  xgroups <- x$input$groups
  matchDesign <- x$input$matchDesign
  twoNets <- x$input$twoNets
  callNetConstr <- x$input$call
  
  count1 <- x$input$countMat1
  count2 <- x$input$countMat2
  countsJoint <- x$input$countsJoint
  
  count_norm1 <- x$input$normCounts1
  count_norm2 <- x$input$normCounts2
  
  if(!twoNets){
    stop("Input contains only a single network.")
  }
  
  #---------------------------------------------------------------------------
  if(!is.null(seed)) set.seed(seed)
  
  # create combined data matrix
  if(assoType == "dissimilarity"){
    n1 <- ncol(count1)
    n2 <- ncol(count2)
    n <- n1 + n2
    nVar <- nrow(count1)
    
    xbind <- cbind(count1, count2)
    
  } else{
    if(jointPrepro){
      n1 <- nrow(count_norm1)
      n2 <- nrow(count_norm2)
      nVar <- ncol(count_norm1)
      
      xbind <- rbind(count_norm1, count_norm2)
      
    } else{
      n1 <- nrow(count1)
      n2 <- nrow(count2)
      nVar <- ncol(count1)
      
      xbind <- rbind(count1, count2)
    }
    
    n <- n1 + n2
  }

  
  #---------------------------------------
  # create matrix with permuted group labels
  
  if(is.null(permGroupMat)){
    
    if(verbose){
      message("Create matrix with permuted group labels ... ", appendLF = FALSE)
    }
    
    perm_group_mat <- get_perm_group_mat(n1 = n1, n2 = n2, n = n, nPerm = nPerm, 
                                         matchDesign = matchDesign)
    
    if(verbose){
      message("Done.")
    }
    
  } else{
    stopifnot(ncol(permGroupMat) == n)
    perm_group_mat <- permGroupMat
  }
  
  #---------------------------------------
  # compute and store association matrices
  
  if(computeAsso){
    
    stopifnot(is.character(fileStoreAssoPerm))
    stopifnot(length(fileStoreAssoPerm) == 1)
    
    if(append){
      fmat <- suppressWarnings(try(fm.open(filenamebase = fileStoreAssoPerm), 
                                   silent = TRUE))
      
      if(class(fmat) == "try-error"){
        fmat <- fm.create(filenamebase = fileStoreAssoPerm, 
                          nrow = (nVar * nPerm), ncol = (2 * nVar))
        
        if(verbose){
          message("New files '", 
                  paste0(fileStoreAssoPerm, ".bmat and "),
                  paste0(fileStoreAssoPerm, ".desc.txt created. "))
        }
        
        startIndex <- 0
        
      } else{
        fmat.tmp <- as.matrix(fmat)
        close(fmat)
        
        fmat <- fm.create(filenamebase = fileStoreAssoPerm, 
                          nrow = nrow(fmat.tmp) + (nVar * nPerm), 
                          ncol = (2 * nVar))
        
        fmat[1:nrow(fmat.tmp), 1:(2 * nVar)] <- fmat.tmp
        
        startIndex <- nrow(fmat.tmp)
      }
      
    } else{
      fmat <- fm.create(filenamebase = fileStoreAssoPerm, 
                        nrow = (nVar * nPerm), ncol = (2 * nVar))
      startIndex <- 0
      
      if(verbose){
        message("Files '", 
                paste0(fileStoreAssoPerm, ".bmat and "),
                paste0(fileStoreAssoPerm, ".desc.txt created. "))
      }
    }
    
    close(fmat)
    
    if(storeCountsPerm){
      stopifnot(is.character(fileStoreCountsPerm))
      stopifnot(length(fileStoreCountsPerm) == 2)

      if(append){
        fmat_counts1 <- suppressWarnings(try(fm.open(filenamebase = fileStoreCountsPerm[1]), 
                                             silent = TRUE))
        
        fmat_counts2 <- suppressWarnings(try(fm.open(filenamebase = fileStoreCountsPerm[2]), 
                                             silent = TRUE))
        
        if(class(fmat_counts1) == "try-error" || class(fmat_counts2) == "try-error"){
          fmat_counts1 <- fm.create(filenamebase = fileStoreCountsPerm[1], 
                                    nrow = (n1 * nPerm), ncol = nVar)
          
          fmat_counts2 <- fm.create(filenamebase = fileStoreCountsPerm[2], 
                                    nrow = (n2 * nPerm), ncol = nVar)
          
          if(verbose){
            message("New files '", 
                    paste0(fileStoreCountsPerm[1], ".bmat, "),
                    paste0(fileStoreCountsPerm[1], ".desc.txt, "),
                    paste0(fileStoreCountsPerm[2], ".bmat, and "),
                    paste0(fileStoreCountsPerm[2], ".desc.txt created. "))
          }
          
          startIndex_counts1 <- startIndex_counts2 <- 0
          
        } else{

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
        
      } else{
        fmat_counts1 <- fm.create(filenamebase = fileStoreCountsPerm[1], 
                                  nrow = (n1 * nPerm), ncol = nVar)
        
        fmat_counts2 <- fm.create(filenamebase = fileStoreCountsPerm[2], 
                                  nrow = (n2 * nPerm), ncol = nVar)
        
        startIndex_counts1 <- startIndex_counts2 <- 0
        
        if(verbose){
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
    
    if(!is.null(seed)){
      seeds <- sample.int(1e8, size = nPerm)
    } else{
      seeds <- NULL
    }
    
    if(cores > 1){
      if(parallel::detectCores() < cores) cores <- parallel::detectCores()
      
      if(verbose){
        message("Starting socket cluster to compute permutation associations in parallel ... ", 
                appendLF = FALSE)
      }
      cl <- makeCluster(cores, outfile = "")
      registerDoSNOW(cl)
      '%do_or_dopar%' <- get('%dopar%')
      
      if(verbose){
        message("Done.")
      }
      
    } else{
      if(verbose){
        message("Compute permutation associations ... ")
      }
      '%do_or_dopar%' <- get('%do%')
    }
    
    if(verbose){
      pb<-txtProgressBar(0, nPerm, style=3)
      
      progress<-function(n){
        setTxtProgressBar(pb,n)
      }
      
      opts <- list(progress=progress)
    } else{
      opts <- list()
    }
    
    if(!is.null(logFile)) cat("", file=logFile, append=FALSE)
    
    p <- NULL
    propsPerm <- foreach(p = 1:nPerm,
                         .packages = c("WGCNA", "SPRING",
                                       "SpiecEasi",
                                       "ccrepe",
                                       "vegan",
                                       "LaplacesDemon",
                                       "propr",
                                       "DESeq2",
                                       "filematrix"),
                         .export = c("calc_association", "cclasso", "gcoda"),
                         .options.snow = opts) %do_or_dopar% {
                         
                           if(!is.null(logFile)){
                             cat(paste("iteration", p,"\n"),
                                 file=logFile, append=TRUE)
                           }
                           
                           if(verbose) progress(p)
                           
                           if(!is.null(seed)) set.seed(seeds[p])
                           
                           if(assoType == "dissimilarity"){
                             count1.tmp <- xbind[ ,which(perm_group_mat[p, ] == 1)]
                             count2.tmp <- xbind[ ,which(perm_group_mat[p, ] == 2)]
                           } else{
                             count1.tmp <- xbind[which(perm_group_mat[p, ] == 1), ]
                             count2.tmp <- xbind[which(perm_group_mat[p, ] == 2), ]
                           }
                           
                           if(!jointPrepro){
                             # zero treatment and normalization necessary if
                             # in network construction two count matrices 
                             # were given or dissimilarity network is created
                             
                             suppressMessages(
                               count1.tmp <- zero_treat(countMat = count1.tmp, 
                                                        zeroMethod = zeroMethod,
                                                        zeroParam = zeroPar, 
                                                        needfrac = needfrac,
                                                        needint = needint, 
                                                        verbose = FALSE))
                             
                             suppressMessages(
                               count2.tmp <- zero_treat(countMat = count2.tmp, 
                                                        zeroMethod = zeroMethod,
                                                        zeroParam = zeroPar, 
                                                        needfrac = needfrac,
                                                        needint = needint, 
                                                        verbose = FALSE))
                             
                             suppressMessages(
                               count1.tmp <- norm_counts(countMat = count1.tmp, 
                                                         normMethod = normMethod,
                                                         normParam = normPar, 
                                                         zeroMethod = zeroMethod,
                                                         needfrac = needfrac, 
                                                         verbose = FALSE))
                             
                             suppressMessages(
                               count2.tmp <- norm_counts(countMat = count2.tmp, 
                                                         normMethod = normMethod,
                                                         normParam = normPar, 
                                                         zeroMethod = zeroMethod,
                                                         needfrac = needfrac, 
                                                         verbose = FALSE))
                           }
                           
                           if(storeCountsPerm){
                             fmat_counts1 <- fm.open(filenamebase = fileStoreCountsPerm[1])
                             fmat_counts2 <- fm.open(filenamebase = fileStoreCountsPerm[2])
                             
                             fmat_counts1[startIndex_counts1 + (p-1) * n1 + (1:n1), 
                                          1:nVar] <- count1.tmp
                             
                             fmat_counts2[startIndex_counts2 + (p-1) * n2 + (1:n2), 
                                          1:nVar] <- count2.tmp
                             
                             close(fmat_counts1)
                             close(fmat_counts2)
                           }
                           
                           assoMat1.tmp <- calc_association(count1.tmp,
                                                            measure = measure,
                                                            measurePar = measurePar,
                                                            verbose = FALSE)
                           assoMat2.tmp <- calc_association(count2.tmp,
                                                            measure = measure,
                                                            measurePar = measurePar,
                                                            verbose = FALSE)
                           
                           fmat <- fm.open(filenamebase = fileStoreAssoPerm)
                           
                           fmat[startIndex + (p-1) * nVar + (1:nVar), 
                                1:nVar] <- assoMat1.tmp
                           fmat[startIndex + (p-1) * nVar + (1:nVar), 
                                nVar + (1:nVar)] <- assoMat2.tmp
                           
                           close(fmat)
                         }

    if(verbose){
      close(pb)
      message("Done.")
    }
    
    
    if(cores > 1){
      if(verbose){
        message("Stopping socket cluster ... ", appendLF = FALSE)
      }
      snow::stopCluster(cl)
      if(verbose){
        message("Done.")
      }
    } 
    
  }
  
  invisible(perm_group_mat)
}