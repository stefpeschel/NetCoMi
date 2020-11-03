#' @title gCoda: conditional dependence network inference for compositional data
#'
#' @description A parallelized implementation of the gCoda approach 
#'   \cite{(Fang et al., 2017)}, published on GitHub \cite{(Fang, 2016)}.
#'
#' @param x numeric matrix (\emph{n}x\emph{p}) with samples in rows and OTUs/taxa in
#'   columns.
#' @param counts logical indicating whether x constains counts or fractions.
#'   Defaults to \code{FALSE} meaning that x contains fractions so that rows
#'   sum up to 1.
#' @param pseudo numeric value giving a pseudo count, which is added to all
#'   counts if \code{counts = TRUE}. Default is 0.5.
#' @param lambda.min.ratio numeric value specifying lambda(max) / lambda(min).
#'   Defaults to 1e-4.
#' @param nlambda numberic value (integer) giving the of tuning parameters.
#'   Defaults to 15.
#' @param ebic.gamma numeric value specifying the gamma value of EBIC.
#'   Defaults to 0.5.
#' @param cores integer indicating the number of CPU cores used for computation. 
#'   Defaults to 1L. For \code{cores} > 1L, \code{\link[foreach]{foreach}} is 
#'   used for parallel execution.
#' @param verbose logical indicating whether a progress indicator is shown 
#' (\code{TRUE} by default). 
#'
#' @return A list containing the following elements:
#' \tabular{ll}{
#'   \code{lambda}\tab lambda sequence for compuation of EBIC score\cr
#'   \code{nloglik}\tab negative log likelihood for lambda sequence\cr
#'   \code{df}\tab number of edges for lambda sequence\cr
#'   \code{path}\tab sparse pattern for lambda sequence\cr
#'   \code{icov}\tab inverse covariance matrix for lambda sequence\cr
#'   \code{ebic.score}\tab EBIC score for lambda sequence\cr
#'   \code{refit}\tab sparse pattern with best EBIC score\cr
#'   \code{opt.icov}\tab inverse covariance matrix with best EBIC score\cr
#'   \code{opt.lambda}\tab lambda with best EBIC score}
#'
#' @author
#'   Fang Huaying, Peking University (R-Code and documentation)\cr
#'   Stefanie Peschel (Parts of the documentation; Parallelization)
#' @references
#'   \insertRef{fang2016gcodaGithub}{NetCoMi}\cr
#'   \insertRef{fang2017gcoda}{NetCoMi}
#'
#' @export

gcoda <- function(x, counts = F, pseudo = 0.5, lambda.min.ratio = 1e-4,
                  nlambda = 15, ebic.gamma = 0.5, cores = 1L, verbose = TRUE) {

  # Counts or fractions?
  if(counts) {
    x <- x + pseudo
    x <- x / Matrix::rowSums(x)
  }
  n <- nrow(x)
  p <- ncol(x)
  # Log transformation for compositional data
  S <- var(log(x) - Matrix::rowMeans(log(x)))
  # Generate lambda via lambda.min.ratio and nlambda
  lambda.max <- max(max(S - diag(p)), -min(S - diag(p)))
  lambda.min <- lambda.min.ratio * lambda.max
  lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))

  # Compute solution paths for gcoda
  icov <- diag(p)
  if(cores > 1){
    cl <- snow::makeCluster(cores)
    #snow::clusterExport(cl, c("S", "icov", "lambda"), envir = environment())
    doSNOW::registerDoSNOW(cl)
    '%do_or_dopar%' <- get('%dopar%')
    if(verbose) message("Start parallel foreach loop ...")
  } else{
    '%do_or_dopar%' <- get('%do%')
  }
  
  if(verbose){
    pb <- utils::txtProgressBar(0,nlambda,style=3)
    progress <- function(n){
      setTxtProgressBar(pb,n)
    }

    opts <- list(progress=progress)
  } else{
    opts <- list()
  }

  fit.tmp <- foreach(i = 1:nlambda, .export = c("gcoda_sub", "obj_gcoda"), 
                     .noexport = "x", .options.snow = opts) %do_or_dopar% {
                      
                      if(verbose) progress(i)
                      
                      
                      outlist <- list()

                      out.gcoda <- gcoda_sub(A = S, iSig = icov, lambda = lambda[i])
                      icov <- out.gcoda$iSig
                      
                      outlist[["nloglik"]] <- out.gcoda$nloglik
                      outlist[["icov"]] <- icov
                      outlist[["path"]] <- 0 + (abs(icov) > 1e-6)
                      diag(outlist[["path"]]) <- 0
                      outlist[["df"]] <- sum(outlist[["path"]]) / 2
                      
                      outlist
                     }
  
  if(verbose) close(pb)

  if(verbose & cores > 1){
    message("Stopping socket cluster ... ", appendLF = FALSE)
  }
  
  if(cores > 1) snow::stopCluster(cl)
  
  if(verbose & cores > 1){
    message("Done.")
  }
  
  
  # Store fit result for gcoda via a series of lambda
  fit <- list()
  fit$lambda <- lambda
  fit$nloglik <- rep(0, nlambda)
  fit$df <- rep(0, nlambda)
  fit$path <- list()
  fit$icov <- list()
  
  for(i in 1:nlambda){
    fit$nloglik[i] <- fit.tmp[[i]]$nloglik
    fit$df[i] <- fit.tmp[[i]]$df
    fit$path[[i]] <- fit.tmp[[i]]$path
    fit$icov[[i]] <- fit.tmp[[i]]$icov
  }

  #####################################
  # Continue if edge density is too small
  imax <- nlambda + 15
  emin <- p * (p - 1) / 2 * 0.618

  if(fit$df[i] <= emin && i <= imax && lambda[i] > 1e-6 && verbose){
    message("\nEdge density too small. Continuing ...")
  }
  
  while(fit$df[i] <= emin && i <= imax && lambda[i] > 1e-6) {
    if(verbose) message("nlambda = ", i, " ... Going on ...\r",appendLF=FALSE)
    lambda <- c(lambda, lambda[i]/2)
    i <- i + 1
    fit$lambda[i] <- lambda[i]
    out.gcoda <- gcoda_sub(A = S, iSig = icov, lambda = lambda[i])
    icov <- out.gcoda$iSig
    fit$nloglik[i] <- out.gcoda$nloglik
    fit$icov[[i]] <- icov
    fit$path[[i]] <- 0 + (abs(icov) > 1e-6)
    diag(fit$path[[i]]) <- 0
    fit$df[i] <- sum(fit$path[[i]]) / 2
  }
  
  if(verbose) message("")
  # Compute EBIC score for lambda selection
  fit$ebic.score <- n * fit$nloglik + log(n) * fit$df +
    4 * ebic.gamma * log(p) * fit$df
  fit$opt.index <- which.min(fit$ebic.score)
  fit$refit <- fit$path[[fit$opt.index]]
  fit$opt.icov <- fit$icov[[fit$opt.index]]
  fit$opt.lambda <- fit$lambda[fit$opt.index]
  return(fit)
}


#-------------------------------------------------------------------------------
# Optimization for gcoda with given lambda
gcoda_sub <- function(A, iSig = NULL, lambda = 0.1, tol_err = 1e-4,
                      k_max = 100) {
  p <- ncol(A)
  if(is.null(iSig)) {
    iSig <- diag(p)
  }

  err <- 1
  k <- 0
  fval_cur <- Inf

  while(err > tol_err && k < k_max) {
    iSig_O <- Matrix::rowSums(iSig)
    iS_iSig <- 1 / sum(iSig_O)
    iSig_O2 <- iSig_O * iS_iSig
    A_iSig_O2 <- Matrix::rowSums(A * rep(iSig_O2, each = p))
    A2 <- A - A_iSig_O2 - rep(A_iSig_O2, each = p) +
      sum(iSig_O2 * A_iSig_O2) + iS_iSig
    iSig2 <- huge::huge(x = A2, lambda = lambda, method = "glasso", 
                        verbose = FALSE)$icov[[1]]

    fval_new <- obj_gcoda(iSig = iSig2, A = A, lambda = lambda)
    xerr <- max(abs(iSig2 - iSig) / (abs(iSig2) + 1))
    err <- min(xerr, abs(fval_cur - fval_new)/(abs(fval_new) + 1))

    k <- k + 1
    iSig <- iSig2
    fval_cur <- fval_new
  }
  nloglik <- fval_cur - lambda * sum(abs(iSig))

  if(k >= k_max) {
    cat("WARNING of gcoda_sub:\n", "\tMaximum Iteration:", k_max,
        "&& Relative error:", err, "!\n")
  }

  return(list(iSig = iSig, nloglik = nloglik))
}
#----------------------------------------
# Objective function value of gcoda (negative log likelihood + penalty)
obj_gcoda <- function(iSig, A, lambda) {
  p <- ncol(A)
  iSig_O <- Matrix::rowSums(iSig)
  S_iSig <- sum(iSig_O)
  nloglik <- - log(det(iSig)) + sum(iSig * A) + log(S_iSig) -
    sum(iSig_O * Matrix::rowSums(A * rep(iSig_O, each = p))) / S_iSig
  pen <- lambda * sum(abs(iSig))
  return(nloglik + pen)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
