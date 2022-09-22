#' @title CCLasso: Correlation inference of Composition data through Lasso method
#'
#' @description Implementation of the CCLasso approach 
#'   \cite{(Fang et al., 2015)}, which is published on GitHub 
#'   \cite{(Fang, 2016)}. The function is extended by a progress message.
#'
#' @param x numeric matrix (\emph{n}x\emph{p}) with samples in rows and 
#'   OTUs/taxa in columns.
#' @param counts logical indicating whether x constains counts or fractions.
#'   Defaults to \code{FALSE} meaning that x contains fractions so that rows
#'   sum up to 1.
#' @param pseudo numeric value giving a pseudo count, which is added to all
#'   counts if \code{counts = TRUE}. Default is 0.5.
#' @param sig numeric matrix giving an initial covariance matrix. If \code{NULL}
#'  (default), \code{diag(rep(1, p))} is used.
#' @param lams numeric vector specifying the tuning parameter sequences. Default
#'   is \code{10^(seq(0, -8, by = -0.01))}.
#' @param K numeric value (integer) giving the folds of crossvalidation.
#'   Defaults to 3.
#' @param kmax numeric value (integer) specifying the maximum iteration for
#'   augmented lagrangian method. Default is 5000.
#' @param verbose logical indicating whether a progress indicator is shown 
#' (\code{TRUE} by default).
#'
#' @return A list containing the following elements:
#' \tabular{ll}{
#'   \code{cov.w}\tab Covariance estimation\cr
#'   \code{cor.w}\tab Correlation estimation\cr
#'   \code{lam}\tab Final tuning parameter}
#'
#' @author
#'   Fang Huaying, Peking University (R code)\cr
#'   Stefanie Peschel (documentation)
#' @references
#'   \insertRef{fang2015cclasso}{NetCoMi}\cr\cr
#'   \insertRef{fang2016cclassoGithub}{NetCoMi}
#' @import Matrix
#' @export

cclasso <- function(x, counts = F, pseudo = 0.5, sig = NULL,
                    lams = 10^(seq(0, -8, by = -0.01)),
                    K = 3, kmax = 5000, verbose = TRUE) {
  
  # data dimension
  p <- ncol(x);
  n <- nrow(x);
  # Counts or Fractions?
  if (counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  # log transformation
  xlog <- log(x);
  # use all data
  vx <- stats::var(xlog);
  # initial value
  res <- list();
  if (is.null(sig)) {
    res$sig <- diag(rep(1, p));
  }
  else {
    res$sig <- sig;
  }
  #-----------------------------------------------------------------------------
  # weight diagonal for loss
  rmean.vx <- rowMeans(vx);
  wd <- 1 / diag(vx - rmean.vx - rep(rmean.vx, each = p) + mean(rmean.vx));
  wd2 <- sqrt(wd);
  #-----------------------------------------------------------------------------
  # preparation for update sigma in augmented lagrange method
  rho <- 1; # needed
  u.f <- eigen(diag(rep(1, p)) - 1 / p)$vectors; # needed
  wd.u <- (t(u.f) %*% (wd * u.f))[-p, -p];
  diag(wd.u) <- diag(wd.u) + rho;
  wd.u.eig <- eigen(wd.u);
  d0.wd <- 2 / outer(wd.u.eig$values, wd.u.eig$values, "+"); # needed
  u0.wd <- wd.u.eig$vectors; # needed
  #-----------------------------------------------------------------------------
  n_lam <- length(lams);
  tol.zero <- 1e-8;
  
  #-----------------------------------------------------------------------------
  # cross validation
  if (n_lam == 1) {
    lamA <- lams[1];
  }
  else {
    tol.loss <- 1e-6;
    loss.old <- Inf;
    k.loss <- n_lam;
    n.b <- floor(n / K);
    # loss <- rep(0, n_lam);
    for (i in 1:n_lam) {
      
      loss.cur <- 0;
      for (k in 1:K) {
        # testing data and training data
        itest <- (n.b * (k-1) + 1):(n.b * k);
        vxk <- stats::var(xlog[itest, ]);
        #
        vx2k <- stats::var(xlog[-itest, ]);
        # for training
        res <- cclasso.sub(vx = vx2k, wd = wd, lam = lams[i],
                           u.f = u.f, u0.wd = u0.wd, d0.wd = d0.wd,
                           sig = res$sig, rho = rho, kmax = kmax);
        # loss cumulation
        res$sig[abs(res$sig) <= tol.zero] <- 0;
        dsig <- res$sig - vxk;
        rmean.dsig <- rowMeans(dsig);
        half.loss <- dsig * rep(wd2, each = p) -
          outer(rmean.dsig, wd2, "*") +
          rep(  (mean(rmean.dsig) - rmean.dsig) * wd2, each = p);
        # loss[i] <- loss[i] + base::norm(half.loss, "F")^2;
        loss.cur <- loss.cur + base::norm(half.loss, "F")^2;
      }
      
      thresh <- tol.loss * max(loss.cur, loss.old, 1)
      if (verbose) {
        
        message("current loss diff.: ", 
                sprintf("%.6f", round(loss.cur - loss.old, 6)), 
                " (breaks if >= ", sprintf("%.6f", round(thresh, 6)),
                ")\r", appendLF=FALSE)
      }
      
      if (loss.cur - loss.old >= thresh) {
        k.loss <- i - 1;
        break;
      }
      else {
        loss.old <- loss.cur;
      }
    }
    # select lambda
    # k.loss <- which.min(loss);
    lamA <- lams[k.loss];
    if (k.loss == 1 || k.loss == n_lam) {
      cat("Warning:", "Tuning (", lamA ,") on boundary!\n");
    }
  }
  
  if (verbose) message("")
  #-----------------------------------------------------------------------------
  
  res <- cclasso.sub(vx = vx, wd = wd, lam = lamA,
                     u.f = u.f, u0.wd = u0.wd, d0.wd = d0.wd,
                     sig = res$sig, rho = rho, kmax = kmax);

  res$sig[abs(res$sig) <= tol.zero] <- 0
  
  #---------
  # Edited by Stefanie Peschel:
  eigres <- eigen(res$sig)$values
  
  if (typeof(eigres) == "complex") {
    eigres <- Re(eigres)
  }
  #---------
    
  if (min(eigres) <= tol.zero) {
    sig.sparse <- abs(res$sig) > tol.zero
    diag(res$sig) <- diag(res$sig) * sign(diag(res$sig))
    res$sig <- as.matrix(Matrix::nearPD(res$sig)$mat) * sig.sparse
  }
  
  # get correlation matrix from covariance matrix
  Is <- sqrt(1 / diag(res$sig));
  cor.w <- Is * res$sig * rep(Is, each = p);
  # remove too small correlation values
  cor.w[abs(cor.w) <= 1e-6] <- 0;
  #
  return(list(cov.w = res$sig, cor.w = cor.w, lam = lamA));
}


# cclasso for only one lambda
cclasso.sub <- function(vx, wd, lam, u.f, u0.wd, d0.wd, sig = NULL,
                        rho = 1, kmax = 5000, x.tol = 1e-6) {
  p <- ncol(vx);
  # initial value
  lam.rho <- lam / rho;
  if (is.null(sig)) {
    sig <- diag(rep(1, p));
  }
  sig2 <- sig;
  LAM <- matrix(0, p, p);
  # loop start
  k <- 0;
  err <- 1;
  while(err > x.tol && k < kmax) {
    # update sigma
    x.sig <- t(u.f) %*% ((sig2  -  vx) - LAM  / rho) %*% u.f;
    x.sig[-p,-p] <- u0.wd %*% ((t(u0.wd) %*% x.sig[-p, -p] %*% u0.wd) *
                                 d0.wd * rho) %*% t(u0.wd);
    sig.new <- vx + u.f %*% x.sig %*% t(u.f);
    # update sigma2
    A <- LAM / rho + sig.new;
    sig2.new <- (A > lam.rho) * (A - lam.rho) +
      (A < -lam.rho) * (A + lam.rho);
    diag(sig2.new) <- diag(A);
    # update Lambda
    LAM <- LAM + rho * (sig.new - sig2.new);
    # calculate error
    err <- max( base::norm(sig.new - sig, "F") / max(1, base::norm(sig)),
                base::norm(sig2.new - sig2, "F") / max(1, base::norm(sig2)));
    # update current value
    sig <- sig.new;
    sig2 <- sig2.new;
    k <- k + 1;
  }
  #
  if (k >= kmax) {
    cat("Warning:", "Maximum ", kmax,
        "iteration while relative error is", err, "\n");
  }
  #
  return(list(sig = sig, k = k));
}


