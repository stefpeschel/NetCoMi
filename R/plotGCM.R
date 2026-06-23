#' @title Plot a Graphlet Correlation Matrix Heatmap
#'
#' @description Plot the Graphlet Correlation Matrix (GCM) heatmap from an
#'   object returned by \code{\link{netAnalyze}}, \code{\link{calcGCM}}, or
#'   \code{\link{calcGCD}}.
#'
#' @param x object of class \code{microNetProps} returned by
#'   \code{\link{netAnalyze}}, object of class \code{GCM} returned by
#'   \code{\link{calcGCM}}, or object of class \code{GCD} returned by
#'   \code{\link{calcGCD}}.
#' @param gcmLCC logical. If \code{FALSE} (default), the GCM heatmap is plotted
#'   for the whole network. If \code{TRUE}, the GCM heatmap is plotted for the
#'   largest connected component. Only used if \code{x} is of class
#'   \code{microNetProps}.
#' @param titles optional character vector with plot titles. For one network,
#'   a single title is used. For two networks, three titles are used for the
#'   first GCM, second GCM, and the difference heatmap.
#' @param mar numeric vector specifying the plot margins passed to
#'   \code{\link{plotHeat}}. If \code{NULL}, suitable defaults are used.
#' @param sigTextCex numeric defining the size of the explanatory significance
#'   code text. Default is 1.
#' @param ... further arguments passed to \code{\link{plotHeat}}.
#'
#' @return Invisibly returns a list with the arguments used by
#'   \code{\link{plotHeat}}.
#'
#' @examples
#' data("amgut1.filt")
#'
#' # Single-network workflow
#' amgut_net <- netConstruct(amgut1.filt, measure = "pearson",
#'                           filtTax = "highestVar",
#'                           filtTaxPar = list(highestVar = 50),
#'                           zeroMethod = "pseudoZO", normMethod = "clr",
#'                           sparsMethod = "threshold", thresh = 0.4)
#'
#' amgut_props <- netAnalyze(amgut_net, clustMethod = "cluster_fast_greedy",
#'                           hubPar = "eigenvector")
#'
#' # Object returned by netAnalyze()
#' plotGCM(amgut_props)
#'
#' # Object returned by calcGCM()
#' gcm <- calcGCM(amgut_net$adjaMat1)
#' plotGCM(gcm)
#'
#' # Two-network workflow
#' data("amgut2.filt.phy")
#' 
#' amgut_season_yes <- phyloseq::subset_samples(amgut2.filt.phy,
#'                                       SEASONAL_ALLERGIES == "yes")
#' amgut_season_no <- phyloseq::subset_samples(amgut2.filt.phy,
#'                                      SEASONAL_ALLERGIES == "no")
#' 
#' n_yes <- phyloseq::nsamples(amgut_season_yes)
#' 
#' amgut_net2 <- netConstruct(data = amgut_season_yes,
#'                            data2 = amgut_season_no,
#'                            measure = "pearson",
#'                            filtSamp = "highestFreq",
#'                            filtSampPar = list(highestFreq = n_yes),
#'                            filtTax = "highestVar",
#'                            filtTaxPar = list(highestVar = 30),
#'                            zeroMethod = "pseudoZO", normMethod = "clr")
#' 
#' # Object returned by calcGCD()
#' gcd <- calcGCD(amgut_net2$adjaMat1, amgut_net2$adjaMat2)
#' plotGCM(gcd)
#'
#' @seealso \code{\link{netAnalyze}}, \code{\link{calcGCM}},
#'   \code{\link{calcGCD}},
#'   \code{\link{testGCM}}, \code{\link{plotHeat}}
#'
#' @importFrom graphics layout par text
#' @export

plotGCM <- function(x,
                    gcmLCC = FALSE,
                    titles = NULL,
                    mar = NULL,
                    sigTextCex = 1,
                    ...) {
  if (!(inherits(x, "microNetProps") || inherits(x, "GCM") ||
        inherits(x, "GCD"))) {
    stop("\"x\" must be of class \"microNetProps\" returned by netAnalyze() ",
         "or class \"GCM\" or \"GCD\" returned by calcGCM() or calcGCD().")
  }
  
  if (!is.logical(gcmLCC) || length(gcmLCC) != 1L || is.na(gcmLCC)) {
    stop("\"gcmLCC\" must be a logical value.")
  }
  
  if (inherits(x, "microNetProps") && is.null(x$graphlet)) {
    stop("No graphlet properties found in \"x\". Please run netAnalyze() ",
         "with graphlet = TRUE.")
  }
  
  if (!is.null(titles) && !is.character(titles)) {
    stop("\"titles\" must be a character vector.")
  }
  
  if (!is.null(mar) && !(is.numeric(mar) && length(mar) == 4L)) {
    stop("\"mar\" must be a numeric vector with 4 elements.")
  }
  
  if (!(is.numeric(sigTextCex) && length(sigTextCex) == 1L &&
        sigTextCex > 0)) {
    stop("\"sigTextCex\" must be a positive numeric value.")
  }
  
  graphlet <- .getGCMPlotData(x, gcmLCC)
  twoNets <- isTRUE(graphlet$twoNets)
  out <- list()
  
  if (twoNets) {
    if (is.null(titles)) {
      titles <- c("GCM1", "GCM2", "Differences between GCMs (GCM1-GCM2)")
    } else if (length(titles) != 3L) {
      stop("\"titles\" must have length 3 for two networks.")
    }
    
    if (is.null(mar)) mar <- c(0, 0, 2, 0)
    
    opar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar), add = TRUE)
    
    mat <- matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE)
    graphics::layout(mat = mat)
    
    out$gcm1 <- .plotGCMHeat(mat = graphlet$gcm1,
                             pmat = graphlet$pAdjust1,
                             title = titles[1],
                             colorLim = c(-1, 1),
                             mar = mar,
                             ...)
    .addOrbRect()
    
    out$gcm2 <- .plotGCMHeat(mat = graphlet$gcm2,
                             pmat = graphlet$pAdjust2,
                             title = titles[2],
                             colorLim = c(-1, 1),
                             mar = mar,
                             ...)
    .addOrbRect()
    
    gcmtest <- testGCM(.asGCM(graphlet$gcm1, graphlet$ocount1),
                       .asGCM(graphlet$gcm2, graphlet$ocount2),
                       verbose = FALSE)
    
    out$diff <- .plotGCMHeat(mat = gcmtest$diff,
                             pmat = gcmtest$pAdjustDiff,
                             title = titles[3],
                             mar = mar,
                             ...)
    .addOrbRect()
    
    graphics::text(14, 1,
                   "Significance codes:\n***: 0.001;  **: 0.01;  *: 0.05",
                   xpd = NA, pos = 4, cex = sigTextCex)
  } else {
    if (is.null(titles)) {
      titles <- "GCM"
    } else if (length(titles) != 1L) {
      stop("\"titles\" must have length 1 for one network.")
    }
    
    if (is.null(mar)) mar <- c(2, 0, 2, 0)
    
    out$gcm <- .plotGCMHeat(mat = graphlet$gcm1,
                            pmat = graphlet$pAdjust1,
                            title = titles,
                            colorLim = c(-1, 1),
                            mar = mar,
                            ...)
    .addOrbRect()
    
    graphics::text(6, -0.2,
                   "Significance codes:  ***: 0.001;  **: 0.01;  *: 0.05",
                   xpd = NA, cex = sigTextCex)
  }
  
  invisible(out)
}


#' @keywords internal
.getGCMPlotData <- function(x, gcmLCC) {
  if (inherits(x, "GCM")) {
    gcmtest <- testGCM(x, verbose = FALSE)
    
    return(list(ocount1 = x$ocount,
                gcm1 = x$gcm,
                pAdjust1 = gcmtest$pAdjust1,
                twoNets = FALSE))
  }
  
  if (inherits(x, "GCD")) {
    gcmtest <- testGCM(x, verbose = FALSE)
    
    return(list(ocount1 = x$ocount1,
                ocount2 = x$ocount2,
                gcm1 = x$gcm1,
                gcm2 = x$gcm2,
                pAdjust1 = gcmtest$pAdjust1,
                pAdjust2 = gcmtest$pAdjust2,
                twoNets = TRUE))
  }
  
  graphlet <- if (gcmLCC) x$graphletLCC else x$graphlet
  graphlet$twoNets <- isTRUE(x$input$twoNets)
  graphlet
}


#' @keywords internal
.asGCM <- function(gcm, ocount) {
  out <- list(gcm = gcm, ocount = ocount)
  class(out) <- "GCM"
  out
}


#' @keywords internal
.plotGCMHeat <- function(mat, pmat, title, mar, colorLim = NULL, ...) {
  args <- list(...)
  defaults <- list(mat = mat,
                   pmat = pmat,
                   type = "mixed",
                   title = title,
                   colorLim = colorLim,
                   mar = mar)
  defaults <- defaults[!vapply(defaults, is.null, logical(1))]
  do.call(plotHeat, utils::modifyList(defaults, args))
}
