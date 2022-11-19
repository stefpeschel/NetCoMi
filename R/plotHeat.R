#' @title Create a heatmap with p-values
#' 
#' @description A function to draw heatmaps with the option to use p-values 
#'   or significance codes as cell text. It allows to draw a mixed heatmap with 
#'   different cell text (values, p-values, or significance code) in the lower 
#'   and upper triangle. The function \code{\link[corrplot]{corrplot}} is used 
#'   for plotting the heatmap.
#' 
#' @param mat numeric matrix with values to be plotted.
#' @param pmat optional matrix with p-values.
#' @param type character defining the type of the heatmap. Possible values are:
#'   \describe{
#'   \item{\code{"full"}}{Default. The cell text specified via \code{textUpp} is 
#'   used for the whole heatmap.}
#'   \item{\code{"mixed"}}{Different cell text is used for the upper and lower
#'   triangle. The upper triangle is specified via \code{textUpp} and the lower 
#'   triangle via \code{textLow}.}
#'   \item{\code{"upper"}}{Only the upper triangle is plotted. The text is 
#'   specified via \code{textUpp}.}
#'   \item{\code{"lower"}}{Only the lower triangle is plotted. The text is 
#'   specified via \code{textLow}.}
#'   }
#' @param textUpp character specifying the cell text either for the 
#'   full heatmap (if \code{type} is "full") or for the upper triangle (if 
#'   \code{type} is "mixed" or "upper"). Default is "mat". Possible values are:
#'   \describe{
#'   \item{\code{"mat"}}{Cells contain the values in the matrix given by 
#'   \code{mat}}
#'   \item{\code{"sigmat"}}{Same as "mat" but insignificant values (and cells) 
#'   are blank.}
#'   \item{\code{"pmat"}}{Cells contain the p-values given by \code{p-mat}.}
#'   \item{\code{"code"}}{Cells contain significance codes corresponding to the 
#'   p-values given by \code{p-mat}. The following coding is used: 
#'   "***: 0.001;  **: 0.01;  *: 0.05".}
#'   \item{\code{"none"}}{No cell text is plotted.}
#'   }
#' @param textLow same as \code{textUpp} but for the lower triangle 
#'   (if \code{type} is "mixed" or "lower"). Default is "code".
#' @param methUpp character specifying how values are represented in the 
#'   full heatmap (if \code{type} is "full") or in the upper triangle (if 
#'   \code{type} is "mixed" or "upper"). 
#'   Possible values are: "circle", "square", "ellipse", "number", 
#'   "shade", "color" (default), "pie". The method is passed to the 
#'   \code{method} argument of \code{\link[corrplot]{corrplot}}.
#' @param methLow same es \code{methUpp} but for the lower triangle.
#' @param diag logical. If \code{TRUE} (default), the diagonal is printed. 
#'   If \code{FALSE} and \code{type} is "full" or "mixed", the diagonal cells 
#'   are white. If \code{FALSE} and \code{type} is "upper" or "lower", only the 
#'   non-diagonal cells are printed.
#' @param title character giving the title.
#' @param mar vector specifying the plot margins. 
#'   See \code{\link[graphics]{par}}. Default is c(0, 0, 1, 0).
#' @param labPos character defining the label position. Possible values are:
#'   "lt"(left and top, default), 
#'   "ld"(left and diagonal; \code{type} must be "lower"), 
#'   "td"(top and diagonal; \code{type} must be "upper"), 
#'   "d"(diagonal only), "n"(no labels). Passed to 
#'   \code{\link[corrplot]{corrplot}} argument \code{tl.pos}.
#' @param labCol label color. Default is "gray40". Passed to 
#'   \code{\link[corrplot]{corrplot}} argument \code{tl.col}.
#' @param labCex numeric defining the label size. Default is 1.1. Passed to 
#'   \code{\link[corrplot]{corrplot}} argument \code{tl.cex}.
#' @param textCol color of the cell text (values, p-values, and code). Default 
#'   is "black".
#' @param textCex numeric defining the text size. Default is 1. Currently only 
#'   works for types "mat" and "code".
#' @param textFont numeric defining the text font. Default is 1. Currently only 
#'   works for type "mat". 
#' @param digits integer defining the number of decimal places used for 
#'   matrix values and p-values.
#' @param legendPos position of the color legend. Possible values are: 
#'   "r"(right; default), "b"(bottom), "n"(no legend).
#' @param colorPal character specifying the color palette used for cell coloring 
#'   if \code{color} is not set. Available are the sequential and diverging 
#'   color palettes from \code{\link[RColorBrewer]{RColorBrewer}}:
#'   \describe{
#'   \item{Sequential:}{"Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", 
#'   "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", 
#'   "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"}
#'   \item{Diverging:}{"BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", 
#'   "RdYlGn", "Spectral"}
#'   }
#'   By default, "RdBu" is used if the first value of \code{colorLim} is 
#'   negative and "YlOrRd" otherwise.
#' @param addWhite logical. If \code{TRUE}, white is added to the color palette.
#'   (first element for sequential palettes and middle element for diverging 
#'   palettes). For a diverging palette, \code{nCol} should be set to an odd 
#'   number so that the middle color is white. 
#' @param nCol integer defining the number of colors to which the color palette 
#'   should be interpolated. Default is 51L. 
#'   \code{\link[grDevices]{colorRamp}} is used for color interpolation.
#' @param colorLim numeric vector with two values defining the color limits. 
#'   The first element of the color vector is assigned to the lower limit and 
#'   the last element of the color vector to the upper limit. Default is 
#'   c(0,1) if the values of \code{mat} are in [0,1], c(-1,1) if the values are 
#'   in [-1,1], and the minimum and maximum values otherwise.
#' @param revCol logical. If \code{TRUE}, the reversed color vector is used. 
#'   Default is \code{FALSE}. Ignored if \code{color} is given.
#' @param color an optional vector with colors used for cell coloring.
#' @param bg background color of the cells. Default is "white".
#' @param argsUpp optional list of arguments passed to 
#'   \code{\link[corrplot]{corrplot}}. Arguments set within \code{plotHeat()} 
#'   are overwritten by arguments in the list. Used for the 
#'   full heatmap if \code{type} is "full" and for the upper triangle if 
#'   \code{type} is "mixed" or "upper".
#' @param argsLow same as \code{argsUpp} but for the lower triangle 
#'   (if \code{type} is "mixed" or "lower").
#' 
#' @return Invisible list with two elements \code{argsUpper} and 
#'   \code{argsLower} containing the \code{\link[corrplot]{corrplot}} 
#'   arguments used for the upper and lower triangle of the heatmap.
#'   
#' @examples 
#' # Load data sets from American Gut Project (from SpiecEasi package)
#' data("amgut2.filt.phy")
#' 
#' # Split data into two groups: with and without seasonal allergies
#' amgut_season_yes <- phyloseq::subset_samples(amgut2.filt.phy, 
#'                                       SEASONAL_ALLERGIES == "yes")
#' amgut_season_no <- phyloseq::subset_samples(amgut2.filt.phy, 
#'                                      SEASONAL_ALLERGIES == "no")
#' 
#' # Sample sizes
#' phyloseq::nsamples(amgut_season_yes)
#' phyloseq::nsamples(amgut_season_no)
#' 
#' # Make sample sizes equal to ensure comparability
#' n_yes <- phyloseq::nsamples(amgut_season_yes)
#' 
#' amgut_season_no <- phyloseq::subset_samples(amgut_season_no, X.SampleID %in% 
#'                                      get_variable(amgut_season_no, 
#'                                      "X.SampleID")[1:n_yes])
#' 
#' # Network construction
#' amgut_net <- netConstruct(data = amgut_season_yes,
#'                           data2 = amgut_season_no,
#'                           measure = "pearson",
#'                           filtTax = "highestVar",
#'                           filtTaxPar = list(highestVar = 50),
#'                           zeroMethod = "pseudoZO", 
#'                           normMethod = "clr",
#'                           sparsMethod = "thresh",
#'                           thresh = 0.4,
#'                           seed = 123456)
#' 
#' # Estimated and sparsified associations of group 1
#' plotHeat(amgut_net$assoEst1, textUpp = "none", labCex = 0.6)
#' plotHeat(amgut_net$assoMat1, textUpp = "none", labCex = 0.6)
#' 
#' # Compute graphlet correlation matrices and perform significance tests
#' adja1 <- amgut_net$adjaMat1
#' adja2 <- amgut_net$adjaMat2
#' 
#' gcm1 <- calcGCM(adja1)
#' gcm2 <- calcGCM(adja2)
#' 
#' gcmtest <- testGCM(obj1 = gcm1, obj2 = gcm2)
#' 
#' # Mixed heatmap of GCM1 and significance codes
#' plotHeat(mat = gcmtest$gcm1, 
#'          pmat = gcmtest$pAdjust1,
#'          type = "mixed", 
#'          textLow = "code")
#' 
#' # Mixed heatmap of GCM2 and p-values (diagonal disabled)
#' plotHeat(mat = gcmtest$gcm1, 
#'          pmat = gcmtest$pAdjust1,
#'          diag = FALSE,
#'          type = "mixed", 
#'          textLow = "pmat")
#' 
#' # Mixed heatmap of differences (GCM1 - GCM2) and significance codes
#' plotHeat(mat = gcmtest$diff, 
#'          pmat = gcmtest$pAdjustDiff,
#'          type = "mixed", 
#'          textLow = "code",
#'          title = "Differences between GCMs (GCM1 - GCM2)",
#'          mar = c(0, 0, 2, 0))
#' 
#' # Heatmap of differences (insignificant values are blank)
#' plotHeat(mat = gcmtest$diff, 
#'          pmat = gcmtest$pAdjustDiff,
#'          type = "full", 
#'          textUpp = "sigmat")
#' 
#' # Same as before but with higher significance level
#' plotHeat(mat = gcmtest$diff, 
#'          pmat = gcmtest$pAdjustDiff,
#'          type = "full", 
#'          textUpp = "sigmat",
#'          argsUpp = list(sig.level = 0.1))
#' 
#' # Heatmap of absolute differences
#' # (different position of labels and legend)
#' plotHeat(mat = gcmtest$absDiff, 
#'          type = "full", 
#'          labPos = "d",
#'          legendPos = "b")
#' 
#' # Mixed heatmap of absolute differences 
#' # (different methods, text options, and color palette)
#' plotHeat(mat = gcmtest$absDiff,
#'          type = "mixed",
#'          textLow = "mat",
#'          methUpp = "number",
#'          methLow = "circle",
#'          labCol = "black",
#'          textCol = "gray50",
#'          textCex = 1.3,
#'          textFont = 2,
#'          digits = 1L,
#'          colorLim = range(gcmtest$absDiff),
#'          colorPal = "Blues",
#'          nCol = 21L,
#'          bg = "darkorange",
#'          addWhite = FALSE)
#' 
#' # Mixed heatmap of differences 
#' # (different methods, text options, and color palette)
#' plotHeat(mat = gcmtest$diff,
#'          type = "mixed",
#'          textLow = "none",
#'          methUpp = "number",
#'          methLow = "pie",
#'          textCex = 1.3,
#'          textFont = 2,
#'          digits = 1L,
#'          colorLim = range(gcmtest$diff),
#'          colorPal = "PiYG",
#'          nCol = 21L,
#'          bg = "gray80")
#' 
#' # Heatmap of differences with given color vector
#' plotHeat(mat = gcmtest$diff, 
#'          nCol = 21L,
#'          color = grDevices::colorRampPalette(c("blue", "white", "orange"))(31))
#'   
#' 
#' @importFrom corrplot corrplot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export

plotHeat <- function(mat, 
                     pmat = NULL,
                     type = "full",
                     textUpp = "mat",
                     textLow = "code",
                     methUpp = "color",
                     methLow = "color",
                     diag = TRUE,
                     title = "",
                     mar = c(0, 0, 1, 0),
                     labPos = "lt",
                     labCol = "gray40",
                     labCex = 1.1,
                     textCol = "black",
                     textCex = 1,
                     textFont = 1,
                     digits = 2L,
                     legendPos = "r",
                     colorPal = NULL,
                     addWhite = TRUE,
                     nCol = 51L,
                     colorLim = NULL,
                     revCol = FALSE,
                     color = NULL,
                     bg = "white",
                     argsUpp = NULL,
                     argsLow = NULL) {
  
  # Check input arguments
  argsIn <- as.list(environment())
  
  argsOut <- .checkArgsPlotHeat(argsIn)
  
  for (i in 1:length(argsOut)) {
    assign(names(argsOut)[i], argsOut[[i]])
  }
  
  #-----------------------------------------------------------------------------
  # Colors
  
  # Sequential and diverging palettes from RColorBrewer function
  sequential <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", 
                  "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
                  "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")
  
  diverging <- c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", 
                 "RdYlGn", "Spectral")
  
  # Set color limits
  if (is.null(colorLim)) {
    mrange <- range(mat)
    
    if (mrange[1] >= -1 & mrange[2] <= 1) {
      if (min(mat) < 0) {
        colorLim <- c(-1, 1)
        
      } else {
        colorLim <-  c(0, 1)
      }
    } else {
      colorLim <- mrange
    }
  }
  
  if (is.null(color)) {
    
    # Use diverging colors if mat in [-1, 1] and sequential colors otherwise
    if (is.null(colorPal)) {
      if (colorLim[1] < 0) {
        colorPal <- "RdBu"
        
      } else {
        colorPal <-  "YlOrRd"
      }
    }
    
    if (colorPal %in% diverging) {
      # Get colors of chosen palette
      color <- RColorBrewer::brewer.pal(n = 11, name = colorPal)
      
      # Replace middle color by white
      if (addWhite) color[6] <- "#FFFFFF"
      
      # Interpolate color vector to nCol colors
      color <- grDevices::colorRampPalette(color)(nCol)
      
    } else {
      
      # Get colors of chosen palette
      color <- RColorBrewer::brewer.pal(n = 9, name = colorPal)
      
      # Replace first color by white
      if (addWhite) color[1] <- "#FFFFFF"
      
      # Interpolate color vector to nCol colors
      color <- grDevices::colorRampPalette(color)(nCol)
    }
    
    # Revert color vector
    if (revCol) color <- rev(color)
  }
  
  #-----------------------------------------------------------------------------
  # Set arguments for the different text options
  argsValues <- list(corr = mat,
                     addCoef.col = textCol,
                     number.cex = textCex,
                     number.font = textFont,
                     insig = "n",
                     number.digits = digits)
  
  if (type == "mixed" && !diag) {
    diag(argsValues$corr) <- 0
    argsValues$p.mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    dimnames(argsValues$p.mat) <- dimnames(mat)
    diag(argsValues$p.mat) <- 1
    argsValues$sig.level <- 0.9
    argsValues$insig <- "blank"
  }
  
  argsSigVals <- list(corr = mat,
                      p.mat = pmat, 
                      insig ='blank',
                      addCoef.col = textCol, 
                      number.cex = textCex,
                      number.font = textFont,
                      number.digits = digits,
                      sig.level = 0.05)
  
  argsPvals <- list(corr = mat,
                    p.mat = pmat,
                    insig = 'p-value',
                    sig.level = -1,
                    pch.col = textCol,
                    pch.cex = textCex,
                    number.cex = textCex,
                    number.font = textFont,
                    number.digits = digits)
  
  argsCode <- list(corr = mat,
                   p.mat = pmat,
                   sig.level = c(0.001, 0.01, 0.05),
                   insig = "label_sig",
                   pch.cex = textCex,
                   pch.col = textCol)
  
  argsNone <- list(corr = mat)
  
  # Arguments added in all cases
  argsAdd <- list(order = "original",
                  tl.pos = labPos,
                  tl.col = labCol,
                  tl.cex = labCex,
                  cl.pos = legendPos,
                  col = color,
                  bg = bg,
                  is.corr = FALSE,
                  col.lim = colorLim,
                  title = title, 
                  mar = mar)
  
  # If TRUE, the upper triangle is plotted first
  upperFirst <- TRUE
  
  #-----------------------------------------------------------------------------
  # Arguments for the upper triangle
  
  argsUpper <- NULL
  
  if (type %in% c("mixed", "full", "upper")) {
    argsUpper <- switch (textUpp,
                         "mat" = argsValues, 
                         "pmat" = argsPvals, 
                         "code" = argsCode, 
                         "sigmat" = argsSigVals, 
                         "none" = argsNone)
    
    argsUpper$method <- methUpp
    
    argsUpper <- append(argsUpper, argsAdd)
    
    argsUpper$type <- switch(type,
                             "mixed" = "upper",
                             "full" = "full",
                             "upper" = "upper")
    
    if (type %in% c("full", "upper") && !diag) {
      argsUpper$diag <- FALSE
      
    } else if (type == "mixed") {
      if (textUpp == "mat") {
        argsUpper$add <- TRUE
        upperFirst <- FALSE
      }
    }
    
    # Append and replace arguments by user-defined arguments
    if (!is.null(argsUpp)) {
      argsUpper <- append(argsUpper, 
                          argsUpp[!names(argsUpp) %in% names(argsUpper)])
      argsUpper[names(argsUpp)] <- argsUpp
    }
    
  }
  
  #-----------------------------------------------------------------------------
  # Arguments for the lower triangle
  
  argsLower <- NULL
  
  if (type %in% c("mixed", "lower")) {
    argsLower <- switch (textLow,
                         "mat" = argsValues, 
                         "pmat" = argsPvals, 
                         "code" = argsCode, 
                         "sigmat" = argsSigVals, 
                         "none" = argsNone)
    
    argsLower$method <- methLow
    
    argsLower <- append(argsLower, argsAdd)
    
    argsLower$type <- "lower"
    
    if (type == "mixed") {
      if (textUpp != "mat") {
        argsLower$add <- TRUE
      }
      
    } else if (!diag) {
      argsLower$diag <- FALSE
    }
    
    # Append and replace arguments by user-defined arguments
    if (!is.null(argsLow)) {
      argsLower <- append(argsLower, 
                          argsLow[!names(argsLow) %in% names(argsLower)])
      argsLower[names(argsLow)] <- argsLow
    }
  }
  
  #-----------------------------------------------------------------------------
  # Plot the heatmaps
  if (type %in% c("mixed", "full", "upper")) {
    argsFirst <- if (upperFirst) argsUpper else argsLower
    
    .suppress_warnings(base::do.call(corrplot, args = argsFirst), 
                       startsWith, 
                       "col.lim interval ")
    
  }
  
  if (type %in% c("mixed", "lower")) {
    argsSec <- if (upperFirst) argsLower else argsUpper
    
    .suppress_warnings(base::do.call(corrplot, args = argsSec), 
                       startsWith, 
                       "col.lim interval ")
  }
  
  output <- list(argsUpper = argsUpper, argsLower = argsLower)
  
  invisible(output)
}





