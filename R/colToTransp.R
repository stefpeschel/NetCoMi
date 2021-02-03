#' Adding transparency to a color
#'
#' @param col color vector specified similar to the \code{col} argument in
#'   \code{\link[grDevices]{col2rgb}}
#' @param percent numeric between 0 and 100 giving the level of transparency. 
#'   Defaults to 50.
#'   
#' @examples 
#' # Excepts hexadecimal strings, written colors, or numbers as input
#' colToTransp("#FF0000FF", 50)
#' colToTransp("black", 50)
#' colToTransp(2)
#' 
#' # Different shades of red
#' r80 <- colToTransp("red", 80)
#' r50 <- colToTransp("red", 50)
#' r20 <- colToTransp("red", 20)
#' 
#' barplot(rep(5, 4), col=c("red", r20, r50, r80), names.arg = 1:4)
#' 
#' # Vector as input
#' rain_transp <- colToTransp(rainbow(5), 50)
#' 
#' barplot(rep(5, 5), col = rain_transp, names.arg = 1:5)
#'
#' @importFrom grDevices col2rgb rgb
#' @export

colToTransp <- function(col, percent = 50){
  stopifnot(is.numeric(percent) && percent <= 100 && percent >= 0)
  
  rgbVal <- grDevices::col2rgb(col)
  
  transpcol <- grDevices::rgb(red = rgbVal[1,], green = rgbVal[2,], 
                              blue = rgbVal[3,], alpha = (100-percent)*255/100, 
                              maxColorValue = 255)
  
  names(transpcol) <- colnames(rgbVal)
  
  return(transpcol)
}

