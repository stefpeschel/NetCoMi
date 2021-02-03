#' @title Generate a Spring Embedded Layout
#'
#' @description The Fruchterman Reingold layout is used to determine node
#'   positions for a given adjacency matrix. R code from
#'   \code{\link[qgraph]{qgraph}} is used in this function.
#'
#' @param adjaMat adjacency matrix
#' @param repulse.rad positive numeric value indicating the strength of
#'   repulsive forces in the "spring" layout. Nodes are placed closer together
#'   for smaller values and further apart for higher values.
#'
#' @importFrom qgraph qgraph.layout.fruchtermanreingold

springLayout <- function(adjaMat, repulse.rad = NULL){
  E <- list()
  input <- adjaMat
  nNodes <- nrow(input)
  incl <- upper.tri(input, diag = TRUE)
  directed <- matrix(FALSE, nNodes, nNodes)
  E$from <- numeric(0)
  E$to <- numeric(0)
  E$weight <- numeric(0)
  E$from <- rep(1:nrow(input), times = nrow(input))
  E$to <- rep(1:nrow(input), each = nrow(input))
  E$weight <- c(input)
  E$from <- E$from[c(incl)]
  E$to <- E$to[c(incl)]
  E$weight <- E$weight[c(incl)]
  keep <- abs(E$weight) > 0
  E$from <- E$from[keep]
  E$to <- E$to[keep]
  E$weight <- E$weight[keep]

  layout <- qgraph.layout.fruchtermanreingold(cbind(E$from, E$to),
                                             abs(E$weight/max(abs(E$weight)))^2,
                                             nNodes, rotation = NULL,
                                             layout.control = 0.5, niter = NULL,
                                             max.delta = NULL,
                                             area = NULL, cool.exp = NULL,
                                             repulse.rad = repulse.rad,
                                             init = NULL, constraints = NULL)
  return(layout)
}



# taken from qgraph
num2hex <- function(x) {
  hex <- unlist(strsplit("0123456789ABCDEF", split = ""))
  
  return(paste(hex[(x - x%%16)/16 + 1], hex[x%%16 + 1], sep = ""))
}




first_unequal_element <- function(x,y){
  stopiter <- FALSE
  lx <- length(x)
  ly <- length(y)
  i <- 0
  
  while(stopiter == FALSE){
    i <- i + 1
    if(x[i] != y[i]) stopiter <- TRUE
    if(x[i] == "[" & y[i] == "[") stopiter <- TRUE
  }
  
  return(i)
}



