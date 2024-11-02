#' @keywords internal
.calcJaccard <- function(group1, group2, sigTest, greater) {
  N <- length(union(group1, group2))
  C <- length(intersect(group1, group2))
  jacc <- ifelse(N>0, C/N, 0)
  p.greater <- NA
  p.less <- NA
  
  if (sigTest) {
    
    if (C==0) {
      s = 0
    } else {
      s <- numeric(C)
      for (x in 0:(C-1)) {
        s[x+1] <- choose(N,x) * 2^(N-x)
      }
    }
    
    # probability for jaccard index being greater than or equal to given value
    p.greater <- 1-(sum(s) / 3^N)
    
    s <- numeric(C+1)
    for (x in 0:(C)) {
      s[x+1] <- choose(N,x) * 2^(N-x)
    }
    # probability for jaccard index being less than or equal to given value
    p.less <- sum(s) / 3^N
    
  }
  
  # p is the propability that the jaccard value differs from what would be
  # expected at random
  output <- c(jacc = jacc, p.greater = p.greater, p.less = p.less)
  
  return(output)
}
