get_perm_group_mat <- function(n1, n2, n, nPerm, matchDesign){
  
  if(is.null(matchDesign)){
    
    groupvec <- c(rep(1, n1), rep(2, n2))
    
    perm_group_mat <- t(sapply(1:(nPerm + 10), function(x) sample(groupvec)))
    
    perm_group_mat <- unique.matrix(perm_group_mat)
    
    while (nrow(perm_group_mat) < nPerm) {
      perm_group_mat <- rbind(perm_group_mat, sample(groupvec))
      perm_group_mat <- unique.matrix(perm_group_mat)
    }
    
  } else{
    # size of matching sets
    set_size <- sum(matchDesign)
    # number of sets
    n_sets <- n / set_size
    # group vector corresponding to matching design
    groups_orig <- rep(c(rep(1, matchDesign[1]), 
                         rep(2, matchDesign[2])), 
                       n_sets)
    
    perm_func <- function(x){
      groups_perm <- NULL
      
      for(i in 1:n_sets){
        groups_perm <- c(groups_perm, 
                         sample(c(rep(1, matchDesign[1]),
                                  rep(2, matchDesign[2]))))
      }
      
      # shuffled groups belonging to original group 1 
      groups_perm1 <- groups_perm[groups_orig == 1]
      
      # shuffled groups belonging to original group 2
      groups_perm2 <- groups_perm[groups_orig == 2]
      
      c(groups_perm1, groups_perm2)
    }
    
    perm_group_mat <- t(sapply(1:(nPerm + 10), perm_func))
    
    perm_group_mat <- unique.matrix(perm_group_mat)
    
    while (nrow(perm_group_mat) < nPerm) {
      groups_perm <- NULL
      
      for(i in 1:n_sets){
        groups_perm <- c(groups_perm, 
                         sample(c(rep(1, matchDesign[1]),
                                  rep(2, matchDesign[2]))))
      }
      
      # shuffled groups belonging to original group 1 
      groups_perm1 <- groups_perm[groups_orig == 1]
      
      # shuffled groups belonging to original group 2
      groups_perm2 <- groups_perm[groups_orig == 2]
      
      perm_group_mat <- rbind(perm_group_mat, c(groups_perm1, groups_perm2))
      perm_group_mat <- unique.matrix(perm_group_mat)
    }
    
  }
  
  perm_group_mat <- perm_group_mat[1:nPerm, ]
  
  return(perm_group_mat)
}

