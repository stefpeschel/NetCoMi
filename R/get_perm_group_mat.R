get_perm_group_mat <- function(n1, n2, n, nPerm, matchDesign) {
  if (is.null(matchDesign)) {
    
    # group vector: n1 samples in group1 and n2 samples in group 2
    groupvec <- c(rep(1, n1), rep(2, n2))
    
    # maximum number of permutations
    maxcomb <- choose(n, n1)
    
    if (maxcomb < nPerm) {
      tmp <- ifelse(maxcomb < 500, " (not recommended).", ".")
      stop(paste0("Possible number of permutations is smaller than 'nPerm'. ",
                  "Either increase sample size or set 'nPerm' to ", 
                  maxcomb, tmp))
    }
    
    if (maxcomb == nPerm) {
      # inverted group vector
      groupvec_inv <- c(rep(2, n1), rep(1, n2))
      
      # matrix with all possible permutations / combinations
      # (each column gives indices of samples in group 1)
      combmat <- utils::combn(n, n1)
      
      # matrix with permuted group labels (each column belongs to a sample)
      perm_group_mat <- matrix(NA, nrow = maxcomb, ncol = n)
      
      for (i in 1:maxcomb) {
        perm_group_mat[i, combmat[,i]] <- 1
        perm_group_mat[i, -combmat[,i]] <- 2
      }
      
    } else {
      perm_group_mat <- t(sapply(1:(nPerm + 10), function(x) sample(groupvec)))
      
      perm_group_mat <- unique.matrix(perm_group_mat)
      
      iter <- 0
      
      while (nrow(perm_group_mat) < nPerm) {
        perm_group_mat <- rbind(perm_group_mat, sample(groupvec))
        perm_group_mat <- unique.matrix(perm_group_mat)
        iter <- iter + 1
        if (iter == 100000) {
          stop(paste0("Unique matrix with permuted group labels ",
                      "could not be computed within 100000 iterations."))
        }
      }
    }
    
  } else {
    # size of matching sets
    set_size <- sum(matchDesign)
    
    # number of sets
    n_sets <- n / set_size
    
    # possible combinations within each set
    possib <- factorial(matchDesign[1] + matchDesign[2]) / 
      (factorial(matchDesign[1]) * factorial(matchDesign[2]))
    
    # maximum number of permutations
    maxcomb <- possib^n_sets
    
    if (maxcomb < nPerm) {
      tmp <- ifelse(maxcomb < 500, " (not recommended).", ".")
      stop(paste0("Possible number of permutations for this match design ", 
                  "is smaller than 'nPerm'. ",
                  "Either increase sample size or set 'nPerm' to ", 
                  maxcomb, tmp))
      
    } else {
      # group vector corresponding to matching design
      groups_orig <- rep(c(rep(1, matchDesign[1]), 
                           rep(2, matchDesign[2])), 
                         n_sets)
      
      perm_func <- function(x) {
        groups_perm <- NULL
        
        for (i in 1:n_sets) {
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
      
      iter <- 0
      
      while (nrow(perm_group_mat) < nPerm) {
        groups_perm <- NULL
        
        for (i in 1:n_sets) {
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
        
        iter <- iter + 1
        if (iter == 100000) {
          stop(paste0("Unique matrix with permuted group labels ",
                      "could not be computed within 100000 iterations."))
        }
      }
    }
  }
  
  perm_group_mat <- perm_group_mat[1:nPerm, ]
  
  return(perm_group_mat)
}

