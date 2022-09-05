#' @title Edit labels
#'
#' @description Function for editing node labels, i.e., shortening to a certain
#'   length and removing unwanted characters.\cr\cr
#'   The function is used by NetCoMi's plot functions 
#'   \code{\link{plot.microNetProps}} and \code{\link{plot.diffnet}}.
#'   
#' @details Consider a vector with three bacteria names: "Streptococcus1", 
#'   "Streptococcus2", and "Streptomyces".\cr\cr
#'   
#'   \code{shortenLabels = "simple"} with \code{labelLength = 6}
#'    leads to shortened labels: 
#'   "Strept", "Strept", and "Strept", which are not distinguishable. \cr\cr
#'   
#'   \code{shortenLabels = "intelligent"} with 
#'   \code{labelPattern = c(5, "'", 3)} leads to shortened labels:
#'   "Strep'coc", "Strep'coc", are "Strep'myc", where the first two are not 
#'   distinguishable.
#'   
#'   \code{shortenLabels = "intelligent"} with 
#'   \code{labelPattern = c(5, "'", 3, "'", 3)} leads to shortened labels:
#'   "Strep'coc'1", "Strep'coc'2", and "Strep'myc", from which the original
#'   labels can be inferred.\cr\cr
#'   
#'   The intelligent approach is as follows:\cr\cr 
#'   First, labels are shortened to the 
#'   defined length (argument \code{labelLength}). The \code{labelPattern} is 
#'   then applied to all duplicated labels. For each group of duplicates, 
#'   the third label part starts at the letter where two or more labels are 
#'   different for the first time. 
#'   The five-part pattern (if given) applies if a group of 
#'   duplicates consists of more than two labels and if the shortened labels are 
#'   not unique after applying the three-part pattern. Then, the fifth part 
#'   starts at the letter where all labels are different for the first time.
#'   \cr\cr  A message is printed if the returned labels are not unique.
#'
#' @param x character vector with node labels.
#' @param shortenLabels character indicating how to shorten the labels. 
#'   Available options are:
#'   \describe{
#'   \item{\code{"intelligent"}}{Elements of \code{charToRm} are removed,
#'   labels are shortened to length \code{labelLength}, and duplicates are
#'   removed using \code{labelPattern}.}
#'   \item{\code{"simple"}}{Elements of \code{charToRm} are  removed and labels
#'   are shortened to length \code{labelLength}.}
#'   \item{\code{"none"}}{Labels are not shortened.} }
#' @param labelLength integer defining the length to which labels shall
#'   be shortened if \code{shortenLabels} is used. Defaults to 6.
#' @param labelPattern vector of three or five elements, which is used if 
#'   argument \code{shortenLabels} is set to \code{"intelligent"}. 
#'   If cutting a label to length \code{labelLength} leads to duplicates, 
#'   the label is shortened according to \code{labelPattern}, 
#'   where the first entry gives the length of the first part, 
#'   the second entry is used a separator, and the third entry
#'   is the length of the third part. If \code{labelPattern} has five elements 
#'   and the shortened labels are still not unique, 
#'   the fourth element serves as further separator, and the fifth element gives
#'   the length of the last label part. Defaults to c(4, "'", 3, "'", 3). 
#'   See details for an example.
#' @param addBrack logical indicating whether to add a closing square bracket. 
#'   If \code{TRUE}, a "]" is added if the first part contains a "[".
#' @param charToRm character vector giving one or more patterns to remove from 
#'   the labels.
#' @param verbose logical. If \code{TRUE}, the function is allowed to return 
#'   messages.
#'
#' @return Character vector with edited labels.
#'
#' @examples
#' labels <- c("Salmonella", 
#'             "Clostridium", "Clostridiales(O)", 
#'             "Ruminococcus", "Ruminococcaceae(F)", 
#'             "Enterobacteriaceae", "Enterococcaceae",
#'             "[Bacillus] alkalinitrilicus",
#'             "[Bacillus] alkalisediminis",
#'             "[Bacillus] oceani")
#' 
#' # Use the "simple" method to shorten labels
#' editLabels(labels, shortenLabels = "simple", labelLength = 6)
#' # -> Original labels cannot be inferred from shortened labels
#' 
#' # Use the "intelligent" method to shorten labels with three-part pattern
#' editLabels(labels, shortenLabels = "intelligent", labelLength = 6,
#'            labelPattern = c(6, "'", 4))
#' # -> [Bacillus] alkalinitrilicus and [Bacillus] alkalisediminis not 
#' #    distinguishable
#' 
#' # Use the "intelligent" method to shorten labels with five-part pattern
#' editLabels(labels, shortenLabels = "intelligent", labelLength = 6,
#'            labelPattern = c(6, "'", 3, "'", 3))
#' 
#' # Same as before but no brackets are added
#' editLabels(labels, shortenLabels = "intelligent", labelLength = 6, 
#'            addBrack = FALSE, labelPattern = c(6, "'", 3, "'", 3))
#' 
#' # Remove character pattern(s) (can also be a vector with multiple patterns)
#' labels <- c("g__Faecalibacterium", "g__Clostridium", "g__Eubacterium", 
#'             "g__Bifidobacterium", "g__Bacteroides")
#'             
#' editLabels(labels, charToRm = "g__")
#' 
#' @export

editLabels <- function(x, 
                       shortenLabels = c("intelligent", "simple", "none"),
                       labelLength = 6, 
                       labelPattern = NULL,
                       addBrack = TRUE,
                       charToRm = NULL,
                       verbose = TRUE) {
  labels <- x
  
  shortenLabels <- match.arg(shortenLabels)
  
  stopifnot(is.numeric(labelLength) & labelLength >= 0)
  
  stopifnot(is.logical(addBrack))
  
  stopifnot(is.logical(verbose))
  
  # Define label pattern
  if (is.null(labelPattern)) {
    labelPattern <- c(4, "'", 3, "'", 3)
  } else {
    stopifnot(length(labelPattern) %in% c(3, 5))
  }
  
  # Intelligent approach
  if (shortenLabels == "intelligent") {
    
    for (char in charToRm) {
      labels <- gsub(char, "", labels)
    }
    
    labels[labels == ""] <- "-"
    shortlabels <- substring(labels, 1,labelLength)
    
    dupli <- which(duplicated(shortlabels))
    
    # find duplicates in variable names
    while(length(dupli) > 0) {
      lpat <- length(labelPattern)
      
      ind <- which(shortlabels == shortlabels[dupli[1]])
      dupnames <- strsplit(labels[ind], "")
      
      # Turn two consecutive numbers into double-digit numbers
      for (i in seq_along(dupnames)) {
        dupnames[[i]] <- .sing2doubDigit(dupnames[[i]])
      }
      
      # Make length of duplicate names equal
      lvec <- unlist(lapply(dupnames, length))
      l <- max(lvec)
      if (min(lvec) != max(lvec)) {
        for (i in 1:length(dupnames)) {
          if (lvec[i] < l) {
            dupnames[[i]] <- c(dupnames[[i]], rep(" ", l-lvec[i]))
          }
        }
      }
      
      pos <- .firstUnequalElement(dupnames)
      first_unequal <- pos$first_unequal
      all_unequal <- pos$all_unequal
      if (first_unequal == all_unequal) lpat <- 3
      
      cut1 <- as.numeric(labelPattern[1])
      cut2 <- as.numeric(labelPattern[3])
      cut2 <- min(cut2, l - first_unequal + 1)
      
      if (lpat == 5) {
        cut3 <- as.numeric(labelPattern[5])
        cut3 <- min(cut3, l - all_unequal + 1)
      } else {
        cut3 <- 0
      }
      
      
      if (first_unequal > length(dupnames[[1]])) {
        shortlabels[ind] <- substring(labels[ind], 1, cut1)
        
      } else {
        for (k in 1:length(dupnames)) {
          str1 <- paste(dupnames[[k]][1:cut1], collapse = "")
          str2 <- labelPattern[2]
          str3 <- paste(dupnames[[k]][first_unequal:(first_unequal+cut2-1)], 
                        collapse = "")
          
          if (lpat == 5) {
            str4 <- labelPattern[4]
            str5 <- paste(dupnames[[k]][all_unequal:(all_unequal+cut3-1)], 
                          collapse = "")
          }
          
          if (addBrack) {
            # Add "]" if str1 contains "["
            if (grepl("\\[", str1)) {
              if (grep("\\]", dupnames[[k]]) < first_unequal) {
                str1 <- paste0(str1, "]")
              }
            }
          }
          
          # Ignore str3 if empty
          str3spl <- strsplit(str3, "")[[1]]
          
          if (all(grepl(" ", str3spl))) {
            shortlabels[ind[k]] <- str1
            
          } else {
            if (lpat == 5) {
              
              # Ignore str5 if empty
              str5spl <- strsplit(str5, "")[[1]]
              
              if (!all(grepl(" ", str5spl))) {
                shortlabels[ind[k]] <- paste0(str1, str2, str3, str4, str5)
              } else {
                shortlabels[ind[k]] <- paste0(str1, str2, str3)
              }
              
            } else {
              shortlabels[ind[k]] <- paste0(str1, str2, str3)
            }
          }
          
        }
      }
      
      dupli <- dupli[!dupli %in% ind]
    }
    
    if (any(duplicated(shortlabels)) & verbose) {
      message("Shortened labels could not be made unique.")
    }
    
    labels <- shortlabels
    
    # Simple approach
  } else if (shortenLabels == "simple") {
    
    for (char in charToRm) {
      labels <- gsub(char, "", labels)
    }
    
    labels[labels == ""] <- "-"
    labels <- substr(labels, 1, labelLength)
    
    # Only unwanted characters removed, no shortening
  } else if (!is.null(charToRm)) {
    
    for (char in charToRm) {
      labels <- gsub(char, "", labels)
    }
    
  }
  
  return(labels)
}


