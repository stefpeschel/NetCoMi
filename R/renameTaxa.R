#' @title Rename taxa
#'
#' @description Function for renaming taxa in a taxonomic table, which can be 
#'   given as matrix or phyloseq object. \cr\cr
#'   It comes with functionality for making unknown 
#'   and unclassified taxa unique and substituting them by the next higher known
#'   taxonomic level, e.g., an unknown genus "g__" can automatically be renamed 
#'   to "1_Streptococcaceae(F)".
#'   User-defined patterns determine the format of known and substituted names. 
#'   Unknown names (e.g., NAs) and unclassified taxa can be 
#'   handled separately. Duplicated names within one or more chosen ranks can 
#'   also be made unique by numbering them consecutively.
#'
#' @param taxtab taxonomic table (matrix containing the taxonomic names; columns 
#'   must be taxonomic ranks) or phyloseq object.
#' @param pat character specifying the pattern of new taxonomic names if the 
#'   current name is KNOWN.
#'   See the examples and default value for a demo. 
#'   Possible space holders are:
#'   \describe{
#'   \item{\code{<name>}}{Taxonomic name (either the original or replaced one)}
#'   \item{\code{<rank>}}{Taxonomic rank in lower case}
#'   \item{\code{<Rank>}}{Taxonomic rank with first letter in upper case}
#'   \item{\code{<r>}}{Abbreviated taxonomic rank in lower case}
#'   \item{\code{<R>}}{Abbreviated taxonomic rank in upper case}
#'   }
#' @param substPat character specifying the pattern of new taxonomic names if the 
#'   current name is UNKNOWN. The current name is substituted by the next higher 
#'   existing name. 
#'   Possible space holders (in addition to that of \code{pat}):
#'   \describe{
#'   \item{\code{<subst_name>}}{Substituted taxonomic name (next higher existing 
#'   name)}
#'   \item{\code{<subst_rank>}}{Taxonomic rank of substitute name in lower case}
#'   \item{\code{<subst_Rank>}}{Taxonomic rank of substitute name with first letter 
#'   in upper case}
#'   \item{\code{<subst_r>}}{Abbreviated taxonomic rank of substitute name in 
#'   lower case}
#'   \item{\code{<subst_R>}}{Abbreviated taxonomic rank of substitute name in 
#'   upper case}
#'   }
#' @param unknown character vector giving the labels of unknown taxa, without 
#'   leading rank label (e.g., "g_" or "g__" for genus level). If 
#'   \code{numUnknown = TRUE}, unknown names are replaced by a number.
#' @param numUnknown logical. If \code{TRUE}, a number is assigned to all 
#'   unknown taxonomic names (defined by \code{unknown}) to make them unique.
#' @param unclass character vector giving the label of unclassified taxa, 
#'   without leading rank label (e.g., "g_" or "g__" for genus level). If 
#'   \code{numUnclass = TRUE}, a number is added to the names of unclassified 
#'   taxa. Note that unclassified taxa and unknown taxa get a separate numbering 
#'   if \code{unclass} is set. To replace all unknown and unclassified taxa by 
#'   numbers, add "unclassified" (or the appropriate counterpart) to 
#'   \code{unknown} and set \code{unclass} to \code{NULL}.
#' @param numUnclass logical. If \code{TRUE}, a number is assigned to all 
#'   unclassified taxa (defined by \code{unclass}) to make them unique. 
#'   The pattern is defined via \code{numUnclassPat}.
#' @param numUnclassPat character defining the pattern used for numbering 
#'   unclassified taxa. Must include a space holder for the name ("<name>") and
#'   one for the number ("<num>"). Default is "<name><num>" resulting e.g., in
#'   "unclassified1".
#' @param numDupli character vector giving the ranks that should be made unique 
#'   by adding a number. Elements must match column names. The pattern is 
#'   defined via \code{numDupliPat}.
#' @param numDupliPat character defining the pattern used for numbering 
#'   duplicated names (if \code{numDupli} is given). Must include a space holder 
#'   for the name ("<name>") and one for the number ("<num>"). 
#'   Default is "<name><num>" resulting e.g., in "Ruminococcus1".
#' @param ranks character vector giving rank names used for renaming the 
#'   taxa. If \code{NULL}, the functions tries to automatically set rank names
#'   based on common usage.
#' @param ranksAbb character vector giving abbreviated rank names, which are 
#'   directly used for the place holders <r>, <subst_r>, <R>, and <subst_R> 
#'   (the former two in lower case and the latter two in upper case). 
#'   If \code{NULL}, the first letter of the rank names is used.
#' @param ignoreCols numeric vector with columns to be ignored. Names remain
#'   unchanged for these columns. Columns containing \code{NA}s are ignored 
#'   automatically only if \code{ignoreCols = NULL}. 
#'   Note: length of \code{ranks} and \code{ranksAbb} must 
#'   match the number of non-ignored columns.
#' @return Renamed taxonomic table (matrix or phyloseq object, depending on the 
#'   input).
#'
#' @examples
#' #--- Load and edit data -----------------------------------------------------
#' 
#' library(phyloseq)
#' data("GlobalPatterns")
#' global <- subset_taxa(GlobalPatterns, Kingdom == "Bacteria")
#' taxtab <- global@tax_table@.Data[1:10, ]
#' 
#' # Add some unclassified taxa
#' taxtab[c(2,3,5), "Species"] <- "unclassified"
#' taxtab[c(2,3), "Genus"] <- "unclassified"
#' taxtab[2, "Family"] <- "unclassified"
#' 
#' # Add some blanks
#' taxtab[7, "Genus"] <- " "
#' taxtab[7:9, "Species"] <- " "
#' 
#' # Add taxon that is unclassified up to Kingdom
#' taxtab[9, ] <- "unclassified"
#' taxtab[9, 1] <- "Unclassified"
#' 
#' # Add row names
#' rownames(taxtab) <- paste0("OTU", 1:nrow(taxtab))
#' 
#' print(taxtab)
#' 
#' #--- Example 1 (default setting) --------------------------------------------
#' 
#' # Example 1 (default setting)
#' # - Known names are replaced by "<r>_<name>"
#' # - Unknown names are replaced by "<r>_<name>_<subst_r>_<subst_name>"
#' # - Unclassified taxa have separate numbering
#' # - Ranks are taken from column names
#' # - e.g., unknown genus -> "g_1_f_Streptococcaceae"
#' 
#' renamed1 <- renameTaxa(taxtab)
#' renamed1
#' 
#' #--- Example 2 --------------------------------------------------------------
#' # - Use phyloseq object (subset of class clostridia to decrease runtime)
#' 
#' global_sub <- subset_taxa(global, Class == "Clostridia")
#' 
#' renamed2 <- renameTaxa(global_sub)
#' tax_table(renamed2)[1:5, ]
#' 
#' #--- Example 3 --------------------------------------------------------------
#' # - Known names remain unchanged
#' # - Substituted names are indicated by their rank in brackets
#' # - Pattern for numbering unclassified taxa changed
#' # - e.g., unknown genus -> "Streptococcaceae (F)"
#' # - Note: Numbering of unknowns is not shown because "<name>" is not 
#' #   included in "substPat"
#' 
#' renamed3 <- renameTaxa(taxtab, numUnclassPat = "<name>_<num>",
#'                          pat = "<name>", 
#'                          substPat = "<subst_name> (<subst_R>)")
#' renamed3
#' 
#' #--- Example 4 --------------------------------------------------------------
#' # - Same as before but numbering shown for unknown names
#' # - e.g., unknown genus -> "1 Streptococcaceae (F)"
#' 
#' renamed4 <- renameTaxa(taxtab, numUnclassPat = "<name>_<num>",
#'                          pat = "<name>", 
#'                          substPat = "<name> <subst_name> (<subst_R>)")
#' renamed4
#' 
#' #--- Example 5 --------------------------------------------------------------
#' # - Same numbering for unkown names and unclassified taxa
#' # - e.g., unknown genus -> "1_Streptococcaceae(F)"
#' # - Note: We get a warning here because "Unclassified" (with capital U) 
#' #   are not included in "unknown" but occur in the data
#' 
#' renamed5 <- renameTaxa(taxtab, unclass = NULL,
#'                          unknown = c(NA, " ", "unclassified"), 
#'                          pat = "<name>", 
#'                          substPat = "<name>_<subst_name>(<subst_R>)")
#' renamed5
#' 
#' #--- Example 6 --------------------------------------------------------------
#' # - Same as before, but OTU9 is now renamed correctly
#' 
#' renamed6 <- renameTaxa(taxtab, unclass = NULL,
#'                          unknown = c(NA, " ", "unclassified", "Unclassified"),
#'                          pat = "<name>", 
#'                          substPat = "<name>_<subst_name>(<subst_R>)")
#' renamed6
#' 
#' #--- Example 7 --------------------------------------------------------------
#' # - Add "(<Rank>: unknown)" to unknown names
#' # - e.g., unknown genus -> "1 Streptococcaceae (Genus: unknown)"
#' 
#' renamed7 <- renameTaxa(taxtab, unclass = NULL,
#'                          unknown = c(NA, " ", "unclassified", "Unclassified"),
#'                          pat = "<name>", 
#'                          substPat = "<name> <subst_name> (<Rank>: unknown)")
#' renamed7
#' 
#' #--- Example 8 --------------------------------------------------------------
#' # - Do not substitute unknowns and unclassified taxa by higher ranks
#' # - e.g., unknown genus -> "1"
#' 
#' renamed8 <- renameTaxa(taxtab, 
#'                          pat = "<name>", substPat = "<name>")
#' renamed8
#' 
#' #--- Example 9 --------------------------------------------------------------
#' # - Error if ranks cannot be automatically determined 
#' #   from column names or taxonomic names
#' 
#' taxtab_noranks <- taxtab
#' colnames(taxtab_noranks) <- paste0("Rank", 1:ncol(taxtab))
#' head(taxtab_noranks)
#' 
#' \dontrun{
#' renamed9 <- renameTaxa(taxtab_noranks, 
#'                          pat = "<name>", 
#'                          substPat = "<name>_<subst_name>(<subst_R>)")
#' }
#' 
#' # Ranks can either be given via "ranks" ... 
#' (ranks <- colnames(taxtab))
#' 
#' renamed9 <- renameTaxa(taxtab_noranks, 
#'                          pat = "<name>", 
#'                          substPat = "<name>_<subst_name>(<subst_R>)",
#'                          ranks = ranks)
#' renamed9
#' 
#' # ... or "ranksAbb" (we now use the lower case within "substPat")
#' (ranks <- substr(colnames(taxtab), 1, 1))
#' 
#' renamed9 <- renameTaxa(taxtab_noranks, 
#'                          pat = "<name>", 
#'                          substPat = "<name>_<subst_name>(<subst_r>)",
#'                          ranksAbb = ranks)
#' renamed9
#' 
#' #--- Example 10 -------------------------------------------------------------
#' # - Make names of ranks "Family" and "Order" unique by adding numbers to 
#' #   duplicated names
#' 
#' renamed10 <- renameTaxa(taxtab, 
#'                           pat = "<name>", 
#'                           substPat = "<name>_<subst_name>(<subst_R>)",
#'                           numDupli = c("Family", "Order"))
#' renamed10
#' 
#' any(duplicated(renamed10[, "Family"]))
#' any(duplicated(renamed10[, "Order"]))
#'
#' @importFrom phyloseq tax_table
#' @export

renameTaxa <- function(taxtab, 
                       pat = "<r>_<name>",
                       substPat = "<r>_<name>_<subst_r>_<subst_name>",
                       unknown = c(NA, "", " ", "__"),
                       numUnknown = TRUE,
                       unclass = c("unclassified", "Unclassified"),
                       numUnclass = TRUE,
                       numUnclassPat = "<name><num>",
                       numDupli = NULL,
                       numDupliPat = "<name><num>",
                       ranks = NULL, ranksAbb = NULL, 
                       ignoreCols = NULL) {
  
  stopifnot(is.character(pat))
  stopifnot(is.character(substPat))
  stopifnot(is.logical(numUnknown))
  stopifnot(is.logical(numUnclass))
  
  if (inherits(taxtab, "phyloseq")) {
    tax <- as.matrix(taxtab@tax_table@.Data)
  } else {
    tax <- as.matrix(taxtab)
  }
  
  # Input check for numUnclass and numDupli
  
  if (numUnclass) {
    if (!grepl("<name>", numUnclassPat)) {
      stop('Argument "numUnclassPat" must contain a space holder "<name>".')
    }
    if (!grepl("<num>", numUnclassPat)) {
      stop('Argument "numUnclassPat" must contain a space holder "<num>".')
    }
  }
  
  if (!is.null(numDupli)) {
    if (!all(numDupli %in% colnames(tax))) {
      stop('Ranks given with "numDupli" must match column names of ' ,
           'taxonomic table.')
    }
    
    if (!grepl("<name>", numDupliPat)) {
      stop('Argument "numDupliPat" must contain a space holder "<name>".')
    }
    
    if (!grepl("<num>", numDupliPat)) {
      stop('Argument "numDupliPat" must contain a space holder "<num>".')
    }
  }
  
  # Remove columns to ignore
  if (!is.null(ignoreCols)) {
    if (is.character(ignoreCols)) {
      ignoreCols <- which(colnames(tax) %in% ignoreCols)
    }
    
  } else {
    # Ignore columns that contain NAs only
    ignore.tmp <- which(apply(tax, 2, function(x) all(is.na(x))))
    
    if (length(ignore.tmp) == 1) {
      message(paste0("Column ", ignore.tmp, 
                     " contains NAs only and is ignored."))
      ignoreCols <- ignore.tmp
      
    } else if (length(ignore.tmp) > 1) {
      message(paste0("Columns ", ignore.tmp, 
                     " contain NAs only and are ignored."))
      ignoreCols <- ignore.tmp
    }
  }
  
  tax_orig <- tax
  
  if (!is.null(ignoreCols)) {
    tax <- tax_orig[, -ignoreCols]
    
  } 
  
  # Define ranks
  nranks <- ncol(tax)
  
  ranks.tmp <- c("Life", "Domain", "Kingdom", "Phylum", "Class", "Order", 
                 "Family", "Genus", "Species")
  ranksAbb.tmp <- substr(tolower(ranks.tmp), 1, 1)
  
  missRanks <- missRanksAbb <- FALSE
  
  if (is.null(ranks)) {
    if (any(ranks.tmp %in% colnames(tax))) {
      ranks <- colnames(tax)
      
    } else {
      
      repl <- 0
      ranks <- rep(NA, ncol(tax))
      
      for (r in seq_along(ranks.tmp)) {
        sel <- grep(paste0(ranksAbb.tmp, "_")[r], t(tax)) %% nranks
        
        if (length(sel) != 0) {
          sel <- sel[1]
          if (sel == 0) sel <- nranks
          ranks[sel] <- ranks.tmp[r]
          repl <- repl + 1
        }
      }
      
      if (repl != nranks) {
        missRanks <- TRUE
        
      } 
    }
    
  } else {
    if (length(ranks) != nranks) {
      stop("Length of 'ranks' does not match number of columns")
    }
  }
  
  
  if (is.null(ranksAbb)) {
    if (missRanks) {
      missRanksAbb <- TRUE
    } else {
      ranksAbb <- substr(tolower(ranks), 1, 1)
    }
    
  } else {
    if (length(ranksAbb) != nranks) {
      stop("Length of 'ranksAbb' does not match number of columns")
    }
  }
  
  # Stop if full rank names are not given but needed
  if (missRanks && (grepl("<rank>", pat) | 
                    grepl("<Rank>", pat))) {
    stop("Ranks could not be determined but are needed for the given pattern. ", 
         "Please provide 'ranks'")
  }
  
  if (missRanks && (grepl("<rank>", substPat) | 
                    grepl("<Rank>", substPat) | 
                    grepl("<subst_rank>", substPat) | 
                    grepl("<subst_Rank>", substPat))) {
    stop("Ranks could not be determined but are needed for the given ", 
         "substitute pattern. Please provide 'ranks'")
  }
  
  # Stop if abbreviated rank names are not given but needed
  if (missRanksAbb && (grepl("<r>", pat) | 
                       grepl("<R>", pat))) {
    stop("Abbreviated ranks could not be determined but are needed for the ", 
         "given pattern. Please provide 'ranks' or 'ranksAbb'.")
  }
  
  if (missRanksAbb && (grepl("<r>", substPat) | 
                       grepl("<R>", substPat) | 
                       grepl("<subst_r>", substPat) | 
                       grepl("<subst_R>", substPat))) {
    stop("Ranks could not be determined but are needed for the given ", 
         "substitute pattern. Please provide 'ranks' or 'ranksAbb'.")
  }
  
  
  for (r in seq_along(ranksAbb.tmp)) {
    tax <- gsub(paste0(ranksAbb.tmp[r], "__"), "", tax)
    tax <- gsub(paste0(ranksAbb.tmp[r], "_"), "", tax)
  }
  
  #-----------------------------------------------------------------------------
  ### Make labels unique
  
  unknownMat <- matrix(data = FALSE, nrow = nrow(tax), ncol  = nranks)
  dimnames(unknownMat) <- dimnames(tax)
  
  # Add number to unclassified taxa
  if (!is.null(unclass)) {
    
    if (numUnclass) {
      # Split pattern
      pat.tmp <- unlist(strsplit(numUnclassPat, "<"))
      pat.tmp <- unlist(strsplit(pat.tmp, ">"))
      numUnclassPat <- pat.tmp
      
      for (c in 1:nranks) {
        # Detect unclassified taxa
        isunclass <- tax[, c] %in% unclass
        
        if (sum(isunclass) > 0) {
          # Replace names according to pattern
          tmp <- matrix(rep(numUnclassPat, 
                            each = sum(isunclass)), nrow = sum(isunclass))
          tmp[tmp == "name"] <- tax[isunclass, c]
          tmp[tmp == "num"] <- 1:sum(isunclass)
          
          tax[isunclass, c] <- apply(tmp, 1, paste, collapse="")
          
          unknownMat[isunclass, c] <- TRUE
        }
      }
    }
    
  } else {
    if (any(grepl("unclass", tax, ignore.case = TRUE))) {
      unclassname <- unique(tax[grep("unclass", tax, ignore.case = TRUE)])
      
      if (!all(unclassname %in% unknown)) {
        unclassname <- unclassname[!unclassname %in% unknown]
        warning(paste0('Taxonomic table contains unclassified taxa. ',
                       'Consider adding "', unclassname, 
                       '" to argument "unknown".'))
      }
    }
  }
  
  # Replace unknowns by number
  if (numUnknown) {
    for (c in 1:nranks) {
      
      # Add a number if the rank is unknown
      isunknown <- tax[, c] %in% unknown
      
      if (sum(isunknown) > 0) {
        tax[isunknown, c] <- 1:sum(isunknown)
        unknownMat[isunknown, c] <- TRUE
      }
    }
  }
  
  # Make duplicates unique
  if (!is.null(numDupli)) {
    # Split pattern
    pat.tmp <- unlist(strsplit(numDupliPat, "<"))
    pat.tmp <- unlist(strsplit(pat.tmp, ">"))
    numDupliPat <- pat.tmp
    
    for (c in seq_along(numDupli)) {
      tax.tmp <- tax[, numDupli[c]]
      dupli <- table(tax.tmp)[table(tax.tmp) > 1]
      
      if (length(dupli) > 0) {
        for (d in seq_along(dupli)) {
          # Replace names according to pattern
          tmp <- matrix(rep(numDupliPat, each = dupli[d]), nrow = dupli[d])
          tmp[tmp == "name"] <- names(dupli[d])
          tmp[tmp == "num"] <- 1:dupli[d]
          
          tax.tmp[tax.tmp == names(dupli[d])] <- apply(tmp, 1, paste, 
                                                       collapse="")
        }
        
        tax[, numDupli[c]] <- tax.tmp
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  
  pat.tmp <- unlist(strsplit(pat, "<"))
  pat.tmp <- unlist(strsplit(pat.tmp, ">"))
  pat <- pat.tmp
  
  pat.tmp <- unlist(strsplit(substPat, "<"))
  pat.tmp <- unlist(strsplit(pat.tmp, ">"))
  substPat <- pat.tmp
  
  for (i in 1:nrow(tax)) {
    
    for (r in nranks:1) {
      name <- tax[i, r]
      isunknown <- unknownMat[i, r]
      usesubst <- NA
      
      if (isunknown) {
        usesubst <- TRUE
        
        # Which ranks are known?
        knownr <- which(!unknownMat[i, ])
        
        if (length(knownr) == 0) {
          # If all ranks are unknown, take name of the highest rank
          sr <- 1
          if (r == 1) usesubst <- FALSE
          
        } else {
          # Select last known rank
          sr <- knownr[length(knownr)]
        }
        
      } else {
        usesubst <- FALSE
      }
      
      # Assign the right pattern
      if (usesubst) {
        patout <- substPat
      } else {
        patout <- pat
      }
      
      
      patout[patout == "r"] <- tolower(ranksAbb[r])
      patout[patout == "R"] <- toupper(ranksAbb[r])
      patout[patout == "rank"] <- tolower(ranks[r])
      patout[patout == "Rank"] <- paste0(toupper(substr(ranks[r], 1, 1)), 
                                         substr(ranks[r], 2, 
                                                nchar(ranks[r])))
      patout[patout == "name"] <- name
      
      if (usesubst) {
        patout[patout == "subst_r"] <- tolower(ranksAbb[sr])
        patout[patout == "subst_R"] <- toupper(ranksAbb[sr])
        patout[patout == "subst_rank"] <- tolower(ranks[sr])
        patout[patout == "subst_Rank"] <- paste0(toupper(substr(ranks[sr], 1, 1)), 
                                                 substr(ranks[sr], 2, 
                                                        nchar(ranks[sr])))
        patout[patout == "subst_name"] <- tax[i, sr]
      }
      
      taxname <- paste(patout, collapse = "")
      
      tax[i, r] <- taxname
    }
  }
  
  if (!is.null(ignoreCols)) {
    taxout <- tax_orig
    taxout[, (1:ncol(tax_orig))[-ignoreCols]] <- tax
  } else {
    taxout <- tax
  }
  
  if (inherits(taxtab, "phyloseq")) {
    phyloseq::tax_table(taxtab) <- phyloseq::tax_table(taxout)
    taxout <- taxtab
  } 
  
  return(taxout)
}




