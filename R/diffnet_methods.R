#' @title Summary Method for Objects of Class diffnet
#'
#' @description The main results returned by \code{\link{diffnet}} are printed
#'   in a compact summary.
#'
#' @param object object of class \code{diffnet} returned by \code{\link{diffnet}}.
#' @param digits integer giving the number of decimal places to which numeric
#'   summaries are rounded. Defaults to 3L.
#' @param digitsPval integer giving the number of decimal places to which the
#'   p-values are rounded. Values at or below the displayed precision are shown
#'   in scientific notation. Defaults to 6L.
#' @param ... not used.
#'
#' @seealso \code{\link{diffnet}}
#'
#' @method summary diffnet
#' @rdname summarize.diffnet
#' @export
summary.diffnet <- function(object,
                            digits = 3L,
                            digitsPval = 6L,
                            ...) {

  if (!inherits(object, "diffnet")) {
    stop('"object" must be of class "diffnet".')
  }

  stopifnot(is.numeric(digits))
  digits <- as.integer(digits)

  stopifnot(is.numeric(digitsPval))
  digitsPval <- as.integer(digitsPval)

  asso1 <- object$assoMat1
  asso2 <- object$assoMat2

  alpha <- if (!is.null(object$alpha)) object$alpha else NA_real_
  lfdrThresh <- if (!is.null(object$lfdrThresh)) object$lfdrThresh else NA_real_
  adjust <- if (!is.null(object$adjusted)) object$adjusted else "none"

  nTestsTotal <- if (!is.null(object$nTestsTotal)) object$nTestsTotal else
    if (!is.null(object$pvalsVec)) length(object$pvalsVec) else NA_integer_
  nTestsAdjusted <- if (!is.null(object$nTestsAdjusted)) object$nTestsAdjusted else
    nTestsTotal

  meta <- data.frame(
    value = c(object$diffMethod,
              adjust,
              as.character(isTRUE(object$adjustNonzeroOnly)),
              nTestsTotal,
              nTestsAdjusted,
              nTestsTotal - nTestsAdjusted,
              alpha,
              lfdrThresh),
    row.names = c("Method",
                  "Adjustment",
                  "Adjust non-zero only",
                  "Pairwise tests total",
                  "Pairwise tests adjusted",
                  "Pairwise tests excluded",
                  "Alpha",
                  "lfdr threshold"),
    check.names = FALSE
  )
  names(meta) <- ""

  assoSummary <- rbind(
    "Network 1" = .summaryDiffnetNumeric(asso1[lower.tri(asso1)], digits),
    "Network 2" = .summaryDiffnetNumeric(asso2[lower.tri(asso2)], digits)
  )

  counts <- data.frame(
    "Network 1" = c(.countDiffnetEdges(asso1),
                    .countDiffnetSingletons(asso1)),
    "Network 2" = c(.countDiffnetEdges(asso2),
                    .countDiffnetSingletons(asso2)),
    row.names = c("Non-zero edges", "Singletons"),
    check.names = FALSE
  )

  diffVec <- abs(asso1 - asso2)[lower.tri(asso1)]
  diffCounts <- data.frame(
    count = c(sum(diffVec != 0),
              if (!is.null(object$diffMat)) .countDiffnetEdges(object$diffMat) else NA,
              if (!is.null(object$diffAdjustMat)) {
                .countDiffnetEdges(object$diffAdjustMat)
              } else {
                NA
              }),
    row.names = c("Non-zero association differences",
                  "Significant before adjustment",
                  "Significant after adjustment")
  )

  pvalSummary <- NULL
  if (!is.null(object$pvalsVec)) {
    pvalSummary <- rbind(
      "Raw p-values" = .summaryDiffnetPvals(object$pvalsVec,
                                            alpha,
                                            digitsPval),
      "Adjusted p-values" = .summaryDiffnetPvals(
        object$pAdjustVec,
        if (identical(adjust, "lfdr")) lfdrThresh else alpha,
        digitsPval
      )
    )
  }

  extra <- NULL
  if (identical(object$diffMethod, "threshold")) {
    extra <- data.frame(value = object$diffThresh,
                        row.names = "Difference threshold")
  } else if (identical(object$diffMethod, "discordant") &&
             !is.null(object$diffProbs)) {
    extra <- data.frame(t(.summaryDiffnetNumeric(
      object$diffProbs[lower.tri(object$diffProbs)], digits
    )))
    rownames(extra) <- "Differential probabilities"
  }

  out <- list(meta = meta,
              assoSummary = assoSummary,
              counts = counts,
              diffCounts = diffCounts,
              pvalSummary = pvalSummary,
              extra = extra,
              call = object$call)

  class(out) <- "summary.diffnet"

  out
}

#' @title Print Summary for objects of class \code{diffnet}
#'
#' @param x object of class \code{summary.diffnet} returned by
#'   \code{\link{summary.diffnet}}.
#' @param ... not used.
#'
#' @method print summary.diffnet
#' @rdname summarize.diffnet
#' @export
print.summary.diffnet <- function(x, ...) {
  cat("\nDifferential Network Summary\n")
  cat("----------------------------")

  cat("\nCALL: \n")
  dput(x$call, control = NULL)

  cat("\n______________________________")
  cat("\nMethod and tests\n")
  cat("````````````````\n")
  .printDiffnetKeyValue(x$meta)

  cat("\n______________________________")
  cat("\nAssociation summaries\n")
  cat("`````````````````````\n")
  print(x$assoSummary)

  cat("\n______________________________")
  cat("\nNetwork counts\n")
  cat("``````````````\n")
  print(x$counts)

  cat("\n______________________________")
  cat("\nDifferential associations\n")
  cat("````````````````````````\n")
  print(x$diffCounts)

  if (!is.null(x$pvalSummary)) {
    cat("\n______________________________")
    cat("\nP-value summaries\n")
    cat("`````````````````\n")
    print(x$pvalSummary)
  }

  if (!is.null(x$extra)) {
    cat("\n______________________________")
    cat("\nMethod-specific summary\n")
    cat("```````````````````````\n")
    print(x$extra)
  }

  invisible(x)
}

#' @keywords internal
.summaryDiffnetNumeric <- function(x, digits) {
  round(c(min = min(x, na.rm = TRUE),
          q1 = unname(stats::quantile(x, 0.25, na.rm = TRUE)),
          median = stats::median(x, na.rm = TRUE),
          mean = mean(x, na.rm = TRUE),
          q3 = unname(stats::quantile(x, 0.75, na.rm = TRUE)),
          max = max(x, na.rm = TRUE)), digits)
}

#' @keywords internal
.summaryDiffnetPvals <- function(x, threshold, digitsPval) {
  x <- x[!is.na(x)]

  if (!length(x)) {
    return(data.frame(n = 0,
                      min = NA_character_,
                      q1 = NA_character_,
                      median = NA_character_,
                      mean = NA_character_,
                      q3 = NA_character_,
                      max = NA_character_,
                      significant = 0,
                      check.names = FALSE))
  }

  data.frame(n = length(x),
             min = .formatPval(min(x, na.rm = TRUE), digitsPval),
             q1 = .formatPval(unname(stats::quantile(x, 0.25, na.rm = TRUE)),
                              digitsPval),
             median = .formatPval(stats::median(x, na.rm = TRUE), digitsPval),
             mean = .formatPval(mean(x, na.rm = TRUE), digitsPval),
             q3 = .formatPval(unname(stats::quantile(x, 0.75, na.rm = TRUE)),
                              digitsPval),
             max = .formatPval(max(x, na.rm = TRUE), digitsPval),
             significant = sum(x <= threshold),
             check.names = FALSE)
}

#' @keywords internal
.countDiffnetEdges <- function(mat) {
  sum(mat[lower.tri(mat)] != 0, na.rm = TRUE)
}

#' @keywords internal
.countDiffnetSingletons <- function(mat) {
  matTmp <- mat
  diag(matTmp) <- 0
  sum(rowSums(matTmp != 0, na.rm = TRUE) == 0)
}

#' @keywords internal
.printDiffnetKeyValue <- function(x) {
  labels <- rownames(x)
  values <- as.character(x[[1]])
  labelWidth <- max(nchar(labels))
  valueWidth <- max(nchar(values), na.rm = TRUE)
  
  cat(paste0(format(labels, width = labelWidth, justify = "left"),
             " ",
             format(values, width = valueWidth, justify = "right")),
      sep = "\n")
  cat("\n")
}
