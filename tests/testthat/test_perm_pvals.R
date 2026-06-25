test_that(".calcPermPval returns active, empirical, and refined p-values", {
  res <- NetCoMi:::.calcPermPval(
    tstat = c(a = 3, b = 4),
    tstatPerm = matrix(c(1, 3, 5, 4, 4, 2), nrow = 3),
    nPerm = 3,
    exact = FALSE
  )

  expect_named(res, c("pval", "empirical", "refined"))
  expect_equal(res$pval, c(a = 0.75, b = 0.75))
  expect_equal(res$empirical, c(a = 0.75, b = 0.75))
  expect_null(res$refined)
})

test_that(".calcPermPval matrix empirical p-values match column-wise calls", {
  tstat <- c(a = 3, b = 4)
  tstatPerm <- matrix(c(1, 3, 5, 4, 4, 2), nrow = 3)

  matrix_res <- NetCoMi:::.calcPermPval(
    tstat = tstat,
    tstatPerm = tstatPerm,
    nPerm = 3,
    exact = FALSE
  )

  col_res <- vapply(seq_along(tstat), function(i) {
    NetCoMi:::.calcPermPval(
      tstat = tstat[i],
      tstatPerm = tstatPerm[, i],
      nPerm = 3,
      exact = FALSE
    )$pval
  }, numeric(1))
  names(col_res) <- names(tstat)

  expect_equal(matrix_res$empirical, col_res)
})

test_that(".calcPermPval does not refine exact permutation tests", {
  res <- NetCoMi:::.calcPermPval(
    tstat = 3,
    tstatPerm = c(1, 3, 5),
    nPerm = 3,
    exact = TRUE,
    refinePermPvals = TRUE
  )

  expect_equal(res$pval, 2 / 3)
  expect_equal(res$empirical, 2 / 3)
  expect_null(res$refined)
})
