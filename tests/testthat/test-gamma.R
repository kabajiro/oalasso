# G.4 / D.7 -- lambda/gamma grid construction.
#
# Hand derivation of the paired gammas (S&E / rejoinder / fOAL convention,
# gamma = 2*(gamma.factor - delta + 1), gcf = 2), for the default delta grid
# in ascending (grid) order:
#   delta: -10    -5    -2    -1   -0.75  -0.5  -0.25  0.25  0.49
#   gamma:  26    16    10     8    7.5    7     6.5   5.5   5.02
# The fAL transcription bug 2*(gcf - delta) + 1 would instead give
#   25, 15, 9, 7, 6.5, 6, 5.5, 4.5, 4.02  (always smaller by 1) and must
# never appear.

test_that("paired mode: gamma = 2*(gamma.factor - delta + 1) over the default grid", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  deltas <- default_deltas()
  expect_equal(fit$path$delta, deltas)
  expect_equal(fit$path$gamma, 2 * (2 - deltas + 1))
  expect_equal(fit$path$gamma, c(26, 16, 10, 8, 7.5, 7, 6.5, 5.5, 5.02))
  expect_identical(fit$info$gamma.mode, "paired")
  # lambda.n = n^delta with n = nrow(data) (kills legacy scoping bug #2)
  expect_equal(fit$path$lambda.n, nrow(d)^deltas)
})

test_that("gamma.factor propagates into the pairing", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, gamma.factor = 3,
             lambda = c(-1, -0.5, 0.25))
  expect_equal(fit$path$gamma, 2 * (3 - c(-1, -0.5, 0.25) + 1))
})

test_that("the fAL transcription-bug formula never appears in the grid", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  fal <- 2 * (2 - fit$path$delta) + 1
  expect_true(all(abs(fit$path$gamma - fal) > 0.5))
})

test_that("scalar gamma is fixed and crossed with the full lambda grid", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, gamma = 2.5)
  deltas <- default_deltas()
  expect_equal(fit$path$delta, deltas)
  expect_equal(fit$path$gamma, rep(2.5, length(deltas)))
  expect_identical(fit$info$gamma.mode, "fixed")
  # the selected point reports the fixed gamma
  expect_equal(unname(fit$lambda[["gamma"]]), 2.5)
})

test_that("the selected point is a row of the grid", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  expect_true(unname(fit$lambda[["delta"]]) %in% fit$path$delta)
  row <- fit$path[fit$path$selected, , drop = FALSE]
  expect_identical(nrow(row), 1L)
  expect_equal(row$delta, unname(fit$lambda[["delta"]]))
  expect_equal(row$gamma, unname(fit$lambda[["gamma"]]))
  expect_equal(row$wamd, fit$criterion.value)
})
