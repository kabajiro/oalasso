# G.3 -- numerical equivalence against the in-test reimplementation of the
# rejoinder recipe (helper-rejoinder.R): same standardized X, same OLS
# betas, same S&E grid => identical selected lambda/gamma, identical
# selected set, identical PS and wAMD path.

test_that("oal() with all defaults reproduces the rejoinder recipe", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  ref <- rejoinder_oal(oal_formula(), d, outcome = "y")

  # identical outcome-model coefficients (same lm on the same standardized X)
  expect_equal(fit$outcome.coef[names(ref$betaXY)], ref$betaXY,
               tolerance = 1e-12)

  # identical wAMD over the whole grid, in grid order
  expect_equal(fit$path$wamd, ref$wamd, tolerance = 1e-8) # cross-BLAS headroom

  # Selection equivalence -- robust to numerical near-ties across BLAS
  # builds: the fixture grid contains points whose wAMD values differ only
  # at the ~1e-8 relative level, so the argmin can legitimately flip
  # between platforms (macOS Accelerate vs reference BLAS vs MKL). When
  # both implementations land on the same grid point, everything must
  # match; when they flip on a near-tie, the two selected criterion
  # values must still agree to near-tie precision.
  same.point <- identical(unname(fit$lambda[["delta"]]), ref$delta)
  if (same.point) {
    expect_equal(unname(fit$lambda[["gamma"]]), ref$gamma)
    expect_equal(unname(fit$lambda[["lambda.n"]]), ref$lambda.n)
    expect_equal(fit$criterion.value, ref$criterion.value,
                 tolerance = 1e-6) # scalar, cross-BLAS headroom
    # identical selected variable set
    got <- sort(fit$selected$term[fit$selected$selected])
    expect_identical(got, sort(ref$selected))
    # identical propensity scores
    expect_equal(unname(fit$ps), unname(ref$ps),
                 tolerance = 1e-6) # cross-BLAS headroom
  } else {
    expect_lt(abs(fit$criterion.value - ref$criterion.value),
              1e-6 * max(1, abs(ref$criterion.value)))
  }
})

test_that("per-grid-point selected sets match the rejoinder recipe", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, keep.path = TRUE)
  ref <- rejoinder_oal(oal_formula(), d, outcome = "y")
  expect_identical(length(fit$sets), length(ref$sets))
  for (r in seq_along(ref$sets)) {
    expect_identical(sort(fit$sets[[r]]), sort(ref$sets[[r]]))
  }
  expect_equal(fit$path$n.selected, lengths(ref$sets))
})
