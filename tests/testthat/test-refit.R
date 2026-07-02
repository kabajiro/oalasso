# G.10 / D.11 -- the Schnitzer post-selection refit variant.

test_that("refit = TRUE: PS equals the unpenalized glm refit on the selected set", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, refit = TRUE)
  sel <- fit$selected$term[fit$selected$selected]
  expect_gt(length(sel), 0L)

  # Unpenalized logistic fitted probabilities are invariant to the affine
  # standardization of the columns, so the refit can be reproduced on the
  # ORIGINAL-scale covariate matrix.
  g <- stats::glm(fit$treat ~ fit$covs[, sel, drop = FALSE],
                  family = stats::binomial())
  ps_ref <- pmin(pmax(stats::fitted(g), fit$clip[1]), fit$clip[2])
  expect_equal(unname(fit$ps), unname(ps_ref), tolerance = 1e-6)
})

test_that("refit PS differs generically from the penalized PS", {
  d <- make_oal_data()
  fit_pen <- oal(oal_formula(), data = d, outcome = ~ y, lambda = 0.25,
                 gamma = 2, refit = FALSE)
  fit_ref <- oal(oal_formula(), data = d, outcome = ~ y, lambda = 0.25,
                 gamma = 2, refit = TRUE)
  # same single grid point, same selection stage; the penalized PS is
  # shrunk while the refit is not
  expect_false(isTRUE(all.equal(unname(fit_pen$ps), unname(fit_ref$ps),
                                tolerance = 1e-6)))
})

test_that("refit convergence is recorded per grid point", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, refit = TRUE)
  expect_true(is.logical(fit$path$refit.converged))
  expect_false(anyNA(fit$path$refit.converged))

  fit0 <- oal(oal_formula(), data = d, outcome = ~ y, refit = FALSE)
  expect_true(all(is.na(fit0$path$refit.converged)))
})

test_that("refit = TRUE is provenance-labeled", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25),
             refit = TRUE)
  expect_match(fit$provenance, "post-selection refit", fixed = TRUE)
})

test_that("refit = TRUE wAMD is computed from the refit PS", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, refit = TRUE)
  res <- oal_wamd(fit$ps, fit$treat, std_covs(fit), coef = fit$outcome.coef,
                  estimand = fit$estimand, clip = fit$clip)
  expect_equal(fit$criterion.value, res$total, tolerance = 1e-10)
})
