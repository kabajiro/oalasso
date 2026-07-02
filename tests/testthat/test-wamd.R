# G.5 / D.12 -- the exported wAMD criterion, hand-checked on the 6-row
# fixture (full derivation in helper-fixtures.R), plus the IPW weight
# formulas and a WeightIt cross-check.

test_that("oal_wamd ATE matches the hand-computed value (with clipping)", {
  fx <- make_wamd_fixture()
  res <- oal_wamd(fx$ps, fx$treat, fx$covs, coef = fx$coef,
                  estimand = "ATE")
  # From the derivation header:
  #   diff_x1 = 112/59, diff_x2 = 5/413
  #   contributions = (2*112/59, 4*5/413) = (1568/413, 20/413)
  #   total = 1588/413
  expect_type(res, "list")
  expect_named(res, c("total", "by.covariate"))
  expect_equal(res$total, 1588 / 413, tolerance = 1e-12)
  expect_equal(unname(res$by.covariate), c(1568 / 413, 20 / 413),
               tolerance = 1e-12)
  expect_equal(sum(res$by.covariate), res$total, tolerance = 1e-12)
})

test_that("oal_wamd ATT matches the hand-computed value (with clipping)", {
  fx <- make_wamd_fixture()
  res <- oal_wamd(fx$ps, fx$treat, fx$covs, coef = fx$coef,
                  estimand = "ATT")
  # From the derivation header:
  #   diff_x1 = 392/401, diff_x2 = 386/1203
  #   contributions = (784/401, 1544/1203) = (2352/1203, 1544/1203)
  #   total = 3896/1203
  expect_equal(res$total, 3896 / 1203, tolerance = 1e-12)
  expect_equal(unname(res$by.covariate), c(784 / 401, 1544 / 1203),
               tolerance = 1e-12)
})

test_that("oal_wamd clips the PS BEFORE forming the IPW weights", {
  fx <- make_wamd_fixture()
  clipped <- pmin(pmax(fx$ps, 0.01), 0.99)
  for (est in c("ATE", "ATT")) {
    raw <- oal_wamd(fx$ps, fx$treat, fx$covs, coef = fx$coef, estimand = est)
    pre <- oal_wamd(clipped, fx$treat, fx$covs, coef = fx$coef, estimand = est)
    expect_equal(raw$total, pre$total, tolerance = 1e-12)
  }
  # unclipped weights would give a different number: recompute without clip
  e <- fx$ps
  w <- fx$treat / e + (1 - fx$treat) / (1 - e)
  d1 <- abs(sum(w * fx$covs[, 1] * fx$treat) / sum(w * fx$treat) -
            sum(w * fx$covs[, 1] * (1 - fx$treat)) / sum(w * (1 - fx$treat)))
  d2 <- abs(sum(w * fx$covs[, 2] * fx$treat) / sum(w * fx$treat) -
            sum(w * fx$covs[, 2] * (1 - fx$treat)) / sum(w * (1 - fx$treat)))
  unclipped <- 2 * d1 + 4 * d2
  clip_res <- oal_wamd(fx$ps, fx$treat, fx$covs, coef = fx$coef,
                       estimand = "ATE")
  expect_false(isTRUE(all.equal(clip_res$total, unclipped,
                                tolerance = 1e-8)))
})

test_that("wAMD uses |coef|: sign of the outcome coefficients is irrelevant", {
  fx <- make_wamd_fixture()
  res1 <- oal_wamd(fx$ps, fx$treat, fx$covs, coef = fx$coef,
                   estimand = "ATE")
  res2 <- oal_wamd(fx$ps, fx$treat, fx$covs, coef = abs(fx$coef),
                   estimand = "ATE")
  expect_equal(res1$total, res2$total, tolerance = 1e-12)
})

test_that("$weights implements the D.12 IPW formulas at the fitted estimand", {
  d <- make_oal_data()
  fit_ate <- oal(oal_formula(), data = d, outcome = ~ y,
                 lambda = c(-1, 0.25), estimand = "ATE")
  expect_equal(unname(fit_ate$weights),
               unname(ref_weights(fit_ate$ps, fit_ate$treat, "ATE")),
               tolerance = 1e-12)
  fit_att <- oal(oal_formula(), data = d, outcome = ~ y,
                 lambda = c(-1, 0.25), estimand = "ATT")
  expect_equal(unname(fit_att$weights),
               unname(ref_weights(fit_att$ps, fit_att$treat, "ATT")),
               tolerance = 1e-12)
  expect_equal(weights(fit_ate), fit_ate$weights)
})

test_that("IPW weights agree with WeightIt::get_w_from_ps", {
  skip_if_not_installed("WeightIt")
  d <- make_oal_data()
  for (est in c("ATE", "ATT")) {
    fit <- oal(oal_formula(), data = d, outcome = ~ y,
               lambda = c(-1, 0.25), estimand = est)
    w_ref <- WeightIt::get_w_from_ps(unname(fit$ps), treat = fit$treat,
                                     estimand = est)
    expect_equal(unname(fit$weights), unname(as.numeric(w_ref)),
                 tolerance = 1e-10)
  }
})

test_that("criterion.value equals oal_wamd on the standardized design", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  expect_identical(fit$criterion, "wamd")
  res <- oal_wamd(fit$ps, fit$treat, std_covs(fit), coef = fit$outcome.coef,
                  estimand = fit$estimand, clip = fit$clip)
  expect_equal(fit$criterion.value, res$total, tolerance = 1e-10)
})
