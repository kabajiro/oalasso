# G.9 / D.5 -- the outcome.coef extension hook: validation, override
# semantics, zero -> hard exclusion, provenance.

hook_coef <- function() {
  c(Xc1 = 0.9, Xc2 = 0.9, Xp1 = 0.7, Xp2 = 0.7,
    Xi1 = -0.05, Xi2 = -0.05, Xs1 = 0.04, Xs2 = 0.06)
}

test_that("unnamed outcome.coef errors", {
  d <- make_oal_data()
  expect_error(oal(oal_formula(), data = d,
                   outcome.coef = unname(hook_coef())), "name")
})

test_that("outcome.coef with missing names errors listing the expected names", {
  d <- make_oal_data()
  bv <- hook_coef()[-8]                     # drop Xs2
  # D.5: the error lists the expected names (the model-matrix columns)
  expect_error(oal(oal_formula(), data = d, outcome.coef = bv), "Xs2")
})

test_that("outcome.coef of the wrong length errors", {
  d <- make_oal_data()
  expect_error(oal(oal_formula(), data = d,
                   outcome.coef = c(Xc1 = 0.9, Xc2 = 0.9)), "Xp1|name|length")
})

test_that("non-finite outcome.coef entries error", {
  d <- make_oal_data()
  bv <- hook_coef(); bv["Xp1"] <- NA_real_
  expect_error(oal(oal_formula(), data = d, outcome.coef = bv),
               "finite|NA|[Mm]issing")
  bv["Xp1"] <- Inf
  expect_error(oal(oal_formula(), data = d, outcome.coef = bv), "finite")
})

test_that("all-zero outcome.coef errors", {
  d <- make_oal_data()
  bv <- hook_coef(); bv[] <- 0
  expect_error(oal(oal_formula(), data = d, outcome.coef = bv), "zero")
})

test_that("a valid hook overrides the internal outcome model (no `outcome` needed)", {
  d <- make_oal_data()
  bv <- hook_coef()
  fit <- oal(oal_formula(), data = d, outcome.coef = bv,
             lambda = c(-1, 0.25))
  expect_s3_class(fit, "oal")
  expect_identical(fit$info$gamma.mode, "user-coef")
  expect_equal(fit$outcome.coef[names(bv)], bv, tolerance = 1e-15)
  expect_match(fit$outcome.model$label, "user")

  # supplied alongside `outcome`, the outcome variable is unused: the fit
  # is identical to the fit without it
  fit2 <- oal(oal_formula(), data = d, outcome = ~ y, outcome.coef = bv,
              lambda = c(-1, 0.25))
  expect_equal(fit2$ps, fit$ps, tolerance = 1e-12)

  # and it genuinely overrides: differs from the internal-model fit
  fit_int <- oal(oal_formula(), data = d, outcome = ~ y,
                 lambda = c(-1, 0.25))
  expect_false(isTRUE(all.equal(unname(fit$outcome.coef),
                                unname(fit_int$outcome.coef),
                                tolerance = 1e-6)))
})

test_that("scrambled outcome.coef names are accepted and aligned to the columns", {
  d <- make_oal_data()
  bv <- hook_coef()
  fit1 <- oal(oal_formula(), data = d, outcome.coef = bv,
              lambda = c(-1, 0.25))
  fit2 <- oal(oal_formula(), data = d, outcome.coef = rev(bv),
              lambda = c(-1, 0.25))
  expect_equal(fit2$ps, fit1$ps, tolerance = 1e-12)
})

test_that("zero entries are hard exclusions across the whole grid", {
  d <- make_oal_data()
  bv <- hook_coef(); bv["Xs1"] <- 0
  fit <- oal(oal_formula(), data = d, outcome.coef = bv, keep.path = TRUE)
  srow <- fit$selected[fit$selected$term == "Xs1", ]
  expect_false(srow$selected)
  expect_identical(srow$coef, 0)
  expect_identical(srow$role, "excluded")
  expect_true(is.infinite(srow$penalty.factor))
  # exactly zero at EVERY grid point (glmnet Inf -> exclude, no epsilon)
  expect_true(all(fit$coef.path["Xs1", ] == 0))
  expect_false(any(vapply(fit$sets, function(s) "Xs1" %in% s, logical(1))))
})

test_that("hook fits carry the experimental provenance label", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome.coef = hook_coef(),
             lambda = c(-1, 0.25))
  expect_match(fit$provenance, "user-supplied outcome coefficients",
               fixed = TRUE)
  expect_match(fit$provenance, "screening use only", fixed = TRUE)
  expect_output(print(fit), "screening use only", fixed = TRUE)
})
