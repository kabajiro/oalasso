# G.1 -- THE core correctness test (ported from the rejoinder repo's
# tests/test-glmnet.R pattern).  For every grid-point type, the returned
# penalized coefficients on the standardized design must satisfy the KKT
# subgradient conditions of the TARGET objective
#     -ll(alpha) + lambda_n * sum_j |b_j|^(-gamma) |alpha_j|
# on the TOTAL-log-likelihood scale:
#     |grad_j| <= lambda_n * pen_j            for alpha_j  = 0
#     grad_j = -lambda_n * pen_j * sign(...)  for alpha_j != 0
# This is the acceptance criterion for the lambda-scale correction
# s = mean(pen_active) * lambda_n / n (D.9): whatever constant satisfies
# these conditions is correct, whatever fails is wrong -- independent of
# glmnet's internal rescaling conventions.
#
# Tolerances are calibrated in helper-fixtures.R (observed <= 4.3e-4 with
# glmnet 4.1.10 defaults; wrong scale constants violate at >= 1e-1).

test_that("plain OAL, paired mode, single grid point satisfies KKT", {
  d <- make_oal_data()
  # single delta = 0.25 -> paired gamma = 2*(2 - 0.25 + 1) = 5.5
  fit <- oal(oal_formula(), data = d, outcome = ~ y, lambda = 0.25,
             refit = FALSE)
  expect_equal(unname(fit$lambda[["gamma"]]), 5.5)
  lam_n <- unname(fit$lambda[["lambda.n"]])
  expect_equal(lam_n, nrow(d)^0.25)

  Xs <- std_covs(fit)
  b <- fit$outcome.coef[colnames(Xs)]
  pen <- abs(b)^(-unname(fit$lambda[["gamma"]]))
  a <- fit$coefficients.std
  k <- kkt_max_violation(a[1], a[-1][colnames(Xs)], Xs, fit$treat, lam_n, pen)
  expect_lt(k$coef, KKT_TOL)
  expect_lt(k$intercept, KKT_INT_TOL)
})

test_that("plain OAL, fixed-gamma mode, selected point satisfies KKT", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, gamma = 2,
             refit = FALSE)
  lam_n <- unname(fit$lambda[["lambda.n"]])
  Xs <- std_covs(fit)
  pen <- abs(fit$outcome.coef[colnames(Xs)])^(-2)
  a <- fit$coefficients.std
  k <- kkt_max_violation(a[1], a[-1][colnames(Xs)], Xs, fit$treat, lam_n, pen)
  expect_lt(k$coef, KKT_TOL)
  expect_lt(k$intercept, KKT_INT_TOL)
})

test_that("GOAL grid point (lambda2 > 0) satisfies KKT on the augmented objective", {
  d <- make_oal_data()
  lam2 <- 0.5
  # single grid point: delta = 0.25, fixed gamma = 2, lambda2 = 0.5
  fit <- oal(oal_formula(), data = d, outcome = ~ y, method = "goal",
             lambda = 0.25, gamma = 2, lambda2 = lam2, refit = FALSE)
  expect_equal(unname(fit$lambda[["lambda2"]]), lam2)
  lam_n <- unname(fit$lambda[["lambda.n"]])

  Xs <- std_covs(fit)
  p <- ncol(Xs)
  pen <- abs(fit$outcome.coef[colnames(Xs)])^(-2)

  # GOAL raw penalty constant is on the AUGMENTED-n base (Balde 2025
  # supplement: adaptive.lasso(lambda = n.q^(il)) with n.q = n + q, q = p
  # augmentation rows) -- NOT n^delta.
  expect_equal(lam_n, (nrow(d) + p)^0.25)

  # Zou-Hastie augmentation exactly as oal_funs.R (D.10): the augmented-
  # data lasso solution is the returned standardized coefficient vector
  # UNDONE by the (1 + lambda2) rescale, which applies to ALL coefficients
  # INCLUDING the intercept (Balde 2025 supplement:
  # expit(cbind(1, X) %*% (1 + lambda2) * coef)).
  Xa <- rbind(Xs, sqrt(lam2) * diag(p))
  Aa <- c(fit$treat, rep(0, p))
  a <- fit$coefficients.std / (1 + lam2)
  alpha_aug <- a[-1][colnames(Xs)]
  k <- kkt_max_violation(a[1], alpha_aug, Xa, Aa, lam_n, pen)
  expect_lt(k$coef, KKT_TOL)
  expect_lt(k$intercept, KKT_INT_TOL)
})

test_that("user outcome.coef with a zero entry: hard exclusion + KKT on the rest", {
  d <- make_oal_data()
  bv <- c(Xc1 = 0.9, Xc2 = 0.9, Xp1 = 0.7, Xp2 = 0.7,
          Xi1 = -0.05, Xi2 = -0.05, Xs1 = 0, Xs2 = 0.06)
  fit <- oal(oal_formula(), data = d, outcome.coef = bv, lambda = 0.25,
             gamma = 2, refit = FALSE)
  lam_n <- unname(fit$lambda[["lambda.n"]])
  Xs <- std_covs(fit)
  pen <- abs(bv[colnames(Xs)])^(-2)      # Xs1 -> Inf (glmnet exclude)
  expect_true(is.infinite(pen[["Xs1"]]))

  a <- fit$coefficients.std
  alpha <- a[-1][colnames(Xs)]
  # hard drop: exactly zero, at the exact solution
  expect_identical(unname(alpha[["Xs1"]]), 0)
  # KKT over the finite-penalty coordinates (the target objective fixes
  # the excluded coordinate at zero). The s-correction uses the effective
  # mean (sum(pen[finite]) + #Inf)/p: glmnet converts Inf factors to
  # `exclude` and internally resets them to 1 BEFORE the sum-to-nvars
  # rescaling, so each excluded column contributes 1 to the mean on the
  # full matrix; a finite-only mean fails these conditions at >= 1e-1.
  k <- kkt_max_violation(a[1], alpha, Xs, fit$treat, lam_n, pen)
  expect_lt(k$coef, KKT_TOL)
  expect_lt(k$intercept, KKT_INT_TOL)
})

test_that("GOAL with lambda2 = 0 satisfies the PLAIN (unaugmented) KKT conditions", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, method = "goal",
             lambda = 0.25, gamma = 2, lambda2 = 0, refit = FALSE)
  lam_n <- unname(fit$lambda[["lambda.n"]])
  Xs <- std_covs(fit)
  pen <- abs(fit$outcome.coef[colnames(Xs)])^(-2)
  a <- fit$coefficients.std
  # no augmentation may run at lambda2 = 0 (D.10): the plain-objective
  # KKT conditions must hold on the ORIGINAL n rows
  k <- kkt_max_violation(a[1], a[-1][colnames(Xs)], Xs, fit$treat, lam_n, pen)
  expect_lt(k$coef, KKT_TOL)
  expect_lt(k$intercept, KKT_INT_TOL)
})
