# G.2 -- glmnet drift guard.  Pure-glmnet tests (no oalasso code) of the
# two upstream behaviors the exact-objective correction (D.9) relies on:
#   (i)  penalty factors are internally rescaled to sum to nvars, so
#        s = mean(pen) * lambda_n / n imposes the TARGET objective
#        -ll + lambda_n * sum(pen_j |alpha_j|) exactly (KKT-verified);
#   (ii) an Inf penalty factor is converted to an exclusion: the
#        coefficient is EXACTLY zero.
# If these tests fail after a glmnet upgrade, the s-correction (and the
# Inf -> exclude handling of user outcome.coef zeros) must be revisited.
#
# NOTE (probed, glmnet 4.1.10): when an Inf penalty factor is PRESENT in
# the supplied vector, glmnet resets it to 1 internally and the sum-to-
# nvars rescale then runs over all nvars entries, so the exact constant
# becomes (sum(pen_finite) + n_inf)/nvars * lambda_n / n -- NOT
# mean(pen_finite) * lambda_n / n.  Fitting on the REDUCED matrix
# (finite-penalty columns only) restores the mean(pen_active) formula.
# The KKT tests in test-kkt.R are the acceptance criterion for whichever
# route the package takes.

glmnet_info <- function(msg) {
  sprintf("glmnet %s: %s", as.character(utils::packageVersion("glmnet")), msg)
}

test_that("s = mean(pen)*lambda_n/n imposes the exact OAL objective (KKT + pinned solution)", {
  skip_if_not_installed("glmnet")
  d <- make_oal_data()
  X <- stats::model.matrix(oal_formula(), d)[, -1]
  A <- d$A
  n <- nrow(X)
  Xs <- scale(X, TRUE, TRUE)
  b <- stats::coef(stats::lm(Y ~ ., data = data.frame(Y = d$y, A = A,
                                                      as.data.frame(Xs))))[colnames(Xs)]

  gam <- 2
  lam_n <- n^0.25
  pen <- abs(b)^(-gam)
  fit <- glmnet::glmnet(x = Xs, y = A, family = "binomial", alpha = 1,
                        standardize = FALSE, intercept = TRUE,
                        penalty.factor = pen)
  s <- mean(pen) * lam_n / n
  cf <- as.numeric(stats::coef(fit, s = s, exact = TRUE, x = Xs, y = A,
                               penalty.factor = pen))

  k <- kkt_max_violation(cf[1], cf[-1], Xs, A, lam_n, pen)
  expect_true(k$coef < KKT_TOL,
              info = glmnet_info("sum-to-nvars rescaling / coef(exact = TRUE) assumption drifted (KKT violated)"))
  expect_true(k$intercept < KKT_INT_TOL,
              info = glmnet_info("unpenalized-intercept stationarity drifted"))

  # Hard-coded exact solution for this fixture at (delta, gamma) = (0.25, 2),
  # precomputed with thresh = 1e-12 (default-thresh solution agrees to
  # 2.8e-8).  Order: intercept, Xc1, Xc2, Xp1, Xp2, Xi1, Xi2, Xs1, Xs2.
  pinned <- c(0.096835877263, 0.183505790942, 0.320052820281,
              0.000000000000, -0.115261558970, 0.000000000000,
              0.000000000000, 0.000000000000, 0.000000000000)
  expect_true(max(abs(cf - pinned)) < 1e-5,
              info = glmnet_info("exact-solution regression pin drifted"))
})

test_that("Inf penalty factors are converted to exclusions (coefficient exactly zero)", {
  skip_if_not_installed("glmnet")
  d <- make_oal_data()
  X <- stats::model.matrix(oal_formula(), d)[, -1]
  A <- d$A
  n <- nrow(X)
  Xs <- scale(X, TRUE, TRUE)
  b <- stats::coef(stats::lm(Y ~ ., data = data.frame(Y = d$y, A = A,
                                                      as.data.frame(Xs))))[colnames(Xs)]
  lam_n <- n^0.25
  pen <- abs(b)^(-2)
  pen[3] <- Inf                       # exclude the 3rd column (Xp1)
  fin <- is.finite(pen)

  fit <- glmnet::glmnet(x = Xs, y = A, family = "binomial", alpha = 1,
                        standardize = FALSE, intercept = TRUE,
                        penalty.factor = pen)
  cf <- as.numeric(stats::coef(fit, s = mean(pen[fin]) * lam_n / n,
                               exact = TRUE, x = Xs, y = A,
                               penalty.factor = pen))
  expect_identical(cf[1 + 3], 0,
                   info = glmnet_info("Inf penalty factor no longer converts to exclude"))

  # Reduced-matrix route: dropping the excluded column restores the
  # mean(pen_active) formula exactly (KKT on the target objective).
  fit2 <- glmnet::glmnet(x = Xs[, fin], y = A, family = "binomial",
                         alpha = 1, standardize = FALSE, intercept = TRUE,
                         penalty.factor = pen[fin])
  s2 <- mean(pen[fin]) * lam_n / n
  cf2 <- as.numeric(stats::coef(fit2, s = s2, exact = TRUE, x = Xs[, fin],
                                y = A, penalty.factor = pen[fin]))
  k2 <- kkt_max_violation(cf2[1], cf2[-1], Xs[, fin], A, lam_n, pen[fin])
  expect_true(k2$coef < KKT_TOL,
              info = glmnet_info("mean(pen_active) formula on the reduced matrix drifted"))
})
