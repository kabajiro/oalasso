# ==========================================================================
# Shared fixtures and reference helpers for the oalasso test suite
# ==========================================================================
#
# (a) make_oal_data(): deterministic S&E-style simulated dataset
#     (Shortreed & Ertefaie 2017 scenario layout): 2 confounders (Xc*,
#     predict A and y), 2 pure outcome predictors (Xp*), 2 instruments
#     (Xi*, predict A only), 2 noise covariates (Xs*).  The seed is fixed
#     INSIDE the helper so every test sees the same data.
#
#     Realized properties with seed = 20260702, n = 300 (verified
#     independently of the package, R 4.5.0 / glmnet 4.1.10):
#       * arms: 157 treated / 143 control (>= 2 per arm);
#       * standardized-scale OLS betas of lm(y ~ A + Xs):
#           Xc1  0.9577, Xc2  0.9120, Xp1 0.7258, Xp2 0.7303,
#           Xi1 -0.0487, Xi2 -0.0531, Xs1 0.0416, Xs2 0.0661
#         -- all nonzero, so the internal degenerate-beta error never
#         fires, and min |b| = 0.0416 keeps |b|^-gamma finite up to the
#         largest paired gamma (gamma = 26 at delta = -10:
#         0.0416^-26 ~ 8.7e35 < .Machine$double.xmax);
#       * unpenalized glm PS range (0.0077, 0.9950), i.e. a handful of
#         units sit at the default clip bounds -- well under the 5%
#         near-positivity threshold.
#
# (b) make_wamd_fixture(): a 6-row HAND-COMPUTABLE fixture for the exported
#     oal_wamd() criterion.  Full derivation below; nothing was copied
#     from running the package.
#
# --------------------------------------------------------------------------
# 6-ROW wAMD FIXTURE -- FULL HAND DERIVATION (D.12 formulas)
# --------------------------------------------------------------------------
# Units u1..u6:
#   A   = (1, 1, 1, 0, 0, 0)
#   x1  = (1, 2, 3, 3, 2, 1)
#   x2  = (0, 1, 1, 0, 0, 1)
#   ps  = (0.5, 0.8, 0.001, 0.5, 0.2, 0.998)     (raw candidate PS)
#   b   = (x1 = 2, x2 = -4)                      (signed coef; wAMD uses |b|)
#
# Clipping to [0.01, 0.99] BEFORE the IPW weights (D.12):
#   e = (0.5, 0.8, 0.01, 0.5, 0.2, 0.99)   (u3: 0.001 -> 0.01; u6: 0.998 -> 0.99)
#
# --- ATE: w_i = A_i/e_i + (1-A_i)/(1-e_i) ---
#   treated  (u1,u2,u3): 1/e     = (2, 5/4, 100),  sum = 413/4  (= 103.25)
#   control  (u4,u5,u6): 1/(1-e) = (2, 5/4, 100),  sum = 413/4
#   x1: treated weighted sum = 2*1 + 5/4*2 + 100*3 = 609/2
#       control weighted sum = 2*3 + 5/4*2 + 100*1 = 217/2
#       diff_x1 = (609/2 - 217/2) / (413/4) = 196/(413/4) = 784/413 = 112/59
#   x2: treated weighted sum = 0 + 5/4 + 100 = 405/4
#       control weighted sum = 0 + 0   + 100 = 100
#       diff_x2 = (405/4 - 400/4) / (413/4) = (5/4)/(413/4) = 5/413
#   contributions: |b_x1|*diff = 2*112/59 = 224/59 = 1568/413
#                  |b_x2|*diff = 4*5/413  = 20/413
#   total_ATE = 1588/413 = 3.845036319612591
#
# --- ATT: w_i = A_i + (1-A_i) e_i/(1-e_i) ---
#   treated: (1, 1, 1), sum = 3
#   control: e/(1-e) = (1, 1/4, 99), sum = 401/4  (= 100.25)
#   x1: treated mean = (1+2+3)/3 = 2
#       control weighted mean = (1*3 + 1/4*2 + 99*1)/(401/4) = (205/2)/(401/4)
#                             = 410/401
#       diff_x1 = 2 - 410/401 = 392/401
#   x2: treated mean = 2/3
#       control weighted mean = (0 + 0 + 99)/(401/4) = 396/401
#       diff_x2 = 396/401 - 2/3 = (1188 - 802)/1203 = 386/1203
#   contributions: 2*392/401 = 784/401 = 2352/1203;  4*386/1203 = 1544/1203
#   total_ATT = 3896/1203 = 3.238570241064007
# ==========================================================================

default_deltas <- function() c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)

make_oal_data <- function(n = 300, seed = 20260702) {
  set.seed(seed)
  Xc1 <- rnorm(n); Xc2 <- rnorm(n)
  Xp1 <- rnorm(n); Xp2 <- rnorm(n)
  Xi1 <- rnorm(n); Xi2 <- rnorm(n)
  Xs1 <- rnorm(n); Xs2 <- rnorm(n)
  lp <- 0.6 * Xc1 + 0.6 * Xc2 + 0.9 * Xi1 + 0.9 * Xi2
  A <- rbinom(n, 1L, plogis(lp))
  y <- 2 + A + 0.9 * Xc1 + 0.9 * Xc2 + 0.7 * Xp1 + 0.7 * Xp2 + rnorm(n)
  d <- data.frame(A, y, Xc1, Xc2, Xp1, Xp2, Xi1, Xi2, Xs1, Xs2)
  rownames(d) <- paste0("id", seq_len(n))
  d
}

oal_formula <- function() {
  stats::as.formula("A ~ Xc1 + Xc2 + Xp1 + Xp2 + Xi1 + Xi2 + Xs1 + Xs2")
}

make_wamd_fixture <- function() {
  list(
    treat = c(1, 1, 1, 0, 0, 0),
    covs  = cbind(x1 = c(1, 2, 3, 3, 2, 1),
                  x2 = c(0, 1, 1, 0, 0, 1)),
    ps    = c(0.5, 0.8, 0.001, 0.5, 0.2, 0.998),
    coef  = c(x1 = 2, x2 = -4)
  )
}

# Small dataset with a factor covariate (dummy expansion + predict xlevels).
# seed 20260703, n = 150: realized arms 80 treated / 70 control.
make_factor_data <- function(n = 150, seed = 20260703) {
  set.seed(seed)
  x1 <- rnorm(n)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  A <- rbinom(n, 1L, plogis(0.7 * x1 + 0.5 * (g == "b") - 0.5 * (g == "c")))
  y <- 1 + A + x1 + (g == "b") - (g == "c") + rnorm(n)
  data.frame(A, y, x1, g)
}

# Near-positivity fixture: very strong confounders push the PS onto the
# clip bounds (used for the D.15 diagnostic; the test skips if the fitted
# model does not actually put > 5% of units at a bound).
make_extreme_data <- function(n = 300, seed = 20260704) {
  set.seed(seed)
  Xc1 <- rnorm(n); Xc2 <- rnorm(n); Xs1 <- rnorm(n)
  A <- rbinom(n, 1L, plogis(2.5 * Xc1 + 2.5 * Xc2))
  y <- 1 + A + 2 * Xc1 + 2 * Xc2 + rnorm(n)
  data.frame(A, y, Xc1, Xc2, Xs1)
}

# --------------------------------------------------------------------------
# Reference IPW weights (D.12), re-derived from the paper, base R only.
# --------------------------------------------------------------------------
ref_weights <- function(ps, treat, estimand) {
  if (estimand == "ATT") ifelse(treat == 1, 1, ps / (1 - ps))
  else                   ifelse(treat == 1, 1 / ps, 1 / (1 - ps))
}

# Reference wAMD (D.12): clipped-PS IPW weights; per-column weighted mean
# difference; total = sum over ALL columns of |coef_j| * diff_j.
ref_wamd <- function(ps, treat, covs, coef, estimand = "ATE",
                     clip = c(0.01, 0.99)) {
  e <- pmin(pmax(ps, clip[1]), clip[2])
  w <- ref_weights(e, treat, estimand)
  diffs <- vapply(seq_len(ncol(covs)), function(j) {
    abs(sum(w * covs[, j] * treat) / sum(w * treat) -
        sum(w * covs[, j] * (1 - treat)) / sum(w * (1 - treat)))
  }, numeric(1))
  by <- abs(coef) * diffs
  list(total = sum(by), by.covariate = by)
}

# --------------------------------------------------------------------------
# KKT checker for the TARGET OAL objective on the total-log-likelihood scale
#     f(a0, alpha) = -sum_i [ A_i eta_i - log(1 + exp(eta_i)) ]
#                    + lambda.n * sum_j pen_j |alpha_j|
# with eta = a0 + X alpha.  Subgradient conditions:
#   alpha_j != 0 :  grad_j + lambda.n * pen_j * sign(alpha_j) = 0
#   alpha_j  = 0 :  |grad_j| <= lambda.n * pen_j
#   pen_j = Inf  :  alpha_j must be exactly 0 (glmnet Inf -> exclude)
#   intercept    :  sum_i (plogis(eta_i) - A_i) = 0
# where grad_j = sum_i (plogis(eta_i) - A_i) X_ij is the gradient of the
# NEGATIVE total log-likelihood.  Violations are normalized by
# max(lambda.n * pen_j, 1) so huge instrument penalties do not mask real
# stationarity errors on retained coefficients.
#
# Tolerance calibration (glmnet 4.1.10, default thresh, this fixture,
# verified before the package existed): max relative violation across the
# full S&E grid was 4.3e-4; a wrong lambda-scale constant (e.g. omitting
# the mean(pen) rescale or the 1/n conversion) violates at >= 1e-1.
# --------------------------------------------------------------------------
KKT_TOL <- 2e-3
KKT_INT_TOL <- 1e-6

kkt_max_violation <- function(intercept, alpha, X, A, lambda.n, pen) {
  eta <- drop(intercept + X %*% alpha)
  pr <- stats::plogis(eta)
  grad <- drop(crossprod(X, pr - A))
  viol <- numeric(length(alpha))
  for (j in seq_along(alpha)) {
    if (!is.finite(pen[j])) {           # excluded column: must be exactly 0
      viol[j] <- abs(alpha[j])
      next
    }
    thr <- lambda.n * pen[j]
    viol[j] <- if (alpha[j] != 0) {
      abs(grad[j] + thr * sign(alpha[j])) / max(thr, 1)
    } else {
      max(0, abs(grad[j]) - thr) / max(thr, 1)
    }
  }
  list(coef = max(viol), intercept = abs(sum(pr - A)) / length(A))
}

# Reconstruct the standardized design actually optimized (D.3) from the
# returned object: training-statistics scaling stored in info$center/scale.
std_covs <- function(fit) {
  scale(fit$covs, center = fit$info$center, scale = fit$info$scale)
}
