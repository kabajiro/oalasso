# ==========================================================================
# In-test reference implementation of the OAL-rejoinder glmnet recipe
# (jmiahjones/OAL-rejoinder-2022, R/oal_funs.R + R/simulate_oal.R pattern,
# transcribed from the implementation dossier -- nothing downloaded).
#
# Steps, verbatim from the recipe:
#   1. X <- model.matrix(formula)[, -1];  Xs <- scale(X, TRUE, TRUE)
#      (scale ONCE with training statistics, ALL columns incl. dummies).
#   2. betaXY <- coef(lm(Y ~ A + Xs))[covariates]  (full sample, both arms;
#      treatment and intercept coefficients discarded).
#   3. S&E grids: delta in c(-10,-5,-2,-1,-0.75,-0.5,-0.25,0.25,0.49),
#      lambda_n = n^delta paired with gamma = 2*(gcf - delta + 1), gcf = 2.
#   4. Per grid row: pen = |betaXY|^(-gamma);
#      glmnet(x = Xs, y = A, family = "binomial", alpha = 1,
#             standardize = FALSE, intercept = TRUE, penalty.factor = pen);
#      evaluate the TARGET objective  -ll + lambda_n * sum(pen_j |alpha_j|)
#      exactly at  s = mean(pen) * lambda_n / n  via
#      coef(fit, s = s, exact = TRUE, x = Xs, y = A, penalty.factor = pen)
#      (the KKT-verified lambda-scale correction: undoes glmnet's
#      sum-to-nvars penalty-factor rescaling and the total -> mean
#      log-likelihood conversion).
#   5. PS = plogis(intercept + Xs alpha), clipped to [0.01, 0.99]; ATE IPW
#      weights; wAMD = sum_j |betaXY_j| * |weighted mean difference of
#      standardized column j|; select the first minimum over the grid.
# ==========================================================================

rejoinder_oal <- function(formula, data, outcome,
                          lambda_vec = c(-10, -5, -2, -1, -0.75,
                                         -0.5, -0.25, 0.25, 0.49),
                          gcf = 2, clip = c(0.01, 0.99)) {
  stopifnot(requireNamespace("glmnet", quietly = TRUE))
  X <- stats::model.matrix(formula, data)[, -1, drop = FALSE]
  A <- data[[deparse(formula[[2L]])]]
  y <- data[[outcome]]
  n <- nrow(X)
  p <- ncol(X)

  Xs <- scale(X, TRUE, TRUE)
  lm.Y <- stats::lm(Y ~ ., data = data.frame(Y = y, A = A,
                                             as.data.frame(Xs)))
  betaXY <- stats::coef(lm.Y)[colnames(Xs)]

  gamma_vals <- 2 * (gcf - lambda_vec + 1)

  G <- length(lambda_vec)
  wamd_vec <- numeric(G)
  ps_list <- vector("list", G)
  set_list <- vector("list", G)

  for (r in seq_len(G)) {
    lam_n <- n^lambda_vec[r]
    gam <- gamma_vals[r]
    pen <- abs(betaXY)^(-gam)
    fit <- glmnet::glmnet(x = Xs, y = A, family = "binomial", alpha = 1,
                          standardize = FALSE, intercept = TRUE,
                          penalty.factor = pen)
    s <- mean(pen) * lam_n / n
    cf <- as.numeric(stats::coef(fit, s = s, exact = TRUE, x = Xs, y = A,
                                 penalty.factor = pen))
    e <- stats::plogis(drop(cbind(1, Xs) %*% cf))
    e <- pmin(pmax(e, clip[1]), clip[2])
    w <- A / e + (1 - A) / (1 - e)          # ATE weights (S&E default)
    diffs <- vapply(seq_len(p), function(j) {
      abs(sum(w * Xs[, j] * A) / sum(w * A) -
          sum(w * Xs[, j] * (1 - A)) / sum(w * (1 - A)))
    }, numeric(1))
    wamd_vec[r] <- sum(abs(betaXY) * diffs)
    ps_list[[r]] <- e
    set_list[[r]] <- colnames(Xs)[cf[-1] != 0]
  }

  best <- which(wamd_vec <= min(wamd_vec) * (1 + 1e-9))[1L]
  list(delta = lambda_vec[best],
       gamma = gamma_vals[best],
       lambda.n = n^lambda_vec[best],
       wamd = wamd_vec,
       criterion.value = wamd_vec[best],
       ps = ps_list[[best]],
       selected = set_list[[best]],
       sets = set_list,
       betaXY = betaXY)
}
