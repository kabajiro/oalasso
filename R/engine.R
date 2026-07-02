# engine.R -- the estimation engine (D.4-D.11): weight (outcome) model,
# penalty factors, tuning grid, glmnet paths with the exact lambda-scale
# correction, GOAL augmentation, and the optional post-selection refit.
# Nothing here is exported.

# ---------------------------------------------------------------------------
# Weight (outcome) model (D.4): lm(y ~ A + Xs) for gaussian() or
# glm(y ~ A + Xs, binomial()) on the FULL sample, both arms pooled (the
# S&E/rejoinder convention). `b` = covariate coefficients only; the treatment
# and intercept coefficients are discarded and the treatment is never
# penalized anywhere.
#
# Degenerate-beta policy (internal model): any NA or exactly-zero b_j is an
# ERROR naming the columns. OLS betas are almost surely nonzero, so an exact
# zero/NA signals a data pathology (rank deficiency / aliasing) that would
# silently corrupt the penalty-scale correction.
.fit_outcome_model <- function(y, treat, Xs, family) {
  XA <- cbind("(Intercept)" = 1, A = as.numeric(treat), Xs)
  fit <- suppressWarnings(stats::glm.fit(x = XA, y = y, family = family))
  cf <- fit$coefficients
  b <- cf[-(1:2)]
  names(b) <- colnames(Xs)
  bad <- names(b)[!is.finite(b) | b == 0]
  if (length(bad)) {
    stop(sprintf(paste0("Degenerate outcome-model coefficient(s) (NA or exactly zero) for column(s) %s: ",
                        "rank-deficient or aliased outcome model; remove collinear columns or supply ",
                        "`outcome.coef`."),
                 paste0('"', bad, '"', collapse = ", ")), call. = FALSE)
  }

  ## compact coefficient table (estimate + SE from the final IRLS/OLS qr)
  se <- rep(NA_real_, length(cf))
  ok <- tryCatch({
    Rm <- qr.R(fit$qr)
    piv <- fit$qr$pivot
    dispersion <- if (identical(family$family, "gaussian")) {
      sum(fit$residuals^2 * fit$weights) / fit$df.residual
    } else 1
    se[piv] <- sqrt(diag(chol2inv(Rm)) * dispersion)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) se <- rep(NA_real_, length(cf))
  tab <- data.frame(term = c("(Intercept)", "A", colnames(Xs)),
                    estimate = as.numeric(cf),
                    se = se,
                    row.names = NULL)

  label <- if (identical(family$family, "binomial")) "glm(y ~ A + X, binomial)" else "lm(y ~ A + X)"
  list(b = b, table = tab, label = label, fit = fit)
}

# ---------------------------------------------------------------------------
# Tuning grid (D.7). `deltas` are exponents (lambda_n = n^delta); gamma is
# either PAIRED, gamma_r = 2 * (gamma.factor - delta_r + 1) -- the
# S&E/rejoinder/fOAL convention, derivable from lambda_n * n^(gamma/2 - 1)
# = n^gcf -- or a FIXED scalar crossed with all deltas (Schnitzer-style).
# fAL's `2*(gcf - delta) + 1` is a transcription bug and is never implemented.
# Grid order (defines tie-breaking, D.13): lambda2 ascending (0 first) OUTER
# x delta ascending INNER.
#
# lambda.n exponent base (Balde 2025 supplement, verified 2026-07-02): plain
# OAL rows (lambda2 = 0) use lambda.n = n^delta; GOAL-augmented rows
# (lambda2 > 0) use lambda.n = (n + q)^delta with q = p augmentation rows --
# Balde's GOAL loop calls adaptive.lasso(lambda = n.q^(il), ...) with
# n.q = n + q.  lambda2 = 0 grid points deliberately stay on the plain-OAL
# base n^delta so that GOAL nests plain OAL exactly (D.10).
.make_grid <- function(deltas, lambda2s, gamma, gamma.factor, n, p) {
  g <- expand.grid(delta = deltas, lambda2 = lambda2s,
                   KEEP.OUT.ATTRS = FALSE)   # delta varies fastest = inner
  g <- g[, c("lambda2", "delta")]
  g$lambda.n <- ifelse(g$lambda2 > 0, (n + p)^g$delta, n^g$delta)
  g$gamma <- if (is.null(gamma)) 2 * (gamma.factor - g$delta + 1) else gamma
  rownames(g) <- NULL
  g
}

# ---------------------------------------------------------------------------
# The exact lambda-scale correction (D.9, THE core correctness requirement).
#
# Target objective (total log-likelihood scale):
#     -l(alpha) + lambda_n * sum_j pen_j * |alpha_j|,  pen_j = |b_j|^(-gamma).
# glmnet minimizes  -(1/n_eff) l(alpha) + s * sum_j pf~_j |alpha_j|  where the
# penalty factors are internally rescaled to sum to nvars. Recipe per the
# Shortreed/Ertefaie rejoinder code (jmiahjones/OAL-rejoinder-2022,
# R/oal_funs.R + tests/test-glmnet.R, KKT-verified): evaluate at
#     s = mean(pen) * lambda_n / n_eff
# via coef(fit, s = s, exact = TRUE, x = , y = , penalty.factor = ), which
# undoes the sum-to-nvars rescaling and converts total -> per-observation
# loss. n_eff = n for plain OAL; with GOAL augmentation 1/n is replaced by
# d[1]/sum(d) over the augmented sample (unit weights -> 1/(n + p)).
#
# GOAL raw penalty constant (Balde 2025 supplement, verified 2026-07-02):
# on augmented grid points (lambda2 > 0) the target objective is the
# AUGMENTED-sample binomial likelihood penalized by
#     lambda.n * sum_j pen_j |alpha_j|  with  lambda.n = (n + q)^delta,
# q = p augmentation rows (Balde: adaptive.lasso(lambda = n.q^(il)),
# n.q = n + q), so glmnet needs s = mean(pen) * (n + q)^delta / (n + q):
# the SAME .penalty_s() formula with lambda.n already on the (n + q) base
# (computed in .make_grid) and n.eff = n + q. Plain OAL keeps n^delta / n.
# The KKT test on the augmented objective pins this normalization.
#
# Infinite penalty factors (user `outcome.coef` zeros): glmnet converts
# pen_j = Inf to `exclude` and internally RESETS those factors to 1 BEFORE the
# sum-to-nvars rescale, so the KKT-correct constant generalizes to
#     s = (sum(pen[finite]) + #Inf) / nvars * lambda_n / n_eff,
# which reduces exactly to mean(pen) * lambda_n / n_eff when all factors are
# finite (KKT-verified against glmnet 4.1.10; the KKT test in the package test
# suite is the acceptance criterion).
.penalty_s <- function(pen, lambda.n, n.eff) {
  fin <- is.finite(pen)
  ((sum(pen[fin]) + sum(!fin)) / length(pen)) * lambda.n / n.eff
}

# ---------------------------------------------------------------------------
# Unpenalized post-selection refit (D.11, Schnitzer variant): logistic
# glm.fit on the selected columns of the SAME standardized matrix. Returns
# full-length standardized coefficients (zeros for unselected columns), the
# convergence flag, and any glm.fit warning messages (collected so oal() can
# emit ONE consolidated warning; no silent fallback).
.refit_glm <- function(Xs, A, selected) {
  p <- ncol(Xs)
  xm <- cbind(1, Xs[, selected, drop = FALSE])
  warn <- character(0)
  rf <- withCallingHandlers(
    stats::glm.fit(x = xm, y = A, family = stats::binomial()),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  cf <- rf$coefficients
  if (anyNA(cf)) {
    bad <- colnames(Xs)[selected][which(is.na(cf[-1]))]
    stop(sprintf(paste0("The post-selection refit produced NA coefficient(s) for column(s) %s ",
                        "(collinear selected columns). Remove collinear columns or use refit = FALSE."),
                 paste0('"', bad, '"', collapse = ", ")), call. = FALSE)
  }
  out <- numeric(p + 1L)
  out[1L] <- cf[1L]
  out[1L + which(selected)] <- cf[-1L]
  list(cf = out, converged = isTRUE(rf$converged), warnings = warn, fit = rf)
}

# ---------------------------------------------------------------------------
# Back-transform standardized-scale coefficients (intercept first) to the
# original covariate scale (D.3):
#   alpha_j    = alpha_std_j / scale_j
#   intercept  = intercept_std - sum(alpha_std * center / scale).
.back_transform <- function(cf.std, center, scale) {
  a <- cf.std[-1L] / scale
  c(cf.std[1L] - sum(cf.std[-1L] * center / scale), a)
}

# ---------------------------------------------------------------------------
# The grid engine (D.8-D.12): one glmnet() path per distinct
# (lambda2, gamma) -- in paired mode one per grid row, in fixed/user-gamma
# mode ONE path reused across the lambda grid -- each grid point evaluated
# EXACTLY at its corrected s. Returns per-grid-point results plus the
# standardized coefficient matrix.
.run_grid <- function(Xs, A, b, grid, refit, clip, estimand, keep.fits,
                      verbose) {
  n <- nrow(Xs)
  p <- ncol(Xs)
  G <- nrow(grid)
  A <- as.numeric(A)

  cs <- matrix(0, nrow = p + 1L, ncol = G,
               dimnames = list(c("(Intercept)", colnames(Xs)), NULL))
  s.used <- numeric(G)
  wamd <- numeric(G)
  n.selected <- integer(G)
  refit.conv <- rep(NA, G)
  sets <- vector("list", G)
  ps.mat <- matrix(NA_real_, nrow = n, ncol = G)
  fits <- if (keep.fits) list() else NULL
  refit.warn <- character(0)

  ## lambda2 = 0 grid points MUST run the plain-OAL code path: zero augmented
  ## rows would still contribute intercept terms to the binomial likelihood
  ## and break exact nesting (D.10).
  for (l2 in unique(grid$lambda2)) {
    rows.l2 <- which(grid$lambda2 == l2)
    if (l2 > 0) {
      aug <- sqrt(l2) * diag(p)
      colnames(aug) <- colnames(Xs)
      xx <- rbind(Xs, aug)              # Zou-Hastie augmentation (D.10)
      yy <- c(A, rep(0, p))
      n.eff <- n + p
    } else {
      xx <- Xs
      yy <- A
      n.eff <- n
    }
    for (gam in unique(grid$gamma[rows.l2])) {
      rows <- rows.l2[grid$gamma[rows.l2] == gam]
      pen <- abs(b)^(-gam)              # penalty factors (D.6)
      path <- glmnet::glmnet(x = xx, y = yy, family = "binomial", alpha = 1,
                             standardize = FALSE, intercept = TRUE,
                             penalty.factor = pen)
      if (keep.fits) {
        fits[[sprintf("lambda2=%s,gamma=%s", format(l2), format(gam))]] <- path
      }
      for (r in rows) {
        s <- .penalty_s(pen, grid$lambda.n[r], n.eff)
        s.used[r] <- s
        ## exact-objective evaluation (D.9): resupply x, y, penalty.factor.
        ## glmnet's lambda.interp() emits an uninformative "collapsing to
        ## unique 'x' values" warning when a degenerate path has a single
        ## lambda; with exact = TRUE the solution is re-solved at s (never
        ## interpolated), so only that specific warning is muffled.
        cf <- withCallingHandlers(
          as.numeric(stats::coef(path, s = s, exact = TRUE,
                                 x = xx, y = yy, penalty.factor = pen)),
          warning = function(w) {
            if (grepl("collapsing to unique 'x' values",
                      conditionMessage(w), fixed = TRUE)) {
              invokeRestart("muffleWarning")
            }
          })
        ## elastic-net rescale (D.10): ALL coefficients INCLUDING the
        ## intercept (Balde 2025 supplement, verified 2026-07-02:
        ## PS = expit(cbind(1, X) %*% (1 + lambda2) * coef))
        if (l2 > 0) cf <- (1 + l2) * cf
        sel <- cf[-1L] != 0

        if (refit) {
          rf <- .refit_glm(Xs, A, sel)
          cf.use <- rf$cf
          refit.conv[r] <- rf$converged
          refit.warn <- c(refit.warn, rf$warnings)
        } else {
          cf.use <- cf
        }

        cs[, r] <- cf.use
        sets[[r]] <- colnames(Xs)[sel]
        n.selected[r] <- sum(sel)
        ## per-grid-point PS from the ORIGINAL n rows only, clipped BEFORE
        ## the wAMD IPW weights (D.11)
        eta <- cf.use[1L] + as.numeric(Xs %*% cf.use[-1L])
        ps <- .clip_ps(stats::plogis(eta), clip)
        ps.mat[, r] <- ps
        ## wAMD on the SAME standardized matrix with |b| weights (D.12);
        ## oal_wamd() is the single source of truth for the criterion
        wamd[r] <- oal_wamd(ps, A, Xs, b, estimand = estimand, clip = clip)$total
        .vmsg(verbose,
              sprintf("grid %d/%d: lambda2 = %s, delta = %s, gamma = %s -> %d selected, wAMD = %s",
                      r, G, .fmt(l2), .fmt(grid$delta[r]), .fmt(gam),
                      n.selected[r], .fmt(wamd[r])))
      }
    }
  }

  if (refit) {
    bad <- which(!refit.conv)
    if (length(bad) || length(refit.warn)) {
      msg <- character(0)
      if (length(bad)) {
        msg <- sprintf("the unpenalized refit did not converge at %d of %d grid point(s) (see path$refit.converged)",
                       length(bad), G)
      }
      if (length(refit.warn)) {
        msg <- c(msg, sprintf("glm.fit reported: %s",
                              paste(unique(refit.warn), collapse = "; ")))
      }
      warning(sprintf("refit = TRUE: %s. No silent fallback is applied; inspect the path before use.",
                      paste(msg, collapse = "; ")), call. = FALSE)
    }
  }

  list(coef.std = cs, s = s.used, wamd = wamd, n.selected = n.selected,
       refit.converged = refit.conv, sets = sets, ps = ps.mat, fits = fits)
}
