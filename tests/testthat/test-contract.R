# G.11 / G.7 / B.2 -- the returned-object contract, including the exact
# psAve `ps.append` contract (numeric, length n, named by rownames(data),
# strictly inside (0, 1)), internal consistency of the two coefficient
# scales, path bookkeeping, and tie-breaking.

test_that("$ps satisfies the psAve/MatchIt/WeightIt handoff contract", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  ps <- fit$ps
  expect_true(is.numeric(ps))
  expect_length(ps, nrow(d))
  expect_identical(names(ps), rownames(d))
  expect_false(anyNA(ps))
  # clipped to the default [0.01, 0.99] (= psave's clip, so re-clipping
  # is a no-op) and therefore strictly inside (0, 1)
  expect_true(all(ps >= 0.01 & ps <= 0.99))
  expect_true(all(ps > 0 & ps < 1))
})

test_that("a custom clip is respected", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25),
             clip = c(0.05, 0.95))
  expect_true(all(fit$ps >= 0.05 & fit$ps <= 0.95))
  expect_equal(fit$clip, c(0.05, 0.95))
})

test_that("$ps round-trips into psAve::psave(ps.append = )", {
  skip_if_not_installed("psAve")
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  ma <- psAve::psave(oal_formula(), data = d, outcome = ~ y,
                     ps.methods = "glm", prog.methods = "glm",
                     ps.append = fit$ps)
  expect_s3_class(ma, "psave")
})

test_that("both coefficient scales are stored and mutually consistent (D.3)", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  p <- ncol(fit$covs)
  expect_length(fit$coefficients, p + 1L)
  expect_length(fit$coefficients.std, p + 1L)
  expect_identical(names(fit$coefficients)[1], "(Intercept)")

  # back-transform (D.3): alpha_j = alpha_std_j / scale_j;
  # intercept = intercept_std - sum(alpha_std * center / scale)
  ctr <- fit$info$center[colnames(fit$covs)]
  scl <- fit$info$scale[colnames(fit$covs)]
  a_std <- fit$coefficients.std[-1][colnames(fit$covs)]
  expect_equal(unname(fit$coefficients[-1][colnames(fit$covs)]),
               unname(a_std / scl), tolerance = 1e-10)
  expect_equal(unname(fit$coefficients[1]),
               unname(fit$coefficients.std[1] - sum(a_std * ctr / scl)),
               tolerance = 1e-10)

  # the original-scale coefficients regenerate $ps
  eta <- drop(cbind(1, fit$covs) %*%
                fit$coefficients[c("(Intercept)", colnames(fit$covs))])
  expect_equal(unname(fit$ps),
               unname(pmin(pmax(stats::plogis(eta), fit$clip[1]),
                           fit$clip[2])),
               tolerance = 1e-10)
})

test_that("path/sets/coef.path bookkeeping is coherent", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  G <- length(default_deltas())
  expect_s3_class(fit$path, "data.frame")
  expect_true(all(c("lambda2", "delta", "lambda.n", "gamma", "s", "wamd",
                    "n.selected", "refit.converged", "selected")
                  %in% names(fit$path)))
  expect_identical(nrow(fit$path), G)
  expect_identical(sum(fit$path$selected), 1L)
  expect_length(fit$sets, G)
  expect_identical(dim(fit$coef.path), c(ncol(fit$covs) + 1L, G))
  # sets vs coef.path agree
  for (r in seq_len(G)) {
    nz <- rownames(fit$coef.path)[-1][fit$coef.path[-1, r] != 0]
    expect_identical(sort(nz), sort(fit$sets[[r]]))
  }
})

test_that("keep.path = FALSE drops path, sets and coef.path", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25),
             keep.path = FALSE)
  expect_null(fit$path)
  expect_null(fit$sets)
  expect_null(fit$coef.path)
  # the deliverable is unaffected
  expect_length(fit$ps, nrow(d))
})

test_that("selection is the documented first-minimum rule (1e-9 relative tolerance)", {
  d <- make_oal_data()
  for (fit in list(oal(oal_formula(), data = d, outcome = ~ y),
                   oal(oal_formula(), data = d, outcome = ~ y,
                       method = "goal", lambda2 = c(0, 0.5, 10)))) {
    w <- fit$path$wamd
    rule <- which(w <= min(w) * (1 + 1e-9))[1L]
    expect_identical(which(fit$path$selected), rule)
  }
})

test_that("ties prefer the first grid row (smallest delta)", {
  d <- make_oal_data()
  # with a FIXED gamma, delta = -10 and -5 give lambda_n = n^-10 and n^-5:
  # both penalties are numerically zero, the two solutions coincide, and
  # the wAMDs tie to ~1e-11 relative (probed with glmnet 4.1.10); the
  # documented rule must return the FIRST row, delta = -10.
  fit <- oal(oal_formula(), data = d, outcome = ~ y, gamma = 2,
             lambda = c(-10, -5))
  w <- fit$path$wamd
  skip_if(abs(w[2] - w[1]) > 1e-9 * abs(w[1]),
          "fixture did not produce a tie on this platform")
  expect_equal(unname(fit$lambda[["delta"]]), -10)
})

test_that("info and scalar bookkeeping fields are populated", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  expect_identical(as.numeric(fit$info$n), as.numeric(nrow(d)))
  expect_identical(as.numeric(fit$info$p), as.numeric(ncol(fit$covs)))
  expect_identical(as.numeric(fit$info$n.treated), as.numeric(sum(d$A == 1)))
  expect_identical(fit$estimand, "ATE")
  expect_identical(fit$method, "oal")
  expect_false(fit$refit)
  expect_identical(fit$criterion, "wamd")
  expect_identical(names(fit$lambda),
                   c("delta", "lambda.n", "gamma", "lambda2"))
  expect_equal(unname(fit$lambda[["lambda.n"]]),
               nrow(d)^unname(fit$lambda[["delta"]]))
  expect_identical(as.character(fit$info$glmnet.version),
                   as.character(utils::packageVersion("glmnet")))
  expect_true(is.call(fit$call))
})

test_that("the balance table has the psAve layout", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  expect_s3_class(fit$balance, "data.frame")
  expect_true(all(c("smd.un", "smd.wt", "ks.un", "ks.wt")
                  %in% names(fit$balance)))
  expect_identical(nrow(fit$balance), ncol(fit$covs))
})

test_that("keep.fits controls the $fits slot", {
  d <- make_oal_data()
  fit0 <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25))
  expect_null(fit0$fits)
  fit1 <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25),
              keep.fits = TRUE)
  expect_true(is.list(fit1$fits))
})
