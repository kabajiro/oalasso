# G.6 / D.10 -- GOAL: exact nesting at lambda2 = 0, grid layout, and
# provenance labeling.  The correctness of the augmentation itself
# (Zou-Hastie trick + (1 + lambda2) rescale) is established by the
# augmented-objective KKT test in test-kkt.R.

test_that("method = 'goal' with lambda2 = 0 reproduces method = 'oal' exactly", {
  d <- make_oal_data()
  fit_oal  <- oal(oal_formula(), data = d, outcome = ~ y, method = "oal")
  fit_goal <- oal(oal_formula(), data = d, outcome = ~ y, method = "goal",
                  lambda2 = 0)
  # lambda2 = 0 grid points MUST run the plain-OAL code path (no
  # augmentation), so these are equal to machine precision, not merely
  # statistically close.
  expect_equal(fit_goal$ps, fit_oal$ps, tolerance = 1e-12)
  expect_equal(fit_goal$coefficients, fit_oal$coefficients,
               tolerance = 1e-12)
  expect_equal(fit_goal$coefficients.std, fit_oal$coefficients.std,
               tolerance = 1e-12)
  expect_equal(fit_goal$weights, fit_oal$weights, tolerance = 1e-12)
  expect_equal(fit_goal$criterion.value, fit_oal$criterion.value,
               tolerance = 1e-12)
  expect_equal(fit_goal$path$wamd, fit_oal$path$wamd, tolerance = 1e-12)
  expect_equal(unname(fit_goal$lambda[["delta"]]),
               unname(fit_oal$lambda[["delta"]]))
  expect_equal(unname(fit_goal$lambda[["gamma"]]),
               unname(fit_oal$lambda[["gamma"]]))
  # the selected lambda2 is reported as 0 (goal) vs NA (oal)
  expect_equal(unname(fit_goal$lambda[["lambda2"]]), 0)
  expect_true(is.na(fit_oal$lambda[["lambda2"]]))
})

test_that("GOAL default lambda2 grid: layout, order and selection bookkeeping", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, method = "goal")
  deltas <- default_deltas()
  # the author's published grid (Balde 2025 supplement, "taken from Zou
  # and Hastie (2005)"): 11 lambda2 values -> 11 x 9 = 99 grid rows
  l2 <- c(0, 10^c(-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1))
  expect_length(l2, 11L)
  # grid order defines tie-breaking (D.7): lambda2 ascending (0 first)
  # OUTER x delta ascending INNER
  expect_identical(nrow(fit$path), length(l2) * length(deltas))
  expect_identical(nrow(fit$path), 99L)
  expect_equal(fit$path$lambda2, rep(l2, each = length(deltas)))
  expect_equal(fit$path$delta, rep(deltas, times = length(l2)))
  expect_identical(sum(fit$path$selected), 1L)
  expect_true(unname(fit$lambda[["lambda2"]]) %in% l2)
  row <- fit$path[fit$path$selected, , drop = FALSE]
  expect_equal(row$lambda2, unname(fit$lambda[["lambda2"]]))
  expect_equal(row$wamd, fit$criterion.value)
  # lambda.n exponent base (Balde 2025 supplement): n^delta on the
  # lambda2 = 0 (plain OAL) rows, (n + q)^delta with q = p augmentation
  # rows on the augmented rows
  n <- nrow(d)
  p <- ncol(fit$covs)
  z <- fit$path$lambda2 == 0
  expect_equal(fit$path$lambda.n[z], n^fit$path$delta[z])
  expect_equal(fit$path$lambda.n[!z], (n + p)^fit$path$delta[!z])
})

test_that("a positive lambda2 changes the fit relative to plain OAL", {
  d <- make_oal_data()
  fit0 <- oal(oal_formula(), data = d, outcome = ~ y, lambda = 0.25,
              gamma = 2)
  fit2 <- oal(oal_formula(), data = d, outcome = ~ y, method = "goal",
              lambda = 0.25, gamma = 2, lambda2 = 10)
  expect_false(isTRUE(all.equal(unname(fit0$coefficients.std),
                                unname(fit2$coefficients.std),
                                tolerance = 1e-6)))
})

test_that("GOAL provenance: author's published grid vs user-specified lambda2 grid", {
  d <- make_oal_data()
  fit_def <- oal(oal_formula(), data = d, outcome = ~ y, method = "goal",
                 lambda = c(-1, 0.25))
  expect_match(fit_def$provenance,
               "GOAL (Balde, Yang & Lefebvre 2023, doi:10.1111/biom.13683)",
               fixed = TRUE)
  expect_match(fit_def$provenance,
               "author's published grid (Balde 2025 supplement)",
               fixed = TRUE)
  expect_output(print(fit_def), "author's published grid", fixed = TRUE)

  fit_usr <- oal(oal_formula(), data = d, outcome = ~ y, method = "goal",
                 lambda = c(-1, 0.25), lambda2 = c(0, 0.5))
  expect_match(fit_usr$provenance, "user-specified lambda2 grid",
               fixed = TRUE)
})

test_that("selection matches Balde's nested lambda2 rule, including on ties", {
  d <- make_oal_data()

  # Balde's nested rule (2025 supplement): per lambda2, which.min of wAMD
  # over the (delta, gamma) grid; then which.min over the per-lambda2
  # minima (first minimum at both levels).  The package applies the same
  # rule with the documented 1e-9 relative tie tolerance at both levels
  # (it reduces to Balde's exact which.min under exact ties); the same
  # tolerance is used here so floating-point ulp noise cannot flip the
  # comparison.
  nested_pick <- function(path, tol = 1e-9) {
    fm <- function(x) which(x <= min(x) + tol * max(1, abs(min(x))))[1L]
    blocks <- unique(path$lambda2)
    per <- vapply(blocks,
                  function(l2) min(path$wamd[path$lambda2 == l2]),
                  numeric(1))
    l2.star <- blocks[fm(per)]
    rows <- which(path$lambda2 == l2.star)
    rows[fm(path$wamd[rows])]
  }

  # (a) generic fixture: flat first-min == nested rule
  fit <- oal(oal_formula(), data = d, outcome = ~ y, method = "goal",
             lambda2 = c(0, 0.5, 10))
  expect_identical(which(fit$path$selected), nested_pick(fit$path))

  # (b) constructed EXACT tie across the whole grid: near-zero user
  # outcome coefficients make every penalty explode, so no covariate is
  # ever selected; the PS is constant at every grid point and the wAMD
  # (weighted mean differences with constant per-arm weights) is the same
  # number everywhere.  Both the flat rule and Balde's nested rule must
  # return the FIRST grid row: smallest lambda2 (= 0), then smallest
  # delta.
  bv <- c(Xc1 = 1e-6, Xc2 = 1e-6, Xp1 = 1e-6, Xp2 = 1e-6,
          Xi1 = 1e-6, Xi2 = 1e-6, Xs1 = 1e-6, Xs2 = 1e-6)
  fit_tie <- oal(oal_formula(), data = d, outcome.coef = bv,
                 method = "goal", lambda = c(0.25, 0.49), gamma = 2,
                 lambda2 = c(0, 0.5))
  skip_if(diff(range(fit_tie$path$wamd)) >
            1e-9 * max(1, abs(min(fit_tie$path$wamd))),
          "fixture did not produce a grid-wide tie on this platform")
  expect_identical(which(fit_tie$path$selected), 1L)
  expect_identical(which(fit_tie$path$selected), nested_pick(fit_tie$path))
  expect_equal(unname(fit_tie$lambda[["lambda2"]]), 0)
  expect_equal(unname(fit_tie$lambda[["delta"]]), 0.25)
})
