# G.11 -- oal_match / oal_weight round-trips and the cobalt bal.tab method
# (suite pattern, mirroring psAve's psave_match / psave_weight tests).

test_that("oal_match returns a matchit object with distance = fit$ps", {
  skip_if_not_installed("MatchIt")
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  m <- oal_match(fit)
  expect_s3_class(m, "matchit")
  expect_equal(unname(as.numeric(m$distance)), unname(fit$ps),
               tolerance = 1e-12)
})

test_that("oal_match forwards arguments verbatim", {
  skip_if_not_installed("MatchIt")
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  m <- oal_match(fit, method = "nearest", caliper = 0.5)
  expect_s3_class(m, "matchit")
  expect_false(is.null(m$caliper))
  expect_equal(unname(as.numeric(m$distance)), unname(fit$ps),
               tolerance = 1e-12)
})

test_that("oal_weight returns a weightit object with ps = fit$ps", {
  skip_if_not_installed("WeightIt")
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  w <- oal_weight(fit)
  expect_s3_class(w, "weightit")
  expect_equal(unname(as.numeric(w$ps)), unname(fit$ps), tolerance = 1e-12)
  expect_identical(w$estimand, fit$estimand)
  # WeightIt's ATE weights from this ps must equal the stored weights
  expect_equal(unname(as.numeric(w$weights)), unname(fit$weights),
               tolerance = 1e-8)
})

test_that("oal_weight honors an overriding estimand", {
  skip_if_not_installed("WeightIt")
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, estimand = "ATE")
  w <- oal_weight(fit, estimand = "ATT")
  expect_identical(w$estimand, "ATT")
  expect_equal(unname(as.numeric(w$ps)), unname(fit$ps), tolerance = 1e-12)
})

test_that("cobalt::bal.tab dispatches on oal objects", {
  skip_if_not_installed("cobalt")
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  bt <- cobalt::bal.tab(fit)
  expect_s3_class(bt, "bal.tab")
})
