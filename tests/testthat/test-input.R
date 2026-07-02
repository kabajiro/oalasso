# G.8 / D.1-D.4 -- input validation, treatment coercion, NA policy,
# degenerate designs, argument guards.  Error expectations match SHORT
# fragments only, never full messages.

# a cheap 2-point grid for tests that need a full fit
fast_args <- function() list(lambda = c(-1, 0.25))

## --- treatment coercion (D.1, MatchIt convention) --------------------------

test_that("numeric 0/1 treatment is kept and stored as integer 0/1", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25))
  expect_true(is.integer(fit$treat) || all(fit$treat %in% c(0, 1)))
  expect_identical(as.integer(fit$treat), as.integer(d$A))
})

test_that("two-level factor treatment: second level is treated = 1", {
  d <- make_oal_data()
  d2 <- d
  d2$A <- factor(ifelse(d$A == 1, "trt", "ctrl"), levels = c("ctrl", "trt"))
  fit  <- oal(oal_formula(), data = d,  outcome = ~ y, lambda = c(-1, 0.25))
  fit2 <- oal(oal_formula(), data = d2, outcome = ~ y, lambda = c(-1, 0.25))
  expect_identical(as.integer(fit2$treat), as.integer(d$A))
  expect_equal(fit2$ps, fit$ps, tolerance = 1e-12)
})

test_that("reversed factor levels flip the treated arm", {
  d <- make_oal_data()
  d2 <- d
  d2$A <- factor(ifelse(d$A == 1, "trt", "ctrl"), levels = c("trt", "ctrl"))
  fit2 <- oal(oal_formula(), data = d2, outcome = ~ y, lambda = c(-1, 0.25))
  # second level in factor order is "ctrl" -> treated = original controls
  expect_identical(as.integer(fit2$treat), as.integer(d$A == 0))
})

test_that("logical treatment is coerced to 0/1", {
  d <- make_oal_data()
  d2 <- d
  d2$A <- d$A == 1
  fit2 <- oal(oal_formula(), data = d2, outcome = ~ y, lambda = c(-1, 0.25))
  expect_identical(as.integer(fit2$treat), as.integer(d$A))
})

test_that("character treatment uses factor-level order (second level treated)", {
  d <- make_oal_data()
  d2 <- d
  d2$A <- ifelse(d$A == 1, "b_trt", "a_ctl")   # alphabetical: a_ctl < b_trt
  fit2 <- oal(oal_formula(), data = d2, outcome = ~ y, lambda = c(-1, 0.25))
  expect_identical(as.integer(fit2$treat), as.integer(d$A))
})

test_that("non-0/1 numeric treatment errors", {
  d <- make_oal_data()
  d$A <- d$A + 1                                # 1/2 coding
  expect_error(oal(oal_formula(), data = d, outcome = ~ y), "binary")
})

test_that("more than two treatment levels errors", {
  d <- make_oal_data()
  d$A <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  expect_error(oal(oal_formula(), data = d, outcome = ~ y), "two")
})

## --- NA policy: error, never drop (D.1) ------------------------------------

test_that("NA in the treatment errors, never drops", {
  d <- make_oal_data()
  d$A[5] <- NA
  expect_error(oal(oal_formula(), data = d, outcome = ~ y), "[Mm]issing|NA")
})

test_that("NA in a covariate errors and names the variable", {
  d <- make_oal_data()
  d$Xc1[7] <- NA
  expect_error(oal(oal_formula(), data = d, outcome = ~ y), "Xc1")
})

test_that("NA in the outcome errors and names the variable", {
  d <- make_oal_data()
  d$y[3] <- NA
  expect_error(oal(oal_formula(), data = d, outcome = ~ y), "y")
})

test_that("NA in an unused column is harmless", {
  d <- make_oal_data()
  d$junk <- NA_real_
  fit <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25))
  expect_s3_class(fit, "oal")
})

test_that("fewer than 2 units per arm errors", {
  d <- make_oal_data()
  d$A <- 0
  d$A[1] <- 1
  expect_error(oal(oal_formula(), data = d, outcome = ~ y), "2")
})

## --- outcome specification (D.1) --------------------------------------------

test_that("two-sided outcome formula gives the fixed teaching error", {
  d <- make_oal_data()
  expect_error(oal(oal_formula(), data = d, outcome = y ~ Xc1 + Xc2),
               "supply outcome = ~ y", fixed = TRUE)
})

test_that("one-sided outcome with more than one variable errors", {
  d <- make_oal_data()
  expect_error(oal(oal_formula(), data = d, outcome = ~ y + Xc1), "one")
})

test_that("missing outcome without outcome.coef errors", {
  d <- make_oal_data()
  expect_error(oal(oal_formula(), data = d), "outcome")
})

## --- degenerate designs (D.2, D.4) ------------------------------------------

test_that("zero-variance covariate errors naming the column", {
  d <- make_oal_data()
  d$konst <- 1
  f <- stats::update(oal_formula(), . ~ . + konst)
  expect_error(oal(f, data = d, outcome = ~ y), "konst")
})

test_that("aliased (duplicated) covariate gives the degenerate-beta error", {
  d <- make_oal_data()
  d$Xc1dup <- d$Xc1
  f <- stats::update(oal_formula(), . ~ . + Xc1dup)
  expect_error(oal(f, data = d, outcome = ~ y), "aliased|degenerate|zero")
})

## --- argument guards (B.1, D.18) ---------------------------------------------

test_that("suspected raw penalty values in `lambda` warn about exponents", {
  d <- make_oal_data()
  expect_warning(oal(oal_formula(), data = d, outcome = ~ y,
                     lambda = c(5, 10)), "exponent")
})

test_that("family = binomial() emits the beyond-validation warning", {
  d <- make_oal_data()
  d$ybin <- as.integer(d$y > stats::median(d$y))
  expect_warning(
    fit <- oal(oal_formula(), data = d, outcome = ~ ybin,
               family = binomial(), lambda = c(-1, 0.25)),
    "beyond"
  )
  expect_s3_class(fit, "oal")
})

test_that("unused ... arguments warn", {
  d <- make_oal_data()
  expect_warning(oal(oal_formula(), data = d, outcome = ~ y,
                     lambda = c(-1, 0.25), bogus_argument = 1),
                 "[Uu]nused|[Ii]gnored|reserved")
})

test_that("lambda2 supplied with method = 'oal' errors", {
  d <- make_oal_data()
  expect_error(oal(oal_formula(), data = d, outcome = ~ y,
                   method = "oal", lambda2 = c(0, 0.5)), "lambda2")
})

test_that("outcome-leakage guard runs even when outcome.coef is supplied", {
  d <- make_oal_data()
  cn <- colnames(model.matrix(A ~ Xc1 + Xc2 + y, data = d))[-1]
  b <- stats::setNames(rep(1, length(cn)), cn)
  expect_error(
    oal(A ~ Xc1 + Xc2 + y, data = d, outcome = ~ y, outcome.coef = b),
    "must not appear"
  )
  # without `outcome`, the hook alone cannot detect leakage (documented):
  fit <- oal(A ~ Xc1 + Xc2 + y, data = d, outcome.coef = b)
  expect_s3_class(fit, "oal")
})
