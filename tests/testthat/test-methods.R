# G.13/G.14 (fragment form) / B.3 -- print/summary output fragments,
# provenance strings (D.19, verbatim), coef(), predict(), plot smoke test,
# and the near-positivity diagnostic (D.15).  Full error messages are
# never matched -- short fragments only; provenance strings are spec'd
# verbatim and matched fixed.

test_that("print shows provenance, retained/excluded lists, rationale and handoff", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  expect_output(print(fit),
                "OAL (Shortreed & Ertefaie 2017, doi:10.1111/biom.12679)",
                fixed = TRUE)
  expect_output(print(fit), "Retained")
  expect_output(print(fit), "Excluded")
  # the FIXED instrument-exclusion one-liner (Brookhart 2006; Myers 2011)
  expect_output(print(fit), "Brookhart")
  expect_output(print(fit), "bias amplification")
  # literal next-call handoff
  expect_output(print(fit), "oal_match", fixed = TRUE)
  expect_output(print(fit), "matchit", fixed = TRUE)
  # provenance is also stored verbatim
  expect_identical(fit$provenance,
                   "OAL (Shortreed & Ertefaie 2017, doi:10.1111/biom.12679)")
})

test_that("summary prints the wAMD selection and balance material", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  s <- summary(fit)
  expect_output(print(s),
                "OAL (Shortreed & Ertefaie 2017, doi:10.1111/biom.12679)",
                fixed = TRUE)
  expect_output(print(s), "wAMD")
})

test_that("binomial and refit provenance suffixes appear (D.19)", {
  d <- make_oal_data()
  d$ybin <- as.integer(d$y > stats::median(d$y))
  fit_b <- suppressWarnings(
    oal(oal_formula(), data = d, outcome = ~ ybin, family = binomial(),
        lambda = c(-1, 0.25))
  )
  expect_match(fit_b$provenance,
               "[binomial outcome model: beyond the validated simulations]",
               fixed = TRUE)
  fit_r <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25),
               refit = TRUE)
  expect_output(print(fit_r), "post-selection refit", fixed = TRUE)
})

test_that("coef() returns both scales", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  expect_equal(coef(fit), fit$coefficients)
  expect_equal(coef(fit, scale = "original"), fit$coefficients)
  expect_equal(coef(fit, scale = "standardized"), fit$coefficients.std)
})

test_that("predict(fit) reproduces $ps without keep.fits, in both refit modes", {
  d <- make_oal_data()
  for (rf in c(FALSE, TRUE)) {
    fit <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25),
               refit = rf, keep.fits = FALSE)
    expect_equal(unname(predict(fit)), unname(fit$ps), tolerance = 1e-12)
    # newdata = the original rows must give the same scores
    expect_equal(unname(predict(fit, newdata = d[1:7, ])),
                 unname(fit$ps[1:7]), tolerance = 1e-10)
  }
})

test_that("predict type = 'link' is the unclipped linear predictor", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  lk <- predict(fit, type = "link")
  expect_equal(unname(pmin(pmax(stats::plogis(lk), fit$clip[1]),
                           fit$clip[2])),
               unname(fit$ps), tolerance = 1e-10)
  # hand-computed link on new rows: original-scale coefficients
  eta <- drop(cbind(1, fit$covs[3:5, , drop = FALSE]) %*%
                fit$coefficients[c("(Intercept)", colnames(fit$covs))])
  expect_equal(unname(predict(fit, newdata = d[3:5, ], type = "link")),
               unname(eta), tolerance = 1e-10)
})

test_that("predict errors on NA in newdata", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y, lambda = c(-1, 0.25))
  nd <- d[1:5, ]
  nd$Xc1[2] <- NA
  expect_error(predict(fit, newdata = nd), "[Mm]issing|NA")
})

test_that("factor covariates: fit works, predict guards factor levels", {
  df <- make_factor_data()
  fit <- oal(A ~ x1 + g, data = df, outcome = ~ y, lambda = c(-1, 0.25))
  expect_length(fit$ps, nrow(df))
  expect_true(all(fit$ps > 0 & fit$ps < 1))
  # dummies are standardized like any column (D.3), so there are 3 or more
  # model-matrix columns
  expect_gte(ncol(fit$covs), 3L)
  nd <- df[1:4, ]
  nd$g <- factor(c("a", "b", "c", "zz"),
                 levels = c("a", "b", "c", "zz"))    # unseen level
  expect_error(predict(fit, newdata = nd), "level|factor")
})

test_that("plot.oal runs for type = 'wamd'", {
  d <- make_oal_data()
  fit <- oal(oal_formula(), data = d, outcome = ~ y)
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({
    grDevices::dev.off()
    unlink(pdf_file)
  })
  expect_no_error(plot(fit, type = "wamd"))
})

test_that("near-positivity diagnostic is printed when > 5% sit at a clip bound", {
  d <- make_extreme_data()
  fit <- oal(A ~ Xc1 + Xc2 + Xs1, data = d, outcome = ~ y)
  at_bound <- mean(fit$ps <= fit$clip[1] | fit$ps >= fit$clip[2])
  skip_if(at_bound <= 0.05,
          "fixture did not trigger the clip-bound diagnostic")
  expect_output(print(fit), "clipping bound")
  # ...and the mild fixture must NOT trigger it
  d2 <- make_oal_data()
  fit2 <- oal(oal_formula(), data = d2, outcome = ~ y)
  if (mean(fit2$ps <= fit2$clip[1] | fit2$ps >= fit2$clip[2]) <= 0.05) {
    out <- paste(utils::capture.output(print(fit2)), collapse = "\n")
    expect_false(grepl("clipping bound", out, fixed = TRUE))
  }
})
