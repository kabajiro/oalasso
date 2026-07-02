# methods.R -- S3 methods for oal objects: print, summary, print.summary,
# coef, predict, weights (B.3), plus the provenance rendering.

# Render the retained/excluded teaching lists for print.oal.
.render_selection <- function(sel, digits) {
  fmt1 <- function(rows) {
    if (!nrow(rows)) return("  (none)")
    paste0("  ", rows$term, "  (|b| = ",
           format(abs(rows$outcome.coef), digits = digits), ")",
           collapse = "\n")
  }
  ret <- sel[sel$selected, , drop = FALSE]
  exc <- sel[!sel$selected, , drop = FALSE]
  list(retained = fmt1(ret), excluded = fmt1(exc), n.excluded = nrow(exc))
}

# The FIXED instrument-exclusion one-liner (B.3, exact wording).
.instrument_note <- function() {
  paste0("Excluded covariates showed no outcome association; excluding instruments ",
         "and noise variables from a propensity score model improves precision and ",
         "avoids bias amplification (Brookhart et al. 2006; Myers et al. 2011).")
}

#' Print an oal object
#'
#' Prints a one-screen teaching summary of a fitted [oal()] object: the
#' provenance label first (which method, which tuning provenance), the sample,
#' the selected tuning point and its wAMD, the retained
#' (outcome-related) versus excluded (outcome-unrelated) covariates with their
#' outcome-coefficient magnitudes, the fixed instrument-exclusion rationale,
#' and then the **literal next call** -- echoing the formula and data name
#' from your own [oal()] call -- that hands the score to
#' [MatchIt::matchit()] or [WeightIt::weightit()]. A near-positivity
#' diagnostic is appended when more than 5\% of the propensity scores sit at a
#' clipping bound.
#'
#' @param x An `oal` object.
#' @param digits Number of significant digits to print. Default 3.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @seealso [oal()], [summary.oal()]
#' @export
print.oal <- function(x, digits = 3, ...) {
  sub <- substitute(x)
  obj <- if (is.name(sub)) as.character(sub) else "fit"

  cat(x$provenance, "\n", sep = "")
  cat("An oal object (outcome-adaptive lasso propensity score)\n")
  cat(sprintf(" - estimand: %s;  refit: %s\n", x$estimand,
              if (isTRUE(x$refit)) "TRUE (post-selection glm)" else "FALSE (penalized fit)"))
  cat(sprintf(" - sample:   %d units (%d treated, %d control); %d covariate column(s)\n",
              x$info$n, x$info$n.treated, x$info$n - x$info$n.treated, x$info$p))
  l <- x$lambda
  cat(sprintf(" - selected: delta = %s (lambda_n = %s), gamma = %s%s\n",
              format(l[["delta"]], digits = digits),
              format(l[["lambda.n"]], digits = digits),
              format(l[["gamma"]], digits = digits),
              if (!is.na(l[["lambda2"]]))
                sprintf(", lambda2 = %s", format(l[["lambda2"]], digits = digits))
              else ""))
  cat(sprintf(" - wAMD at the selection: %s\n",
              format(x$criterion.value, digits = digits)))

  r <- .render_selection(x$selected, digits)
  cat("\nRetained (outcome-related):\n", r$retained, "\n", sep = "")
  cat("Excluded (outcome-unrelated):\n", r$excluded, "\n", sep = "")
  if (r$n.excluded > 0L) {
    cat("\n", .instrument_note(), "\n", sep = "")
  }

  ## the literal next call, echoing the user's own formula/data symbols
  f.txt <- tryCatch(deparse1(x$call$formula), error = function(e) NULL)
  if (is.null(f.txt) || !nzchar(f.txt) || f.txt == "NULL") f.txt <- deparse1(x$formula)
  d.txt <- tryCatch(deparse1(x$call$data), error = function(e) "data")
  cat("\nNext:\n")
  cat(sprintf("  MatchIt::matchit(%s, data = %s, distance = %s$ps)\n",
              f.txt, d.txt, obj))
  cat(sprintf("    or: oal_match(%s)\n", obj))
  cat(sprintf("  WeightIt::weightit(%s, data = %s, ps = %s$ps, estimand = \"%s\")\n",
              f.txt, d.txt, obj, x$estimand))
  cat(sprintf("    or: oal_weight(%s)\n", obj))
  cat(sprintf("  ($ps is numeric, named by rownames(%s), strictly inside (0, 1);\n", d.txt))
  cat(sprintf("   it also satisfies psAve::psave(ps.append = cbind(oal = %s$ps)).)\n", obj))

  share <- .clip_share(x$ps, x$clip)
  if (share > 5) {
    cat("\n", .clip_note(share), "\n", sep = "")
  }
  invisible(x)
}

#' Summarize an oal object
#'
#' Produces (a) the full per-covariate selection table (`selected`), (b) a
#' path summary (the best grid row per `lambda2` value, when the path was
#' kept), (c) the outcome (weight) model coefficient table, (d) the full
#' balance table (unweighted vs. weighted SMD and KS, with a `*` marker at
#' weighted SMD > 0.1), and (e) the near-positivity clip diagnostic.
#'
#' @param object An `oal` object.
#' @param ... Ignored.
#' @param x A `summary.oal` object.
#' @param digits Number of significant digits to print. Default 3.
#'
#' @return For `summary.oal()`, an object of class `"summary.oal"`: a list
#'   with elements `call`, `provenance`, `method`, `estimand`, `refit`,
#'   `lambda`, `criterion`, `criterion.value`, `selected`, `path.summary`,
#'   `outcome.model`, `balance`, `clip`, `clip.share`, and `nn`.
#'   `print.summary.oal()` returns `x` invisibly.
#' @seealso [oal()], [print.oal()], [bal.tab.oal()]
#' @export
summary.oal <- function(object, ...) {
  path.summary <- NULL
  if (!is.null(object$path)) {
    sp <- split(object$path, object$path$lambda2)
    if (!length(sp)) sp <- list(object$path)   # lambda2 all NA (method = "oal")
    best <- do.call(rbind, lapply(sp, function(d) d[.first_min(d$wamd), , drop = FALSE]))
    rownames(best) <- NULL
    path.summary <- best
  }
  out <- list(call = object$call,
              provenance = object$provenance,
              method = object$method,
              estimand = object$estimand,
              refit = object$refit,
              lambda = object$lambda,
              criterion = object$criterion,
              criterion.value = object$criterion.value,
              selected = object$selected,
              path.summary = path.summary,
              outcome.model = object$outcome.model,
              balance = object$balance,
              clip = object$clip,
              clip.share = .clip_share(object$ps, object$clip),
              nn = c(control = sum(object$treat == 0L),
                     treated = sum(object$treat == 1L)))
  class(out) <- "summary.oal"
  out
}

#' @rdname summary.oal
#' @export
print.summary.oal <- function(x, digits = 3, ...) {
  cat(x$provenance, "\n", sep = "")
  cat("Summary of an oal fit\n")
  cat(sprintf("Call: %s\n\n", deparse1(x$call)))
  cat(sprintf("Estimand: %s;  criterion: %s (weighted absolute mean difference);  refit: %s\n",
              x$estimand, x$criterion, x$refit))
  cat(sprintf("Sample: %d treated, %d control\n",
              x$nn[["treated"]], x$nn[["control"]]))
  l <- x$lambda
  cat(sprintf("Selected: delta = %s (lambda_n = %s), gamma = %s%s;  wAMD = %s\n\n",
              format(l[["delta"]], digits = digits),
              format(l[["lambda.n"]], digits = digits),
              format(l[["gamma"]], digits = digits),
              if (!is.na(l[["lambda2"]]))
                sprintf(", lambda2 = %s", format(l[["lambda2"]], digits = digits))
              else "",
              format(x$criterion.value, digits = digits)))

  cat("Selection table (outcome coefficients on the standardized scale):\n")
  sel <- x$selected
  sel$outcome.coef <- round(sel$outcome.coef, digits)
  sel$penalty.factor <- signif(sel$penalty.factor, digits)
  sel$coef <- round(sel$coef, digits)
  print(sel)
  if (any(!sel$selected)) {
    cat("\n", .instrument_note(), "\n", sep = "")
  }

  if (!is.null(x$path.summary)) {
    cat("\nBest grid row per lambda2 (by wAMD):\n")
    ps <- x$path.summary
    for (nm in c("lambda.n", "s", "wamd")) ps[[nm]] <- signif(ps[[nm]], digits)
    print(ps)
  }

  cat(sprintf("\nOutcome (weight) model: %s\n", x$outcome.model$label))
  om <- x$outcome.model$coefficients
  for (nm in intersect(c("estimate", "se"), names(om))) {
    om[[nm]] <- signif(om[[nm]], digits)
  }
  print(om)

  cat("\nBalance (original-scale covariates):\n")
  b <- round(x$balance, digits)
  b$` ` <- ifelse(x$balance$smd.wt > 0.1, "*", "")
  print(b)
  cat("---\n'*' = weighted SMD > 0.1\n")

  if (x$clip.share > 5) {
    cat("\n", .clip_note(x$clip.share), "\n", sep = "")
  }
  invisible(x)
}

#' Extract propensity score model coefficients
#'
#' Returns the selected outcome-adaptive lasso propensity score coefficients
#' (intercept first), either on the original covariate scale (default; the
#' scale used by [predict.oal()]) or on the standardized scale that the
#' penalized objective was actually optimized on (the scale of the adaptive
#' weights and of the wAMD -- methods-paper material).
#'
#' @param object An `oal` object.
#' @param scale `"original"` (default) or `"standardized"`.
#' @param ... Ignored.
#'
#' @return A named numeric vector of length p + 1.
#' @seealso [oal()], [predict.oal()]
#' @export
coef.oal <- function(object, scale = c("original", "standardized"), ...) {
  scale <- match.arg(scale)
  if (scale == "original") object$coefficients else object$coefficients.std
}

#' Predict outcome-adaptive lasso propensity scores for new data
#'
#' Computes propensity scores (or the linear predictor) for new observations
#' from the stored terms, factor levels, and ORIGINAL-scale coefficients of an
#' [oal()] fit. Works without `keep.fits = TRUE` and for both `refit` modes;
#' the training standardization is never re-estimated (the original-scale
#' coefficients already absorb it). Missing values in `newdata` are an error.
#'
#' @param object An `oal` object.
#' @param newdata A data frame containing the variables of the propensity
#'   score formula, or `NULL` (default) for the in-sample values.
#' @param type `"ps"` (default): propensity scores clipped to the fit's
#'   `clip`; `"link"`: the unclipped linear predictor.
#' @param ... Ignored.
#'
#' @return A numeric vector with one value per row, named by rownames. For
#'   `newdata = NULL` and `type = "ps"` this is exactly `object$ps`.
#' @seealso [oal()], [coef.oal()]
#' @export
predict.oal <- function(object, newdata = NULL, type = c("ps", "link"), ...) {
  type <- match.arg(type)
  cf <- object$coefficients
  if (is.null(newdata)) {
    if (type == "ps") return(object$ps)
    eta <- cf[1L] + as.numeric(object$covs %*% cf[-1L])
    names(eta) <- rownames(object$data)
    return(eta)
  }
  newdata <- as.data.frame(newdata)
  des <- .build_design(object$terms, newdata, xlev = object$xlevels,
                       contrasts = object$info$contrasts, what = "`newdata`")
  X <- des$X
  if (!setequal(colnames(X), names(cf)[-1L])) {
    stop(paste0("The model matrix built from `newdata` does not reproduce the training columns; ",
                "supply the same variables (and factor levels) as at fit time."),
         call. = FALSE)
  }
  X <- X[, names(cf)[-1L], drop = FALSE]
  eta <- cf[1L] + as.numeric(X %*% cf[-1L])
  names(eta) <- rownames(newdata)
  if (type == "link") return(eta)
  .clip_ps(stats::plogis(eta), object$clip)
}

#' Extract the inverse-probability weights of an oal fit
#'
#' Returns `object$weights`: the IPW weights at the *fitted* estimand implied
#' by the clipped propensity score (ATE:
#' \eqn{A/e + (1-A)/(1-e)}; ATT: \eqn{A + (1-A)\,e/(1-e)}). For other
#' estimands use `WeightIt::get_w_from_ps(object$ps, object$treat,
#' estimand = ...)` or [oal_weight()].
#'
#' @param object An `oal` object.
#' @param ... Ignored.
#'
#' @return The numeric vector `object$weights`, named by `rownames(data)`.
#' @seealso [oal()], [oal_weight()], [oal_wamd()]
#' @export
weights.oal <- function(object, ...) {
  object$weights
}
