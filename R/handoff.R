# handoff.R -- thin, guarded pass-throughs to MatchIt and WeightIt (B.3).
# They add NO arguments of their own and reuse the formula and data stored in
# the oal object, which eliminates the row-misalignment hazard of retyping
# `data =` by construction. The explicit two-step call remains canonical.

#' Match on an outcome-adaptive lasso propensity score
#'
#' Convenience pass-through to [MatchIt::matchit()]: matches on the propensity
#' score of an [oal()] fit, reusing the formula and data stored in the object.
#' Equivalent to the canonical explicit call
#' \preformatted{MatchIt::matchit(<formula>, data = <data>, distance = fit$ps, ...)}
#' but with no opportunity for row misalignment between the two steps. All
#' `...` arguments are forwarded verbatim; the return value is an ordinary
#' `matchit` object, so the full \pkg{MatchIt}/\pkg{cobalt} toolkit applies.
#'
#' @param object An `oal` object.
#' @param ... Arguments forwarded verbatim to [MatchIt::matchit()] (e.g.,
#'   `method`, `caliper`, `ratio`, `replace`).
#'
#' @return A `matchit` object; see [MatchIt::matchit()].
#' @seealso [oal()], [oal_weight()], [MatchIt::matchit()]
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' data("lalonde", package = "MatchIt")
#' fit <- oal(treat ~ age + educ + married + re74, data = lalonde,
#'            outcome = ~ re78)
#' m <- oal_match(fit, method = "nearest", caliper = 0.2)
#' @export
oal_match <- function(object, ...) {
  if (!inherits(object, "oal")) {
    stop("`object` must be an oal object (the result of oal()).", call. = FALSE)
  }
  .require_pkg("MatchIt", "for oal_match()")
  MatchIt::matchit(formula = object$formula, data = object$data,
                   distance = object$ps, ...)
}

#' Weight by an outcome-adaptive lasso propensity score
#'
#' Convenience pass-through to [WeightIt::weightit()]: constructs balancing
#' weights from the propensity score of an [oal()] fit, reusing the stored
#' formula and data. Equivalent to the canonical explicit call
#' \preformatted{WeightIt::weightit(<formula>, data = <data>, ps = fit$ps, estimand = <estimand>, ...)}
#' All `...` arguments are forwarded verbatim; the return value is an ordinary
#' `weightit` object.
#'
#' @param object An `oal` object.
#' @param estimand The estimand passed to [WeightIt::weightit()]; defaults to
#'   the estimand of the fit (note that the wAMD *selection* already used the
#'   fitted estimand's weights).
#' @param ... Arguments forwarded verbatim to [WeightIt::weightit()].
#'
#' @return A `weightit` object; see [WeightIt::weightit()].
#' @seealso [oal()], [oal_match()], [WeightIt::weightit()],
#'   [WeightIt::get_w_from_ps()]
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE) && requireNamespace("WeightIt", quietly = TRUE)
#' data("lalonde", package = "MatchIt")
#' fit <- oal(treat ~ age + educ + married + re74, data = lalonde,
#'            outcome = ~ re78)
#' w <- oal_weight(fit)
#' @export
oal_weight <- function(object, estimand = object$estimand, ...) {
  if (!inherits(object, "oal")) {
    stop("`object` must be an oal object (the result of oal()).", call. = FALSE)
  }
  .require_pkg("WeightIt", "for oal_weight()")
  WeightIt::weightit(formula = object$formula, data = object$data,
                     ps = object$ps, estimand = estimand, ...)
}
