# baltab.R -- display-oriented balance table (D.14) and the bal.tab method
# for oal objects, registered on cobalt's generic so that
# cobalt::bal.tab(fit) "just works" (suite pattern).

# Display-oriented balance table (B.2 `balance`, D.14): computed on the
# ORIGINAL-scale covariates with cobalt native conventions (binary columns
# standardized by sqrt(p(1-p))), absolute values; columns smd.un, smd.wt,
# ks.un, ks.wt -- layout identical to psAve::psave()$balance. The wAMD
# *selection criterion* is native code and never delegated to cobalt (D.12);
# this table is display only.
.balance_table <- function(covs, treat, W, s.d.denom, bin.vars) {
  smd.un <- cobalt::col_w_smd(covs, treat = treat, weights = NULL, std = TRUE,
                              s.d.denom = s.d.denom, abs = TRUE,
                              bin.vars = bin.vars)
  smd.wt <- cobalt::col_w_smd(covs, treat = treat, weights = W, std = TRUE,
                              s.d.denom = s.d.denom, abs = TRUE,
                              bin.vars = bin.vars)
  ks.un <- cobalt::col_w_ks(covs, treat = treat, weights = NULL,
                            bin.vars = bin.vars)
  ks.wt <- cobalt::col_w_ks(covs, treat = treat, weights = W,
                            bin.vars = bin.vars)
  data.frame(smd.un = unname(smd.un), smd.wt = unname(smd.wt),
             ks.un = unname(ks.un), ks.wt = unname(ks.wt),
             row.names = colnames(covs))
}

#' Balance tables for oal objects
#'
#' A method for [cobalt::bal.tab()]: assesses balance on the (original-scale)
#' covariates of an [oal()] fit under the inverse-probability weights implied
#' by the outcome-adaptive lasso propensity score, which is supplied as a
#' `distance` measure.
#'
#' The call delegates to the default \pkg{cobalt} machinery as
#' `cobalt::bal.tab(<covariates>, treat = x$treat, weights = x$weights,
#' s.d.denom = <by estimand>, distance = data.frame(ps = x$ps), ...)`, so all
#' the usual \pkg{cobalt} arguments (`un`, `stats`, `thresholds`, ...) are
#' available and display conventions are \pkg{cobalt}'s own. The *selection
#' criterion* inside [oal()] is the papers' wAMD instead, computed natively on
#' the standardized covariates; see [oal_wamd()].
#'
#' @param x An `oal` object.
#' @param ... Further arguments passed on to [cobalt::bal.tab()] (e.g.,
#'   `un = TRUE`, `thresholds = c(m = 0.1)`).
#'
#' @return A `bal.tab` object; see [cobalt::bal.tab()].
#'
#' @references
#' Shortreed SM, Ertefaie A (2017). Outcome-adaptive lasso: variable selection
#' for causal inference. *Biometrics*, 73(4), 1111-1122.
#' \doi{10.1111/biom.12679}
#'
#' @seealso [oal()], [cobalt::bal.tab()], [plot.oal()]
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' data("lalonde", package = "MatchIt")
#' fit <- oal(treat ~ age + educ + married + re74, data = lalonde,
#'            outcome = ~ re78)
#' cobalt::bal.tab(fit, un = TRUE)
#' @exportS3Method cobalt::bal.tab
bal.tab.oal <- function(x, ...) {
  cobalt::bal.tab(as.data.frame(x$covs),
                  treat = x$treat,
                  weights = x$weights,
                  s.d.denom = .sd_denom(x$estimand),
                  distance = data.frame(ps = as.numeric(x$ps)),
                  ...)
}
