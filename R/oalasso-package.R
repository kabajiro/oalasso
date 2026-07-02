# oalasso-package.R -- package-level documentation and import directives.

#' oalasso: Outcome-Adaptive Lasso Propensity Scores
#'
#' Estimates propensity scores by the outcome-adaptive lasso (OAL) of
#' Shortreed and Ertefaie (2017) and the generalized outcome-adaptive lasso
#' (GOAL) of Balde, Yang and Lefebvre (2023), using \pkg{glmnet} with an exact
#' penalty-scale correction so that the published objectives and tuning grids
#' are reproduced, and tuning by the papers' weighted absolute mean difference
#' (wAMD) balance criterion. The resulting score is designed to be supplied
#' directly to [MatchIt::matchit()] as a distance measure, to
#' [WeightIt::weightit()] as a propensity score, or to `psAve::psave()` as an
#' appended candidate.
#'
#' The single estimation function is [oal()]. Its result hands off to the
#' existing ecosystem: [oal_match()] / [oal_weight()] (or the equivalent
#' explicit [MatchIt::matchit()] / [WeightIt::weightit()] calls) and
#' [cobalt::bal.tab()] (which has a method for `oal` objects). [oal_wamd()]
#' exposes the selection criterion for methods research and testing.
#'
#' ## Design notes
#' \pkg{glmnet} sits in `Imports` -- the one principled departure from the
#' suite's engines-in-Suggests rule -- because it is the sole solver of the
#' estimator's objective, and the exact-scale correction depends on its
#' penalty-factor rescaling, its infinite-penalty-to-exclusion conversion, and
#' its `coef(..., exact = TRUE)` re-supply semantics (floor `>= 4.1-2`). The
#' archived \pkg{lqa} package of the original reference implementation is not
#' a dependency. \pkg{cobalt} powers the display balance table and
#' [bal.tab.oal()]; the wAMD criterion itself is native code.
#'
#' Note the name collision: the GOAL implemented here (`method = "goal"`) is
#' the generalized outcome-adaptive lasso of Balde, Yang and Lefebvre (2023)
#' for binary treatments, not the unrelated "GOAL" generalized propensity
#' score software for continuous exposures by Gao and colleagues.
#'
#' @references
#' Shortreed SM, Ertefaie A (2017). Outcome-adaptive lasso: variable selection
#' for causal inference. *Biometrics*, 73(4), 1111-1122.
#' \doi{10.1111/biom.12679}
#'
#' Balde I, Yang Y, Lefebvre G (2023). Reader reaction to "Outcome-adaptive
#' lasso: variable selection for causal inference" by Shortreed and Ertefaie
#' (2017). *Biometrics*, 79(1), 514-520. \doi{10.1111/biom.13683}
#'
#' Jones B, Ertefaie A, Shortreed SM (2023). Rejoinder to reader reaction "On
#' the use of the outcome-adaptive lasso for propensity score estimation".
#' *Biometrics*, 79(1), 521-525. \doi{10.1111/biom.13681}
#'
#' Balde I (2025). Oracle properties of the generalized outcome-adaptive
#' lasso. *Statistics & Probability Letters*.
#' \doi{10.1016/j.spl.2025.110379}
#'
#' Zou H (2006). The adaptive lasso and its oracle properties. *Journal of the
#' American Statistical Association*, 101(476), 1418-1429.
#' \doi{10.1198/016214506000000735}
#'
#' Brookhart MA, Schneeweiss S, Rothman KJ, Glynn RJ, Avorn J, Sturmer T
#' (2006). Variable selection for propensity score models. *American Journal
#' of Epidemiology*, 163(12), 1149-1156. \doi{10.1093/aje/kwj149}
#'
#' Myers JA, Rassen JA, Gagne JJ, Huybrechts KF, Schneeweiss S, Rothman KJ,
#' Joffe MM, Glynn RJ (2011). Effects of adjusting for instrumental variables
#' on bias and precision of effect estimates. *American Journal of
#' Epidemiology*, 174(11), 1213-1222. \doi{10.1093/aje/kwr364}
#'
#' @importFrom glmnet glmnet
#' @importFrom stats gaussian binomial coef predict weights plogis
#' @importFrom graphics plot
#' @keywords internal
"_PACKAGE"
