# wamd.R -- the wAMD selection criterion (D.12) and its IPW weight formulas,
# implemented ONCE in the exported oal_wamd(); grid selection (D.13).

# Estimand-specific IPW weights from a (clipped) propensity score (D.12):
#   ATE: w_i = A_i / e_i + (1 - A_i) / (1 - e_i)
#   ATT: w_i = A_i + (1 - A_i) * e_i / (1 - e_i)
.ipw_from_ps <- function(ps, treat, estimand) {
  w <- numeric(length(ps))
  t1 <- treat == 1L
  if (estimand == "ATT") {
    w[t1] <- 1
    w[!t1] <- ps[!t1] / (1 - ps[!t1])
  } else {
    w[t1] <- 1 / ps[t1]
    w[!t1] <- 1 / (1 - ps[!t1])
  }
  w
}

#' Weighted absolute mean difference (wAMD) of Shortreed and Ertefaie (2017)
#'
#' Scores an arbitrary candidate propensity score on the exact wAMD balance
#' criterion used by [oal()] to select its tuning parameters. This function is
#' the package's **single source of truth** for the criterion: the internal
#' grid search calls it verbatim, and exporting it lets any candidate score --
#' from any package -- be scored on the same yardstick (the mirror of
#' `psAve::psave_criteria()`).
#'
#' @details
#' The propensity scores are first clipped to `clip`; the inverse-probability
#' weights at `estimand` are then (with \eqn{e_i} the clipped score)
#' \deqn{ATE:\; w_i = A_i/e_i + (1 - A_i)/(1 - e_i), \qquad
#'       ATT:\; w_i = A_i + (1 - A_i)\, e_i/(1 - e_i).}
#' For each covariate column \eqn{j} the weighted absolute mean difference is
#' \deqn{d_j = \left| \frac{\sum_i w_i X_{ij} A_i}{\sum_i w_i A_i}
#'   - \frac{\sum_i w_i X_{ij} (1 - A_i)}{\sum_i w_i (1 - A_i)} \right|,}
#' and the criterion is \eqn{\mathrm{wAMD} = \sum_j |\beta_j|\, d_j}, summed
#' over **all** columns of `covs` (including any the propensity score model
#' excluded), with \eqn{\beta_j} the outcome-model coefficients supplied in
#' `coef` (their absolute values are taken internally).
#'
#' Inside [oal()], `covs` is the **standardized** model matrix and `coef` the
#' outcome coefficients on that same standardized scale (the wAMD, like the
#' adaptive weights, is scale-dependent). The formula is implemented natively
#' -- never delegated to \pkg{cobalt} -- so that it reproduces the published
#' criterion exactly.
#'
#' @param ps Numeric vector of propensity scores, strictly inside (0, 1)
#'   before clipping.
#' @param treat Treatment vector: numeric 0/1, logical, or a two-level
#'   factor/character (second level = treated), as in [oal()].
#' @param covs Numeric matrix of covariates, one row per unit; no missing
#'   values.
#' @param coef Numeric vector of outcome-model coefficients, one per column of
#'   `covs` (matched by name when both are named, else by position). Absolute
#'   values are used as the balance weights.
#' @param estimand `"ATE"` (default) or `"ATT"`; selects the IPW weight
#'   formula above.
#' @param clip Length-2 numeric; the scores are clipped to
#'   `[clip[1], clip[2]]` **before** the weights are formed (default
#'   `c(0.01, 0.99)`).
#'
#' @return A list with elements
#' \describe{
#'   \item{`total`}{the wAMD, \eqn{\sum_j |\beta_j| d_j};}
#'   \item{`by.covariate`}{named numeric vector of the per-covariate
#'     contributions \eqn{|\beta_j| d_j} (they sum to `total`).}
#' }
#'
#' @references
#' Shortreed SM, Ertefaie A (2017). Outcome-adaptive lasso: variable selection
#' for causal inference. *Biometrics*, 73(4), 1111-1122.
#' \doi{10.1111/biom.12679}
#'
#' @seealso [oal()], [weights.oal()]
#' @examples
#' set.seed(7)
#' X <- matrix(rnorm(200), 100, 2, dimnames = list(NULL, c("x1", "x2")))
#' A <- rbinom(100, 1, plogis(X[, 1]))
#' ps <- plogis(0.8 * X[, 1])
#' oal_wamd(ps, A, X, coef = c(x1 = 0.5, x2 = 0.1), estimand = "ATE")
#' @export
oal_wamd <- function(ps, treat, covs, coef, estimand = c("ATE", "ATT"),
                     clip = c(0.01, 0.99)) {
  estimand <- match.arg(estimand)
  if (!is.numeric(clip) || length(clip) != 2L || anyNA(clip) ||
      clip[1L] <= 0 || clip[2L] >= 1 || clip[1L] >= clip[2L]) {
    stop("`clip` must be a length-2 numeric vector with 0 < clip[1] < clip[2] < 1.",
         call. = FALSE)
  }
  if (!is.matrix(covs) || !is.numeric(covs)) {
    covs <- as.matrix(covs)
    if (!is.numeric(covs)) {
      stop("`covs` must be a numeric matrix (one row per unit).", call. = FALSE)
    }
  }
  n <- nrow(covs)
  if (!is.numeric(ps) || length(ps) != n || anyNA(ps) ||
      any(ps <= 0) || any(ps >= 1)) {
    stop(sprintf("`ps` must be a numeric vector of %d propensity scores strictly inside (0, 1).", n),
         call. = FALSE)
  }
  if (length(treat) != n) {
    stop(sprintf("`treat` must have one value per row of `covs` (%d rows).", n),
         call. = FALSE)
  }
  a <- .coerce_treat(treat, name = "treatment (`treat`)")
  if (anyNA(covs)) {
    stop("`covs` contains missing values; oalasso requires complete cases.",
         call. = FALSE)
  }
  if (!is.numeric(coef) || length(coef) != ncol(covs) || any(!is.finite(coef))) {
    stop(sprintf("`coef` must be a finite numeric vector with one value per column of `covs` (%d columns).",
                 ncol(covs)), call. = FALSE)
  }
  if (!is.null(names(coef)) && !is.null(colnames(covs))) {
    if (!setequal(names(coef), colnames(covs))) {
      stop("The names of `coef` must match colnames(covs).", call. = FALSE)
    }
    coef <- coef[colnames(covs)]
  }

  ## PS clipped BEFORE the IPW weights (D.12)
  e <- .clip_ps(ps, clip)
  w <- .ipw_from_ps(e, a, estimand)
  w1 <- w * (a == 1L)
  w0 <- w * (a == 0L)
  m1 <- as.numeric(crossprod(covs, w1)) / sum(w1)
  m0 <- as.numeric(crossprod(covs, w0)) / sum(w0)
  d <- abs(m1 - m0)
  contrib <- abs(coef) * d
  names(contrib) <- colnames(covs)
  list(total = sum(contrib), by.covariate = contrib)
}
