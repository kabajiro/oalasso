# oal.R -- the main (and only) estimation function: validation and
# orchestration only. The heavy lifting lives in input.R, engine.R, wamd.R.

#' Outcome-adaptive lasso propensity scores
#'
#' `oal()` estimates a propensity score by the outcome-adaptive lasso (OAL) of
#' Shortreed and Ertefaie (2017) or its generalized elastic-net form (GOAL) of
#' Balde, Yang and Lefebvre (2023): an adaptive lasso on the treatment
#' (propensity score) log-likelihood whose penalty weights come from an
#' *outcome* regression, so that covariates unrelated to the outcome --
#' instruments and noise variables -- are excluded from the propensity score
#' model. Tuning parameters are selected by the papers' weighted absolute mean
#' difference (wAMD) balance criterion. The result is deliberately modest: a
#' numeric score vector designed to be handed to [MatchIt::matchit()] as
#' `distance`, to [WeightIt::weightit()] as `ps`, or to `psAve::psave()` as an
#' appended candidate.
#'
#' @param formula A two-sided formula `treat ~ x1 + x2 + ...`, exactly as in
#'   [MatchIt::matchit()]. The right-hand side defines the propensity score
#'   covariates; the same covariates (plus the treatment) form the outcome
#'   model that supplies the adaptive penalty weights.
#' @param data A data frame containing the variables in `formula` and
#'   `outcome`. Complete cases in all used variables are REQUIRED; any missing
#'   value is an error naming the offending variables, never a silent row
#'   drop.
#' @param outcome A one-sided formula `~ y` naming the outcome variable
#'   (exactly one variable). Required unless `outcome.coef` is supplied, in
#'   which case it is unused. A two-sided formula is an error: OAL derives its
#'   penalty weights from an outcome model on the *same* covariates as the
#'   propensity score model, so there is no separate outcome design to
#'   specify.
#' @param method `"oal"` (default): the outcome-adaptive lasso; `"goal"`: the
#'   generalized outcome-adaptive lasso, which adds an elastic-net
#'   \eqn{\lambda_2 \sum_j \alpha_j^2} term (grouping effect and numerical
#'   stability under correlated covariates and near-positivity violations).
#' @param estimand `"ATE"` (default) or `"ATT"`; determines the
#'   inverse-probability weights used inside the wAMD criterion and returned
#'   in `weights`. The ATE default is Shortreed and Ertefaie's wAMD weighting
#'   and deliberately differs from `psAve::psave()`'s ATT default.
#' @param family The outcome (weight) model family: `gaussian()` (default;
#'   the specification validated by the papers' simulations) or `binomial()`.
#'   `binomial()` emits a one-time warning: binomial outcome models are beyond
#'   the simulations validated by Shortreed and Ertefaie (2017); GLM theory
#'   per Balde (2025).
#' @param lambda Numeric vector of EXPONENTS \eqn{\delta}; the penalties
#'   actually applied are \eqn{\lambda_n = n^{\delta}}. The default is
#'   Shortreed and Ertefaie's grid
#'   `c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)`. A guard warns when
#'   `max(lambda) > 3`: such values look like raw penalty values, not
#'   exponents (convert a raw grid via `delta = log(lambda_n)/log(n)`).
#' @param gamma `NULL` (default) pairs \eqn{\gamma} with each \eqn{\delta} as
#'   \eqn{\gamma = 2(\code{gamma.factor} - \delta + 1)} (the
#'   Shortreed-Ertefaie/rejoinder convention); a single number (e.g. `2.5`)
#'   fixes \eqn{\gamma} and crosses it with the whole `lambda` grid (the
#'   Schnitzer-style fixed-\eqn{\gamma} mode). No vector form is accepted.
#' @param gamma.factor The convergence factor `gcf` in the pairing formula
#'   (default `2`); ignored when `gamma` is given.
#' @param lambda2 Numeric grid of elastic-net ridge constants for
#'   `method = "goal"` only (supplying it otherwise is an error). Default
#'   `c(0, 10^c(-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1))` -- the
#'   author's published grid (Balde 2025 supplement, "taken from Zou and
#'   Hastie (2005)"), verified against the official GOAL code on 2026-07-02;
#'   `0` is always a member and nests plain OAL exactly. `lambda2` is
#'   selected jointly with `(lambda, gamma)` by the wAMD.
#' @param outcome.coef The extension hook: a NAMED numeric vector of outcome
#'   coefficients \eqn{b_j} over the **standardized** model-matrix columns
#'   (names must exactly cover the columns), overriding the internal outcome
#'   model. **The values are interpreted on the standardized scale of the
#'   internal design matrix** -- adaptive weights are scale-dependent. All
#'   values must be finite and not all zero. Zeros ARE allowed and mean a
#'   hard drop: \eqn{|0|^{-\gamma} = \infty} and glmnet natively converts an
#'   infinite penalty factor to an exclusion. No epsilon or cap is ever
#'   applied. Using this hook flips the provenance label to
#'   "screening use only" (see Details) and sets
#'   `info$gamma.mode = "user-coef"`. If `outcome` is also supplied, it is
#'   still validated and the outcome-leakage guard still runs (the outcome
#'   variable may not appear among the PS covariates); if `outcome` is
#'   omitted, the package **cannot** detect outcome leakage in `formula` --
#'   ensuring the coefficients were derived without post-treatment or outcome
#'   information is then entirely the caller's responsibility.
#' @param refit `FALSE` (default): the propensity score comes from the
#'   penalized fit itself (the Shortreed-Ertefaie/rejoinder convention).
#'   `TRUE`: at every grid point an *unpenalized* logistic regression is refit
#'   on the selected covariates, and both the propensity score and the wAMD
#'   come from the refit (the Schnitzer et al. 2025 variant, adapted from a
#'   longitudinal setting); per-grid-point convergence is recorded in
#'   `path$refit.converged` and any non-convergence triggers one warning --
#'   never a silent fallback.
#' @param clip Length-2 numeric: propensity scores are clipped to
#'   `[clip[1], clip[2]]` BEFORE the wAMD IPW weights are formed and in the
#'   returned `ps` (default `c(0.01, 0.99)`, equal to `psAve::psave()`'s
#'   default so that psave's re-clipping is a no-op). The default is a
#'   Shortreed-Ertefaie-lineage safety choice; Balde's official reference
#'   code runs UNCLIPPED weights (`1/e`, `1/(1 - e)` with no truncation).
#'   To reproduce that behavior, effectively disable clipping with
#'   `clip = c(1e-12, 1 - 1e-12)`.
#' @param keep.path If `TRUE` (default), the full per-grid-point results are
#'   stored in `path`, `sets`, and `coef.path`.
#' @param keep.fits If `TRUE`, the fitted \pkg{glmnet} paths, the outcome
#'   model, and (with `refit = TRUE`) the selected refit are retained in
#'   `fits`. Default `FALSE`.
#' @param verbose If `TRUE`, progress messages report the grid evaluation and
#'   the selected point.
#' @param ... Reserved for future use; supplying unused arguments triggers a
#'   warning.
#'
#' @details
#' ## Objective and exact glmnet evaluation
#' OAL solves, on the TOTAL log-likelihood scale,
#' \deqn{\hat\alpha = \arg\min_\alpha\; -\ell(\alpha; A, X_s)
#'   + \lambda_n \sum_j |b_j|^{-\gamma} |\alpha_j|,}
#' where \eqn{\ell} is the binomial log-likelihood of treatment on the
#' standardized covariates and \eqn{b_j} are the outcome-model coefficients.
#' \pkg{glmnet} minimizes the *per-observation* loss and internally rescales
#' penalty factors to sum to the number of variables, so each grid point is
#' evaluated **exactly** at
#' \deqn{s = \overline{\mathrm{pen}} \cdot \lambda_n / n, \qquad
#'   \mathrm{pen}_j = |b_j|^{-\gamma},}
#' via `coef(fit, s = s, exact = TRUE, x = , y = , penalty.factor = )`, which
#' re-solves at that penalty (the KKT-verified recipe of the
#' Shortreed/Ertefaie rejoinder code). When `outcome.coef` contains zeros,
#' some \eqn{\mathrm{pen}_j = \infty}: glmnet converts those to exclusions and
#' internally resets their factor to 1 *before* the sum-to-nvars rescale, so
#' the constant generalizes to
#' \eqn{s = (\sum_{j: \mathrm{finite}} \mathrm{pen}_j + \#\{\mathrm{pen}_j =
#' \infty\})/p \cdot \lambda_n / n} -- identical to
#' \eqn{\overline{\mathrm{pen}}\,\lambda_n/n} when all factors are finite.
#' This correction is always on; it is what makes the \eqn{n^\delta} grid
#' carry its published meaning.
#'
#' ## Standardization protocol (fixed, not an argument)
#' The model matrix `X <- model.matrix(formula, data)[, -1]` (default
#' treatment contrasts) is standardized ONCE, `Xs <- scale(X, TRUE, TRUE)`,
#' all columns including dummies (the rejoinder convention). The outcome model
#' AND glmnet are both fit on this SAME `Xs` with
#' `glmnet(standardize = FALSE)`. Per-subset scaling is structurally
#' impossible. Coefficients are back-transformed for output:
#' \eqn{\alpha_j = \alpha^{std}_j / s_j} and
#' \eqn{\alpha_0 = \alpha^{std}_0 - \sum_j \alpha^{std}_j c_j / s_j}, with
#' centers \eqn{c} and scales \eqn{s} stored in `info`.
#'
#' ## Weight (outcome) model and the degenerate-beta policy
#' By default \eqn{b} comes from `lm(y ~ A + Xs)` (or
#' `glm(y ~ A + Xs, binomial)`) on the FULL sample, both arms pooled -- the
#' Shortreed-Ertefaie/rejoinder convention (not control-arm-only). Only the
#' covariate coefficients are used; the treatment coefficient is discarded and
#' the treatment is never penalized anywhere. Any `NA` or *exactly zero*
#' internal coefficient is an ERROR naming the columns: OLS betas are almost
#' surely nonzero, so an exact zero or `NA` signals rank deficiency or
#' aliasing that would silently corrupt the penalty-scale correction. Remove
#' the collinear columns or supply `outcome.coef`.
#'
#' ## Tuning grid and the \eqn{\gamma} pairing
#' `lambda` supplies exponents \eqn{\delta} with
#' \eqn{\lambda_n = n^{\delta}}, `n = nrow(model.matrix)` always. With
#' `gamma = NULL` each \eqn{\delta} is paired with
#' \eqn{\gamma = 2(\mathrm{gcf} - \delta + 1)} (default
#' \eqn{\mathrm{gcf} = 2}), the formula derivable from
#' \eqn{\lambda_n\, n^{\gamma/2 - 1} = n^{\mathrm{gcf}}}; there are no
#' \eqn{\lambda \times \gamma} cross-products in paired mode. A scalar
#' `gamma` is crossed with all \eqn{\delta} (Schnitzer's \eqn{\gamma = 2.5}
#' is reachable this way; convert their raw \eqn{\lambda} grid via
#' \eqn{\delta = \log \lambda_n / \log n}).
#'
#' ## GOAL
#' For \eqn{\lambda_2 > 0} the elastic-net term is solved by the Zou-Hastie
#' augmentation: \eqn{X_{aug} = \mathrm{rbind}(X_s, \sqrt{\lambda_2} I_p)}
#' with \eqn{p} pseudo-responses 0, followed by rescaling ALL coefficients
#' INCLUDING the intercept by \eqn{(1 + \lambda_2)} (Balde 2025 supplement:
#' \eqn{PS = \mathrm{expit}(\mathrm{cbind}(1, X)\,(1 + \lambda_2)\,\hat\alpha)});
#' the propensity score uses the original \eqn{n} rows only. On augmented
#' grid points the raw penalty constant is
#' \eqn{\lambda_n = (n + q)^{\delta}} with \eqn{q = p} augmentation rows
#' (Balde's `adaptive.lasso(lambda = n.q^(il))`, `n.q = n + q`), and the
#' penalty-scale constant uses \eqn{1/(n + q)} in place of \eqn{1/n}, i.e.
#' \eqn{s = \overline{\mathrm{pen}} (n + q)^{\delta} / (n + q)}.
#' \eqn{\lambda_2 = 0} grid points run the plain-OAL code path (no
#' augmentation, \eqn{\lambda_n = n^{\delta}}), so `method = "goal"` with
#' `lambda2 = 0` nests `method = "oal"` exactly. (Balde's script literally
#' augments with \eqn{q} zero rows and the \eqn{(n+q)^\delta} base even at
#' \eqn{\lambda_2 = 0}; the zero rows still shift the intercept, so exact
#' nesting is deliberately preferred here.) \eqn{\lambda_2} is selected
#' jointly with \eqn{(\delta, \gamma)} by the wAMD: the flat first-minimum
#' over the \eqn{\lambda_2}-outer \eqn{\times} \eqn{\delta}-inner grid,
#' which is equivalent to Balde's nested rule (per-\eqn{\lambda_2} minima
#' over \eqn{(\delta, \gamma)}, then the first minimum over
#' \eqn{\lambda_2}).
#'
#' ## wAMD criterion and selection
#' Every grid point is scored by [oal_wamd()]: propensity scores clipped to
#' `clip` first (note that Balde's reference implementation is unclipped;
#' see the `clip` argument), IPW weights at `estimand`, per-column weighted
#' absolute mean differences on the standardized matrix summed with weights
#' \eqn{|b_j|} over ALL columns (including excluded ones). The selected grid point is the
#' FIRST minimum (relative tolerance `1e-9`) in the grid order
#' \eqn{\lambda_2} ascending (0 first) outer \eqn{\times} \eqn{\delta}
#' ascending inner -- ties therefore prefer plain OAL and then the smallest
#' penalty. The wAMD is a balance heuristic: no published proof guarantees it
#' lands \eqn{\lambda_n} inside the window required by the theory
#' (\eqn{\lambda_n/\sqrt{n} \to 0},
#' \eqn{\lambda_n n^{\gamma/2 - 1} \to \infty}, \eqn{\gamma > 1}), and it
#' inherits the outcome model's misspecification vulnerability.
#'
#' ## Why exclude instruments?
#' Covariates that predict treatment but not outcome should be excluded from
#' a propensity score model: they increase the variance of the effect estimate
#' without decreasing bias (Brookhart et al. 2006) and amplify the bias from
#' any unmeasured confounding (Myers et al. 2011). OAL's penalty implements
#' exactly this exclusion.
#'
#' ## Provenance labels
#' `print()` and `summary()` always print a provenance line first:
#' `"OAL (Shortreed & Ertefaie 2017, doi:10.1111/biom.12679)"` or
#' `"GOAL (Balde, Yang & Lefebvre 2023, doi:10.1111/biom.13683)"` with a
#' clause stating whether the `lambda2` grid is the author's published grid
#' (Balde 2025 supplement) or user-specified. Supplying `outcome.coef` appends
#' `"user-supplied outcome coefficients — screening use only; no oracle
#' property"`; `family = binomial()` and `refit = TRUE` append their own
#' flags.
#'
#' ## Handoff contract
#' `$ps` is a numeric vector of length `n`, named by `rownames(data)`,
#' strictly inside (0, 1) after clipping -- it satisfies
#' `psAve::psave(ps.append = )`, `MatchIt::matchit(distance = )`, and
#' `WeightIt::weightit(ps = )` by construction.
#'
#' There is no `seed` argument: the whole pipeline is deterministic.
#'
#' @return An object of class `"oal"`: a list with components
#' \describe{
#'   \item{`ps`}{numeric(n), named by `rownames(data)`: the propensity score,
#'     clipped to `clip`, strictly in (0, 1) -- **the deliverable**.}
#'   \item{`weights`}{numeric(n): the IPW weights at `estimand` implied by
#'     `ps`.}
#'   \item{`coefficients`}{named numeric(p+1): selected propensity score
#'     coefficients on the ORIGINAL covariate scale (intercept first).}
#'   \item{`coefficients.std`}{the same on the standardized scale actually
#'     optimized.}
#'   \item{`selected`}{data frame, one row per model-matrix column: `term`,
#'     `outcome.coef` (signed \eqn{b_j}, standardized scale),
#'     `penalty.factor`, `coef` (original scale), `selected` (logical), and
#'     `role` (`"retained"`/`"excluded"`).}
#'   \item{`lambda`}{named numeric: `delta`, `lambda.n` (\eqn{= n^\delta};
#'     \eqn{= (n + p)^\delta} on GOAL grid points with \eqn{\lambda_2 > 0}),
#'     `gamma`, `lambda2` (`NA` for `method = "oal"`) at the selected point.}
#'   \item{`criterion`, `criterion.value`}{`"wamd"` and its value at the
#'     selection.}
#'   \item{`path`}{data frame (or `NULL` without `keep.path`), one row per
#'     grid point in grid order: `lambda2`, `delta`, `lambda.n`, `gamma`, `s`
#'     (the exact glmnet `s` used), `wamd`, `n.selected`, `refit.converged`
#'     (`NA` unless `refit`), `selected`.}
#'   \item{`sets`}{list (or `NULL`): the selected variable set at every grid
#'     point.}
#'   \item{`coef.path`}{(p+1) x G matrix (or `NULL`): propensity score
#'     coefficients across the grid, original scale.}
#'   \item{`outcome.coef`}{named numeric(p): the \eqn{b_j} actually used
#'     (internal or user), standardized scale.}
#'   \item{`outcome.model`}{compact list: coefficient table, `family`, and a
#'     `label` (`"lm(y ~ A + X)"`, `"glm(y ~ A + X, binomial)"`, or
#'     `"user-supplied"`); the full fit object only under `keep.fits`.}
#'   \item{`balance`}{data frame per ORIGINAL covariate column: `smd.un`,
#'     `smd.wt`, `ks.un`, `ks.wt` via \pkg{cobalt} -- identical layout to
#'     `psAve::psave()$balance`.}
#'   \item{`provenance`}{the fixed provenance label, printed first by
#'     `print()`/`summary()`.}
#'   \item{`treat`}{integer(n) 0/1 treatment as coerced.}
#'   \item{`covs`}{n x p numeric original-scale model matrix with
#'     `attr(, "bin.vars")`.}
#'   \item{`estimand`, `method`, `refit`, `clip`}{as resolved.}
#'   \item{`formula`, `data`, `terms`, `xlevels`}{stored to power
#'     [oal_match()], [oal_weight()] and [predict.oal()].}
#'   \item{`fits`}{list (or `NULL`): `glmnet` (one path per
#'     \eqn{(\lambda_2, \gamma)}), `outcome`, `refit` -- iff `keep.fits`.}
#'   \item{`info`}{list: `n`, `n.treated`, `p`, `grid.size`, `family`,
#'     `gamma.mode` (`"paired"`/`"fixed"`/`"user-coef"`), `center`, `scale`,
#'     `contrasts`, `glmnet.version`, `oalasso.version`.}
#'   \item{`call`}{the matched call.}
#' }
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
#' Schnitzer ME, Talbot D, Liu Y, Berger C, Wang G, O'Loughlin J, Sylvestre
#' M-P, Ertefaie A (2025). Outcome-adaptive LASSO for longitudinal data.
#' *Statistics in Medicine* (arXiv:2410.08283).
#'
#' @seealso [oal_match()], [oal_weight()], [oal_wamd()], [predict.oal()],
#'   [plot.oal()], [bal.tab.oal()]
#' @examplesIf requireNamespace("MatchIt", quietly = TRUE)
#' data("lalonde", package = "MatchIt")
#'
#' ## 1) Matching
#' fit <- oal(treat ~ age + educ + race + married + nodegree + re74 + re75,
#'            data = lalonde, outcome = ~ re78)
#' fit
#' m <- oal_match(fit)   # = MatchIt::matchit(f, lalonde, distance = fit$ps)
#'
#' ## 2) Weighting (requires WeightIt)
#' # w <- oal_weight(fit)  # = WeightIt::weightit(f, lalonde, ps = fit$ps, "ATE")
#'
#' ## 3) psAve composition (requires psAve); cbind() labels the candidate
#' ##    while keeping the row names aligned with rownames(data)
#' # ma <- psAve::psave(treat ~ age + educ + race + married + nodegree +
#' #                      re74 + re75,
#' #                    data = lalonde, outcome = ~ re78,
#' #                    ps.append = cbind(oal = fit$ps))
#' @export
oal <- function(formula,
                data,
                outcome,
                method       = c("oal", "goal"),
                estimand     = c("ATE", "ATT"),
                family       = gaussian(),
                lambda       = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
                gamma        = NULL,
                gamma.factor = 2,
                lambda2      = c(0, 10^c(-2, -1.5, -1, -0.75, -0.5,
                                         -0.25, 0, 0.25, 0.5, 1)),
                outcome.coef = NULL,
                refit        = FALSE,
                clip         = c(0.01, 0.99),
                keep.path    = TRUE,
                keep.fits    = FALSE,
                verbose      = FALSE,
                ...) {
  cl <- match.call()
  method <- match.arg(method)
  estimand <- match.arg(estimand)

  dots <- list(...)
  if (length(dots)) {
    warning(sprintf("Ignoring unused argument(s): %s.",
                    paste0("`", names(dots), "`", collapse = ", ")),
            call. = FALSE)
  }

  ## --- scalar argument validation --------------------------------------------
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) family <- family()
  if (!inherits(family, "family") ||
      !family$family %in% c("gaussian", "binomial")) {
    stop("`family` must be gaussian() or binomial().", call. = FALSE)
  }
  if (!is.numeric(lambda) || !length(lambda) || any(!is.finite(lambda))) {
    stop("`lambda` must be a non-empty numeric vector of finite exponents delta (lambda_n = n^delta).",
         call. = FALSE)
  }
  if (max(lambda) > 3) {
    warning(paste0("max(lambda) > 3: these look like raw penalty values; `lambda` takes EXPONENTS ",
                   "delta with lambda_n = n^delta -- supply exponents delta ",
                   "(convert a raw grid via delta = log(lambda_n)/log(n))."),
            call. = FALSE)
  }
  if (!is.null(gamma)) {
    if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) ||
        gamma <= 0) {
      stop(paste0("`gamma` must be NULL (paired with each delta as 2*(gamma.factor - delta + 1)) ",
                  "or a single positive number (fixed gamma crossed with the lambda grid). ",
                  "No vector form is accepted."), call. = FALSE)
    }
  }
  if (!is.numeric(gamma.factor) || length(gamma.factor) != 1L ||
      !is.finite(gamma.factor)) {
    stop("`gamma.factor` must be a single finite number.", call. = FALSE)
  }
  if (method == "oal") {
    if (!missing(lambda2)) {
      stop("`lambda2` applies only to method = \"goal\"; it is the elastic-net ridge constant of GOAL.",
           call. = FALSE)
    }
  } else {
    if (!is.numeric(lambda2) || !length(lambda2) || any(!is.finite(lambda2)) ||
        any(lambda2 < 0)) {
      stop("`lambda2` must be a non-empty numeric vector of finite values >= 0 (0 nests plain OAL).",
           call. = FALSE)
    }
  }
  if (!is.numeric(clip) || length(clip) != 2L || anyNA(clip) ||
      clip[1L] <= 0 || clip[2L] >= 1 || clip[1L] >= clip[2L]) {
    stop("`clip` must be a length-2 numeric vector with 0 < clip[1] < clip[2] < 1.",
         call. = FALSE)
  }
  for (nm in c("refit", "keep.path", "keep.fits", "verbose")) {
    v <- get(nm, inherits = FALSE)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
    }
  }

  ## --- inputs (D.1-D.3) -------------------------------------------------------
  user.coef <- !is.null(outcome.coef)
  if (missing(outcome)) outcome <- NULL
  inp <- .process_inputs(formula, data, outcome, family,
                         need.outcome = !user.coef)
  n <- inp$n
  p <- inp$p
  treat <- inp$treat
  Xs <- inp$Xs
  .vmsg(verbose, sprintf("n = %d units (%d treated, %d control); p = %d model-matrix column(s).",
                         n, sum(treat == 1L), sum(treat == 0L), p))

  ## --- outcome (weight) model or user hook (D.4-D.5) --------------------------
  if (user.coef) {
    b <- .validate_outcome_coef(outcome.coef, colnames(Xs))
    om <- list(b = b,
               table = data.frame(term = names(b), estimate = as.numeric(b),
                                  row.names = NULL),
               label = "user-supplied", fit = NULL)
    gamma.mode <- "user-coef"
  } else {
    if (identical(family$family, "binomial")) {
      warning(paste0("binomial outcome models are beyond the simulations validated by ",
                     "Shortreed & Ertefaie (2017); GLM theory per Balde (2025, ",
                     "doi:10.1016/j.spl.2025.110379)."), call. = FALSE)
    }
    om <- .fit_outcome_model(inp$y, treat, Xs, family)
    b <- om$b
    gamma.mode <- if (is.null(gamma)) "paired" else "fixed"
  }

  ## --- tuning grid (D.7) -------------------------------------------------------
  deltas <- sort(unique(as.numeric(lambda)))
  l2s <- if (method == "goal") sort(unique(as.numeric(lambda2))) else 0
  grid <- .make_grid(deltas, l2s, gamma, gamma.factor, n, p)
  G <- nrow(grid)
  .vmsg(verbose, sprintf("Grid: %d point(s) (%d delta x %d lambda2); gamma mode: %s.",
                         G, length(deltas), length(l2s), gamma.mode))

  ## --- engine (D.8-D.12) -------------------------------------------------------
  eng <- .run_grid(Xs, treat, b, grid, refit = refit, clip = clip,
                   estimand = estimand, keep.fits = keep.fits,
                   verbose = verbose)

  ## --- selection (D.13): first minimum in grid order ---------------------------
  idx <- .first_min(eng$wamd)
  .vmsg(verbose, sprintf("Selected grid point %d: delta = %s, gamma = %s%s; wAMD = %s (%d covariate(s) retained).",
                         idx, .fmt(grid$delta[idx]), .fmt(grid$gamma[idx]),
                         if (method == "goal") sprintf(", lambda2 = %s", .fmt(grid$lambda2[idx])) else "",
                         .fmt(eng$wamd[idx]), eng$n.selected[idx]))

  ## --- assembly (B.2) ----------------------------------------------------------
  ps <- eng$ps[, idx]
  names(ps) <- inp$rn
  W <- .ipw_from_ps(ps, treat, estimand)
  names(W) <- inp$rn

  cf.std <- eng$coef.std[, idx]
  cf.orig <- .back_transform(cf.std, inp$center, inp$scale)
  names(cf.orig) <- names(cf.std)

  sel.log <- cf.std[-1L] != 0
  pen.sel <- abs(b)^(-grid$gamma[idx])
  selected <- data.frame(term = colnames(Xs),
                         outcome.coef = as.numeric(b),
                         penalty.factor = as.numeric(pen.sel),
                         coef = as.numeric(cf.orig[-1L]),
                         selected = as.logical(sel.log),
                         role = ifelse(sel.log, "retained", "excluded"),
                         row.names = NULL)

  lambda.sel <- c(delta = grid$delta[idx],
                  lambda.n = grid$lambda.n[idx],
                  gamma = grid$gamma[idx],
                  lambda2 = if (method == "goal") grid$lambda2[idx] else NA_real_)

  path <- NULL
  sets <- NULL
  coef.path <- NULL
  if (keep.path) {
    path <- data.frame(lambda2 = if (method == "goal") grid$lambda2 else NA_real_,
                       delta = grid$delta,
                       lambda.n = grid$lambda.n,
                       gamma = grid$gamma,
                       s = eng$s,
                       wamd = eng$wamd,
                       n.selected = eng$n.selected,
                       refit.converged = if (refit) eng$refit.converged else NA,
                       selected = seq_len(G) == idx,
                       row.names = NULL)
    sets <- eng$sets
    coef.path <- apply(eng$coef.std, 2L, .back_transform,
                       center = inp$center, scale = inp$scale)
    rownames(coef.path) <- names(cf.std)
    colnames(coef.path) <- sprintf("lambda2=%s,delta=%s",
                                   format(path$lambda2), format(path$delta))
  }

  balance <- .balance_table(inp$X, treat, W, .sd_denom(estimand), inp$bin.vars)

  covs.out <- inp$X
  attr(covs.out, "bin.vars") <- inp$bin.vars

  fits <- NULL
  if (keep.fits) {
    fits <- list(glmnet = eng$fits,
                 outcome = om$fit,
                 refit = if (refit) .refit_glm(Xs, treat, sel.log)$fit else NULL)
  }

  info <- list(n = n,
               n.treated = sum(treat == 1L),
               p = p,
               grid.size = G,
               family = if (user.coef) NA_character_ else family$family,
               gamma.mode = gamma.mode,
               center = inp$center,
               scale = inp$scale,
               contrasts = inp$contrasts,
               glmnet.version = tryCatch(as.character(utils::packageVersion("glmnet")),
                                         error = function(e) NA_character_),
               oalasso.version = tryCatch(as.character(utils::packageVersion("oalasso")),
                                          error = function(e) NA_character_))

  provenance <- .provenance(method,
                            lambda2.user = (method == "goal" && !missing(lambda2)),
                            user.coef = user.coef,
                            family = if (user.coef) "gaussian" else family$family,
                            refit = refit)

  structure(list(ps = ps,
                 weights = W,
                 coefficients = cf.orig,
                 coefficients.std = cf.std,
                 selected = selected,
                 lambda = lambda.sel,
                 criterion = "wamd",
                 criterion.value = eng$wamd[idx],
                 path = path,
                 sets = sets,
                 coef.path = coef.path,
                 outcome.coef = b,
                 outcome.model = list(coefficients = om$table,
                                      family = if (user.coef) NA_character_ else family$family,
                                      label = om$label),
                 balance = balance,
                 provenance = provenance,
                 treat = treat,
                 covs = covs.out,
                 estimand = estimand,
                 method = method,
                 refit = refit,
                 clip = clip,
                 formula = formula,
                 data = inp$data,
                 terms = inp$terms,
                 xlevels = inp$xlevels,
                 fits = fits,
                 info = info,
                 call = cl),
            class = "oal")
}
