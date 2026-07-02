# Outcome-adaptive lasso propensity scores

`oal()` estimates a propensity score by the outcome-adaptive lasso (OAL)
of Shortreed and Ertefaie (2017) or its generalized elastic-net form
(GOAL) of Balde, Yang and Lefebvre (2023): an adaptive lasso on the
treatment (propensity score) log-likelihood whose penalty weights come
from an *outcome* regression, so that covariates unrelated to the
outcome – instruments and noise variables – are excluded from the
propensity score model. Tuning parameters are selected by the papers'
weighted absolute mean difference (wAMD) balance criterion. The result
is deliberately modest: a numeric score vector designed to be handed to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
as `distance`, to
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
as `ps`, or to
[`psAve::psave()`](https://kabajiro.github.io/psAve/reference/psave.html)
as an appended candidate.

## Usage

``` r
oal(
  formula,
  data,
  outcome,
  method = c("oal", "goal"),
  estimand = c("ATE", "ATT"),
  family = gaussian(),
  lambda = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
  gamma = NULL,
  gamma.factor = 2,
  lambda2 = c(0, 10^c(-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1)),
  outcome.coef = NULL,
  refit = FALSE,
  clip = c(0.01, 0.99),
  keep.path = TRUE,
  keep.fits = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- formula:

  A two-sided formula `treat ~ x1 + x2 + ...`, exactly as in
  [`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html).
  The right-hand side defines the propensity score covariates; the same
  covariates (plus the treatment) form the outcome model that supplies
  the adaptive penalty weights.

- data:

  A data frame containing the variables in `formula` and `outcome`.
  Complete cases in all used variables are REQUIRED; any missing value
  is an error naming the offending variables, never a silent row drop.

- outcome:

  A one-sided formula `~ y` naming the outcome variable (exactly one
  variable). Required unless `outcome.coef` is supplied, in which case
  it is unused. A two-sided formula is an error: OAL derives its penalty
  weights from an outcome model on the *same* covariates as the
  propensity score model, so there is no separate outcome design to
  specify.

- method:

  `"oal"` (default): the outcome-adaptive lasso; `"goal"`: the
  generalized outcome-adaptive lasso, which adds an elastic-net
  \\\lambda_2 \sum_j \alpha_j^2\\ term (grouping effect and numerical
  stability under correlated covariates and near-positivity violations).

- estimand:

  `"ATE"` (default) or `"ATT"`; determines the inverse-probability
  weights used inside the wAMD criterion and returned in `weights`. The
  ATE default is Shortreed and Ertefaie's wAMD weighting and
  deliberately differs from
  [`psAve::psave()`](https://kabajiro.github.io/psAve/reference/psave.html)'s
  ATT default.

- family:

  The outcome (weight) model family:
  [`gaussian()`](https://rdrr.io/r/stats/family.html) (default; the
  specification validated by the papers' simulations) or
  [`binomial()`](https://rdrr.io/r/stats/family.html).
  [`binomial()`](https://rdrr.io/r/stats/family.html) emits a one-time
  warning: binomial outcome models are beyond the simulations validated
  by Shortreed and Ertefaie (2017); GLM theory per Balde (2025).

- lambda:

  Numeric vector of EXPONENTS \\\delta\\; the penalties actually applied
  are \\\lambda_n = n^{\delta}\\. The default is Shortreed and
  Ertefaie's grid `c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)`.
  A guard warns when `max(lambda) > 3`: such values look like raw
  penalty values, not exponents (convert a raw grid via
  `delta = log(lambda_n)/log(n)`).

- gamma:

  `NULL` (default) pairs \\\gamma\\ with each \\\delta\\ as \\\gamma =
  2(\code{gamma.factor} - \delta + 1)\\ (the
  Shortreed-Ertefaie/rejoinder convention); a single number (e.g. `2.5`)
  fixes \\\gamma\\ and crosses it with the whole `lambda` grid (the
  Schnitzer-style fixed-\\\gamma\\ mode). No vector form is accepted.

- gamma.factor:

  The convergence factor `gcf` in the pairing formula (default `2`);
  ignored when `gamma` is given.

- lambda2:

  Numeric grid of elastic-net ridge constants for `method = "goal"` only
  (supplying it otherwise is an error). Default
  `c(0, 10^c(-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1))` – the
  author's published grid (Balde 2025 supplement, "taken from Zou and
  Hastie (2005)"), verified against the official GOAL code on
  2026-07-02; `0` is always a member and nests plain OAL exactly.
  `lambda2` is selected jointly with `(lambda, gamma)` by the wAMD.

- outcome.coef:

  The extension hook: a NAMED numeric vector of outcome coefficients
  \\b_j\\ over the **standardized** model-matrix columns (names must
  exactly cover the columns), overriding the internal outcome model.
  **The values are interpreted on the standardized scale of the internal
  design matrix** – adaptive weights are scale-dependent. All values
  must be finite and not all zero. Zeros ARE allowed and mean a hard
  drop: \\\|0\|^{-\gamma} = \infty\\ and glmnet natively converts an
  infinite penalty factor to an exclusion. No epsilon or cap is ever
  applied. Using this hook flips the provenance label to "screening use
  only" (see Details) and sets `info$gamma.mode = "user-coef"`. If
  `outcome` is also supplied, it is still validated and the
  outcome-leakage guard still runs (the outcome variable may not appear
  among the PS covariates); if `outcome` is omitted, the package
  **cannot** detect outcome leakage in `formula` – ensuring the
  coefficients were derived without post-treatment or outcome
  information is then entirely the caller's responsibility.

- refit:

  `FALSE` (default): the propensity score comes from the penalized fit
  itself (the Shortreed-Ertefaie/rejoinder convention). `TRUE`: at every
  grid point an *unpenalized* logistic regression is refit on the
  selected covariates, and both the propensity score and the wAMD come
  from the refit (the Schnitzer et al. 2025 variant, adapted from a
  longitudinal setting); per-grid-point convergence is recorded in
  `path$refit.converged` and any non-convergence triggers one warning –
  never a silent fallback.

- clip:

  Length-2 numeric: propensity scores are clipped to
  `[clip[1], clip[2]]` BEFORE the wAMD IPW weights are formed and in the
  returned `ps` (default `c(0.01, 0.99)`, equal to
  [`psAve::psave()`](https://kabajiro.github.io/psAve/reference/psave.html)'s
  default so that psave's re-clipping is a no-op). The default is a
  Shortreed-Ertefaie-lineage safety choice; Balde's official reference
  code runs UNCLIPPED weights (`1/e`, `1/(1 - e)` with no truncation).
  To reproduce that behavior, effectively disable clipping with
  `clip = c(1e-12, 1 - 1e-12)`.

- keep.path:

  If `TRUE` (default), the full per-grid-point results are stored in
  `path`, `sets`, and `coef.path`.

- keep.fits:

  If `TRUE`, the fitted glmnet paths, the outcome model, and (with
  `refit = TRUE`) the selected refit are retained in `fits`. Default
  `FALSE`.

- verbose:

  If `TRUE`, progress messages report the grid evaluation and the
  selected point.

- ...:

  Reserved for future use; supplying unused arguments triggers a
  warning.

## Value

An object of class `"oal"`: a list with components

- `ps`:

  numeric(n), named by `rownames(data)`: the propensity score, clipped
  to `clip`, strictly in (0, 1) – **the deliverable**.

- `weights`:

  numeric(n): the IPW weights at `estimand` implied by `ps`.

- `coefficients`:

  named numeric(p+1): selected propensity score coefficients on the
  ORIGINAL covariate scale (intercept first).

- `coefficients.std`:

  the same on the standardized scale actually optimized.

- `selected`:

  data frame, one row per model-matrix column: `term`, `outcome.coef`
  (signed \\b_j\\, standardized scale), `penalty.factor`, `coef`
  (original scale), `selected` (logical), and `role`
  (`"retained"`/`"excluded"`).

- `lambda`:

  named numeric: `delta`, `lambda.n` (\\= n^\delta\\; \\= (n +
  p)^\delta\\ on GOAL grid points with \\\lambda_2 \> 0\\), `gamma`,
  `lambda2` (`NA` for `method = "oal"`) at the selected point.

- `criterion`, `criterion.value`:

  `"wamd"` and its value at the selection.

- `path`:

  data frame (or `NULL` without `keep.path`), one row per grid point in
  grid order: `lambda2`, `delta`, `lambda.n`, `gamma`, `s` (the exact
  glmnet `s` used), `wamd`, `n.selected`, `refit.converged` (`NA` unless
  `refit`), `selected`.

- `sets`:

  list (or `NULL`): the selected variable set at every grid point.

- `coef.path`:

  (p+1) x G matrix (or `NULL`): propensity score coefficients across the
  grid, original scale.

- `outcome.coef`:

  named numeric(p): the \\b_j\\ actually used (internal or user),
  standardized scale.

- `outcome.model`:

  compact list: coefficient table, `family`, and a `label`
  (`"lm(y ~ A + X)"`, `"glm(y ~ A + X, binomial)"`, or
  `"user-supplied"`); the full fit object only under `keep.fits`.

- `balance`:

  data frame per ORIGINAL covariate column: `smd.un`, `smd.wt`, `ks.un`,
  `ks.wt` via cobalt – identical layout to `psAve::psave()$balance`.

- `provenance`:

  the fixed provenance label, printed first by
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html).

- `treat`:

  integer(n) 0/1 treatment as coerced.

- `covs`:

  n x p numeric original-scale model matrix with `attr(, "bin.vars")`.

- `estimand`, `method`, `refit`, `clip`:

  as resolved.

- `formula`, `data`, `terms`, `xlevels`:

  stored to power
  [`oal_match()`](https://kabajiro.github.io/oalasso/reference/oal_match.md),
  [`oal_weight()`](https://kabajiro.github.io/oalasso/reference/oal_weight.md)
  and
  [`predict.oal()`](https://kabajiro.github.io/oalasso/reference/predict.oal.md).

- `fits`:

  list (or `NULL`): `glmnet` (one path per \\(\lambda_2, \gamma)\\),
  `outcome`, `refit` – iff `keep.fits`.

- `info`:

  list: `n`, `n.treated`, `p`, `grid.size`, `family`, `gamma.mode`
  (`"paired"`/`"fixed"`/`"user-coef"`), `center`, `scale`, `contrasts`,
  `glmnet.version`, `oalasso.version`.

- `call`:

  the matched call.

## Details

### Objective and exact glmnet evaluation

OAL solves, on the TOTAL log-likelihood scale, \$\$\hat\alpha =
\arg\min\_\alpha\\ -\ell(\alpha; A, X_s) + \lambda_n \sum_j
\|b_j\|^{-\gamma} \|\alpha_j\|,\$\$ where \\\ell\\ is the binomial
log-likelihood of treatment on the standardized covariates and \\b_j\\
are the outcome-model coefficients. glmnet minimizes the
*per-observation* loss and internally rescales penalty factors to sum to
the number of variables, so each grid point is evaluated **exactly** at
\$\$s = \overline{\mathrm{pen}} \cdot \lambda_n / n, \qquad
\mathrm{pen}\_j = \|b_j\|^{-\gamma},\$\$ via
`coef(fit, s = s, exact = TRUE, x = , y = , penalty.factor = )`, which
re-solves at that penalty (the KKT-verified recipe of the
Shortreed/Ertefaie rejoinder code). When `outcome.coef` contains zeros,
some \\\mathrm{pen}\_j = \infty\\: glmnet converts those to exclusions
and internally resets their factor to 1 *before* the sum-to-nvars
rescale, so the constant generalizes to \\s = (\sum\_{j:
\mathrm{finite}} \mathrm{pen}\_j + \\\\\mathrm{pen}\_j = \infty\\)/p
\cdot \lambda_n / n\\ – identical to
\\\overline{\mathrm{pen}}\\\lambda_n/n\\ when all factors are finite.
This correction is always on; it is what makes the \\n^\delta\\ grid
carry its published meaning.

### Standardization protocol (fixed, not an argument)

The model matrix `X <- model.matrix(formula, data)[, -1]` (default
treatment contrasts) is standardized ONCE, `Xs <- scale(X, TRUE, TRUE)`,
all columns including dummies (the rejoinder convention). The outcome
model AND glmnet are both fit on this SAME `Xs` with
`glmnet(standardize = FALSE)`. Per-subset scaling is structurally
impossible. Coefficients are back-transformed for output: \\\alpha_j =
\alpha^{std}\_j / s_j\\ and \\\alpha_0 = \alpha^{std}\_0 - \sum_j
\alpha^{std}\_j c_j / s_j\\, with centers \\c\\ and scales \\s\\ stored
in `info`.

### Weight (outcome) model and the degenerate-beta policy

By default \\b\\ comes from `lm(y ~ A + Xs)` (or
`glm(y ~ A + Xs, binomial)`) on the FULL sample, both arms pooled – the
Shortreed-Ertefaie/rejoinder convention (not control-arm-only). Only the
covariate coefficients are used; the treatment coefficient is discarded
and the treatment is never penalized anywhere. Any `NA` or *exactly
zero* internal coefficient is an ERROR naming the columns: OLS betas are
almost surely nonzero, so an exact zero or `NA` signals rank deficiency
or aliasing that would silently corrupt the penalty-scale correction.
Remove the collinear columns or supply `outcome.coef`.

### Tuning grid and the \\\gamma\\ pairing

`lambda` supplies exponents \\\delta\\ with \\\lambda_n = n^{\delta}\\,
`n = nrow(model.matrix)` always. With `gamma = NULL` each \\\delta\\ is
paired with \\\gamma = 2(\mathrm{gcf} - \delta + 1)\\ (default
\\\mathrm{gcf} = 2\\), the formula derivable from \\\lambda_n\\
n^{\gamma/2 - 1} = n^{\mathrm{gcf}}\\; there are no \\\lambda \times
\gamma\\ cross-products in paired mode. A scalar `gamma` is crossed with
all \\\delta\\ (Schnitzer's \\\gamma = 2.5\\ is reachable this way;
convert their raw \\\lambda\\ grid via \\\delta = \log \lambda_n / \log
n\\).

### GOAL

For \\\lambda_2 \> 0\\ the elastic-net term is solved by the Zou-Hastie
augmentation: \\X\_{aug} = \mathrm{rbind}(X_s, \sqrt{\lambda_2} I_p)\\
with \\p\\ pseudo-responses 0, followed by rescaling ALL coefficients
INCLUDING the intercept by \\(1 + \lambda_2)\\ (Balde 2025 supplement:
\\PS = \mathrm{expit}(\mathrm{cbind}(1, X)\\(1 +
\lambda_2)\\\hat\alpha)\\); the propensity score uses the original \\n\\
rows only. On augmented grid points the raw penalty constant is
\\\lambda_n = (n + q)^{\delta}\\ with \\q = p\\ augmentation rows
(Balde's `adaptive.lasso(lambda = n.q^(il))`, `n.q = n + q`), and the
penalty-scale constant uses \\1/(n + q)\\ in place of \\1/n\\, i.e. \\s
= \overline{\mathrm{pen}} (n + q)^{\delta} / (n + q)\\. \\\lambda_2 =
0\\ grid points run the plain-OAL code path (no augmentation,
\\\lambda_n = n^{\delta}\\), so `method = "goal"` with `lambda2 = 0`
nests `method = "oal"` exactly. (Balde's script literally augments with
\\q\\ zero rows and the \\(n+q)^\delta\\ base even at \\\lambda_2 = 0\\;
the zero rows still shift the intercept, so exact nesting is
deliberately preferred here.) \\\lambda_2\\ is selected jointly with
\\(\delta, \gamma)\\ by the wAMD: the flat first-minimum over the
\\\lambda_2\\-outer \\\times\\ \\\delta\\-inner grid, which is
equivalent to Balde's nested rule (per-\\\lambda_2\\ minima over
\\(\delta, \gamma)\\, then the first minimum over \\\lambda_2\\).

### wAMD criterion and selection

Every grid point is scored by
[`oal_wamd()`](https://kabajiro.github.io/oalasso/reference/oal_wamd.md):
propensity scores clipped to `clip` first (note that Balde's reference
implementation is unclipped; see the `clip` argument), IPW weights at
`estimand`, per-column weighted absolute mean differences on the
standardized matrix summed with weights \\\|b_j\|\\ over ALL columns
(including excluded ones). The selected grid point is the FIRST minimum
(relative tolerance `1e-9`) in the grid order \\\lambda_2\\ ascending (0
first) outer \\\times\\ \\\delta\\ ascending inner – ties therefore
prefer plain OAL and then the smallest penalty. The wAMD is a balance
heuristic: no published proof guarantees it lands \\\lambda_n\\ inside
the window required by the theory (\\\lambda_n/\sqrt{n} \to 0\\,
\\\lambda_n n^{\gamma/2 - 1} \to \infty\\, \\\gamma \> 1\\), and it
inherits the outcome model's misspecification vulnerability.

### Why exclude instruments?

Covariates that predict treatment but not outcome should be excluded
from a propensity score model: they increase the variance of the effect
estimate without decreasing bias (Brookhart et al. 2006) and amplify the
bias from any unmeasured confounding (Myers et al. 2011). OAL's penalty
implements exactly this exclusion.

### Provenance labels

[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) always print a
provenance line first:
`"OAL (Shortreed & Ertefaie 2017, doi:10.1111/biom.12679)"` or
`"GOAL (Balde, Yang & Lefebvre 2023, doi:10.1111/biom.13683)"` with a
clause stating whether the `lambda2` grid is the author's published grid
(Balde 2025 supplement) or user-specified. Supplying `outcome.coef`
appends
`"user-supplied outcome coefficients — screening use only; no oracle property"`;
`family = binomial()` and `refit = TRUE` append their own flags.

### Handoff contract

`$ps` is a numeric vector of length `n`, named by `rownames(data)`,
strictly inside (0, 1) after clipping – it satisfies
`psAve::psave(ps.append = )`, `MatchIt::matchit(distance = )`, and
`WeightIt::weightit(ps = )` by construction.

There is no `seed` argument: the whole pipeline is deterministic.

## References

Shortreed SM, Ertefaie A (2017). Outcome-adaptive lasso: variable
selection for causal inference. *Biometrics*, 73(4), 1111-1122.
[doi:10.1111/biom.12679](https://doi.org/10.1111/biom.12679)

Balde I, Yang Y, Lefebvre G (2023). Reader reaction to "Outcome-adaptive
lasso: variable selection for causal inference" by Shortreed and
Ertefaie (2017). *Biometrics*, 79(1), 514-520.
[doi:10.1111/biom.13683](https://doi.org/10.1111/biom.13683)

Jones B, Ertefaie A, Shortreed SM (2023). Rejoinder to reader reaction
"On the use of the outcome-adaptive lasso for propensity score
estimation". *Biometrics*, 79(1), 521-525.
[doi:10.1111/biom.13681](https://doi.org/10.1111/biom.13681)

Balde I (2025). Oracle properties of the generalized outcome-adaptive
lasso. *Statistics & Probability Letters*.
[doi:10.1016/j.spl.2025.110379](https://doi.org/10.1016/j.spl.2025.110379)

Zou H (2006). The adaptive lasso and its oracle properties. *Journal of
the American Statistical Association*, 101(476), 1418-1429.
[doi:10.1198/016214506000000735](https://doi.org/10.1198/016214506000000735)

Brookhart MA, Schneeweiss S, Rothman KJ, Glynn RJ, Avorn J, Sturmer T
(2006). Variable selection for propensity score models. *American
Journal of Epidemiology*, 163(12), 1149-1156.
[doi:10.1093/aje/kwj149](https://doi.org/10.1093/aje/kwj149)

Myers JA, Rassen JA, Gagne JJ, Huybrechts KF, Schneeweiss S, Rothman KJ,
Joffe MM, Glynn RJ (2011). Effects of adjusting for instrumental
variables on bias and precision of effect estimates. *American Journal
of Epidemiology*, 174(11), 1213-1222.
[doi:10.1093/aje/kwr364](https://doi.org/10.1093/aje/kwr364)

Schnitzer ME, Talbot D, Liu Y, Berger C, Wang G, O'Loughlin J, Sylvestre
M-P, Ertefaie A (2025). Outcome-adaptive LASSO for longitudinal data.
*Statistics in Medicine* (arXiv:2410.08283).

## See also

[`oal_match()`](https://kabajiro.github.io/oalasso/reference/oal_match.md),
[`oal_weight()`](https://kabajiro.github.io/oalasso/reference/oal_weight.md),
[`oal_wamd()`](https://kabajiro.github.io/oalasso/reference/oal_wamd.md),
[`predict.oal()`](https://kabajiro.github.io/oalasso/reference/predict.oal.md),
[`plot.oal()`](https://kabajiro.github.io/oalasso/reference/plot.oal.md),
[`bal.tab.oal()`](https://kabajiro.github.io/oalasso/reference/bal.tab.oal.md)

## Examples

``` r
data("lalonde", package = "MatchIt")

## 1) Matching
fit <- oal(treat ~ age + educ + race + married + nodegree + re74 + re75,
           data = lalonde, outcome = ~ re78)
fit
#> OAL (Shortreed & Ertefaie 2017, doi:10.1111/biom.12679)
#> An oal object (outcome-adaptive lasso propensity score)
#>  - estimand: ATE;  refit: FALSE (penalized fit)
#>  - sample:   614 units (185 treated, 429 control); 8 covariate column(s)
#>  - selected: delta = -10 (lambda_n = 1.31e-28), gamma = 26
#>  - wAMD at the selection: 877
#> 
#> Retained (outcome-related):
#>   age  (|b| =  128)
#>   educ  (|b| = 1062)
#>   racehispan  (|b| =  560)
#>   racewhite  (|b| =  621)
#>   married  (|b| =  201)
#>   nodegree  (|b| =  126)
#>   re74  (|b| = 1920)
#>   re75  (|b| =  763)
#> Excluded (outcome-unrelated):
#>   (none)
#> 
#> Next:
#>   MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde, distance = x$ps)
#>     or: oal_match(x)
#>   WeightIt::weightit(treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde, ps = x$ps, estimand = "ATE")
#>     or: oal_weight(x)
#>   ($ps is numeric, named by rownames(lalonde), strictly inside (0, 1);
#>    it also satisfies psAve::psave(ps.append = cbind(oal = x$ps)).)
m <- oal_match(fit)   # = MatchIt::matchit(f, lalonde, distance = fit$ps)

## 2) Weighting (requires WeightIt)
# w <- oal_weight(fit)  # = WeightIt::weightit(f, lalonde, ps = fit$ps, "ATE")

## 3) psAve composition (requires psAve); cbind() labels the candidate
##    while keeping the row names aligned with rownames(data)
# ma <- psAve::psave(treat ~ age + educ + race + married + nodegree +
#                      re74 + re75,
#                    data = lalonde, outcome = ~ re78,
#                    ps.append = cbind(oal = fit$ps))
```
