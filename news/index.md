# Changelog

## oalasso 1.0.0

Initial release, implementing the outcome-adaptive lasso of Shortreed &
Ertefaie (2017), *Biometrics* 73(4):1111–1122, <doi:10.1111/biom.12679>,
and the generalized outcome-adaptive lasso (GOAL) of Baldé, Yang &
Lefebvre (2023), *Biometrics* 79(1):514–520, <doi:10.1111/biom.13683>.

- [`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md):
  outcome-adaptive lasso propensity scores. Adaptive penalty weights
  `|b_j|^-gamma` from an unpenalized outcome model `Y ~ A + X` on the
  full sample ([`gaussian()`](https://rdrr.io/r/stats/family.html)
  validated; [`binomial()`](https://rdrr.io/r/stats/family.html) allowed
  with a one-time warning, GLM theory per Baldé 2025,
  <doi:10.1016/j.spl.2025.110379>); tuning over Shortreed & Ertefaie’s
  exponent grid `lambda_n = n^delta` with the paired penalty exponent
  `gamma = 2*(gamma.factor - delta + 1)`; selection by the weighted
  absolute mean difference (wAMD) with `|b_j|` weights, ATE (default) or
  ATT weighting; ties take the first minimum in documented grid order.
- Exact reproduction of the published objective on `glmnet`: every grid
  point is evaluated at `s = mean(pen) * lambda_n / n` with
  `coef(..., exact = TRUE)`, undoing glmnet’s internal sum-to-nvars
  penalty-factor rescaling and its per-observation loss scale (the
  rejoinder recipe of Jones, Ertefaie & Shortreed 2023,
  <doi:10.1111/biom.13681>). The KKT conditions of the target objective
  are the acceptance criterion in the test suite, with a drift guard
  against the installed glmnet version.
- `method = "goal"`: Baldé’s elastic-net generalization via the
  Zou–Hastie data augmentation with the `(1 + lambda2)` rescale of ALL
  coefficients including the intercept (Baldé 2025 supplement:
  `expit(cbind(1, X) %*% (1 + lambda2) * coef)`); on augmented grid
  points (`lambda2 > 0`) the raw penalty constant is `(n + q)^delta`
  with `q = p` augmentation rows, matching the author’s
  `adaptive.lasso(lambda = n.q^(il))` with `n.q = n + q`. `lambda2` is
  selected jointly with `(delta, gamma)` by the wAMD — the flat
  first-minimum, equivalent to Baldé’s nested
  per-`lambda2`-then-over-`lambda2` rule — and `lambda2 = 0` exactly
  nests plain OAL (no augmentation on that code path). The default
  `lambda2` grid
  `c(0, 10^c(-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1))` is the
  author’s published grid (Baldé 2025 supplement, “taken from Zou and
  Hastie (2005)”), verified against the official GOAL code on
  2026-07-02.
- Baldé’s reference code runs UNCLIPPED IPW weights; the default
  `clip = c(0.01, 0.99)` is a Shortreed–Ertefaie-lineage safety choice,
  and `clip = c(1e-12, 1 - 1e-12)` effectively disables clipping to
  reproduce the reference behavior.
- Fixed standardization protocol: covariates are standardized once
  (training statistics, all columns including dummies); the outcome
  model and glmnet fit on the same matrix with `standardize = FALSE`;
  coefficients returned on both scales.
- Options with provenance labels printed on every fit: `refit = TRUE`
  (post-selection unpenalized refit, Schnitzer et al. 2025,
  <doi:10.1002/sim.70316>, adapted from a longitudinal setting;
  convergence recorded per grid point), scalar `gamma` (Schnitzer-style
  fixed exponent crossed with the delta grid), and the `outcome.coef`
  hook (user-supplied weights on the standardized scale — screening use
  only; zeros become hard exclusions via glmnet’s Inf-to-exclude).
- Ecosystem integration: `fit$ps` (numeric, named by `rownames(data)`,
  clipped to `clip = c(0.01, 0.99)`) drops into
  `MatchIt::matchit(distance = )`, `WeightIt::weightit(ps = )`, and
  `psAve::psave(ps.append = )`;
  [`oal_match()`](https://kabajiro.github.io/oalasso/reference/oal_match.md)
  /
  [`oal_weight()`](https://kabajiro.github.io/oalasso/reference/oal_weight.md)
  wrappers reuse the stored formula and data;
  [`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
  works on `oal` objects.
- Full S3 suite: [`print()`](https://rdrr.io/r/base/print.html)
  (provenance line first; retained vs excluded covariates with the
  instrument-exclusion rationale of Brookhart et al. 2006 and Myers et
  al. 2011; the literal next call; near-positivity diagnostic when \> 5%
  of scores sit at a clip bound),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (`"wamd"`,
  `"coef"`, `"balance"`), [`coef()`](https://rdrr.io/r/stats/coef.html)
  (both scales), [`predict()`](https://rdrr.io/r/stats/predict.html)
  (works without `keep.fits`),
  [`weights()`](https://rdrr.io/r/stats/weights.html).
- Exported criterion:
  [`oal_wamd()`](https://kabajiro.github.io/oalasso/reference/oal_wamd.md)
  scores any candidate propensity score on the exact S&E wAMD formula
  (single source of truth; never delegated to cobalt).
- Policies: any `NA` in used variables is an error naming the variables
  (never a silent drop); zero-variance and aliased/rank-deficient
  covariates error; fully deterministic pipeline (no seed argument);
  English-only messages.
- Documented divergences from the author’s legacy `fAL`/`fOAL` code
  (which is superseded, not reproduced): the `fAL` gamma transcription
  bug, the `varlist` wAMD bug, global-`n` scoping, per-subset
  [`scale()`](https://rdrr.io/r/base/scale.html), and the archived `lqa`
  dependency. See
  [`vignette("method-details", package = "oalasso")`](https://kabajiro.github.io/oalasso/articles/method-details.md).
- Three vignettes: Getting Started (matching and weighting workflows,
  reading the output), Using oalasso with psAve (the `ps.append`
  composition and its caveats), and Method details and provenance (all
  formulas, the exact penalty-scale correction, provenance table,
  relation to other software, honest nonlinear-extension notes).
