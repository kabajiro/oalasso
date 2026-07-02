# Method details and provenance

This vignette documents the method implemented by
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) at the
level of formulas, explains the one nontrivial piece of numerical
machinery (the exact penalty-scale correction for `glmnet`), records the
provenance of every default (which choices are Shortreed & Ertefaie’s,
which are Schnitzer’s, which are the package’s), documents the
differences from the author’s own legacy code, situates `oalasso` among
related software, and closes with an honest section on nonlinear
extensions. Throughout, “S&E” is Shortreed & Ertefaie (2017) and “the
rejoinder” is Jones, Ertefaie & Shortreed (2023), whose companion code
is the reference implementation this package reproduces.

## Notation

For subjects $`i = 1, \dots, n`$: covariates $`X_i`$ (a $`p`$-column
numeric matrix after full dummy expansion of factors), binary treatment
$`A_i \in \{0, 1\}`$, and outcome $`Y_i`$. The propensity score model is
logistic,
$`e(X_i) = \Pr(A_i = 1 \mid X_i) = \operatorname{expit}(\alpha_0 + X_i^\top \alpha)`$.
Covariates are conceptually partitioned into confounders (predict $`A`$
and $`Y`$), pure outcome predictors (predict $`Y`$ only), instruments
(predict $`A`$ only), and noise. All fitting is done on the standardized
matrix $`X^s`$ (below); coefficients are back-transformed to the
original scale for output.

## The OAL estimator

### Step 1: the outcome (weight) model

An unpenalized regression of the outcome on treatment and *all* PS-model
covariates, fit on the full sample (both arms pooled — the S&E/rejoinder
convention):

``` math
Y_i = \eta A_i + \sum_{j=1}^{p} \beta_j X^s_{ij} + \varepsilon_i
\quad\text{(via } \texttt{lm} \text{ for } \texttt{gaussian()}\text{, } \texttt{glm} \text{ for } \texttt{binomial()}\text{)}.
```

Only the covariate coefficients $`\hat\beta_1, \dots, \hat\beta_p`$ are
used; the treatment and intercept coefficients are discarded, and the
treatment is never penalized anywhere. A
[`binomial()`](https://rdrr.io/r/stats/family.html) outcome family is
allowed with a one-time warning: S&E’s validating simulations used
continuous outcomes; the GLM extension of the theory is Baldé (2025,
<doi:10.1016/j.spl.2025.110379>).

Because OLS coefficients are almost surely nonzero, an exactly-zero or
`NA` $`\hat\beta_j`$ from the internal model signals a data pathology
(rank deficiency, aliased columns) that would silently corrupt the
penalty scaling below —
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) therefore
stops with an error naming the offending columns rather than continuing.

### Step 2: the adaptively penalized treatment model

``` math
\hat\alpha(\lambda_n, \gamma) \;=\; \arg\min_\alpha \;
-\ell_n(\alpha; A, X^s) \;+\; \lambda_n \sum_{j=1}^{p} \hat w_j\, |\alpha_j|,
\qquad
\hat w_j = |\hat\beta_j|^{-\gamma},
```

where $`\ell_n`$ is the *total* (not per-observation) binomial
log-likelihood. The mechanism: $`\hat\beta_j`$ is
$`\sqrt n`$-consistent, so for outcome-related covariates
$`\hat w_j = O_p(1)`$, while for outcome-unrelated covariates
(instruments, noise) $`|\hat\beta_j| = O_p(n^{-1/2})`$ and the penalty
weight explodes at rate $`n^{\gamma/2}`$ — the outcome model decides
which coefficients get an exploding penalty in the treatment model.
S&E’s Theorem 1 gives the oracle-type result under the rate conditions

``` math
\lambda_n / \sqrt{n} \to 0,
\qquad
\lambda_n\, n^{\gamma/2 - 1} \to \infty,
\qquad
\gamma > 1 :
```

outcome-unrelated coefficients are set to zero with probability tending
to one, and the retained coefficients are $`\sqrt n`$-asymptotically
normal with the oracle variance. Note the second condition is *heavier*
than the standard adaptive-lasso rate
($`\lambda_n n^{(\gamma-1)/2} \to \infty`$): an instrument has a
genuinely nonzero treatment-model coefficient, so the treatment
likelihood pulls it in, and only a heavier penalty overrides that pull.
$`\gamma > 1`$ is exactly what makes the window between the two
conditions nonempty. These guarantees are conditional on the linear
outcome model being correctly specified — see the final section.

### Step 3: the tuning grid

`lambda` supplies **exponents** $`\delta`$, with
$`\lambda_n = n^\delta`$ and $`n`$ always the number of rows of the
model matrix. The default is S&E’s grid:

``` math
\delta \in \{-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49\}.
```

Each $`\delta`$ is **paired** with a penalty exponent (never crossed, so
there is no silent grid explosion):

\$\$ \gamma\_\delta = 2\\(g\_{cf} - \delta + 1), \qquad g\_{cf} = 2
\text{ (the \texttt{gamma.factor} argument)}, \$\$

which is the value solving
$`\lambda_n \, n^{\gamma/2 - 1} = n^{g_{cf}}`$ — every grid point sits
on the same divergence rate for the second oracle condition while
$`\lambda_n`$ itself ranges from negligible to nearly $`\sqrt n`$. This
is the S&E / rejoinder convention. Supplying a scalar `gamma`
(e.g. `gamma = 2.5`) instead fixes $`\gamma`$ and crosses it with the
whole $`\delta`$ grid — Schnitzer et al.’s (2025) configuration; their
raw lambda grid converts to exponents via
$`\delta = \log(\lambda_n)/\log(n)`$. A guard warns if `lambda` values
look like raw penalties rather than exponents (`max(lambda) > 3`).

### Step 4: the wAMD selection criterion

For each grid point, the fitted scores are clipped to `clip` (default
$`[0.01, 0.99]`$) and converted to IPW weights at the requested
estimand:

``` math
\text{ATE:}\quad w_i = \frac{A_i}{e_i} + \frac{1 - A_i}{1 - e_i},
\qquad\qquad
\text{ATT:}\quad w_i = A_i + (1 - A_i)\,\frac{e_i}{1 - e_i}.
```

The weighted absolute mean difference is then, over **all** $`p`$
standardized columns (including excluded ones),

``` math
\mathrm{wAMD} \;=\; \sum_{j=1}^{p} |\hat\beta_j| \cdot
\left| \frac{\sum_i w_i X^s_{ij} A_i}{\sum_i w_i A_i}
     - \frac{\sum_i w_i X^s_{ij} (1 - A_i)}{\sum_i w_i (1 - A_i)} \right|,
```

and the selected grid point is the argmin — the first minimum within a
$`10^{-9}`$ relative tolerance, in the documented grid order (`lambda2`
ascending with 0 first, then $`\delta`$ ascending), so ties prefer plain
OAL and then the smallest penalty. This formula is implemented once,
natively, in the exported
[`oal_wamd()`](https://kabajiro.github.io/oalasso/reference/oal_wamd.md);
it is never delegated to `cobalt` (formula fidelity), while the
*display* balance table (`$balance`, `bal.tab()`) uses cobalt’s native
conventions, exactly as in `psAve`.

Known issues of the criterion, stated plainly: wAMD is a balance
*heuristic* — there is no theorem that it lands $`\lambda_n`$ inside the
oracle window; it reuses the same $`|\hat\beta_j|`$ as the penalty and
so inherits the same outcome-model-misspecification vulnerability; and
it is unstable under near-positivity violations, where a few huge
weights dominate the sums — the failure mode the rejoinder identifies in
correlated designs, and the motivation for GOAL.
[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
therefore warn when more than 5% of the scores sit at a clipping bound.

### Step 5: the propensity score

With the default `refit = FALSE`, the score is $`\operatorname{expit}`$
of the *penalized* coefficients (S&E/rejoinder). With `refit = TRUE`
(the Schnitzer variant, adapted from their longitudinal setting and
provenance-labeled as such), an unpenalized logistic regression is refit
on the selected covariates at every grid point, and both the wAMD and
the final score come from the refit; per-grid-point convergence is
recorded in `path$refit.converged`. All scores are clipped to `clip`
before the wAMD and in the returned `$ps`.

## GOAL: the generalized outcome-adaptive lasso

`method = "goal"` (Baldé, Yang & Lefebvre 2023) adds a ridge term to the
objective:

``` math
\hat\alpha(\lambda_n, \gamma, \lambda_2) \;=\; \arg\min_\alpha \;
-\ell_n(\alpha) \;+\; \lambda_n \sum_j \hat w_j |\alpha_j| \;+\; \lambda_2 \sum_j \alpha_j^2 .
```

The $`\ell_2`$ term supplies an elastic-net grouping effect and
numerical stability, so correlated confounders are selected together
rather than arbitrarily, and it behaves like propensity score truncation
under near-positivity — the rejoinder’s own remedy for OAL’s failure
mode. Baldé (2025) proves the oracle property under the OAL rate
conditions plus $`\lambda_2/\sqrt n \to 0`$.

The solver uses the Zou–Hastie data-augmentation trick, exactly as in
the rejoinder code and as described by Baldé (2025): for
$`\lambda_2 > 0`$, augment

``` math
X^s_{\text{aug}} = \begin{pmatrix} X^s \\ \sqrt{\lambda_2}\, I_p \end{pmatrix},
\qquad
A_{\text{aug}} = (A^\top, 0_p^\top)^\top,
```

fit the lasso on the augmented data, then rescale **all** coefficients —
*including* the intercept — by $`(1 + \lambda_2)`$ (Baldé’s official
code: `expit(cbind(1, X) %*% (1 + lambda2) * coef)`); scores are
computed from the original $`n`$ rows only. On augmented grid points the
raw penalty constant uses the **augmented** sample size,
$`\lambda_n = (n + q)^\delta`$ with $`q = p`$ augmentation rows (the
author’s `adaptive.lasso(lambda = n.q^(il))`, `n.q = n + q`). Grid
points with $`\lambda_2 = 0`$ run the *plain OAL code path* — no
augmentation, $`\lambda_n = n^\delta`$ — because even zero-valued
augmented rows would contribute intercept terms to the binomial
likelihood and break exact nesting. (Baldé’s script literally augments
with $`q`$ zero rows and keeps the $`(n+q)^\delta`$ base even at
$`\lambda_2 = 0`$; `oalasso` deliberately prefers exact nesting there.)
`lambda2` is selected jointly with $`(\delta, \gamma)`$ by the same
wAMD: the flat first-minimum over the $`\lambda_2`$-outer
$`\times`$$`\delta`$-inner grid, equivalent to Baldé’s nested rule
(per-$`\lambda_2`$ minima over $`(\delta, \gamma)`$, then the first
minimum over $`\lambda_2`$).

**Provenance.** The default grid
`lambda2 = c(0, 10^c(-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1))`
is the **author’s published grid** (Baldé 2025 supplement, “taken from
Zou and Hastie (2005)”), verified against the official lqa-based GOAL
code on 2026-07-02; the printed provenance line says so. Supplying your
own `lambda2` switches the label to “user-specified lambda2 grid”.

**Clipping.** Baldé’s reference code runs *unclipped* IPW weights
($`1/e`$ and $`1/(1-e)`$ with no truncation). The `oalasso` default
`clip = c(0.01, 0.99)` is a Shortreed–Ertefaie-lineage safety choice
(the suite constant shared with `psAve`); to reproduce the reference
behavior, effectively disable clipping with
`clip = c(1e-12, 1 - 1e-12)`.

## Standardization protocol (fixed, not an argument)

$`X^s = \texttt{scale}(X, \texttt{TRUE}, \texttt{TRUE})`$ is computed
**once**, with training statistics, over **all** columns including
dummies (the rejoinder convention). The outcome model and every `glmnet`
fit use the *same* $`X^s`$ with `standardize = FALSE`. This matters
because adaptive penalty weights are scale-dependent: letting glmnet
standardize internally, or scaling different subsets separately, changes
$`|\hat\beta_j|^{-\gamma}`$ and hence the selections. Output
coefficients are back-transformed ($`\alpha_j = \alpha^s_j / s_j`$;
intercept adjusted by $`-\sum_j \alpha^s_j\, c_j / s_j`$), and both
scales are returned (`coefficients`, `coefficients.std`); the centers
and scales are stored in `info`.

## The exact penalty-scale correction (why `glmnet` alone is not enough)

This is the core correctness requirement of the package, and the reason
a naive `glmnet(penalty.factor = ...)` translation of OAL does **not**
reproduce the published method.

Two scale mismatches stand between glmnet and the OAL objective:

1.  **Per-observation vs total loss.** glmnet minimizes
    $`-\tfrac{1}{n}\ell_n(\alpha) + \lambda \sum_j \tilde{v}_j |\alpha_j|`$
    — the *average* log-likelihood — while OAL’s
    $`\lambda_n = n^\delta`$ grid is defined against the *total*
    log-likelihood. Ignoring this changes the effective penalty by a
    factor of $`n`$, and the $`n^\delta`$ grid loses its theoretical
    meaning.
2.  **Penalty-factor rescaling.** glmnet internally rescales the
    supplied penalty factors to sum to `nvars`:
    $`\tilde v_j = v_j / \overline{v}`$ where
    $`\overline v = \tfrac1p \sum_j v_j`$. The documented consequence
    (“the penalty factors are internally rescaled to sum to nvars, and
    the lambda sequence will reflect this change”) is that the lambda
    you pass is silently reinterpreted. With adaptive weights this is
    not benign: one near-zero $`\hat\beta_j`$ inflates $`\overline v`$
    enormously and thereby *relaxes* the penalty on every other
    covariate.

Both are undone at once by evaluating each grid point at

``` math
s \;=\; \overline{v}_{\text{eff}} \cdot \frac{\lambda_n}{n},
\qquad
\overline{v}_{\text{eff}} = \frac{\sum_{j:\,v_j < \infty} v_j \;+\; \#\{j : v_j = \infty\}}{p},
```

via
`coef(fit, s = s, exact = TRUE, x = Xs, y = A, penalty.factor = pen)` —
`exact = TRUE` re-solves at exactly this $`s`$ rather than interpolating
along the path, and requires re-supplying the data and penalty factors.
The effective mean $`\overline{v}_{\text{eff}}`$ reflects glmnet’s
actual behavior with infinite factors: an infinite $`v_j`$
($`\hat b_j = 0`$ under the user hook below) is converted to `exclude` —
a hard drop — and then *internally reset to 1 before the sum-to-`nvars`
rescaling*, so each excluded column contributes 1 (not nothing) to the
rescaling mean on the full matrix. (The plain
$`\operatorname{mean}\{v_j : v_j < \infty\}`$ would be correct only on a
reduced matrix from which excluded columns are physically removed;
`oalasso` keeps the full matrix and uses $`\overline{v}_{\text{eff}}`$,
which reduces exactly to $`\operatorname{mean}(v)`$ when all factors are
finite.) Under GOAL augmentation both constants switch to the augmented
sample size: $`\lambda_n = (n + q)^\delta`$ and the $`1/n`$ becomes the
augmented per-observation weight $`1/(n + q)`$, $`q = p`$ augmentation
rows.

The correction is **always on** — it is not a knob — and its acceptance
criterion is not the formula but the *Karush–Kuhn–Tucker conditions of
the target objective*: the test suite verifies, at the exact solution of
every grid-point type, that
$`|\partial \ell / \partial \alpha_j| \le \lambda_n \hat w_j`$ for zero
coefficients with equality for nonzero ones. This recipe and its KKT
verification are ported from the rejoinder’s companion code
(`jmiahjones/OAL-rejoinder-2022`, `tests/test-glmnet.R`); a drift-guard
test pins the behavior against the installed glmnet version.

A minimal illustration of the machinery on simulated data (`x5`, `x6`
are instruments):

``` r

library(oalasso)

set.seed(20260702)
n <- 300
X <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, paste0("x", 1:6)))
A <- rbinom(n, 1, plogis(X %*% c(1, 1, 0, 0, 1, 1)))
Y <- as.numeric(2 * A + X %*% c(0.6, 0.6, 0.6, 0.6, 0, 0) + rnorm(n))
sim <- data.frame(A, Y, X)

fit <- oal(A ~ x1 + x2 + x3 + x4 + x5 + x6, data = sim, outcome = ~ Y)

## the path records, per grid point: lambda2, delta, lambda.n = n^delta,
## the paired gamma = 2*(2 - delta + 1), the exact glmnet s actually used,
## the wAMD, and the selected-set size
fit$path[, c("delta", "lambda.n", "gamma", "s", "wamd", "n.selected")]
#>    delta     lambda.n gamma            s       wamd n.selected
#> 1 -10.00 1.693509e-25 26.00 9.012153e+19 0.07692437          4
#> 2  -5.00 4.115226e-13 16.00 8.995054e+13 0.07691585          4
#> 3  -2.00 1.111111e-05 10.00 2.386841e+10 0.07691109          4
#> 4  -1.00 3.333333e-03  8.00 1.553156e+09 0.08201828          4
#> 5  -0.75 1.387264e-02  7.50 7.853034e+08 0.10329190          4
#> 6  -0.50 5.773503e-02  7.00 3.972424e+08 0.16404414          3
#> 7  -0.25 2.402811e-01  6.50 2.010352e+08 0.32596912          3
#> 8   0.25 4.161791e+00  5.50 5.155981e+07 0.70276946          0
#> 9   0.49 1.636023e+01  5.02 2.684984e+07 0.70276946          0
all.equal(fit$path$gamma, 2 * (2 - fit$path$delta + 1))
#> [1] TRUE
```

## Provenance of every default

| Setting | Default | Source |
|----|----|----|
| Penalty objective, adaptive weights $`\hat w_j = \lvert\hat\beta_j\rvert^{-\gamma}`$ | — | S&E (2017) |
| Outcome model $`Y \sim A + X`$, full sample, OLS/GLM | `family = gaussian()` | S&E (2017) / rejoinder code |
| $`\lambda`$ grid (exponents $`\delta`$) | `c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)` | S&E (2017) |
| $`\gamma`$ pairing $`2(g_{cf} - \delta + 1)`$, $`g_{cf} = 2`$ | `gamma = NULL`, `gamma.factor = 2` | S&E (2017) / rejoinder code |
| Fixed scalar $`\gamma`$ crossed with the $`\delta`$ grid | opt-in (`gamma = 2.5`) | Schnitzer et al. (2025) |
| Tuning criterion: wAMD with $`\lvert\hat\beta_j\rvert`$ weights | the only criterion | S&E (2017) |
| Estimand for the wAMD weights | `"ATE"` | S&E (2017) (note: `psave()` defaults to ATT) |
| PS from the penalized fit | `refit = FALSE` | S&E (2017) / rejoinder code |
| Post-selection unpenalized refit | opt-in (`refit = TRUE`) | Schnitzer et al. (2025), adapted from a longitudinal setting |
| GOAL $`\ell_2`$ term, augmentation + $`(1+\lambda_2)`$ rescale (all coefficients incl. intercept) | `method = "goal"` | Baldé et al. (2023); Baldé (2025) supplement |
| GOAL `lambda2` grid | `c(0, 10^c(-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1))` | **Author’s published grid** (Baldé 2025 supplement, “taken from Zou and Hastie (2005)”; verified against the official code 2026-07-02) |
| GOAL penalty base $`(n+q)^\delta`$ on augmented grid points | fixed | Baldé (2025) supplement (`n.q = n + q`) |
| Standardize once, `glmnet(standardize = FALSE)` | fixed protocol | Rejoinder code convention |
| Exact penalty-scale correction $`s = \overline v_{\text{active}} \lambda_n / n`$ | always on | Rejoinder code (KKT-verified) |
| Clipping constants | `clip = c(0.01, 0.99)` | Package choice, S&E lineage (suite constant shared with `psAve`; same constants in the legacy code). **Baldé’s reference runs unclipped** — use `clip = c(1e-12, 1 - 1e-12)` to reproduce it |
| Ties: first minimum, $`\lambda_2`$ then $`\delta`$ ascending | fixed | Package choice (consistent with Schnitzer et al.’s smallest-lambda tie note and psAve’s tolerant first-minimum) |
| `NA` anywhere in used variables is an error | fixed | Package/suite policy (never silent dropping) |

## Differences from the author’s legacy code

`oalasso` supersedes the OAL implementations (`fAL`, `fOAL`) in the
author’s earlier simulation code, which were built on the archived `lqa`
package. Four defects of that code are deliberately not reproduced, so
outputs *will* differ where those code paths were exercised:

1.  **The `fAL` gamma formula.** `fAL` computed
    `gamma = 2*(gcf - delta) + 1`, which disagrees with its own comment
    and gives $`\gamma`$ smaller by 1 than the S&E convention — a
    transcription bug. Only the correct pairing
    $`\gamma = 2(g_{cf} - \delta + 1)`$ (as in `fOAL` and the rejoinder)
    is implemented; the buggy formula appears nowhere in the package,
    and a regression test asserts its absence from the grid.
2.  **The `varlist` wAMD bug.** One `fOAL` branch passed the *treatment*
    variable name as the covariate list to its wAMD function, so the
    criterion was computed over the treatment column alone and the
    lambda selection in that branch was meaningless.
    [`oal_wamd()`](https://kabajiro.github.io/oalasso/reference/oal_wamd.md)
    always sums over all $`p`$ covariate columns.
3.  **Global-`n` scoping.** In the legacy code, the `n` in
    `lambda = n^delta` was never defined inside the functions — it
    resolved from the calling (global) environment, silently using
    whatever `n` the surrounding script happened to hold. In `oalasso`,
    $`n`$ is always `nrow(model.matrix)`, by construction.
4.  **Per-subset [`scale()`](https://rdrr.io/r/base/scale.html).** The
    legacy code standardized fitting and evaluation subsets separately,
    by different means and SDs — which changes the adaptive weights and
    the scores. The fixed protocol above (scale once, training
    statistics) makes this unreproducible by design.

## Relation to other software

- **`lqa`** (Ulbricht 2010) — the package S&E’s original implementation
  was built on — is archived and removed from CRAN; `oalasso` does not
  depend on or reference it as a dependency. The `glmnet` route with the
  exact-scale correction replaces it.
- **`jmiahjones/OAL-rejoinder-2022`** (GitHub) — the rejoinder’s
  companion code — is the reference implementation: its `oal_funs.R`
  recipe (standardization, pairing, $`s`$-correction, GOAL augmentation)
  and its KKT test are what `oalasso` reproduces and tests against.
- **Schnitzer-Biostats-Lab/Longitudinal-outcome-adaptive-LASSO**
  (GitHub; Schnitzer et al. 2025, <doi:10.1002/sim.70316>) — the
  longitudinal OAL pipeline from which the `refit = TRUE` variant and
  the fixed-$`\gamma`$ mode are adapted. Their SE-scaled balance
  criterion is *not* implemented in v1 (its validation is
  longitudinal-only, and its bootstrap would break the package’s
  determinism).
- **GOAL disambiguation.** Baldé’s GOAL has no public repository; the
  official code lives in the Wiley supplement of Baldé et al. (2023).
  The GitHub repository `QianGao-SXMU/GOAL` (Gao et al.) is a
  *different* method — a generalized propensity score for continuous
  exposures — that happens to share the acronym. `method = "goal"` here
  is Baldé’s.
- **`ctmle`** (CRAN) implements collaborative TMLE (van der Laan &
  Gruber 2010, <doi:10.2202/1557-4679.1181>; scalable versions Ju et
  al. 2019, <doi:10.1177/0962280217729845>) — the target-driven,
  learner-agnostic parent framework of outcome-aware nuisance selection.
- **`hal9001`** (CRAN) implements the highly adaptive lasso, on which
  the outcome-highly-adaptive lasso (oHAL; Ju, Benkeser & van der Laan
  2020, <doi:10.1111/biom.13121>) is built — see the next section.
- A Python OAL implementation exists at
  `tom-beer/Outcome-Adaptive-LASSO` (GitHub).

## A verified scale difference against the historical `lqa`-based references

During development we ran `oalasso` head-to-head against the *official*
GOAL/OAL reference code from the supplement of Baldé (2025) on identical
data, with the archived `lqa` engine installed in an isolated library
(the full report ships with the package sources under `verification/`).
The outcome-model coefficients agree to machine precision, but the
selected $`(\lambda, \gamma)`$ can differ. The root cause is
reproducible and instructive: `lqa.default()` is called with its default
`standardize = TRUE`, which rescales the (already SD-standardized)
columns to unit *Euclidean norm* internally. On SD-standardized data
that norm is exactly $`\sqrt{n-1}`$, so **for plain OAL** the reference
code effectively solves the published objective with the L1 constant
inflated to $`\lambda_n\sqrt{n-1}`$. Bridging `glmnet` with that factor
reproduces the `lqa` OAL grid points and the reference’s exact OAL
selection (the residual gap traces to `lqa`’s loose `conv.eps = 1e-3`;
tightening it collapses the difference to $`\sim 10^{-6}`$). **For GOAL
with $`\lambda_2 > 0`$ the story is messier and a single scale factor
does not suffice:** `lqa` standardizes the *augmented* design (realized
column norms $`\sqrt{n - 1 + \lambda_2(1 - 1/(n+q))}`$, with nonzero
augmented-row column means interacting with the intercept), and its
modified PIRLS treats the augmentation rows as weight-1/response-0
*ridge* terms inside each IRLS step — structurally different from the
binomial pseudo-row objective a `glmnet` augmentation solves. GOAL grid
points therefore agree only approximately with the reference solver even
after bridging; selections on our test data agreed at $`\lambda_2 = 0`$
and differed in PS by $`\sim 4\times 10^{-2}`$ at the selected point
overall (full details in the equivalence report).

`oalasso` deliberately solves the objective *as printed in the papers* —
that is what the KKT acceptance tests pin down — rather than reproducing
the historical solver’s scale. Because $`\lambda_n = n^{\delta}`$ and
$`\sqrt{n-1}\approx n^{1/2}`$, users who want selections on the
historical `lqa` scale can simply shift the exponent grid:

``` r

oal(A ~ ., data = dat, outcome = ~ y,
    lambda = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49) + 0.5)
```

In practice the wAMD criterion re-selects the tuning adaptively, so
downstream propensity scores are typically close; but simulation
comparisons against numbers published with `lqa`-based code should keep
this scale difference in mind.

## Beyond linearity: honest notes on nonlinear extensions

### The linear limitation

Everything above — the penalty weights, the wAMD weights, and both
published oracle theorems — rests on the *linear* outcome model. Under
misspecification, $`\hat\beta`$ converges to best-linear-projection
coefficients, and the set it zeroes need not be the outcome-unrelated
set: a confounder whose outcome effect is purely nonlinear can have a
zero projection coefficient and be wrongly dropped, while a null
covariate correlated with a strong predictor can be wrongly kept.
Selection consistency for the confounder set is therefore not guaranteed
under misspecification, and to our knowledge no dedicated published
study isolates OAL’s selection behavior in that regime (a literature
reading, not a theorem); the rejoinder explicitly used a correctly
specified linear model and pointed to flexible cross-fitted methods as
the way forward. Hand-coding basis terms (splines, interactions) into
`formula` and `outcome` is the within-scope remedy, since both models
then share the expanded columns.

For genuinely nonparametric outcome-adaptive selection, use the
published tools: the **outcome-highly-adaptive lasso** (Ju, Benkeser &
van der Laan 2020, <doi:10.1111/biom.13121>), which selects over HAL’s
spline-interaction basis via `hal9001` and pairs with TMLE, and the
**C-TMLE** framework (`ctmle`; van der Laan & Gruber 2010), which
selects nuisance models by their value for the target parameter with any
learner. Those methods carry their own theory; `oalasso` does not
attempt to reimplement them.

### The `outcome.coef` hook

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md)’s single
extension seam is the `outcome.coef` argument: a named numeric vector
$`b_j`$ over the standardized model-matrix columns that replaces the
internal outcome model in *both* the penalty
($`\hat w_j = |b_j|^{-\gamma}`$) and the wAMD weights. **The values are
interpreted on the standardized-**$`X^s`$\*\* scale\*\* — adaptive
weights are scale-dependent, so coefficients from an unstandardized fit
will not mean what you think. Zeros are allowed here (unlike in the
internal model): $`|0|^{-\gamma} = \infty`$ becomes a hard `exclude` in
glmnet. Using the hook switches the printed provenance to an explicit
experimental label — screening use only.

The reference construction (reproducing the internal weights by hand)
and the ridge variant for $`p`$ close to $`n`$ (the natural analogue of
the author’s legacy `shrink = TRUE` branch):

``` r

## Reference construction: what oal() computes internally
f  <- A ~ x1 + x2 + x3 + x4 + x5 + x6
X  <- model.matrix(f, sim)[, -1, drop = FALSE]
Xs <- scale(X, TRUE, TRUE)                     # the package's fixed protocol
b  <- coef(lm(sim$Y ~ sim$A + Xs))[-(1:2)]     # drop intercept and treatment
names(b) <- colnames(Xs)
fit.hook <- oal(f, data = sim, outcome.coef = b)   # same selections as outcome = ~ Y
                                                   # (provenance label differs)

## Ridge outcome model for p close to n
set.seed(1)                                    # cv.glmnet folds are random
cv <- glmnet::cv.glmnet(x = cbind(A = sim$A, Xs), y = sim$Y, alpha = 0,
                        penalty.factor = c(0, rep(1, ncol(Xs))))
b  <- as.numeric(coef(cv, s = "lambda.min"))[-(1:2)]
names(b) <- colnames(Xs)
fit.ridge <- oal(f, data = sim, outcome.coef = b)
```

### Experimental: importance-based weights (a recipe, not a feature)

A tempting idea is to make the weights nonlinear-aware by replacing
$`|\hat\beta_j|`$ with a tree-ensemble variable importance. The package
deliberately does **not** export this: no published, peer-reviewed
method plugs tree importance into an OAL-style adaptive penalty, and the
theory does not transfer. The following framing, reproduced verbatim
from the package’s theory review, is the ceiling of what may be claimed
for such a construction:

> 1.  Claim it as a **screening / dimension-reduction** procedure
>     targeting outcome-relevant covariates (X_C ∪ X_P), not an oracle
>     selector.
> 2.  State the target guarantee as a **sure screening property** (Fan &
>     Lv 2008) — no false negatives w.p.→1 — with a second-stage refit
>     for false positives.
> 3.  Do **not** claim the OAL/GOAL oracle property (√n asymptotic
>     normality with oracle variance): tree importance is not a
>     √n-consistent signed coefficient, so the Zou/Shortreed–Ertefaie
>     proof does not transfer.
> 4.  Do **not** claim λ_n rate conditions (λ_n/√n→0, λ_n n^{γ/2−1}→∞)
>     confer selection consistency here — those depend on
>     \|β̃ⱼ\|=O_p(n^{−1/2}) under the null, which importance measures do
>     not satisfy.
> 5.  Explicitly warn that **impurity/permutation importance is biased**
>     with mixed-scale (Strobl 2007) and correlated (Strobl 2008)
>     covariates, which can let an instrument masquerade as relevant.
> 6.  **Default to conditional permutation importance and
>     unbiased/honest trees** to mitigate those biases; document this as
>     a requirement, not an option.
> 7.  Offer **model-X knockoffs** (Candès et al. 2018) as the
>     finite-sample **FDR-controlled** confounder-selection mode,
>     flagging its dependence on an accurate covariate model.
> 8.  Keep the **instrument-exclusion rationale** (Brookhart 2006; Myers
>     2011; Pearl 2011): beneficial under possible unmeasured
>     confounding and IPW estimation; note nonparametric outcome
>     modeling does not remove Z-bias or positivity concerns.
> 9.  Recommend pairing selection with a **doubly robust estimator
>     (AIPW/TMLE) + cross-fitting**, so misspecification/selection error
>     degrades gracefully rather than breaking a fragile oracle claim.
> 10. Summarize honestly: “empirically competitive,
>     screening-consistent, FDR-controllable confounder selection for
>     nonlinear settings” — **not** “oracle-consistent causal variable
>     selection.”

Under that framing, the recipe below feeds *conditional* permutation
importance from a conditional-inference forest (`party::cforest` with
`cforest_unbiased()` — the Strobl-recommended configuration; see points
5–6 above) into the `outcome.coef` hook. It is intentionally left
unevaluated: run it only with the framing above in mind, and pair the
result with a doubly robust downstream estimator (point 9).

``` r

## EXPERIMENTAL — screening use only, per the ten points above. Requires `party`.
library(party)

f  <- A ~ x1 + x2 + x3 + x4 + x5 + x6
X  <- model.matrix(f, sim)[, -1, drop = FALSE]
Xs <- scale(X, TRUE, TRUE)                 # same standardization as oal()
df <- data.frame(Y = sim$Y, Xs, check.names = FALSE)

set.seed(20260702)                         # forests are stochastic
cf <- cforest(reformulate(colnames(Xs), response = "Y"), data = df,
              controls = cforest_unbiased(ntree = 1000,
                                          mtry  = ceiling(sqrt(ncol(Xs)))))
vi <- varimp(cf, conditional = TRUE)       # conditional permutation importance
vi <- pmax(vi[colnames(Xs)], 0)            # negative importance -> 0 = hard exclusion

fit.imp <- oal(f, data = sim, outcome.coef = vi)
fit.imp    # provenance line flags the user-supplied weights; screening use only
```

Notes on the recipe. The names of the vector must exactly cover the
standardized model-matrix columns (factor covariates must be expanded
via `model.matrix` first, as above). Importance is nonnegative, which is
fine: only $`|b_j|`$ enters the penalty and the wAMD. Zeroed importances
become hard exclusions. And the choice of importance measure is
load-bearing, not cosmetic — plain impurity or marginal permutation
importance can hand an instrument a small penalty precisely when it is
correlated with a confounder, defeating the method’s purpose (Strobl et
al. 2007, 2008).

## References

Baldé, I. (2025). The oracle property of the generalized
outcome-adaptive lasso. *Statistics & Probability Letters*, 221, 110379.
<doi:10.1016/j.spl.2025.110379>

Baldé, I., Yang, Y. A., & Lefebvre, G. (2023). Reader reaction to
“Outcome-adaptive lasso: Variable selection for causal inference” by
Shortreed and Ertefaie (2017). *Biometrics*, 79(1), 514–520.
<doi:10.1111/biom.13683>

Brookhart, M. A., Schneeweiss, S., Rothman, K. J., Glynn, R. J., Avorn,
J., & Stürmer, T. (2006). Variable selection for propensity score
models. *American Journal of Epidemiology*, 163(12), 1149–1156.
<doi:10.1093/aje/kwj149>

Candès, E., Fan, Y., Janson, L., & Lv, J. (2018). Panning for gold:
‘model-X’ knockoffs for high dimensional controlled variable selection.
*Journal of the Royal Statistical Society: Series B*, 80(3), 551–577.
<doi:10.1111/rssb.12265>

Fan, J., & Lv, J. (2008). Sure independence screening for ultrahigh
dimensional feature space. *Journal of the Royal Statistical Society:
Series B*, 70(5), 849–911. <doi:10.1111/j.1467-9868.2008.00674.x>

Jones, J., Ertefaie, A., & Shortreed, S. M. (2023). Rejoinder to “Reader
reaction to ‘Outcome-adaptive lasso: Variable selection for causal
inference’ by Shortreed and Ertefaie (2017)”. *Biometrics*, 79(1),
521–525. <doi:10.1111/biom.13681>

Ju, C., Benkeser, D., & van der Laan, M. J. (2020). Robust inference on
the average treatment effect using the outcome highly adaptive lasso.
*Biometrics*, 76(1), 109–118. <doi:10.1111/biom.13121>

Ju, C., Gruber, S., Lendle, S. D., Chambaz, A., Franklin, J. M., Wyss,
R., Schneeweiss, S., & van der Laan, M. J. (2019). Scalable
collaborative targeted learning for high-dimensional data. *Statistical
Methods in Medical Research*, 28(2), 532–554.
<doi:10.1177/0962280217729845>

Myers, J. A., Rassen, J. A., Gagne, J. J., Huybrechts, K. F.,
Schneeweiss, S., Rothman, K. J., Joffe, M. M., & Glynn, R. J. (2011).
Effects of adjusting for instrumental variables on bias and precision of
effect estimates. *American Journal of Epidemiology*, 174(11),
1213–1222. <doi:10.1093/aje/kwr364>

Pearl, J. (2011). Invited commentary: Understanding bias amplification.
*American Journal of Epidemiology*, 174(11), 1223–1227.
<doi:10.1093/aje/kwr352>

Schnitzer, M. E., Talbot, D., Liu, Y., Berger, C., Wang, G., O’Loughlin,
J., Sylvestre, M.-P., & Ertefaie, A. (2025). Adaptive sparsening and
smoothing of the treatment model for longitudinal causal inference using
outcome-adaptive LASSO and marginal fused LASSO. *Statistics in
Medicine*. <doi:10.1002/sim.70316>

Shortreed, S. M., & Ertefaie, A. (2017). Outcome-adaptive lasso:
Variable selection for causal inference. *Biometrics*, 73(4), 1111–1122.
<doi:10.1111/biom.12679>

Strobl, C., Boulesteix, A.-L., Zeileis, A., & Hothorn, T. (2007). Bias
in random forest variable importance measures: Illustrations, sources
and a solution. *BMC Bioinformatics*, 8, 25.
<doi:10.1186/1471-2105-8-25>

Strobl, C., Boulesteix, A.-L., Kneib, T., Augustin, T., & Zeileis, A.
(2008). Conditional variable importance for random forests. *BMC
Bioinformatics*, 9, 307. <doi:10.1186/1471-2105-9-307>

van der Laan, M. J., & Gruber, S. (2010). Collaborative double robust
targeted maximum likelihood estimation. *The International Journal of
Biostatistics*, 6(1), Article 17. <doi:10.2202/1557-4679.1181>

Zou, H. (2006). The adaptive lasso and its oracle properties. *Journal
of the American Statistical Association*, 101(476), 1418–1429.
<doi:10.1198/016214506000000735>
