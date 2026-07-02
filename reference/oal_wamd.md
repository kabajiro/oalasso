# Weighted absolute mean difference (wAMD) of Shortreed and Ertefaie (2017)

Scores an arbitrary candidate propensity score on the exact wAMD balance
criterion used by
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) to select
its tuning parameters. This function is the package's **single source of
truth** for the criterion: the internal grid search calls it verbatim,
and exporting it lets any candidate score – from any package – be scored
on the same yardstick (the mirror of
[`psAve::psave_criteria()`](https://kabajiro.github.io/psAve/reference/psave_criteria.html)).

## Usage

``` r
oal_wamd(
  ps,
  treat,
  covs,
  coef,
  estimand = c("ATE", "ATT"),
  clip = c(0.01, 0.99)
)
```

## Arguments

- ps:

  Numeric vector of propensity scores, strictly inside (0, 1) before
  clipping.

- treat:

  Treatment vector: numeric 0/1, logical, or a two-level
  factor/character (second level = treated), as in
  [`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md).

- covs:

  Numeric matrix of covariates, one row per unit; no missing values.

- coef:

  Numeric vector of outcome-model coefficients, one per column of `covs`
  (matched by name when both are named, else by position). Absolute
  values are used as the balance weights.

- estimand:

  `"ATE"` (default) or `"ATT"`; selects the IPW weight formula above.

- clip:

  Length-2 numeric; the scores are clipped to `[clip[1], clip[2]]`
  **before** the weights are formed (default `c(0.01, 0.99)`).

## Value

A list with elements

- `total`:

  the wAMD, \\\sum_j \|\beta_j\| d_j\\;

- `by.covariate`:

  named numeric vector of the per-covariate contributions \\\|\beta_j\|
  d_j\\ (they sum to `total`).

## Details

The propensity scores are first clipped to `clip`; the
inverse-probability weights at `estimand` are then (with \\e_i\\ the
clipped score) \$\$ATE:\\ w_i = A_i/e_i + (1 - A_i)/(1 - e_i), \qquad
ATT:\\ w_i = A_i + (1 - A_i)\\ e_i/(1 - e_i).\$\$ For each covariate
column \\j\\ the weighted absolute mean difference is \$\$d_j = \left\|
\frac{\sum_i w_i X\_{ij} A_i}{\sum_i w_i A_i} - \frac{\sum_i w_i X\_{ij}
(1 - A_i)}{\sum_i w_i (1 - A_i)} \right\|,\$\$ and the criterion is
\\\mathrm{wAMD} = \sum_j \|\beta_j\|\\ d_j\\, summed over **all**
columns of `covs` (including any the propensity score model excluded),
with \\\beta_j\\ the outcome-model coefficients supplied in `coef`
(their absolute values are taken internally).

Inside [`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
`covs` is the **standardized** model matrix and `coef` the outcome
coefficients on that same standardized scale (the wAMD, like the
adaptive weights, is scale-dependent). The formula is implemented
natively – never delegated to cobalt – so that it reproduces the
published criterion exactly.

## References

Shortreed SM, Ertefaie A (2017). Outcome-adaptive lasso: variable
selection for causal inference. *Biometrics*, 73(4), 1111-1122.
[doi:10.1111/biom.12679](https://doi.org/10.1111/biom.12679)

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`weights.oal()`](https://kabajiro.github.io/oalasso/reference/weights.oal.md)

## Examples

``` r
set.seed(7)
X <- matrix(rnorm(200), 100, 2, dimnames = list(NULL, c("x1", "x2")))
A <- rbinom(100, 1, plogis(X[, 1]))
ps <- plogis(0.8 * X[, 1])
oal_wamd(ps, A, X, coef = c(x1 = 0.5, x2 = 0.1), estimand = "ATE")
#> $total
#> [1] 0.2089005
#> 
#> $by.covariate
#>          x1          x2 
#> 0.207586166 0.001314295 
#> 
```
