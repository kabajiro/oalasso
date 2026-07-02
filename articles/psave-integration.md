# Using oalasso with psAve

## Two packages, two kinds of propensity score model uncertainty

`oalasso` and [psAve](https://github.com/kabajiro/psAve) are companion
packages that address *different* sources of uncertainty in a propensity
score (PS) analysis:

- **Which covariates?** The outcome-adaptive lasso (Shortreed &
  Ertefaie 2017) selects the covariates that belong in the PS model —
  confounders and outcome predictors — and excludes instruments and
  noise variables, whose inclusion inflates variance and can amplify
  unmeasured-confounding bias (Brookhart et al. 2006; Myers et
  al. 2011).
- **Which model form?**
  [`psAve::psave()`](https://rdrr.io/pkg/psAve/man/psave.html) (Kabata,
  Stuart & Shintani 2024) averages candidate PS models of different
  functional forms (logistic regression, CART, random forest, gradient
  boosting, …), with convex mixing weights selected by balance on the
  prognostic score.

The two compose cleanly because both packages speak the same minimal
contract: a propensity score is a plain numeric vector of length $`n`$,
named by `rownames(data)`, strictly inside (0, 1). `psave()` accepts
such a vector as an extra candidate via its `ps.append` argument, so the
OAL score can compete with — and be averaged with — the flexible
learners:

``` r

fit <- oal(treat ~ x1 + x2 + ..., data = dat, outcome = ~ y, estimand = "ATT")
ma  <- psAve::psave(treat ~ x1 + x2 + ..., data = dat, outcome = ~ y,
                    ps.append = cbind(oal = fit$ps))
```

Contract compliance is by construction:
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) errors on
missing data rather than dropping rows (so no reordering or shrinking
can misalign rows), its default `clip = c(0.01, 0.99)` equals
`psave()`’s default clip (so psave’s re-clipping of the appended column
is a no-op), and `fit$ps` carries the rownames of `data`, which
`psave()` checks. Passing the bare vector (`ps.append = fit$ps`) also
works; the candidate is then labeled `"append"` in psave’s output.
Wrapping it as a one-column matrix, `cbind(oal = fit$ps)`, keeps the
rownames for psave’s alignment check *and* gives the candidate a
readable label.

## When this composition helps

Use it when both kinds of uncertainty are live at once:

- **Suspected instruments.** Some covariates plausibly predict treatment
  assignment but not the outcome (e.g., prescriber preference, calendar
  or access variables). None of psave’s default candidate learners
  performs variable selection against instruments; the OAL candidate
  contributes exactly that.
- **Model-form uncertainty.** The covariate–treatment relationship may
  be nonlinear or involve interactions that OAL’s logistic-linear PS
  model misses; the flexible candidates contribute that. (OAL’s
  *selection* is also linear-in-covariates on the outcome side — see the
  method-details vignette for honest notes on this limitation.)

The averaging step then adjudicates empirically: if the instrument-free
linear score balances the prognostic score better, it earns weight; if
the flexible learners capture something the linear model misses, they
do.

## A worked example

We use the `lalonde` data. Two alignment choices deserve a comment.
First, [`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md)’s
default estimand is ATE (Shortreed & Ertefaie’s wAMD convention) while
`psave()`’s is ATT, so for a coherent pipeline set them explicitly to
the same estimand — here ATT. Second, give both functions the same
`outcome` variable: OAL uses it to build its penalty weights, psave to
build the prognostic score its criterion balances.

``` r

library(oalasso)
data("lalonde", package = "MatchIt")

fit <- oal(treat ~ age + educ + race + married + nodegree + re74 + re75,
           data = lalonde, outcome = ~ re78, estimand = "ATT")
fit$lambda            # the selected (delta, gamma) tuning point
#>         delta      lambda.n         gamma       lambda2 
#> -1.000000e+01  1.313156e-28  2.600000e+01            NA
```

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) is
deterministic; the seed below is for `psave()`, whose default learner
menu includes stochastic learners. To keep this vignette light we
restrict psave to its two deterministic learners — with the full default
menu (`ranger`, `xgboost`) the call is the same:

``` r

set.seed(1234)
ma <- psAve::psave(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   data = lalonde, outcome = ~ re78,
                   ps.methods   = c("glm", "rpart"),
                   prog.methods = c("glm", "rpart"),
                   ps.append    = cbind(oal = fit$ps))
ma
#> A psave object (model-averaged propensity score)
#>  - estimand:  ATT
#>  - criterion: prog (weighted ASMD of the model-averaged prognostic score)
#>  - sample:    614 units (185 treated, 429 control)
#> 
#> lambda (PS mixing weights):
#>   glm    0.000  |                    |
#>   rpart  0.000  |                    |
#>   oal    1.000  |====================|
#> 
#> gamma (prognostic mixing weights):
#>   glm    0.050  |=                   |
#>   rpart  0.950  |=================== |
#> 
#> Criterion value at selected lambda: 0.0387
#> 
#> Balance preview (worst covariates + prognostic score):
#>          smd.un smd.wt ks.un ks.wt
#> age       0.309  0.119 0.158 0.308
#> married   0.826  0.047 0.324 0.019
#> nodegree  0.245  0.041 0.111 0.018
#> prog      0.426  0.039 0.395 0.088
#> 
#> Next:
#>   MatchIt::matchit(treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde, distance = x$ps)
#>     or: psave_match(x)
#>   WeightIt::weightit(treat ~ age + educ + race + married + nodegree + re74 + re75, data = lalonde, ps = x$ps, estimand = "ATT")
#>     or: psave_weight(x)
```

The printed mixing weights $`\lambda`$ now include an `oal` column: the
share of the averaged score contributed by the outcome-adaptive lasso
candidate. A weight of 0 means the criterion preferred the base learners
(appended candidates are placed last on psave’s grid, so exact ties
favor the base candidates — a documented psAve rule); a large weight
means instrument exclusion paid off in prognostic-score balance.
Everything downstream is unchanged psAve usage: `psave_match(ma)`,
`psave_weight(ma)`, `cobalt::bal.tab(ma)`.

## The criterion-interaction caveat

One honest paragraph before you interpret the mixing weights. An OAL
candidate pairs *naturally* with psave’s default `criterion = "prog"`:
both are outcome-oriented — OAL leaves instruments unbalanced *on
purpose* because they cannot cause confounding, and the prognostic-score
criterion does not charge it for that, because instruments contribute
nothing to the prognostic score. The other criteria can mis-rank an
instrument-free candidate for reasons that have nothing to do with bias.
Covariate-balance criteria (`"smd"`, `"ks"`) treat every covariate
equally, so they penalize the OAL candidate for residual imbalance on
the very instruments it deliberately declined to balance — imbalance
that is harmless (instruments do not confound) and whose “repair” is
what Myers et al. (2011) warn amplifies unmeasured-confounding bias. The
prediction criterion (`"logloss"`) is worse still: instruments are, by
definition, strong treatment predictors, so a score that excludes them
*must* predict treatment less accurately — log loss systematically
rewards exactly the covariates OAL exists to remove. If you append an
OAL candidate, keep `criterion = "prog"` (the default); with `"smd"`,
`"ks"`, or `"logloss"`, a low weight on the `oal` column is not evidence
against instrument exclusion.

## Validation status

Label this composition honestly in your reporting: it is an **extension
beyond both packages’ validated sets**. The simulations of Kabata,
Stuart & Shintani (2024) did not include an outcome-adaptive-lasso
candidate among the averaged models, and Shortreed & Ertefaie (2017) did
not study model averaging over their score. Each component is
implemented exactly as published (each package’s provenance labels say
precisely what is validated and what is not), but the *combination* has
no dedicated published evaluation — treat it as a principled sensitivity
analysis rather than a validated named method, and cite it as such
(e.g., “an OAL propensity score (Shortreed & Ertefaie 2017), included as
a candidate in prognostic-score-based model averaging (Kabata et
al. 2024)”).

## References

Brookhart, M. A., Schneeweiss, S., Rothman, K. J., Glynn, R. J., Avorn,
J., & Stürmer, T. (2006). Variable selection for propensity score
models. *American Journal of Epidemiology*, 163(12), 1149–1156.
<doi:10.1093/aje/kwj149>

Kabata, D., Stuart, E. A., & Shintani, A. (2024). Prognostic score-based
model averaging approach for propensity score estimation. *BMC Medical
Research Methodology*, 24, 228. <doi:10.1186/s12874-024-02350-y>

Myers, J. A., Rassen, J. A., Gagne, J. J., Huybrechts, K. F.,
Schneeweiss, S., Rothman, K. J., Joffe, M. M., & Glynn, R. J. (2011).
Effects of adjusting for instrumental variables on bias and precision of
effect estimates. *American Journal of Epidemiology*, 174(11),
1213–1222. <doi:10.1093/aje/kwr364>

Shortreed, S. M., & Ertefaie, A. (2017). Outcome-adaptive lasso:
Variable selection for causal inference. *Biometrics*, 73(4), 1111–1122.
<doi:10.1111/biom.12679>
