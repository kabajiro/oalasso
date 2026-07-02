# oalasso: Outcome-Adaptive Lasso Propensity Scores

Estimates propensity scores by the outcome-adaptive lasso (OAL) of
Shortreed and Ertefaie (2017) and the generalized outcome-adaptive lasso
(GOAL) of Balde, Yang and Lefebvre (2023), using glmnet with an exact
penalty-scale correction so that the published objectives and tuning
grids are reproduced, and tuning by the papers' weighted absolute mean
difference (wAMD) balance criterion. The resulting score is designed to
be supplied directly to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
as a distance measure, to
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
as a propensity score, or to
[`psAve::psave()`](https://rdrr.io/pkg/psAve/man/psave.html) as an
appended candidate.

## Details

The single estimation function is
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md). Its
result hands off to the existing ecosystem:
[`oal_match()`](https://kabajiro.github.io/oalasso/reference/oal_match.md)
/
[`oal_weight()`](https://kabajiro.github.io/oalasso/reference/oal_weight.md)
(or the equivalent explicit
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
/
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html)
calls) and
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
(which has a method for `oal` objects).
[`oal_wamd()`](https://kabajiro.github.io/oalasso/reference/oal_wamd.md)
exposes the selection criterion for methods research and testing.

### Design notes

glmnet sits in `Imports` – the one principled departure from the suite's
engines-in-Suggests rule – because it is the sole solver of the
estimator's objective, and the exact-scale correction depends on its
penalty-factor rescaling, its infinite-penalty-to-exclusion conversion,
and its `coef(..., exact = TRUE)` re-supply semantics (floor
`>= 4.1-2`). The archived lqa package of the original reference
implementation is not a dependency. cobalt powers the display balance
table and
[`bal.tab.oal()`](https://kabajiro.github.io/oalasso/reference/bal.tab.oal.md);
the wAMD criterion itself is native code.

Note the name collision: the GOAL implemented here (`method = "goal"`)
is the generalized outcome-adaptive lasso of Balde, Yang and Lefebvre
(2023) for binary treatments, not the unrelated "GOAL" generalized
propensity score software for continuous exposures by Gao and
colleagues.

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

## See also

Useful links:

- <https://github.com/kabajiro/oalasso>

- Report bugs at <https://github.com/kabajiro/oalasso/issues>

## Author

**Maintainer**: Daijiro Kabata <daijiro.kabata@port.kobe-u.ac.jp>
\[copyright holder\]

Authors:

- Daijiro Kabata <daijiro.kabata@port.kobe-u.ac.jp> \[copyright holder\]
