# Balance tables for oal objects

A method for
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html):
assesses balance on the (original-scale) covariates of an
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) fit under
the inverse-probability weights implied by the outcome-adaptive lasso
propensity score, which is supplied as a `distance` measure.

## Usage

``` r
# S3 method for class 'oal'
bal.tab(x, ...)
```

## Arguments

- x:

  An `oal` object.

- ...:

  Further arguments passed on to
  [`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
  (e.g., `un = TRUE`, `thresholds = c(m = 0.1)`).

## Value

A `bal.tab` object; see
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html).

## Details

The call delegates to the default cobalt machinery as
`cobalt::bal.tab(<covariates>, treat = x$treat, weights = x$weights, s.d.denom = <by estimand>, distance = data.frame(ps = x$ps), ...)`,
so all the usual cobalt arguments (`un`, `stats`, `thresholds`, ...) are
available and display conventions are cobalt's own. The *selection
criterion* inside
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) is the
papers' wAMD instead, computed natively on the standardized covariates;
see
[`oal_wamd()`](https://kabajiro.github.io/oalasso/reference/oal_wamd.md).

## References

Shortreed SM, Ertefaie A (2017). Outcome-adaptive lasso: variable
selection for causal inference. *Biometrics*, 73(4), 1111-1122.
[doi:10.1111/biom.12679](https://doi.org/10.1111/biom.12679)

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html),
[`plot.oal()`](https://kabajiro.github.io/oalasso/reference/plot.oal.md)

## Examples

``` r
data("lalonde", package = "MatchIt")
fit <- oal(treat ~ age + educ + married + re74, data = lalonde,
           outcome = ~ re78)
cobalt::bal.tab(fit, un = TRUE)
#> Balance Measures
#>             Type Diff.Un Diff.Adj
#> ps      Distance  0.8366  -0.1688
#> age      Contin. -0.2419  -0.1062
#> educ     Contin.  0.0448   0.1938
#> married   Binary -0.3236   0.0340
#> re74     Contin. -0.5958   0.5336
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    409.07   50.16
```
