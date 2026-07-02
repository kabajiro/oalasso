# Weight by an outcome-adaptive lasso propensity score

Convenience pass-through to
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html):
constructs balancing weights from the propensity score of an
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) fit,
reusing the stored formula and data. Equivalent to the canonical
explicit call

    WeightIt::weightit(<formula>, data = <data>, ps = fit$ps, estimand = <estimand>, ...)

All `...` arguments are forwarded verbatim; the return value is an
ordinary `weightit` object.

## Usage

``` r
oal_weight(object, estimand = object$estimand, ...)
```

## Arguments

- object:

  An `oal` object.

- estimand:

  The estimand passed to
  [`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html);
  defaults to the estimand of the fit (note that the wAMD *selection*
  already used the fitted estimand's weights).

- ...:

  Arguments forwarded verbatim to
  [`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html).

## Value

A `weightit` object; see
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html).

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`oal_match()`](https://kabajiro.github.io/oalasso/reference/oal_match.md),
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html),
[`WeightIt::get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.html)

## Examples

``` r
data("lalonde", package = "MatchIt")
fit <- oal(treat ~ age + educ + married + re74, data = lalonde,
           outcome = ~ re78)
w <- oal_weight(fit)
```
