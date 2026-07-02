# Match on an outcome-adaptive lasso propensity score

Convenience pass-through to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html):
matches on the propensity score of an
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) fit,
reusing the formula and data stored in the object. Equivalent to the
canonical explicit call

    MatchIt::matchit(<formula>, data = <data>, distance = fit$ps, ...)

but with no opportunity for row misalignment between the two steps. All
`...` arguments are forwarded verbatim; the return value is an ordinary
`matchit` object, so the full MatchIt/cobalt toolkit applies.

## Usage

``` r
oal_match(object, ...)
```

## Arguments

- object:

  An `oal` object.

- ...:

  Arguments forwarded verbatim to
  [`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
  (e.g., `method`, `caliper`, `ratio`, `replace`).

## Value

A `matchit` object; see
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html).

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`oal_weight()`](https://kabajiro.github.io/oalasso/reference/oal_weight.md),
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)

## Examples

``` r
data("lalonde", package = "MatchIt")
fit <- oal(treat ~ age + educ + married + re74, data = lalonde,
           outcome = ~ re78)
m <- oal_match(fit, method = "nearest", caliper = 0.2)
```
