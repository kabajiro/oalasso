# Extract propensity score model coefficients

Returns the selected outcome-adaptive lasso propensity score
coefficients (intercept first), either on the original covariate scale
(default; the scale used by
[`predict.oal()`](https://kabajiro.github.io/oalasso/reference/predict.oal.md))
or on the standardized scale that the penalized objective was actually
optimized on (the scale of the adaptive weights and of the wAMD –
methods-paper material).

## Usage

``` r
# S3 method for class 'oal'
coef(object, scale = c("original", "standardized"), ...)
```

## Arguments

- object:

  An `oal` object.

- scale:

  `"original"` (default) or `"standardized"`.

- ...:

  Ignored.

## Value

A named numeric vector of length p + 1.

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`predict.oal()`](https://kabajiro.github.io/oalasso/reference/predict.oal.md)
