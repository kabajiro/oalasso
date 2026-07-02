# Predict outcome-adaptive lasso propensity scores for new data

Computes propensity scores (or the linear predictor) for new
observations from the stored terms, factor levels, and ORIGINAL-scale
coefficients of an
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) fit.
Works without `keep.fits = TRUE` and for both `refit` modes; the
training standardization is never re-estimated (the original-scale
coefficients already absorb it). Missing values in `newdata` are an
error.

## Usage

``` r
# S3 method for class 'oal'
predict(object, newdata = NULL, type = c("ps", "link"), ...)
```

## Arguments

- object:

  An `oal` object.

- newdata:

  A data frame containing the variables of the propensity score formula,
  or `NULL` (default) for the in-sample values.

- type:

  `"ps"` (default): propensity scores clipped to the fit's `clip`;
  `"link"`: the unclipped linear predictor.

- ...:

  Ignored.

## Value

A numeric vector with one value per row, named by rownames. For
`newdata = NULL` and `type = "ps"` this is exactly `object$ps`.

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`coef.oal()`](https://kabajiro.github.io/oalasso/reference/coef.oal.md)
