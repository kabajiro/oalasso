# Extract the inverse-probability weights of an oal fit

Returns `object$weights`: the IPW weights at the *fitted* estimand
implied by the clipped propensity score (ATE: \\A/e + (1-A)/(1-e)\\;
ATT: \\A + (1-A)\\e/(1-e)\\). For other estimands use
`WeightIt::get_w_from_ps(object$ps, object$treat, estimand = ...)` or
[`oal_weight()`](https://kabajiro.github.io/oalasso/reference/oal_weight.md).

## Usage

``` r
# S3 method for class 'oal'
weights(object, ...)
```

## Arguments

- object:

  An `oal` object.

- ...:

  Ignored.

## Value

The numeric vector `object$weights`, named by `rownames(data)`.

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`oal_weight()`](https://kabajiro.github.io/oalasso/reference/oal_weight.md),
[`oal_wamd()`](https://kabajiro.github.io/oalasso/reference/oal_wamd.md)
