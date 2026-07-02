# Summarize an oal object

Produces (a) the full per-covariate selection table (`selected`), (b) a
path summary (the best grid row per `lambda2` value, when the path was
kept), (c) the outcome (weight) model coefficient table, (d) the full
balance table (unweighted vs. weighted SMD and KS, with a `*` marker at
weighted SMD \> 0.1), and (e) the near-positivity clip diagnostic.

## Usage

``` r
# S3 method for class 'oal'
summary(object, ...)

# S3 method for class 'summary.oal'
print(x, digits = 3, ...)
```

## Arguments

- object:

  An `oal` object.

- ...:

  Ignored.

- x:

  A `summary.oal` object.

- digits:

  Number of significant digits to print. Default 3.

## Value

For `summary.oal()`, an object of class `"summary.oal"`: a list with
elements `call`, `provenance`, `method`, `estimand`, `refit`, `lambda`,
`criterion`, `criterion.value`, `selected`, `path.summary`,
`outcome.model`, `balance`, `clip`, `clip.share`, and `nn`.
`print.summary.oal()` returns `x` invisibly.

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`print.oal()`](https://kabajiro.github.io/oalasso/reference/print.oal.md),
[`bal.tab.oal()`](https://kabajiro.github.io/oalasso/reference/bal.tab.oal.md)
