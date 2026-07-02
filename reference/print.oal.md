# Print an oal object

Prints a one-screen teaching summary of a fitted
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) object:
the provenance label first (which method, which tuning provenance), the
sample, the selected tuning point and its wAMD, the retained
(outcome-related) versus excluded (outcome-unrelated) covariates with
their outcome-coefficient magnitudes, the fixed instrument-exclusion
rationale, and then the **literal next call** – echoing the formula and
data name from your own
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) call –
that hands the score to
[`MatchIt::matchit()`](https://kosukeimai.github.io/MatchIt/reference/matchit.html)
or
[`WeightIt::weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.html).
A near-positivity diagnostic is appended when more than 5\\ clipping
bound.

## Usage

``` r
# S3 method for class 'oal'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An `oal` object.

- digits:

  Number of significant digits to print. Default 3.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`summary.oal()`](https://kabajiro.github.io/oalasso/reference/summary.oal.md)
