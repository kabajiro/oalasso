# Plot an oal object

Three diagnostic displays for an
[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md) fit:

- `"wamd"`:

  the wAMD selection criterion against the exponent \\\delta\\
  (`lambda`), one curve per `lambda2` value for `method = "goal"`; the
  selected point is marked. Requires `keep.path = TRUE`.

- `"coef"`:

  the original-scale propensity score coefficient paths against
  \\\delta\\ at the selected `lambda2`, one line per covariate;
  covariates driven to zero (instruments, noise) are visible as lines
  absorbed into the axis. Requires `keep.path = TRUE`.

- `"balance"`:

  a Love plot of covariate balance before/after weighting, via
  [`cobalt::love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.html)
  (dispatched through
  [`bal.tab.oal()`](https://kabajiro.github.io/oalasso/reference/bal.tab.oal.md);
  cobalt is an Import, so always available).

## Usage

``` r
# S3 method for class 'oal'
plot(x, type = c("wamd", "coef", "balance"), ...)
```

## Arguments

- x:

  An `oal` object.

- type:

  One of `"wamd"` (default), `"coef"`, `"balance"`.

- ...:

  For `"balance"`, further arguments to
  [`cobalt::love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.html)
  (e.g., `thresholds = 0.1`); otherwise further graphical parameters
  passed to the base plotting calls.

## Value

For `"balance"`, the `ggplot` object from
[`cobalt::love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.html)
(invisibly, after printing); otherwise `x`, invisibly.

## See also

[`oal()`](https://kabajiro.github.io/oalasso/reference/oal.md),
[`bal.tab.oal()`](https://kabajiro.github.io/oalasso/reference/bal.tab.oal.md),
[`cobalt::love.plot()`](https://ngreifer.github.io/cobalt/reference/love.plot.html)
