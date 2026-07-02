# plot.R -- plot.oal: "wamd" (criterion path), "coef" (coefficient paths),
# "balance" (Love plot via cobalt). Base graphics except the Love plot.

#' Plot an oal object
#'
#' Three diagnostic displays for an [oal()] fit:
#' \describe{
#'   \item{`"wamd"`}{the wAMD selection criterion against the exponent
#'     \eqn{\delta} (`lambda`), one curve per `lambda2` value for
#'     `method = "goal"`; the selected point is marked. Requires
#'     `keep.path = TRUE`.}
#'   \item{`"coef"`}{the original-scale propensity score coefficient paths
#'     against \eqn{\delta} at the selected `lambda2`, one line per covariate;
#'     covariates driven to zero (instruments, noise) are visible as lines
#'     absorbed into the axis. Requires `keep.path = TRUE`.}
#'   \item{`"balance"`}{a Love plot of covariate balance before/after
#'     weighting, via [cobalt::love.plot()] (dispatched through
#'     [bal.tab.oal()]; \pkg{cobalt} is an Import, so always available).}
#' }
#'
#' @param x An `oal` object.
#' @param type One of `"wamd"` (default), `"coef"`, `"balance"`.
#' @param ... For `"balance"`, further arguments to [cobalt::love.plot()]
#'   (e.g., `thresholds = 0.1`); otherwise further graphical parameters passed
#'   to the base plotting calls.
#'
#' @return For `"balance"`, the `ggplot` object from [cobalt::love.plot()]
#'   (invisibly, after printing); otherwise `x`, invisibly.
#' @seealso [oal()], [bal.tab.oal()], [cobalt::love.plot()]
#' @export
plot.oal <- function(x, type = c("wamd", "coef", "balance"), ...) {
  type <- match.arg(type)
  switch(type,
         wamd = .plot_wamd(x, ...),
         coef = .plot_coef(x, ...),
         balance = .plot_balance(x, ...))
}

.need_path <- function(x) {
  if (is.null(x$path)) {
    stop("No tuning path is stored (keep.path = FALSE). Re-run oal() with keep.path = TRUE.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.plot_wamd <- function(x, ...) {
  .need_path(x)
  path <- x$path
  l2 <- path$lambda2
  groups <- if (all(is.na(l2))) list(seq_len(nrow(path))) else
    split(seq_len(nrow(path)), l2)
  graphics::plot(NA, xlim = range(path$delta), ylim = range(path$wamd),
                 xlab = expression(delta ~ "(lambda exponent;" ~ lambda[n] == n^delta * ")"),
                 ylab = "wAMD", ...)
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    o <- order(path$delta[idx])
    graphics::lines(path$delta[idx][o], path$wamd[idx][o], col = k, lty = k)
  }
  sel <- which(path$selected)
  graphics::points(path$delta[sel], path$wamd[sel], pch = 19, col = "firebrick")
  if (length(groups) > 1L) {
    graphics::legend("topright",
                     legend = sprintf("lambda2 = %s", names(groups)),
                     col = seq_along(groups), lty = seq_along(groups),
                     bty = "n", cex = 0.8)
  }
  invisible(x)
}

.plot_coef <- function(x, ...) {
  .need_path(x)
  path <- x$path
  ## coefficient paths at the SELECTED lambda2 (all rows for method = "oal")
  rows <- if (all(is.na(path$lambda2))) seq_len(nrow(path)) else
    which(path$lambda2 == x$lambda[["lambda2"]])
  o <- rows[order(path$delta[rows])]
  B <- x$coef.path[-1L, o, drop = FALSE]   # drop the intercept row
  d <- path$delta[o]
  graphics::matplot(d, t(B), type = "l", lty = 1, col = seq_len(nrow(B)),
                    xlab = expression(delta ~ "(lambda exponent;" ~ lambda[n] == n^delta * ")"),
                    ylab = "coefficient (original scale)", ...)
  graphics::abline(h = 0, col = "grey70")
  graphics::abline(v = x$lambda[["delta"]], col = "firebrick", lty = 3)
  graphics::legend("topleft", legend = rownames(B), col = seq_len(nrow(B)),
                   lty = 1, bty = "n", cex = 0.7)
  invisible(x)
}

.plot_balance <- function(x, ...) {
  p <- cobalt::love.plot(x, ...)
  print(p)
  invisible(p)
}
