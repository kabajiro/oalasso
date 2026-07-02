# utils.R -- internal utilities: messaging, dependency guards, tie-breaking,
# small helpers. Nothing here is exported.

# Verbose messaging helper. All progress output goes through here so that
# `verbose = FALSE` (the default) is completely silent.
.vmsg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

# Guard for Suggests packages: every conditional dependency error names the
# missing package and the exact install.packages() call (psAve convention).
.require_pkg <- function(pkg, purpose) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(paste0("Package \"%s\" is required %s but is not installed.\n",
                        "Install it with: install.packages(\"%s\")"),
                 pkg, purpose, pkg),
         call. = FALSE)
  }
  invisible(TRUE)
}

# First minimum with a small relative numerical tolerance (D.13). The
# documented rule is "the FIRST grid row attaining the minimum" in the grid
# order lambda2-ascending (0 first) outer x delta-ascending inner, so exact
# ties prefer lambda2 = 0 (plain OAL) and then the smallest penalty.
#
# Equivalence to Balde's NESTED rule (2025 supplement, verified 2026-07-02):
# Balde minimizes wAMD over (delta, gamma) per lambda2 (which.min = first
# min), then takes the first lambda2 attaining the smallest per-lambda2
# minimum. Because the flat grid enumerates lambda2 blocks in ascending
# order with delta ascending inside each block, the flat first-minimum
# equals the nested first-minimum: a lambda2 block contains a grid point
# within tolerance of the global minimum iff its block minimum does, so
# both rules land on the same block and, inside it, on the same first
# delta. Pinned by the constructed-tie test in test-goal.R.
.first_min <- function(x, tol = 1e-9) {
  vmin <- min(x)
  which(x <= vmin + tol * max(1, abs(vmin)))[1L]
}

# TRUE for columns taking exactly two distinct values (used by the cobalt
# balance table; identical to psAve's detection).
.detect_bin <- function(mat) {
  apply(mat, 2L, function(x) length(unique(x)) == 2L)
}

# cobalt s.d.denom implied by the estimand (cobalt's own convention): the
# treated-group SD for the ATT, the pooled SD for the ATE. Used only by the
# display-oriented balance table and bal.tab.oal(); the wAMD selection
# criterion itself involves no standardization (D.12).
.sd_denom <- function(estimand) {
  if (identical(estimand, "ATT")) "treated" else "pooled"
}

# Clip a propensity score vector to [clip[1], clip[2]].
.clip_ps <- function(ps, clip) {
  pmin(pmax(ps, clip[1L]), clip[2L])
}

# Percentage of propensity scores sitting AT either clipping bound (D.15).
.clip_share <- function(ps, clip) {
  100 * mean(ps <= clip[1L] | ps >= clip[2L])
}

# The near-positivity diagnostic message (D.15, exact wording).
.clip_note <- function(share) {
  sprintf(paste0("%.1f%% of propensity scores are at the clipping bound; ",
                 "near-positivity violations destabilize wAMD selection. ",
                 "Consider method = 'goal', whose L2 term stabilizes the weights ",
                 "(Jones, Ertefaie & Shortreed 2023, doi:10.1111/biom.13681)."),
          share)
}

# Provenance label (D.19, exact strings; printed first by print()/summary()).
.provenance <- function(method, lambda2.user = FALSE, user.coef = FALSE,
                        family = "gaussian", refit = FALSE) {
  base <- if (identical(method, "goal")) {
    paste0("GOAL (Balde, Yang & Lefebvre 2023, doi:10.1111/biom.13683) -- ",
           if (isTRUE(lambda2.user)) "user-specified lambda2 grid"
           else "lambda2 grid: author's published grid (Balde 2025 supplement)")
  } else {
    "OAL (Shortreed & Ertefaie 2017, doi:10.1111/biom.12679)"
  }
  if (isTRUE(user.coef)) {
    base <- paste0(base,
                   "; user-supplied outcome coefficients -- screening use only; no oracle property")
  }
  if (identical(family, "binomial")) {
    base <- paste0(base, " [binomial outcome model: beyond the validated simulations]")
  }
  if (isTRUE(refit)) {
    base <- paste0(base,
                   " [post-selection refit: Schnitzer et al. 2025 component, adapted from a longitudinal setting]")
  }
  base
}

# Format a numeric scalar for messages.
.fmt <- function(x, digits = 4) format(x, digits = digits)
