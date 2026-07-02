# input.R -- input validation, treatment coercion, model matrix, and the
# fixed standardization protocol (D.1-D.3). Nothing here is exported.

# Coerce a treatment vector to integer 0/1 per MatchIt convention (D.1):
# numeric 0/1 kept; logical -> 0/1; two-level factor/character -> the SECOND
# level (in factor-level order) is treated = 1; anything else errors.
# (Identical to psAve's .coerce_treat -- suite convention.)
.coerce_treat <- function(a, name = "treatment") {
  if (is.numeric(a)) {
    u <- unique(a[!is.na(a)])
    if (!all(u %in% c(0, 1))) {
      stop(sprintf("The %s variable must be binary: numeric 0/1, logical, or a two-level factor/character.",
                   name), call. = FALSE)
    }
    out <- as.integer(a)
  } else if (is.logical(a)) {
    out <- as.integer(a)
  } else if (is.factor(a) || is.character(a)) {
    f <- droplevels(as.factor(a))
    if (nlevels(f) != 2L) {
      stop(sprintf("The %s variable must have exactly two levels; found %d (%s).",
                   name, nlevels(f), paste0('"', levels(f), '"', collapse = ", ")),
           call. = FALSE)
    }
    out <- as.integer(f == levels(f)[2L])
  } else {
    stop(sprintf("The %s variable must be binary: numeric 0/1, logical, or a two-level factor/character.",
                 name), call. = FALSE)
  }
  out
}

# Error (never drop) on missing values, naming the offending variables (D.1).
.check_no_na <- function(data, vars, what) {
  vars <- vars[vars %in% names(data)]
  if (!length(vars)) return(invisible(TRUE))
  bad <- vars[vapply(vars, function(v) anyNA(data[[v]]), logical(1L))]
  if (length(bad)) {
    stop(sprintf(paste0("Missing values found in %s variable(s): %s.\n",
                        "oalasso requires complete cases in all used variables and never drops rows silently; ",
                        "handle missing data (e.g., by imputation) before calling oal()."),
                 what, paste0('"', bad, '"', collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

# Model matrix for the PS model (D.2): default treatment contrasts, intercept
# column dropped -- X <- model.matrix(formula, data)[, -1, drop = FALSE].
# (This is the expansion the OAL logistic PS model itself uses; the cobalt
# full-dummy expansion would create exactly collinear columns and trip the
# degenerate-beta error of D.4 for every factor.) Stores terms/xlevels/
# contrasts for predict().
.build_design <- function(tt, data, xlev = NULL, contrasts = NULL,
                          what = "formula") {
  mf <- stats::model.frame(tt, data = data, na.action = stats::na.pass,
                           xlev = xlev, drop.unused.levels = is.null(xlev))
  if (anyNA(mf)) {
    stop(sprintf(paste0("Missing values found among the %s variables. ",
                        "oalasso requires complete cases and never drops rows silently."),
                 what), call. = FALSE)
  }
  for (j in seq_along(mf)) {
    if (is.character(mf[[j]])) mf[[j]] <- factor(mf[[j]])
  }
  X <- stats::model.matrix(tt, data = mf, contrasts.arg = contrasts)
  contr <- attr(X, "contrasts")
  X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]
  if (ncol(X) == 0L) {
    stop("`formula` must contain at least one covariate on the right-hand side.",
         call. = FALSE)
  }
  storage.mode(X) <- "double"
  list(X = X,
       xlev = stats::.getXlevels(tt, mf),
       contrasts = contr)
}

# Full input processing for oal() (D.1-D.3). Returns everything the
# orchestrator needs; errors early with actionable messages.
.process_inputs <- function(formula, data, outcome, family, need.outcome) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula of the form `treat ~ x1 + x2 + ...`.",
         call. = FALSE)
  }
  if (missing(data) || is.null(data)) {
    stop("`data` must be supplied as a data.frame.", call. = FALSE)
  }
  data <- as.data.frame(data)
  if (nrow(data) == 0L) stop("`data` has zero rows.", call. = FALSE)
  rn <- rownames(data)

  ## --- treatment (D.1) ------------------------------------------------------
  tname <- deparse1(formula[[2L]])
  a <- eval(formula[[2L]], data, environment(formula))
  if (length(a) != nrow(data)) {
    stop(sprintf("The treatment variable `%s` has length %d; expected %d (one per row of `data`).",
                 tname, length(a), nrow(data)), call. = FALSE)
  }
  if (anyNA(a)) {
    stop(sprintf(paste0("Missing values found in the treatment variable `%s`. ",
                        "oalasso requires complete cases and never drops rows silently."),
                 tname), call. = FALSE)
  }
  treat <- .coerce_treat(a, name = sprintf("treatment (`%s`)", tname))
  if (sum(treat == 1L) < 2L || sum(treat == 0L) < 2L) {
    stop(sprintf("At least 2 units are required in each treatment arm; found %d treated and %d control.",
                 sum(treat == 1L), sum(treat == 0L)), call. = FALSE)
  }

  ## --- covariates / model matrix (D.2) --------------------------------------
  tt <- stats::delete.response(stats::terms(formula, data = data))
  if (length(attr(tt, "term.labels")) == 0L) {
    stop("`formula` must contain at least one covariate on the right-hand side.",
         call. = FALSE)
  }
  .check_no_na(data, all.vars(tt), "covariate")
  des <- .build_design(tt, data, what = "covariate (formula right-hand side)")
  X <- des$X
  n <- nrow(X)  # n is ALWAYS nrow(X); n^delta never resolves from any other scope (D.2)
  if (n != nrow(data)) {
    stop("Internal error: model matrix rows do not match `data` rows.", call. = FALSE)
  }
  p <- ncol(X)
  if (p < 2L) {
    stop(paste0("The propensity score model must contain at least two covariate columns: ",
                "glmnet requires >= 2 columns and outcome-adaptive selection among a single ",
                "covariate is not meaningful."), call. = FALSE)
  }
  const <- colnames(X)[apply(X, 2L, function(col) length(unique(col)) == 1L)]
  if (length(const)) {
    stop(sprintf(paste0("Zero-variance covariate column(s): %s. ",
                        "Constant columns cannot be standardized or selected; remove them from `formula`."),
                 paste0('"', const, '"', collapse = ", ")), call. = FALSE)
  }

  ## --- standardization protocol (D.3, fixed) --------------------------------
  ## scale(X, TRUE, TRUE) ONCE with training statistics, ALL columns including
  ## dummies; the weight model AND glmnet are both fit on this SAME matrix with
  ## glmnet(standardize = FALSE). Never per-subset scaling.
  Xs <- scale(X, center = TRUE, scale = TRUE)
  center <- attr(Xs, "scaled:center")
  scl <- attr(Xs, "scaled:scale")
  attr(Xs, "scaled:center") <- NULL
  attr(Xs, "scaled:scale") <- NULL

  ## --- outcome (D.1) ---------------------------------------------------------
  y <- NULL
  outcome.name <- NA_character_
  if (need.outcome && is.null(outcome)) {
    stop(paste0("`outcome` is required (unless `outcome.coef` is supplied): OAL derives its ",
                "penalty weights from an outcome model on the same covariates as the PS model; ",
                "supply outcome = ~ y."), call. = FALSE)
  }
  ## Parse and validate `outcome` whenever it is SUPPLIED, even when the
  ## penalty weights come from `outcome.coef` -- the outcome-leakage guard
  ## must still run (a supplied outcome variable may not appear among the PS
  ## covariates). When `outcome.coef` is used WITHOUT `outcome`, no leakage
  ## check is possible; this is documented in ?oal.
  if (!is.null(outcome)) {
    if (!inherits(outcome, "formula")) {
      stop("`outcome` must be a one-sided formula `~ y` naming the outcome variable.",
           call. = FALSE)
    }
    if (length(outcome) == 3L) {
      stop(paste0("OAL derives its penalty weights from an outcome model on the same covariates ",
                  "as the PS model; supply outcome = ~ y."), call. = FALSE)
    }
    yexpr <- outcome[[2L]]
    yvars <- all.vars(yexpr)
    if (length(yvars) != 1L) {
      stop("`outcome` must be a one-sided formula `~ y` with exactly one variable.",
           call. = FALSE)
    }
    outcome.name <- deparse1(yexpr)
    .check_no_na(data, yvars, "outcome")
    y <- eval(yexpr, data, environment(outcome))
    if (length(y) != n) {
      stop(sprintf("The outcome variable `%s` has length %d; expected %d (one per row of `data`).",
                   outcome.name, length(y), n), call. = FALSE)
    }
    if (anyNA(y)) {
      stop(sprintf(paste0("Missing values found in the outcome variable `%s`. ",
                          "oalasso requires complete cases and never drops rows silently."),
                   outcome.name), call. = FALSE)
    }
    if (any(yvars %in% all.vars(tt))) {
      stop(sprintf("The outcome variable `%s` must not appear among the propensity score covariates in `formula`.",
                   outcome.name), call. = FALSE)
    }
    if (identical(family$family, "binomial")) {
      y <- as.numeric(.coerce_treat(y, name = sprintf("binary outcome (`%s`)", outcome.name)))
    } else {
      if (!is.numeric(y)) {
        stop(sprintf("With family = gaussian(), the outcome variable `%s` must be numeric.",
                     outcome.name), call. = FALSE)
      }
      y <- as.numeric(y)
    }
  }

  list(data = data, n = n, p = p, rn = rn,
       treat = treat, treat.name = tname,
       X = X, Xs = Xs, center = center, scale = scl,
       bin.vars = .detect_bin(X),
       terms = tt, xlevels = des$xlev, contrasts = des$contrasts,
       y = y, outcome.name = outcome.name)
}

# Validate a user-supplied `outcome.coef` vector (D.5): named numeric whose
# names exactly cover the standardized model-matrix columns; all finite; not
# all zero. Zeros ARE allowed here: |0|^-gamma = Inf, which glmnet natively
# converts to `exclude` (a hard drop). Returned reordered to colnames(Xs).
.validate_outcome_coef <- function(oc, cn) {
  expected <- paste0('"', cn, '"', collapse = ", ")
  if (!is.numeric(oc) || is.null(names(oc))) {
    stop(sprintf(paste0("`outcome.coef` must be a NAMED numeric vector over the standardized ",
                        "model-matrix columns; expected names: %s."), expected),
         call. = FALSE)
  }
  nm <- names(oc)
  if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) ||
      length(oc) != length(cn) || !setequal(nm, cn)) {
    stop(sprintf(paste0("The names of `outcome.coef` must exactly cover the model-matrix ",
                        "columns (one coefficient per column, no duplicates); expected names: %s."),
                 expected), call. = FALSE)
  }
  if (any(!is.finite(oc))) {
    stop("All values of `outcome.coef` must be finite (zeros are allowed and cause a hard drop of the column).",
         call. = FALSE)
  }
  if (all(oc == 0)) {
    stop("`outcome.coef` must not be all zero: at least one covariate must carry a finite adaptive penalty.",
         call. = FALSE)
  }
  oc <- oc[cn]
  storage.mode(oc) <- "double"
  oc
}
