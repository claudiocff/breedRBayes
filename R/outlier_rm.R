#' MAD score of a numeric vector about its median
#'
#' Robust standardised distance `|x - median(x)| / (constant * MAD)`, the basis of
#' the Hampel / median-absolute-deviation outlier rule. When the MAD is zero (more
#' than half the values are identical, so the scaled MAD collapses) the spread is
#' re-estimated from the mean absolute deviation about the median, itself rescaled
#' to be consistent for the normal. If both are zero (a constant vector) every
#' score is `0` — nothing is flagged.
#'
#' @param x Numeric vector (may contain `NA`, which score as `NA`).
#' @param constant Consistency constant making the scaled MAD an unbiased
#'   estimator of the standard deviation under normality (`1.4826`).
#' @return A numeric vector of robust z-scores, the same length as `x`.
#' @keywords internal
.mad_score <- function(x, constant = 1.4826) {
  med <- stats::median(x, na.rm = TRUE)
  s   <- stats::mad(x, center = med, constant = constant, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {                 # >50% ties collapse the MAD
    s <- constant * mean(abs(x - med), na.rm = TRUE)  # fall back to a rescaled MeanAD
  }
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))  # constant vector: nothing to flag
  abs(x - med) / s
}

#' Grouping factor from a list of columns (or a single "all" stratum)
#' @keywords internal
.outlier_group <- function(group_cols, n) {
  if (length(group_cols)) interaction(group_cols, drop = TRUE, sep = ":")
  else factor(rep("all", n))
}

#' MAD outlier flags computed within strata
#'
#' Scores `y` with [.mad_score()] separately inside each level of `grp` and flags
#' values beyond `k`. Strata with fewer than `min_n` finite values are left
#' unscored (score `0`, never flagged) because the MAD is unreliable for tiny
#' groups. `NA` values score `NA` and are never flagged.
#'
#' @param y Numeric vector.
#' @param grp Grouping factor, same length as `y`.
#' @param k Threshold in robust standard deviations.
#' @param min_n Minimum finite observations per stratum to attempt scoring.
#' @return A list with the logical `flag` and numeric `score`, both length `y`.
#' @keywords internal
.mad_flag <- function(y, grp, k, min_n) {
  score <- rep(NA_real_, length(y))
  for (lv in levels(grp)) {
    idx <- which(grp == lv)
    if (!length(idx)) next
    if (sum(is.finite(y[idx])) < min_n) score[idx] <- 0     # too few to estimate a MAD
    else                                score[idx] <- .mad_score(y[idx])
  }
  flag <- is.finite(score) & score > k                      # NA -> not flagged
  list(flag = flag, score = score)
}

#' Robust MAD-based outlier detection and removal
#'
#' Flags observations that lie more than `k` robust standard deviations from the
#' median, using the median absolute deviation (MAD) rule
#' `|value - median| / (1.4826 * MAD) > k`. Unlike mean/standard-deviation
#' screening, the median and MAD are not themselves dragged around by the outliers
#' they are meant to catch, which makes the rule well suited to breeding-trial
#' data.
#'
#' `outlier_rm()` is generic. Screen a raw response with the [formula][stats::formula]
#' / vector method ([outlier_rm.default()]), or screen the **residuals of a fitted
#' model** with [outlier_rm.breedRB_fit()] — the latter accounts for the fixed and
#' random structure already in the model, so what remains to be judged is genuine
#' lack of fit rather than, say, a real environment effect.
#'
#' @param x A `breedRB_fit`, a `formula`, or a numeric vector — see the methods.
#' @param ... Passed to methods.
#' @return The cleaned data (see each method); always carries an `"outliers"`
#'   attribute describing what was flagged.
#' @seealso [outlier_rm.breedRB_fit()], [outlier_rm.default()].
#' @export
outlier_rm <- function(x, ...) {
  UseMethod("outlier_rm")
}

#' Screen a raw response for outliers (formula / vector method)
#'
#' The response and any grouping factors follow the package's formula convention:
#' the left-hand side is the trait to screen and the (optional) right-hand side
#' lists factors defining the strata within which the median and MAD are computed
#' — e.g. `yield ~ env` screens for outliers separately in each environment. A
#' bare numeric vector is also accepted (with `data = NULL`).
#'
#' @param x A formula `response ~ group1 + group2 ...` (screening within the
#'   crossed strata of the right-hand-side factors), a one-sided formula
#'   `response ~ 1` / `~ response` (screen the whole column), or a numeric vector.
#' @param data A `data.frame` holding the variables in `x`. Omit when `x` is a
#'   numeric vector.
#' @param by Optional grouping, as an alternative to putting factors on the
#'   formula's right-hand side: a factor / vector (vector input) or a character
#'   vector of column names (data-frame input). Combined with any RHS factors.
#' @param k Threshold in robust standard deviations. `3` is a common default;
#'   larger is more permissive. Must be positive.
#' @param action What to do with flagged observations: `"na"` sets the response to
#'   `NA` (default; keeps every row, so the design stays aligned), `"remove"`
#'   drops the offending rows, `"flag"` changes nothing but returns the flags.
#' @param min_n Strata with fewer than this many non-missing observations are left
#'   untouched (the MAD is unreliable for tiny groups). Default `5`.
#' @param verbose Logical; print a one-line summary of how many were flagged.
#' @param ... Unused.
#'
#' @return For a `data.frame` input, the data frame after applying `action`
#'   (`"flag"` returns it unchanged); an added logical column
#'   `<response>_outlier` marks the flagged rows when `action != "remove"`. For a
#'   vector input, the cleaned numeric vector, or the logical flag vector when
#'   `action = "flag"`. In all cases the result carries an `"outliers"` attribute:
#'   a list with the logical `flag`, the integer `index` of flagged rows, the
#'   robust `score`s and the `k` used.
#'
#' @examples
#' data(wheat, package = "breedRBayes")
#' # screen yield within each environment, blanking outliers to NA
#' clean <- outlier_rm(yield ~ env, data = wheat, k = 3)
#' # or drop the rows entirely
#' outlier_rm(yield ~ env, data = wheat, k = 4, action = "remove")
#' # plain vector
#' outlier_rm(c(1, 2, 3, 4, 100))
#' @export
outlier_rm.default <- function(x, data = NULL, by = NULL, k = 3,
                               action = c("na", "remove", "flag"),
                               min_n = 5L, verbose = TRUE, ...) {
  action <- match.arg(action)
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("outlier_rm(): `k` must be a single positive number.", call. = FALSE)
  }

  # --- resolve the response vector `y` and the grouping factors ---------------
  resp_name <- NULL
  grp_names <- character(0)

  if (inherits(x, "formula")) {
    if (is.null(data)) {
      stop("outlier_rm(): supply `data` when `x` is a formula.", call. = FALSE)
    }
    tt  <- stats::terms(x, data = data)
    lhs <- attr(tt, "response")
    vars <- as.character(attr(tt, "variables"))[-1]
    if (lhs == 0L) {                         # one-sided ~response : response is the sole term
      if (length(vars) != 1L) {
        stop("outlier_rm(): a one-sided formula must name exactly one response, ",
             "e.g. `~ yield`.", call. = FALSE)
      }
      resp_name <- vars[1]
    } else {
      resp_name <- vars[lhs]
      rhs <- attr(tt, "term.labels")
      grp_names <- rhs[rhs != "1"]           # `response ~ 1` -> no grouping
    }
    y <- data[[resp_name]]
    if (is.null(y)) {
      stop("outlier_rm(): response `", resp_name, "` not found in `data`.", call. = FALSE)
    }
  } else if (is.numeric(x)) {
    y <- x
  } else {
    stop("outlier_rm(): `x` must be a breedRB_fit, a formula or a numeric vector.",
         call. = FALSE)
  }
  if (!is.numeric(y)) {
    stop("outlier_rm(): the response must be numeric.", call. = FALSE)
  }

  # --- assemble the grouping key (RHS factors + `by`), one entry per row ------
  group_cols <- list()
  for (g in grp_names) {
    if (is.null(data[[g]])) {
      stop("outlier_rm(): grouping variable `", g, "` not found in `data`.", call. = FALSE)
    }
    group_cols[[g]] <- data[[g]]
  }
  group_cols <- .resolve_by(by, group_cols, data, length(y))

  grp <- .outlier_group(group_cols, length(y))

  # --- score within each stratum ---------------------------------------------
  ff    <- .mad_flag(y, grp, k, min_n)
  flag  <- ff$flag
  score <- ff$score
  index <- which(flag)

  if (verbose) {
    n_grp <- if (length(group_cols)) nlevels(grp) else 1L
    message("outlier_rm(): flagged ", length(index), " of ", length(y),
            " observation(s) beyond ", k, " MAD",
            if (n_grp > 1L) paste0(" (within ", n_grp, " strata)") else "",
            if (!is.null(resp_name)) paste0(" on `", resp_name, "`") else "", ".")
  }

  info <- list(flag = flag, index = index, score = score, k = k)

  # --- apply the requested action --------------------------------------------
  if (is.null(data)) {                                  # vector interface
    out <- switch(action,
      na     = { z <- y; z[flag] <- NA; z },
      remove = y[!flag],
      flag   = flag)
    attr(out, "outliers") <- info
    return(out)
  }

  out <- data                                           # data-frame interface
  if (action == "remove") {
    out <- out[!flag, , drop = FALSE]                   # `info` still indexes original rows
  } else {
    if (action == "na") out[[resp_name]][flag] <- NA
    out[[paste0(resp_name, "_outlier")]] <- flag        # keep a visible marker
  }
  attr(out, "outliers") <- info
  out
}

#' Add `by=` columns/vector to a list of grouping columns
#' @keywords internal
.resolve_by <- function(by, group_cols, data, n) {
  if (is.null(by)) return(group_cols)
  if (is.character(by) && !is.null(data)) {            # column names
    for (g in by) {
      if (is.null(data[[g]])) {
        stop("outlier_rm(): `by` column `", g, "` not found in the data.", call. = FALSE)
      }
      group_cols[[g]] <- data[[g]]
    }
  } else {                                             # a factor / vector aligned with y
    if (length(by) != n) {
      stop("outlier_rm(): `by` must have the same length as the response.", call. = FALSE)
    }
    group_cols[["by"]] <- by
  }
  group_cols
}

#' Remove outliers from a fitted model by screening its residuals
#'
#' Screens the **residuals** of a fitted [bbglr()] model — observed minus
#' posterior-mean fitted value ([residuals.breedRB_fit()]) — with the robust MAD
#' rule, and returns the model's data frame (`fit$data`) with the offending
#' observations removed (or blanked / flagged). Because the residual already
#' subtracts everything the model explains (the fixed part, the genomic / random
#' effects, and for a random regression the fitted reaction norm), an extreme
#' residual is genuine lack of fit for that record rather than a real design
#' effect — the appropriate quantity to screen before refitting.
#'
#' The usual workflow is to fit once, drop the residual outliers, and refit on the
#' returned data frame:
#' \preformatted{
#'   fit   <- bbglr(yield ~ 1 + env, random = ~ mrk(gen, M), data = dat,
#'                  relmat = list(M = M))
#'   clean <- outlier_rm(fit, k = 3)          # residual-based, returns a data.frame
#'   fit2  <- bbglr(yield ~ 1 + env, random = ~ mrk(gen, M), data = clean,
#'                  relmat = list(M = M))
#' }
#'
#' For a multi-trait fit each trait's residual column is screened separately; a
#' row is flagged when **any** trait residual is an outlier. `action = "na"`
#' blanks only the offending trait cell(s); `action = "remove"` drops the whole
#' row.
#'
#' @param x A `breedRB_fit` from [bbglr()].
#' @param k Threshold in robust standard deviations (default `3`; larger is more
#'   permissive). Must be positive.
#' @param action What to do with a flagged record: `"remove"` (default) drops the
#'   row, `"na"` sets the offending response(s) to `NA` (keeping the row so the
#'   design stays aligned) and adds a `<trait>_outlier` marker column, `"flag"`
#'   changes nothing but adds the marker column(s).
#' @param by Optional grouping for the residual screen: a character vector of
#'   `fit$data` column names, or a factor / vector aligned to the modelled rows.
#'   Residuals are usually mean-zero across the model's strata already, so a global
#'   screen (the default, `NULL`) is normally what you want; use `by` only to allow
#'   a different residual scale per group.
#' @param min_n Strata with fewer than this many finite residuals are left
#'   untouched. Default `5`.
#' @param verbose Logical; print a one-line summary of how many were flagged.
#' @param ... Unused.
#'
#' @return The model's data frame (`fit$data`) after applying `action`: rows
#'   removed for `"remove"`, else the same rows with the response(s) blanked
#'   (`"na"`) or unchanged (`"flag"`) plus a logical `<trait>_outlier` column. The
#'   result carries an `"outliers"` attribute — a list with the logical `flag`
#'   (per modelled row), the integer `index` of flagged rows, the screened
#'   `residuals`, the robust `score`s and the `k` used.
#'
#' @seealso [residuals.breedRB_fit()], [fitted.breedRB_fit()], [model_fit()],
#'   [outlier_rm.default()].
#' @export
outlier_rm.breedRB_fit <- function(x, k = 3,
                                   action = c("remove", "na", "flag"),
                                   by = NULL, min_n = 5L, verbose = TRUE, ...) {
  fit    <- x
  action <- match.arg(action)
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("outlier_rm(): `k` must be a single positive number.", call. = FALSE)
  }

  data   <- fit$data
  traits <- fit$response$traits
  mt     <- isTRUE(fit$response$multitrait)
  res    <- stats::residuals(fit)                       # vector (ST) or matrix (MT)
  resm   <- if (mt) as.matrix(res) else matrix(res, ncol = 1L)
  n      <- nrow(resm)

  # --- grouping for the residual screen (global unless `by` is given) ---------
  group_cols <- .resolve_by(by, list(), data, n)
  grp <- .outlier_group(group_cols, n)

  # --- score each trait's residuals; a row is flagged if any trait is out -----
  flagm  <- matrix(FALSE, nrow = n, ncol = ncol(resm))
  scorem <- matrix(NA_real_, nrow = n, ncol = ncol(resm))
  for (j in seq_len(ncol(resm))) {
    ff <- .mad_flag(resm[, j], grp, k, min_n)
    flagm[, j]  <- ff$flag
    scorem[, j] <- ff$score
  }
  rowflag <- apply(flagm, 1L, any)
  index   <- which(rowflag)

  if (verbose) {
    n_grp <- if (length(group_cols)) nlevels(grp) else 1L
    message("outlier_rm(): flagged ", length(index), " of ", n,
            " modelled record(s) with residuals beyond ", k, " MAD",
            if (n_grp > 1L) paste0(" (within ", n_grp, " strata)") else "",
            if (mt) paste0(" across ", length(traits), " trait(s)") else "", ".")
  }

  info <- list(flag = rowflag, index = index,
               residuals = if (mt) resm else drop(resm),
               score = if (mt) scorem else drop(scorem), k = k)

  # --- return the cleaned data frame ------------------------------------------
  out <- data
  if (action == "remove") {
    out <- out[!rowflag, , drop = FALSE]                 # `info` still indexes modelled rows
  } else {
    for (j in seq_along(traits)) {
      if (action == "na") out[[traits[j]]][flagm[, j]] <- NA
      out[[paste0(traits[j], "_outlier")]] <- flagm[, j]   # per-trait marker
    }
  }
  attr(out, "outliers") <- info
  out
}
