# ---------------------------------------------------------------------------
# bbglr_control(): one place to set the MCMC controls and model hyperparameters,
# so the main bbglr() signature stays short.
# ---------------------------------------------------------------------------

#' Control settings for [bbglr()]
#'
#' Bundles the MCMC controls and model hyperparameters into a single object, so
#' the main [bbglr()] call stays short. Pass it as `control = bbglr_control(...)`;
#' any value left at its default is used as-is.
#'
#' @param nIter,burnIn,thin MCMC controls passed to \pkg{BGLR}: total iterations,
#'   burn-in, and thinning interval.
#' @param nChains Integer number of independent chains (distinct seeds). Enables
#'   Gelman-Rubin diagnostics when `> 1`.
#' @param seed Base random seed; chain `c` uses `seed + c - 1`.
#' @param saveAt Directory for BGLR binary output. `NULL` (default) uses a fresh
#'   temporary directory.
#' @param verbose Logical; print BGLR progress.
#' @param exp_var_rank Genomic low-rank control for GBLUP-style terms
#'   (`mrk()` in GBLUP mode and `vm()`). A number in `(0, 1]` keeps the smallest
#'   set of leading principal components of the genomic relationship that together
#'   explain at least this fraction of its total variance (default `0.99`), giving
#'   a principled, self-scaling low-rank / PC-GBLUP fit that shrinks the design
#'   and speeds up sampling while retaining essentially all of the genetic signal.
#'   When \pkg{RSpectra} is available the required components are found by growing
#'   the top eigenpairs incrementally (the total variance is the free-to-compute
#'   trace of the relationship), so for a low-rank panel the `n x n` matrix is
#'   never formed. Set to `1` or `NA` to keep every (non-null) component. An
#'   explicit integer `rank =` inside `mrk()` overrides this for that term.
#' @return A list of class `breedRB_control`.
#' @examples
#' bbglr_control(nIter = 20000, burnIn = 5000, exp_var_rank = 0.95)
#' @seealso [bbglr()]
#' @export
bbglr_control <- function(nIter = 5000, burnIn = 1000, thin = 5, nChains = 1,
                          seed = 123, saveAt = NULL, verbose = TRUE,
                          exp_var_rank = 0.99) {
  .chk_count <- function(x, nm, min = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min ||
        x != as.integer(x)) {
      stop("`", nm, "` must be a single integer >= ", min, ".", call. = FALSE)
    }
    as.integer(x)
  }
  nIter   <- .chk_count(nIter,   "nIter")
  burnIn  <- .chk_count(burnIn,  "burnIn", min = 0L)
  thin    <- .chk_count(thin,    "thin")
  nChains <- .chk_count(nChains, "nChains")
  seed    <- .chk_count(seed,    "seed", min = -.Machine$integer.max)
  if (burnIn >= nIter) {
    stop("`burnIn` (", burnIn, ") must be smaller than `nIter` (", nIter, ").",
         call. = FALSE)
  }
  if (!is.null(saveAt) && (!is.character(saveAt) || length(saveAt) != 1L)) {
    stop("`saveAt` must be NULL or a single directory path.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("`verbose` must be a single logical.", call. = FALSE)
  }
  # exp_var_rank: NULL/NA = no variance-based truncation; else in (0, 1].
  if (is.null(exp_var_rank)) exp_var_rank <- NA_real_
  if (!is.na(exp_var_rank) &&
      (!is.numeric(exp_var_rank) || length(exp_var_rank) != 1L ||
       exp_var_rank <= 0 || exp_var_rank > 1)) {
    stop("`exp_var_rank` must be NA (keep all) or a single number in (0, 1].",
         call. = FALSE)
  }

  structure(
    list(nIter = nIter, burnIn = burnIn, thin = thin, nChains = nChains,
         seed = seed, saveAt = saveAt, verbose = verbose,
         exp_var_rank = as.numeric(exp_var_rank)),
    class = "breedRB_control"
  )
}

#' @export
print.breedRB_control <- function(x, ...) {
  cat("<breedRB_control>\n")
  cat("  nIter: ", x$nIter, "  burnIn: ", x$burnIn, "  thin: ", x$thin,
      "  nChains: ", x$nChains, "\n", sep = "")
  cat("  seed: ", x$seed, "  verbose: ", x$verbose, "\n", sep = "")
  cat("  exp_var_rank: ",
      if (is.na(x$exp_var_rank)) "NA (keep all components)"
      else format(x$exp_var_rank), "\n", sep = "")
  invisible(x)
}
