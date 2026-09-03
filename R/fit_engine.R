# ---------------------------------------------------------------------------
# bbglr(): main entry point. Parses formulas, builds the BGLR ETA, and runs one
# or more MCMC chains, returning a `breedRB_fit` object.
# ---------------------------------------------------------------------------

#' Fit a Bayesian mixed model with the BGLR engine using asreml-style formulas
#'
#' @param fixed Two-sided formula for the fixed/mean part, e.g. `yield ~ 1 + env`.
#'   Use `cbind(t1, t2) ~ ...` for a multi-trait model.
#' @param random One-sided formula of random terms, e.g. `~ mrk(gen, M)`. Special
#'   functions: `mrk(f, M, method)` genomic effect from a **marker matrix** `M`
#'   that automatically fits GBLUP or RR-BLUP, whichever is cheaper (`method` is
#'   `"auto"`, `"GBLUP"` or `"RRBLUP"`; `"auto"` picks GBLUP when markers >=
#'   genotypes, else RR-BLUP — the two are prediction-equivalent). A `mrk()` fit
#'   can score new genotypes with [predict.breedRB_fit()] and recover marker
#'   effects with [solve_SNP()]; it is the recommended genomic term. `vm(f, K)`
#'   effect of factor `f` with a supplied covariance matrix `K` — use it for a
#'   kernel you already have and cannot rebuild from markers (pedigree A-matrix,
#'   environmic kernel, externally-computed G, or any custom kernel). `leg(x, n)`
#'   / `rr(x, n)` random regression of order `n`; `fa(f, k)` / `rrc(f, k)`
#'   factor-analytic. Bare names are factors; `a:b` is an interaction.
#' @param residual One-sided residual formula. `NULL`/`~ units` = homogeneous;
#'   `~ dsum(~units | env)` = heterogeneous residual variances by `env`.
#' @param data A data frame.
#' @param relmat Named list of matrices referenced by `mrk()` and `vm()`. For
#'   `mrk()` the entry is a marker matrix with genotypes in rows (row names =
#'   genotype IDs) and markers in columns. For `vm()` the entry is a covariance /
#'   kernel matrix (K, not its inverse) — genomic, pedigree, environmic or any
#'   custom kernel — with row/column names matching the factor levels.
#' @param nIter,burnIn,thin MCMC controls passed to BGLR.
#' @param nChains Integer number of independent chains (distinct seeds). Enables
#'   Gelman-Rubin diagnostics when `> 1`.
#' @param seed Base random seed; chain `c` uses `seed + c - 1`.
#' @param saveAt Directory for BGLR binary output. Defaults to a fresh temp dir.
#' @param verbose Logical; print BGLR progress.
#'
#' @return An object of class `breedRB_fit`.
#' @examples
#' \donttest{
#' G <- readRDS(system.file("extdata", "kinship.matrix.rds", package = "breedRBayes"))
#' dat <- read.csv(system.file("extdata", "data_soy.csv", package = "breedRBayes"))
#' fit <- bbglr(yield ~ 1 + env, random = ~ vm(gen, G), data = dat,
#'              relmat = list(G = G), nIter = 2000, burnIn = 500, nChains = 2)
#' }
#' @export
bbglr <- function(fixed, random = NULL, residual = NULL, data,
                  relmat = list(), nIter = 5000, burnIn = 1000, thin = 5,
                  nChains = 1, seed = 123, saveAt = NULL, verbose = TRUE) {

  cl <- match.call()
  RhpcBLASctl::blas_set_num_threads(1)
  data <- droplevels(as.data.frame(data))

  parsed <- parse_model(fixed, random, residual)
  data   <- .filter_to_relmat(parsed, data, relmat)
  built  <- build_eta(parsed, data, relmat)

  # response
  if (parsed$response$multitrait) {
    y <- as.matrix(data[, parsed$response$traits, drop = FALSE])
    # Supply the factor-analytic loadings pattern (lower-triangular identifiability
    # constraint) for any FA term now that the number of traits is known.
    nt <- ncol(y)
    built$ETA <- lapply(built$ETA, function(e) {
      if (!is.null(e$Cov) && identical(e$Cov$type, "FA")) {
        nF <- e$Cov$nF
        if (nF > nt) stop("fa(): number of factors (", nF, ") exceeds number of traits (",
                          nt, ").", call. = FALSE)
        M <- matrix(TRUE, nt, nF)
        M[upper.tri(M)] <- FALSE        # loading estimated only for factor j <= trait i
        e$Cov$M <- M
      }
      e
    })
  } else {
    if (any(vapply(built$ETA, function(e) !is.null(e$Cov), logical(1)))) {
      stop("fa()/rrc() factor-analytic structure requires a multi-trait model ",
           "(use cbind(...) on the left-hand side).", call. = FALSE)
    }
    y <- data[[parsed$response$traits]]
  }

  if (is.null(saveAt)) saveAt <- tempfile(pattern = "breedRBayes_")
  dir.create(saveAt, showWarnings = FALSE, recursive = TRUE)

  chains <- vector("list", nChains)
  paths  <- character(nChains)
  for (ch in seq_len(nChains)) {
    set.seed(seed + ch - 1L)
    prefix <- file.path(saveAt, paste0("chain", ch, "_"))
    paths[ch] <- prefix
    if (parsed$response$multitrait) {
      chains[[ch]] <- BGLR::Multitrait(
        y = y, ETA = built$ETA, nIter = nIter, burnIn = burnIn, thin = thin,
        saveAt = prefix, verbose = verbose)
    } else {
      args <- list(y = y, ETA = built$ETA, nIter = nIter, burnIn = burnIn,
                   thin = thin, saveAt = prefix, verbose = verbose)
      if (!is.null(built$group)) args$groups <- built$group
      chains[[ch]] <- do.call(BGLR::BGLR, args)
    }
  }

  structure(
    list(
      call     = cl,
      parsed   = parsed,
      meta     = built$meta,
      group    = built$group,
      chains   = chains,
      paths    = paths,
      data     = data,
      relmat   = relmat,
      response = parsed$response,
      control  = list(nIter = nIter, burnIn = burnIn, thin = thin,
                      nChains = nChains, seed = seed, saveAt = saveAt)
    ),
    class = "breedRB_fit"
  )
}

#' @export
print.breedRB_fit <- function(x, ...) {
  cat("<breedRB_fit>\n")
  cat("  Call: ", deparse(x$call), "\n", sep = "")
  cat("  Traits: ", paste(x$response$traits, collapse = ", "),
      if (x$response$multitrait) " (multi-trait)" else "", "\n", sep = "")
  cat("  Fixed terms:  ", if (length(x$parsed$fixed)) paste(vapply(x$parsed$fixed,
      `[[`, character(1), "label"), collapse = ", ") else "(intercept only)", "\n", sep = "")
  cat("  Random terms: ", if (length(x$parsed$random)) paste(vapply(x$parsed$random,
      `[[`, character(1), "label"), collapse = ", ") else "(none)", "\n", sep = "")
  cat("  Residual: ", if (isTRUE(x$parsed$residual$hetero))
      paste0("heterogeneous by ", x$parsed$residual$group) else "homogeneous", "\n", sep = "")
  cat("  Chains: ", x$control$nChains, "  nIter: ", x$control$nIter,
      "  burnIn: ", x$control$burnIn, "  thin: ", x$control$thin, "\n", sep = "")
  invisible(x)
}
