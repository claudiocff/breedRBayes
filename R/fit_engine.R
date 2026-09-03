# ---------------------------------------------------------------------------
# bbglr(): main entry point. Parses formulas, builds the BGLR ETA, and runs one
# or more MCMC chains, returning a `breedRB_fit` object.
# ---------------------------------------------------------------------------

#' Fit a Bayesian mixed model with the BGLR engine using asreml-style formulas
#'
#' @param fixed Two-sided formula for the fixed/mean part, e.g. `yield ~ 1 + env`.
#'   Use `cbind(t1, t2) ~ ...` for a multi-trait model.
#' @param random One-sided formula of random terms, e.g. `~ mrk(gen, M)`. Special
#'   functions: `mrk(f, M, method, rank)` genomic effect from a **marker matrix**
#'   `M` that automatically fits GBLUP or RR-BLUP, whichever is cheaper (`method`
#'   is `"auto"`, `"GBLUP"` or `"RRBLUP"`; `"auto"` picks GBLUP when markers >=
#'   genotypes, else RR-BLUP — the two are prediction-equivalent). `rank` (GBLUP
#'   only, default all) keeps just the top-`rank` principal components of the
#'   genomic relationship — a low-rank / PC-GBLUP approximation that shrinks the
#'   design and, when \pkg{RSpectra} is installed, is computed straight from the
#'   markers without ever forming the `n x n` relationship matrix (much cheaper
#'   for many genotypes). A `mrk()` fit
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
#'   genotype IDs) and markers in columns, coded as **0/1/2 allele dosages**.
#'   Missing values are mean-imputed per marker and (near-)constant markers are
#'   dropped automatically. In GBLUP mode `mrk()` then builds the genomic
#'   relationship internally by VanRaden's method 1 (\eqn{G = ZZ'/2\sum p_j q_j})
#'   and adds a tiny diagonal ridge so `G` is always non-singular. For `vm()` the entry is a covariance / kernel matrix
#'   (K, not its inverse) — genomic, pedigree, environmic or any custom kernel —
#'   with row/column names matching the factor levels.
#' @param control A [bbglr_control()] object bundling the MCMC controls (`nIter`,
#'   `burnIn`, `thin`, `nChains`, `seed`, `saveAt`, `verbose`) and model
#'   hyperparameters (`exp_var_rank`, the explained-variance target for the
#'   genomic low-rank rotation).
#' @param ... Individual [bbglr_control()] settings passed directly as a shortcut
#'   (e.g. `nIter = 20000`); they override the corresponding value in `control`.
#'
#' @return An object of class `breedRB_fit`.
#' @examples
#' \donttest{
#' G <- readRDS(system.file("extdata", "kinship.matrix.rds", package = "breedRBayes"))
#' dat <- read.csv(system.file("extdata", "data_soy.csv", package = "breedRBayes"))
#' fit <- bbglr(yield ~ 1 + env, random = ~ vm(gen, G), data = dat,
#'              relmat = list(G = G),
#'              control = bbglr_control(nIter = 2000, burnIn = 500, nChains = 2))
#' }
#' @seealso [bbglr_control()] for the available control settings.
#' @export
bbglr <- function(fixed, random = NULL, residual = NULL, data,
                  relmat = list(), control = bbglr_control(), ...) {

  cl <- match.call()

  # Accept individual control settings passed directly (back-compat + shortcut),
  # merging them over the supplied/defaulted control object.
  dots <- list(...)
  if (length(dots)) {
    ctrl_args <- names(formals(bbglr_control))
    bad <- setdiff(names(dots), ctrl_args)
    if (length(bad)) {
      stop("Unknown argument(s) to bbglr(): ", paste(bad, collapse = ", "),
           ". Set model controls via control = bbglr_control(...).", call. = FALSE)
    }
    control <- do.call(bbglr_control,
                       utils::modifyList(unclass(control)[ctrl_args], dots))
  }
  if (!inherits(control, "breedRB_control")) {
    control <- do.call(bbglr_control, as.list(control))
  }

  nIter   <- control$nIter;   burnIn <- control$burnIn; thin    <- control$thin
  nChains <- control$nChains; seed   <- control$seed;   verbose <- control$verbose
  saveAt  <- control$saveAt
  exp_var <- if (is.na(control$exp_var_rank)) NA_real_ else control$exp_var_rank

  RhpcBLASctl::blas_set_num_threads(1)
  data <- droplevels(as.data.frame(data))

  parsed <- parse_model(fixed, random, residual)
  data   <- .filter_to_relmat(parsed, data, relmat)
  built  <- build_eta(parsed, data, relmat, exp_var = exp_var)

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
  control$saveAt <- saveAt                       # record the resolved output directory

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
      control  = control
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
