# ---------------------------------------------------------------------------
# Formula parsing: turn asreml-style fixed/random/residual formulas into an
# internal, engine-agnostic list of "term specs" consumed by build_eta().
#
# A term spec has the shape:
#   list(label      = <original term label>,
#        role       = "fixed" | "random",
#        components = list(<component spec>, ...))   # interaction components
#
# A component spec is one of:
#   list(type = "factor",    var = <name>, center = <lgl>)
#   list(type = "vm",        var = <name>, relmat = <key>)   # genomic / GBLUP
#   list(type = "mrk",       var = <name>, relmat = <key>, method = <chr>, rank = <int>)
#                                                            # markers, auto GBLUP/RRBLUP
#   list(type = "leg",       var = <name>, order = <int>)    # random regression
#   list(type = "covariate", var = <name>)                   # numeric covariate
# ---------------------------------------------------------------------------

#' Flatten an interaction expression (`a:b:c`) into its components
#' @keywords internal
.flatten_interaction <- function(e) {
  if (is.call(e) && identical(e[[1]], as.name(":"))) {
    c(.flatten_interaction(e[[2]]), .flatten_interaction(e[[3]]))
  } else {
    list(e)
  }
}

#' Classify a single interaction component (symbol or special call)
#' @keywords internal
.classify_component <- function(comp) {
  if (is.symbol(comp)) {
    return(list(type = "factor", var = as.character(comp)))
  }
  if (!is.call(comp)) {
    stop("Unsupported term component: ", deparse(comp), call. = FALSE)
  }
  fname <- as.character(comp[[1]])
  args  <- as.list(comp)[-1]
  argnm <- function(i, default = NULL) {
    if (length(args) >= i) as.character(args[[i]]) else default
  }
  switch(fname,
    vm = list(type = "vm", var = argnm(1),
              relmat = if (length(args) >= 2) as.character(args[[2]]) else NULL),
    mrk = {
      # match named or positional args against mrk(var, M, method, rank)
      mc <- match.call(function(var, M, method = "auto", rank = NULL) NULL, comp)
      list(type = "mrk", var = as.character(mc$var),
           relmat = if (!is.null(mc$M)) as.character(mc$M) else NULL,
           method = if (!is.null(mc$method))
                      match.arg(as.character(mc$method), c("auto", "GBLUP", "RRBLUP"))
                    else "auto",
           rank = if (!is.null(mc$rank)) as.integer(eval(mc$rank)) else NA_integer_)
    },
    leg = ,
    rr  = list(type = "leg", var = argnm(1),
               order = if (length(args) >= 2) as.integer(eval(args[[2]])) else 1L),
    fa  = ,
    rrc = list(type = "fa", var = argnm(1),
               nfac = if (length(args) >= 2) as.integer(eval(args[[2]])) else 1L),
    cov = ,
    lin = list(type = "covariate", var = argnm(1)),
    stop("Unknown special function '", fname, "()' in model formula.", call. = FALSE)
  )
}

#' Parse a single term label (e.g. "env:vm(gen, G)") into a term spec
#' @keywords internal
.parse_term_label <- function(label, role) {
  e <- str2lang(label)
  comps <- lapply(.flatten_interaction(e), .classify_component)
  list(label = label, role = role, components = comps)
}

#' Parse one side (`fixed`/`random`) formula into a list of term specs
#'
#' @return A list with `intercept` (logical) and `terms` (list of term specs).
#' @keywords internal
.parse_side <- function(formula, role) {
  if (is.null(formula)) return(list(intercept = FALSE, terms = list()))
  tt <- terms(formula, keep.order = TRUE)
  labels <- attr(tt, "term.labels")
  intercept <- attr(tt, "intercept") == 1
  terms_list <- lapply(labels, .parse_term_label, role = role)
  list(intercept = intercept, terms = terms_list)
}

#' Parse the response (LHS) of the fixed formula
#'
#' Supports a single trait (`yield ~ ...`) or several traits stacked with
#' `cbind(t1, t2, ...) ~ ...` (triggers a multi-trait fit).
#' @keywords internal
.parse_response <- function(fixed) {
  if (length(fixed) < 3L) stop("`fixed` must be a two-sided formula.", call. = FALSE)
  lhs <- fixed[[2]]
  if (is.call(lhs) && identical(lhs[[1]], as.name("cbind"))) {
    traits <- vapply(as.list(lhs)[-1], as.character, character(1))
    list(traits = traits, multitrait = TRUE)
  } else {
    list(traits = as.character(lhs), multitrait = FALSE)
  }
}

#' Parse a residual specification
#'
#' `~ units` (or `NULL`) means a single homogeneous residual variance.
#' `~ dsum(~units | env)` (or `~ diag(env)`) means heterogeneous residual
#' variances grouped by `env`, mapped to BGLR's `groups=` argument.
#' @keywords internal
.parse_residual <- function(residual) {
  if (is.null(residual)) return(list(hetero = FALSE, group = NULL))
  rhs <- residual[[length(residual)]]
  if (is.call(rhs)) {
    fname <- as.character(rhs[[1]])
    if (fname == "dsum") {
      # dsum(~units | env)  ->  extract grouping var after '|'
      inner <- rhs[[2]]
      cond  <- inner[[length(inner)]]
      if (is.call(cond) && identical(cond[[1]], as.name("|"))) {
        return(list(hetero = TRUE, group = as.character(cond[[3]])))
      }
    }
    if (fname %in% c("diag", "dsum")) {
      return(list(hetero = TRUE, group = as.character(rhs[[2]])))
    }
  }
  list(hetero = FALSE, group = NULL)
}

#' Parse asreml-style fixed / random / residual formulas
#'
#' Turns the three model formulas into an internal, engine-agnostic description
#' used by [bbglr()]. Not usually called directly.
#'
#' @param fixed Two-sided formula for the mean/fixed part, e.g.
#'   `yield ~ 1 + env`. Multi-trait models use `cbind(t1, t2) ~ ...`.
#' @param random One-sided formula for random terms, e.g.
#'   `~ vm(gen, G) + env:vm(gen, G)`. `NULL` for none.
#' @param residual One-sided residual formula. `NULL` (default) or `~ units`
#'   gives homogeneous residuals; `~ dsum(~units | env)` gives heterogeneous
#'   residual variances by `env`.
#'
#' @return A `breedRB_terms` list: `response`, `fixed`, `random`, `residual`.
#' @examples
#' parse_model(yield ~ 1 + env, ~ vm(gen, G))
#' @export
parse_model <- function(fixed, random = NULL, residual = NULL) {
  response <- .parse_response(fixed)
  fixed_side  <- .parse_side(fixed, "fixed")
  random_side <- .parse_side(random, "random")
  structure(
    list(
      response  = response,
      intercept = fixed_side$intercept,
      fixed     = fixed_side$terms,
      random    = random_side$terms,
      residual  = .parse_residual(residual)
    ),
    class = "breedRB_terms"
  )
}
