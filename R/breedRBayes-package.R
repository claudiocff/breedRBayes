#' breedRBayes: ASReml-style formula front-end for the BGLR Bayesian engine
#'
#' Fit Bayesian mixed models with [BGLR::BGLR()] / [BGLR::Multitrait()] using a
#' formula interface inspired by \pkg{asreml} and \pkg{sommer}. Models are
#' described with separate `fixed`, `random` and `residual` formulas and a set of
#' in-formula "special" functions ([vm()], [leg()], [fa()], ...). The package
#' supports genomic prediction (GBLUP), random-regression / reaction-norm,
#' multi-trait and factor-analytic models, and returns posterior chains,
#' variance components, heritabilities (as full posterior distributions) and
#' Markov-chain convergence diagnostics.
#'
#' @section Main entry point:
#' [bbglr()] fits a model. See [varcomp()], [heritability()], [mcmc_diag()] and
#' [solution()] for post-processing.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats terms as.formula model.matrix sd median quantile var mad pchisq
#'   reformulate prcomp
#' @importFrom utils read.table
## usethis namespace: end
NULL
