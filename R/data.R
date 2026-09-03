# ---------------------------------------------------------------------------
# Package example data: the CIMMYT wheat genomic-prediction dataset.
# ---------------------------------------------------------------------------

#' CIMMYT wheat multi-environment yield trial
#'
#' Grain yield of 599 CIMMYT wheat lines evaluated in four environments, in
#' long format ready for the \pkg{breedRBayes} formula interface. This is the
#' long-format reshape of the classic `wheat` dataset distributed with
#' \pkg{BGLR} (`wheat.Y`); the yields are the pre-standardised values from that
#' source. The matching marker matrix is [wheat_M] and the two share the same
#' genotype identifiers, so the same data serve the phenotype-only
#' (`random = ~ gen`) and genomic (`random = ~ mrk(gen, wheat_M)`) analyses.
#'
#' @format A data frame with 2396 rows (599 lines x 4 environments) and 3
#'   columns:
#' \describe{
#'   \item{gen}{Line identifier (factor, 599 levels), matching `rownames(wheat_M)`.}
#'   \item{env}{Environment (factor, `E1`..`E4`).}
#'   \item{yield}{Standardised grain yield (one record per line x environment).}
#' }
#' @source Distributed with the \pkg{BGLR} package; originally from CIMMYT's
#'   Global Wheat Program. See Crossa et al. (2010) \emph{Genetics} 186:713--724.
#' @seealso [wheat_M] for the marker matrix.
#' @examples
#' \donttest{
#' data(wheat)
#' fit <- bbglr(yield ~ env, random = ~ gen, data = wheat)
#' heritability(fit, genetic = "gen")
#' }
"wheat"

#' CIMMYT wheat marker matrix
#'
#' Molecular markers for the 599 wheat lines in [wheat]: 1279 DArT markers coded
#' `0/1`. Row names are the line identifiers (the levels of `wheat$gen`), so the
#' matrix can be passed straight to `mrk(gen, wheat_M)` for a genomic (GBLUP /
#' RRBLUP) analysis.
#'
#' @format A numeric matrix, 599 lines (rows) x 1279 markers (columns), entries
#'   in `{0, 1}`. Row names are the line identifiers; column names are `m1`..`m1279`.
#' @source Distributed with the \pkg{BGLR} package (`wheat.X`); CIMMYT Global
#'   Wheat Program. See Crossa et al. (2010) \emph{Genetics} 186:713--724.
#' @seealso [wheat] for the phenotypes.
#' @examples
#' \donttest{
#' data(wheat)
#' fit_g <- bbglr(yield ~ env, random = ~ mrk(gen, wheat_M), data = wheat)
#' }
"wheat_M"
