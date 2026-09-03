# breedRBayes — end-to-end examples
#
# Run from a session with the package loaded (devtools::load_all() or
# library(breedRBayes)). Example data ships in inst/extdata/.

library(breedRBayes)

## --- Data ------------------------------------------------------------------
G   <- readRDS(system.file("extdata", "kinship.matrix.rds", package = "breedRBayes"))
dat <- read.csv(system.file("extdata", "data_soy.csv",      package = "breedRBayes"))

str(dat)
length(unique(dat$env))
length(unique(dat$gen))

## --- 1. GBLUP with two chains ----------------------------------------------
fit <- bbglr(
  yield ~ 1 + env,
  random = ~ vm(gen, G),
  data   = dat,
  relmat = list(G = G),
  nIter  = 6000, burnIn = 2000, nChains = 2
)

varcomp(fit)                       # posterior variance components
heritability(fit)                  # h2 as a full posterior distribution
mcmc_diag(fit, plot = TRUE)        # R-hat, ESS, Geweke + ggplot2 trace plots

solution(fit, term = "gen")        # posterior genomic breeding values
solution(fit, term = "env", type = "fixed")   # fixed-effect estimates

plot_trace(fit)                    # chain trace plots
plot_posterior(heritability(fit))  # heritability posterior density

## --- 2. Random effect without a relationship matrix ------------------------
fit_r <- bbglr(yield ~ env + env:rep, random = ~ gen, residual = ~ units,
               data = dat)
solution(fit_r, term = "gen", type = "random")   # BLUPs, no G matrix needed

## --- 3. Random regression (reaction norm) ----------------------------------
## Needs an environmental covariate column (e.g. env_index) in `dat`.
# fit_rr <- bbglr(
#   yield ~ 1 + env,
#   random = ~ vm(gen, G) + leg(env_index, 2):vm(gen, G),
#   data = dat, relmat = list(G = G),
#   nIter = 6000, burnIn = 2000
# )

## --- 4. Multi-trait and factor-analytic ------------------------------------
# fit_mt <- bbglr(cbind(t1, t2, t3) ~ env, random = ~ vm(gen, G),
#                 data = dat, relmat = list(G = G))
# varcomp(fit_mt)      # per-trait variances + genetic correlations
# heritability(fit_mt) # per-trait h2
#
# fit_fa <- bbglr(cbind(t1, t2, t3) ~ env, random = ~ fa(gen, 1),
#                 data = dat, relmat = list(G = G))
