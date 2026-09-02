# breedRBayes

**An ASReml-style formula front-end for the [BGLR](https://cran.r-project.org/package=BGLR) Bayesian engine.**

`breedRBayes` lets you describe Bayesian mixed models with the familiar
`fixed` / `random` / `residual` triple-formula syntax of ASReml-R and `sommer`,
and fits them with the Gibbs sampler in **BGLR**. It targets plant- and
animal-breeding genetic evaluation and exposes the full Bayesian output:
posterior chains, variance components, heritabilities as complete posterior
distributions, genetic correlations, and Markov-chain convergence diagnostics.

## Features

- **One formula interface, four model classes**
  - **GBLUP / fixed / random** — genomic relationship via `vm(gen, G)`.
  - **Random regression / reaction norm** — `leg(x, order)` / `rr(x, order)`
    Legendre bases.
  - **Multi-trait** — `cbind(t1, t2, ...) ~ ...` with unstructured genetic
    covariance and genetic correlations.
  - **Factor-analytic** — `fa(gen, k)` / `rrc(gen, k)` reduced-rank genetic
    covariance.
- **Multiple chains** for Gelman–Rubin R-hat and other diagnostics.
- **Posterior everything** — `varcomp()`, `heritability()` and `gebv()` return
  full posterior draws, not just point estimates.
- **Diagnostics** — `mcmc_diag()` (R-hat, effective sample size, Geweke) plus
  trace and posterior-density plots.
- **Finlay–Wilkinson** two-step reaction-norm analysis (`FW_bglr()`) built on
  top of the same engine, with a value-added probability-of-superiority layer.

## Installation

```r
# install.packages("remotes")
remotes::install_github("claudiocff/breedRBayes")
```

`breedRBayes` depends on `BGLR`, `coda`, `RhpcBLASctl`, `dplyr`, `tidyr`,
`reshape2`, `tibble` and `ggplot2`.

## The formula interface

Models are described with up to three formulas:

| Argument   | Purpose                          | Example                              |
|------------|----------------------------------|--------------------------------------|
| `fixed`    | mean / fixed effects (LHS = trait) | `yield ~ 1 + env`                  |
| `random`   | random effects                   | `~ vm(gen, G) + env:vm(gen, G)`      |
| `residual` | residual variance structure      | `~ dsum(~units \| env)` (heterogeneous) |

Special functions inside the formulas:

| Special                | Meaning                                              |
|------------------------|------------------------------------------------------|
| `vm(f, K)`             | random effect of factor `f` with covariance `K` (GBLUP). `K` is the **covariance** (relationship) matrix, not its inverse. |
| `leg(x, n)` / `rr(x, n)` | random regression on an orthonormal Legendre basis of order `n` |
| `fa(f, k)` / `rrc(f, k)` | factor-analytic genetic covariance with `k` factors (multi-trait) |
| `a:b`                  | interaction (e.g. `env:vm(gen, G)` for GxE)          |
| `cbind(t1, t2) ~ ...`  | stack traits → multi-trait model                     |

Relationship matrices are supplied through `relmat = list(G = K)`.

## Quick start

```r
library(breedRBayes)

G   <- readRDS(system.file("extdata", "kinship.matrix.rds", package = "breedRBayes"))
dat <- read.csv(system.file("extdata", "data_soy.csv",      package = "breedRBayes"))

## --- GBLUP with two chains --------------------------------------------------
fit <- bbglr(
  yield ~ 1 + env,
  random = ~ vm(gen, G),
  data   = dat,
  relmat = list(G = G),
  nIter  = 6000, burnIn = 2000, nChains = 2
)

varcomp(fit)          # posterior variance components
heritability(fit)     # h2 as a full posterior distribution
mcmc_diag(fit)        # R-hat, effective sample size, Geweke
gebv(fit)             # posterior genomic breeding values

plot_trace(fit)       # chain trace plots
plot_posterior(heritability(fit))
```

### Random regression (reaction norm)

```r
fit_rr <- bbglr(
  yield ~ 1 + env,
  random = ~ vm(gen, G) + leg(env_index, 2):vm(gen, G),
  data = dat, relmat = list(G = G),
  nIter = 6000, burnIn = 2000
)
```

### Multi-trait and genetic correlations

```r
fit_mt <- bbglr(
  cbind(t1, t2, t3) ~ env,
  random = ~ vm(gen, G),
  data = dat, relmat = list(G = G)
)
varcomp(fit_mt)       # per-trait genetic/residual variances + rg for each pair
heritability(fit_mt)  # per-trait h2
```

### Factor-analytic

```r
fit_fa <- bbglr(
  cbind(t1, t2, t3) ~ env,
  random = ~ fa(gen, 1),      # 1-factor genetic covariance
  data = dat, relmat = list(G = G)
)
```

## Finlay–Wilkinson reaction norm

`FW_bglr()` runs the two-step Finlay–Wilkinson analysis (estimate the
environment index, then a genotype-specific intercept + slope random
regression) using `bbglr()` as the engine:

```r
fw  <- FW_bglr(dat, gen = "gen", env = "env", trait = "yield",
               kinship.matrix = G, nIter = 3000, burnIn = 1000, order = 1)
ext <- extr_outs_bglr(fw, gen = "gen")   # posterior chains, variances, G params
ps  <- prob_sup_bglr(ext, int = 0.2)     # probability-of-superiority summaries
```

## Example data

Example datasets ship in `inst/extdata/` (soybean phenotypes, a genomic
relationship matrix, and environmental-covariate matrices). Reach them with
`system.file("extdata", "<file>", package = "breedRBayes")`. Runnable end-to-end
scripts live in `inst/scripts/`.

## Documentation

Every exported function is documented; see `?bbglr`, `?varcomp`,
`?heritability`, `?mcmc_diag`, `?gebv` and `?FW_bglr`. A full PDF reference
manual (`breedRBayes-manual.pdf`) is included at the repository root.

## License

MIT © breedRBayes authors. See [LICENSE](LICENSE).
