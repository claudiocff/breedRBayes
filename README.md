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
  - **GBLUP / fixed / random** — genomic relationship via `vm(gen, G)`, or a raw
    marker matrix via `mrk(gen, M)` (auto GBLUP/RR-BLUP; `solve_SNP()` for effects).
  - **Random regression / reaction norm** — `leg(x, order)` / `rr(x, order)`
    Legendre bases.
  - **Multi-trait** — `cbind(t1, t2, ...) ~ ...` with unstructured genetic
    covariance and genetic correlations.
  - **Factor-analytic** — `fa(gen, k)` / `rrc(gen, k)` reduced-rank genetic
    covariance.
- **Multiple chains** for Gelman–Rubin R-hat and other diagnostics.
- **Posterior everything** — `varcomp()`, `heritability()` and `solution()` return
  full posterior draws, not just point estimates.
- **Diagnostics** — `mcmc_diag()` (R-hat, effective sample size, Geweke) plus
  trace and posterior-density plots.

## Installation

```r
# install.packages("remotes")
remotes::install_github("claudiocff/breedRBayes")
```

`breedRBayes` depends on `BGLR`, `coda`, `RhpcBLASctl` and `ggplot2`.

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
| `mrk(f, M, method)`    | genomic effect of `f` from a **marker matrix** `M` (genotypes in rows). Auto-selects GBLUP or RR-BLUP (`method = "auto"`, `"GBLUP"`, `"RRBLUP"`). |
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
mcmc_diag(fit, plot = TRUE)          # R-hat, ESS, Geweke + ggplot2 trace plots
solution(fit, term = "gen")          # posterior genomic breeding values

plot_trace(fit)       # chain trace plots
plot_posterior(heritability(fit))
```

### Extracting solutions for any term

`solution()` extracts posterior solutions for any term — genomic `vm()` BLUPs, a
plain random factor (fitted without a relationship matrix), random-regression
coefficients, or fixed effects:

```r
fit <- bbglr(yield ~ env + env:rep, random = ~ gen, residual = ~ units, data = dat)

solution(fit, term = "gen", type = "random")   # random BLUPs, no G matrix needed
solution(fit, term = "env", type = "fixed")    # fixed-effect estimates
```

`type` is optional and, when given, is checked against the term's actual role.
Set `add_mu = TRUE` to add the model intercept to every draw, shifting the
solutions onto the overall-mean scale (`BLUP + mu`):

```r
solution(fit, term = "gen", type = "random", add_mu = TRUE)
```

### Marker matrix: automatic GBLUP / RR-BLUP

`vm(gen, G)` fits GBLUP from a ready-made relationship matrix. When you have the
**marker matrix** instead, `mrk(gen, M)` builds the model directly and picks the
cheaper of the two mathematically equivalent parameterizations:

- **markers ≥ genotypes** → **GBLUP** (an `n × n` genomic relationship, fitted in
  its principal-component basis),
- **genotypes > markers** → **RR-BLUP** (estimate the `p` marker effects directly).

Both give the same breeding values and the same heritability; only the compute
cost differs. Pass `method = "GBLUP"` or `"RRBLUP"` to force one.

```r
fit <- bbglr(yield ~ 1 + env, random = ~ mrk(gen, M), data = dat, relmat = list(M = M))
solution(fit, term = "mrk(gen, M)", type = "random")   # genomic breeding values
```

`solve_SNP()` returns per-marker (allele-substitution) effects on the
centred-marker scale, such that `GEBV = Mc %*% b`. For an RR-BLUP fit the effects
are read directly; for a GBLUP fit they are **back-solved** from the breeding
values (`b = Mcᵀ (Mc Mcᵀ)⁻¹ u`, applied to every posterior draw). GBLUP needs the
marker matrix passed in; RR-BLUP does not:

```r
snp <- solve_SNP(fit, M)                 # GBLUP fit: M required to back-solve
head(snp[order(-abs(snp$effect)), ])     # largest-effect markers
```

### Probability of selection

`pr()` gives the posterior probability that each level ranks in the top
fraction. Every MCMC draw is ranked independently and the best
`round(threshold × n)` levels are flagged; the reported probability is how often
a level is flagged — the Bayesian "probability of being in the top X%":

```r
pr(fit, term = "gen", type = "random", threshold = 0.20)  # P(top 20%)
pr(fit, term = "gen", type = "random", threshold = 0.10, higher = FALSE)  # P(bottom 10%)
```

With `pair = TRUE` it returns the pairwise table of `P(A > B)` for every pair of
levels (`A` is the member with the larger posterior mean):

```r
pr(fit, term = "gen", type = "random", pair = TRUE)   # columns: A, B, prob = P(A > B)
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

## Example data

Example datasets ship in `inst/extdata/` (soybean phenotypes, a genomic
relationship matrix, and environmental-covariate matrices). Reach them with
`system.file("extdata", "<file>", package = "breedRBayes")`. Runnable end-to-end
scripts live in `inst/scripts/`.

## Documentation

Every exported function is documented; see `?bbglr`, `?varcomp`,
`?heritability`, `?mcmc_diag` and `?solution`. A full PDF reference
manual (`breedRBayes-manual.pdf`) is included at the repository root.

## License

MIT © breedRBayes authors. See [LICENSE](LICENSE).
