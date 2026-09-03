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
  - **GBLUP / fixed / random** — a raw marker matrix via `mrk(gen, M)` (auto
    GBLUP/RR-BLUP; self-scaling low-rank PC-GBLUP via `exp_var_rank`;
    `predict()` / `predict_pr()` for new genotypes, `solve_SNP()` for effects),
    or a ready-made relationship matrix via `vm(gen, G)`.
  - **Random regression / reaction norm** — `leg(x, order)` / `rr(x, order)`
    Legendre bases.
  - **Multi-trait** — `cbind(t1, t2, ...) ~ ...` with unstructured genetic
    covariance and genetic correlations.
  - **Factor-analytic** — `fa(gen, k)` / `rrc(gen, k)` reduced-rank genetic
    covariance.
- **One control object** — MCMC settings and model hyperparameters (including the
  genomic low-rank target) live in `bbglr_control()`, keeping the main call short.
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

Marker matrices (`mrk()`) and relationship/kernel matrices (`vm()`) are both
supplied through `relmat = list(M = ..., G = ...)`.

## Quick start

```r
library(breedRBayes)

## A small genomic dataset: 150 genotypes x 500 markers, 130 of them phenotyped.
set.seed(1)
n <- 150; p <- 500
M <- matrix(rbinom(n * p, 2, 0.3), n, p,
            dimnames = list(paste0("G", 1:n), paste0("m", 1:p)))
train <- paste0("G", 1:130)          # phenotyped genotypes
newg  <- paste0("G", 131:150)        # genotyped but NOT phenotyped

g   <- scale(scale(M[train, ], scale = FALSE) %*% rnorm(p)) * 3   # true breeding values
dat <- data.frame(gen = rep(train, each = 3))
dat$yield <- g[dat$gen, 1] + rnorm(nrow(dat), 0, 3)

## --- genomic model with two chains -----------------------------------------
## mrk() takes the marker matrix directly and auto-selects GBLUP or RR-BLUP.
fit <- bbglr(
  yield ~ 1,
  random  = ~ mrk(gen, M),
  data    = dat,
  relmat  = list(M = M[train, ]),
  control = bbglr_control(nIter = 6000, burnIn = 2000, nChains = 2)
)

varcomp(fit)                          # posterior variance components
heritability(fit)                     # h2 as a full posterior distribution
mcmc_diag(fit, plot = TRUE)           # R-hat, ESS, Geweke + ggplot2 trace plots
solution(fit, term = "mrk(gen, M)")   # posterior genomic breeding values

## Probability each genotype ranks in the top 10% (Bayesian selection):
pr(fit, term = "mrk(gen, M)", threshold = 0.10)

## Predict the 20 genotyped-but-unphenotyped lines from their markers:
predict(fit, M[newg, ], add_mu = TRUE)   # ID, prediction, sd, lower, upper

## P(each new line beats the top 10% of the training population):
predict_pr(fit, M[newg, ], threshold = 0.10)   # ID, prediction, sd, prob

## Marker effects (back-solved from a GBLUP fit, read directly from RR-BLUP):
solve_SNP(fit)

plot_trace(fit)       # chain trace plots
plot_posterior(heritability(fit))
```

### Controlling the fit

MCMC settings and model hyperparameters are bundled in `bbglr_control()`, so the
main `bbglr()` call stays focused on the model:

```r
ctrl <- bbglr_control(
  nIter        = 20000,   # total iterations
  burnIn       = 5000,    # burn-in
  thin         = 10,      # thinning interval
  nChains      = 4,       # independent chains (enables R-hat)
  seed         = 1,       # base seed (chain c uses seed + c - 1)
  exp_var_rank = 0.99     # genomic low-rank target (fraction of variance kept)
)
fit <- bbglr(yield ~ 1, random = ~ mrk(gen, M), data = dat, relmat = list(M = M),
             control = ctrl)
```

The individual settings may also be passed straight to `bbglr()` as a shortcut
(`bbglr(..., nIter = 20000)`); they override the corresponding value in
`control`. `exp_var_rank` is described under low-rank GBLUP below.

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

### Genomic models from a marker matrix (recommended)

`mrk(gen, M)` is the recommended way to fit a genomic model: give it the raw
**marker matrix** (genotypes in rows, markers in columns) and it builds the model
directly, picking the cheaper of the two mathematically equivalent
parameterizations:

- **markers ≥ genotypes** → **GBLUP** (an `n × n` genomic relationship built
  internally by **VanRaden's method 1**, `G = ZZ'/2Σpⱼqⱼ`, then fitted in its
  principal-component basis; a small diagonal ridge keeps `G` non-singular),
- **genotypes > markers** → **RR-BLUP** (estimate the `p` marker effects directly).

Both give the same breeding values and the same heritability; only the compute
cost differs. Pass `method = "GBLUP"` or `"RRBLUP"` to force one. The marker
matrix should be coded as **0/1/2** allele dosages; missing values are
mean-imputed per marker and near-constant markers are dropped automatically.

```r
fit <- bbglr(yield ~ 1 + env, random = ~ mrk(gen, M), data = dat, relmat = list(M = M))
solution(fit, term = "mrk(gen, M)", type = "random")   # genomic breeding values
```

**Large genotype panels — low-rank GBLUP.** For many genotypes the `n × n`
relationship and its eigendecomposition dominate the cost, so GBLUP is fitted in
a **low-rank principal-component basis**. By default the rank is chosen to
explain 99% of the genomic variance (`bbglr_control(exp_var_rank = 0.99)`): the
long tail of numerically-tiny components — which carry almost no genetic signal —
is dropped, shrinking the design and speeding up sampling while retaining
essentially all of the information. Lower the target for a more aggressive
approximation, or set `exp_var_rank = NA` to keep every component:

```r
# self-scaling default (99% of variance); tune it through the control object
fit <- bbglr(yield ~ 1, random = ~ mrk(gen, M), data = dat, relmat = list(M = M),
             control = bbglr_control(exp_var_rank = 0.95))
```

The `exp_var_rank` search does **not** need the full spectrum. The total genomic
variance is the trace of the relationship — `Σ rowSums(Mc²)`, free to compute —
so it is the exact denominator up front; with **RSpectra** installed the leading
components are then grown incrementally (`svds`) and the search stops as soon as
their cumulative share crosses the target. For a low-rank panel (a handful of PCs
dominate, as in real genomic data) only a small `k ≪ n` is ever computed and the
`n × n` matrix is *never formed*; it falls back to a single full
eigendecomposition only when the target would need most of the spectrum.

`predict()` and `solve_SNP()` reuse the same cached rotation, so scoring and
marker back-solving stay fast.

Because the marker matrix is retained, a `mrk()` fit can **predict new
genotypes** and **recover marker effects** — a `vm()` fit (below) cannot.

**Predict unobserved genotypes** with `predict()`: pass a marker matrix `M_new`
for genotypes that were not in the training data. Their value is
`Mc_new %*% b` (new markers centred by the *training* means, times the posterior
marker effects), with the full posterior propagated to a credible interval.
`add_mu = TRUE` (the default) puts predictions on the mean/phenotype scale:

```r
pred <- predict(fit, M_new, add_mu = TRUE)   # score genotypes never seen in training
head(pred)                                   # ID, prediction, sd, lower, upper
```

**Probability a new line clears a population bar** with `predict_pr()`: the
posterior probability that a new genotype exceeds a threshold defined by the
**training population** (e.g. its top 10%). It is computed per posterior draw —
in each draw the bar is that quantile of the *training* genotypes and the new
line's predicted value is compared against it — so it uses the whole posterior,
not just the point estimate. A new line with **more marker information** (a
tighter predictive posterior) that sits above the bar earns a higher probability
than an equally-ranked but more uncertain one:

```r
predict_pr(fit, M_new, threshold = 0.10)                 # P(above the training top 10%)
predict_pr(fit, M_new, threshold = 0.10, higher = FALSE) # P(below the training bottom 10%)
```

**Marker effects** with `solve_SNP()`: per-marker allele-substitution effects on
the centred-marker scale (`GEBV = Mc %*% b`). RR-BLUP effects are read directly;
GBLUP effects are back-solved (`b = Mcᵀ (Mc Mcᵀ)⁻¹ u`) per posterior draw, reusing
the training markers held in the fit:

```r
snp <- solve_SNP(fit)                     # no need to re-supply the marker matrix
head(snp[order(-abs(snp$effect)), ])      # largest-effect markers
```

### `vm(f, K)`: bring your own relationship / kernel matrix

Use `vm(f, K)` when you already have a **covariance matrix** and not the raw
markers. `K` can be any positive (semi-)definite kernel indexed by the levels of
`f`:

- a **pedigree** numerator relationship matrix (A),
- an **environmic** kernel (environmental covariance between locations/years,
  e.g. from `EnvRtype`),
- a genomic `G` computed elsewhere (custom method / allele frequencies),
- any other custom covariance kernel.

It fits the same reaction-norm / GBLUP-style model as `mrk()` in its GBLUP mode,
but since it holds no markers it cannot `predict()` new levels or `solve_SNP()`.
For raw marker data, prefer `mrk()`.

```r
fit <- bbglr(yield ~ 1 + env, random = ~ vm(gen, A), data = dat, relmat = list(A = A))  # pedigree
fit <- bbglr(yield ~ 1,       random = ~ vm(env, W), data = dat, relmat = list(W = W))  # environmic
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
