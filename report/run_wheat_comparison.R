## Fit the genomic model ladder on the CIMMYT wheat data and cache everything
## needed for the HTML report (report/wheat_comparison.rda).
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all("."))
data(wheat); data(wheat_M)
rel  <- list(wheat_M = wheat_M)
ctrl <- bbglr_control(nIter = 5000, burnIn = 2000, thin = 2, nChains = 2, verbose = FALSE)
TERM <- "mrk(gen, wheat_M)"                      # genomic grouping (intercept) label

## Environmental index x -------------------------------------------------------
## The BGLR wheat phenotypes are standardised WITHIN each environment (column
## means 0, sds 1), so the classic Finlay-Wilkinson index (the environment mean)
## is flat and carries no gradient. There is, however, strong CROSSOVER GxE
## (env 1 is negatively correlated with envs 2/4/5). We therefore define the
## environmental index as the leading GxE axis: the first right-singular-vector
## (environment score) of the genotype x environment matrix -- the AMMI/FA
## environmental score. It orders environments along the dominant interaction
## pattern and, unlike the mean, has real spread.
Ywide <- reshape(wheat, idvar = "gen", timevar = "env", direction = "wide")
Ymat  <- as.matrix(Ywide[, paste0("yield.", levels(wheat$env))])
rownames(Ymat) <- as.character(Ywide$gen)
sv    <- svd(scale(Ymat, center = TRUE, scale = FALSE))
escore <- sv$v[, 1]
if (sum(escore) < 0) escore <- -escore                 # sign convention
xmap  <- data.frame(env = levels(wheat$env),
                    x = as.numeric(scale(escore)))      # standardised env score
dat   <- merge(wheat, xmap, by = "env")
cat("env index x (leading GxE axis):\n"); print(xmap)
cat("var explained by axis 1:", round(sv$d[1]^2 / sum(sv$d^2), 3), "\n")

## Nested ladder: genomic main effect (GBLUP) + genomic leg() slopes of order 0..4
fit_order <- function(k) {
  slope <- if (k == 0) NULL else sprintf(" + mrk(gen, wheat_M):leg(x, %d)", k)
  rf <- stats::as.formula(paste0("~ mrk(gen, wheat_M)", slope %||% "", " + env"))
  bbglr(yield ~ leg(x), random = rf, residual = ~ units, data = dat,
        relmat = rel, control = ctrl)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
fits <- list()
for (k in 0:4) { cat(sprintf("\n==== fitting order %d ====\n", k)); fits[[k + 1]] <- fit_order(k) }
names(fits) <- paste0("M", 0:4)
labs <- c(M0 = "GBLUP (order 0)", M1 = "RR linear (1)", M2 = "RR quadratic (2)",
          M3 = "RR cubic (3)", M4 = "RR quartic (4)")

## ---- model-fit / DIC ladder -----------------------------------------------
mf  <- lapply(fits, model_fit)
dic <- data.frame(
  model = names(fits), label = labs[names(fits)],
  dic  = vapply(mf, `[[`, numeric(1), "dic"),
  pD   = vapply(mf, `[[`, numeric(1), "pD"),
  r2   = vapply(mf, `[[`, numeric(1), "r2"),
  rmse = vapply(mf, `[[`, numeric(1), "rmse"),
  varE = vapply(mf, `[[`, numeric(1), "varE"),
  row.names = NULL)
cat("\n-- DIC ladder --\n"); print(dic)

## ---- variance components + heritability -----------------------------------
vc <- lapply(fits, varcomp)
h2_main <- lapply(fits, function(f) heritability(f, genetic = TERM)$summary)
h2_coef <- lapply(fits[-1], function(f) {
  tt <- attr(f$meta[[grep(":leg", vapply(f$meta, `[[`, "", "label"), fixed = TRUE)[1]]], "label")
  heritability(f, genetic = grep(":leg", vapply(f$meta, `[[`, "", "label"),
                                 value = TRUE)[1])$summary
})

## ---- genotype merit: GEBV intercepts per model + rank correlations ---------
gebv <- data.frame(gen = solution(fits$M0, term = TERM, type = "random")$effect)
for (m in names(fits)) {
  s <- solution(fits[[m]], term = TERM, type = "random")
  gebv[[m]] <- s$solution[match(gebv$gen, s$effect)]
}
rank_cor <- cor(gebv[, names(fits)], method = "spearman")
cat("\n-- Spearman rank corr of GEBV intercepts across models --\n"); print(round(rank_cor, 3))

## ---- reaction-norm curves (quadratic model, all genotypes) -----------------
best <- dic$model[which.min(dic$dic)]
rr_term <- grep(":leg", vapply(fits[["M2"]]$meta, `[[`, "", "label"), value = TRUE)[1]
curves  <- reaction_norm(fits[["M2"]], term = rr_term, add_fixed_reg = TRUE,
                         leg_basis = FALSE, n_grid = 120, plot = FALSE)

## ---- across-gradient analytics for each RR model ---------------------------
grad <- lapply(fits[-1], function(f) {
  tm <- grep(":leg", vapply(f$meta, `[[`, "", "label"), value = TRUE)[1]
  rr_gradient(f, term = tm, n_grid = 12L, threshold = 0.20, leg_basis = FALSE, plot = FALSE)
})
names(grad) <- names(fits)[-1]

## ---- GxE summary (adaptability / responsiveness / stability) from M2 -------
gxe_tab <- gxe(fits[["M2"]], term = rr_term)

## ---- per-environment genotype ranking (M2 curves at each env x) ------------
env_x   <- xmap[order(xmap$x), ]
rn_env  <- reaction_norm(fits[["M2"]], term = rr_term, add_fixed_reg = TRUE,
                         leg_basis = FALSE, n_grid = 200, plot = FALSE)
gg <- attr(rn_env, "grid")
pred_env <- sapply(env_x$x, function(xv) {
  gi <- which.min(abs(gg - xv))
  v  <- rn_env[abs(rn_env$gradient - gg[gi]) < 1e-9, ]
  stats::setNames(v$value, v$id)
})
colnames(pred_env) <- env_x$env
pred_env <- pred_env[gebv$gen, , drop = FALSE]

## ---- selection PROBABILITIES for the best model ---------------------------
## Fully-Bayesian: pr() ranks EVERY MCMC draw independently and reports the
## fraction of draws in which a genotype lands in the top fraction -> the
## posterior probability of being a top selection (honest ranking uncertainty).
best_fit <- fits[[best]]
pr_sel <- data.frame(
  gen   = gebv$gen,
  gebv  = gebv[[best]])
for (thr in c(0.05, 0.10, 0.20)) {
  p <- pr(best_fit, term = TERM, type = "random", threshold = thr, higher = TRUE)
  pr_sel[[sprintf("p_top%02d", round(100 * thr))]] <-
    p$prob[match(pr_sel$gen, p[[1]])]
}
cat("\n-- selection probabilities (head, ordered by P(top10%)) --\n")
print(utils::head(pr_sel[order(-pr_sel$p_top10), ], 10))

## pairwise P(A beats B) among the top genotypes by overall merit -------------
D <- attr(solution(best_fit, term = TERM, type = "random"), "draws")  # [nDraws x nGeno]
top_ids <- as.character(gebv$gen[head(order(-gebv[[best]]), 12)])
stopifnot(all(top_ids %in% colnames(D)))
Dt <- D[, top_ids, drop = FALSE]
pr_pair <- matrix(0, length(top_ids), length(top_ids), dimnames = list(top_ids, top_ids))
for (i in seq_along(top_ids)) pr_pair[i, ] <- colMeans(Dt[, i] > Dt)  # P(row > col)
diag(pr_pair) <- NA

res <- list(dic = dic, vc = vc, h2_main = h2_main, h2_coef = h2_coef,
            gebv = gebv, rank_cor = rank_cor, curves = as.data.frame(curves),
            grad = grad, gxe = as.data.frame(gxe_tab), gxe_attr = attributes(gxe_tab)[c("term","higher")],
            pred_env = pred_env, xmap = xmap, env_x = env_x, labs = labs,
            pr_sel = pr_sel, pr_pair = pr_pair, top_ids = top_ids,
            best = best, rr_term = rr_term, ctrl = ctrl,
            n_gen = nlevels(wheat$gen), n_env = nlevels(wheat$env),
            n_mrk = ncol(wheat_M), n_obs = nrow(wheat))
saveRDS(res, "report/wheat_comparison.rds")
cat("\nSAVED report/wheat_comparison.rds  best model:", best, "\n")
