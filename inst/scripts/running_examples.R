#Running examples

data <- read.csv("R/data_examples/data_soy.csv")
data |> 
  dplyr::select(env,gen) |> 
  table()

length(unique(data$env))
length(unique(data$gen))

str(data)
check <- sample(unique(data$gen), 5)
data$check <- data$gen%in%check
colnames(data)

library(SpATS)

env <- unique(data$env)
test <- data[data$env == "E01",]
test$block <- as.factor(test$block)
test$gen <- as.factor(test$gen)
env

blue_data <- data.frame()
for(i in env){
  
  temp <- data |> 
    dplyr::filter(env == i) |> 
    dplyr::mutate(gen = factor(gen),
                  block = factor(block))
  
  m0 <- SpATS(response = "yield",
              spatial = ~ SAP(col, row, nseg = c(10,20), degree = 3, pord = 2), 
              genotype.as.random = FALSE,
              genotype = "gen",
              random = ~ block,
              data = temp, 
              control =  list(tolerance = 1e-03))
  
  blue_temp <- predict(m0, which = "gen") |> 
    dplyr::select(gen,predicted.values) |> 
    dplyr::mutate(env = i) |> 
    dplyr::rename(yield = predicted.values)
  
  blue_data <- rbind(blue_data,blue_temp)

  
}

blue_data$check <- blue_data$gen%in%check

W <- as.matrix(readRDS("R/data_examples/W_soy.rds"))


model <- FW_bglr(data = blue_data, 
                  gen = "gen", 
                  env = "env", 
                  trial = NULL, 
                  block = NULL, 
                  check = "check", 
                  trait = "yield",
                  res.het = FALSE, 
                  nIter = 10000, 
                  burnIn = 5000, 
                  verbose = FALSE,
                  saveAt = "~/learning-claudio/R/BGLR_outs/", 
                  kinship.matrix = NULL,
                  W = W,
                  order = 3,
                  seed = 123)





out <- extr_outs_bglr(model,
                      gen = "gen")


out$effective_sizes
results <- prob_sup_bglr(out)
results$across$bayesian_stats
model$fit2$fit$DIC
gen_by_env <- results$within$hpd_gge |> 
  dplyr::select(ENV,ID,Median) |> 
  dplyr::left_join(blue_data,dplyr::join_by(ID == gen, ENV == env)) |> 
  tidyr::drop_na() |> 
  dplyr::group_by(ENV) |> 
  dplyr::reframe(r = cor(yield,Median))
  

mean(gen_by_env$r)


gen_by_env <- results$within$hpd_gge |> 
  dplyr::select(ENV,ID,Median) |> 
  dplyr::left_join(blue_data,dplyr::join_by(ID == gen, ENV == env)) |> 
  tidyr::drop_na()


ggplot(gen_by_env, aes(x = Median, y = yield, color = ENV))+
  geom_point()


out$variances
plot_fw(results)
plot_pair_prob(results)
plot_prob_within(results)

model$INT
model$COEF

Wp <- readRDS("R/data_examples/Wp_soy.rds")
l_new <- predict_env_gradient_pca(model, Wp)
max(l_new,na.rm = TRUE)
min(l_new,na.rm = TRUE)

winners <- prob_best_by_gradient(model, l_new, tie.method = "share")

prob <- reshape2::melt(winners$prob) |> 
  dplyr::rename(genotype = Var1,
                env = Var2,
                prob = value)

prod <- prod_summary_by_gradient(model,l_new)$mean |> 
  reshape2::melt() |> 
  dplyr::rename(genotype = Var1,
                env = Var2,
                yield = value)

# Predict l for new environments (Wp rows named by env)
pred_grad <- data.frame(
  env = rownames(Wp),
  gradient = as.numeric(l_new)  # l_hat for each env
)

# Bring LAT/LON
covamb_pred <- readRDS("R/data_examples/covam_pred.rds")
pred_grad <- pred_grad |>
  dplyr::left_join(covamb_pred, by = "env") |>
  dplyr::select(env, LAT, LON, gradient) |>
  dplyr::distinct() |>
  tidyr::drop_na()

# Scale l to x in [-1,1] using the model’s saved range
rng_env <- model$range
if (diff(rng_env) > 0) {
  x_pred <- 2 * (pred_grad$gradient - rng_env[1]) / diff(rng_env) - 1
} else {
  x_pred <- rep(0, nrow(pred_grad))
}

# Legendre basis (orthonormal), force deg0 = 1
P_pred <- legendre_basis(x_pred, order = model$order, orthonormal = TRUE)
P_pred[, 1] <- 1
colnames(P_pred) <- paste0("deg", 0:model$order)
rownames(P_pred) <- pred_grad$env


pred_grad$leg <- x_pred


# Posterior mean coefficients per genotype 



geno_pred <- prob |>
  dplyr::left_join(pred_grad, by = "env") |>
  dplyr::select(genotype,env, LAT, LON, prob,leg) |> 
  dplyr::left_join(prod, by = c("env","genotype")) |> 
  tidyr::drop_na()



library(sf)
library(terra)
library(gstat)
library(dplyr)
library(ggplot2)

# 1) Ler shapefile e projetar para metros
br_est  <- raster::shapefile("BR_shape_files/BR_UF_2024.shp")
est_sf  <- sf::st_as_sf(br_est)
MS_geo  <- est_sf[est_sf$SIGLA_UF == "MS", ] |> st_make_valid()

crs_target <- 5880  # SIRGAS 2000 / Brazil Albers
MS_proj    <- st_transform(MS_geo, crs_target)

# 2) Construir grade de 5 km cobrindo MS e pegar centróides para interpolação
grid_5km <- st_make_grid(MS_proj, cellsize = 10000, square = FALSE)
grid_5km <- st_intersection(st_sf(id = seq_along(grid_5km), geometry = grid_5km), MS_proj) # recorta
grid_pts <- st_centroid(grid_5km)  # pontos alvo (centros de células)

# 3) Pontos de yield: sf projetado e filtrado dentro de MS
geno_sf <- st_as_sf(geno_pred, coords = c("LON", "LAT"), crs = 4326) |>
  st_transform(crs_target)
geno_sf <- geno_sf[MS_proj, ]

# 4) Interpolar por genótipo com gstat::idw (newdata = centróides da grade)
geno_names <-results$across$perfo |> 
  dplyr::pull(ID)
geno_names <- geno_names[1:30]
geno_names <- c("C054","G087","G088","G089","G090","G100","G114","G120","G170","G177","G187")

# for (g in geno_names) {
#   pts_g <- geno_sf[geno_sf$genotype == g, ]  # preserva classe sf
#   col_name <- paste0("y_", g)
#   
#   if (nrow(pts_g) >= 3) {
#     # 1) Clamp e logit para evitar infinitos
#     eps <- 1e-6
#     p_raw  <- pts_g$prob
#     p_clip <- pmin(pmax(p_raw, eps), 1 - eps)
#     pts_g$logit_p <- log(p_clip / (1 - p_clip))  # logit(p)
#     
#     # 2) Converter para sp (gstat usa 'sp')
#     pts_sp  <- as(pts_g, "Spatial")
#     grid_sp <- as(grid_pts, "Spatial")
#     
#     # 3) Parâmetros iniciais e variograma empírico mais estável
#     #    (use CRS projetado para distâncias em metros)
#     coords <- sf::st_coordinates(pts_g)
#     Dmat   <- sp::spDists(coords, longlat = FALSE)
#     maxd   <- max(Dmat, na.rm = TRUE)
#     width  <- maxd / 12          # largura dos bins
#     cutoff <- 0.8 * maxd         # distância máxima considerada
#     
#     vg <- gstat::variogram(logit_p ~ 1, pts_sp,
#                            cutoff = cutoff, width = width,
#                            cressie = TRUE)  # estimador robusto
#     
#     # Valores iniciais razoáveis a partir da variância da logit_p
#     psill0  <- stats::var(pts_g$logit_p, na.rm = TRUE)
#     nugget0 <- 0.10 * psill0
#     range0  <- maxd / 3
#     
#     # 4) Tentar ajustar Spherical; se falhar, tentar Exponential; se falhar, Gaussian
#     fit_ok <- FALSE
#     vg_models <- list(
#       gstat::vgm(psill = psill0, model = "Sph", range = range0, nugget = nugget0),
#       gstat::vgm(psill = psill0, model = "Exp", range = range0, nugget = nugget0),
#       gstat::vgm(psill = psill0, model = "Gau", range = range0, nugget = nugget0)
#     )
#     
#     vg_fit <- NULL
#     for (vm in vg_models) {
#       vg_fit_try <- try(suppressWarnings(gstat::fit.variogram(vg, vm)), silent = TRUE)
#       bad_fit <- inherits(vg_fit_try, "try-error") ||
#         any(!is.finite(vg_fit_try$psill)) ||
#         any(vg_fit_try$psill < 0) ||
#         !is.finite(sum(vg_fit_try$psill))
#       if (!bad_fit) {
#         vg_fit <- vg_fit_try
#         fit_ok <- TRUE
#         break
#       }
#     }
#     
#     if (!fit_ok) {
#       # Fallback: IDW direto na probabilidade (mantém seu comportamento atual)
#       pred <- gstat::idw(
#         formula   = prob ~ 1,
#         locations = pts_sp,
#         newdata   = grid_sp,
#         idp       = 2
#       )
#       grid_5km[[col_name]] <- pred$var1.pred
#       
#     } else {
#       # 5) Krigagem em logit_p e inversa da logit para probabilidade
#       kr <- gstat::krige(logit_p ~ 1, pts_sp, grid_sp, model = vg_fit)
#       pred_prob <- 1 / (1 + exp(-kr$var1.pred))  # inverse-logit
#       
#       grid_5km[[col_name]] <- pred_prob
#     }
#     
#   } else {
#     grid_5km[[col_name]] <- NA_real_
#   }
# }

for (g in geno_names) {
  pts_g <- geno_sf[geno_sf$genotype == g, ]  # preserva classe sf
  col_name <- paste0("y_", g)                # y_* = produtividade interpolada
  
  if (nrow(pts_g) >= 3) {
    # Converter para sp
    pts_sp  <- as(pts_g, "Spatial")
    grid_sp <- as(grid_pts, "Spatial")
    
    # Distâncias (CRS projetado em metros)
    coords <- sf::st_coordinates(pts_g)
    Dmat   <- sp::spDists(coords, longlat = FALSE)
    maxd   <- max(Dmat, na.rm = TRUE)
    width  <- maxd / 12
    cutoff <- 0.8 * maxd
    
    # Variograma robusto da produtividade
    vg_y <- gstat::variogram(yield ~ 1, pts_sp,
                             cutoff = cutoff, width = width,
                             cressie = TRUE)
    psill0_y  <- stats::var(pts_g$yield, na.rm = TRUE)
    nugget0_y <- 0.10 * psill0_y
    range0_y  <- maxd / 3
    
    # Tenta Sph/Exp/Gau; se falhar, cai para IDW
    fit_ok_y <- FALSE
    vg_models_y <- list(
      gstat::vgm(psill = psill0_y, model = "Sph", range = range0_y, nugget = nugget0_y),
      gstat::vgm(psill = psill0_y, model = "Exp", range = range0_y, nugget = nugget0_y),
      gstat::vgm(psill = psill0_y, model = "Gau", range = range0_y, nugget = nugget0_y)
    )
    vg_fit_y <- NULL
    for (vm in vg_models_y) {
      vg_fit_try <- try(suppressWarnings(gstat::fit.variogram(vg_y, vm)), silent = TRUE)
      bad_fit <- inherits(vg_fit_try, "try-error") ||
        any(!is.finite(vg_fit_try$psill)) ||
        any(vg_fit_try$psill < 0) ||
        !is.finite(sum(vg_fit_try$psill))
      if (!bad_fit) {
        vg_fit_y <- vg_fit_try
        fit_ok_y <- TRUE
        break
      }
    }
    
    if (!fit_ok_y) {
      # Fallback: IDW direto na produtividade
      pred_y <- gstat::idw(
        formula   = yield ~ 1,
        locations = pts_sp,
        newdata   = grid_sp,
        idp       = 2
      )
      grid_5km[[col_name]] <- pred_y$var1.pred
    } else {
      # Krigagem da produtividade
      kr_y <- gstat::krige(yield ~ 1, pts_sp, grid_sp, model = vg_fit_y)
      grid_5km[[col_name]] <- kr_y$var1.pred
    }
    
  } else {
    grid_5km[[col_name]] <- NA_real_
  }
}


# 5) Escolher vencedor por célula (maior yield entre colunas y_* )
y_cols    <- paste0("y_", geno_names)
y_mat     <- as.matrix(st_drop_geometry(grid_5km[, y_cols]))
idx_max   <- max.col(y_mat, ties.method = "first")
has_data  <- rowSums(is.finite(y_mat)) > 0


geno_cols <- grep("^y_", names(grid_5km), value = TRUE)
mat <- as.matrix(st_drop_geometry(grid_5km)[, geno_cols])
has_data <- rowSums(!is.na(mat)) > 0

mat2 <- mat
mat2[is.na(mat2)] <- -Inf  # assim o max.col ignora NAs
idx_max <- max.col(mat2, ties.method = "first")
winner_raw <- ifelse(has_data, geno_cols[idx_max], NA_character_)
grid_5km$winner <- factor(winner_raw, levels = geno_cols)
grid_5km$winner <- gsub("y_","",grid_5km$winner)


ggplot() +
  geom_sf(data = grid_5km, aes(fill = winner), color = "white") +
  geom_sf(data = MS_proj, fill = NA, color = "black", linewidth = 1.2) +
  coord_sf(crs = st_crs(MS_proj), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()



ggplot() +
  geom_sf(data = grid_5km, aes(fill = y_G074), color = "white") +
  geom_sf(data = MS_proj, fill = NA, color = "black", linewidth = 1.2) +
  coord_sf(crs = st_crs(MS_proj), expand = FALSE) +
  viridis::scale_fill_viridis()+
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()





