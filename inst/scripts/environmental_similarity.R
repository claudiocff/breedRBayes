covamb = readRDS("R/data_examples/covamb.rds")
covamb_pred = readRDS("R/data_examples/covam_pred.rds")


matcovamb = subset(covamb, select = c(-DOY, -YYYYMMDD, 
                                      -daysFromStart, -n, -N, -RTA))
matcovamb = matcovamb |> dplyr::reframe(across(LON:soc, mean), .by = 'env')
matcovamb = tibble::column_to_rownames(matcovamb, var = 'env')

# Unobserved environments
covamb_pred = covamb_pred |> 
  dplyr::filter(YEAR != 2000 | !MONTH %in% c('JAN','FEB','MAR')) |> 
  dplyr::filter(YEAR != 2021 | !MONTH %in% c('OCT','NOV','DEC')) |> 
  dplyr::select(-YEAR, -MONTH, - ALT) |> 
  dplyr::reframe(across(ALLSKY_SFC_LW_DWN:soc, mean), .by = c(mun, env)) |>  
  dplyr::rename(T2M_RANGE = T2RANGE, PRECTOT = PRECTOTCORR)

matcovamb = matcovamb[, sort(colnames(matcovamb))]

matcovamb_pred = covamb_pred |> 
  tibble::column_to_rownames('env') |> 
  dplyr::select(-mun)
matcovamb_pred = matcovamb_pred[, sort(colnames(matcovamb_pred))]
matcovamb_pred = as.matrix(matcovamb_pred)


colnames(matcovamb)
colnames(matcovamb_pred)


matcovamb <- matcovamb[,colnames(matcovamb)%in%colnames(matcovamb_pred)]
matcovamb_pred <- matcovamb_pred[,colnames(matcovamb_pred)%in%colnames(matcovamb)]


dim(matcovamb)
dim(matcovamb_pred)


all(colnames(matcovamb) == colnames(matcovamb_pred))

eucl = as.matrix(pdist::pdist(scale(matcovamb), scale(matcovamb_pred)))

dim(eucl)

rownames(eucl) = rownames(matcovamb)
colnames(eucl) = rownames(matcovamb_pred)


dist.env = lapply(apply(eucl, 1, function(x) reshape2::melt(x)), 
                  function(x) x |> 
                    dplyr::rename(dist = value) |> 
                    tibble::rownames_to_column('point') |> 
                    dplyr::left_join(covamb_pred[,c('env', 'LON', 'LAT')], 
                                     by = c('point' = 'env')))
rm(eucl)

br_est = raster::shapefile("BR_shape_files/BR_UF_2024.shp")
est_sf = sf::st_as_sf(br_est)
MS = est_sf[est_sf$SIGLA_UF == "MS",]
MS_pol = br_est[br_est$SIGLA_UF == 'MS',]


obs_window = spatstat.geom::owin(xrange = MS_pol@bbox[1,], yrange = MS_pol@bbox[2,])

cl = parallel::makeCluster(parallel::detectCores(logical = F)-1)
parallel::clusterEvalQ(cl, {library(spatstat)
  library(raster) 
  library(sf)})
parallel::clusterExport(cl, varlist = c('obs_window', 'MS'))


dist.list = parallel::parLapply(cl = cl, dist.env, function(x){
  
  ppp_test = spatstat.geom::ppp(x = x$LON, y = x$LAT, 
                                marks = x$dist, window = obs_window)
  
  powers <- seq(0.1, 5, 0.1)
  mse_result <- NULL
  for(power in powers){
    CV_idw <- spatstat.explore::idw(ppp_test, power=power, at="points")
    mse_result <- c(mse_result, mean((ppp_test$marks-CV_idw)^2, na.rm=T))
  }
  optimal_power <- powers[which.min(mse_result)]
  
  tr_idw = spatstat.explore::idw(ppp_test, power = optimal_power, at = 'pixels')
  
  idw_raster = as(raster::mask(raster::raster(tr_idw, crs = crs(br_est)), MS), 
                  "SpatialPixelsDataFrame")
  
  idw_raster.df = as.data.frame(idw_raster)
})
parallel::stopCluster(cl)
rm(cl)

save(dist.list, file = 'R/data_examples/dist_list.RData')


mindist = apply(do.call(cbind, lapply(dist.list, function(x) x[,-c(2,3)])), 
                1, function(x) x[which.min(x)])



cbind(dist.list$E01[,-1], 
      layer = abs((mindist - min(mindist))/(max(mindist) - min(mindist))-1)) |> 
  ggplot() +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true",
                                    pad_x = unit(0.5,"cm"), pad_y = unit(0.7,"cm"),
                                    style = ggspatial::north_arrow_nautical())+ 
  geom_raster(aes(x = x, y = y , fill = layer), interpolate = T) +
  geom_sf(data = MS, fill = NA, linewidth = 1, color = "black")  + 
  labs(x = 'Longitude', y = 'Latitude', fill = 'Environmental \n similarity') + 
  scale_fill_viridis_c(option = 'turbo', direction = 1, 
                       limits = c(0, 1)) + 
  theme(legend.position = 'right') + 
  geom_point(
    data = data.frame("lat" = c(-22.187, -22.191,-21.434,-22.634,
                                -20.442,-22.078,-22.304,-21.614,
                                -23.065,-21.801,-19.395,-20.931,-20.443),
                      "lon" = c(-52.717, -55.947,-54.438,-54.822,-54.646,
                                -54.789,-53.815,-55.168,-54.190,-54.546,
                                -54.566,-54.961,-54.868)),
    aes(x = lon, y = lat, size = 1.2), show.legend = F, color = 'black',
  ) 
rm(dist.list, dist.env)































