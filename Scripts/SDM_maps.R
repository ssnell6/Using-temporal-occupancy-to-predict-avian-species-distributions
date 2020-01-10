library(tidyverse)
library(glm2)
library(gam)
library(randomForest)
library(dismo)
library(raster)
library(maptools)
library(rgdal)
library(mice)
library(pROC)
library(hydroGOF)
library(tmap)
library(RColorBrewer)

bbs_occ = read.csv("Data/bbs_2001_2015.csv", header=TRUE)
bbs_occ_sub = bbs_occ %>% 
  group_by(aou) %>%
  dplyr::count(stateroute) %>% 
  dplyr::mutate(occ = n/15) 

exp_pres = read.csv("Data/expect_pres.csv", header = TRUE) %>%
  filter(!stateroute %in% bad_rtes$stateroute)
# remove routes where bbs_occ_sub
traits = read.csv("Data/Master_RO_Correlates.csv", header = TRUE)
bsize = read.csv("data/DunningBodySize_old_2008.11.12.csv", header = TRUE)
lat_long = read.csv("Data/latlongs.csv", header = TRUE)
tax_code = read.csv("Data/Tax_AOU_Alpha.csv", header = TRUE)
bad_rtes = read.csv("Data/bad_rtes.csv", heade = TRUE)
bi_env = read.csv("Data/all_env.csv", header = TRUE)
bi_means = bi_env[,c("stateroute","mat.mean", "elev.mean", "map.mean", "ndvi.mean")]
env_bio = read.csv("Data/env_bio.csv", header = TRUE)
env_bio = na.omit(env_bio)
env_bio_sub = env_bio[,c(1, 21:39)]


##### read in raw bbs data for 2016 ####
bbs_new <- read.csv("Data/bbs_2016.csv", header = TRUE) 
bbs_new$presence = 1
bbs_new_exp_pres <- read.csv("Data/expect_pres_2016.csv", header = TRUE)
bbs_new_all <- left_join(bbs_new_exp_pres, bbs_new, by = c("spAOU"="aou", "stateroute" = "stateroute"))
bbs_new_all$presence <- case_when(is.na(bbs_new_all$presence) == TRUE ~ 0, 
                                  bbs_new_all$presence == 1 ~ 1)

all_env = left_join(bi_means, env_bio_sub, by = "stateroute")

#update tax_code Winter Wren
tax_code$AOU_OUT[tax_code$AOU_OUT == 7220] <- 7222
tax_code$AOU_OUT[tax_code$AOU_OUT == 4810] <- 4812
tax_code$AOU_OUT[tax_code$AOU_OUT == 4123] <- 4120

# BBS cleaning
bbs_inc_absence = full_join(bbs_occ_sub, exp_pres, by = c("aou" ="spAOU", "stateroute" = "stateroute")) %>%
  dplyr::select(aou, stateroute, occ)
bbs_inc_absence$occ[is.na(bbs_inc_absence$occ)] <- 0
bbs_inc_absence$presence = 0
bbs_inc_absence$presence[bbs_inc_absence$occ > 0] <- 1
num_occ = bbs_inc_absence %>% group_by(aou) %>% tally(presence) %>% left_join(bbs_inc_absence, ., by = "aou")

# 412 focal species
bbs_final_occ = filter(num_occ,n > 49)
bbs_occ_code = left_join(bbs_final_occ, tax_code, by = c("aou" = "AOU_OUT"))

# 319 focal species
bbs_final_occ_ll = left_join(bbs_occ_code, lat_long, by = "stateroute") %>%
  filter(aou %in% tax_code$AOU_OUT & stateroute %in% bbs_occ_sub$stateroute) 

bbs_final_occ_ll$sp_success = 15 * bbs_final_occ_ll$occ
bbs_final_occ_ll$sp_fail = 15 * (1 - bbs_final_occ_ll$occ) 
bbs_final_occ_ll$presence <- as.numeric(bbs_final_occ_ll$presence)
# temp filter for vis purposes
auc_df = read.csv("Data/auc_df.csv", header = TRUE)

#### change spp here ##### 
sdm_input <- filter(bbs_final_occ_ll, aou == 6280) %>% 
  left_join(all_env, by = "stateroute") 
sdm_notrans <- filter(sdm_input, occ > 0.33| occ == 0) %>% na.omit(.)

# Determine geographic extent of our data using AOU = i
max.lat <- ceiling(max(sdm_input$latitude))
min.lat <- floor(min(sdm_input$latitude))
max.lon <- ceiling(max(sdm_input$longitude))
min.lon <- floor(min(sdm_input$longitude))

glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_input)
gam_pres <- mgcv::gam(presence ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = sdm_input)
rf_occ <- randomForest(sp_success/15 ~elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
rf_pres <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)

max_ind_pres = maxent(sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")], sdm_input$presence)

# predict
pred_glm_occ <- predict(glm_occ,type=c("response"))
pred_glm_pr <- predict(glm_pres,type=c("response"))
pred_gam_occ <- predict(gam_occ,type=c("response"))
pred_gam_pr <- predict(gam_pres,type=c("response"))
pred_rf_occ <- predict(rf_occ,type=c("response"))
pred_rf_pr <- predict(rf_pres,type=c("response"))
max_pred_pres <- predict(max_ind_pres, sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")])

sdm_output = data.frame(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_pres)

##### predict #####
# test.poly <- readShapePoly("Z:/GIS/birds/All/All/Vireo_flavifrons_22705237.shp")
# proj4string(test.poly) <- intl_proj
# 
# all_env_ll <- full_join(all_env, sdm_input[,c("stateroute", "aou")], by = "stateroute") %>%
#   left_join(lat_long, by = "stateroute")
# coordinates(all_env_ll) <- c("longitude", "latitude")
# proj4string(all_env_ll) <- intl_proj
# sporigin = test.poly[test.poly@data$SEASONAL == 1|test.poly@data$SEASONAL == 2|test.poly@data$SEASONAL ==5,]
# routes_inside <- all_env_ll[!is.na(sp::over(all_env_ll, as(sporigin,"SpatialPolygons"))),]
# routes_inside = data.frame(routes_inside) 
#   
# routes_inside$pred_glm_occ <- predict(glm_occ, newdata = routes_inside)
# routes_inside$pred_glm_pr <- predict(glm_pres,newdata = routes_inside)
# routes_inside$pred_gam_occ <- predict(gam_occ,newdata = routes_inside)
# routes_inside$pred_gam_pr <- predict(gam_pres,newdata = routes_inside)
# routes_inside$pred_rf_occ <- predict(rf_occ,newdata = routes_inside)
# routes_inside$pred_rf_pr <- predict(rf_pres,newdata = routes_inside)
# max_pred_pres <- predict(max_ind_pres, routes_inside[,c("elev.mean", "ndvi.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14")])
# routes_inside_pred <- data.frame(routes_inside, max_pred_pres)
# 
# sdm_output = full_join(sdm_input[,c("aou","stateroute","occ","presence","n","PRIMARY_COM_NAME","ALPHA.CODE", "sp_success", "sp_fail")], routes_inside_pred, by = "stateroute") %>%
#   na.omit()

##### multiple imputation #####
max.lat <- ceiling(max(sdm_input$latitude))
min.lat <- floor(min(sdm_input$latitude))
max.lon <- ceiling(max(sdm_input$longitude))
min.lon <- floor(min(sdm_input$longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

mod.r <- SpatialPointsDataFrame(coords = sdm_output[,c("longitude", "latitude")],
                                data = sdm_output[,c("latitude", "longitude", "pred_glm_pr", "pred_glm_occ", "pred_gam_pr", "pred_gam_occ",  "pred_rf_pr", "pred_rf_occ","max_pred_pres")], 
                                proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
r = raster(nrows =22,ncols = 30, geographic.extent, 1) 
projection(r) <- "+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"
p = rasterToPoints(r)
p = data.frame(p)
names(p) = c("longitude", "latitude")
mod.p <- SpatialPointsDataFrame(coords = p, data = p, proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
sdm_output_impute = full_join(sdm_output, p, by = c("latitude", "longitude")) %>%
  select(latitude, longitude, pred_glm_occ, pred_glm_pr, pred_gam_occ, pred_gam_pr, pred_rf_occ, pred_rf_pr,  max_pred_pres)
imputed_Data <- mice(sdm_output_impute, m=5, maxit = 50, method = 'pmm', seed = 500)
completedData <- complete(imputed_Data,1)

# original code
mod.r <- SpatialPointsDataFrame(coords = completedData[,c("longitude", "latitude")],
                                data = completedData, proj4string = CRS("+proj=longlat +datum=WGS84"))
r.r = raster(mod.r, res = 1) # 40x40 km/111 (degrees) * 2 tp eliminate holes, bioclim is 4 km
plot.r <- rasterize(mod.r, r.r)
mask <- readOGR("Z:/GIS/geography/continent.shp") 
projection(mask) <- "+proj=longlat +datum=WGS84"
r_mask <- mask(plot.r, mask)


##### no trans ######
glm_occ_notrans <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_notrans)
glm_pres_notrans <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_notrans)
gam_occ_notrans <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_notrans)
gam_pres_notrans <- mgcv::gam(presence ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = sdm_notrans)
rf_occ_notrans <- randomForest(sp_success/15 ~elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_notrans)
rf_pres_notrans <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_notrans)
max_ind_pres_notrans = maxent(sdm_notrans[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")], sdm_notrans$presence)
# predict
pred_glm_occ_notrans <- predict(glm_occ_notrans,type=c("response"))
pred_glm_pr_notrans <- predict(glm_pres_notrans,type=c("response"))
pred_gam_occ_notrans <- predict(gam_occ_notrans,type=c("response"))
pred_gam_pr_notrans <- predict(gam_pres_notrans,type=c("response"))
pred_rf_occ_notrans <- predict(rf_occ_notrans,type=c("response"))
pred_rf_pr_notrans <- predict(rf_pres_notrans,type=c("response"))
max_pred_pres_notrans <- predict(max_ind_pres_notrans, sdm_notrans[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")])
sdm_output_notrans = data.frame(sdm_notrans, pred_glm_pr_notrans, pred_glm_occ_notrans, pred_gam_pr_notrans, pred_gam_occ_notrans, pred_rf_occ_notrans, pred_rf_pr_notrans, max_pred_pres_notrans) 

#### rmse ####
rmse_occ <- rmse(sdm_output$pred_glm_occ, sdm_output$occ)
rmse_pres <- rmse(sdm_output$pred_glm_pr, as.numeric(sdm_output$presence))
rmse_notrans <- rmse(sdm_output_notrans$pred_glm_pr_notrans,as.numeric(sdm_output_notrans$presence))
rmse_gam <- rmse(as.vector(sdm_output$pred_gam_occ), sdm_output$occ)
rmse_gam_pres <- rmse(as.vector(sdm_output$pred_gam_pr), as.numeric(sdm_output$presence))
rmse_gam_notrans <- rmse(as.vector(sdm_output_notrans$pred_gam_pr_notrans), as.numeric(sdm_output_notrans$presence))
rmse_rf <- rmse(sdm_output$pred_rf_occ, sdm_output$occ)
rmse_rf_pres <- rmse(as.vector(as.numeric(sdm_output$pred_rf_pr)), as.numeric(sdm_output$presence))
rmse_rf_notrans <- rmse(as.vector(as.numeric(sdm_output_notrans$pred_rf_pr_notrans)), sdm_output_notrans$presence)
rmse_me_pres <- rmse(sdm_output$max_pred_pres, as.numeric(sdm_output$presence))
rmse_me_pres_notrans <- rmse(sdm_output_notrans$max_pred_pres_notrans, as.numeric(sdm_output_notrans$presence))


#### multiple imputation ####
mod.core <- SpatialPointsDataFrame(coords = sdm_output_notrans[,c("longitude", "latitude")],
                                   data = sdm_output_notrans[,c("latitude", "longitude", "pred_glm_pr_notrans", "pred_glm_occ_notrans", "pred_gam_pr_notrans", "pred_gam_occ_notrans", "pred_rf_occ_notrans", "pred_rf_pr_notrans", "max_pred_pres_notrans")], proj4string = CRS("+proj=longlat +datum=WGS84"))
r = raster(nrows =22,ncols = 30, geographic.extent, 1) 
projection(r) <- "+proj=longlat +datum=WGS84"
p = rasterToPoints(r)
p = data.frame(p)
names(p) = c("longitude", "latitude")
mod.p <- SpatialPointsDataFrame(coords = p, data = p, proj4string = CRS("+proj=longlat +datum=WGS84"))
sdm_output_impute.core = full_join(sdm_output_notrans, p, by = c("latitude", "longitude")) %>%
  select(latitude, longitude, pred_glm_occ_notrans, pred_glm_pr_notrans, pred_gam_occ_notrans, pred_gam_pr_notrans, pred_rf_occ_notrans, pred_rf_pr_notrans,  max_pred_pres_notrans)
imputed_Data.core <- mice(sdm_output_impute.core, m=5, maxit = 50, method = 'pmm', seed = 500)
completedData.core <- complete(imputed_Data.core,1)

# original code
mod.core <- SpatialPointsDataFrame(coords = completedData.core[,c("longitude", "latitude")],
                   data = completedData.core, proj4string = CRS("+proj=longlat +datum=WGS84"))
r.core = raster(mod.core, res = 1) # 40x40 km/111 (degrees) * 2 tp eliminate holes, bioclim is 4 km
plot.core <- rasterize(mod.core, r.core)
c_mask <- mask(plot.core, mask)

sdm_output$presence <- factor(sdm_output$presence,
                              levels = c(1,0), ordered = TRUE)

##### plots ####
us_sf <- read_sf("Z:/GIS/geography/continent.shp")
us_sf <- st_transform(us_sf, crs = "+proj=longlat +datum=WGS84")

sdm_output$core <- 0
sdm_output$core[sdm_output$stateroute %in% sdm_output_notrans$stateroute] <- 1
routes_sf <- st_as_sf(sdm_output,  coords = c("longitude", "latitude"))
# CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km")
routes_notrans <- st_as_sf(sdm_output[sdm_output$occ >= 0.33,], coords = c("longitude", "latitude"))

us <- tm_shape(us_sf) + tm_borders() + tm_fill(col = "white")

#### need to add in core spp ####
routes_sf$presence_cat <- case_when(routes_sf$occ >= 0.33 ~ "core",
                                    routes_sf$occ < 0.33 & routes_sf$occ > 0 ~ "trans",
                                    routes_sf$occ == 0 ~ "absent")
routes_sf$presence_cat <- factor(routes_sf$presence_cat, levels = c("core", "trans", "absent"), ordered = TRUE)
# point_map <- tm_shape(routes_sf, legend.show = FALSE) + 
#  tm_symbols(size = 0.75, shape="presence_cat", shapes = c(16,16,4), col = 'presence_cat', palette = c(
#    "#008837","#7fbf7b","#7b3294")) + 
#  tm_shape(us_sf) + tm_borders( "black", lwd = 3) +
#  tm_layout("Occurrence", title.size = 2, title.position = c("right","bottom")) 
palette = brewer.pal(5, "PRGn")

scale_fun <- function(r){ 
  min <- min(na.omit(r))
  max <- max(na.omit(r))
  break1 <- min+((max-min)/5)
  break2 <- break1+((max-min)/5)
  break3 <- break2+((max-min)/5)
  break4 <- break3+((max-min)/5)
  pal <- c(min, break1, break2, break3, break4, max)
}

point_map <- tm_shape(routes_sf) + 
  tm_symbols(size = 0.75, shape="presence", shapes = c(16,4), alpha = 0.5, col = "black") + 
  tm_legend(show=FALSE) +
  tm_shape(us_sf) + tm_borders( "black", lwd = 3) + 
  tm_shape(routes_notrans)  + 
  tm_symbols(col = "presence", palette = "-PRGn", size = 0.75, shapes = c(16,4)) + tm_legend(outside = TRUE)+ 
  tm_layout("  Observed \nOccurences", title.size = 2, title.position = c("right","bottom")) +
  tm_layout(main.title = "A") 

point_map_allen <- tm_shape(routes_sf) + 
  tm_symbols(size = 0.75, shape="presence", shapes = c(16,4), alpha = 0.5, col = "black") + 
  tm_legend(show=FALSE) +
  tm_shape(us_sf) + tm_borders( "black", lwd = 3) + 
  tm_shape(routes_notrans)  + 
  tm_symbols(col = "occ", palette = "PRGn", size = 0.75, shapes = c(16,4)) + tm_legend(outside = TRUE)+ 
  tm_layout("  Observed \nOccurences", title.size = 2, title.position = c("right","bottom")) +
  tm_layout(main.title = "A") 

sdm_maxent_pr <- tm_shape(r_mask) + tm_raster("max_pred_pres", palette = palette, style = "cont", breaks=quantile(r_mask$max_pred_pres, probs = seq(0.2,0.8, by = 0.2)) , legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout(paste("RMSE =",signif(rmse_me_pres, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "B") 

sdm_maxent_core <- tm_shape(c_mask) + tm_raster("max_pred_pres_notrans", palette = palette, style = "cont", title = "MaxEnt Core", breaks=quantile(c_mask$max_pred_pres_notrans, probs = seq(0.2,0.8, by = 0.2)),legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout(paste("RMSE =",signif(rmse_me_pres_notrans, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "C") 

sdm_glm_occ <- tm_shape(r_mask) + tm_raster("pred_glm_occ", palette = palette, style = "cont", title = "GLM Occ",breaks=quantile(r_mask$pred_glm_occ, probs = seq(0.2,0.8, by = 0.2)), legend.show = FALSE) + 
  tm_shape(us_sf) + tm_borders( "black", lwd = 3) + 
  tm_layout(paste("RMSE =",signif(rmse_occ, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "D") 
#sdm_occ

sdm_glm_pr <- tm_shape(r_mask) + tm_raster("pred_glm_pr", palette = palette, style = "cont", title = "GLM Pres", breaks=quantile(r_mask$pred_glm_pr, probs = seq(0.2,0.8, by = 0.2)), legend.show = FALSE) + tm_shape(us_sf) + 
  tm_borders(col = "black", lwd = 3) + 
  tm_layout(paste("RMSE =",signif(rmse_pres, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "E") 
#sdm_pr

sdm_glm_core <- tm_shape(c_mask) + tm_raster("pred_glm_pr_notrans", palette = palette, style = "cont", breaks=quantile(c_mask$pred_glm_pr_notrans, probs = seq(0.2,0.8, by = 0.2)),legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + 
  tm_layout(paste("RMSE =",signif(rmse_notrans, 2)), title.size = 2, title.position = c("right","bottom")) +
  tm_layout(main.title = "F") 

sdm_gam_occ <- tm_shape(r_mask) + tm_raster("pred_gam_occ", palette = palette, style = "cont", title = "GAM Occ", breaks=quantile(r_mask$pred_gam_occ, probs = seq(0.2,0.8, by = 0.2)), legend.show = FALSE) + tm_shape(us_sf) + tm_borders( "black", lwd = 3) + tm_layout(paste("RMSE =",signif(rmse_gam, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "G") 
#sdm_occ

sdm_gam_pr <- tm_shape(r_mask) + tm_raster("pred_gam_pr", palette = palette, style = "cont",breaks=quantile(r_mask$pred_gam_pr, probs = seq(0.2,0.8, by = 0.2)), legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + 
  tm_layout(paste("RMSE =",signif(rmse_gam_pres, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "H") 
#sdm_pr

sdm_gam_core <- tm_shape(c_mask) + tm_raster("pred_gam_pr_notrans", palette = palette, style = "cont", breaks=quantile(c_mask$pred_gam_pr_notrans, probs = seq(0.2,0.8, by = 0.2)), legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout(paste("RMSE =",signif(rmse_gam_notrans, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "I") 

# had to do manual bc mean was less than first break.
r_mask$pred_rf_occ <- abs(r_mask$pred_rf_occ)
sdm_rf_occ <- tm_shape(r_mask) + tm_raster("pred_rf_occ", palette = palette, style = "cont", title = "RF Occ", breaks=quantile(r_mask$pred_rf_occ, probs = seq(0.2,0.8, by = 0.2)), legend.show = FALSE) + tm_shape(us_sf) + tm_borders( "black", lwd = 3) + tm_layout(paste("RMSE =",signif(rmse_rf, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "J") 
#sdm_occ

sdm_rf_pr <- tm_shape(r_mask) + tm_raster("pred_rf_pr", palette = palette, style = "cont", title = "RF Pres", breaks=quantile(r_mask$pred_rf_pr, probs = seq(0.2,0.8, by = 0.2)), legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout(paste("RMSE =",signif(rmse_rf_pres, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "K") 
#sdm_pr

# had to do manual bc mean was less than first break.
c_mask$pred_rf_pr_notrans <- abs(c_mask$pred_rf_pr_notrans)
sdm_rf_core <- tm_shape(c_mask) + tm_raster("pred_rf_pr_notrans", palette = palette, style = "cont", title = "RF Core", breaks=scale_fun(as.vector(c_mask$pred_rf_pr_notrans)), legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout(paste("RMSE =",signif(rmse_rf_notrans, 2)), title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "L") 

MaxEnt_plot <- tmap_arrange(point_map, sdm_maxent_pr, sdm_maxent_core, ncol = 1)
fig_glm <- tmap_arrange(sdm_glm_occ, sdm_glm_pr, sdm_glm_core, ncol = 1)
fig_gam <- tmap_arrange(sdm_gam_occ, sdm_gam_pr, sdm_gam_core, ncol = 1)
fig_rf <- tmap_arrange(sdm_rf_occ, sdm_rf_pr, sdm_rf_core, ncol = 1)

# final_fig1 <- tmap_arrange(point_map, sdm_maxent_pr, sdm_maxent_core, sdm_glm_occ, sdm_glm_pr, sdm_glm_core, sdm_gam_occ, sdm_gam_pr, sdm_gam_core, sdm_rf_occ, sdm_rf_pr, sdm_rf_core, nrow = 4, ncol = 3) 
# tmap_save(final_fig1, "Figures/Figure1.pdf", height = 26, width = 20)


allen_fig1 <- tmap_arrange(point_map_allen, sdm_maxent_pr, sdm_maxent_core, sdm_glm_occ, sdm_glm_pr, sdm_glm_core, sdm_gam_occ, sdm_gam_pr, sdm_gam_core, sdm_rf_occ, sdm_rf_pr, sdm_rf_core, nrow = 4, ncol = 3) 
tmap_save(allen_fig1, "Figures/Figure1_allen.pdf", height = 16, width = 20)


#### MAPS #####
auc_df = read.csv("Data/auc_df.csv", header = TRUE)
setwd("Figures/maps/")


mapfun <- function(pdf_name, vec, num){ 
  
  pdf(pdf_name, height = 8, width = 10)
  par(mfrow = c(2, 3)) # makes plots too small
  
  for(i in unique(sp_list)){
    print(i)
    sdm_output = c()
    
    bbs_sub <- filter(bbs_final_occ_ll, aou == i)
    temp <- filter(all_env, stateroute %in% bbs_sub$stateroute)
    sdm_input <- left_join(bbs_sub, temp, by = "stateroute")
    sdm_input = na.omit(sdm_input)
    j = unique(sdm_input$ALPHA.CODE)
    
    # print(length(sdm_input$stateroute))
    if(length(unique(sdm_input$stateroute)) > 40  & length(unique(sdm_input$presence)) >1){
      #   if(levels(as.factor(sdm_input$presence)) > 1){
      glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_input)
      gam_pres <- mgcv::gam(presence ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = sdm_input)
      rf_occ <- randomForest(sp_success/15 ~elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      rf_pres <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      
      max_ind_pres = dismo::maxent(sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")], sdm_input$presence)
      
      # predict
      pred_glm_occ <- predict(glm_occ,type=c("response"))
      pred_glm_pr <- predict(glm_pres,type=c("response"))
      pred_gam_occ <- predict(gam_occ,type=c("response"))
      pred_gam_pr <- predict(gam_pres,type=c("response"))
      pred_rf_occ <- predict(rf_occ,type=c("response"))
      pred_rf_pr <- predict(rf_pres,type=c("response"))
      max_pred_pres <- predict(max_ind_pres, sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")])
      
      sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_pres)  
      
      
      # Determine geographic extent of our data using AOU = i
      max.lat <- ceiling(max(sdm_input$latitude))
      min.lat <- floor(min(sdm_input$latitude))
      max.lon <- ceiling(max(sdm_input$longitude))
      min.lon <- floor(min(sdm_input$longitude))
      
      geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
      
      
      mod.r <- SpatialPointsDataFrame(coords = sdm_output[,c("longitude", "latitude")],
                                      data = sdm_output[,c("latitude", "longitude", "pred_glm_pr", "pred_glm_occ", "pred_gam_pr", "pred_gam_occ",  "pred_rf_pr", "pred_rf_occ","max_pred_pres")], 
                                      proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
      r = raster(mod.r, res = 0.6) # 40x40 km/111 (degrees) * 2 tp eliminate holes
      # bioclim is 4 km
      plot.r = rasterize(mod.r, r)
      
      data(wrld_simpl)
      # Plot the base map
      plot(wrld_simpl, 
           xlim = c(min.lon, max.lon),
           ylim = c(min.lat, max.lat),
           axes = TRUE, 
           col = "grey95", main = paste("SDM plot for ", j, sep=""))
      
      plot(plot.r[[num]], add = TRUE)
      
      # points(x = bbs_sub$longitude, y = bbs_sub$latitude, col = bbs_sub$presence, pch = 20, cex = 0.75)
      # Add the points for individual observation if necessary
      # sdm_input$presence <-droplevels(sdm_input$presence, exclude = c("0"))
      # sdm_input$col = c("black", "white")
      points(x = sdm_input$longitude, y = sdm_input$latitude, col = sdm_input$presence, pch = 20, cex = 0.75)
      
      box()
    }
    
  }
  dev.off()
}
# dev.off()

# using function to generate maps for all methods
mapfun(pdf_name = 'SDM_glm_pres_6280.pdf',vec = "pred_glm_pr", 4)

mapfun(pdf_name = 'SDM_glm_occ_6280.pdf', vec = "pred_glm_occ", 5)

mapfun(pdf_name = 'SDM_gam_pres_6280.pdf',vec = "pred_gam_pr", 6)

mapfun(pdf_name = 'SDM_gam_occ_6280.pdf', vec = "pred_gam_occ", 7)

mapfun(pdf_name = 'SDM_rf_pr_6280.pdf', vec = "pred_rf_pr", 8)

mapfun(pdf_name = 'SDM_rf_occ_6280.pdf', vec = "pred_rf_occ", 9)

mapfun(pdf_name = 'SDM_me_pres_6280.pdf', vec = "pred_me_pres", 10)

setwd("C:/Git/SDMs")

env.data = as.matrix(sdm_input[,c("latitude", "longitude", "elev.mean", "ndvi.mean", "bio.mean.bio1", "bio.mean.bio12")])

env.proj = SpatialPointsDataFrame(coords = sdm_input[,c("longitude", "latitude")],
                                  data = sdm_input[,c("latitude", "longitude", "elev.mean", "ndvi.mean", "bio.mean.bio1", "bio.mean.bio12")], 
                                  proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
r = raster(env.proj)
env.proj.raster = rasterize(env.proj, r)
env.stack = raster::stack(env.proj.raster@data$elev.mean, env.proj.raster@data$ndvi.mean)



# temp filter for vis purposes
sdm_input <- filter(bbs_final_occ_ll, Aou == 6280) %>% left_join(all_env, by = "stateroute") 
sdm_notrans <- filter(sdm_input, occ >= 0.33) 

# Determine geographic extent of our data using AOU = i
max.lat <- ceiling(max(sdm_input$latitude))
min.lat <- floor(min(sdm_input$latitude))
max.lon <- ceiling(max(sdm_input$longitude))
min.lon <- floor(min(sdm_input$longitude))

#   if(levels(as.factor(sdm_input$presence)) > 1){
glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)

# predict
pred_glm_occ <- predict(glm_occ,type=c("response"))
pred_glm_pr <- predict(glm_pres,type=c("response"))

sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ)  

mod.r <- SpatialPointsDataFrame(coords = sdm_output[,c("longitude", "latitude")],
                                data = sdm_output[,c("latitude", "longitude", "pred_glm_pr", "pred_glm_occ")], 
                                proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
r = raster(mod.r, res = 0.6) # 40x40 km/111 (degrees) * 2 tp eliminate holes
# bioclim is 4 km
plot.r = rasterize(mod.r, r)

glm_occ_notrans <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_notrans)
glm_pres_notrans <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_notrans)

# predict
pred_glm_occ_notrans <- predict(glm_occ_notrans,type=c("response"))
pred_glm_pr_notrans <- predict(glm_pres_notrans,type=c("response"))

sdm_output_notrans = cbind(sdm_notrans, pred_glm_pr_notrans, pred_glm_occ_notrans)  

mod.core <- SpatialPointsDataFrame(coords = sdm_output[,c("longitude", "latitude")],
                                   data = sdm_output[,c("latitude", "longitude", "pred_glm_pr", "pred_glm_occ")], 
                                   proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
r.core = raster(mod.core, res = 0.6) # 40x40 km/111 (degrees) * 2 tp eliminate holes
# bioclim is 4 km
plot.core = rasterize(mod.core, r.core)

