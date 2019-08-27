
library(tidyverse)
library(glm2)
library(gam)
library(randomForest)
library(dismo)
library(raster)
library(maptools)
library(pROC)
library(hydroGOF)
library(sf)
library(tmap)

bbs_occ = read.csv("Data/bbs_2001_2015.csv", header=TRUE) %>% filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010)

bbs_occ_sub = bbs_occ %>% 
  dplyr::count(aou, stateroute) %>% 
  filter(n < 16) %>% 
  dplyr::mutate(occ = n/15) 

exp_pres = read.csv("Data/expect_pres.csv", header = TRUE)
traits = read.csv("Data/Master_RO_Correlates.csv", header = TRUE)
bsize = read.csv("data/DunningBodySize_old_2008.11.12.csv", header = TRUE)
lat_long = read.csv("Data/latlongs.csv", header = TRUE)
tax_code = read.csv("Data/Tax_AOU_Alpha.csv", header = TRUE)
bi_env = read.csv("Data/all_env.csv", header = TRUE)
bi_means = bi_env[,c("stateroute","mat.mean", "elev.mean", "map.mean", "ndvi.mean")]
env_bio = read.csv("Data/env_bio.csv", header = TRUE)
env_bio = na.omit(env_bio)
env_bio_sub = env_bio[,c(1, 21:39)]

# env_bio[,c("stateroute","bio.mean.bio4", "bio.mean.bio5", "bio.mean.bio6", "bio.mean.bio13", "bio.mean.bio14")]

all_env = left_join(bi_means, env_bio_sub, by = "stateroute")

#update tax_code Winter Wren
tax_code$AOU_OUT[tax_code$AOU_OUT == 7220] <- 7222
tax_code$AOU_OUT[tax_code$AOU_OUT == 4810] = 4812
tax_code$AOU_OUT[tax_code$AOU_OUT == 4123] = 4120

# BBS cleaning
bbs_inc_absence = full_join(bbs_occ_sub, exp_pres, by = c("aou" ="spAOU", "stateroute" = "stateroute"))
bbs_inc_absence$occ[is.na(bbs_inc_absence$occ)] <- 0
bbs_inc_absence$presence = 0
bbs_inc_absence$presence[bbs_inc_absence$occ > 0] <- 1
num_occ = bbs_inc_absence %>% group_by(aou) %>% tally(presence) %>% left_join(bbs_inc_absence, ., by = "aou")

# 412 focal species
bbs_final_occ = filter(num_occ,n.y > 1)
bbs_occ_code = left_join(bbs_final_occ, tax_code, by = c("aou" = "AOU_OUT"))

# 319 focal species
bbs_final_occ_ll = filter(bbs_occ_code, aou %in% tax_code$AOU_OUT) %>%
  dplyr::select(aou, stateroute, occ, presence, ALPHA.CODE, latitude, longitude)
bbs_final_occ_ll$sp_success = 15 * bbs_final_occ_ll$occ
bbs_final_occ_ll$sp_fail = 15 * (1 - bbs_final_occ_ll$occ) 
bbs_final_occ_ll$presence <- as.numeric(bbs_final_occ_ll$presence)
# temp filter for vis purposes
auc_df = read.csv("Data/auc_df.csv", header = TRUE)

#### change spp here ##### 
sdm_input <- filter(bbs_final_occ_ll, aou == 6280) %>% 
  left_join(all_env, by = "stateroute") %>% 
  na.omit(.)
sdm_notrans <- filter(sdm_input, occ >= 0.33| occ == 0) %>% na.omit(.)

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

sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_pres) 
  

mod.r <- SpatialPointsDataFrame(coords = sdm_output[,c("longitude", "latitude")],
               data = sdm_output[,c("latitude", "longitude", "pred_glm_pr", "pred_glm_occ", "pred_gam_pr", "pred_gam_occ", "pred_rf_occ", "pred_rf_pr", "max_pred_pres")], proj4string = CRS("+proj=longlat +datum=WGS84"))
r = raster(mod.r, res = 0.6) # 40x40 km/111 (degrees) * 2 tp eliminate holes
# bioclim is 4 km
plot.r = rasterize(mod.r, r)

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
sdm_output_notrans = cbind(sdm_notrans, pred_glm_pr_notrans, pred_glm_occ_notrans, pred_gam_pr_notrans, pred_gam_occ_notrans, pred_rf_occ_notrans, pred_rf_pr_notrans, max_pred_pres_notrans) 

mod.core <- SpatialPointsDataFrame(coords = sdm_output_notrans[,c("longitude", "latitude")],
            data = sdm_output_notrans[,c("latitude", "longitude", "pred_glm_pr_notrans", "pred_glm_occ_notrans", "pred_gam_pr_notrans", "pred_gam_occ_notrans", "pred_rf_occ_notrans", "pred_rf_pr_notrans", "max_pred_pres_notrans")], proj4string = CRS("+proj=longlat +datum=WGS84"))
r.core = raster(mod.core, res = 0.6) # 40x40 km/111 (degrees) * 2 tp eliminate holes
# bioclim is 4 km
plot.core = rasterize(mod.core, r.core)

sdm_output$presence <- factor(sdm_output$presence,
                              levels = c(1,0), ordered = TRUE)

us_sf <- read_sf("Z:/GIS/geography/continent.shp")
us_sf <- st_transform(us_sf, crs = "+proj=longlat +datum=WGS84")
# spData::us_states
# read_sf("Z:/GIS/birds/BCR.shp")
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

point_map <- tm_shape(routes_sf) + 
  tm_symbols(size = 0.75, shape="presence", shapes = c(16,4), alpha = 0.5, col = "#5E5E5E") + 
  tm_legend(show=FALSE) +
  tm_shape(us_sf) + tm_borders( "black", lwd = 3) + 
  tm_shape(routes_notrans)  + 
  tm_symbols(col = "presence", palette = "-PRGn", size = 0.75, shapes = c(16,4)) + tm_legend(outside = TRUE)+ 
  tm_layout("  Observed \nOccurences", title.size = 2, title.position = c("right","bottom")) +
  tm_layout(main.title = "A") 

sdm_maxent_pr <- tm_shape(plot.r) + tm_raster("max_pred_pres", palette = "PRGn", style = "cont", title = "MaxEnt Pres", legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout("MaxEnt Pres", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "B") 

sdm_maxent_core <- tm_shape(plot.core) + tm_raster("max_pred_pres_notrans", palette = "PRGn", style = "cont", title = "MaxEnt Core", legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout("MaxEnt Core", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "C") 

sdm_glm_occ <- tm_shape(plot.r) + tm_raster("pred_glm_occ", palette = "PRGn", style = "cont", title = "GLM Occ", legend.show = FALSE) + 
  tm_shape(us_sf) + tm_borders( "black", lwd = 3) + 
  tm_layout("GLM Occ", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "D") 
#sdm_occ

sdm_glm_pr <- tm_shape(plot.r) + tm_raster("pred_glm_pr", palette = "PRGn", style = "cont", title = "GLM Pres", breaks=c(0.05,0.1,0.15 ,0.2,0.25), legend.show = FALSE) + tm_shape(us_sf) + 
  tm_borders(col = "black", lwd = 3) + 
  tm_layout("GLM Pres", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "E") 
#sdm_pr

sdm_glm_core <- tm_shape(plot.core) + tm_raster("pred_glm_pr_notrans", palette = "PRGn", style = "cont", legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + 
  tm_layout("GLM Core", title.size = 2, title.position = c("right","bottom")) +
  tm_layout(main.title = "F") 

sdm_gam_occ <- tm_shape(plot.r) + tm_raster("pred_gam_occ", palette = "PRGn", style = "cont", title = "GAM Occ", legend.show = FALSE) + tm_shape(us_sf) + tm_borders( "black", lwd = 3) + tm_layout("GAM Occ", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "G") 
#sdm_occ

sdm_gam_pr <- tm_shape(plot.r) + tm_raster("pred_gam_pr", palette = "PRGn", style = "cont", breaks=c(0.05,0.1,0.15 ,0.2,0.25), legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + 
  tm_layout("GAM Pres", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "H") 
#sdm_pr

sdm_gam_core <- tm_shape(plot.core) + tm_raster("pred_gam_pr_notrans", palette = "PRGn", style = "cont", legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout("GAM Core", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "I") 

sdm_rf_occ <- tm_shape(plot.r) + tm_raster("pred_rf_occ", palette = "PRGn", style = "cont", title = "RF Occ", legend.show = FALSE) + tm_shape(us_sf) + tm_borders( "black", lwd = 3) + tm_layout("RF Occ", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "J") 
#sdm_occ

sdm_rf_pr <- tm_shape(plot.r) + tm_raster("pred_rf_pr", palette = "PRGn", style = "cont", title = "RF Pres", breaks=c(0.05,0.1,0.15 ,0.2,0.25), legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout("RF Pres", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "K") 
#sdm_pr

sdm_rf_core <- tm_shape(plot.core) + tm_raster("pred_rf_pr_notrans", palette = "PRGn", style = "cont", title = "RF Core", legend.show = FALSE) + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) + tm_layout("RF Core", title.size = 2, title.position = c("right","bottom"), legend.bg.color = "white") +
  tm_layout(main.title = "L") 

MaxEnt_plot <- tmap_arrange(point_map, sdm_maxent_pr, sdm_maxent_core, ncol = 1)
fig_glm <- tmap_arrange(sdm_glm_occ, sdm_glm_pr, sdm_glm_core, ncol = 1)
fig_gam <- tmap_arrange(sdm_gam_occ, sdm_gam_pr, sdm_gam_core, ncol = 1)
fig_rf <- tmap_arrange(sdm_rf_occ, sdm_rf_pr, sdm_rf_core, ncol = 1)

final_fig1 <- tmap_arrange(point_map, sdm_maxent_pr, sdm_maxent_core, sdm_glm_occ, sdm_glm_pr, sdm_glm_core, sdm_gam_occ, sdm_gam_pr, sdm_gam_core, sdm_rf_occ, sdm_rf_pr, sdm_rf_core, nrow = 4, ncol = 3) 
tmap_save(final_fig1, "Figures/Fig1.pdf", height = 16, width = 20)
