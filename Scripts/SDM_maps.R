
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

bbs_occ = read.csv("Data/bbs_sub1.csv", header=TRUE)
bbs_occ_sub = bbs_occ %>% filter(Aou > 2880) %>%
  filter(Aou < 3650 | Aou > 3810) %>%
  filter(Aou < 3900 | Aou > 3910) %>%
  filter(Aou < 4160 | Aou > 4210) %>%
  filter(Aou != 7010)

exp_pres = read.csv("Data/expect_pres.csv", header = TRUE)
exp_pres = exp_pres[,c("stateroute","spAOU")] 
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
bbs_inc_absence = full_join(bbs_occ_sub, exp_pres, by = c("Aou" ="spAOU", "stateroute" = "stateroute"))
bbs_inc_absence$occ[is.na(bbs_inc_absence$occ)] <- 0
bbs_inc_absence$presence = 0
bbs_inc_absence$presence[bbs_inc_absence$occ > 0] <- 1
num_occ = bbs_inc_absence %>% group_by(Aou) %>% tally(presence) %>% left_join(bbs_inc_absence, ., by = "Aou")

# 412 focal species
bbs_final_occ = filter(num_occ,nn > 1)
bbs_occ_code = left_join(bbs_final_occ, tax_code, by = c("Aou" = "AOU_OUT"))

# 319 focal species
bbs_focal_spp = filter(bbs_occ_code, Aou %in% tax_code$AOU_OUT)

bbs_final_occ_ll = left_join(bbs_focal_spp, lat_long, by = "stateroute")
bbs_final_occ_ll = bbs_final_occ_ll[,c("Aou", "stateroute", "occ", "presence", "ALPHA.CODE",
                                       "latitude", "longitude")]
bbs_final_occ_ll$sp_success = 15 * bbs_final_occ_ll$occ
bbs_final_occ_ll$sp_fail = 15 * (1 - bbs_final_occ_ll$occ) 
# temp filter for vis purposes
auc_df = read.csv("Data/auc_df.csv", header = TRUE)

#### change spp here #####
sdm_input <- filter(bbs_final_occ_ll, Aou == 6280) %>% left_join(all_env, by = "stateroute") %>% na.omit(.)
sdm_notrans <- filter(sdm_input, occ > 0.33) %>% na.omit(.)

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
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
r = raster(mod.r, res = 0.6) # 40x40 km/111 (degrees) * 2 tp eliminate holes
# bioclim is 4 km
plot.r = rasterize(mod.r, r)

glm_occ_notrans <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_notrans)
glm_pres_notrans <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_notrans)

# predict
pred_glm_occ_notrans <- predict(glm_occ_notrans,type=c("response"))
pred_glm_pr_notrans <- predict(glm_pres_notrans,type=c("response"))

sdm_output_notrans = cbind(sdm_notrans, pred_glm_pr_notrans, pred_glm_occ_notrans)  

mod.core <- SpatialPointsDataFrame(coords = sdm_output_notrans[,c("longitude", "latitude")],
                                   data = sdm_output_notrans[,c("latitude", "longitude", "pred_glm_pr_notrans", "pred_glm_occ_notrans")], 
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
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
routes_notrans <- st_as_sf(sdm_notrans, coords = c("longitude", "latitude"))

us <- tm_shape(us_sf) + tm_borders() + tm_fill(col = "white")
point_map <- tm_shape(routes_sf) + tm_symbols(size = 0.75, shape="presence", shapes = c(16,4)) + tm_shape(us_sf) + tm_borders( "black", lwd = 3) 
point_map

sdm_occ <- tm_shape(plot.r) + tm_raster("pred_glm_occ", palette = "PRGn", style = "cont") + tm_shape(us_sf) + tm_borders( "black", lwd = 3) 
sdm_occ
sdm_pr <- tm_shape(plot.r) + tm_raster("pred_glm_pr", palette = "PRGn", style = "cont") + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) 
sdm_pr
sdm_core <- tm_shape(plot.core) + tm_raster("pred_glm_pr", palette = "PRGn", style = "cont") + tm_shape(us_sf) + tm_borders(col = "black", lwd = 3) 
sdm_core

tmap_arrange(point_map, sdm_occ, sdm_pr, sdm_core)


