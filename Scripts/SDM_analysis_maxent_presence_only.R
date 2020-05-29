library(tidyverse)
library(glm2)
library(gam)
library(randomForest)
library(dismo)
library(raster)
library(maptools)
library(sf)
library(pROC)
library(hydroGOF)

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
bad_rtes = read.csv("Data/bad_rtes.csv", header = TRUE)
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
# write.csv(bbs_final_occ_ll, "Data/bbs_final_occ_ll.csv", row.names = FALSE)

# Thuiller 2014 source for choosing these vars
# http://worldclim.org/bioclim
# BIO4 = temperature seasonality (intra-annual standard deviation * 100)
# BIO5 = maximum temperature of the warmest month
# BIO6 = minimum temperature of the coldest month
# BIO13 = precipitation of the wettest month 
# BIO14 = precipitation of the driest month

# all_bio.pca <- prcomp(env_bio[,c(21:39)], center = TRUE,scale. = TRUE)
# library(ggbiplot)
# ggbiplot(all_bio.pca)

# you need to do cor.test and not excede 5 vars for power. Need  >= 50 presences
sdm_input_global <- left_join(bbs_final_occ_ll, all_env, by = "stateroute")
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_251') # for 64-bit version

###### main analysis ######
all_env_raster <- stack("Z:/GIS/all_env_maxent_mw.tif")
# all_env_raster <- projectRaster(all_env_raster, crs = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"))

auc_df = c()

sp_list = unique(bbs_final_occ_ll$aou)
#  change back when done w overfit

for(i in sp_list){
  sdm_output = c()
  print(i)
  bbs_sub <- filter(bbs_final_occ_ll, aou == i) 
  bbs_new_sub <- filter(bbs_new, aou == i) 
  bbs_new_sub$pres_2016 <- bbs_new_sub$presence
  sdm_input <- filter(all_env, stateroute %in% bbs_sub$stateroute) %>%
    full_join(bbs_sub, by = "stateroute") 
  sdm_input <- data.frame(sdm_input) %>%
    filter(presence == 1)
  if(length(unique(sdm_input$stateroute)) > 40){
    if(nrow(filter(sdm_input, presence == 1)) > 59){
      ll <- data.frame(lon = sdm_input$longitude, lat = sdm_input$latitude)
      
      ll_spat <- SpatialPoints(ll, proj4string=CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
      ll_spat_laea <- spTransform(ll_spat, CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"))
      
      max_ind_pres = dismo::maxent(all_env_raster, ll_spat_laea)
      max_pred_pres <- predict(max_ind_pres, all_env_raster, progress='text')
      
      max_pred_points <- raster::extract(max_pred_pres, ll_spat)
      sdm_output = cbind(sdm_input, max_pred_points) 
      
      rmse_me_pres <- rmse(sdm_output$max_pred_points, sdm_output$presence)
      auc_df = rbind(auc_df, c(i, rmse_me_pres))
      j = unique(sdm_input$ALPHA.CODE)
    }
  }
  # write.csv(sdm_output, paste("sdm_output_notrans_", i, ".csv",  sep=""), row.names = FALSE)
}


auc_df = data.frame(auc_df)
names(auc_df) = c("AOU","rmse_me_pres")
write.csv(auc_df, "Data/auc_df_ME_only.csv", row.names = FALSE)

