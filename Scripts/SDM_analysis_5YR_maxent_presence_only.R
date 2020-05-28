library(tidyverse)
library(glm2)
library(gam)
library(randomForest)
library(dismo)
library(raster)
library(maptools)
library(pROC)
library(hydroGOF)

bbs_occ = read.csv("Data/bbs_2001_2015.csv", header=TRUE) %>%
  filter(year %in% c(2001, 2004, 2007, 2010, 2013))
# 2006:2010

bbs_occ_sub = bbs_occ %>% 
  dplyr::count(aou, stateroute) %>% 
  dplyr::mutate(occ = n/5) 

exp_pres = read.csv("Data/expect_pres.csv", header = TRUE) %>%
  filter(!stateroute %in% bad_rtes$stateroute)
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


all_env = left_join(bi_means, env_bio_sub, by = "stateroute")

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

bbs_final_occ_ll$sp_success = 5 * bbs_final_occ_ll$occ
bbs_final_occ_ll$sp_fail = 5 * (1 - bbs_final_occ_ll$occ) 
bbs_final_occ_ll$presence <- as.numeric(bbs_final_occ_ll$presence)

# https://www.rdocumentation.org/packages/dismo/versions/1.1-4/topics/maxent
###### SDM analysis ######
auc_df_5 = c()
sp_list = unique(bbs_final_occ_ll$aou)

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
      
      max_ind_pres = dismo::maxent(all_env_raster, ll)
      max_pred_pres <- predict(max_ind_pres, all_env_raster)
      max_pred_points <- raster::extract(max_pred_pres, ll)
      sdm_output = cbind(sdm_input, max_pred_points) 
      
      rmse_me_pres <- rmse(sdm_output$max_pred_points, sdm_output$presence)
      auc_df = rbind(auc_df, c(i, rmse_me_pres))
      j = unique(sdm_input$ALPHA.CODE)
    }
  }
  # write.csv(sdm_output, paste("sdm_output_notrans_", i, ".csv",  sep=""), row.names = FALSE)
}


auc_df = data.frame(auc_df)
names(auc_df) = c("AOU","rmse_me_PO_5")
write.csv(auc_df, "Data/auc_df_ME_only_5.csv", row.names = FALSE)
