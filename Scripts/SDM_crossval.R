library(tidyverse)
library(glm2)
library(gam)
library(randomForest)
library(dismo)
library(raster)
library(maptools)
library(pROC)
library(purrr)
library(hydroGOF)
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_241') # for 64-bit version


bbs_occ = read.csv("Data/bbs_2001_2015.csv", header=TRUE) %>% filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010)

bbs_occ_sub = bbs_occ %>% 
  dplyr::count(aou, stateroute) %>% 
  dplyr::mutate(occ = n/15) 

auc_df <- read.csv("Data/auc_df.csv", header = TRUE)
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
# NOTE - only 1 raster layer, need to rerun to get stack
# all_env_raster <- raster("Z:/GIS/gimms/all_env_maxent.tif")

# read in raw bbs data for 2016
bbs_new <- read.csv("Data/bbs_2016.csv", header = TRUE) 
bbs_new$presence = 1
bbs_new_exp_pres <- read.csv("Data/expect_pres_2016.csv", header = TRUE)
bbs_new_all <- left_join(bbs_new_exp_pres, bbs_new, by = c("spAOU"="aou", "stateroute" = "stateroute"))
bbs_new_all$presence <- case_when(is.na(bbs_new_all$presence) == TRUE ~ 0, 
                                  bbs_new_all$presence == 1 ~ 1)

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
num_occ = bbs_inc_absence %>% 
  group_by(aou) %>% 
  tally(presence) %>% 
  left_join(bbs_inc_absence, ., by = "aou")

# 412 focal species
bbs_final_occ = filter(num_occ,n.y > 1)
bbs_occ_code = left_join(bbs_final_occ, tax_code, by = c("aou" = "AOU_OUT"))

# 319 focal species
bbs_focal_spp = filter(bbs_occ_code, aou %in% tax_code$AOU_OUT & stateroute %in% bbs_occ_sub$stateroute)

bbs_final_occ_ll = left_join(bbs_focal_spp, lat_long, by = "stateroute")
bbs_final_occ_ll = bbs_final_occ_ll[,c("aou", "stateroute", "occ", "presence", "ALPHA.CODE",
                                       "latitude", "longitude")]
bbs_final_occ_ll$sp_success = 15 * bbs_final_occ_ll$occ
bbs_final_occ_ll$sp_fail = 15 * (1 - bbs_final_occ_ll$occ) 


threshfun <- function(pred_vals){ 
 thresh <-  max(pred_vals) * 0.5
}

sp_list = unique(auc_df$AOU)

#### temporal crossval ######
test_df = c()
for(i in sp_list){
  sdm_output = c()
  print(i)
  bbs_sub <- filter(bbs_final_occ_ll, aou == i) # %>% filter(occ <= 0.33333333) RUN FOR EXCL TRANS
  bbs_new_sub <- filter(bbs_new_all, spAOU == i) 
  temp <- filter(all_env, stateroute %in% bbs_sub$stateroute)
  sdm_input <- left_join(bbs_sub, temp, by = "stateroute")
  sdm_input = na.omit(sdm_input)
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
    if(nrow(filter(sdm_input, presence == 1)) > 59){
      
      glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_input)
      gam_pres <- mgcv::gam(presence ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = sdm_input)
      rf_occ <- randomForest(sp_success/15 ~elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      rf_pres <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      
      
      pred_glm_occ <- predict(glm_occ,type=c("response"))
      pred_glm_pr <- predict(glm_pres,type=c("response"))
      pred_gam_occ <- predict(gam_occ,type=c("response"))
      pred_gam_pr <- predict(gam_pres,type=c("response"))
      pred_rf_occ <- predict(rf_occ,type=c("response"))
      pred_rf_pr <- predict(rf_pres,type=c("response"))
      
      ll <- data.frame(lon = sdm_input$longitude, lat = sdm_input$latitude)
      max_ind_pres = dismo::maxent(all_env_raster, ll)
      max_pred_pres <- predict(max_ind_pres, all_env_raster)
      max_pred_points <- raster::extract(max_pred_pres, ll)
      
      threshglm_occ <-threshfun(pred_glm_occ)
      threshglm_pr <-threshfun(pred_glm_pr)
      threshgam_occ <-threshfun(pred_gam_occ)
      threshgam_pr <-threshfun(pred_gam_pr)
      threshrf_occ <-threshfun(pred_rf_occ)
      threshrf_pr <-threshfun(pred_rf_pr)
      threshmax_pres <-threshfun(max_pred_points)
    
      sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_points) 
      pred_2016 <- left_join(sdm_output, bbs_new_sub[c("stateroute", "presence")], by = "stateroute") %>%
        mutate(predicted_glm_occ = ifelse(pred_glm_occ > threshglm_occ, 1, 0),
               predicted_glm_pr = ifelse(pred_glm_pr > threshglm_pr, 1, 0),
               predicted_gam_occ = ifelse(pred_gam_occ > threshgam_occ, 1, 0),
               predicted_gam_pr = ifelse(pred_gam_pr > threshgam_pr, 1, 0),
               predicted_rf_occ = ifelse(pred_rf_occ > threshrf_occ, 1, 0),
               predicted_rf_pr = ifelse(pred_rf_pr > threshrf_pr, 1, 0),
               predicted_max_pres = ifelse(max_pred_points > threshmax_pres, 1, 0)) %>%
        dplyr::select(aou, stateroute, occ, presence.x, latitude, longitude, pred_gam_occ, presence.y, predicted_glm_occ, predicted_glm_pr,predicted_gam_occ, predicted_gam_pr, predicted_rf_occ, predicted_rf_pr, predicted_max_pres)
      test_df = rbind(test_df, pred_2016)
    }
  }
}


# write.csv(test_df, "Data/temporal_crossval_df_5.csv", row.names = FALSE) # wrote _5 for thresh of .5, med for median

##### temporal processing ####
test_df <- read.csv("Data/temporal_crossval_df_5.csv", header = TRUE) 
# to account for species not detected in 2015-2016 but are within the range
# presence.x is 2001-2015, presence.y is 2015-2016
test_df$presence.y[is.na(test_df$presence.y)] = 0 

# CONFUSION MATRIX is PREDICTED X ACTUAL
glmocc <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(glmocc = purrr::map(data, ~{
      dat <- .
      dat %>%
        group_by(predicted_glm_occ, presence.y) %>%
        count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(glmocc)) %>%
  mutate(cats = paste(predicted_glm_occ, presence.y, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_glmocc",
                              cats == "1_0" ~ "pres_abs_glmocc",
                              cats == "0_1" ~ "abs_pres_glmocc",
                              cats == "0_0" ~ "abs_abs_glmocc")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)


glmpr <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(glmpr = purrr::map(data, ~{
  newdat <- .
  newdat %>%
    group_by(predicted_glm_pr, presence.y) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(glmpr)) %>%
  mutate(cats = paste(predicted_glm_pr, presence.y, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_glmpr",
                              cats == "1_0" ~ "pres_abs_glmpr",
                              cats == "0_1" ~ "abs_pres_glmpr",
                              cats == "0_0" ~ "abs_abs_glmpr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

gamocc <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(gamocc = purrr::map(data, ~{
  newdat2 <- .
  newdat2 %>%
    group_by(predicted_gam_occ, presence.y) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(gamocc)) %>%
  mutate(cats = paste(predicted_gam_occ, presence.y, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_gamocc",
                              cats == "1_0" ~ "pres_abs_gamocc",
                              cats == "0_1" ~ "abs_pres_gamocc",
                              cats == "0_0" ~ "abs_abs_gamocc")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

gampr <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(gampr = purrr::map(data, ~{
  newdat3 <- .
  newdat3 %>%
    group_by(predicted_gam_pr, presence.y) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(gampr)) %>%
  mutate(cats = paste(predicted_gam_pr, presence.y, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_gampr",
                              cats == "1_0" ~ "pres_abs_gampr",
                              cats == "0_1" ~ "abs_pres_gampr",
                              cats == "0_0" ~ "abs_abs_gampr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

rfocc <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(rfocc = purrr::map(data, ~{
  newdat4 <- .
  newdat4 %>%
    group_by(predicted_rf_occ, presence.y) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(rfocc)) %>%
  mutate(cats = paste(predicted_rf_occ, presence.y, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_rfocc",
                              cats == "1_0" ~ "pres_abs_rfocc",
                              cats == "0_1" ~ "abs_pres_rfocc",
                              cats == "0_0" ~ "abs_abs_rfocc")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

rfpr <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(rfpres = purrr::map(data, ~{
  newdat5 <- .
  newdat5 %>%
    group_by(predicted_rf_pr, presence.y) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(rfpres)) %>%
  mutate(cats = paste(predicted_rf_pr, presence.y, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_rfpr",
                              cats == "1_0" ~ "pres_abs_rfpr",
                              cats == "0_1" ~ "abs_pres_rfpr",
                              cats == "0_0" ~ "abs_abs_rfpr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

maxpr <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(maxpres = purrr::map(data, ~{
     newdat6 <- .
     newdat6 %>%
       group_by(predicted_max_pres, presence.y) %>%
       count()})) %>%
  dplyr::select(aou, maxpres) %>%
  unnest(cols = c(maxpres)) %>%
  mutate(predicted_max_pres = if_else(is.na(predicted_max_pres)|predicted_max_pres == 0, 0, 1),
         cats = paste(predicted_max_pres, presence.y, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_maxpr",
                              cats == "1_0" ~ "pres_abs_maxpr",
                              cats == "0_1" ~ "abs_pres_maxpr",
                              cats == "0_0" ~ "abs_abs_maxpr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n,  values_fill = list(n = 0))

length <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(length = purrr::map(data, ~{
  newdatlength <- .
  length(newdatlength$stateroute)})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(length))

pres_matrix <- full_join(glmocc, glmpr, by = "aou") %>%
  full_join(gamocc, by = "aou") %>%
  full_join(gampr, by = "aou") %>%
  full_join(rfocc, by = "aou") %>%
  full_join(rfpr, by = "aou") %>%
  full_join(maxpr, by = "aou") %>%
  full_join(length, by = "aou") %>%
  replace_na(list(abs_abs_glmocc=0, abs_pres_glmocc=0,  pres_abs_glmocc=0,  pres_pres_glmocc=0, abs_abs_glmpr=0,   
  abs_pres_glmpr=0, pres_abs_glmpr=0, pres_pres_glmpr=0, abs_abs_gamocc=0, abs_pres_gamocc=0, pres_abs_gamocc=0, pres_pres_gamocc=0, abs_abs_gampr=0, abs_pres_gampr=0, pres_abs_gampr=0, pres_pres_gampr=0, abs_abs_rfocc=0, abs_pres_rfocc=0, pres_abs_rfocc=0, pres_pres_rfocc=0, abs_abs_rfpr=0, abs_pres_rfpr=0, pres_abs_rfpr=0, pres_pres_rfpr=0, abs_abs_maxpr=0, abs_pres_maxpr=0, pres_abs_maxpr=0, pres_pres_maxpr=0))

pres_matrix <- data.frame(pres_matrix)

pres_matrix$sensitivity_glmocc <- (pres_matrix$pres_pres_glmocc/(pres_matrix$pres_pres_glmocc + pres_matrix$abs_pres_glmocc))*100
pres_matrix$specificity_glmocc <- (pres_matrix$abs_abs_glmocc/(pres_matrix$pres_abs_glmocc + pres_matrix$abs_abs_glmocc))*100
pres_matrix$pp_glmocc <- (pres_matrix$pres_pres_glmocc/(pres_matrix$pres_pres_glmocc + pres_matrix$pres_abs_glmocc))*100
pres_matrix$np_glmocc <- (pres_matrix$abs_abs_glmocc/(pres_matrix$abs_abs_glmocc + pres_matrix$abs_pres_glmocc))*100


pres_matrix$sensitivity_glmpr <- (pres_matrix$pres_pres_glmpr/(pres_matrix$pres_pres_glmpr + pres_matrix$abs_pres_glmpr))*100
pres_matrix$specificity_glmpr <- (pres_matrix$abs_abs_glmpr/(pres_matrix$pres_abs_glmpr + pres_matrix$abs_abs_glmpr))*100
pres_matrix$pp_glmpr <- (pres_matrix$pres_pres_glmpr/(pres_matrix$pres_pres_glmpr + pres_matrix$pres_abs_glmpr))*100
pres_matrix$np_glmpr <- (pres_matrix$abs_abs_glmpr/(pres_matrix$abs_abs_glmpr + pres_matrix$abs_pres_glmpr))*100


pres_matrix$sensitivity_gamocc <- (pres_matrix$pres_pres_gamocc/(pres_matrix$pres_pres_gamocc + pres_matrix$abs_pres_gamocc))*100
pres_matrix$specificity_gamocc <- (pres_matrix$abs_abs_gamocc/(pres_matrix$pres_abs_gamocc + pres_matrix$abs_abs_gamocc))*100
pres_matrix$pp_gamocc <- (pres_matrix$pres_pres_gamocc/(pres_matrix$pres_pres_gamocc + pres_matrix$pres_abs_gamocc))*100
pres_matrix$np_gamocc <- (pres_matrix$abs_abs_gamocc/(pres_matrix$abs_abs_gamocc + pres_matrix$abs_pres_gamocc))*100


pres_matrix$sensitivity_gampr <- (pres_matrix$pres_pres_gampr/(pres_matrix$pres_pres_gampr + pres_matrix$abs_pres_gampr))*100
pres_matrix$specificity_gampr <- (pres_matrix$abs_abs_gampr/(pres_matrix$pres_abs_gampr + pres_matrix$abs_abs_gampr))*100
pres_matrix$pp_gampr <- (pres_matrix$pres_pres_gampr/(pres_matrix$pres_pres_gampr + pres_matrix$pres_abs_gampr))*100
pres_matrix$np_gampr <- (pres_matrix$abs_abs_gampr/(pres_matrix$abs_abs_gampr + pres_matrix$abs_pres_gampr))*100


pres_matrix$sensitivity_rfocc <- (pres_matrix$pres_pres_rfocc/(pres_matrix$pres_pres_rfocc + pres_matrix$abs_pres_rfocc))*100
pres_matrix$specificity_rfocc <- (pres_matrix$abs_abs_rfocc/(pres_matrix$pres_abs_rfocc + pres_matrix$abs_abs_rfocc))*100
pres_matrix$pp_rfocc <- (pres_matrix$pres_pres_rfocc/(pres_matrix$pres_pres_rfocc + pres_matrix$pres_abs_rfocc))*100
pres_matrix$np_rfocc <- (pres_matrix$abs_abs_rfocc/(pres_matrix$abs_abs_rfocc + pres_matrix$abs_pres_rfocc))*100


pres_matrix$sensitivity_rfpr <- (pres_matrix$pres_pres_rfpr/(pres_matrix$pres_pres_rfpr + pres_matrix$abs_pres_rfpr))*100
pres_matrix$specificity_rfpr <- (pres_matrix$abs_abs_rfpr/(pres_matrix$pres_abs_rfpr + pres_matrix$abs_abs_rfpr))*100
pres_matrix$pp_rfpr <- (pres_matrix$pres_pres_rfpr/(pres_matrix$pres_pres_rfpr + pres_matrix$pres_abs_rfpr))*100
pres_matrix$np_rfpr <- (pres_matrix$abs_abs_rfpr/(pres_matrix$abs_abs_rfpr + pres_matrix$abs_pres_rfpr))*100


pres_matrix$sensitivity_max <- (pres_matrix$pres_pres_max/(pres_matrix$pres_pres_max + pres_matrix$abs_pres_max))*100
pres_matrix$specificity_max <- (pres_matrix$abs_abs_max/(pres_matrix$pres_abs_max + pres_matrix$abs_abs_max))*100
pres_matrix$pp_max <- (pres_matrix$pres_pres_max/(pres_matrix$pres_pres_max + pres_matrix$pres_abs_max))*100
pres_matrix$np_max <- (pres_matrix$abs_abs_max/(pres_matrix$abs_abs_max + pres_matrix$abs_pres_max))*100

 
pres_matrix_means <- pres_matrix %>%
  filter(!is.nan(np_max)) %>% # removed NaNs, there were 100% predicted presence for these spp.
  filter(!is.nan(np_glmpr)) %>%
  filter(!is.nan(np_rfpr)) %>%
  summarise(n = n(),
  mean_sensitivity_gamocc = mean(sensitivity_gamocc),
  mean_specificity_gamocc = mean(specificity_gamocc),
  mean_pp_gamocc = mean(pp_gamocc),          
  mean_np_gamocc = mean(np_gamocc),   
 
  mean_sensitivity_gampr = mean(sensitivity_gampr),
  mean_specificity_gampr = mean(specificity_gampr),
  mean_pp_gampr = mean(pp_gampr),         
  mean_np_gampr = mean(np_gampr), 
  
  mean_sensitivity_glmocc = mean(sensitivity_glmocc),
  mean_specificity_glmocc = mean(specificity_glmocc),
  mean_pp_glmocc = mean(pp_glmocc),    
  mean_np_glmocc = mean(np_glmocc),
    
  mean_sensitivity_glmpr = mean(sensitivity_glmpr),
  mean_specificity_glmpr = mean(specificity_glmpr),
  mean_pp_glmpr = mean(pp_glmpr),      
  mean_np_glmpr = mean(np_glmpr), 
  
 mean_sensitivity_rfocc = mean(sensitivity_rfocc),
 mean_specificity_rfocc = mean(specificity_rfocc),
 mean_pp_rfocc  = mean(pp_rfocc),        
 mean_np_rfocc = mean(np_rfocc), 

 mean_sensitivity_rfpr = mean(sensitivity_rfpr),
 mean_specificity_rfpr = mean(specificity_rfpr), 
 mean_pp_rfpr = mean(pp_rfpr),        
 mean_np_rfpr = mean(np_rfpr), 
  
 mean_sensitivity_max  = mean(sensitivity_max), 
 mean_specificity_max  = mean(specificity_max), 
 mean_pp_max = mean(pp_max),           
 mean_np_max = mean(np_max),
 
 sd_sensitivity_gamocc = sd(sensitivity_gamocc),
 sd_specificity_gamocc = sd(specificity_gamocc),
 sd_pp_gamocc = sd(pp_gamocc),          
 sd_np_gamocc = sd(np_gamocc),
 
 sd_sensitivity_gampr = sd(sensitivity_gampr),
 sd_specificity_gampr = sd(specificity_gampr),
 sd_pp_gampr = sd(pp_gampr),         
 sd_np_gampr = sd(np_gampr),
 
 sd_sensitivity_glmocc = sd(sensitivity_glmocc),
 sd_specificity_glmocc = sd(specificity_glmocc),
 sd_pp_glmocc = sd(pp_glmocc),    
 sd_np_glmocc = sd(np_glmocc), 
 
 sd_sensitivity_glmpr = sd(sensitivity_glmpr),
 sd_specificity_glmpr = sd(specificity_glmpr),
 sd_pp_glmpr = sd(pp_glmpr),      
 sd_np_glmpr = sd(np_glmpr),
 
 sd_sensitivity_rfocc = sd(sensitivity_rfocc),
 sd_specificity_rfocc = sd(specificity_rfocc),
 sd_pp_rfocc  = sd(pp_rfocc),        
 sd_np_rfocc = sd(np_rfocc),
 
 sd_sensitivity_rfpr = sd(sensitivity_rfpr),
 sd_specificity_rfpr = sd(specificity_rfpr), 
 sd_pp_rfpr = sd(pp_rfpr),        
 sd_np_rfpr = sd(np_rfpr),
 
 sd_sensitivity_max  = sd(sensitivity_max), 
 sd_specificity_max  = sd(specificity_max), 
 sd_pp_max = sd(pp_max),           
 sd_np_max = sd(np_max))
  

pres_matrix_mean <- pivot_longer(pres_matrix_means, mean_sensitivity_gamocc:mean_np_max, "Mod") %>%
  mutate(complement = value,
         mod = substring(Mod, 6)) %>%
  dplyr::select(n, mod, Mod, complement) 
pres_matrix_sd <- pivot_longer(pres_matrix_means, sd_sensitivity_gamocc:sd_np_max, "Mod", "SD") %>%
  mutate(SD = value, 
         mod = substring(Mod, 4)) %>%
  dplyr::select(n, mod, Mod, SD) %>%
  mutate(se = SD/sqrt(n))

pres_matrix_plot2 <-  full_join(pres_matrix_mean, pres_matrix_sd, by = c("n", "mod")) %>%
  separate(col = Mod.x, into = c("Mean","Measure", "Modtype"), sep = "_")  %>%
  mutate(value = 100 - complement)
pres_matrix_plot2$Modtype <- factor(pres_matrix_plot2$Modtype, levels = c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"), ordered = TRUE)

pres_matrix_plot2$Measure <- factor(pres_matrix_plot2$Measure, levels = c("pp","np","sensitivity","specificity") , ordered = TRUE)
#+ scale_color_manual(values=c("#034e7b","#034e7b","steelblue2", "steelblue2","#238b45", "#238b45" ,"purple"), labels=c("rmse_gam", "rmse_gam_pres", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres",  "rmse_me_pres"))

pplot = pres_matrix_plot2 %>%
  ggplot(aes(x = Measure, y = value, group = Modtype)) +   
  geom_bar(aes(fill = factor(Modtype)), position='dodge', stat="identity", width=0.75,color = "white", lwd = 3) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), position = position_dodge(0.75), width = 0.2) +
  theme_classic()+ 
  theme(axis.title.x=element_text(size=54),axis.title.y=element_text(size=54, angle=90)) + xlab(bquote("")) + ylab(bquote("Percent")) +
  scale_fill_manual(values = c("#034e7b","navyblue","steelblue2", "dodgerblue2","#238b45", "darkgreen" ,"purple"),
                    breaks=c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"),
                    labels=c("GAM - Occ","GAM - Pr","GLM - Occ","GLM - Pr","RF - Occ","RF - Pr","MaxEnt - Pr")) +
  scale_x_discrete(labels=c("False \ndiscovery rate","False \nomission rate", "False \nnegative rate", "False \npositive rate")) +
  theme(axis.text.x=element_text(size = 50, colour = "black"),axis.ticks=element_blank(), axis.text.y=element_text(size=50, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=50), legend.key.width=unit(2, "lines"), legend.key.size = unit(2, "cm")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) 
# ggsave("Figures/temp_crossval_5.pdf", width = 30, height = 20)


##### spatial_crossval #######
bbs_env <- left_join(bbs_final_occ_ll, all_env, by = "stateroute")
rmse = c()
sdm_space_cval = c()
for(i in sp_list[124:191]){
  print(i)
  space_sub <- dplyr::filter(bbs_env,  aou == i)
  #Randomly shuffle the data
  sdm_input<-space_sub[sample(nrow(space_sub)),]
  sdm_input = na.omit(sdm_input)
  #Create 10 equally size folds
  folds <- cut(seq(1,nrow(sdm_input)),breaks=10,labels=FALSE)
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
    if(nrow(filter(sdm_input, presence == 1)) > 49){
      #Perform 10 fold cross validation
      for(j in 1:10){
        print(j)
        #Segement your data by fold using the which() function 
        testIndexes <- which(folds==j,arr.ind=TRUE)
        testData <- sdm_input[testIndexes, ]
        trainData <- sdm_input[-testIndexes, ]
        gam_occ_train <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = trainData)
        glm_occ_train <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = trainData)
        glm_pres_train <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = trainData)
        gam_occ_train <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = trainData)
        gam_pres_train <- mgcv::gam(presence ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = trainData)
        rf_occ_train <- randomForest(sp_success/15 ~elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = trainData)
        rf_pres_train <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = trainData)
        
        pred_gam_occ <- predict(gam_occ_train, testData,type= c("response")) 
        pred_gam_pr <- predict(gam_pres_train,  testData,type= c("response"))
        pred_glm_occ <- predict(glm_occ_train, testData,type= c("response"))
        pred_glm_pr <- predict(glm_pres_train, testData,type= c("response"))
        pred_rf_occ <- predict(rf_occ_train, testData,type= c("response"))
        pred_rf_pr <- predict(rf_pres_train, testData,type= c("response"))

        lltrain <- data.frame(lon = trainData$longitude, lat = trainData$latitude)
        lltest <- data.frame(lon = testData$longitude, lat = testData$latitude)
        max_ind_pres = dismo::maxent(all_env_raster, lltrain)
        all_env_points <- raster::extract(all_env_raster, lltest)
        max_pred_points <- predict(max_ind_pres, all_env_points,type= c("response"))
        
        threshglm_occ <-threshfun(pred_glm_occ[!is.na(pred_glm_occ)])
        threshglm_pr <-threshfun(pred_glm_pr[!is.na(pred_glm_pr)])
        threshgam_occ <-threshfun(pred_gam_occ[!is.na(pred_gam_occ)])
        threshgam_pr <-threshfun(pred_gam_pr[!is.na(pred_gam_pr)])
        threshrf_occ <-threshfun(pred_rf_occ[!is.na(pred_rf_occ)])
        threshrf_pr <-threshfun(pred_rf_pr[!is.na(pred_rf_pr)])
        threshmax_pres <-threshfun(max_pred_points[!is.na(max_pred_points)])
      
        sdm_test <- sdm_input[testIndexes, ] %>%  
          mutate(predicted_glm_occ = ifelse(pred_glm_occ > threshglm_occ, 1, 0),
                 predicted_glm_pr = ifelse(pred_glm_pr > threshglm_pr, 1, 0),
                 predicted_gam_occ = ifelse(pred_gam_occ > threshgam_occ, 1, 0),
                 predicted_gam_pr = ifelse(pred_gam_pr > threshgam_pr, 1, 0),
                 predicted_rf_occ = ifelse(pred_rf_occ > threshrf_occ, 1, 0),
                 predicted_rf_pr = ifelse(pred_rf_pr > threshrf_pr, 1, 0),
                 predicted_max_pres = ifelse(max_pred_points > threshmax_pres, 1, 0),
                 j = j) 
        rmse_occ <- rmse(pred_glm_occ, sdm_test$occ)
        rmse_pres <- rmse(pred_glm_pr, sdm_test$presence)
        rmse_gam <- rmse(as.vector(pred_gam_occ), sdm_test$occ)
        rmse_gam_pres <- rmse(as.vector(pred_gam_pr), sdm_test$presence)
        rmse_rf <- rmse(pred_rf_occ, sdm_test$occ)
        rmse_rf_pres <- rmse(as.vector(as.numeric(pred_rf_pr)), sdm_test$presence)
        rmse_me_pres <- rmse(max_pred_points, sdm_test$presence)
        rmse = rbind(rmse, c(i, rmse_occ, rmse_pres, rmse_gam, rmse_gam_pres, rmse_rf, rmse_rf_pres, rmse_me_pres))
        
        sdm_space_cval <- rbind(sdm_space_cval, sdm_test)
      }
    }
  }
}
sdm_space_cval <- data.frame(sdm_space_cval)
# write.csv(sdm_space_cval,"Data/space_cval_5.csv", row.names = FALSE)
rmse <- data.frame(rmse)
names(rmse) <- c("aou","rmse_occ", "rmse_pres", "rmse_gam", "rmse_gam_pres", "rmse_rf", "rmse_rf_pres", "rmse_me_pres")
# write.csv(rmse,"Data/space_cval_rmse.csv", row.names = FALSE)
##### process spatial data #####
sdm_space_cval <- read.csv("Data/space_cval_5.csv", header = TRUE) 

# CONFUSION MATRIX is PREDICTED X ACTUAL
glmocc <- sdm_space_cval %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(glmocc = purrr::map(data, ~{
    dat <- .
    dat %>%
      group_by(predicted_glm_occ, presence) %>%
      count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(glmocc)) %>%
  mutate(cats = paste(predicted_glm_occ, presence, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_glmocc",
                              cats == "1_0" ~ "pres_abs_glmocc",
                              cats == "0_1" ~ "abs_pres_glmocc",
                              cats == "0_0" ~ "abs_abs_glmocc")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)


glmpr <- sdm_space_cval %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(glmpr = purrr::map(data, ~{
    newdat <- .
    newdat %>%
      group_by(predicted_glm_pr, presence) %>%
      count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(glmpr)) %>%
  mutate(cats = paste(predicted_glm_pr, presence, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_glmpr",
                              cats == "1_0" ~ "pres_abs_glmpr",
                              cats == "0_1" ~ "abs_pres_glmpr",
                              cats == "0_0" ~ "abs_abs_glmpr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

gamocc <- sdm_space_cval %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(gamocc = purrr::map(data, ~{
    newdat2 <- .
    newdat2 %>%
      group_by(predicted_gam_occ, presence) %>%
      count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(gamocc)) %>%
  mutate(cats = paste(predicted_gam_occ, presence, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_gamocc",
                              cats == "1_0" ~ "pres_abs_gamocc",
                              cats == "0_1" ~ "abs_pres_gamocc",
                              cats == "0_0" ~ "abs_abs_gamocc")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

gampr <- sdm_space_cval %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(gampr = purrr::map(data, ~{
    newdat3 <- .
    newdat3 %>%
      group_by(predicted_gam_pr, presence) %>%
      count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(gampr)) %>%
  mutate(cats = paste(predicted_gam_pr, presence, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_gampr",
                              cats == "1_0" ~ "pres_abs_gampr",
                              cats == "0_1" ~ "abs_pres_gampr",
                              cats == "0_0" ~ "abs_abs_gampr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

rfocc <- sdm_space_cval %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(rfocc = purrr::map(data, ~{
    newdat4 <- .
    newdat4 %>%
      group_by(predicted_rf_occ, presence) %>%
      count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(rfocc)) %>%
  mutate(cats = paste(predicted_rf_occ, presence, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_rfocc",
                              cats == "1_0" ~ "pres_abs_rfocc",
                              cats == "0_1" ~ "abs_pres_rfocc",
                              cats == "0_0" ~ "abs_abs_rfocc")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

rfpr <- sdm_space_cval %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(rfpres = purrr::map(data, ~{
    newdat5 <- .
    newdat5 %>%
      group_by(predicted_rf_pr, presence) %>%
      count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(rfpres)) %>%
  mutate(cats = paste(predicted_rf_pr, presence, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_rfpr",
                              cats == "1_0" ~ "pres_abs_rfpr",
                              cats == "0_1" ~ "abs_pres_rfpr",
                              cats == "0_0" ~ "abs_abs_rfpr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

maxpr <- sdm_space_cval %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(maxpres = purrr::map(data, ~{
    newdat6 <- .
    newdat6 %>%
      group_by(predicted_max_pres, presence) %>%
      count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(maxpres)) %>%
  mutate(cats = paste(predicted_max_pres, presence, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_maxpr",
                              cats == "1_0" ~ "pres_abs_maxpr",
                              cats == "0_1" ~ "abs_pres_maxpr",
                              cats == "0_0" ~ "abs_abs_maxpr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

length <- sdm_space_cval %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(length = purrr::map(data, ~{
    newdatlength <- .
    length(newdatlength$stateroute)})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(length))

pres_matrix <- full_join(glmocc, glmpr, by = "aou") %>%
  full_join(gamocc, by = "aou") %>%
  full_join(gampr, by = "aou") %>%
  full_join(rfocc, by = "aou") %>%
  full_join(rfpr, by = "aou") %>%
  full_join(maxpr, by = "aou") %>%
  full_join(length, by = "aou") %>%
  replace_na(list(abs_abs_glmocc=0, abs_pres_glmocc=0,  pres_abs_glmocc=0,  pres_pres_glmocc=0, abs_abs_glmpr=0,   
                  abs_pres_glmpr=0, pres_abs_glmpr=0, pres_pres_glmpr=0, abs_abs_gamocc=0, abs_pres_gamocc=0, pres_abs_gamocc=0, pres_pres_gamocc=0, abs_abs_gampr=0, abs_pres_gampr=0, pres_abs_gampr=0, pres_pres_gampr=0, abs_abs_rfocc=0, abs_pres_rfocc=0, pres_abs_rfocc=0, pres_pres_rfocc=0, abs_abs_rfpr=0, abs_pres_rfpr=0, pres_abs_rfpr=0, pres_pres_rfpr=0, abs_abs_maxpr=0, abs_pres_maxpr=0, pres_abs_maxpr=0, pres_pres_maxpr=0))

pres_matrix <- data.frame(pres_matrix)

pres_matrix$sensitivity_glmocc <- (pres_matrix$pres_pres_glmocc/(pres_matrix$pres_pres_glmocc + pres_matrix$abs_pres_glmocc))*100
pres_matrix$specificity_glmocc <- (pres_matrix$abs_abs_glmocc/(pres_matrix$pres_abs_glmocc + pres_matrix$abs_abs_glmocc))*100
pres_matrix$pp_glmocc <- (pres_matrix$pres_pres_glmocc/(pres_matrix$pres_pres_glmocc + pres_matrix$pres_abs_glmocc))*100
pres_matrix$np_glmocc <- (pres_matrix$abs_abs_glmocc/(pres_matrix$abs_abs_glmocc + pres_matrix$abs_pres_glmocc))*100


pres_matrix$sensitivity_glmpr <- (pres_matrix$pres_pres_glmpr/(pres_matrix$pres_pres_glmpr + pres_matrix$abs_pres_glmpr))*100
pres_matrix$specificity_glmpr <- (pres_matrix$abs_abs_glmpr/(pres_matrix$pres_abs_glmpr + pres_matrix$abs_abs_glmpr))*100
pres_matrix$pp_glmpr <- (pres_matrix$pres_pres_glmpr/(pres_matrix$pres_pres_glmpr + pres_matrix$pres_abs_glmpr))*100
pres_matrix$np_glmpr <- (pres_matrix$abs_abs_glmpr/(pres_matrix$abs_abs_glmpr + pres_matrix$abs_pres_glmpr))*100


pres_matrix$sensitivity_gamocc <- (pres_matrix$pres_pres_gamocc/(pres_matrix$pres_pres_gamocc + pres_matrix$abs_pres_gamocc))*100
pres_matrix$specificity_gamocc <- (pres_matrix$abs_abs_gamocc/(pres_matrix$pres_abs_gamocc + pres_matrix$abs_abs_gamocc))*100
pres_matrix$pp_gamocc <- (pres_matrix$pres_pres_gamocc/(pres_matrix$pres_pres_gamocc + pres_matrix$pres_abs_gamocc))*100
pres_matrix$np_gamocc <- (pres_matrix$abs_abs_gamocc/(pres_matrix$abs_abs_gamocc + pres_matrix$abs_pres_gamocc))*100


pres_matrix$sensitivity_gampr <- (pres_matrix$pres_pres_gampr/(pres_matrix$pres_pres_gampr + pres_matrix$abs_pres_gampr))*100
pres_matrix$specificity_gampr <- (pres_matrix$abs_abs_gampr/(pres_matrix$pres_abs_gampr + pres_matrix$abs_abs_gampr))*100
pres_matrix$pp_gampr <- (pres_matrix$pres_pres_gampr/(pres_matrix$pres_pres_gampr + pres_matrix$pres_abs_gampr))*100
pres_matrix$np_gampr <- (pres_matrix$abs_abs_gampr/(pres_matrix$abs_abs_gampr + pres_matrix$abs_pres_gampr))*100


pres_matrix$sensitivity_rfocc <- (pres_matrix$pres_pres_rfocc/(pres_matrix$pres_pres_rfocc + pres_matrix$abs_pres_rfocc))*100
pres_matrix$specificity_rfocc <- (pres_matrix$abs_abs_rfocc/(pres_matrix$pres_abs_rfocc + pres_matrix$abs_abs_rfocc))*100
pres_matrix$pp_rfocc <- (pres_matrix$pres_pres_rfocc/(pres_matrix$pres_pres_rfocc + pres_matrix$pres_abs_rfocc))*100
pres_matrix$np_rfocc <- (pres_matrix$abs_abs_rfocc/(pres_matrix$abs_abs_rfocc + pres_matrix$abs_pres_rfocc))*100


pres_matrix$sensitivity_rfpr <- (pres_matrix$pres_pres_rfpr/(pres_matrix$pres_pres_rfpr + pres_matrix$abs_pres_rfpr))*100
pres_matrix$specificity_rfpr <- (pres_matrix$abs_abs_rfpr/(pres_matrix$pres_abs_rfpr + pres_matrix$abs_abs_rfpr))*100
pres_matrix$pp_rfpr <- (pres_matrix$pres_pres_rfpr/(pres_matrix$pres_pres_rfpr + pres_matrix$pres_abs_rfpr))*100
pres_matrix$np_rfpr <- (pres_matrix$abs_abs_rfpr/(pres_matrix$abs_abs_rfpr + pres_matrix$abs_pres_rfpr))*100


pres_matrix$sensitivity_max <- (pres_matrix$pres_pres_max/(pres_matrix$pres_pres_max + pres_matrix$abs_pres_max))*100
pres_matrix$specificity_max <- (pres_matrix$abs_abs_max/(pres_matrix$pres_abs_max + pres_matrix$abs_abs_max))*100
pres_matrix$pp_max <- (pres_matrix$pres_pres_max/(pres_matrix$pres_pres_max + pres_matrix$pres_abs_max))*100
pres_matrix$np_max <- (pres_matrix$abs_abs_max/(pres_matrix$abs_abs_max + pres_matrix$abs_pres_max))*100


pres_matrix_means <- pres_matrix %>%
  filter(!is.nan(np_max)) %>% # removed NaNs, there were 100% predicted presence for these spp.
  filter(!is.nan(np_glmpr)) %>%
  filter(!is.nan(np_rfpr)) %>%
  summarise(n = n(),
    mean_sensitivity_gamocc = mean(sensitivity_gamocc),
    mean_specificity_gamocc = mean(specificity_gamocc),
    mean_pp_gamocc = mean(pp_gamocc),          
    mean_np_gamocc = mean(np_gamocc),       
    
    mean_sensitivity_gampr = mean(sensitivity_gampr),
    mean_specificity_gampr = mean(specificity_gampr),
    mean_pp_gampr = mean(pp_gampr),         
    mean_np_gampr = mean(np_gampr),         
    
    mean_sensitivity_glmocc = mean(sensitivity_glmocc),
    mean_specificity_glmocc = mean(specificity_glmocc),
    mean_pp_glmocc = mean(pp_glmocc),    
    mean_np_glmocc = mean(np_glmocc),     
    
    mean_sensitivity_glmpr = mean(sensitivity_glmpr),
    mean_specificity_glmpr = mean(specificity_glmpr),
    mean_pp_glmpr = mean(pp_glmpr),      
    mean_np_glmpr = mean(np_glmpr),         
    
    mean_sensitivity_rfocc = mean(sensitivity_rfocc),
    mean_specificity_rfocc = mean(specificity_rfocc),
    mean_pp_rfocc  = mean(pp_rfocc),        
    mean_np_rfocc = mean(np_rfocc),        
    
    mean_sensitivity_rfpr = mean(sensitivity_rfpr),
    mean_specificity_rfpr = mean(specificity_rfpr), 
    mean_pp_rfpr = mean(pp_rfpr),        
    mean_np_rfpr = mean(np_rfpr),           
    
    mean_sensitivity_max  = mean(sensitivity_max), 
    mean_specificity_max  = mean(specificity_max), 
    mean_pp_max = mean(pp_max),           
    mean_np_max = mean(np_max),
    
    
    sd_sensitivity_gamocc = sd(sensitivity_gamocc),
    sd_specificity_gamocc = sd(specificity_gamocc),
    sd_pp_gamocc = sd(pp_gamocc),          
    sd_np_gamocc = sd(np_gamocc),
    
    sd_sensitivity_gampr = sd(sensitivity_gampr),
    sd_specificity_gampr = sd(specificity_gampr),
    sd_pp_gampr = sd(pp_gampr),         
    sd_np_gampr = sd(np_gampr),
    
    sd_sensitivity_glmocc = sd(sensitivity_glmocc),
    sd_specificity_glmocc = sd(specificity_glmocc),
    sd_pp_glmocc = sd(pp_glmocc),    
    sd_np_glmocc = sd(np_glmocc), 
    
    sd_sensitivity_glmpr = sd(sensitivity_glmpr),
    sd_specificity_glmpr = sd(specificity_glmpr),
    sd_pp_glmpr = sd(pp_glmpr),      
    sd_np_glmpr = sd(np_glmpr),
    
    sd_sensitivity_rfocc = sd(sensitivity_rfocc),
    sd_specificity_rfocc = sd(specificity_rfocc),
    sd_pp_rfocc  = sd(pp_rfocc),        
    sd_np_rfocc = sd(np_rfocc),
    
    sd_sensitivity_rfpr = sd(sensitivity_rfpr),
    sd_specificity_rfpr = sd(specificity_rfpr), 
    sd_pp_rfpr = sd(pp_rfpr),        
    sd_np_rfpr = sd(np_rfpr),
    
    sd_sensitivity_max  = sd(sensitivity_max), 
    sd_specificity_max  = sd(specificity_max), 
    sd_pp_max = sd(pp_max),           
    sd_np_max = sd(np_max))      

pres_matrix_plot <- gather(pres_matrix_means, "Mod", "complement", mean_sensitivity_gamocc:mean_np_max)
pres_matrix_mean <- pivot_longer(pres_matrix_means, mean_sensitivity_gamocc:mean_np_max, "Mod")  %>%
  mutate(complement = value,
         mod = substring(Mod, 6)) %>%
  dplyr::select(n, mod, Mod, complement) 
pres_matrix_sd <- pivot_longer(pres_matrix_means, sd_sensitivity_gamocc:sd_np_max, "Mod", "SD") %>%
  mutate(SD = value, 
         mod = substring(Mod, 4)) %>%
  dplyr::select(n, mod, Mod, SD) %>%
  mutate(se = SD/sqrt(n))

pres_spatial_plot2 <-  full_join(pres_matrix_mean, pres_matrix_sd, by = c("n", "mod")) %>%
  separate(col = Mod.x, into = c("Mean","Measure", "Modtype"), sep = "_")  %>%
  mutate(value = 100 - complement)

pres_spatial_plot2$Modtype <- factor(pres_spatial_plot2$Modtype, levels = c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"), ordered = TRUE)

pres_spatial_plot2$Measure <- factor(pres_spatial_plot2$Measure, levels = c("pp","np","sensitivity","specificity") , ordered = TRUE)

splot = pres_spatial_plot2 %>%
  ggplot(aes(x = Measure, y = value, group = Modtype)) +   
  geom_bar(aes(fill = factor(Modtype)), position='dodge', stat="identity", width=0.75,color = "white", lwd = 3) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), position = position_dodge(0.75), width = 0.2) +
  theme_classic()+ theme(axis.title.x=element_text(size=54),axis.title.y=element_text(size=54, angle=90)) + xlab(bquote("")) + ylab(bquote("Percent")) +
  scale_fill_manual(values = c("#034e7b","navyblue","steelblue2", "dodgerblue2","#238b45", "darkgreen" ,"purple"),
                    breaks=c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"),
                    labels=c("GAM - Occ","GAM - Pr","GLM - Occ","GLM - Pr","RF - Occ","RF - Pr","MaxEnt - Pr")) +
  scale_x_discrete(labels=c("False \ndiscovery rate","False \nomission rate", "False \nnegative rate", "False \npositive rate")) +
  theme(axis.text.x=element_text(size = 50, colour = "black"),axis.ticks=element_blank(), axis.text.y=element_text(size=50, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=50), legend.key.width=unit(2, "lines"), legend.key.size = unit(2, "cm")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) 
# ggsave("Figures/spatial_crossval_25.pdf", width = 30, height = 20)


library(cowplot)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
grid <- plot_grid(pplot + theme(legend.position="top"),
                  splot + theme(legend.position="none"),
                  align = 'v',
                  labels = c("A","B"),
                  label_x = 0.12, 
                  label_size = 30,
                  nrow = 2) 
ggsave("Figures/Figure_5.pdf", height = 20, width = 24)



#### global t tests  ######
# kappa(na.omit(sdm_test$presence), na.omit(sdm_test$predicted_pres), cutoff = 0.7)

t_test <- pres_matrix_plot2 %>%
  mutate(occ_pres = substring(Modtype, 4)) 
t_test$occ_pres <- case_when(t_test$occ_pres == "occ" ~"occ",
                             t_test$occ_pres == "pr" ~ "pr",
                             t_test$occ_pres == "cc" ~ "occ",
                             t_test$occ_pres == "r" ~ "pr",
                             t_test$occ_pres == "" ~ "pr")

t_wide <- t_test %>%
  dplyr::select(Mean, Measure, Modtype, occ_pres, value) %>%
  group_by(Measure, occ_pres) %>%
  pivot_wider(names_from = Measure, values_from = value) %>%
  filter(Modtype != "max") 

 
# false disc rate
t.test(t_wide$pp[t_wide$occ_pres =="occ"], t_wide$pp[t_wide$occ_pres =="pr"], paired = TRUE, alternative= "two.sided")
# false omission rate
t.test(t_wide$np[t_wide$occ_pres =="occ"], t_wide$np[t_wide$occ_pres =="pr"], paired = TRUE, alternative= "two.sided")
# false neg rate
t.test(t_wide$sensitivity[t_wide$occ_pres =="occ"], t_wide$sensitivity[t_wide$occ_pres =="pr"], paired = TRUE, alternative= "two.sided")
# false pos rate
t.test(t_wide$specificity[t_wide$occ_pres =="occ"], t_wide$specificity[t_wide$occ_pres =="pr"], paired = TRUE, alternative= "two.sided")


## spatial 
s_test <- pres_spatial_plot2 %>%
  mutate(occ_pres = substring(Modtype, 4)) 
s_test$occ_pres <- case_when(s_test$occ_pres == "occ" ~"occ",
                             s_test$occ_pres == "pr" ~ "pr",
                             s_test$occ_pres == "cc" ~ "occ",
                             s_test$occ_pres == "r" ~ "pr",
                             s_test$occ_pres == "" ~ "pr")

s_wide <- s_test %>%
  dplyr::select(Mean, Measure, Modtype, occ_pres, value) %>%
  group_by(Measure, occ_pres) %>%
  pivot_wider(names_from = Measure, values_from = value) %>%
  filter(Modtype != "max") 

# false disc rate
t.test(s_wide$pp[s_wide$occ_pres =="occ"], s_wide$pp[s_wide$occ_pres =="pr"], paired = TRUE, alternative= "two.sided")
# false omission rate
t.test(s_wide$np[s_wide$occ_pres =="occ"], s_wide$np[s_wide$occ_pres =="pr"], paired = TRUE, alternative= "two.sided")
# false neg rate
t.test(s_wide$sensitivity[s_wide$occ_pres =="occ"], s_wide$sensitivity[s_wide$occ_pres =="pr"], paired = TRUE, alternative= "two.sided")
# false pos rate
t.test(s_wide$specificity[s_wide$occ_pres =="occ"], s_wide$specificity[s_wide$occ_pres =="pr"], paired = TRUE, alternative= "two.sided")

