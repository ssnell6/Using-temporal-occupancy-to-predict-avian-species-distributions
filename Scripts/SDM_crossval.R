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
# thresh <- median(pred_vals)
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
      max_ind_pres = maxent(sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")], sdm_input$presence)
      
      pred_glm_occ <- predict(glm_occ,type=c("response"))
      pred_glm_pr <- predict(glm_pres,type=c("response"))
      pred_gam_occ <- predict(gam_occ,type=c("response"))
      pred_gam_pr <- predict(gam_pres,type=c("response"))
      pred_rf_occ <- predict(rf_occ,type=c("response"))
      pred_rf_pr <- predict(rf_pres,type=c("response"))
      max_pred_pres <- predict(max_ind_pres, sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")])
      
      threshglm_occ <-threshfun(pred_glm_occ)
      threshglm_pr <-threshfun(pred_glm_pr)
      threshgam_occ <-threshfun(pred_gam_occ)
      threshgam_pr <-threshfun(pred_gam_pr)
      threshrf_occ <-threshfun(pred_rf_occ)
      threshrf_pr <-threshfun(pred_rf_pr)
      threshmax_pres <-threshfun(max_pred_pres)
    
      sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_pres) 
      pred_2016 <- left_join(sdm_output, bbs_new_sub[c("stateroute", "presence")], by = "stateroute") %>%
        mutate(predicted_glm_occ = ifelse(pred_glm_occ > threshglm_occ, 1, 0),
               predicted_glm_pr = ifelse(pred_glm_pr > threshglm_pr, 1, 0),
               predicted_gam_occ = ifelse(pred_gam_occ > threshgam_occ, 1, 0),
               predicted_gam_pr = ifelse(pred_gam_pr > threshgam_pr, 1, 0),
               predicted_rf_occ = ifelse(pred_rf_occ > threshrf_occ, 1, 0),
               predicted_rf_pr = ifelse(pred_rf_pr > threshrf_pr, 1, 0),
               predicted_max_pres = ifelse(max_pred_pres > threshmax_pres, 1, 0)) %>%
        dplyr::select(aou, stateroute, occ, presence.x, latitude, longitude, pred_gam_occ, presence.y, predicted_glm_occ, predicted_glm_pr,predicted_gam_occ, predicted_gam_pr, predicted_rf_occ, predicted_rf_pr, predicted_max_pres)
      test_df = rbind(test_df, pred_2016)
    }
  }
}


# write.csv(test_df, "Data/temporal_crossval_df_75.csv", row.names = FALSE) # wrote _5 for thresh of .5, med for median


test_df_25 <- read.csv("Data/temporal_crossval_df_25.csv", header = TRUE) 
test_df_5 <- read.csv("Data/temporal_crossval_df_5.csv", header = TRUE) 
test_df <- read.csv("Data/temporal_crossval_df_75.csv", header = TRUE) 

# data.frame(test_df_75) %>%
#   ggplot(aes(x = pred_gam_occ, group = as.factor(aou))) + geom_density()

  
# to account for species not detected in 2015-2016 but are within the range
test_df$presence.y[is.na(test_df$presence.y)] = 0 

# temp measure to avoid error
#test_df2 <- filter(test_df, !aou %in% c(4090,4900,4950,4860,
      # 7550, 7310, 6120, 3870, 3960,4430,4560,4641,5970,6290,7660))
empty = data.frame(stateroute = 0,occ = 0,presence.x = 0,latitude = 0,longitude = 0,pred_gam_occ = 0,presence.y = 0,predicted_glm_occ = 0,predicted_glm_pr = 0, predicted_gam_occ = 0,predicted_gam_pr = 0,predicted_rf_occ = 0, predicted_rf_pr = 0, predicted_max_pres = 0)

table_rows <- function(dat, col, name, pos1, pos2){
if(nrow(table(dat$col, dat$presence.x)) == 1){
  sub.1 <- bind_rows(dat, empty)
  name = table(sub.1$col, sub.1$presence.x)[pos1, pos2] 
} else{
  name = table(dat$col, dat$presence.x)[pos1, pos2]
  }
}

#### WORK ON THIS 
glmocc <- test_df %>% group_by(aou) %>%
  nest() %>%
  dplyr::mutate(glmocc = purrr::map(data, ~{
      dat <- .
      dat %>%
        group_by(predicted_glm_occ, presence.x) %>%
        count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(glmocc)) %>%
  mutate(cats = paste(predicted_glm_occ, presence.x, sep = "_"),
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
    group_by(predicted_glm_pr, presence.x) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(glmpr)) %>%
  mutate(cats = paste(predicted_glm_pr, presence.x, sep = "_"),
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
    group_by(predicted_gam_occ, presence.x) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(gamocc)) %>%
  mutate(cats = paste(predicted_gam_occ, presence.x, sep = "_"),
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
    group_by(predicted_gam_pr, presence.x) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(gampr)) %>%
  mutate(cats = paste(predicted_gam_pr, presence.x, sep = "_"),
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
    group_by(predicted_rf_occ, presence.x) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(rfocc)) %>%
  mutate(cats = paste(predicted_rf_occ, presence.x, sep = "_"),
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
    group_by(predicted_rf_pr, presence.x) %>%
    count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(rfpres)) %>%
  mutate(cats = paste(predicted_rf_pr, presence.x, sep = "_"),
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
       group_by(predicted_max_pres, presence.x) %>%
       count()})) %>%
  dplyr::select(-data) %>%
  unnest(cols = c(maxpres)) %>%
  mutate(cats = paste(predicted_max_pres, presence.x, sep = "_"),
         cats_cat = case_when(cats == "1_1" ~ "pres_pres_maxpr",
                              cats == "1_0" ~ "pres_abs_maxpr",
                              cats == "0_1" ~ "abs_pres_maxpr",
                              cats == "0_0" ~ "abs_abs_maxpr")) %>%
  dplyr::select(aou, cats_cat, n) %>%
  pivot_wider(names_from = cats_cat, values_from = n)

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
  full_join(length, by = "aou")

pres_matrix <- data.frame(pres_matrix)


pres_matrix$accuracy_gamocc <- ((pres_matrix$pres_pres_gamocc + pres_matrix$abs_abs_gamocc)/pres_matrix$length)*100
pres_matrix$sensitivity_gamocc <- (pres_matrix$pres_pres_gamocc/(pres_matrix$pres_pres_gamocc + pres_matrix$pres_abs_gamocc))*100
pres_matrix$specificity_gamocc <- (pres_matrix$abs_abs_gamocc/(pres_matrix$abs_pres_gamocc + pres_matrix$abs_abs_gamocc))*100
pres_matrix$pp_gamocc <- (pres_matrix$pres_pres_gamocc/(pres_matrix$pres_pres_gamocc + pres_matrix$abs_pres_gamocc))*100
pres_matrix$np_gamocc <- (pres_matrix$abs_abs_gamocc/(pres_matrix$abs_abs_gamocc + pres_matrix$pres_abs_gamocc))*100

pres_matrix$accuracy_gampr <- ((pres_matrix$pres_pres_gampr + pres_matrix$abs_abs_gampr)/pres_matrix$length)*100
pres_matrix$sensitivity_gampr <- (pres_matrix$pres_pres_gampr/(pres_matrix$pres_pres_gampr + pres_matrix$pres_abs_gampr))*100
pres_matrix$specificity_gampr <- (pres_matrix$abs_abs_gampr/(pres_matrix$abs_pres_gampr + pres_matrix$abs_abs_gampr))*100
pres_matrix$pp_gampr <- (pres_matrix$pres_pres_gampr/(pres_matrix$pres_pres_gampr + pres_matrix$abs_pres_gampr))*100
pres_matrix$np_gampr <- (pres_matrix$abs_abs_gampr/(pres_matrix$abs_abs_gampr + pres_matrix$pres_abs_gampr))*100

pres_matrix$accuracy_glmocc <- ((pres_matrix$pres_pres_glmocc + pres_matrix$abs_abs_glmocc)/pres_matrix$length)*100
pres_matrix$sensitivity_glmocc <- (pres_matrix$pres_pres_glmocc/(pres_matrix$pres_pres_glmocc + pres_matrix$pres_abs_glmocc))*100
pres_matrix$specificity_glmocc <- (pres_matrix$abs_abs_glmocc/(pres_matrix$abs_pres_glmocc + pres_matrix$abs_abs_glmocc))*100
pres_matrix$pp_glmocc <- (pres_matrix$pres_pres_glmocc/(pres_matrix$pres_pres_glmocc + pres_matrix$abs_pres_glmocc))*100
pres_matrix$np_glmocc <- (pres_matrix$abs_abs_glmocc/(pres_matrix$abs_abs_glmocc + pres_matrix$pres_abs_glmocc))*100

pres_matrix$accuracy_glmpr <- ((pres_matrix$pres_pres_glmpr + pres_matrix$abs_abs_glmpr)/pres_matrix$length)*100
pres_matrix$sensitivity_glmpr <- (pres_matrix$pres_pres_glmpr/(pres_matrix$pres_pres_glmpr + pres_matrix$pres_abs_glmpr))*100
pres_matrix$specificity_glmpr <- (pres_matrix$abs_abs_glmpr/(pres_matrix$abs_pres_glmpr + pres_matrix$abs_abs_glmpr))*100
pres_matrix$pp_glmpr <- (pres_matrix$pres_pres_glmpr/(pres_matrix$pres_pres_glmpr + pres_matrix$abs_pres_glmpr))*100
pres_matrix$np_glmpr <- (pres_matrix$abs_abs_glmpr/(pres_matrix$abs_abs_glmpr + pres_matrix$pres_abs_glmpr))*100

pres_matrix$accuracy_rfocc <- ((pres_matrix$pres_pres_rfocc + pres_matrix$abs_abs_rfocc)/pres_matrix$length)*100
pres_matrix$sensitivity_rfocc <- (pres_matrix$pres_pres_rfocc/(pres_matrix$pres_pres_rfocc + pres_matrix$pres_abs_rfocc))*100
pres_matrix$specificity_rfocc <- (pres_matrix$abs_abs_rfocc/(pres_matrix$abs_pres_rfocc + pres_matrix$abs_abs_rfocc))*100
pres_matrix$pp_rfocc <- (pres_matrix$pres_pres_rfocc/(pres_matrix$pres_pres_rfocc + pres_matrix$abs_pres_rfocc))*100
pres_matrix$np_rfocc <- (pres_matrix$abs_abs_rfocc/(pres_matrix$abs_abs_rfocc + pres_matrix$pres_abs_rfocc))*100

pres_matrix$accuracy_rfpr <- ((pres_matrix$pres_pres_rfpr + pres_matrix$abs_abs_rfpr)/pres_matrix$length)*100
pres_matrix$sensitivity_rfpr <- (pres_matrix$pres_pres_rfpr/(pres_matrix$pres_pres_rfpr + pres_matrix$pres_abs_rfpr))*100
pres_matrix$specificity_rfpr <- (pres_matrix$abs_abs_rfpr/(pres_matrix$abs_pres_rfpr + pres_matrix$abs_abs_rfpr))*100
pres_matrix$pp_rfpr <- (pres_matrix$pres_pres_rfpr/(pres_matrix$pres_pres_rfpr + pres_matrix$abs_pres_rfpr))*100
pres_matrix$np_rfpr <- (pres_matrix$abs_abs_rfpr/(pres_matrix$abs_abs_rfpr + pres_matrix$pres_abs_rfpr))*100

pres_matrix$accuracy_max <- ((pres_matrix$pres_pres_max + pres_matrix$abs_abs_max)/pres_matrix$length)*100
pres_matrix$sensitivity_max <- (pres_matrix$pres_pres_max/(pres_matrix$pres_pres_max + pres_matrix$pres_abs_max))*100
pres_matrix$specificity_max <- (pres_matrix$abs_abs_max/(pres_matrix$abs_pres_max + pres_matrix$abs_abs_max))*100
pres_matrix$pp_max <- (pres_matrix$pres_pres_max/(pres_matrix$pres_pres_max + pres_matrix$abs_pres_max))*100
pres_matrix$np_max <- (pres_matrix$abs_abs_max/(pres_matrix$abs_abs_max + pres_matrix$pres_abs_max))*100
 
pres_matrix_means <- pres_matrix %>%
  na.omit() %>%
  summarise(
  mean_accuracy_gamocc = mean(accuracy_gamocc),
  mean_sensitivity_gamocc = mean(sensitivity_gamocc),
  mean_specificity_gamocc = mean(specificity_gamocc),
  mean_pp_gamocc = mean(pp_gamocc),          
  mean_np_gamocc = mean(np_gamocc),       
  mean_accuracy_gampr = mean(accuracy_gampr),  
  mean_sensitivity_gampr = mean(sensitivity_gampr),
  mean_specificity_gampr = mean(specificity_gampr),
  mean_pp_gampr = mean(pp_gampr),         
  mean_np_gampr = mean(np_gampr),         
  mean_accuracy_glmocc = mean(accuracy_glmocc),   
 mean_sensitivity_glmocc = mean(sensitivity_glmocc),
 mean_specificity_glmocc = mean(specificity_glmocc),
 mean_pp_glmocc = mean(pp_glmocc),    
 mean_np_glmocc = mean(np_glmocc),     
 mean_accuracy_glmpr = mean(accuracy_glmpr),   
 mean_sensitivity_glmpr = mean(sensitivity_glmpr),
 mean_specificity_glmpr = mean(specificity_glmpr),
 mean_pp_glmpr = mean(pp_glmpr),      
 mean_np_glmpr = mean(np_glmpr),         
 mean_accuracy_rfocc = mean(accuracy_rfocc),  
 mean_sensitivity_rfocc = mean(sensitivity_rfocc),
 mean_specificity_rfocc = mean(specificity_rfocc),
 mean_pp_rfocc  = mean(pp_rfocc),        
 mean_np_rfocc = mean(np_rfocc),        
 mean_accuracy_rfpr = mean(accuracy_rfpr),   
 mean_sensitivity_rfpr = mean(sensitivity_rfpr),
 mean_specificity_rfpr = mean(specificity_rfpr), 
 mean_pp_rfpr = mean(pp_rfpr),        
 mean_np_rfpr = mean(np_rfpr),           
 mean_accuracy_max  = mean(accuracy_max),    
 mean_sensitivity_max  = mean(sensitivity_max), 
 mean_specificity_max  = mean(specificity_max), 
 mean_pp_max = mean(pp_max),           
 mean_np_max = mean(np_max))       
 
pres_matrix_plot <- gather(pres_matrix_means, "Mod", "inverse", mean_accuracy_gamocc:mean_np_max)
pres_matrix_plot2 <-  separate(data = pres_matrix_plot, col = Mod, into = c("Mean","Measure", "Modtype"), sep = "_")   %>% mutate(value = 100-inverse)
pres_matrix_plot2$Modtype <- factor(pres_matrix_plot2$Modtype, levels = c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"), ordered = TRUE)

pres_matrix_plot2$Measure <- factor(pres_matrix_plot2$Measure, levels = c("pp","np","sensitivity","specificity", "accuracy") , ordered = TRUE)
#+ scale_color_manual(values=c("#034e7b","#034e7b","steelblue2", "steelblue2","#238b45", "#238b45" ,"purple"), labels=c("rmse_gam", "rmse_gam_pres", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres",  "rmse_me_pres"))


pplot = filter(pres_matrix_plot2, Measure != "accuracy") %>%
  ggplot(aes(x = Measure, y = value)) +   
  geom_bar(aes(fill = factor(Modtype)), position = "dodge", stat="identity", color = "white", lwd = 3) +
  theme_classic()+ theme(axis.title.x=element_text(size=54),axis.title.y=element_text(size=54, angle=90)) + xlab(bquote("")) + ylab(bquote("Percent")) +
  scale_fill_manual(values = c("#034e7b","navyblue","steelblue2", "dodgerblue2","#238b45", "darkgreen" ,"purple"),
                    breaks=c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"),
                    labels=c("GAM - Occ","GAM - Pr","GLM - Occ","GLM - Pr","RF - Occ","RF - Pr","MaxEnt - Pr")) +
  scale_x_discrete(labels=c("False \ndiscovery rate","False \nomission rate", "False \nnegative rate", "False \npositive rate")) +
  theme(axis.text.x=element_text(size = 50, colour = "black"),axis.ticks=element_blank(), axis.text.y=element_text(size=50, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=50), legend.key.width=unit(2, "lines"), legend.key.size = unit(2, "cm")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) 
ggsave("Figures/temp_crossval_75.pdf", width = 30, height = 20)

####### dummy data #####
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

# read in raw bbs data for 2016
bbs_new <- read.csv("Data/bbs_2015_on.csv", header = TRUE) %>%
  filter(Year == 2016)
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
bbs_final_occ_ll = bbs_final_occ_ll[,c("aou", "stateroute", "occ", "presence", "ALPHA.CODE",
                                       "latitude", "longitude")]
bbs_final_occ_ll$sp_success = 15 * bbs_final_occ_ll$occ
bbs_final_occ_ll$sp_fail = 15 * (1 - bbs_final_occ_ll$occ) 
# write.csv(bbs_occ_code, "Data/final_focal_spp.csv", row.names = FALSE)



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


# test = cor(na.omit(sdm_input_global))
# corrplot(test)
bbs_final_occ_ll <- read.csv("Data/test_val_wide.csv", header = TRUE) %>%
  gather( "Aou","p_o", AOU1:AOU2_occ) 

bbs_final_occ_ll$sp_success = case_when(bbs_final_occ_ll$Aou == "AOU1_occ"|bbs_final_occ_ll$Aou == "AOU2_occ" ~ 15 * bbs_final_occ_ll$p_o)
bbs_final_occ_ll$sp_fail = case_when(bbs_final_occ_ll$Aou == "AOU1_occ"|bbs_final_occ_ll$Aou == "AOU2_occ" ~ 15 * (1 - bbs_final_occ_ll$p_o))

setwd("Data/sdm_dfs/")
pdf('AUC_Curves.pdf', height = 8, width = 10)
par(mfrow = c(2, 3))
auc_df = c()

sp_list = unique(bbs_final_occ_ll$Aou)

for(i in sp_list){
  sdm_output = c()
  
  sdm_input <- filter(bbs_final_occ_ll, Aou == "AOU2"| Aou == "AOU2_occ") # %>% filter(occ <= 0.33333333) RUN FOR EXCL TRANS
  occ = filter(sdm_input, Aou == "AOU2_occ")   
  gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = occ)
  presence = filter(sdm_input, Aou == "AOU2")
  gam_pres <- mgcv::gam(p_o ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = presence)
  
  pred_gam_occ <- predict(gam_occ,type=c("response"))
  pred_gam_pr <- predict(gam_pres,type=c("response"))
  
  
  # sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_pres) 
  #pred_2016 <- left_join(sdm_output, bbs_new_sub[c("stateroute", "pres_2016")], by = "stateroute")
  #thresh <- max(pred_2016$pred_gam_occ) * 0.7
  #gam_rescale <- filter(pred_2016, pred_gam_occ > thresh)
  
  
  rmse_gam <- rmse(as.vector(pred_gam_occ), occ$p_o)
  
  rmse_gam_pres <- rmse(as.vector(pred_gam_pr), occ$p_o)
  
  
  auc_df = rbind(auc_df, c(i, rmse_occ, auc, rmse_pres, auc_pres, rmse_gam, auc_gam, rmse_gam_pres, auc_gam_pres,  rmse_rf, auc_rf, rmse_rf_pres, auc_rf_pres, rmse_me_pres, auc_me_pres, pred_2016))
  j = unique(sdm_input$ALPHA.CODE)
  plot = plot(auc, main = paste("AUC Curve for ", j, ".csv", sep=""))
}
dev.off()

##### spatial_crossval #######
bbs_env <- left_join(bbs_final_occ_ll, all_env, by = "stateroute")
rmse = c()
sdm_space_cval = c()
for(i in sp_list){
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
        max_train = dismo::maxent(trainData[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")], trainData$presence)
        
        pred_gam_occ <- predict(gam_occ_train, testData,type= c("response")) 
        pred_gam_pr <- predict(gam_pres_train,  testData,type= c("response"))
        pred_glm_occ <- predict(glm_occ_train, testData,type= c("response"))
        pred_glm_pr <- predict(glm_pres_train, testData,type= c("response"))
        pred_rf_occ <- predict(rf_occ_train, testData,type= c("response"))
        pred_rf_pr <- predict(rf_pres_train, testData,type= c("response"))
        max_pred_pres <- predict(max_train, as.matrix(testData[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")]), type= c("response"))
        
        threshglm_occ <-threshfun(pred_glm_occ[!is.na(pred_glm_occ)])
        threshglm_pr <-threshfun(pred_glm_pr[!is.na(pred_glm_pr)])
        threshgam_occ <-threshfun(pred_gam_occ[!is.na(pred_gam_occ)])
        threshgam_pr <-threshfun(pred_gam_pr[!is.na(pred_gam_pr)])
        threshrf_occ <-threshfun(pred_rf_occ[!is.na(pred_rf_occ)])
        threshrf_pr <-threshfun(pred_rf_pr[!is.na(pred_rf_pr)])
        threshmax_pres <-threshfun(max_pred_pres[!is.na(max_pred_pres)])
      
        sdm_test <- sdm_input[testIndexes, ] %>%  
          mutate(predicted_glm_occ = ifelse(pred_glm_occ > threshglm_occ, 1, 0),
                 predicted_glm_pr = ifelse(pred_glm_pr > threshglm_pr, 1, 0),
                 predicted_gam_occ = ifelse(pred_gam_occ > threshgam_occ, 1, 0),
                 predicted_gam_pr = ifelse(pred_gam_pr > threshgam_pr, 1, 0),
                 predicted_rf_occ = ifelse(pred_rf_occ > threshrf_occ, 1, 0),
                 predicted_rf_pr = ifelse(pred_rf_pr > threshrf_pr, 1, 0),
                 predicted_max_pres = ifelse(max_pred_pres > threshmax_pres, 1, 0),
                 j = j) 
        rmse_occ <- rmse(pred_glm_occ, sdm_test$occ)
        rmse_pres <- rmse(pred_glm_pr, sdm_test$presence)
        rmse_gam <- rmse(as.vector(pred_gam_occ), sdm_test$occ)
        rmse_gam_pres <- rmse(as.vector(pred_gam_pr), sdm_test$presence)
        rmse_rf <- rmse(pred_rf_occ, sdm_test$occ)
        rmse_rf_pres <- rmse(as.vector(as.numeric(pred_rf_pr)), sdm_test$presence)
        rmse_me_pres <- rmse(max_pred_pres, sdm_test$presence)
        rmse = rbind(rmse, c(i, rmse_occ, rmse_pres, rmse_gam, rmse_gam_pres, rmse_rf, rmse_rf_pres, rmse_me_pres))
        
        sdm_space_cval <- rbind(sdm_space_cval, sdm_test)
      }
    }
  }
}
sdm_space_cval <- data.frame(sdm_space_cval)
# write.csv(sdm_space_cval,"Data/space_cval_5.csv", row.names = FALSE)
sdm_space_cval <- read.csv("Data/space_cval_5.csv", header = TRUE) 
rmse <- data.frame(rmse)
names(rmse) <- c("aou","rmse_occ", "rmse_pres", "rmse_gam", "rmse_gam_pres", "rmse_rf", "rmse_rf_pres", "rmse_me_pres")
# write.csv(rmse,"Data/space_cval_rmse.csv", row.names = FALSE)
# temp measure to avoid error
empty = data.frame(stateroute = 0,occ = 0,presence.x = 0,latitude = 0,longitude = 0,pred_gam_occ = 0,presence.y = 0,predicted_glm_occ = 0,predicted_glm_pr = 0, predicted_gam_occ = 0,predicted_gam_pr = 0,predicted_rf_occ = 0, predicted_rf_pr = 0, predicted_max_pres = 0)

pres_spatial <- sdm_space_cval %>% group_by(aou) %>%
  na.omit() %>%
  dplyr::select(stateroute, aou, occ, presence, latitude, longitude, predicted_glm_occ, predicted_glm_pr, predicted_gam_occ, predicted_gam_pr, predicted_rf_occ, predicted_rf_pr, predicted_max_pres) %>%
  nest() %>%
  dplyr::mutate(pres_pres_glmocc = purrr::map_dbl(data, ~{
    dat <- .
    if(nrow(table(dat$predicted_glm_occ, dat$presence)) == 1){
      sub.1 <- bind_rows(dat, empty)
      pres_pres_glmocc = table(sub.1$predicted_glm_occ, sub.1$presence)[2,2] 
    } else{
      pres_pres_glmocc = table(dat$predicted_glm_occ, dat$presence)[2,2]
    }}),
    pres_abs_glmocc = purrr::map_dbl(data, ~{
      newdat <- .
      table(newdat$predicted_glm_occ, newdat$presence)[1,2]}),
    abs_abs_glmocc = purrr::map_dbl(data, ~{
      newdat2 <- .
      table(newdat2$predicted_glm_occ, newdat2$presence)[1,1]}),
    abs_pres_glmocc = purrr::map_dbl(data, ~{
      newdat3 <- .
      table(newdat3$predicted_glm_occ, newdat3$presence)[2,1]}),
    
    pres_pres_glmpr = purrr::map_dbl(data, ~{
      newdat4 <- .
      if(nrow(table(newdat4$predicted_glm_pr, newdat4$presence)) == 1){
        sub.1 <- bind_rows(newdat4, empty)
        pres_pres_glmpr = table(sub.1$predicted_glm_pr, sub.1$presence)[2,2] 
      } else{
        pres_pres_glmpr = table(newdat4$predicted_glm_pr, newdat4$presence)[2,2]
      }}),
    pres_abs_glmpr = purrr::map_dbl(data, ~{
      newdat5 <- .
      table(newdat5$predicted_glm_pr, newdat5$presence)[1,2]}),
    abs_abs_glmpr = purrr::map_dbl(data, ~{
      newdat6 <- .
      table(newdat6$predicted_glm_pr, newdat6$presence)[1,1]}),
    abs_pres_glmpr = purrr::map_dbl(data, ~{
      newdat7 <- .
      if(nrow(table(newdat7$predicted_glm_pr, newdat7$presence)) == 1){
        sub.1 <- bind_rows(newdat7, empty)
        abs_pres_glmpr = table(sub.1$predicted_glm_pr, sub.1$presence)[2,1] 
      } else{
        abs_pres_glmpr = table(newdat7$predicted_glm_pr, newdat7$presence)[2,1]
      }}),
    
    pres_pres_gamocc = purrr::map_dbl(data, ~{
      newdat8 <- .
      table(newdat8$predicted_gam_occ, newdat8$presence)[2,2]}),
    pres_abs_gamocc = purrr::map_dbl(data, ~{
      newdat9 <- .
      table(newdat9$predicted_gam_occ, newdat9$presence)[1,2]}),
    abs_abs_gamocc = purrr::map_dbl(data, ~{
      newdat10 <- .
      table(newdat10$predicted_gam_occ, newdat10$presence)[1,1]}),
    abs_pres_gamocc = purrr::map_dbl(data, ~{
      newdat11 <- .
      table(newdat11$predicted_gam_occ, newdat11$presence)[2,1]}),
    
    pres_pres_gampr = purrr::map_dbl(data, ~{
      newdat12 <- .
      if(nrow(table(newdat12$predicted_gam_pr, newdat12$presence)) == 1){
        sub.1 <- bind_rows(newdat12, empty)
        abs_pres_glmpr = table(sub.1$predicted_gam_pr, sub.1$presence)[2,2] 
      } else{
        abs_pres_glmpr = table(newdat12$predicted_gam_pr, newdat12$presence)[2,2]
      }}),
    pres_abs_gampr = purrr::map_dbl(data, ~{
      newdat13 <- .
      table(newdat13$predicted_gam_pr, newdat13$presence)[1,2]}),
    abs_abs_gampr = purrr::map_dbl(data, ~{
      newdat14 <- .
      table(newdat14$predicted_gam_pr, newdat14$presence)[1,1]}),
    abs_pres_gampr = purrr::map_dbl(data, ~{
      newdat15 <- .
        if(nrow(table(newdat15$predicted_gam_pr, newdat15$presence)) == 1){
          sub.1 <- bind_rows(newdat15, empty)
          abs_pres_glmpr = table(sub.1$predicted_gam_pr, sub.1$presence)[2,1] 
        } else{
          abs_pres_glmpr = table(newdat15$predicted_gam_pr, newdat15$presence)[2,1]
        }}),
    
    pres_pres_rfocc = purrr::map_dbl(data, ~{
      newdat16 <- .
      table(newdat16$predicted_rf_occ, newdat16$presence)[2,2]}),
    pres_abs_rfocc = purrr::map_dbl(data, ~{
      newdat17 <- .
      table(newdat17$predicted_rf_occ, newdat17$presence)[1,2]}),
    abs_abs_rfocc = purrr::map_dbl(data, ~{
      newdat18 <- .
      table(newdat18$predicted_rf_occ, newdat18$presence)[1,1]}),
    abs_pres_rfocc = purrr::map_dbl(data, ~{
      newdat19 <- .
      table(newdat19$predicted_rf_occ, newdat19$presence)[2,1]}),
    
    pres_pres_rfpres = purrr::map_dbl(data, ~{
      newdat20 <- .
      if(nrow(table(newdat20$predicted_rf_pr, newdat20$presence)) == 1){
        sub.1 <- bind_rows(newdat20, empty)
        pres_pres_rfpres = table(sub.1$predicted_rf_pr, sub.1$presence)[2,2] 
      } else{
        pres_pres_rfpres = table(newdat20$predicted_rf_pr, newdat20$presence)[2,2]
      }}),
    pres_abs_rfpres = purrr::map_dbl(data, ~{
      newdat21 <- .
      table(newdat21$predicted_rf_pr, newdat21$presence)[1,2]}),
    abs_abs_rfpres = purrr::map_dbl(data, ~{
      newdat22 <- .
      table(newdat22$predicted_rf_pr, newdat22$presence)[1,1]}),
    abs_pres_rfpres = purrr::map_dbl(data, ~{
      newdat23 <- .
      if(nrow(table(newdat23$predicted_rf_pr, newdat23$presence)) == 1){
        sub.1 <- bind_rows(newdat23, empty)
        abs_pres_rfpres = table(sub.1$predicted_rf_pr, sub.1$presence)[2,1] 
      } else{
        abs_pres_rfpres = table(newdat23$predicted_rf_pr, newdat23$presence)[2,1]
      }}),
    pres_pres_max = purrr::map_dbl(data, ~{
      newdat24 <- .
      if(nrow(table(newdat24$predicted_max_pres, newdat24$presence)) == 1){
        sub.1 <- bind_rows(newdat24, empty)
        pres_pres_max = table(sub.1$predicted_max_pres, sub.1$presence)[2,2]
      } else{
        pres_pres_max = table(newdat24$predicted_max_pres, newdat24$presence)[2,2]
      }}),
    pres_abs_max = purrr::map_dbl(data, ~{
      newdat25 <- .
      table(newdat25$predicted_max_pres, newdat25$presence)[1,2]}),
    abs_abs_max = purrr::map_dbl(data, ~{
      newdat26 <- .
      table(newdat26$predicted_max_pres, newdat26$presence)[1,1]}),
    abs_pres_max = purrr::map_dbl(data, ~{
      newdat27 <- .
      if(nrow(table(newdat27$predicted_max_pres, newdat27$presence)) == 1){
        sub.1 <- bind_rows(newdat27, empty)
        abs_pres_max = table(sub.1$predicted_max_pres, sub.1$presence)[2,1]
      } else{
        abs_pres_max = table(newdat27$predicted_max_pres, newdat27$presence)[2,1]
      }}),
    length = purrr::map_dbl(data, ~{
      newdatlength <- .
      length(newdatlength$stateroute)})) %>%
  dplyr::select(-data) 


pres_spatial$accuracy_gamocc <- ((pres_spatial$pres_pres_gamocc + pres_spatial$abs_abs_gamocc)/pres_spatial$length)*100
pres_spatial$sensitivity_gamocc <- (pres_spatial$pres_pres_gamocc/(pres_spatial$pres_pres_gamocc + pres_spatial$pres_abs_gamocc))*100
pres_spatial$specificity_gamocc <- (pres_spatial$abs_abs_gamocc/(pres_spatial$abs_pres_gamocc + pres_spatial$abs_abs_gamocc))*100
pres_spatial$pp_gamocc <- (pres_spatial$pres_pres_gamocc/(pres_spatial$pres_pres_gamocc + pres_spatial$abs_pres_gamocc))*100
pres_spatial$np_gamocc <- (pres_spatial$abs_abs_gamocc/(pres_spatial$abs_abs_gamocc + pres_spatial$pres_abs_gamocc))*100

pres_spatial$accuracy_gampr <- ((pres_spatial$pres_pres_gampr + pres_spatial$abs_abs_gampr)/pres_spatial$length)*100
pres_spatial$sensitivity_gampr <- (pres_spatial$pres_pres_gampr/(pres_spatial$pres_pres_gampr + pres_spatial$pres_abs_gampr))*100
pres_spatial$specificity_gampr <- (pres_spatial$abs_abs_gampr/(pres_spatial$abs_pres_gampr + pres_spatial$abs_abs_gampr))*100
pres_spatial$pp_gampr <- (pres_spatial$pres_pres_gampr/(pres_spatial$pres_pres_gampr + pres_spatial$abs_pres_gampr))*100
pres_spatial$np_gampr <- (pres_spatial$abs_abs_gampr/(pres_spatial$abs_abs_gampr + pres_spatial$pres_abs_gampr))*100

pres_spatial$accuracy_glmocc <- ((pres_spatial$pres_pres_glmocc + pres_spatial$abs_abs_glmocc)/pres_spatial$length)*100
pres_spatial$sensitivity_glmocc <- (pres_spatial$pres_pres_glmocc/(pres_spatial$pres_pres_glmocc + pres_spatial$pres_abs_glmocc))*100
pres_spatial$specificity_glmocc <- (pres_spatial$abs_abs_glmocc/(pres_spatial$abs_pres_glmocc + pres_spatial$abs_abs_glmocc))*100
pres_spatial$pp_glmocc <- (pres_spatial$pres_pres_glmocc/(pres_spatial$pres_pres_glmocc + pres_spatial$abs_pres_glmocc))*100
pres_spatial$np_glmocc <- (pres_spatial$abs_abs_glmocc/(pres_spatial$abs_abs_glmocc + pres_spatial$pres_abs_glmocc))*100

pres_spatial$accuracy_glmpr <- ((pres_spatial$pres_pres_glmpr + pres_spatial$abs_abs_glmpr)/pres_spatial$length)*100
pres_spatial$sensitivity_glmpr <- (pres_spatial$pres_pres_glmpr/(pres_spatial$pres_pres_glmpr + pres_spatial$pres_abs_glmpr))*100
pres_spatial$specificity_glmpr <- (pres_spatial$abs_abs_glmpr/(pres_spatial$abs_pres_glmpr + pres_spatial$abs_abs_glmpr))*100
pres_spatial$pp_glmpr <- (pres_spatial$pres_pres_glmpr/(pres_spatial$pres_pres_glmpr + pres_spatial$abs_pres_glmpr))*100
pres_spatial$np_glmpr <- (pres_spatial$abs_abs_glmpr/(pres_spatial$abs_abs_glmpr + pres_spatial$pres_abs_glmpr))*100

pres_spatial$accuracy_rfocc <- ((pres_spatial$pres_pres_rfocc + pres_spatial$abs_abs_rfocc)/pres_spatial$length)*100
pres_spatial$sensitivity_rfocc <- (pres_spatial$pres_pres_rfocc/(pres_spatial$pres_pres_rfocc + pres_spatial$pres_abs_rfocc))*100
pres_spatial$specificity_rfocc <- (pres_spatial$abs_abs_rfocc/(pres_spatial$abs_pres_rfocc + pres_spatial$abs_abs_rfocc))*100
pres_spatial$pp_rfocc <- (pres_spatial$pres_pres_rfocc/(pres_spatial$pres_pres_rfocc + pres_spatial$abs_pres_rfocc))*100
pres_spatial$np_rfocc <- (pres_spatial$abs_abs_rfocc/(pres_spatial$abs_abs_rfocc + pres_spatial$pres_abs_rfocc))*100

pres_spatial$accuracy_rfpr <- ((pres_spatial$pres_pres_rfpres + pres_spatial$abs_abs_rfpres)/pres_spatial$length)*100
pres_spatial$sensitivity_rfpr <- (pres_spatial$pres_pres_rfpres/(pres_spatial$pres_pres_rfpres + pres_spatial$pres_abs_rfpres))*100
pres_spatial$specificity_rfpr <- (pres_spatial$abs_abs_rfpres/(pres_spatial$abs_pres_rfpres + pres_spatial$abs_abs_rfpres))*100
pres_spatial$pp_rfpr <- (pres_spatial$pres_pres_rfpres/(pres_spatial$pres_pres_rfpres + pres_spatial$abs_pres_rfpres))*100
pres_spatial$np_rfpr <- (pres_spatial$abs_abs_rfpres/(pres_spatial$abs_abs_rfpres + pres_spatial$pres_abs_rfpres))*100

pres_spatial$accuracy_max <- ((pres_spatial$pres_pres_max + pres_spatial$abs_abs_max)/pres_spatial$length)*100
pres_spatial$sensitivity_max <- (pres_spatial$pres_pres_max/(pres_spatial$pres_pres_max + pres_spatial$pres_abs_max))*100
pres_spatial$specificity_max <- (pres_spatial$abs_abs_max/(pres_spatial$abs_pres_max + pres_spatial$abs_abs_max))*100
pres_spatial$pp_max <- (pres_spatial$pres_pres_max/(pres_spatial$pres_pres_max + pres_spatial$abs_pres_max))*100
pres_spatial$np_max <- (pres_spatial$abs_abs_max/(pres_spatial$abs_abs_max + pres_spatial$pres_abs_max))*100

# kappa(test_df$presence, test_df$predicted_pres, cutoff = 0.7)
# t.test(test_df$predicted_pres, test_df$presence, paired = TRUE, alternative= "two.sided")
pres_spatial_plot <- gather(pres_spatial, "Mod", "value", accuracy_gamocc:np_max)
pres_spatial_plot2 <-  separate(data = pres_spatial_plot, col = Mod, into = c("Measure", "Modtype"), sep = "_")
pres_spatial_plot2$Modtype <- factor(pres_spatial_plot2$Modtype, levels = c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max", ordered = TRUE))

pres_spatial <- data.frame(pres_spatial)
pres_spatial_means <- pres_spatial %>%
  summarise(
    mean_accuracy_gamocc = mean(accuracy_gamocc),
    mean_sensitivity_gamocc = mean(sensitivity_gamocc),
    mean_specificity_gamocc = mean(specificity_gamocc),
    mean_pp_gamocc = mean(pp_gamocc),          
    mean_np_gamocc = mean(np_gamocc),       
    mean_accuracy_gampr = mean(accuracy_gampr),  
    mean_sensitivity_gampr = mean(sensitivity_gampr),
    mean_specificity_gampr = mean(specificity_gampr),
    mean_pp_gampr = mean(pp_gampr),         
    mean_np_gampr = mean(np_gampr),         
    mean_accuracy_glmocc = mean(accuracy_glmocc),   
    mean_sensitivity_glmocc = mean(sensitivity_glmocc),
    mean_specificity_glmocc = mean(specificity_glmocc),
    mean_pp_glmocc = mean(pp_glmocc),    
    mean_np_glmocc = mean(np_glmocc),     
    mean_accuracy_glmpr = mean(accuracy_glmpr),   
    mean_sensitivity_glmpr = mean(sensitivity_glmpr),
    mean_specificity_glmpr = mean(specificity_glmpr),
    mean_pp_glmpr = mean(pp_glmpr),      
    mean_np_glmpr = mean(np_glmpr),         
    mean_accuracy_rfocc = mean(accuracy_rfocc),  
    mean_sensitivity_rfocc = mean(sensitivity_rfocc),
    mean_specificity_rfocc = mean(specificity_rfocc),
    mean_pp_rfocc  = mean(pp_rfocc),        
    mean_np_rfocc = mean(np_rfocc),        
    mean_accuracy_rfpr = mean(accuracy_rfpr),   
    mean_sensitivity_rfpr = mean(sensitivity_rfpr),
    mean_specificity_rfpr = mean(specificity_rfpr), 
    mean_pp_rfpr = mean(pp_rfpr),        
    mean_np_rfpr = mean(np_rfpr),           
    mean_accuracy_max  = mean(accuracy_max),    
    mean_sensitivity_max  = mean(sensitivity_max), 
    mean_specificity_max  = mean(specificity_max), 
    mean_pp_max = mean(pp_max),           
    mean_np_max = mean(np_max))       

pres_spatial_plot <- gather(pres_spatial_means, "Mod", "inverse", mean_accuracy_gamocc:mean_np_max)
pres_spatial_plot2 <-  separate(data = pres_spatial_plot, col = Mod, into = c("Mean","Measure", "Modtype"), sep = "_")  %>%
  mutate(value = 100 - inverse)
pres_spatial_plot2$Modtype <- factor(pres_spatial_plot2$Modtype, levels = c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max", ordered = TRUE))
pres_spatial_plot2$Measure <- factor(pres_spatial_plot2$Measure, levels = c("pp","np","sensitivity","specificity", "accuracy") , ordered = TRUE)

splot = filter(pres_spatial_plot2, Measure != "accuracy") %>%
  ggplot(aes(x = Measure, y = value)) +   
  geom_bar(aes(fill = factor(Modtype)), position = "dodge", stat="identity",  color = "white", lwd = 3) +
  theme_classic()+ theme(axis.title.x=element_text(size=54),axis.title.y=element_text(size=54, angle=90)) + xlab(bquote("")) + ylab(bquote("Percent")) +
  scale_fill_manual(values = c("#034e7b","navyblue","steelblue2", "dodgerblue2","#238b45", "darkgreen" ,"purple"),
                    breaks=c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"),
                    labels=c("GAM - Occ","GAM - Pr","GLM - Occ","GLM - Pr","RF - Occ","RF - Pr","MaxEnt - Pr")) +
  scale_x_discrete(labels=c("False \ndiscovery rate","False \nomission rate", "False \nnegative rate", "False \npositive rate")) +
  theme(axis.text.x=element_text(size = 50, colour = "black"),axis.ticks=element_blank(), axis.text.y=element_text(size=50, colour = "black")) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=50), legend.key.width=unit(2, "lines"), legend.key.size = unit(2, "cm")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) 
ggsave("Figures/spatial_crossval.pdf", width = 30, height = 20)


library(cowplot)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
grid <- plot_grid(pplot + theme(legend.position="top"),
                  splot + theme(legend.position="none"),
                  align = 'v',
                  labels = c("A","B"),
                  label_x = 0.12, 
                  label_size = 30,
                  nrow = 2) 
ggsave("Figures/xval_plot.pdf", height = 20, width = 24)

glm_acc <- ggplot(pres_spatial, aes(x = accuracy_glmocc/100, y = accuracy_glmpr/100)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Spatial Occ Accuracy")) + ylab(bquote("Spatial Pres Accuracy"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(shape=16, aes(size = pres_spatial$length)) + 
  scale_y_continuous(limit = c(0, 1)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) +
    theme(legend.title=element_blank(), legend.text=element_text(size=15)) +
  annotate("text", x = 0.8, y = 0.01, label = "GLM", size = 10)

gam_acc <- ggplot(pres_spatial, aes(x = accuracy_gamocc/100, y = accuracy_gampr/100)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Spatial Occ Accuracy")) + ylab(bquote("Spatial Pres Accuracy"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(shape=16, aes(size = pres_spatial$length)) + 
  scale_y_continuous(limit = c(0, 1)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15)) +
  annotate("text", x = 0.8, y = 0.01, label = "GAM", size = 10)

rf_acc <- ggplot(pres_spatial, aes(x = accuracy_rfocc/100, y = accuracy_rfpr/100)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Spatial Occ Accuracy")) + ylab(bquote("Spatial Pres Accuracy"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(shape=16, aes(size = pres_spatial$length)) + 
  scale_y_continuous(limit = c(0, 1)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15)) +
  annotate("text", x = 0.8, y = 0.01, label = "RF", size = 10)
#ggsave("Figures/acc_occ_rf.pdf", height = 10, width = 24)

library(cowplot)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
plot_grid(glm_acc + theme(legend.position="none"),
          gam_acc + theme(legend.position="none"),
          rf_acc + theme(legend.position="none"),
          align = 'hv',
          labels = c("A","B", "C"),
          label_x = 0.22, 
          label_size = 30,
          nrow = 1, 
          scale = 0.9) 
ggsave("Figures/acc_occ_pres.pdf", height = 10, width = 24)





glm_acc <- ggplot(pres_matrix, aes(x = accuracy_glmocc/100, y = accuracy_glmpr/100)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Temporal Occ Accuracy")) + ylab(bquote("Temporal Pres Accuracy"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(shape=16, aes(size = pres_matrix$length)) + 
  scale_y_continuous(limit = c(0, 1)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15)) +
  annotate("text", x = 0.8, y = 0.01, label = "GLM", size = 10)

gam_acc <- ggplot(pres_matrix, aes(x = accuracy_gamocc/100, y = accuracy_gampr/100)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Temporal Occ Accuracy")) + ylab(bquote("Temporal Pres Accuracy"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(shape=16, aes(size = pres_matrix$length)) + 
  scale_y_continuous(limit = c(0, 1)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15))  +
  annotate("text", x = 0.8, y = 0.01, label = "GAM", size = 10)

rf_acc <- ggplot(pres_matrix, aes(x = accuracy_rfocc/100, y = accuracy_rfpr/100)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Temporal Occ Accuracy")) + ylab(bquote("Temporal Pres Accuracy"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(shape=16, aes(size = pres_matrix$length)) + 
  scale_y_continuous(limit = c(0, 1)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15))  +
  annotate("text", x = 0.8, y = 0.01, label = "RF", size = 10)
# ggsave("Figures/acc_occ_rf.pdf", height = 10, width = 24)

library(cowplot)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
plot_grid(glm_acc + theme(legend.position="none"),
          gam_acc + theme(legend.position="none"),
          rf_acc + theme(legend.position="none"),
          align = 'hv',
          labels = c("A","B", "C"),
          label_x = 0.22, 
          label_size = 30,
          nrow = 1, 
          scale = 0.9) 
ggsave("Figures/acc_occ_pres_temp.pdf", height = 10, width = 24)

pres_spatial$diff <- pres_spatial$accuracy_glmocc - pres_spatial$accuracy_glmpr
sub2 <- filter(pres_spatial, diff < 0) %>%
  left_join(bbs_final_occ_ll, by = c("aou")) %>%
  select(aou, stateroute) %>%
  group_by(aou) %>%
  summarise(n = n_distinct(stateroute))
hist <- filter(sub2, aou == 6270)
hist(hist$occ)


#### test  ######
test<- c()
for(i in na.omit(unique(sdm_test$Aou))){
  pres_matrix <- sdm_test %>% group_by(Aou) %>%
    filter(Aou == i) %>%
    na.omit() 
  if(sum(pres_matrix$predicted_pres) > 0){
    pres_pres = table(pres_matrix$predicted_pres, pres_matrix$presence)[2,2]
      pres_abs = table(pres_matrix$predicted_pres, pres_matrix$presence)[1,2]
      abs_abs = table(pres_matrix$predicted_pres, pres_matrix$presence)[1,1]
      abs_pres = table(pres_matrix$predicted_pres, pres_matrix$presence)[2,1] 
      length = length(pres_matrix$pred_gam_test)
  test = rbind(test, c(i, length, pres_abs, pres_pres, abs_abs, abs_pres))
  }
}
pres_matrix <- data.frame(test)
names(pres_matrix) <- c("Aou", "length", "pres_abs", "pres_pres", "abs_abs", "abs_pres")
pres_matrix$pres_pres <- unlist(pres_matrix$pres_pres)
pres_matrix$abs_pres <- unlist(pres_matrix$abs_pres)
pres_matrix$abs_abs <- unlist(pres_matrix$abs_abs)
pres_matrix$pres_abs <- unlist(pres_matrix$pres_abs)
pres_matrix$length <- unlist(pres_matrix$length)
pres_matrix$truepos <- pres_matrix$pres_pres/pres_matrix$length
pres_matrix$falsepos <- pres_matrix$pres_abs/pres_matrix$length
pres_matrix$accuracy <- ((pres_matrix$pres_pres + pres_matrix$abs_abs)/pres_matrix$length)*100
pres_matrix$sensitivity <- (pres_matrix$pres_pres/(pres_matrix$pres_pres + pres_matrix$abs_abs))*100
pres_matrix$specificity <- (pres_matrix$abs_abs/(pres_matrix$pres_pres + pres_matrix$abs_abs))*100
# 23 species had 0 predicted presences
  
kappa(na.omit(sdm_test$presence), na.omit(sdm_test$predicted_pres), cutoff = 0.7)
t.test(sdm_test$predicted_pres, sdm_test$presence, paired = TRUE, alternative= "two.sided")
table(sdm_space_cval$presence, sdm_space_cval$predicted_pres)


