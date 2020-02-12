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

##### temporal processing ####
test_df <- read.csv("Data/temporal_crossval_df_25.csv", header = TRUE) 
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
  dplyr::select(-data) %>%
  unnest(cols = c(maxpres)) %>%
  mutate(cats = paste(predicted_max_pres, presence.y, sep = "_"),
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
  summarise(
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
 mean_np_max = mean(np_max))       
 
pres_matrix_plot <- gather(pres_matrix_means, "Mod", "complement", mean_sensitivity_gamocc:mean_np_max)
pres_matrix_plot2 <-  separate(data = pres_matrix_plot, col = Mod, into = c("Mean","Measure", "Modtype"), sep = "_")  %>%
  mutate(value = 100 - complement)
pres_matrix_plot2$Modtype <- factor(pres_matrix_plot2$Modtype, levels = c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"), ordered = TRUE)

pres_matrix_plot2$Measure <- factor(pres_matrix_plot2$Measure, levels = c("pp","np","sensitivity","specificity") , ordered = TRUE)
#+ scale_color_manual(values=c("#034e7b","#034e7b","steelblue2", "steelblue2","#238b45", "#238b45" ,"purple"), labels=c("rmse_gam", "rmse_gam_pres", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres",  "rmse_me_pres"))

pplot = pres_matrix_plot2 %>%
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
# ggsave("Figures/temp_crossval_5.pdf", width = 30, height = 20)

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
# write.csv(sdm_space_cval,"Data/space_cval_75.csv", row.names = FALSE)
rmse <- data.frame(rmse)
names(rmse) <- c("aou","rmse_occ", "rmse_pres", "rmse_gam", "rmse_gam_pres", "rmse_rf", "rmse_rf_pres", "rmse_me_pres")
# write.csv(rmse,"Data/space_cval_rmse.csv", row.names = FALSE)
##### process spatial data #####
sdm_space_cval <- read.csv("Data/space_cval_75.csv", header = TRUE) 

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
  summarise(
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
    mean_np_max = mean(np_max))       

pres_matrix_plot <- gather(pres_matrix_means, "Mod", "complement", mean_sensitivity_gamocc:mean_np_max)
pres_spatial_plot2 <-  separate(data = pres_matrix_plot, col = Mod, into = c("Mean","Measure", "Modtype"), sep = "_")  %>%
  mutate(value = 100 - complement)
pres_spatial_plot2$Modtype <- factor(pres_spatial_plot2$Modtype, levels = c("gamocc","gampr","glmocc","glmpr","rfocc","rfpr","max"), ordered = TRUE)

pres_spatial_plot2$Measure <- factor(pres_spatial_plot2$Measure, levels = c("pp","np","sensitivity","specificity") , ordered = TRUE)

splot = pres_spatial_plot2 %>%
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
ggsave("Figures/xval_plot_75.pdf", height = 20, width = 24)



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





