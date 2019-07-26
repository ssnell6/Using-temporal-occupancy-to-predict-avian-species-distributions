library(tidyverse)
library(glm2)
library(gam)
library(randomForest)
library(dismo)
library(raster)
library(maptools)
library(pROC)
library(purrr)

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
bbs_new_all$pres_2016 <- case_when(is.na(bbs_new_all$presence) == TRUE ~ 0, 
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
num_occ = bbs_inc_absence %>% 
  group_by(Aou) %>% 
  tally(presence) %>% 
  left_join(bbs_inc_absence, ., by = "Aou")

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

sp_list = unique(bbs_final_occ_ll$Aou)
test_df = c()
for(i in sp_list){
  sdm_output = c()
  print(i)
  bbs_sub <- filter(bbs_final_occ_ll, Aou == i) # %>% filter(occ <= 0.33333333) RUN FOR EXCL TRANS
  bbs_new_sub <- filter(bbs_new_all, spAOU == i) 
  temp <- filter(all_env, stateroute %in% bbs_sub$stateroute)
  sdm_input <- left_join(bbs_sub, temp, by = "stateroute")
  sdm_input = na.omit(sdm_input)
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
    if(nrow(filter(sdm_input, presence == 1)) > 49){
      
      gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_input)
      
      pred_gam_occ <- predict(gam_occ,type=c("response"))
      
      thresh <- max(pred_gam_occ) * 0.7 # can test 0.5-0.9
      
      sdm_output = cbind(sdm_input, pred_gam_occ) 
      pred_2016 <- left_join(sdm_output, bbs_new_sub[c("stateroute", "pres_2016")], by = "stateroute") %>%
        mutate(predicted_pres = ifelse(pred_gam_occ > thresh, 1, 0)) %>%
        dplyr::select(Aou, stateroute, occ, presence, latitude, longitude, pred_gam_occ, pres_2016, predicted_pres)
      test_df = rbind(test_df, pred_2016)
    }
  }
}

test_df <- data.frame(test_df)

# write.csv(test_df, "temporal_crossval_df.csv", row.names = FALSE)
test_df <- read.csv("Data/temporal_crossval_df.csv", header = TRUE)
# to account for species not detected in 2015-2016 but are within the range
test_df$pres_2016[is.na(test_df$pres_2016)] = 0 

pres_matrix <- test_df %>% group_by(Aou) %>%
  nest() %>%
  mutate(pres_pres = purrr::map(data, ~{
    data <- .
    table(data$predicted_pres, data$presence)[2,2]}),
    pres_abs = purrr::map(data, ~{
      newdat <- .
      table(newdat$predicted_pres, newdat$presence)[1,2]}), 
    abs_abs = purrr::map(data, ~{
      newdat2 <- .
      table(newdat2$predicted_pres, newdat2$presence)[1,1]}), 
    abs_pres = purrr::map(data, ~{
      newdat3 <- .
      table(newdat3$predicted_pres, newdat3$presence)[2,1]}), 
    length = purrr::map(data, ~{
      newdatlength <- .
      length = length(newdatlength$pred_gam_occ)})) %>%
  dplyr::select(Aou, length, pres_pres, pres_abs, abs_abs, abs_pres) 
pres_matrix <- data.frame(pres_matrix)
pres_matrix$pres_pres <- unlist(pres_matrix$pres_pres)
pres_matrix$abs_pres <- unlist(pres_matrix$abs_pres)
pres_matrix$abs_abs <- unlist(pres_matrix$abs_abs)
pres_matrix$pres_abs <- unlist(pres_matrix$pres_abs)
pres_matrix$length <- unlist(pres_matrix$length)
pres_matrix$true_pos <- pres_matrix$pres_pres/pres_matrix$length
pres_matrix$false_pos <- pres_matrix$pres_abs/pres_matrix$length
pres_matrix$accuracy <- ((pres_matrix$pres_pres + pres_matrix$abs_abs)/pres_matrix$length)*100
pres_matrix$sensitivity <- (pres_matrix$pres_pres/(pres_matrix$pres_pres + pres_matrix$abs_abs))*100
pres_matrix$specificity <- (pres_matrix$abs_abs/(pres_matrix$pres_pres + pres_matrix$abs_abs))*100

kappa(test_df$presence, test_df$predicted_pres, cutoff = 0.7)
t.test(test_df$predicted_pres, test_df$presence, paired = TRUE, alternative= "two.sided")

pplot = ggplot(test_df, aes(x = pres_2016, y = presence)) + geom_point() + theme_classic()+ theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)

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
bbs_final_occ_ll = bbs_final_occ_ll[,c("Aou", "stateroute", "occ", "presence", "ALPHA.CODE",
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
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_202') # for 64-bit version


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
}
write.csv(sdm_output, paste("sdm_output_notrans_", i, ".csv",  sep=""), row.names = FALSE)
}

dev.off()

##### spatial_crossval #######
# code from https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation
sdm_space_cval = c()
for(i in sp_list){
  print(i)
  space_sub <- filter(bbs_final_occ_ll,  Aou == i)
  #Randomly shuffle the data
  space_sub<-space_sub[sample(nrow(space_sub)),]
  #Create 10 equally size folds
  folds <- cut(seq(1,nrow(space_sub)),breaks=10,labels=FALSE)
  
  temp <- filter(all_env, stateroute %in% space_sub$stateroute)
  sdm_input <- left_join(space_sub, temp, by = "stateroute")
  sdm_input = na.omit(sdm_input)
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
    if(nrow(filter(sdm_input, presence == 1)) > 49){
      #Perform 10 fold cross validation
      for(j in 1:10){
        #Segement your data by fold using the which() function 
        testIndexes <- which(folds==j,arr.ind=TRUE)
        testData <- sdm_input[testIndexes, ]
        trainData <- sdm_input[-testIndexes, ]
        gam_train <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = trainData)
        pred_gam_test <- predict(gam_train, testData) 
        thresh <- max(pred_gam_test) * 0.7
        sdm_test <- sdm_input[testIndexes, ] %>%  
          mutate(predicted_pres = ifelse(pred_gam_test > thresh, 1, 0)) %>%
          cbind(., pred_gam_test)
        sdm_space_cval <- rbind(sdm_space_cval, sdm_test)
      }
    }
  }
}

# write.csv(sdm_space_cval,"Data/space_cval.csv", row.names = FALSE)
sdm_space_cval <- read.csv("Z:/hurlbertlab/Snell/space_cval.csv", header = TRUE) 

pres_matrix <- sdm_space_cval %>% group_by(Aou) %>%
  na.omit() %>%
  nest() %>%
  mutate(pres_pres = purrr::map(data, ~{
    data1 <- .
    table(data1$predicted_pres, data1$presence)[2,2]}),
    pres_abs = purrr::map(data, ~{
      newdat <- .
      table(newdat$predicted_pres, newdat$presence)[1,2]}), 
    abs_abs = purrr::map(data, ~{
      newdat2 <- .
      table(newdat2$predicted_pres, newdat2$presence)[1,1]}), 
    abs_pres = purrr::map(data, ~{
      newdat3 <- .
      table(newdat3$predicted_pres, newdat3$presence)[2,1]}), 
    length = purrr::map(data, ~{
      newdatlength <- .
      length = length(newdatlength$pred_gam_test)})) %>%
  dplyr::select(Aou, length, pres_pres, pres_abs, abs_abs, abs_pres) 
pres_matrix <- data.frame(pres_matrix)
pres_matrix$pres_pres <- unlist(pres_matrix$pres_pres)
pres_matrix$abs_pres <- unlist(pres_matrix$abs_pres)
pres_matrix$abs_abs <- unlist(pres_matrix$abs_abs)
pres_matrix$pres_abs <- unlist(pres_matrix$pres_abs)
pres_matrix$length <- unlist(pres_matrix$length)
pres_matrix$true_pos <- pres_matrix$pres_pres/pres_matrix$length
pres_matrix$false_pos <- pres_matrix$pres_abs/pres_matrix$length
pres_matrix$accuracy <- ((pres_matrix$pres_pres + pres_matrix$abs_abs)/pres_matrix$length)*100
pres_matrix$sensitivity <- (pres_matrix$pres_pres/(pres_matrix$pres_pres + pres_matrix$abs_abs))*100
pres_matrix$specificity <- (pres_matrix$abs_abs/(pres_matrix$pres_pres + pres_matrix$abs_abs))*100

kappa(test_df$presence, test_df$predicted_pres, cutoff = 0.7)
t.test(test_df$predicted_pres, test_df$presence, paired = TRUE, alternative= "two.sided")


pplot = ggplot(test_df, aes(x = pres_2016, y = presence)) + geom_point() + theme_classic()+ theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)

t.test(test_df$pres_2016, test_df$presence, paired = TRUE, alternative= "two.sided")


