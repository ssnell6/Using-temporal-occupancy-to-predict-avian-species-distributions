library(tidyverse)
library(glm2)
library(gam)
library(randomForest)
library(dismo)
library(raster)
library(maptools)
library(pROC)
library(hydroGOF)

bbs_occ = read.csv("Data/bbs_2001_2015.csv", header=TRUE) 

bbs_occ_sub = bbs_occ %>% 
  dplyr::count(aou, stateroute) %>% 
  filter(n < 6) %>%
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
  select(aou, stateroute, occ)
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


###### SDM analysis ######
auc_df_5 = c()
sp_list = unique(bbs_final_occ_ll$aou)

for(i in sp_list){
  sdm_output = c()
  print(i)
  bbs_sub <- filter(bbs_final_occ_ll, aou == i) 
  sdm_input <- filter(all_env, stateroute %in% bbs_sub$stateroute) %>%
    full_join(bbs_sub, by = "stateroute")
  sdm_input <- data.frame(sdm_input)
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
  if(nrow(filter(sdm_input, presence == 1)) > 59){
    glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
    glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
    gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_input)
    gam_pres <- mgcv::gam(presence ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = sdm_input)
    rf_occ <- randomForest(sp_success/5 ~elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
    rf_pres <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
   
    max_ind_pres = maxent(sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")], sdm_input$presence)
    
    # predict
    pred_glm_occ <- predict(glm_occ,type=c("response"))
    pred_glm_pr <- predict(glm_pres,type=c("response"))
    # validation across time
    # need to get predict fct working, currently super low. Need to tie to stateroute somehow. Need to set threshold.
    # pred_glm_occ_val <- predict(sdm_input$presence, glm_occ)
      # pred_glm_occ[pred_glm_occ > 0.3]
 
    pred_gam_occ <- predict(gam_occ,type=c("response"))
    pred_gam_pr <- predict(gam_pres,type=c("response"))
    pred_rf_occ <- predict(rf_occ,type=c("response"))
    pred_rf_pr <- predict(rf_pres,type=c("response"))
    max_pred_pres <- predict(max_ind_pres, sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")])
    
    sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_pres) 
    
 rmse_occ <- rmse(sdm_output$pred_glm_occ, sdm_output$occ)
 rmse_pres <- rmse(sdm_output$pred_glm_pr, sdm_output$presence)
 rmse_gam <- rmse(as.vector(sdm_output$pred_gam_occ), sdm_output$occ)
 rmse_gam_pres <- rmse(as.vector(sdm_output$pred_gam_pr), sdm_output$presence)
 rmse_rf <- rmse(sdm_output$pred_rf_occ, sdm_output$occ)
 rmse_rf_pres <- rmse(as.vector(as.numeric(sdm_output$pred_rf_pr)), sdm_output$presence)
 rmse_me_pres <- rmse(sdm_output$max_pred_pres, sdm_output$presence)
 auc_df_5 = rbind(auc_df_5, c(i, rmse_occ, rmse_pres, rmse_gam, rmse_gam_pres, rmse_rf, rmse_rf_pres, rmse_me_pres))
 j = unique(sdm_input$ALPHA.CODE)
  }
  }
}

auc_df_5 = data.frame(auc_df_5)
names(auc_df_5) = c("AOU", "rmse_occ", "rmse_pres","rmse_gam", "rmse_gam_pres", "rmse_rf", "rmse_rf_pres","rmse_me_pres")
# write.csv(auc_df_5, "Data/auc_df_5.csv", row.names = FALSE)


##### no trans #####
auc_df_notrans_5 = c()
sp_list = unique(bbs_final_occ_ll$aou)

for(i in sp_list){
  sdm_output = c()
  print(i)
  bbs_sub <- filter(bbs_final_occ_ll, aou == i) %>% 
    mutate(excl_pres = ifelse(occ <= 0.33, 0, 1))
  sdm_input <- filter(all_env, stateroute %in% bbs_sub$stateroute) %>%
    full_join(bbs_sub, by = "stateroute")
  sdm_input <- data.frame(sdm_input)
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
    if(nrow(filter(sdm_input, presence == 1)) > 59){
      glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      glm_pres <- glm(excl_pres ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_input)
      gam_pres <- mgcv::gam(excl_pres ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = sdm_input)
      rf_occ <- randomForest(sp_success/5 ~elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      rf_pres <- randomForest(excl_pres ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
      
      max_ind_pres = maxent(sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")], sdm_input$excl_pres)
      
      # predict
      pred_glm_occ <- predict(glm_occ,type=c("response"))
      pred_glm_pr <- predict(glm_pres,type=c("response"))
      pred_gam_occ <- predict(gam_occ,type=c("response"))
      pred_gam_pr <- predict(gam_pres,type=c("response"))
      pred_rf_occ <- predict(rf_occ,type=c("response"))
      pred_rf_pr <- predict(rf_pres,type=c("response"))
      max_pred_pres <- predict(max_ind_pres, sdm_input[,c("elev.mean", "bio.mean.bio4","bio.mean.bio5","bio.mean.bio6","bio.mean.bio13","bio.mean.bio14", "ndvi.mean")])
      
      sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_pres) 

      rmse_occ_notrans <- rmse(sdm_output$pred_glm_occ, sdm_output$occ)
      rmse_pres_notrans <- rmse(sdm_output$pred_glm_pr, sdm_output$excl_pres)
      rmse_gam_notrans <- rmse(as.vector(sdm_output$pred_gam_occ), sdm_output$occ)
      rmse_gam_pres_notrans <- rmse(as.vector(sdm_output$pred_gam_pr), sdm_output$excl_pres)
      rmse_rf_notrans <- rmse(sdm_output$pred_rf_occ, sdm_output$occ)
      rmse_rf_pres_notrans <- rmse(as.vector(as.numeric(sdm_output$pred_rf_pr)), sdm_output$excl_pres)
      rmse_me_pres_notrans <- rmse(sdm_output$max_pred_pres, sdm_output$excl_pres)
      auc_df_notrans_5 = rbind(auc_df_notrans_5, c(i, rmse_occ_notrans, rmse_pres_notrans, rmse_gam_notrans, rmse_gam_pres_notrans, rmse_rf_notrans, rmse_rf_pres_notrans, rmse_me_pres_notrans))
      j = unique(sdm_input$ALPHA.CODE)
    }
  }
  # write.csv(sdm_output, paste("sdm_output_notrans_", i, ".csv",  sep=""), row.names = FALSE)
}

auc_df_notrans_5 = data.frame(auc_df_notrans_5)
names(auc_df_notrans_5) = c("AOU", "rmse_occ_notrans", "rmse_pres_notrans","rmse_gam_notrans", "rmse_gam_pres_notrans", "rmse_rf_notrans", "rmse_rf_pres_notrans","rmse_me_pres_notrans")
# write.csv(auc_df_notrans_5, "Data/auc_df_notrans_5.csv", row.names = FALSE)

bbs_final_occ_ll$presence <- factor(bbs_final_occ_ll$presence,levels = c('1','0'), ordered = TRUE)

###### global plots ####
auc_df <- read.csv("Data/auc_df_5.csv", header = TRUE)
auc_df_notrans <- read.csv("Data/auc_df_notrans_5.csv", header = TRUE)

num_pres = bbs_final_occ_ll %>%
  group_by(aou) %>% 
  dplyr::filter(., presence == "1") %>%
  dplyr::summarise(n_pres = n_distinct(stateroute))  

num_routes = bbs_final_occ_ll %>% group_by(aou) %>% 
  dplyr::summarise(n = n_distinct(stateroute))  %>%
  dplyr::filter(., n >= 40)

auc_df_merge = left_join(auc_df, num_pres, by = c("AOU" = "aou"))
auc_df_traits = left_join(auc_df_merge, traits, by = "AOU") %>%
  left_join(., tax_code, by = c("AOU" = "AOU_OUT"))

# plot GLM occ v pres 
#  + geom_label(data = auc_df_traits, aes(x = AUC, y = AUC_pres, label = ALPHA.CODE))
r1 = ggplot(auc_df_traits, aes(x = rmse_occ, y = rmse_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("GLM RMSE")) + ylab(bquote("Pres GLM RMSE"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = auc_df_traits$n_pres)) + scale_y_continuous(limit = c(0, 0.5)) + scale_x_continuous(limit = c(0, .5)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) 
# ggsave("Figures/Occ_Pres_labelled.pdf", height = 8, width = 12)
  
r2 = ggplot(auc_df_traits, aes(x = rmse_gam, y = rmse_gam_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("GAM RMSE")) + ylab(bquote("Pres GAM RMSE"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(shape=16, aes(size = auc_df_traits$n_pres))+  scale_y_continuous(limit = c(0, .5)) + scale_x_continuous(limit = c(0, .5))  + 
    theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
    guides(colour = guide_legend(override.aes = list(shape = 15))) +
    theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) 
#  ggsave("Figures/Occ_numPres_RF.pdf", height = 8, width = 12)

r3 = ggplot(auc_df_traits, aes(x = rmse_rf, y = rmse_rf_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("RF RMSE")) + ylab(bquote("Pres RF RMSE"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + scale_y_continuous(limit = c(0, .5)) + scale_x_continuous(limit = c(0, .5)) + geom_point(shape=16, aes(size = auc_df_traits$n_pres)) + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) + guides(colour = guide_legend(override.aes = list(shape = 15))) + theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) 
#  ggsave("Figures/Occ_numPres_gam.pdf", height = 8, width = 12)

# density plot
rmse_plot = gather(auc_df, occ_mod, occ_rmse, c("rmse_occ", "rmse_gam",  "rmse_rf")) %>%
  gather(pres_mod, pres_rmse, "rmse_pres","rmse_gam_pres","rmse_rf_pres",  "rmse_me_pres")

r4 = ggplot(rmse_plot) + geom_density(lwd = 1.5, lty = 2, aes(occ_rmse, color = occ_mod))  + geom_density(lwd = 1.5, aes(pres_rmse, color = pres_mod)) + theme_classic() + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) + theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + scale_color_manual(values=c("#034e7b","#034e7b", "purple", "steelblue2", "steelblue2","#238b45", "#238b45"), labels=c("rmse_gam", "rmse_gam_pres",  "rmse_me_pres", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres")) + xlab("RMSE") + ylab("Density") + guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_blank()) +theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm"))  


legend <- ggplot(rmse_plot, aes(pres_rmse, pres_mod)) + geom_line(lwd = 1.5, aes(color = pres_mod)) + scale_color_manual(values=c("#034e7b", "purple", "steelblue2", "#238b45"), labels=c("GAM",  "MaxEnt", "GLM",  "RF")) +theme(legend.spacing.x = unit(0.5, 'cm'), legend.key.size = unit(1, "cm"), legend.title = element_blank())


library(cowplot)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
grid <- plot_grid(r1 + theme(legend.position="none"),
          r2 + theme(legend.position="none"),
          r3 + theme(legend.position="none"),
          r4 + theme(legend.position="none"),
          align = 'hv',
          labels = c("A","B", "C", "D"),
          label_x = 0.22, 
          label_size = 30,
          nrow = 2) 


ggsave("Figures/rmse_plot_5.pdf", height = 10, width = 14)

# pres:pres
pres_pres <- left_join(auc_df_notrans, auc_df, by = "AOU") %>%
  left_join(num_pres, by = c("AOU" = "aou"))
# plot GLM occ v pres 
r1 = ggplot(pres_pres, aes(x = rmse_pres_notrans, y = rmse_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + ylab(bquote("Pres GLM RMSE")) + xlab(bquote("GLM No Transients"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = pres_pres$n))  + scale_y_continuous(limit = c(0, 0.5)) + scale_x_continuous(limit = c(0, .5)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) 
# ggsave("Figures/Occ_Pres_labelled.pdf", height = 8, width = 12)


r2 =  ggplot(pres_pres, aes(x = rmse_gam_pres_notrans, y = rmse_gam_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + ylab(bquote("Pres GAM RMSE")) + xlab(bquote("GAM No Transients"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = pres_pres$n))  + scale_y_continuous(limit = c(0, 0.5)) + scale_x_continuous(limit = c(0, .5)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) 
#  ggsave("Figures/Occ_numPres_RF.pdf", height = 8, width = 12)

r3 =  ggplot(pres_pres, aes(x = rmse_rf_pres_notrans, y = rmse_rf_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + ylab(bquote("Pres RF RMSE")) + xlab(bquote("RF No Transients"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = pres_pres$n))  + scale_y_continuous(limit = c(0, 0.5)) + scale_x_continuous(limit = c(0, .5)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) 
#  ggsave("Figures/Occ_numPres_gam.pdf", height = 8, width = 12)

r4 =  ggplot(pres_pres, aes(x = rmse_me_pres_notrans, y = rmse_me_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + ylab(bquote("Pres ME RMSE")) + xlab(bquote("ME No Transients"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = pres_pres$n)) + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) 


# density plot
rmse_plot_pres = gather(pres_pres, pres_mod, pres_rmse, c("rmse_pres","rmse_gam_pres", "rmse_rf_pres", "rmse_me_pres")) %>%
  gather(notrans_mod, notrans_rmse, c("rmse_pres_notrans","rmse_gam_pres_notrans", "rmse_rf_pres_notrans", "rmse_me_pres_notrans"))

r5 = ggplot(rmse_plot_pres) + geom_density(lwd = 1.5, aes(pres_rmse, color = pres_mod)) + geom_density(lwd = 1.5, lty = 4,aes(notrans_rmse, color = notrans_mod)) + theme_classic() + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) + theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + scale_color_manual(values=c("#034e7b","#034e7b", "purple","purple", "steelblue2", "steelblue2","#238b45", "#238b45"), labels=c("rmse_gam_pres","rmse_gam_pres_notrans",   "rmse_me_pres","rmse_me_pres_notrans", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres")) +  xlab("RMSE") + ylab("Density") + guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_blank()) 

# scale_color_manual(values=c("#034e7b","#034e7b", "purple", "steelblue2", "steelblue2","#238b45", "#238b45"), labels=c("rmse_gam", "rmse_gam_pres",  "rmse_me_pres", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres"))

legend <- ggplot(rmse_plot_pres, aes(rmse, mod)) + geom_line(lwd = 1.5, aes(color = mod)) + scale_color_manual(values=c("#006d2c","#ce1256","purple", "#9F84BD", "#66c2a4", "plum", "springgreen2", "#df65b0"), labels=c("GAM pres","GAM pres notrans",  "MaxEnt pres","MaxEnt pres notrans", "GLM pres","GLM pres notrans",  "RF pres", "RF pres notrans")) 
   

library(cowplot)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
plot_grid(r1 + theme(legend.position="none"),
          r2 + theme(legend.position="none"),
          r3 + theme(legend.position="none"),
          r4 + theme(legend.position="none"),
          r5 + theme(legend.position="none"),
          align = 'hv',
          labels = c("A","B", "C", "D", "E"),
          label_x = 0.22, 
          label_size = 30,
          nrow = 2, 
          scale = 0.9) 
ggsave("Figures/rmse_pres_pres_5.pdf", height = 14, width = 24)
