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
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_211') # for 64-bit version

###### main analysis ######
# test = cor(na.omit(sdm_input_global))
# corrplot(test)

auc_df = c()

  sp_list = unique(comp_plot$AOU)
    # unique(bbs_final_occ_ll$aou) change back when done w overfit

for(i in sp_list){
  sdm_output = c()
  print(i)
  bbs_sub <- filter(bbs_final_occ_ll, aou == i) 
  bbs_new_sub <- filter(bbs_new, aou == i) 
  bbs_new_sub$pres_2016 <- bbs_new_sub$presence
  sdm_input <- filter(all_env, stateroute %in% bbs_sub$stateroute) %>%
    full_join(bbs_sub, by = "stateroute") 
  sdm_input <- data.frame(sdm_input)
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
  if(nrow(filter(sdm_input, presence == 1)) > 59){
    glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
    glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
    gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_input)
    gam_pres <- mgcv::gam(presence ~   s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14), family = binomial(link = logit), data = sdm_input)
    rf_occ <- randomForest(sp_success/15 ~elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
    rf_pres <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = binomial(link = logit), data = sdm_input)
    
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
    
    sdm_output = cbind(sdm_input, pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr) 
    pred_2016 <- left_join(sdm_output, bbs_new_sub[c("stateroute", "pres_2016")], by = "stateroute")
    thresh <- max(pred_2016$pred_gam_occ) * 0.7
    gam_rescale <- filter(pred_2016, pred_gam_occ > thresh)
    
 rmse_occ <- rmse(sdm_output$pred_glm_occ, sdm_output$occ)
 rmse_pres <- rmse(sdm_output$pred_glm_pr, sdm_output$presence)
 rmse_gam <- rmse(as.vector(sdm_output$pred_gam_occ), sdm_output$occ)
 rmse_gam_pres <- rmse(as.vector(sdm_output$pred_gam_pr), sdm_output$presence)
 rmse_rf <- rmse(sdm_output$pred_rf_occ, sdm_output$occ)
 rmse_rf_pres <- rmse(as.vector(as.numeric(sdm_output$pred_rf_pr)), sdm_output$presence)
 auc_df = rbind(auc_df, c(i, rmse_occ, rmse_pres, rmse_gam, rmse_gam_pres, rmse_rf, rmse_rf_pres))
 j = unique(sdm_input$ALPHA.CODE)
  }
  }
# write.csv(sdm_output, paste("sdm_output_notrans_", i, ".csv",  sep=""), row.names = FALSE)
}


auc_df = data.frame(auc_df)
names(auc_df) = c("AOU", "rmse_occ", "rmse_pres","rmse_gam", "rmse_gam_pres", "rmse_rf", "rmse_rf_pres","rmse_me_pres")
# write.csv(auc_df, "Data/auc_df.csv", row.names = FALSE)

###### global plots ####
me <- read.csv("Data/auc_df_ME_only.csv", header = TRUE)
#  mutate(ME_PB = rmse_me_pres)
auc_df <- read.csv("Data/auc_df.csv", header = TRUE) %>%
  dplyr::select(-rmse_me_pres) %>%
  left_join(me, by = "AOU")

ggplot(
  data = auc_df, aes(x = rmse_me_pres, y = ME_PB)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "blue", lwd = 1.5) +
  xlab("old ME RMSE") + 
  ylab("new ME RMSE")

# auc_df_notrans <- read.csv("Data/auc_df_notrans.csv", header = TRUE)

num_pres = bbs_final_occ_ll %>%
  group_by(aou) %>% 
  dplyr::filter(., presence == "1") %>%
  dplyr::summarise(n_pres = n_distinct(stateroute))  

num_abs = bbs_final_occ_ll %>%
  group_by(aou) %>% 
  dplyr::filter(., presence == "0") %>%
  dplyr::summarise(n_abs = n_distinct(stateroute))  


num_routes = bbs_final_occ_ll %>% group_by(aou) %>% 
  dplyr::summarise(n = n_distinct(stateroute)) 

auc_df_merge = left_join(auc_df, num_pres, by = c("AOU" = "aou")) %>%
  left_join(num_routes,  by = c("AOU" = "aou")) %>%
  left_join(num_abs, by = c("AOU" = "aou"))

auc_df_traits = left_join(auc_df_merge, traits, by = "AOU") %>%
  left_join(., tax_code, by = c("AOU" = "AOU_OUT")) %>%
  filter(!AOU %in% c(3250, 3390))

auc_df_merge$glm_diff <- auc_df_merge$rmse_occ - auc_df_merge$rmse_pres
auc_df_merge$gam_diff <- auc_df_merge$rmse_gam - auc_df_merge$rmse_gam_pres
auc_df_merge$rf_diff <- auc_df_merge$rmse_rf - auc_df_merge$rmse_rf_pres
auc_df_merge$RO <- auc_df_merge$n_pres/auc_df_merge$n


# plot GLM occ v pres 
#  + geom_label(data = auc_df_traits, aes(x = AUC, y = AUC_pres, label = ALPHA.CODE))
r1 = ggplot(auc_df_traits, aes(x = rmse_occ, y = rmse_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -1),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Occupancy RMSE")) + ylab(bquote("Presence RMSE"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = auc_df_traits$n_pres)) + scale_y_continuous(limit = c(0, 0.5)) + scale_x_continuous(limit = c(0, .5)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm")) +
  annotate("text", x = 0.45, y = 0.02, label = "GLM", size = 10) 
# ggsave("Figures/Occ_Pres_labelled.pdf", height = 8, width = 12)
  
r2 = ggplot(auc_df_traits, aes(x = rmse_gam, y = rmse_gam_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -1),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Occupancy RMSE")) + ylab(bquote("Presence RMSE"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(shape=16, aes(size = auc_df_traits$n_pres))+  scale_y_continuous(limit = c(0, .5)) + scale_x_continuous(limit = c(0, .5))  +
    theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
    guides(colour = guide_legend(override.aes = list(shape = 15))) +
    theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm"))  +
  annotate("text", x = 0.45, y = 0.01, label = "GAM", size = 10)
#  ggsave("Figures/Occ_numPres_RF.pdf", height = 8, width = 12)

r3 = ggplot(auc_df_traits, aes(x = rmse_rf, y = rmse_rf_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -1),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + xlab(bquote("Occupancy RMSE")) + ylab(bquote("Presence RMSE"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + scale_y_continuous(limit = c(0, .5)) + scale_x_continuous(limit = c(0, .5)) + geom_point(shape=16, aes(size = auc_df_traits$n_pres)) + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) + guides(colour = guide_legend(override.aes = list(shape = 15))) + theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) + theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm"))  +
  annotate("text", x = 0.45, y = 0.01, label = "RF", size = 10)
#  ggsave("Figures/Occ_numPres_gam.pdf", height = 8, width = 12)

# density plot
rmse_plot = gather(auc_df, occ_mod, occ_rmse, c("rmse_occ", "rmse_gam",  "rmse_rf")) %>%
  gather(pres_mod, pres_rmse, "rmse_pres","rmse_gam_pres","rmse_rf_pres",  "rmse_me_pres")

r4 = ggplot(rmse_plot) + geom_density(lwd = 1.5, lty = 2, aes(occ_rmse, color = occ_mod))  + geom_density(lwd = 1.5, aes(pres_rmse, color = pres_mod)) + theme_classic() + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) + theme(axis.title.x=element_text(size=34, vjust = -1),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + scale_color_manual(values=c("#034e7b","#034e7b", "purple", "steelblue2", "steelblue2","#238b45", "#238b45"), labels=c("rmse_gam", "rmse_gam_pres",  "rmse_me_pres", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres")) + xlab("RMSE") + ylab("Density") + guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_blank()) +theme(plot.margin=unit(c(1.2,1.2,1.2,1.2),"cm"))  


legend <- ggplot(rmse_plot, aes(pres_rmse, pres_mod)) + geom_line(lwd = 1.5, aes(color = pres_mod)) + scale_color_manual(values=c("#034e7b", "purple", "steelblue2", "#238b45"), labels=c("GAM",  "MaxEnt", "GLM",  "RF")) +theme(legend.spacing.x = unit(0.5, 'cm'), legend.key.size = unit(1, "cm"), legend.title = element_blank())

rmse_plot_sub <- filter(rmse_plot, pres_mod == "rmse_pres" | pres_mod == "rmse_gam_pres")
rmse_plot_sub$pres_mod[rmse_plot_sub$pres_mod == "rmse_pres"] <- "Temporal Occupancy"
rmse_plot_sub$pres_mod[rmse_plot_sub$pres_mod == "rmse_gam_pres"] <- "Presence"
rmse_plot_sub$pres_mod <- factor(rmse_plot_sub$pres_mod, levels = c("Presence", "Temporal Occupancy"), ordered = TRUE)
legend_bw <- ggplot(rmse_plot_sub, aes(pres_rmse, pres_mod)) + geom_line(lwd = 1.5, aes(lty = pres_mod)) +theme(legend.spacing.x = unit(0.5, 'cm'), legend.key.size = unit(1, "cm"), legend.title = element_blank())


library(cowplot)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
grid <- plot_grid(r1 + theme(legend.position="top"),
          r2 + theme(legend.position="none"),
          r3 + theme(legend.position="none"),
          r4 + theme(legend.position="none"),
          align = 'hv',
          labels = c("A","B", "C", "D"),
          label_x = 0.22, 
          label_size = 30,
          nrow = 2) 


ggsave("Figures/Figure_2.pdf", height = 10, width = 14)

# t test for each occ:pres
t.test(auc_df_traits$rmse_occ, auc_df_traits$rmse_pres)
t.test(auc_df_traits$rmse_gam, auc_df_traits$rmse_gam_pres)
t.test(auc_df_traits$rmse_rf, auc_df_traits$rmse_rf_pres)
# pres:pres
pres_pres1 <- left_join(auc_df[,c("AOU","rmse_pres","rmse_gam_pres", "rmse_rf_pres", "rmse_me_pres")], by = "AOU") %>%
  left_join(num_pres, by = c("AOU" = "aou"))


pres_pres1$glm_diff <-  pres_pres1$rmse_pres_notrans - pres_pres1$rmse_pres
pres_pres1$gam_diff <-  pres_pres1$rmse_gam_pres_notrans - pres_pres1$rmse_gam_pres
pres_pres1$rf_diff <-  pres_pres1$rmse_rf_pres_notrans - pres_pres1$rmse_rf_pres
# pres_pres1$me_diff <-  pres_pres1$rmse_me_pres_notrans - pres_pres1$rmse_me_pres

# plot GLM occ v pres 
r1 = ggplot(pres_pres1, aes(x = rmse_pres_notrans, y = rmse_pres)) +
  theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -1),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + 
  ylab(bquote("Presence RMSE")) + xlab(bquote("Presence No Transients"))+ 
  geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = pres_pres1$n_pres))  + 
  scale_y_continuous(limit = c(0, 0.5)) + scale_x_continuous(limit = c(0, .5)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines"))  +
  annotate("text", x = 0.45, y = 0.02, label = "GLM", size = 10)
# ggsave("Figures/Occ_Pres_labelled.pdf", height = 8, width = 12)

r2 =  ggplot(pres_pres1, aes(x = rmse_gam_pres_notrans, y = rmse_gam_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -1),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + ylab(bquote("Presence RMSE")) + xlab(bquote("Presence No Transients"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = pres_pres1$n_pres))  + scale_y_continuous(limit = c(0, 0.5)) + scale_x_continuous(limit = c(0, .5)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) +annotate("text", x = 0.45, y = 0.02, label = "GAM", size = 10)
#  ggsave("Figures/Occ_numPres_RF.pdf", height = 8, width = 12)

r3 =  ggplot(pres_pres1, aes(x = rmse_rf_pres_notrans, y = rmse_rf_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -1),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + ylab(bquote("Presence RMSE")) + xlab(bquote("Presence No Transients"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = pres_pres1$n_pres))  + scale_y_continuous(limit = c(0, 0.5)) + scale_x_continuous(limit = c(0, .5)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) + annotate("text", x = 0.45, y = 0.02, label = "RF", size = 10)
#  ggsave("Figures/Occ_numPres_gam.pdf", height = 8, width = 12)

r4 =  ggplot(pres_pres1, aes(x = rmse_me_pres_notrans, y = rmse_me_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -4),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + ylab(bquote("Pres ME RMSE")) + xlab(bquote("ME No Transients"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  +
  geom_point(shape=16, aes(size = pres_pres1$n_pres)) + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines"))


# density plot
rmse_plot_pres = gather(pres_pres1, pres_mod, pres_rmse, c("rmse_pres","rmse_gam_pres", "rmse_rf_pres", "rmse_me_pres")) %>%
  gather(notrans_mod, notrans_rmse, c("rmse_pres_notrans","rmse_gam_pres_notrans", "rmse_rf_pres_notrans", "rmse_me_pres_notrans"))

r5 = ggplot(rmse_plot_pres) + geom_density(lwd = 1.5, aes(pres_rmse, color = pres_mod)) + geom_density(lwd = 1.5, lty = 2,aes(notrans_rmse, color = notrans_mod)) + theme_classic() + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) + theme(axis.title.x=element_text(size=34, vjust = -1),axis.title.y=element_text(size=34, angle=90, vjust = 5)) + scale_color_manual(values=c("#034e7b","#034e7b", "purple","purple", "steelblue2", "steelblue2","#238b45", "#238b45"), labels=c("rmse_gam_pres","rmse_gam_pres_notrans",   "rmse_me_pres","rmse_me_pres_notrans", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres")) +  xlab("RMSE") + ylab("Density") + guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_blank()) 

# scale_color_manual(values=c("#034e7b","#034e7b", "purple", "steelblue2", "steelblue2","#238b45", "#238b45"), labels=c("rmse_gam", "rmse_gam_pres",  "rmse_me_pres", "rmse_occ", "rmse_pres", "rmse_rf", "rmse_rf_pres"))

legend <- ggplot(rmse_plot_pres, aes(rmse, mod)) + geom_line(lwd = 1.5, aes(color = mod)) + scale_color_manual(values=c("#006d2c","#ce1256","purple", "#9F84BD", "#66c2a4", "plum", "springgreen2", "#df65b0"), labels=c("GAM pres","GAM pres notrans",  "MaxEnt pres","MaxEnt pres notrans", "GLM pres","GLM pres notrans",  "RF pres", "RF pres notrans")) 
   

library(cowplot)
theme_set(theme_cowplot(font_size=24,font_family = "URWHelvetica"))
plot_grid(r1 + theme(legend.position="top"),
          r2 + theme(legend.position="none"),
          r3 + theme(legend.position="none"),
          r4 + theme(legend.position="none"),
          r5 + theme(legend.position="none"),
          align = 'hv',
          labels = c("A","B", "C", "D", "E"),
          label_x = 0.02, 
          label_size = 30,
          nrow = 2, 
          scale = 0.9) 
# ggsave("Figures/rmse_pres_pres.pdf", height = 14, width = 20)

# experimental figure 3
##### stacked bar chart
auc_df_plot <- gather(auc_df_merge, "glm_mod", "val", c(rmse_occ, rmse_pres, glm_diff))
ggplot(data=auc_df_plot, aes(factor(AOU), y=val, fill=factor(glm_mod))) + 
  geom_bar(stat = "identity") + theme_classic() 

#### diff vs RO
auc_df_merge$sign <- ifelse(auc_df_merge$glm_diff >= 0, "pos", "neg")
glm <- ggplot(auc_df_merge, aes(x = RO, y = glm_diff)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -2),axis.title.y=element_text(size=34, angle=90, vjust = 3)) + xlab(bquote("RO")) + ylab(bquote("TO - Presence RMSE")) + 
  geom_point(shape=16, size = 3, aes(color = sign)) + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines"))  +
  scale_color_manual(values = c("steelblue2", "purple4")) +
  annotate("text", x = 0.9, y = -0.3, label = "GLM", size = 10)
auc_df_merge$sign <- ifelse(auc_df_merge$gam_diff >= 0, "pos", "neg")
gam <- ggplot(auc_df_merge, aes(x = RO, y = gam_diff)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -2),axis.title.y=element_text(size=34, angle=90, vjust = 3)) + xlab(bquote("RO")) + ylab(bquote("TO - Presence RMSE")) + 
  geom_point(shape=16,size = 3, aes(color = sign)) + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines"))  +
  scale_color_manual(values = c("#034e7b", "purple4")) +
  annotate("text", x = 0.9, y = -0.4, label = "GAM", size = 10)
auc_df_merge$sign <- ifelse(auc_df_merge$rf_diff >= 0, "pos", "neg")
rf <- ggplot(auc_df_merge, aes(x = RO, y = rf_diff)) +theme_classic()+ theme(axis.title.x=element_text(size=34, vjust = -2),axis.title.y=element_text(size=34, angle=90, vjust = 3)) + xlab(bquote("RO")) + ylab(bquote("TO - Presence RMSE")) + 
  geom_point(shape=16, size = 3,aes(color = sign)) + theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) + 
  scale_color_manual(values = c("#238b45", "purple4")) +
  annotate("text", x = 0.9, y = -0.4, label = "RF", size = 10)
theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
plot_grid(glm + theme(legend.position="none"),
          gam + theme(legend.position="none"),
          rf + theme(legend.position="none"),
          align = 'hv',
          labels = c("A","B", "C"),
          label_x = 0.22, 
          label_size = 30,
          nrow = 1, 
          scale = 0.8) 
# ggsave("Figures/RO_v_diff.pdf", height = 10, width = 14)

open_holes <- filter(ro_plot, rmse_gam_pres < 1.0e-5)

ro_plot <- gather(auc_df_merge, "mod", "diff", c(glm_diff, gam_diff, rf_diff))
ro_plot$mod[ro_plot$mod == "gam_diff"] = "GAM"
ro_plot$mod[ro_plot$mod == "glm_diff"] = "GLM"
ro_plot$mod[ro_plot$mod == "rf_diff"] = "RF"

ro_plot %>%
  filter(!AOU %in% c(3250, 3390)) %>%
  ggplot(aes(x = RO, y = diff, group = mod)) +theme_classic()+ 
  theme(axis.title.x=element_text(size=34),axis.title.y=element_text(size=34, angle=90)) + xlab(bquote("Range Occupancy")) + 
  ylab(expression(Delta~"RMSE (Occupancy - Presence)")) + 
  geom_point(size = 6, aes(color = mod, shape = mod)) + 
  geom_point(data = open_holes, aes(color = mod, shape = mod)) +
  theme(axis.text.x=element_text(size = 30),axis.ticks=element_blank(), axis.text.y=element_text(size=35)) +
  geom_hline(yintercept = 0, lty = 2, lwd =1.5, color = "black") +
  theme(legend.title=element_blank(), legend.text=element_text(size=30), legend.position = c(0.1,0.9), legend.key.width=unit(4, "lines")) + 
  scale_color_manual(values = c("#034e7b","steelblue2","#238b45")) 
ggsave("Figures/Figure_4.pdf", height = 10, width = 14)
