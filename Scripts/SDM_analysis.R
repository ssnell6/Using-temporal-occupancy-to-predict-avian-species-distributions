library(dplyr)
library(glm2)
library(gam)
library(randomForest)
library(ggplot2)
library(dismo)
library(raster)
library(maptools)
library(pROC)

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
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_201') # for 64-bit version


# test = cor(na.omit(sdm_input_global))
# corrplot(test)

setwd("Data/sdm_dfs/")
pdf('AUC_Curves.pdf', height = 8, width = 10)
par(mfrow = c(2, 3))
auc_df = c()

sp_list = unique(bbs_final_occ_ll$Aou)

for(i in sp_list){
  sdm_output = c()
  bbs_sub <- filter(bbs_final_occ_ll, Aou == i)
  temp <- filter(all_env, stateroute %in% bbs_sub$stateroute)
  sdm_input <- left_join(bbs_sub, temp, by = "stateroute")
  sdm_input = na.omit(sdm_input)
  # print(length(sdm_input$stateroute))
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
 #   if(levels(as.factor(sdm_input$presence)) > 1){
    glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio2 + bio.mean.bio3 + bio.mean.bio4 + bio.mean.bio7 + bio.mean.bio8 + bio.mean.bio9 + bio.mean.bio10 + bio.mean.bio11 + bio.mean.bio12 + bio.mean.bio13 + bio.mean.bio14 +bio.mean.bio15+ bio.mean.bio16 +bio.mean.bio17 +bio.mean.bio18 +bio.mean.bio19, family = binomial(link = logit), data = sdm_input)
    glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio2 + bio.mean.bio3 + bio.mean.bio4 + bio.mean.bio7 + bio.mean.bio8 + bio.mean.bio9 + bio.mean.bio10 + bio.mean.bio11 + bio.mean.bio12 + bio.mean.bio13 + bio.mean.bio14 +bio.mean.bio15+ bio.mean.bio16 +bio.mean.bio17 +bio.mean.bio18 +bio.mean.bio19, family = binomial(link = logit), data = sdm_input)
    gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio1) + s(bio.mean.bio2) + s(bio.mean.bio3) + s(bio.mean.bio4) + s(bio.mean.bio7) + s(bio.mean.bio8) + s(bio.mean.bio9) + s(bio.mean.bio10) + s(bio.mean.bio11) + s(bio.mean.bio12) + s(bio.mean.bio13) + s(bio.mean.bio14) + s(bio.mean.bio15) + s(bio.mean.bio16) + s(bio.mean.bio17) + s(bio.mean.bio18) + s(bio.mean.bio19), family = binomial(link = logit), data = sdm_input)
    gam_pres <- mgcv::gam(presence ~  s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio1) + s(bio.mean.bio2) + s(bio.mean.bio3) + s(bio.mean.bio4) + s(bio.mean.bio7) + s(bio.mean.bio8) + s(bio.mean.bio9) + s(bio.mean.bio10) + s(bio.mean.bio11) + s(bio.mean.bio12) + s(bio.mean.bio13) + s(bio.mean.bio14) + s(bio.mean.bio15) + s(bio.mean.bio16) + s(bio.mean.bio17) + s(bio.mean.bio18) + s(bio.mean.bio19), family = binomial(link = logit), data = sdm_input)
    rf_occ <- randomForest(sp_success ~elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio2 + bio.mean.bio3 + bio.mean.bio4 + bio.mean.bio7 + bio.mean.bio8 + bio.mean.bio9 + bio.mean.bio10 + bio.mean.bio11 + bio.mean.bio12 + bio.mean.bio13 + bio.mean.bio14 +bio.mean.bio15+ bio.mean.bio16 +bio.mean.bio17 +bio.mean.bio18 +bio.mean.bio19, family = binomial(link = logit), data = sdm_input)
    rf_pres <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio2 + bio.mean.bio3 + bio.mean.bio4 + bio.mean.bio7 + bio.mean.bio8 + bio.mean.bio9 + bio.mean.bio10 + bio.mean.bio11 + bio.mean.bio12 + bio.mean.bio13 + bio.mean.bio14 +bio.mean.bio15+ bio.mean.bio16 +bio.mean.bio17 +bio.mean.bio18 +bio.mean.bio19, family = binomial(link = logit), data = sdm_input)
   
    max_pres = SpatialPointsDataFrame(coords = sdm_input[,c("longitude", "latitude")],
     data = sdm_input[,c("longitude", "latitude", "presence")], 
     proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
   # max_ind_occ = dismo::maxent(sdm_input[,c("bio.mean.bio1", "elev.mean", "bio.mean.bio2", "ndvi.mean")], sdm_input$presence)
    max_ind_pres = maxent(sdm_input[,c("bio.mean.bio1", "elev.mean", "bio.mean.bio2","bio.mean.bio3","bio.mean.bio4","bio.mean.bio7","bio.mean.bio8","bio.mean.bio9","bio.mean.bio10","bio.mean.bio11","bio.mean.bio12","bio.mean.bio13","bio.mean.bio14","bio.mean.bio15","bio.mean.bio16","bio.mean.bio17","bio.mean.bio18","bio.mean.bio19", "ndvi.mean")], sdm_input$presence)
    
    # predict
    pred_glm_occ <- predict(glm_occ,type=c("response"))
    pred_glm_pr <- predict(glm_pres,type=c("response"))
    pred_gam_occ <- predict(gam_occ,type=c("response"))
    pred_gam_pr <- predict(gam_pres,type=c("response"))
    pred_rf_occ <- predict(rf_occ,type=c("response"))
    pred_rf_pr <- predict(rf_pres,type=c("response"))
    max_pred_pres <- predict(max_ind_pres, sdm_input[,c("bio.mean.bio1", "elev.mean", "bio.mean.bio2","bio.mean.bio3","bio.mean.bio4","bio.mean.bio7","bio.mean.bio8","bio.mean.bio9","bio.mean.bio10","bio.mean.bio11","bio.mean.bio12","bio.mean.bio13","bio.mean.bio14","bio.mean.bio15","bio.mean.bio16","bio.mean.bio17","bio.mean.bio18","bio.mean.bio19", "ndvi.mean")])
   # max_pred_occ <- predict(max_ind_occ, sdm_input[,c("bio.mean.bio1", "elev.mean", "bio.mean.bio2", "ndvi.mean")])
    
    sdm_output = cbind(sdm_input,pred_glm_pr, pred_glm_occ, pred_gam_pr, pred_gam_occ, pred_rf_occ, pred_rf_pr, max_pred_pres) 
 
 roccurve <- roc(sdm_output$occ ~ sdm_output$pred_glm_occ)
 auc =  roc(sdm_output$occ ~ sdm_output$pred_glm_occ)$auc[1]
 
 rocpres <- roc(sdm_output$presence ~ sdm_output$pred_glm_pr)
 auc_pres =  roc(sdm_output$presence ~ sdm_output$pred_glm_pr)$auc[1]
 
 roccurve_gam <- roc(sdm_output$occ ~ sdm_output$pred_gam_occ)
 auc_gam =  roc(sdm_output$occ ~ sdm_output$pred_gam_occ)$auc[1]
 
 rocpres_gam <- roc(sdm_output$presence ~ sdm_output$pred_gam_pr)
 auc_pres_gam =  roc(sdm_output$presence ~ sdm_output$pred_gam_pr)$auc[1]
 
 
 roccurve_rf <- roc(sdm_output$occ ~ sdm_output$pred_rf_occ)
 auc_rf =  roc(sdm_output$occ ~ sdm_output$pred_rf_occ)$auc[1]
 
 rocpres_rf <- roc(sdm_output$presence ~ sdm_output$pred_rf_pr)
 auc_pres_rf =  roc(sdm_output$presence ~ sdm_output$pred_rf_pr)$auc[1]
 
 rocpres_me <- roc(sdm_output$presence ~ sdm_output$max_pred_pres)
 auc_pres_me =  roc(sdm_output$presence ~ sdm_output$max_pred_pres)$auc[1]
 
 auc_df = rbind(auc_df, c(i, auc, auc_pres, auc_rf, auc_gam, auc_pres_gam, auc_pres_rf, auc_pres_me))
 j = unique(sdm_input$ALPHA.CODE)
 plot = plot(roccurve, main = paste("AUC Curve for ", j, ".csv",   
                                     sep=""))
  }
 write.csv(sdm_output, paste("sdm_output_", i, ".csv",   
                            sep=""), row.names = FALSE)
}

dev.off()

setwd("C:/Git/SDMs")
auc_df = data.frame(auc_df)
names(auc_df) = c("AOU", "AUC", "AUC_pres","AUC_RF", "AUC_gam", "AUC_gam_pres",  "AUC_RF_pres", "AUC_me_pres")
# write.csv(auc_df, "Data/auc_df.csv", row.names = FALSE)
test = dplyr::filter(auc_df, AUC > 0.75 & AUC < 1.0)


bbs_final_occ_ll$presence <- factor(bbs_final_occ_ll$presence,levels = c('1','0'), ordered = TRUE)

#### MAPS #####
auc_df = read.csv("Data/auc_df.csv", header = TRUE)
setwd("Figures/maps/")
# temp filter for vis purposes
sp_list_bigauc = filter(bbs_final_occ_ll, Aou %in% test$AOU)

mapfun <- function(pdf_name, vec, num){ 

  pdf(pdf_name, height = 8, width = 10)
  par(mfrow = c(2, 3)) # makes plots too small

for(i in unique(sp_list)){
  print(i)
  sdm_output = c()
  
  bbs_sub <- filter(bbs_final_occ_ll, Aou == i)
  temp <- filter(all_env, stateroute %in% bbs_sub$stateroute)
  sdm_input <- left_join(bbs_sub, temp, by = "stateroute")
  sdm_input = na.omit(sdm_input)
  j = unique(sdm_input$ALPHA.CODE)
  
  # print(length(sdm_input$stateroute))
  if(length(unique(sdm_input$stateroute)) > 40  & length(unique(sdm_input$presence)) >1){
    #   if(levels(as.factor(sdm_input$presence)) > 1){
    glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio12, family = binomial(link = logit), data = sdm_input)
    glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio12, family = binomial(link = logit), data = sdm_input)
    
    rf_occ <- randomForest(sp_success ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio12, family = binomial(link = logit), data = sdm_input)
    rf_pres <- randomForest(presence ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio12, family = binomial(link = logit), data = sdm_input)
    
    # max_pres = max_env = SpatialPointsDataFrame(coords = sdm_input[,c("longitude", "latitude")],
     #  data = sdm_input[,c("longitude", "latitude", "presence")], 
     #  proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
   #  max_ind_pres = maxent(sdm_input[,c("bio.mean.bio1", "elev.mean", "bio.mean.bio2", "ndvi.mean")], sdm_input$presence)
    
    pred_glm_occ <- predict(glm_occ,type=c("response"))
    pred_glm_pr <- predict(glm_pres,type=c("response"))
    
    pred_rf_occ <- predict(rf_occ,type=c("response"))
    pred_rf_pr <- predict(rf_pres,type=c("response"))
    
    sdm_output = cbind(sdm_input,pred_glm_pr, pred_glm_occ, pred_rf_occ, pred_rf_pr) 
  
    # Determine geographic extent of our data using AOU = i
  max.lat <- ceiling(max(sdm_input$latitude))
  min.lat <- floor(min(sdm_input$latitude))
  max.lon <- ceiling(max(sdm_input$longitude))
  min.lon <- floor(min(sdm_input$longitude))

  geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))


  mod.r <- SpatialPointsDataFrame(coords = sdm_output[,c("longitude", "latitude")],
   data = sdm_output[,c("latitude", "longitude", "pred_glm_pr", "pred_glm_occ", "pred_rf_occ")], 
   proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
  r = raster(mod.r, res = 0.6) # 40x40 km/111 (degrees) * 2 tp eliminate holes
  # bioclim is 4 km
  plot.r = rasterize(mod.r, r)

  data(wrld_simpl)
  # Plot the base map
  plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95", main = paste("SDM plot for ", j, sep=""))

  plot(plot.r[[num]], add = TRUE)

  # Add the points for individual observation if necessary
  sdm_input$presence <-droplevels(sdm_input$presence, exclude = c("0"))
  # sdm_input$col = c("black", "white")
  points(x = sdm_input$longitude, y = sdm_input$latitude, col = sdm_input$presence, pch = 20, cex = 0.75)

  box()
  }

}
  dev.off()
}
# dev.off()

# there is still something weird about the dev.off(), need to run dev off or restart R to get running
mapfun(pdf_name = 'SDM_glm_pres_maps.pdf',vec = "pred_glm_pr", 4)

mapfun(pdf_name = 'SDM_glm_occ_maps.pdf', vec = "pred_glm_occ", 5)

# mapfun(pdf_name = 'SDM_rf_pres_maps.pdf', vec = "pred_rf_pr", mod = plot.r$pred_rf_pr)

mapfun(pdf_name = 'SDM_rf_occ_maps.pdf', vec = "pred_rf_occ", 6)

mapfun(pdf_name = 'SDM_me_pres_maps.pdf', vec = "pred_me_pres", mod = plot.r$pred_me_pres)

setwd("C:/Git/SDMs")

env.data = as.matrix(sdm_input[,c("latitude", "longitude", "elev.mean", "ndvi.mean", "bio.mean.bio1", "bio.mean.bio12")])

env.proj = SpatialPointsDataFrame(coords = sdm_input[,c("longitude", "latitude")],
                                  data = sdm_input[,c("latitude", "longitude", "elev.mean", "ndvi.mean", "bio.mean.bio1", "bio.mean.bio12")], 
                                  proj4string = CRS("+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"))
r = raster(env.proj)
env.proj.raster = rasterize(env.proj, r)
env.stack = raster::stack(env.proj.raster@data$elev.mean, env.proj.raster@data$ndvi.mean)

auc_df_traits$diff = auc_df_traits$AUC - auc_df_traits$AUC_pres
occ_sub = subset(auc_df_traits,diff >= 0)
pres_sub = subset(auc_df_traits,diff < 0)

num_pres = bbs_final_occ_ll %>%
  group_by(Aou) %>% 
  filter(., presence == "1") %>%
  summarise(n = n_distinct(stateroute))  

num_routes = bbs_final_occ_ll %>% group_by(Aou) %>% 
  summarise(n = n_distinct(stateroute))  %>%
  filter(., n >= 40)

auc_df_merge = left_join(auc_df, num_pres, by = c("AOU" = "Aou"))
auc_df_traits = left_join(auc_df_merge, traits, by = "AOU") %>%
  left_join(., tax_code, by = c("AOU" = "AOU_OUT"))
# plot GLM occ v pres 
#  + geom_label(data = auc_df_traits, aes(x = AUC, y = AUC_pres, label = ALPHA.CODE))
r1 = ggplot(auc_df_traits, aes(x = AUC, y = AUC_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + xlab(bquote("Occupancy AUC")) + ylab(bquote("Presence AUC"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5)  + 
  geom_point(shape=16, aes(size = auc_df_traits$n))  + geom_smooth(method='lm', se=FALSE, col="blue",linetype="longdash", lwd =2.5) + scale_y_continuous(limit = c(.5, 1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) + scale_x_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(axis.text.x=element_text(size = 32),axis.ticks=element_blank(), axis.text.y=element_text(size=32)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines")) 
# ggsave("Figures/Occ_Pres_labelled.pdf", height = 8, width = 12)
  
  
r2 = ggplot(auc_df_traits, aes(x = AUC_RF, y = AUC_RF_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + xlab(bquote("Occupancy AUC")) + ylab(bquote("Presence AUC"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + scale_y_continuous(limit = c(.5, 1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) + scale_x_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))  + 
    geom_point(shape=16, aes(size = auc_df_traits$n))  + geom_smooth(method='lm', se=FALSE, col="blue",linetype="longdash", lwd =2.5) + theme(axis.text.x=element_text(size = 32),axis.ticks=element_blank(), axis.text.y=element_text(size=32)) +
    guides(colour = guide_legend(override.aes = list(shape = 15))) +
    theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines"))
#  ggsave("Figures/Occ_numPres_RF.pdf", height = 8, width = 12)

r2 = ggplot(auc_df_traits, aes(x = AUC_gam, y = AUC_gam_pres)) +theme_classic()+ theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + xlab(bquote("Occupancy AUC")) + ylab(bquote("Presence AUC"))+ geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + scale_y_continuous(limit = c(.5, 1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) + scale_x_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))  + 
  geom_point(shape=16, aes(size = auc_df_traits$n))  + geom_smooth(method='lm', se=FALSE, col="blue",linetype="longdash", lwd =2.5) + theme(axis.text.x=element_text(size = 32),axis.ticks=element_blank(), axis.text.y=element_text(size=32)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines"))
#  ggsave("Figures/Occ_numPres_gam.pdf", height = 8, width = 12)


# density plot
auc_plot = gather(auc_df, mod, AUC, AUC:AUC_me_pres)
ggplot(auc_plot, aes(AUC, color = mod)) + geom_density(lwd = 1.5)  + theme_classic() + theme(axis.text.x=element_text(size = 32),axis.ticks=element_blank(), axis.text.y=element_text(size=32)) + theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + guides(colour = guide_legend(override.aes = list(shape = 15))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=15), legend.position = c(0.1,0.9), legend.key.width=unit(2, "lines"))
#  ggsave("Figures/density_mod_comp.pdf", height = 9, width = 12)
