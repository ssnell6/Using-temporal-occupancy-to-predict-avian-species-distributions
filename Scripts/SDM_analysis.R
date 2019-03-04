library(dplyr)
library(glm2)

bbs_occ = read.csv("data/bbs_sub1.csv", header=TRUE)
exp_pres = read.csv("data/expect_pres.csv", header = TRUE)
exp_pres = exp_pres[,c("stateroute","spAOU")]
traits = read.csv("data/Master_RO_Correlates.csv", header = TRUE)
lat_long = read.csv("data/latlongs.csv", header = TRUE)
bi_env = read.csv("data/all_env.csv", header = TRUE)
bi_means = bi_env[,c("stateroute","mat.mean", "elev.mean", "map.mean", "ndvi.mean")]
env_bio = read.csv("data/env_bio.csv", header = TRUE)
env_bio = na.omit(env_bio)
env_bio_sub = env_bio[,c(1, 21:39)]
# env_bio[,c("stateroute","bio.mean.bio4", "bio.mean.bio5", "bio.mean.bio6", "bio.mean.bio13", "bio.mean.bio14")]

all_env = left_join(bi_means, env_bio_sub, by = "stateroute")

# BBS cleaning
bbs_inc_absence = full_join(bbs_occ, exp_pres, by = c("Aou" ="spAOU", "stateroute" = "stateroute"))
bbs_inc_absence$occ[is.na(bbs_inc_absence$occ)] <- 0
bbs_inc_absence$presence = 0
bbs_inc_absence$presence[bbs_inc_absence$occ > 0] <- 1
num_occ = bbs_inc_absence %>% group_by(Aou) %>% tally(presence) %>% left_join(bbs_inc_absence, ., by = "Aou")
# 586 focal species
bbs_final_occ = filter(num_occ,nn > 1)
bbs_final_occ_ll = left_join(bbs_final_occ, lat_long, by = "stateroute")
bbs_final_occ_ll = bbs_final_occ_ll[,c("Aou", "stateroute", "occ", "presence", 
                                       "latitude", "longitude")]
  
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
test = cor(na.omit(sdm_input_global))

corrplot(test)

# the XY coordinates of species data
coords <- CRS("+proj=longlat +datum=WGS84")

# running mods on each spp
# need to add in presence absence!
library(pROC)
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
  if(length(unique(sdm_input$stateroute)) > 40){
  glm_occ <- glm2(occ ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio12, family = gaussian(link = "identity"), data = sdm_input)
  glm_pres <- glm2(presence ~ elev.mean + ndvi.mean + bio.mean.bio4 + bio.mean.bio5 + bio.mean.bio6 + bio.mean.bio13 + bio.mean.bio14, family = gaussian(link = "identity"), data = sdm_input)
  predocc <- predict(glm_occ,type=c("response"))
  predpr <- predict(glm_pres,type=c("response"))
 sdm_output = cbind(sdm_input,predpr, predocc)
 
 roccurve <- roc(sdm_output$occ ~ sdm_output$predocc)
 auc =  roc(sdm_output$occ ~ sdm_output$predocc)$auc[1]
 
# rocpres <- roc(sdm_output$presence ~ sdm_output$predpr)
# auc_pres =  roc(sdm_output$presence ~ sdm_output$predpr)$auc[1]
 
 auc_df = rbind(auc_df, c(i, auc, auc_pres))
 plot = plot(roccurve, main = paste("AUC Curve for AOU", i, ".csv",   
                                     sep=""))
 # lines(rocpres, col = "red")
 }
 write.csv(sdm_output, paste("sdm_output_", i, ".csv",   
                            sep=""), row.names = FALSE)
}

dev.off()

auc_df = data.frame(auc_df)
write.csv(auc_df, "auc_df.csv", row.names = FALSE)
names(auc_df) = c("AOU", "AUC")
test = dplyr::filter(auc_df, AUC > 0.75 & AUC < 1.0)
