library(dplyr)
library(glm2)
library(ggplot2)
library(dismo)
library(maptools)

bbs_occ = read.csv("data/bbs_sub1.csv", header=TRUE)
bbs_occ_sub = bbs_occ %>% filter(Aou > 2880) %>%
  filter(Aou < 3650 | Aou > 3810) %>%
  filter(Aou < 3900 | Aou > 3910) %>%
  filter(Aou < 4160 | Aou > 4210) %>%
  filter(Aou != 7010)

exp_pres = read.csv("data/expect_pres.csv", header = TRUE)
exp_pres = exp_pres[,c("stateroute","spAOU")]
traits = read.csv("data/Master_RO_Correlates.csv", header = TRUE)
lat_long = read.csv("data/latlongs.csv", header = TRUE)
tax_code = read.csv("data/Tax_AOU_Alpha.csv", header = TRUE)
bi_env = read.csv("data/all_env.csv", header = TRUE)
bi_means = bi_env[,c("stateroute","mat.mean", "elev.mean", "map.mean", "ndvi.mean")]
env_bio = read.csv("data/env_bio.csv", header = TRUE)
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
bi_focal_spp = filter(bbs_occ_code, Aou %in% exp_pres$spAOU)
  
bbs_final_occ_ll = left_join(bi_focal_spp , lat_long, by = "stateroute")
bbs_final_occ_ll = bbs_final_occ_ll[,c("Aou", "stateroute", "occ", "presence", "ALPHA.CODE",
                                       "latitude", "longitude")]
bbs_final_occ_ll$sp_success = 15 * bbs_final_occ_ll$occ
bbs_final_occ_ll$sp_fail = 15 * (1 - bbs_final_occ_ll$occ) 

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
 #   if(levels(as.factor(sdm_input$presence)) > 1){
  glm_occ <- glm(cbind(sp_success, sp_fail) ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio12, family = binomial(link = logit), data = sdm_input)
  glm_pres <- glm(presence ~ elev.mean + ndvi.mean +bio.mean.bio1 + bio.mean.bio12, family = binomial(link = logit), data = sdm_input)
  predocc <- predict(glm_occ,type=c("response"))
  predpr <- predict(glm_pres,type=c("response"))
 sdm_output = cbind(sdm_input,predpr, predocc)
 
 roccurve <- roc(sdm_output$occ ~ sdm_output$predocc)
 auc =  roc(sdm_output$occ ~ sdm_output$predocc)$auc[1]
 
 rocpres <- roc(sdm_output$presence ~ sdm_output$predpr)
 auc_pres =  roc(sdm_output$presence ~ sdm_output$predpr)$auc[1]
 
 auc_df = rbind(auc_df, c(i, auc, auc_pres))
 j = unique(sdm_input$ALPHA.CODE)
 plot = plot(roccurve, main = paste("AUC Curve for ", j, ".csv",   
                                     sep=""))
 # lines(rocpres, col = "red")
#    }
  }
 write.csv(sdm_output, paste("sdm_output_", i, ".csv",   
                            sep=""), row.names = FALSE)
}

dev.off()

setwd("C:/Git/SDMs")
auc_df = data.frame(auc_df)
names(auc_df) = c("AOU", "AUC", "pres_AUC")
# write.csv(auc_df, "Data/auc_df.csv", row.names = FALSE)
test = dplyr::filter(auc_df, AUC > 0.75 & AUC < 1.0)

#### Mapping for spp ####
# the XY coordinates of species data
coords <- CRS("+proj=longlat +datum=WGS84")

# Determine geographic extent of our data using AOU = i
max.lat <- ceiling(max(sdm_input$latitude))
min.lat <- floor(min(sdm_input$latitude))
max.lon <- ceiling(max(sdm_input$longitude))
min.lon <- floor(min(sdm_input$longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
env.data = as.matrix(sdm_input[,c("latitude", "longitude", "elev.mean", "ndvi.mean", "bio.mean.bio1", "bio.mean.bio12")])
env.data = raster(env.data)

data(wrld_simpl)
# Plot the base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

# Add the points for individual observation
points(x = sdm_input$longitude, 
       y = sdm_input$latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)
box()

# Build species distribution model
obs.data = sdm_input[,c("longitude","latitude")]

rbrick = brick(env.data)
plot(rbrick)

# plot GLM occ v pres
r1 = ggplot(auc_df, aes(x = AUC, y = pres_AUC)) +theme_classic()+ theme(axis.title.x=element_text(size=36, vjust = 2),axis.title.y=element_text(size=36, angle=90, vjust = 2)) + xlab(bquote("Occupancy AUC")) + ylab(bquote("Presence AUC"))+
  geom_abline(intercept = 0, slope = 1, col = "black", lwd = 1.5) + geom_point(cex =4, shape=24)+geom_smooth(method='lm', se=FALSE, col="blue",linetype="longdash", lwd =2.5) +
  theme(axis.text.x=element_text(size = 32),axis.ticks=element_blank(), axis.text.y=element_text(size=32))+ scale_colour_manual("", values=c("#dd1c77","#2ca25f","dark gray"))+
  guides(colour = guide_legend(override.aes = list(shape = 15)))+theme(legend.title=element_blank(), legend.text=element_text(size=36), legend.position = c(0.8,0.2), legend.key.width=unit(2, "lines"), legend.key.height =unit(3, "lines")) 
ggsave("Figures/Occ_Pres.pdf", height = 8, width = 12)


