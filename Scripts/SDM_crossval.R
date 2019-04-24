sp_list = unique(bbs_final_occ_ll$Aou)
test_df = c()
for(i in sp_list){
  sdm_output = c()
  print(i)
  bbs_sub <- filter(bbs_final_occ_ll, Aou == i) # %>% filter(occ <= 0.33333333) RUN FOR EXCL TRANS
  bbs_new_sub <- filter(bbs_new, aou == i) 
  bbs_new_sub$pres_2016 <- bbs_new_sub$presence
  temp <- filter(all_env, stateroute %in% bbs_sub$stateroute)
  sdm_input <- left_join(bbs_sub, temp, by = "stateroute")
  sdm_input = na.omit(sdm_input)
  if(length(unique(sdm_input$stateroute)) > 40 & length(unique(sdm_input$presence)) >1){
    if(nrow(filter(sdm_input, presence == 1)) > 49){
     
      gam_occ <- mgcv::gam(cbind(sp_success, sp_fail) ~ s(elev.mean) + s(ndvi.mean) + s(bio.mean.bio4) + s(bio.mean.bio5) + s(bio.mean.bio6) + s(bio.mean.bio13) + s(bio.mean.bio14) , family = binomial(link = logit), data = sdm_input)
   
      pred_gam_occ <- predict(gam_occ,type=c("response"))
     
      sdm_output = cbind(sdm_input, pred_gam_occ) 
      pred_2016 <- left_join(sdm_output, bbs_new_sub[c("stateroute", "pres_2016")], by = "stateroute")
      thresh <- max(pred_2016$pred_gam_occ) * 0.7
      gam_rescale <- filter(pred_2016, pred_gam_occ > thresh) %>%
        dplyr::select(Aou, stateroute,occ, presence, latitude, longitude, pred_gam_occ, pres_2016)
      gam_rescale$rescale <- 1 # discretize the data
      # get absences back in, make a table of the zeroes and ones
      # nest by species
      table(sdm_output$pres_2016, sdm_output$presence)
      test_df = rbind(test_df, gam_rescale)
    }
  }
}

