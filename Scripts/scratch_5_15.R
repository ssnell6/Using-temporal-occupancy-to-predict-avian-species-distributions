# 
# 
# join_pres <- full_join(num_pres_15, num_pres_5, by = "aou") %>%
#   full_join(num_routes, by = "aou") %>%
#   mutate(five_year_prop = n_pres.y/n,
#          fift_year_prop = n_pres.x/n)
# write.csv(join_pres, "data/five_fifteen_diff.csv", row.names = FALSE)

auc_df_15 <- read.csv("Data/auc_df.csv", header = TRUE) %>%
  mutate(glm_diff = rmse_occ - rmse_pres,
         gam_diff = rmse_gam - rmse_gam_pres,
         rf_diff = rmse_rf - rmse_rf_pres) %>%
  filter(gam_diff > 0)
auc_df_5 <- read.csv("Data/auc_df_5.csv", header = TRUE) %>%
  mutate(glm_diff_5 = rmse_occ - rmse_pres,
         gam_diff_5 = rmse_gam - rmse_gam_pres,
         rf_diff_5 = rmse_rf - rmse_rf_pres)
comp_plot <- left_join(auc_df_15, auc_df_5, by = "AOU") %>%
  left_join(tax_code, by = c("AOU" = "AOU_OUT")) %>%
  mutate(delta_glm_rmse = glm_diff - glm_diff_5,
         delta_gam_rmse = gam_diff - gam_diff_5,
         delta_rf_rmse = rf_diff - rf_diff_5) 

ggplot(comp_plot, aes(x = delta_rpres, y = delta_rmse)) + geom_text(aes(label = ALPHA.CODE)) + geom_hline(yintercept = 0) 

ggplot(comp_plot, aes(x = rmse_occ.x, y = rmse_occ.y)) + geom_point() + geom_hline(yintercept = 0)
cor(comp_plot$rmse_occ.x, comp_plot$rmse_occ.y)
