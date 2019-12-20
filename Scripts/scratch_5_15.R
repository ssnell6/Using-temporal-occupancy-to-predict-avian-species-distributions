# 
# 
# join_pres <- full_join(num_pres_15, num_pres_5, by = "aou") %>%
#   full_join(num_routes, by = "aou") %>%
#   mutate(five_year_prop = n_pres.y/n,
#          fift_year_prop = n_pres.x/n)
# write.csv(join_pres, "data/five_fifteen_diff.csv", row.names = FALSE)

auc_df_15 <- read.csv("Data/auc_df.csv", header = TRUE)
auc_df_5 <- read.csv("Data/auc_df_5.csv", header = TRUE)
join_pres <- read.csv("Data/five_fifteen_diff.csv", header = TRUE)

comp_plot <- full_join(auc_df_15, auc_df_5, by = "AOU") %>%
  full_join(join_pres, by = c("AOU" = "aou")) %>%
  mutate(delta_rmse = rmse_occ.x - rmse_occ.y,
         delta_rpres = fifteen_pres - five_pres)

ggplot(comp_plot, aes(x = delta_rpres, y = delta_rmse)) + geom_point() + geom_hline(yintercept = 0)

ggplot(comp_plot, aes(x = rmse_occ.x, y = rmse_occ.y)) + geom_point() + geom_hline(yintercept = 0)
cor(comp_plot$rmse_occ.x, comp_plot$rmse_occ.y)
