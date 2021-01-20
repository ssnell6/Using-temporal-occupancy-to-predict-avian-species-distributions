### Supplemental information: data table and route map

library(tidyverse)
library(tmap)
library(sf)

### Data table

# Read in data
spec_info <- read.csv("Data/master_RO_correlates.csv", stringsAsFactors = F)
bird_taxo <- read.csv("Data/Bird_Taxonomy.csv", stringsAsFactors = F)

auc_df <- read.csv("Data/auc_df.csv", stringsAsFactors = F)
auc_df_5 <- read.csv("Data/auc_df_5.csv", stringsAsFactors = F)
pres_matrix <- read.csv("Data/pres_matrix.csv", stringsAsFactors = F)
pres_spat <- read.csv("Data/pres_spat.csv", stringsAsFactors = F)
bbs_final_occ_ll <- read.csv("Data/bbs_final_occ_ll.csv", stringsAsFactors = F)

# n routes per species

num_routes <- bbs_final_occ_ll %>% 
  group_by(aou) %>% 
  dplyr::summarise(n_routes = n_distinct(stateroute)) 

# Join to make results table w/ spec, traits, model output

model_out <- auc_df %>%
  right_join(auc_df_5, by= c("AOU"), suffix = c("_15yr", "_5yr")) %>%
  left_join(select(pres_matrix, c(aou,sensitivity_glmocc:np_max)), by = c("AOU"="aou")) %>%
  left_join(select(pres_spat, c(aou,sensitivity_glmocc:np_max)), by= c("AOU"="aou"))

spec_trait <- spec_info %>%
  left_join(bird_taxo, by = c("AOU" = "AOU_OUT")) %>%
  select(AOU, CommonName, CRC_SCI_NAME, migclass) %>%
  distinct()

res_table <- spec_trait %>%
  left_join(num_routes, by = c("AOU" = "aou")) %>%
  right_join(model_out)
write.csv(res_table, "Data/species_model_output_table.csv", row.names = F)

## Map of BBS routes used in analysis

bbs <- read.csv("Data/bbs_2001_2015.csv", stringsAsFactors = F)
bbs_latlon <- read.csv("Data/bbs_routes_20170712.csv", stringsAsFactors = F) %>%
  mutate(stateroute = statenum*1000 + route)

bbs_sf <- bbs %>%
  select(stateroute) %>%
  distinct() %>%
  left_join(bbs_latlon) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) 

na_map <- read_sf("Data/maps/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 == "USA" | sr_adm0_a3 == "CAN", iso_3166_2 != "US-HI")

bbs_map <- tm_shape(na_map) + tm_polygons(col = "gray") + 
  tm_shape(bbs_sf) + tm_dots(col = "black", size = 0.05) + tm_layout("BBS routes",
                        legend.title.size = 1)


tmap_save(bbs_map, "Figures/bbs_route_map.pdf", units = "in", height = 5, width = 8)

