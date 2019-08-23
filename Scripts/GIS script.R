library(maps)
library(rgdal)
library(shapefiles)
library(maptools)
library(raster)
library(rgeos)
library(gtools)
library(sp)
library(tidyr)
library(dplyr)

# read in lat long data
bbs_routes = read.csv("Data/latlongs.csv",header =TRUE)
AOU = read.csv("Data/Bird_Taxonomy.csv", header = TRUE) # taxonomy data
auc_df <- read.csv("Data/auc_df.csv", header = TRUE)
# read in bird range shps
shapefile_path = 'Z:/GIS/birds/All/All/'
# on mac shapefile_path = '/Volumes/hurlbertlab/GIS/birds/All/All'
all_spp_list = list.files(shapefile_path)

bbs_final_occ_ll = read.csv("Data/bbs_final_occ_ll.csv", header = TRUE) # %>% filter(Year == 2016)
  # read.csv("Data/final_focal_spp.csv", header = TRUE)
bbs_final_names.1 = left_join(bbs_final_occ_ll, AOU[,c("AOU_OUT", "CRC_SCI_NAME")], by = c("aou" = "AOU_OUT"))
bbs_final_names.1$focalcat = gsub(" ", "_",bbs_final_names.1$CRC_SCI_NAME)
bbs_final_names.2 = bbs_final_names.1[-grep("_spuh", bbs_final_names.1$focalcat),] 
bbs_final_names.3 = bbs_final_names.2[-grep("/", bbs_final_names.2$focalcat),] 
# bbs_final_names.4 = bbs_final_names.3[-grep("", bbs_final_names.3$focalcat),] 

bbs_final_names = bbs_final_names.3 %>% filter(aou %in% auc_df$AOU) %>%
  unique() 

# data cleaning 
bbs_final_names$focalcat[bbs_final_names$focalcat =="Oreothlypis_peregrina"] = "Vermivora_peregrina"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Vermivora_pinus"] = "Vermivora_cyanoptera"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Stellula_calliope"] = "Selasphorus_calliope"

bbs_final_names$focalcat = gsub('Setophaga_', 'Dendroica_', bbs_final_names$focalcat)

bbs_final_names$focalcat[bbs_final_names$focalcat =="Dendroica_ruticilla"] = "Setophaga_ruticilla"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Picoides_nuttallii"] = "Dryobates_nuttallii"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Cardellina_canadensis"] = "Wilsonia_canadensis"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Geothlypis_philadelphia"] = "Oporornis_philadelphia"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Oreothlypis_ruficapilla"] = "Vermivora_ruficapilla"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Oreothlypis_celata"] = "Vermivora_celata"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Cardellina_pusilla"] = "Wilsonia_pusilla"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Oreothlypis_virginiae"] = "Vermivora_virginiae"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Poecile_hudsonica"] = "Parus_hudsonicus"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Pica_hudsonia"] = "Pica_pica"

bbs_final_names$focalcat = gsub('Poecile_', 'Parus_', bbs_final_names$focalcat)

bbs_final_names$focalcat[bbs_final_names$focalcat =="Dendroica_citrina"] = "Wilsonia_citrina"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Geothlypis_formosus"] = "Oporornis_formosus"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Oreothlypis_luciae"] = "Vermivora_luciae"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Geothlypis_tolmiei"] = "Oporornis_tolmiei"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Calcarius_mccownii"] = "Rhynchophanes_mccownii"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Picoides_villosus"] = "Leuconotopicus_villosus"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Picoides_pubescens"] = "Dryobates_pubescens"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Picoides_dorsalis"] = "Picoides_tridactylus"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Picoides_scalaris"] = "Dryobates_scalaris"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Picoides_albolarvatus"] = "Leuconotopicus_albolarvatus"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Picoides_borealis"] = "Leuconotopicus_borealis"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Aimophila_cassinii"] = "Peucaea_cassinii"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Aimophila_aestivalis"] = "Peucaea_aestivalis"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Aimophila_botterii"] = "Peucaea_botterii"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Aimophila_carpalis"] = "Peucaea_carpalis"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Oreothlypis_crissalis"] = "Vermivora_crissalis"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Ixoreus_naevius"] = "Zoothera_naevia"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Pipilo_fuscus"] = "Melozone_fuscus"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Pipilo_crissalis"] = "Melozone_crissalis"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Pipilo_aberti"] = "Melozone_aberti"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Seiurus_noveboracensis"] = "Parkesia_noveboracensis"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Seiurus_motacilla"] = "Parkesia_motacilla"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Picoides_arizonae"] = "Leuconotopicus_arizonae"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Dryocopus_pileatus"] = "Hylatomus_pileatus"

bbs_final_names$focalcat[bbs_final_names$focalcat =="Carduelis_hornemanni"] = "Carduelis_flammea"

# test to see which species d/n have matching names
all_spp_list2 = gsub("(.*)_.*", "\\1", all_spp_list) 
all_spp_list2 = unique(all_spp_list2)
all_spp_list2 = data.frame(all_spp_list2)
match = anti_join(bbs_final_names, all_spp_list2, by = c("focalcat" = "all_spp_list2")) 
match = unique(match$focalcat)

# for loop to select a genus_spp from pairwise table, read in shp, subset to permanent habitat, plot focal distribution
filesoutput = c()

intl_proj = CRS("+proj=longlat +datum=WGS84")
sp_proj = CRS("+proj=laea +lat_0=40 +lon_0=-100 +units=km")

######## Calculating expected presences for each species using whole range #####
bbs_routes$latitude = abs(bbs_routes$latitude)
bbs_routes = dplyr::select(bbs_routes, -c(statenum, route))
routes = unique(bbs_routes$stateroute)

coordinates(bbs_routes) <- c("longitude", "latitude")
proj4string(bbs_routes) <- intl_proj

# dropping non-intersecting polygons
focal_spp1 = filter(bbs_final_names, !focalcat %in% match) %>%  
  filter(., focalcat != "Columba_livia") %>%  
  filter(., focalcat != "Circus_cyaneus") %>%  
  filter(., focalcat != "Buteo_nitidus") %>%  
  filter(., focalcat != "Trogon_elegans")  %>%  
  filter(., focalcat != "Aphelocoma_ultramarina") %>%
  filter(., focalcat != "Phylloscopus_borealis") %>%
  filter(., focalcat != "Streptopelia_decaocto") %>%
  filter(., focalcat != "Falco_rusticolus")

focal_spp <- focal_spp1 %>% group_by(stateroute, aou) %>%
  filter(., unique(stateroute) > 40)

sp_list = unique(focal_spp$focalcat)

expect_pres = c()

file_names = dir('sp_routes/')

for (sp in sp_list){
  print(sp)
  
  focalAOU = subset(bbs_final_names, focalcat == sp)
  spAOU = unique(focalAOU$aou)
  
  t1 = all_spp_list[grep(sp, all_spp_list)]
  t2 = t1[grep('.shp', t1)]
  t3 = strsplit(t2, ".shp")
  
  test.poly <- readShapePoly(paste(shapefile_path, t3, ".shp", sep = "")) # reads in species-specific shapefile
  proj4string(test.poly) <- intl_proj
  sporigin = test.poly[test.poly@data$SEASONAL == 1|test.poly@data$SEASONAL == 2|test.poly@data$SEASONAL ==5,]
  
  plot(sporigin)

  routes_inside <- bbs_routes[!is.na(sp::over(bbs_routes,as(sporigin,"SpatialPolygons"))),]
  plot(routes_inside, add = T)
  routes_inside = data.frame(routes_inside)
  routes_inside = cbind(routes_inside, spAOU)
  expect_pres=rbind(expect_pres, routes_inside)
}


expect_pres = data.frame(expect_pres)
expect_pres = dplyr::select(expect_pres, -optional)
# write.csv(expect_pres, "Data/expect_pres_2016.csv", row.names = FALSE)
# write.csv(expect_pres, "Data/expect_pres.csv", row.names = FALSE)