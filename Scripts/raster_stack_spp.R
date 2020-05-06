## going to have to run on LL
# get each spp range 
# raster stack all the env layers
# output as raster

library(raster)
library(tidyverse)
library(rgdal)
library(sp)
library(maptools)
library(gimms)

# reading in bioclim, TO, RO, and lat long data for BBS
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5, # 2.5 minutes of a degree
                        path = "data/")

crop_extent <- c(-180, 180, -60, 90)
bioclim.sub <- subset(bioclim.data, c("bio4", "bio5", "bio6", "bio13", "bio14"))

elev <- raster("Z:/GIS/DEM/sdat_10003_1_20170424_102000103.tif")

gimms_files <- list.files("\\\\BioArk\\HurlbertLab\\GIS\\gimms\\")
gimms_df <- data.frame(file_name = gimms_files[-1], year = as.numeric(substr(gimms_files[-1], 15, 18)))
files <- filter(gimms_df, year %in% 2000:2014)
setwd("\\\\BioArk\\HurlbertLab\\GIS\\gimms\\")
gimms_jan <- rasterizeGimms(as.character(files$file_name)[1])
gimms_jul <- rasterizeGimms(as.character(files$file_name)[2])
gimms_breeding <- stack(c(gimms_jan[[11:12]], gimms_jul[[1:4]]))
gimms_reclass <- reclassify(gimms_breeding, cbind(0, NA))
gimms_mean <- mean(gimms_reclass, na.rm=FALSE)
# gimms_crop <- crop(gimms_mean, extent(elev))
elev_agg <- aggregate(elev, fact=10, fun=mean)

gimms_crop <- raster::crop(gimms_mean, extent(elev_agg))
gimms_reproj <- raster(vals=values(gimms_crop),ext=extent(elev_agg),crs=crs(elev_agg),
              nrows=dim(elev_agg)[1],ncols=dim(elev_agg)[2])

bio_crop <- crop(bioclim.sub, extent(elev_agg))
# not working here
bio_resample <- raster::resample(bio_crop, elev_agg, method = "nbg")
gimms_reproj <- raster(vals=values(bio_crop),ext=extent(elev_agg),crs=crs(elev_agg),
                       nrows=dim(elev_agg)[1],ncols=dim(elev_agg)[2])
elev_ndvi <- stack(gimms_reproj, elev_agg)

setwd("C:/Git/SDMs")
# read in lat long data
auc_df <- read.csv("Data/auc_df.csv", header = TRUE)
bbs_routes = read.csv("Data/latlongs.csv",header =TRUE)
AOU = read.csv("Data/Bird_Taxonomy.csv", header = TRUE) %>% 
  filter(AOU_OUT %in% auc_df$AOU)
bbs_final_names = read.csv("Data/bbs_final_occ_ll.csv", header = TRUE) %>% 
  filter(aou %in% auc_df$AOU) %>%
  left_join(AOU[,c("AOU_OUT", "CRC_SCI_NAME")], by = c("aou" = "AOU_OUT")) %>%
  mutate(focalcat = gsub(" ", "_", CRC_SCI_NAME))

# read in bird range shps
shapefile_path = "Z:/GIS/birds/All/All/"
# '/proj/hurlbertlab/ssnell/bird_range_shps'
# on mac shapefile_path = '/Volumes/hurlbertlab/GIS/birds/All/All'
all_spp_list = list.files(shapefile_path)

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
# write.csv(focal_spp, "Data/spp_ranges_for_ll.csv", row.names = F)
sp_list = unique(focal_spp$focalcat)

intl_proj = CRS("+proj=longlat +datum=WGS84")
sp_proj = CRS("+proj=laea +lat_0=40 +lon_0=-100 +units=km")

expect_pres = c()

file_names = dir('sp_routes/')

for (sp in sp_list){
  print(sp)
  
  focalAOU = subset(bbs_final_names, focalcat == sp)
  spAOU = unique(focalAOU$aou)
  
  t1 = all_spp_list[grep(sp, all_spp_list)]
  t2 = t1[grep('.shp', t1)]
  t3 = strsplit(t2, ".shp")
  test.poly <- st_read(paste(shapefile_path, t3, ".shp", sep = "")) %>%
    filter(SEASONAL %in% c(1,2,5)) %>%
    select(SCINAME, geometry)
  sporigin <- st_transform(test.poly, sp_proj) 
  
  ggplot() +
    geom_sf(data = sporigin) +
    # env vars to stack and extract
}