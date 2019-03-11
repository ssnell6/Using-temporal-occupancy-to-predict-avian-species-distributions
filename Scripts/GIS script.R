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
# read in bird range shps
shapefile_path = 'Z:/GIS/birds/All/All/'
# on mac shapefile_path = '/Volumes/hurlbertlab/GIS/birds/All/All'
all_spp_list = list.files(shapefile_path)

bbs_final_occ_ll = read.csv("Data/final_focal_spp.csv", header = TRUE)
bbs_final_names = left_join(bbs_final_occ_ll, AOU[,c("AOU_OUT", "SCI_NAME")], by = c("Aou" = "AOU_OUT"))
bbs_final_names$focalcat = gsub(" ", "_",bbs_final_names$SCI_NAME)

# test to see which species d/n have matching names
all_spp_list2 = gsub("(.*)_.*", "\\1", all_spp_list) 
all_spp_list2 = unique(all_spp_list2)
all_spp_list2 = data.frame(all_spp_list2)
match = anti_join(bbs_final_names, all_spp_list2, by = c("focalcat" = "all_spp_list2")) 
match = unique(match)

# for loop to select a genus_spp from pairwise table, read in shp, subset to permanent habitat, plot focal distribution
filesoutput = c()
# dropping non-intersecting polygons
focal_spp = unique(bbs_final_names$focalcat)

intl_proj = CRS("+proj=longlat +datum=WGS84")
sp_proj = CRS("+proj=laea +lat_0=40 +lon_0=-100 +units=km")

######## Calculating expected presences for each species using whole range #####
bbs_routes$latitude = abs(bbs_routes$latitude)
bbs_routes = dplyr::select(bbs_routes, -c(statenum, route))
routes = unique(bbs_routes$stateroute)

coordinates(bbs_routes) <- c("longitude", "latitude")
proj4string(bbs_routes) <- intl_proj

expect_pres = c()

file_names = dir('sp_routes/')
setwd("sp_routes/")

for (sp in focal_spp){
  print(sp)
  
  focalAOU = subset(bbs_final_names, focalcat == sp)
  spAOU = unique(focalAOU$Aou)
  
  t1 = all_spp_list[grep(sp, all_spp_list)]
  t2 = t1[grep('.shp', t1)]
  t3 = strsplit(t2, ".shp")
  
  test.poly <- readShapePoly(paste(shapefile_path, t3, ".shp", sep = "")) # reads in species-specific shapefile
  proj4string(test.poly) <- intl_proj
  sporigin = test.poly[test.poly@data$SEASONAL == 1|test.poly@data$SEASONAL == 2|test.poly@data$SEASONAL ==5,]
  
  plot(sporigin)
  # source: http://gis.stackexchange.com/questions/63793/how-to-overlay-a-polygon-over-spatialpointsdataframe-and-preserving-the-spdf-dat
  routes_inside <- bbs_routes[!is.na(sp::over(bbs_routes,as(sporigin,"SpatialPolygons"))),]
  plot(routes_inside, add = T)
  routes_inside = data.frame(routes_inside)
  routes_inside = cbind(routes_inside, spAOU)
  expect_pres=rbind(expect_pres, routes_inside)
}


expect_pres = data.frame(expect_pres)
  
expect_pres = dplyr::select(expect_pres, -optional)
