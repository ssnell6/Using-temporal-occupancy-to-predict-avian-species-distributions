library(dismo)
library(rgdal)
library(sp)
library(raster)
library(maptools)
library(dplyr)
library(gimms)

# reading in bioclim, TO, RO, and lat long data for BBS
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5, # 2.5 minutes of a degree
                        path = "data/")

bbs_occ = read.csv("data/bbs_sub1.csv", header=TRUE)
traits = read.csv("data/Master_RO_Correlates.csv", header = TRUE)
lat_long = read.csv("data/latlong_rtes.csv", header = TRUE)

bbs_rtes = left_join(bbs_occ, lat_long[,c("latitude", "longitude", "stateroute")], by = "stateroute")

# From BI

# Define projection to be used throughout analysis
prj.string <- "+proj=laea +lat_0=45.235 +lon_0=-106.675 +units=km"
# Makes routes into a spatialPointsDataframe
coordinates(lat_long)=c('longitude','latitude')
projection(lat_long) = CRS("+proj=longlat +ellps=WGS84")

# bc.model <- bioclim(x = bioclim.data, p = bbs_rtes)

# Transforms routes to an equal-area projection - see previously defined prj.string
routes.laea = spTransform(lat_long, CRS(prj.string))

# A function that draws a circle of radius r around a point: p (x,y)
RADIUS = 40

make.cir = function(p,r){
  points=c()
  for(i in 1:360){
    theta = i*2*pi/360
    y = p[2] + r*cos(theta)
    x = p[1] + r*sin(theta)
    points = rbind(points,c(x,y))
  }
  points=rbind(points,points[1,])
  circle=Polygon(points,hole=F)
  circle
}

#Draw circles around all routes
circs = sapply(1:nrow(routes.laea), function(x){
  circ =  make.cir(routes.laea@coords[x,],RADIUS)
  circ = Polygons(list(circ),ID=routes.laea@data$stateroute[x])
}
)
circs.sp = SpatialPolygons(circs, proj4string=CRS(prj.string))

# Check that circle locations look right
plot(circs.sp)

# extract bioclim stack
bioclim.data2 = projectRaster(bioclim.data, crs = prj.string) 

bio.point = raster::extract(bioclim.data, routes)
bio.mean = raster::extract(bioclim.data, circs.sp, fun = mean, na.rm=T)
bio.var = raster::extract(bioclim.data, circs.sp, fun = var, na.rm=T)

env_bio = data.frame(stateroute = names(circs.sp), bio.point = bio.point, bio.mean = bio.mean, bio.var = bio.var)
write.csv(env_bio, "data/env_bio.csv", row.names = FALSE)



prj.string = "+proj=longlat +datum=WGS84"
#### for maps ####
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

gimms_crop <- crop(gimms_mean, elev)
bio_crop <- crop(bioclim.sub, elev)

# RADIUS = 0.1
# mod.p$rowname <- 1:660
# #Draw circles around all routes
# circs = sapply(1:nrow(mod.p), function(x){
#   circ =  make.cir(mod.p@coords[x,],RADIUS)
#   circ = Polygons(list(circ),ID=mod.p$rowname[x])
# }
# )
# circs.sp = SpatialPolygons(circs, proj4string=CRS(prj.string))


elev_points <- raster::extract(elev, mod.p, na.rm=T)
gimms_points <- raster::extract(gimms_crop, mod.p, na.rm=T)
bio_points <- raster::extract(bio_crop, mod.p, na.rm=T)
env_pred = data.frame(lat = mod.p$latitude, long = mod.p$longitude, elev_point = elev_points, bio_point = bio_points, gimms_point = gimms_points)
write.csv(env_pred, "Data/predicted_env_lat_longs.csv", row.names = FALSE)

r = raster(nrows = 22, ncols = 30, geographic.extent, 1) 
projection(r) <- "+proj=longlat +datum=WGS84"
p = rasterToPoints(r)
p = data.frame(p)
names(p) = c("longitude", "latitude")
mod.p <- SpatialPointsDataFrame(coords = p, data = p, proj4string = CRS("+proj=longlat +datum=WGS84"))


bioclim.data

