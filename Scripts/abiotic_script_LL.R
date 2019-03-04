library(dismo)
library(rgdal)
library(sp)
library(raster)
# library(maptools)
library(dplyr)

# reading in bioclim, TO, RO, and lat long data for BBS
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5, # 2.5 minutes of a degree
                        path = "data/")

bbs_occ = read.csv("data/bbs_sub1.csv", header=TRUE)
traits = read.csv("data/Master_RO_Correlates.csv", header = TRUE)
lat_long = read.csv("data/latlong_rtes.csv", header = TRUE)
env_bio = read.csv("data/env_bio.csv", header = TRUE)

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

# extract bioclim stack
bioclim.data2 = projectRaster(bioclim.data, crs = prj.string) 

bio.point = raster::extract(bioclim.data, lat_long)
bio.mean = raster::extract(bioclim.data, circs.sp, fun = mean, na.rm=T)
bio.var = raster::extract(bioclim.data, circs.sp, fun = var, na.rm=T)

env_bio = data.frame(stateroute = names(circs.sp), bio.point = bio.point, bio.mean = bio.mean, bio.var = bio.var)

write.csv(env_bio, "data/env_bio.csv", row.names = FALSE)

# evhdf <- data.frame(routes, unlist(evh.point.ext), nrow=5652, byrow=T)