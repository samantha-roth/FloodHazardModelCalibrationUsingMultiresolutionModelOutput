#Selinsgrove specific
rm(list=ls())
#plot Lisflood calibration preds
library(raster)
library(RColorBrewer)
library(ggmap)
library(terra)
library(rgdal)

#coords in 50m bounds
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mCoordsin50m.RData")

#compare to the observed predictions
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/noisytruth.RData")


RunTrue.10m= raster("/storage/work/svr5482/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/Extent/RunTrue_1.asc")
#set crs
crs(RunTrue.10m)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
values(RunTrue.10m)[coordsIwantInds]<- noisy2.3truevals.10m
truevals.10m <- raster::extract(RunTrue.10m,coords.10m)

nPCs=18
load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/meanCalPredMR.RData",sep=""))
nPCs=8
load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/cal/",nPCs,"PCs/.25/3cm/meanCalPredHR.RData",sep=""))  


Run_homML_mean<- RunTrue.10m
values(Run_homML_mean)<- meanCalPredMR


Run_homGP_mean<- RunTrue.10m
values(Run_homGP_mean)<- meanCalPredHR


#GET MAP OF SELINSGROVE

max.y.10m<- max(coords.10m[,2])
min.y.10m<- min(coords.10m[,2])
max.x.10m<- max(coords.10m[,1])
min.x.10m<- min(coords.10m[,1])

#Produce map of Selinsgrove
points<- coords.10m
v <- vect(points, crs="+proj=utm +zone=18 +datum=WGS84  +units=m")
v
y <- terra::project(v, "+proj=longlat +datum=WGS84")
y

lonlat <- geom(y)[, c("x", "y")]
head(lonlat, 3)

max.lat.10m<- max(lonlat[,2])
min.lat.10m<- min(lonlat[,2])
max.lon.10m<- max(lonlat[,1])
min.lon.10m<- min(lonlat[,1])


selingsgrove_box <- c(
  left = min.lon.10m,
  bottom = min.lat.10m,
  right = max.lon.10m,
  top = max.lat.10m
)

selinsgrove_stamen <- get_stamenmap(
  bbox = selingsgrove_box,
  maptype = "toner",
  zoom = 15
)

selinsgrove_ggmap<- ggmap(selinsgrove_stamen)

################################################################################
################################################################################
#good method,  now I just need to restrcit raster layer
#r <- raster(system.file("external/test.grd", package="raster"))

#plot true vals >0 

# just to make it reproducible with ggmap we have to transform to wgs84
RunTrue.10mG0<- RunTrue.10m
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
vals.10m<- raster::extract(RunTrue.10m,coords.10m)
ZeroInds<- which(vals.10m==0)
vals.10mG0<- vals.10m
vals.10mG0[ZeroInds]<- NA
values(RunTrue.10mG0)<- vals.10mG0
vals.10mincm<- as.integer(round(100*vals.10mG0))
values(RunTrue.10mG0)<- vals.10mincm #convert to cm since rtp can only deal with integer values

r <- projectRaster(RunTrue.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

dataType(r)
dataType(r)="INT4S"

rtp <- rasterToPolygons(r)

selinsgrove_ggmap + 
  geom_polygon(data = rtp, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.5)  + 
  scale_fill_gradientn("Flood Height (cm)", colors = heat.colors(500)) 



################################################################################
#HomML 

crs(Run_homML_mean)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homML_mean.10mG0<- Run_homML_mean
Run_homML_meanvals.10m<- raster::extract(Run_homML_mean,coords.10m)
ZeroInds<- which(Run_homML_meanvals.10m==0)
Run_homML_meanvals.10mG0<- Run_homML_meanvals.10m
Run_homML_meanvals.10mG0[ZeroInds]<- NA
values(Run_homML_mean.10mG0)<- Run_homML_meanvals.10mG0
Run_homML_meanvals.10mincm<- as.integer(round(100*Run_homML_meanvals.10mG0))
values(Run_homML_mean.10mG0)<- Run_homML_meanvals.10mincm #convert to cm since rtp can only deal with integer values

r.homML <- projectRaster(Run_homML_mean.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homML)="INT4S"

rtp.homML <- rasterToPolygons(r.homML)


################################################################################
#HomGP10

crs(Run_homGP_mean)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homGP_mean.10mG0<- Run_homGP_mean
Run_homGP_meanvals.10m<- raster::extract(Run_homGP_mean,coords.10m)
ZeroInds<- which(Run_homGP_meanvals.10m==0)
Run_homGP_meanvals.10mG0<- Run_homGP_meanvals.10m
Run_homGP_meanvals.10mG0[ZeroInds]<- NA
values(Run_homGP_mean.10mG0)<- Run_homGP_meanvals.10mG0
Run_homGP_meanvals.10mincm<- as.integer(round(100*Run_homGP_meanvals.10mG0))
values(Run_homGP_mean.10mG0)<- Run_homGP_meanvals.10mincm #convert to cm since rtp can only deal with integer values

r.homGP <- projectRaster(Run_homGP_mean.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homGP)="INT4S"

rtp.homGP <- rasterToPolygons(r.homGP)

################################################################################

max(rtp$RunTrue_1)
max(rtp.homML$RunTrue_1)
max(rtp.homGP$RunTrue_1)
max.m<-round(max(c(max(rtp$RunTrue_1),max(rtp.homML$RunTrue_1),max(rtp.homGP$RunTrue_1)))/100)

rtp$RunTrue_1<- rtp$RunTrue_1/100
rtp.homML$RunTrue_1<- rtp.homML$RunTrue_1/100
rtp.homGP$RunTrue_1<- rtp.homGP$RunTrue_1/100
################################################################################

cuts=seq(0,max.m,by=.01)
pal <- colorRampPalette(c("white","blue"))
plot(RunTrue.10m, breaks=cuts, col = pal(length(cuts)+1)) #plot with defined breaks

#par(mfrow=c(3,2))

jpeg(filename = "/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/plots/Observation1.jpeg",
     width = 800, height = 400) #height was 500
selinsgrove_ggmap + 
  geom_polygon(data = rtp, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Flood Height (m)", colors = pal(max.m*100), limits= c(0,max.m)) +
  labs(title= expression(theta^"*1"*": simulated observation")) +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()

################################################################################

jpeg(filename = "/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/plots/HomMR_meanpred.jpeg",
     width = 800, height = 400) #height was 500
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homML, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homML$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  scale_fill_gradientn("Flood Height (m)", colors = pal(max.m*100), limits= c(0,max.m)) +
  labs(title= expression(theta^"*1"*": MR calibrated projection"))+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24))
dev.off()


################################################################################

jpeg(filename = "/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/plots/HomGP10_meanpred.jpeg",
     width = 800, height = 400) #height was 500
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homGP, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homGP$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  scale_fill_gradientn("Flood Height (m)", colors = pal(max.m*100), limits= c(0,max.m)) +
  labs(title= expression(theta^"*1"*": HR calibrated projection"))+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24))
dev.off()


################################################################################
#plot differences in mean predictions

################################################################################

#now do this for the mean predictions- truth

Run_homML_diff<- RunTrue.10m
values(Run_homML_diff)<- meanCalPredMR - truevals.10m

Run_homGP_diff<- RunTrue.10m
values(Run_homGP_diff)<- meanCalPredHR - truevals.10m


################################################################################
#HomML 

crs(Run_homML_diff)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homML_diff.10mG0<- Run_homML_diff
Run_homML_diffvals.10m<- raster::extract(Run_homML_diff,coords.10m)
ZeroInds<- which(Run_homML_diffvals.10m==0)
Run_homML_diffvals.10mG0<- Run_homML_diffvals.10m
Run_homML_diffvals.10mG0[ZeroInds]<- NA
values(Run_homML_diff.10mG0)<- Run_homML_diffvals.10mG0
Run_homML_diffvals.10mincm<- as.integer(round(100*Run_homML_diffvals.10mG0))
values(Run_homML_diff.10mG0)<- Run_homML_diffvals.10mincm #convert to cm since rtp can only deal with integer values

r.homMLdiff <- projectRaster(Run_homML_diff.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homMLdiff)="INT4S"

rtp.homMLdiff <- rasterToPolygons(r.homMLdiff)
rtp.homMLdiff$RunTrue_1<- rtp.homMLdiff$RunTrue_1/100

################################################################################
#HomGP10

crs(Run_homGP_diff)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
Run_homGP_diff.10mG0<- Run_homGP_diff
Run_homGP_diffvals.10m<- raster::extract(Run_homGP_diff,coords.10m)
ZeroInds<- which(Run_homGP_diffvals.10m==0)
Run_homGP_diffvals.10mG0<- Run_homGP_diffvals.10m
Run_homGP_diffvals.10mG0[ZeroInds]<- NA
values(Run_homGP_diff.10mG0)<- Run_homGP_diffvals.10mG0
Run_homGP_diffvals.10mincm<- as.integer(round(100*Run_homGP_diffvals.10mG0))
values(Run_homGP_diff.10mG0)<- Run_homGP_diffvals.10mincm #convert to cm since rtp can only deal with integer values

r.homGPdiff <- projectRaster(Run_homGP_diff.10mG0, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
dataType(r.homGPdiff)="INT4S"

rtp.homGPdiff <- rasterToPolygons(r.homGPdiff)
rtp.homGPdiff$RunTrue_1<- rtp.homGPdiff$RunTrue_1/100
################################################################################

max.resid<- max(c(max( meanCalPredMR - truevals.10m),max(meanCalPredHR - truevals.10m)))

min.resid<- min(c(min( meanCalPredMR - truevals.10m),min( meanCalPredHR - truevals.10m)))

#par(mfrow=c(3,2))

#maxmax<- round(max(c(max(rtp.homMLdiff$RunTrue_1),max(rtp.homGPdiff$RunTrue_1))),1)
#minmin<- round(min(c(min(rtp.homMLdiff$RunTrue_1),min(rtp.homGPdiff$RunTrue_1))),1)

absmax<-round(max(c(abs(min.resid),max.resid)),2)

################################################################################
jpeg(filename = "/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/plots/HomML_diffpred2011.jpeg",
     width = 800, height = 400) #height was 600
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homMLdiff, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homMLdiff$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Difference (m)", colors = terrain.colors(100*(2*absmax)), limits= c(-absmax,absmax)) +
  labs(title= expression(theta^"*1"*": MR calibrated projection - observation"))+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))

dev.off()


################################################################################

################################################################################
jpeg(filename = "/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/plots/HomGP10_diffpred2011.jpeg",
     width = 800, height = 400) #height was 600
selinsgrove_ggmap + 
  geom_polygon(data = rtp.homGPdiff, 
               aes(x = long, y = lat, group = group, 
                   fill = rep(rtp.homGPdiff$RunTrue_1, each= 5)), 
               #fill=layer),
               size = 0, 
               alpha = 0.7)  + 
  #scale_fill_gradientn("Flood Height (cm)", colors = rev(heat.colors(800)), limits= c(0,900)) +
  scale_fill_gradientn("Difference (m)", colors = terrain.colors(100*(2*absmax)), limits= c(-absmax,absmax)) +
  labs(title= expression(theta^"*1"*": HR calibrated projection - observation")) +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()

################################################################################

