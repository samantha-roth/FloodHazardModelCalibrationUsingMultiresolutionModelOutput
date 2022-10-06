#nCh and RWE only
#50 10m, 200 50m
#downscale 50m to 10m using bilinear interpolation 

rm(list=ls())
library(boot); library(raster)
library(ggplot2); library(viridis)
library(fields); library(akima)

nRuns10m<- 100
nRuns50m<- 400
res.e<-10
res.c<-50

RunTrue.10m= raster("/storage/work/svr5482/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/Extent/RunTrue_1.asc")
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
truevals.10m <- extract(RunTrue.10m,coords.10m)
summary(coords.10m)
miny<- min(coords.10m[,2])
maxy<- max(coords.10m[,2])
minx<- min(coords.10m[,1])
maxx<- max(coords.10m[,1])

#plot(RunTrue.10m)

#load pseudo observations (run 5)
RunTrue.50m= raster("/storage/work/svr5482/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/Extent/RunTrue_1.asc")
#crs(RunTrue.50m)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.50m <- xyFromCell(RunTrue.50m,1:ncell(RunTrue.50m))
truevals.50m<- raster::extract(RunTrue.50m,coords.50m)

min.x50<- min(coords.50m[,1])
min.y50<- min(coords.50m[,2])
max.x50<- max(coords.50m[,1])
max.y50<- max(coords.50m[,2])

xIndsIwant1<- which(coords.10m[,1]>min.x50)
xIndsIwant2<- which(coords.10m[,1]<max.x50)
yIndsIwant1<- which(coords.10m[,2]>min.y50)
yIndsIwant2<- which(coords.10m[,2]<max.y50)

yIndsIwant<- intersect(yIndsIwant1,yIndsIwant2)
xIndsIwant<- intersect(xIndsIwant1,xIndsIwant2)
coordsIwantInds<- intersect(yIndsIwant,xIndsIwant)
coordsIwant.10m<- coords.10m[coordsIwantInds,]
Run.10m.sized<- rasterFromCells(RunTrue.10m,coordsIwantInds)
values(Run.10m.sized)<- truevals.10m[coordsIwantInds]
#plot(Run.10m.smaller)

save(coordsIwant.10m,coordsIwantInds,file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mCoordsin50m.RData")

################################################################################

load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mCoordsin50m.RData")

#what about predicting at the 10m resolution using 50m resolution model output?
#plot(RunTrue.50m)
#plot(RunTrue.10m.sized)

#identical(parVals10m,parVals50m[1:nRuns10m,])

y10m<- unique(coordsIwant.10m[,2])
x10m<- unique(coordsIwant.10m[,1])
y50m<- unique(coords.50m[,2])
x50m<- unique(coords.50m[,1])

nx10m<- length(x10m)
ny10m<- length(y10m)
nx50m<- length(x50m)
ny50m<- length(y50m)

################################################################################
##try bilinear interpolation with akima package

##library(akima)

##Example
##data(akima760)
### interpolate at the diagonal of the grid [0,8]x[0,10]
##akima.bil <- bilinear(akima760$x,akima760$y,akima760$z,
##                      seq(0,8,length=50), seq(0,10,length=50))
##plot(sqrt(akima.bil$x^2+akima.bil$y^2), akima.bil$z, type="l")


##x.inds<- which(coords.10m[,1]>min.x50)
##y.inds<- which(coords.10m[,2]<max.y50)
##coordstointerpolate<- coords.10m[intersect(x.inds,y.inds),]

#z1= matrix(truevals.50m, nrow= ny50m, ncol= nx50m, byrow= TRUE)
##need to switch x and y axes
##order of true vals: max y in coords (x in matrix), min x in coords (y in matrix) 
##to min y in coords (x in matrix), max x in coords (y in matrix) 

#z= matrix(NA,nrow= ny50m, ncol= nx50m)
#for(i in 1:nx50m){
#  z[,i]<- rev(z1[,i])
#}

##now the coordinates are ordered min x, min y to max x max y for bilinear interpolation

#test<- bilinear(x= rev(y50m), y= x50m, z= z, 
#                x0= rev(coordsIwant.10m[,2]), y0= coordsIwant.10m[,1])

#z2= matrix(test$z, nrow= ny10m, ncol= nx10m, byrow= TRUE)

#translate back to original coordinates
#downscaled.z<- matrix(NA,nrow= ny10m, ncol= nx10m)
#for(i in 1:nx10m){
#  downscaled.z[,i]<- rev(z2[,i])
#}

#downscaled.z.vec<- c(t(downscaled.z))

#now we're back to max y, min x to min y, max x (top left to bottom right)

#RunTrue.50m.downscaled<- Run.10m.sized
#values(RunTrue.50m.downscaled)<- downscaled.z.vec
#plot(RunTrue.50m.downscaled)
#plot(RunTrue.50m)


################################################################################

res.e<-10
res.c<-50
nLoc10m<- 126791
nLoc10min50m<- 125050
nLoc50m<- 5146


if (dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c") == F){
  dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c")}

if (dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/downscale") == F){
  dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/downscale")}

if(dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c") ==F){
  dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c")}

if(dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale") ==F){
  dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale")}

#create matrix of interpolated 50m output

Runs50mDownscaled<- matrix(NA, nrow= nRuns50m, ncol= nLoc10min50m)
st1 <- Sys.time()
for(i in 1:nRuns50m){
  run50m<- raster(paste("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs50m/nCh_RWE/100e400c/Extent/Run_",i,".asc",sep=""))
  vals.50m<- extract(run50m,coords.50m)
  z1= matrix(vals.50m, nrow= ny50m, ncol= nx50m, byrow= TRUE)
  z= matrix(NA,nrow= ny50m, ncol= nx50m)
  for(j in 1:nx50m){ z[,j]<- rev(z1[,j]) }
  test<- bilinear(x= rev(y50m), y= x50m, z= z, 
                  x0= rev(coordsIwant.10m[,2]), y0= coordsIwant.10m[,1])
  z2= matrix(test$z, nrow= ny10m, ncol= nx10m, byrow= TRUE)
  #translate back to original coordinates
  downscaled.z<- matrix(NA,nrow= ny10m, ncol= nx10m)
  for(k in 1:nx10m){downscaled.z[,k]<- rev(z2[,k])}
  downscaled.z.vec<- c(t(downscaled.z))
  Runs50mDownscaled[i,]<- downscaled.z.vec
}

save(Runs50mDownscaled,file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/Runs50mDownscaled.RData")
en1 <- Sys.time()
time= en1-st1
save(time,file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/time_downscale.RData")