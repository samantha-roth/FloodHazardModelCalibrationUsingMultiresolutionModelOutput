# compile calibrated predictions


rm(list=ls())
library(raster)

nLoc10m<- 126791

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs10m/nCh_RWE/100e400c/MRMCMC100/Extent")

sumCalPredMR<- rep(0,nLoc10m)

nRuns10m=100

for(i in 1:nRuns10m){
  print(i)
  Run10m<- raster(paste("Run_",i,".asc",sep=""))
  #Run10m<- raster(paste("Run_",i,".asc",sep=""))
  coords10m <- xyFromCell(Run10m,1:ncell(Run10m))
  vals10m <- raster::extract(Run10m,coords10m)
  sumCalPredMR<- sumCalPredMR + vals10m
} 
meanCalPredMR<- sumCalPredMR/nRuns10m
print("MR done")

nPCs=18
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
save(meanCalPredMR,file="meanCalPredMR.RData")

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs10m/nCh_RWE/100e400c/HRMCMC100/Extent")

sumCalPredHR<- rep(0,nLoc10m)

nRuns10m=100

for(i in 1:nRuns10m){
  print(i)
  Run10m<- raster(paste("Run_",i,".asc",sep=""))
  #Run10m<- raster(paste("Run_",i,".asc",sep=""))
  coords10m <- xyFromCell(Run10m,1:ncell(Run10m))
  vals10m <- raster::extract(Run10m,coords10m)
  sumCalPredHR<- sumCalPredHR + vals10m
} 

meanCalPredHR<- sumCalPredHR/nRuns10m
print("HR done")

nPCs=8
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/cal/",nPCs,"PCs/.25/3cm",sep=""))  
save(meanCalPredHR,file="meanCalPredHR.RData")

################################################################################
nPCs=18
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("meanCalPredMR.RData")

nPCs=8
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/cal/",nPCs,"PCs/.25/3cm",sep=""))  
load("meanCalPredHR.RData")

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


################################################################################
#Calculate RMSE

homML_mean<- meanCalPredMR
homGP_mean<- meanCalPredHR

homML_RMSE<- sqrt(mean((homML_mean - truevals.10m)^2))

homGP_RMSE<- sqrt(mean((homGP_mean - truevals.10m)^2))

homML_RMSE; homGP_RMSE

#PBIAS

homML_PBias<- 100 * (sum( homML_mean - truevals.10m ) / sum( truevals.10m ))  

homGP_PBias<- 100 * (sum( homGP_mean - truevals.10m ) / sum( truevals.10m ))  

homML_PBias

homGP_PBias 

#VERY SLIGHTLY smaller RMSE from HR than MR
#

################################################################################

#Calculate Fit and Correctness
res.e<- 10

mGridArea<- res.e^2
rGridArea<- res.e^2

Ar<- length(which(truevals.10m>0))*rGridArea


Am<- length(which(homML_mean>0))*mGridArea
Arm<- length(intersect(which(truevals.10m>0),which(homML_mean>0)))*rGridArea
Fvals_homML<- Arm/(Am + Ar - Arm)
Cvals_homML<- Arm/Ar



Am<- length(which(homGP_mean>0))*mGridArea
Arm<- length(intersect(which(truevals.10m>0),which(homGP_mean>0)))*rGridArea
Fvals_homGP<- Arm/(Am + Ar - Arm)
Cvals_homGP<- Arm/Ar


Fvals_homML
Fvals_homGP


Cvals_homML
Cvals_homGP

save(homML_RMSE, homGP_RMSE,Fvals_homML, Fvals_homGP, Cvals_homML, Cvals_homGP,
     file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/ED_F_Cvals.RData")

