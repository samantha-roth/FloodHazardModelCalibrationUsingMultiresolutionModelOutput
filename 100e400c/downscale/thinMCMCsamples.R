rm(list=ls())

nPCs=18
##setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25/3cm")
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_rw_MH_homMR_sqexp_nug.RData")
load("vatprops_rw_MH_homMR_sqexp_nug.RData")
load("vat_time_rw_MH_homMR_sqexp_nug.RData")

res1<- res
ptFinal1<- ptFinal

setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_rw_MH_homMR_sqexp_nug2.RData")
load("vatprops_rw_MH_homMR_sqexp_nug2.RData")
load("vat_time_rw_MH_homMR_sqexp_nug2.RData")

res2<- res
ptFinal2<- ptFinal

res<- rbind(res1,res2); rm(res1); rm(res2)

nThin<- 100
nSamples<- 3e5

thinInds<- seq(from=nSamples/nThin, to= nSamples, by= nSamples/nThin)
resThin<- res[thinInds,]
plot(density(resThin[,1])); plot(density(resThin[,2]))


#set.seed(29)
#inds2<-sample(1:nrow(res),100)
#resThin<- res[inds2,]

plot(density(resThin[,1])); plot(density(resThin[,2]))

save(resThin,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/resThin.RData",sep=""))


nPCs=18
load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/resThin.RData",sep=""))
plot(density(resThin[,1])); plot(density(resThin[,2]))

