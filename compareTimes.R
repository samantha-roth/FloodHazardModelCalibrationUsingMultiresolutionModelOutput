#add times

#MR calibration 50e200c
nPCs=15
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/50e200c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_time_rw_MH_homMR_sqexp_nug.RData")
calTime_MR50e200c<- (ptFinal)/(60^2)

#HR calibration 50e200c
nPCs=7
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/50e200c/just10m/cal/",nPCs,"PCs/.25/3cm",sep=""))  
load("vattime_rw_MH_homGP10_sqexp_nug.RData")
calTime_HR50e200c<- (ptFinal*.75)/(60^2)


#MR calibration 100e400c

nPCs=18
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_time_rw_MH_homMR_sqexp_nug.RData")
ptFinal1<- ptFinal
load("vat_time_rw_MH_homMR_sqexp_nug2.RData")
ptFinal2<- ptFinal
ptFinal<- ptFinal1+ptFinal2
calTime_MR100e400c<- (ptFinal)/(60^2)

#HR calibration 100e400c
nPCs=8
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/cal/",nPCs,"PCs/.25/3cm",sep=""))  
load("vattime_rw_MH_homGP10_sqexp_nug.RData")
calTime_HR100e400c<- (ptFinal*.75)/(60^2)


#########################################LOHI###################################

nPCs=13
#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/downscale/cal/.25/3cm")
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_time_rw_MH_homMR_sqexp_nug.RData")
calTime_MR50e200clohi<- (ptFinal)/(60^2)

nPCs=5
#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/just10m/cal/.25/3cm")
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/just10m/cal/",nPCs,"PCs/.25/3cm",sep=""))  
load("vattime_rw_MH_homGP10_sqexp_nug.RData")
calTime_HR50e200clohi<- (ptFinal*.75)/(60^2)

nPCs=16
#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/100e400c/downscale/cal/.25/3cm")
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_time_rw_MH_homMR_sqexp_nug.RData")
calTime_MR100e400clohi<- (ptFinal)/(60^2)

nPCs=7
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/100e400c/just10m/cal/",nPCs,"PCs/.25/3cm",sep=""))  
load("vattime_rw_MH_homGP10_sqexp_nug.RData")
calTime_HR100e400clohi<- (ptFinal*.75)/(60^2)

calTime_MR50e200c
calTime_HR50e200c

calTime_MR100e400c
calTime_HR100e400c

calTime_MR50e200clohi
calTime_HR50e200clohi

calTime_MR100e400clohi
calTime_HR100e400clohi


################################################################################
#time for downscaling

#100e400c
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/time_downscale.RData")
time

#50e200c
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/50e200c/downscale/time_downscale.RData")
time

################################################################################
#time for emulation

#100e400cMR regular
nPCs=18
emulateTimesMR100e400c<- rep(NA,nPCs)
for(k in 1:nPCs){
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/time.homMLSqExpnug.PC",k,".RData",sep=""))
  print(time)
  emulateTimesMR100e400c[k]<-as.numeric(time)
}

totalTimesMR100e400c<- sum(emulateTimesMR100e400c)

#100e400cHR regular
nPCs=8
emulateTimesHR100e400c<- rep(NA,nPCs)
for(k in 1:nPCs){
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/",nPCs,"PCs/time.homGPSqExp10nug.PC",k,".RData",sep=""))
  print(time)
  emulateTimesHR100e400c[k]<-as.numeric(time)
}

totalTimesHR100e400c<- sum(emulateTimesHR100e400c)


#100e400cMR lohi
nPCs=16
emulateTimesMR100e400clohi<- rep(NA,nPCs)
for(k in 1:nPCs){
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/100e400c/downscale/",nPCs,"PCs/time.homMLSqExpnug.PC",k,".RData",sep=""))
  print(time)
  emulateTimesMR100e400clohi[k]<-as.numeric(time)
}

totalTimesMR100e400clohi<- sum(emulateTimesMR100e400clohi)

#100e400cHR lohi
nPCs=7
emulateTimesHR100e400clohi<- rep(NA,nPCs)
for(k in 1:nPCs){
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/100e400c/just10m/",nPCs,"PCs/time.homGPSqExp10nug.PC",k,".RData",sep=""))
  print(time)
  emulateTimesHR100e400clohi[k]<-as.numeric(time)
}

totalTimesHR100e400clohi<- sum(emulateTimesHR100e400clohi)

################################################################################

#50e200cMR regular
nPCs=15
emulateTimesMR50e200c<- rep(NA,nPCs)
for(k in 1:nPCs){
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/50e200c/downscale/",nPCs,"PCs/time.homMLSqExpnug.PC",k,".RData",sep=""))
  print(time)
  emulateTimesMR50e200c[k]<-as.numeric(time)
}

totalTimesMR50e200c<- sum(emulateTimesMR50e200c)


#50e200cHR regular
nPCs=7
emulateTimesHR50e200c<- rep(NA,nPCs)
for(k in 1:nPCs){
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/50e200c/just10m/",nPCs,"PCs/time.homGPSqExp10nug.PC",k,".RData",sep=""))
  print(time)
  emulateTimesHR50e200c[k]<-as.numeric(time)
}

totalTimesHR50e200c<- sum(emulateTimesHR50e200c)


#50e200cMR lohi
nPCs=13
emulateTimesMR50e200clohi<- rep(NA,nPCs)
for(k in 1:nPCs){
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/downscale/",nPCs,"PCs/time.homMLSqExpnug.PC",k,".RData",sep=""))
  print(time)
  emulateTimesMR50e200clohi[k]<-as.numeric(time)
}

totalTimesMR50e200clohi<- sum(emulateTimesMR50e200clohi)

#50e200cHR lohi
nPCs=5
emulateTimesHR50e200clohi<- rep(NA,nPCs)
for(k in 1:nPCs){
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/just10m/",nPCs,"PCs/time.homGPSqExp10nug.PC",k,".RData",sep=""))
  print(time)
  emulateTimesHR50e200clohi[k]<-as.numeric(time)
}

totalTimesHR50e200clohi<- sum(emulateTimesHR50e200clohi)



totalTimesMR100e400c
totalTimesMR100e400clohi
totalTimesMR50e200c
totalTimesMR50e200clohi

totalTimesHR100e400c
totalTimesHR100e400clohi
totalTimesHR50e200c
totalTimesHR50e200clohi

totalTimesMR100e400c/18
totalTimesMR100e400clohi/16
totalTimesMR50e200c/15
totalTimesMR50e200clohi/13

totalTimesHR100e400c/8
totalTimesHR100e400clohi/7
totalTimesHR50e200c/7
totalTimesHR50e200clohi/5

load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/100e400c/downscale/time_getPCs.RData")