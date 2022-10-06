#Compare MCMC results for downscale
rm(list=ls())

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal")
load("aao_rw_MH_homMR.RData")
load("time_aao_rw_MH_homMR.RData")
res_aao_rw_MH<- res; rm(res)
time_aao_rw_MH<- ptFinal; rm(ptFinal)

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.5")
load("vat_rw_MH_homMR.RData")
load("vat_time_rw_MH_homMR.RData")
load("vatprops_rw_MH_homMR.RData")
res_vat_rw_MH.5<- res; rm(res)
time_vat_rw_MH.5<- ptFinal; rm(ptFinal)

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25")
load("vat_rw_MH_homMR.RData")
load("vat_time_rw_MH_homMR.RData")
load("vatprops_rw_MH_homMR.RData")
res_vat_rw_MH.25<- res; rm(res)
time_vat_rw_MH.25<- ptFinal; rm(ptFinal)

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.1")
load("vat_rw_MH_homMR.RData")
load("vat_time_rw_MH_homMR.RData")
load("vatprops_rw_MH_homMR.RData")
res_vat_rw_MH.1<- res; rm(res)
time_vat_rw_MH.1<- ptFinal; rm(ptFinal)

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.01")
load("vat_rw_MH_homMR.RData")
load("vat_time_rw_MH_homMR.RData")
load("vatprops_rw_MH_homMR.RData")
res_vat_rw_MH.01<- res; rm(res)
time_vat_rw_MH.01<- ptFinal; rm(ptFinal)

nPars=2

for(i in 1:(nPars+1)){
  print(paste0("AAO first 100k variable ",i, " number unique: " ,length(unique(res_aao_rw_MH[1:1e5,i]))))
  print(paste0("VAT.25 first 100k variable ",i, " number unique: " ,length(unique(res_vat_rw_MH.25[1:1e5,i]))))
  print(paste0("VAT.1 first 100k variable ",i, " number unique: " ,length(unique(res_vat_rw_MH.1[1:1e5,i]))))
  print(paste0("VAT.01 first 100k variable ",i, " number unique: " ,length(unique(res_vat_rw_MH.01[1:1e5,i]))))
}

library(batchmeans)

for(i in 1:(nPars+1)){
  print(i)
  print(paste("AAO first 100k ESS is ",ess(res_aao_rw_MH[1:1e5,i]), sep=""))
  print(paste("VAT.25 first 100k ESS is ",ess(res_vat_rw_MH.25[1:1e5,i]), sep=""))
  print(paste("VAT.1 first 100k ESS is ",ess(res_vat_rw_MH.1[1:1e5,i]), sep=""))
  print(paste("VAT.01 first 100k ESS is ",ess(res_vat_rw_MH.01[1:1e5,i]), sep=""))
}

for(i in 1:(nPars+1)){
  print(i)
  ess(res_aao_rw_MH[1:1e5,i],verbose=TRUE)
  ess(res_vat_rw_MH.25[1:1e5,i],verbose=TRUE)
  ess(res_vat_rw_MH.1[1:1e5,i],verbose=TRUE)
  ess(res_vat_rw_MH.01[1:1e5,i],verbose=TRUE)
}

for(i in 1:(nPars+1)){
  print(i)
  print(paste("AAO first 100k ESS is ",ess(res_aao_rw_MH[1:1e5,i],imse=FALSE), sep=""))
  print(paste("VAT.25 first 100k ESS is ",ess(res_vat_rw_MH.25[1:1e5,i],imse=FALSE), sep=""))
  print(paste("VAT.1 first 100k ESS is ",ess(res_vat_rw_MH.1[1:1e5,i],imse=FALSE), sep=""))
  print(paste("VAT.01 first 100k ESS is ",ess(res_vat_rw_MH.01[1:1e5,i],imse=FALSE), sep=""))
}


#effective sample sizes are slightly larger for the MR approach
#except for sigma2
#all at least 1000

for(i in 1:(nPars+1)){
  print(i)
  print(paste("AAO first 100k est is ",bm(res_aao_rw_MH[1:1e5,i])[1]," and the MCSE is ", bm(res_aao_rw_MH[1:1e5,i])[2], sep=""))
  print(paste("VAT.25 first 100k est is ",bm(res_vat_rw_MH.25[1:1e5,i])[1]," and the MCSE is ", bm(res_vat_rw_MH.25[1:1e5,i])[2], sep=""))
  print(paste("VAT.1 first 100k est is ",bm(res_vat_rw_MH.1[1:1e5,i])[1]," and the MCSE is ", bm(res_vat_rw_MH.1[1:1e5,i])[2], sep=""))
  print(paste("VAT.01 first 100k est is ",bm(res_vat_rw_MH.01[1:1e5,i])[1]," and the MCSE is ", bm(res_vat_rw_MH.01[1:1e5,i])[2], sep=""))
  }

for(i in 1:(nPars+1)){
  print(i)
  acf(res_MR[,i],main= paste("HomMR Variable ",i,sep=""))
  acf(res_MR10[,i],main= paste("HomMR10 Variable ",i,sep=""))
  acf(res_HR[,i],main= paste("HomGP Variable ",i,sep=""))
}


for(i in 1:(nPars+1)){
  print(i)
  hist(res_aao_rw_MH[1:1e5,i])
  hist(res_vat_rw_MH.25[1:1e5,i])
  hist(res_vat_rw_MH.1[1:1e5,i])
  hist(res_vat_rw_MH.01[1:1e5,i])
}