
#compare by CVRMSE

rm(list=ls())

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/RMSE_SqExpnugCV.RData",sep=""))

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/RMSE_SqExpnugCV.RData",sep=""))


plot(density(RMSE_MR_CV-RMSE_10m_CV))

summary(RMSE_MR_CV-RMSE_10m_CV)
