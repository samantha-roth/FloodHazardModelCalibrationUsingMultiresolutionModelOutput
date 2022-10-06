load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/100e400c/just10m/RMSE_SqExpnugExtreme.RData",sep=""))

RMSE_HR_100e400clohi<- RMSE_10m_Extreme

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/100e400c/downscale/RMSE_SqExpnugExtreme.RData",sep=""))

RMSE_MR_100e400clohi<- RMSE_MR_Extreme


load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/just10m/RMSE_SqExpnugExtreme.RData",sep=""))

RMSE_HR_50e200clohi<- RMSE_10m_Extreme

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/downscale/RMSE_SqExpnugExtreme.RData",sep=""))

RMSE_MR_50e200clohi<- RMSE_MR_Extreme


##############################

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/RMSE_SqExpnugCV.RData",sep=""))

RMSE_HR_100e400c<- RMSE_10m_CV

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/RMSE_SqExpnugCV.RData",sep=""))

RMSE_MR_100e400c<- RMSE_MR_CV


load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/50e200c/just10m/RMSE_SqExpnugCV.RData",sep=""))

RMSE_HR_50e200c<- RMSE_10m_CV

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/50e200c/downscale/RMSE_SqExpnugCV.RData",sep=""))

RMSE_MR_50e200c<- RMSE_MR_CV


summary(RMSE_HR_100e400c); summary(RMSE_HR_100e400clohi)

summary(RMSE_MR_100e400c); summary(RMSE_MR_100e400clohi)

summary(RMSE_HR_50e200c); summary(RMSE_HR_50e200clohi)

summary(RMSE_MR_50e200c); summary(RMSE_MR_50e200clohi)