#compare calibrate results for each parameter
#realistic setting
rm(list=ls())
library(ggplot2)

nPCs=8
#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/cal/.25/3cm")
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/just10m/cal/",nPCs,"PCs/.25/3cm",sep=""))  
load("vat_rw_MH_homGP10_sqexp_nug.RData")

res_HR<- res[1:3e5,]; rm(res)

nPCs=18
#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25/3cm")
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_rw_MH_homMR_sqexp_nug.RData")

res_MR<- res; rm(res)

nCh.df<- data.frame("value"= c(res_HR[,1],res_MR[,1]),
                    "source"= c(rep("HR",nrow(res_HR)),rep("MR",nrow(res_MR))))

RWE.df<- data.frame("value"= c(res_HR[,2],res_MR[,2]),
                    "source"= c(rep("HR",nrow(res_HR)),rep("MR",nrow(res_MR))))


jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/posterdensity_nCh_HRvsMR_100e400c.jpeg",sep=""),
     width = 800, height = 600)
ggplot(nCh.df, aes(x=value, fill=source)) +
  geom_density(alpha=.25)+
  geom_vline(xintercept=0.0305, linetype="dashed", color = "red") +
  ggtitle("Channel roughness- 100e400c")+
  #geom_text(x=.031, y=650, label=expression("n"["ch"]^"*1"*"=.0305"))+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()

jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/posterdensity_RWE_HRvsMR_100e400c.jpeg",sep=""),
     width = 800, height = 600)
ggplot(RWE.df, aes(x=value, fill=source)) +
  geom_density(alpha=.25)+
  geom_vline(xintercept=1, linetype="dashed", color = "red") +
  ggtitle("River width error- 100e400c")+
  #geom_text(x=.99, y=18, label=expression("RWE"^"*1"*"=1"))+
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()

mean(res_HR[,1]); mean(res_MR[,1])
mean(res_HR[,2]); mean(res_MR[,2])