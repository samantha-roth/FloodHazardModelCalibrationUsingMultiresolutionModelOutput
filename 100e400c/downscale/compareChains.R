#compare results from four different Markov chains

rm(list=ls())

nPCs=18

setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_rw_MH_homMR_sqexp_nug.RData")
load("vatprops_rw_MH_homMR_sqexp_nug.RData")
load("vat_time_rw_MH_homMR_sqexp_nug.RData")

resC1<- res[1:100000,]; rm(res)


setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain2",sep=""))
load("vat_rw_MH_homMR_sqexp_nug.RData")
load("vatprops_rw_MH_homMR_sqexp_nug.RData")
load("vat_time_rw_MH_homMR_sqexp_nug.RData")

resC2<- res; rm(res)

setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain3",sep=""))
load("vat_rw_MH_homMR_sqexp_nug.RData")
load("vatprops_rw_MH_homMR_sqexp_nug.RData")
load("vat_time_rw_MH_homMR_sqexp_nug.RData")

resC3<- res; rm(res)

setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))
load("vat_rw_MH_homMR_sqexp_nug.RData")
load("vatprops_rw_MH_homMR_sqexp_nug.RData")
load("vat_time_rw_MH_homMR_sqexp_nug.RData")

resC4<- res; rm(res)

plot(1:nrow(resC1),resC1[,1],type="l")
lines(1:nrow(resC2),resC2[,1],col="red")
lines(1:nrow(resC3),resC3[,1],col="blue")
lines(1:nrow(resC4),resC4[,1],col="green")

plot(1:nrow(resC1),resC1[,2],type="l")
lines(1:nrow(resC2),resC2[,2],col="red")
lines(1:nrow(resC3),resC3[,2],col="blue")
lines(1:nrow(resC4),resC4[,2],col="green")

plot(1:nrow(resC1),resC1[,3],type="l")
lines(1:nrow(resC2),resC2[,3],col="red")
lines(1:nrow(resC3),resC3[,3],col="blue")
lines(1:nrow(resC4),resC4[,3],col="green")

#all look the same

apply(resC1, 2, mean); apply(resC2, 2, mean); apply(resC3, 2, mean); apply(resC4, 2, mean)

apply(resC1, 2, sd); apply(resC2, 2, sd); apply(resC3, 2,sd); apply(resC4, 2, sd)

apply(resC1, 2, summary); apply(resC2, 2, summary); apply(resC3, 2,summary); apply(resC4, 2, summary)

#plot(density(resC1[,1])); lines(density(resC2[,1]),col="red"); lines(1:nrow(resC3),resC3[,1],col="blue")


nch_post<- data.frame(values = c(resC1[,1],
                                 resC2[,1],
                                 resC3[,1],
                                 resC4[,1]),
                      group = c(rep("C1", nrow(resC1)),
                                rep("C2", nrow(resC2)),
                                rep("C3", nrow(resC3)),
                                rep("C4", nrow(resC4))))

rwe_post<- data.frame(values = c(resC1[,2],
                                 resC2[,2],
                                 resC3[,2],
                                 resC4[,2]),
                      group = c(rep("C1", nrow(resC1)),
                                rep("C2", nrow(resC2)),
                                rep("C3", nrow(resC3)),
                                rep("C4", nrow(resC4))))


library(ggplot2)
ggplot(nch_post, aes(x = values, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.25, bins = 50)

ggplot(rwe_post, aes(x = values, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.25, bins = 50)