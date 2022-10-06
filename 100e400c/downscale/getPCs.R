#disaggregate then get PCs with all model output at 10m resolution
rm(list=ls())

#library(raster)
library(boot)
library(ggplot2)
library(viridis)
library(stats)

nRuns10m<- 100
nRuns50m<- 400
res.e<-10
res.c<-50
nLoc10m<- 126791
nLoc10min50m<- 125050
nLoc50m<- 5146

p = nRuns10m+ nRuns50m
n = nLoc10min50m

#RunTrue.10m= raster("/storage/work/svr5482/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/Extent/RunTrue_1.asc")
#coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
#rm(RunTrue.10m)

##load the coordinates where 10m runs are within the bounds of the 50m runs
#load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mCoordsin50m.RData")

#setwd("/storage/work/svr5482/FloodingModelCalibrationProject/04-Spatial_Stats_Samantha/Outputs10m/nCh_RWE/100e400c/Extent")

#all10min50mVals<- matrix(NA, nrow= nRuns10m, ncol= nLoc10min50m)
#for(i in 1:nRuns10m){
#  run<- raster(paste("Run_",i,".asc",sep=""))
#  vals<- extract(run,coords.10m)
#  all10min50mVals[i,]<- vals[coordsIwantInds]
#}

#save(all10min50mVals,file= "/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/all10min50mVals.RData")
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/all10min50mVals.RData")

st1 <- Sys.time()

#load downscaled 50m values
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/Runs50mDownscaled.RData")

##concatenate the matrix of 10m preds with the matrix of 50m preds downscaled via bilinear interpolation to 10m
preds10m50mat10m<- rbind(all10min50mVals,Runs50mDownscaled)

colSumPreds<- colSums(preds10m50mat10m)

whichAllZero<- which(colSumPreds==0)
save(whichAllZero,file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/whichAllZero.RData")

preds10m50mat10m.NZ<- preds10m50mat10m[,-whichAllZero]

preds10m50mat10m.NZ.C<- scale(preds10m50mat10m.NZ,scale=FALSE) #center the data

meanPreds10m50mat10m.NZ<- colMeans(preds10m50mat10m.NZ)

#SVD.preds10m50mat10m.NZ.C<- svd(preds10m50mat10m.NZ.C)

#X= UDV'= ULA'
#the columns of V are the right singular vectors aka the vectors of PC loadings

mat.Y1=t(preds10m50mat10m.NZ.C) #transpose the centered observations

svd.Y=svd(mat.Y1) #apply SVD 

#Determine how many eigenvectors to use by proportion of variation explained
explained.vars=cumsum(svd.Y$d)/sum(svd.Y$d) #the proportion of explained variance

plot(1:(length(explained.vars)),explained.vars)
plot(1:50,explained.vars[1:50])
plot(1:100,explained.vars[1:100])
#Calculate proportion of explained variance

prop.var.g90<- which(explained.vars>.90)
prop.var.g95<- which(explained.vars>.95)
prop.var.g99<- which(explained.vars>.99)
explained.vars[explained.vars[1]] 
explained.vars[50]
explained.vars[30]
explained.vars[20]
explained.vars[10]

################################################################################

#set up so that I can try with different numbers of PCs ranging from 26 to 50
#explain at least 95% of variance

#for(nPCs in 20:30){
nPCs=18
J_y= nPCs
  
  ##my way
  #K_y= matrix(NA, nrow= nrow(SVD.preds10m50mat10m.NZ.C$v), ncol= J_y)
  #compute K_y
  #for(i in 1:J_y){
  #  K_y[,i]<- SVD.preds10m50mat10m.NZ.C$d[i]*SVD.preds10m50mat10m.NZ.C$v[,i]
  #}
  #WC way- difference is that they are orthonormalized
  K_y=-svd.Y$u%*%diag(svd.Y$d)[,1:J_y]/sqrt(p) #orthonormalized Principal component
  
  #get the matrix Y_R (uncentered since I have a mean function)
  
  Y_R<- t(solve(t(K_y)%*%K_y)%*%t(K_y)%*%t(preds10m50mat10m.NZ.C))

  
  if(dir.exists(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs",sep=""))==F)
  {dir.create(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs",sep=""))}
  
  
  save(Y_R,K_y,J_y,meanPreds10m50mat10m.NZ,
       file= paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))
  
#}

en1 <- Sys.time()
time= en1-st1
save(time,file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/time_getPCs.RData")

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))
