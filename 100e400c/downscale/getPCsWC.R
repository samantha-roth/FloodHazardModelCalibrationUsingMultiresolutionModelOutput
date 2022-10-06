#get PCs using exact WC code with my data
#downscaled

rm(list=ls())

source("/storage/work/svr5482/FloodingModelCalibrationProject/WC2014_code/function/functions_reduced_d.r")

library(raster)
library(boot)
library(ggplot2)
library(viridis)
library(stats)
library(MASS)
library(Matrix)
library(fields)
library(mvtnorm)
library(spam)
library(mcmc)

nRuns10m<- 100
nRuns50m<- 400
res.e<-10
res.c<-50
nLoc10m<- 126791
nLoc10min50m<- 125050
nLoc50m<- 5146

st1 <- Sys.time()


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

#load downscaled 50m values
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/Runs50mDownscaled.RData")

##concatenate the matrix of 10m preds with the matrix of 50m preds downscaled via bilinear interpolation to 10m
preds10m50mat10m<- rbind(all10min50mVals,Runs50mDownscaled)

colSumPreds<- colSums(preds10m50mat10m)

whichAllZero<- which(colSumPreds==0)
save(whichAllZero,file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/whichAllZero.RData")

preds10m50mat10m.NZ<- preds10m50mat10m[,-whichAllZero]
mat.Y<- preds10m50mat10m.NZ
#########################################################################
# Defining the dimensionality of matrices                               #
# n : The number of spatial observations in UVic outputs                #
# p : The number of parameter settings in UVic outputs                  #
# J : The dimensionality of the reduced space                           #
#########################################################################
n=ncol(mat.Y)
p=nrow(mat.Y)
J_y= nPCs= J.eta=18 #explain 95% of variation

################################################
# Finding pricinpal components of mat.Y        #
# K_eta : basis matrix                         #
################################################

#mat.Y=matrix(Y,p,n) mat.Y is a matrix with computer model outputs. column: different spatial locations row: parameter settings
mat.Y1=t(scale(mat.Y,scale=FALSE)) #Center mat.Y
#svd.Y=svd(mat.Y1, LINPACK = FALSE) #apply SVD
svd.Y=svd(mat.Y1) #apply SVD #I removed LINPCACK argument
K_y=K.eta=-svd.Y$u%*%diag(svd.Y$d)[,1:J.eta]/sqrt(p) #orthonormalized Principal component
explained.var=sum(svd.Y$d[1:J.eta])/sum(svd.Y$d) #the proportion of explained variance
Y=Y.centered=as.vector(mat.Y1) #return the outcomes

explained.var

#change mat.Y as a centered one
mean.Y=apply(mat.Y,2,mean)
mat.Y=matrix(Y,p,n,byrow=T)

###################################################################
# Project Y matrix onto a reduced space                           #
###################################################################
Y_R= mat.Y.red=t(solve(t(K.eta)%*%K.eta)%*%t(K.eta)%*%t(mat.Y))


if(dir.exists(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs",sep=""))==F)
{dir.create(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs",sep=""))}

save(Y_R,K_y,J_y,mean.Y,
     file= paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))

en1 <- Sys.time()
time= en1-st1
save(time,file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/time_getPCs.RData")

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))

explained.vars=cumsum(svd.Y$d)/sum(svd.Y$d)

#look at amount of explained variaton for other PCs:
prop.var.g50<- which(explained.vars>.50)
prop.var.g63<- which(explained.vars>.63)
prop.var.g75<- which(explained.vars>.75)
prop.var.g90<- which(explained.vars>.90)

atleast50<- prop.var.g50[1]; atleast63<- prop.var.g63[1]
atleast75<- prop.var.g75[1]; atleast90<- prop.var.g90[1]
save(atleast50,atleast63,atleast75,atleast90,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/emulatorPCsExplainHomMLSqExpnug.RData",sep=""))
