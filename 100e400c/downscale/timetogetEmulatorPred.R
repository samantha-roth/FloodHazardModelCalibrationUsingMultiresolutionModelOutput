#nCh and RWE only
#calibrate using HomML approach
#all at once using a more elegant metropolis-hastings algorithm
#proposal in X = logit((parameter-a)/(b-a)) space

rm(list=ls())

library(mvtnorm); library(invgamma); library(spam)
library(foreach); library(doParallel)
library(emulator); library(dplyr)

#functions I changed from functions in https://github.com/jcken95/sml-athena

setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires")
source("code/myFunctions.R")
source("spatial/changPC2014/10m50m/code/nCh_RWE/100e400c/downscale/calibrate/calFuncs.R")


nRuns10m<- 100
nRuns50m<- 400
nRuns<- nRuns10m+nRuns50m
res.e<-10
res.c<-50
nLoc10m<- 126791
nLoc10min50m<- 125050
nLoc50m<- 5146
nPars=2
#change
nPCs=18
nug=TRUE

#what covariance function is being used?
cf="sqexp"

#load parameters
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/100e400c/nChRWE_lhs_samples_allU.RData")

parVals10m<- data.frame("run"=1:nRuns10m, "n_ch"= samp.E[,1], "rwe"= samp.E[,2])
parVals50m<- data.frame("run"=1:nRuns50m, "n_ch"= samp.C[,1], "rwe"= samp.C[,2])

#no nugget
if(nug==FALSE) emulatorParsLoc<- paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLSqExp.PC",sep="")

#with nugget
if(nug==TRUE) emulatorParsLoc<- paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLSqExpnug.PC",sep="")
################################################################################
#SIMULATED OBSERVATION INFO
parsTrue<- c(.0305,1)

#load the simulated observation with errors added
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/noisytruth.RData")

#load the locations where the prediction is always zero
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/whichAllZero.RData")

noisytruth<- noisy2.3truevals.10m[-whichAllZero]

rm(noisy2truevals.10m); rm(noisy3truevals.10m); rm(noisy4truevals.10m)
################################################################################

#load predictions in PC space
load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))

invtKK<- chol2inv(chol(crossprod(K_y)))
#Z_R<- t(solve(t(K_y)%*%K_y)%*%t(K_y)%*%t(noisytruth-mean.Y))
Z_R<- as.numeric(solve(t(K_y)%*%K_y)%*%t(K_y)%*%(noisytruth-mean.Y))

#prepare data
x.true<- data.frame("n_ch"= parsTrue[1],"rwe"= parsTrue[2])

x.c <- parVals50m[,-1]
x.e <- parVals10m[,-1]

##training data
x.e.std <- scale(x.e)
x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))

n.ml<- n.e<- nRuns10m
n.c<- nRuns50m

x.c.rev<- as.matrix(x.c[n.c:1,])
x.e.rev<- as.matrix(x.e[n.e:1,])

x.c.std.rev <- x.c.std[n.c:1,]
x.e.std.rev <- x.e.std[n.ml:1,]

m.h <- cbind(cbind(1, x.c.rev), rbind(matrix(0, ncol = (ncol(x.e)+1), nrow = n.c - n.e), cbind(1, x.e.rev)))

var.beta <- diag(1, ncol(m.h))

var.bc <- var.beta[1:(ncol(x.e)+1),1:(ncol(x.e)+1)]
var.be <- var.beta[(ncol(x.e)+2):(ncol(m.h)), (ncol(x.e)+2):(ncol(m.h))]

H.c <- cbind(1, x.c.rev)
H.e <- cbind(1, as.matrix(x.c[n.e:1,]))

################################################################################
##load emulator parameters

#calculated design cov mat here to so we don't have to iterate over it
precmat<- array(NA,c(nRuns,nRuns,nPCs))
rho_V.bc_tH.c<- array(NA,c(nPars+1,nRuns50m,nPCs))
rho2_V.bcPlusV.be_tH.e<- array(NA,c(nPars+1,nRuns10m,nPCs))

for(k in 1:nPCs){
  #load emulator pars
  load(paste(emulatorParsLoc,k,".RData",sep=""))
  
  
  designmat <- matrix(0, ncol=ncol(m.h), nrow=nRuns)
  designmat[1:nRuns50m,1:(ncol(x.e)+1)] <- cbind(1, x.c.rev)
  designmat[-(1:nRuns50m), (1:(ncol(x.e)+1))] <- pars$rho*cbind(1, x.e.rev)
  designmat[-(1:nRuns50m), -(1:(ncol(x.e)+1))] <- cbind(1, x.e.rev)
  
  #squared exponential
  #dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta, designmat, power=2)
  
  #exponential
  if(nug==FALSE) dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta, designmat, covfunc=cf)
  
  #exponential with nugget
  if(nug==TRUE) dat.covmat <- dataCovmatHomNug(x.c.std.rev, x.e.std.rev, pars, var.beta, designmat, covfunc=cf)
  
  precmat.ml <- chol2inv(chol(dat.covmat))
  rm(dat.covmat)
  
  precmat[,,k]<- precmat.ml
  
  
  rho_V.bc_tH.c[,,k]<- pars$rho*tcrossprod(var.bc,H.c)
  
  rho2_V.bcPlusV.be_tH.e[,,k]<- tcrossprod(pars$rho^2*var.bc+var.be,H.e) #(pars$rho^2*var.bc + var.be) %*% t(H.e)
}

y.e2mat<- apply(Y_R[1:nRuns10m,],2,rev)
y.c2mat<- apply(Y_R[(nRuns10m+1):(nRuns10m+nRuns50m),],2,rev)


################################################################################
nPCs=18
setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_rw_MH_homMR_sqexp_nug.RData")
load("vatprops_rw_MH_homMR_sqexp_nug.RData")
load("vat_time_rw_MH_homMR_sqexp_nug.RData")

res1<- res
ptFinal1<- ptFinal

#setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))
load("vat_rw_MH_homMR_sqexp_nug2.RData")
load("vatprops_rw_MH_homMR_sqexp_nug2.RData")
load("vat_time_rw_MH_homMR_sqexp_nug2.RData")

res2<- res
ptFinal2<- ptFinal

res<- rbind(res1,res2)


resThin1000<- res[seq(300,3e5,by=300),]

pt<- proc.time()

for(i in 1:1000){
  x.test<- matrix(res[i,1:2], nrow=1, ncol=nPars)
  pred_mean_cov(x.test)
}

ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3] # End Time to be used in Effective Samples per Second Calculation
