#get predictions for the untested model runs using HomMLExp emulator

rm(list=ls())
library(foreach)
library(doParallel)
library(emulator)

source("/storage/work/svr5482/FloodingModelCalibrationProject/multires/code/myFunctions.R")

nRuns10m<- 100
nRuns50m<- 400
res.e<-10
res.c<-50
nLoc10m<- 126791
nLoc10min50m<- 125050
nLoc50m<- 5146
nPars=2
nPCs=18

st <- Sys.time()

#load parameters
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/100e400c/nChRWE_lhs_samples_allU.RData")
parVals10m<- data.frame("run"=1:nRuns10m, "n_ch"= samp.E[,1], "rwe"= samp.E[,2])
parVals50m<- data.frame("run"=1:nRuns50m, "n_ch"= samp.C[,1], "rwe"= samp.C[,2])

#load test parameters
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/100e400c/nChRWE_lhs_Eextra.RData")
parVals10mExtra<- samp.Eextra

#divide parameter values into training and testing
x.test<- parVals10mExtra

x.c <- parVals50m[,-1]
x.e <- parVals10m[,-1]

#load PCs from both 10m and 50m obs
load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))

##setup parallel backend to use many processors
#cores=detectCores()
#cl <- makeCluster(cores[1]-1) # -1 not to overload system
#registerDoParallel(cl)


#foreach(k = 1:nPCs)%dopar%{
for(k in 1:nPCs){
  
  Y_R1E<- Y_R[1:nRuns10m,k]
  Y_R1C<- Y_R[(nRuns10m+1):(nRuns50m+nRuns10m),k]
  
  #divide expensive model output into training and testing
  #y.test <- Y_R1E[test.inds]
  y.c<- Y_R1C
  y.e<- Y_R1E
  
  ##training data
  x.e.std <- scale(x.e)
  x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
  
  ##scale validation data accordingly
  x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
  
  
  ## need to get the data into shape!
  n.ml <- length(y.e)
  
  y.c2 <- rev(y.c)
  y.e2 <- rev(y.e)
  n.c <- length(y.c2); n.e <- length(y.e2)
  x.c.std.rev <- x.c.std[n.c:1,]
  x.e.std.rev <- x.e.std[n.ml:1,]
  
  #pairs(cbind(y.c2, x.c.std))
  #pairs(cbind(y.e2, x.e.std))
  ## reverse everything to get the order correct ...
  n.c <- length(y.c2); n.e <- length(y.e2)
  
  ## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"
  
  x.c.rev<- as.matrix(x.c[n.c:1,])
  x.e.rev<- as.matrix(x.e[n.e:1,])
  
  m.h <- cbind(cbind(1, x.c.rev), rbind(matrix(0, ncol = (nPars+1), nrow = n.c - n.e), cbind(1, x.e.rev)))
  #tail(m.h)
  var.beta <- diag(1, ncol(m.h))
  
  #get emulator parameter estimates
  #load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLExp.PC",k,".RData",sep=""))
  
  #with nugget
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLExpnug.PC",k,".RData",sep=""))
  
  designmat <- matrix(0, ncol=ncol(m.h), nrow=length(c(y.c, y.e)))
  designmat[1:length(y.c),1:(nPars+1)] <- cbind(1, x.c.rev)
  designmat[-(1:length(y.c)), (1:(nPars+1))] <- pars$rho*cbind(1, x.e.rev)
  designmat[-(1:length(y.c)), -(1:(nPars+1))] <- cbind(1, x.e.rev)
  
  #squared exponential
  #dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta, designmat, power=2)
  #exponential
  #dat.covmat <- dataCovmatHom(x.c.std.rev, x.e.std.rev, pars, var.beta, designmat, covfunc="sqexp")
  #exponential with nugget
  dat.covmat <- dataCovmatHomNug(x.c.std.rev, x.e.std.rev, pars, var.beta, designmat, covfunc="sqexp")
  
  precmat.ml <- chol2inv(chol(dat.covmat))
  rm(dat.covmat)
  
  H.c <- cbind(1, x.c.rev)
  H.e <- cbind(1, as.matrix(x.c[n.e:1,]))
  H <- cbind(1, as.matrix(x.test))
  
  var.bc <- var.beta[1:(nPars+1),1:(nPars+1)]
  var.be <- var.beta[(nPars+2):(ncol(m.h)), (nPars+2):(ncol(m.h))]
  
  #squared exponential
  #crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be, power= 2)
  #exponential
  crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be, covfunc="sqexp")
  
  #assume a priori E(y|x) = 0 for any x
  mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))
  
  #squared exponential
  #var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,power=2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,power=2)
  #exponential
  var.ml1 <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,covfunc="sqexp") + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,covfunc="sqexp")
  
  var.ml <- var.ml1 + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
  
  #no nugget
  #var.ml.full <- var.ml
  
  #nugget
  var.ml.full <- var.ml+ diag(pars$m_nugget,nrow(var.ml))
  
  #MSE.homml <- mean((mean.ml-y.test)^2)
  
  #squared exponential
  #save(mean.ml,var.ml.full,MSE.homml,var.ml.full,
  #     file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/downscale/",nPCs,"PCs/emulatorPredsHomML.PC",k,"Extra600.RData",sep=""))
  
  #exponential
  #save(mean.ml,var.ml.full,
  #     file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/emulatorPredsHomMLExp.PC",k,"Extra.RData",sep=""))
  
  #exponential with nugget
  save(mean.ml,var.ml.full,
       file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/emulatorPredsHomMLExpnug.PC",k,"Extra.RData",sep=""))
  
}

en <- Sys.time() - st
time= en

#save(time,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/time_emulatorPredsHomMLExpExtra.RData",sep=""))

#nugget
save(time,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/time_emulatorPredsHomMLExpnugExtra.RData",sep=""))

#load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/emulatorPredsHomMLExp.PC",k,"Extra.RData",sep=""))


