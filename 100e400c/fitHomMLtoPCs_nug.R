#fit a Gaussian process emulator to each column of c(Y_R,YC_R) in parallel
rm(list=ls())

library(rstan)
library(boot)
library("R.matlab")
library(foreach)
library(doParallel)
source("/storage/work/svr5482/FloodingModelCalibrationProject/multires/code/myFunctions.R")
setwd("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main")
#source("GPfunctionsOptim.R")
#source("hetGPfunctions.R")
#rstan_options(auto_write = TRUE)
#homSML<- stan_model("hom-SML.stan",auto_write = TRUE)
#homML<- stan_model("hom-ML.stan")
#homMLExp<- stan_model("hom-MLExp.stan")
#homMLMat<- stan_model("hom-MLMatern.stan")

##################################with nugget###################################
#homMLExpnug<- stan_model("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main/hom-MLExp-nugget.stan")
homMLSqExpnug<- stan_model("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main/hom-ML-nugget.stan")
#homMLMat3_2nug<- stan_model("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main/hom-MLMat3_2-nugget.stan")
#homMLMat5_2nug<- stan_model("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main/hom-MLMat5_2-nugget.stan")
################################################################################

nRuns10m<- 100
nRuns50m<- 400
res.e<-10
res.c<-50
nLoc10m<- 126791
nLoc10min50m<- 125050
nLoc50m<- 5146
nPars=2
nPCs=18

#load parameters
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/100e400c/nChRWE_lhs_samples_allU.RData")

parVals10m<- data.frame("run"=1:nRuns10m, "n_ch"= samp.E[,1], "rwe"= samp.E[,2])
parVals50m<- data.frame("run"=1:nRuns50m, "n_ch"= samp.C[,1], "rwe"= samp.C[,2])

#parsTrue<- c(.0305,1)

#x.test<- data.frame("n_ch"= parsTrue[1],"rwe"= parsTrue[2])

x.c <- parVals50m[,-1]
x.e <- parVals10m[,-1]

################################################################################


#setup parallel backend to use many processors
cores=detectCores()
cl <- parallel::makeCluster(cores[1]-1) # -1 not to overload system
registerDoParallel(cl)


foreach(k = 1:nPCs)%dopar%{
  
  #load PCs from both 10m and 50m obs
  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))
  
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
  #x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale")) 
  
  
  ## need to get the data into shape!
  n.ml <- length(y.e)
  
  y.c2 <- rev(y.c)
  y.e2 <- rev(y.e)
  n.c <- length(y.c2)
  n.e <- length(y.e2)
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
  
  ml.data <- list(
    
    ## data ##
    
    m_p = (nPars+1),
    N = n.e + n.c, 
    K = nPars,
    n_c = n.c, 
    n_e = n.e,
    x_e = x.e.std.rev, 
    x_c = x.c.std.rev,
    m_H = m.h,
    y_c = y.c2, 
    y_e = y.e2,
    
    ## priors ##
    
    m_beta_m = rep(0, (nPars+1)), 
    m_beta_s = var.beta[1:(nPars+1),1:(nPars+1)],
    
    ##########################
    #for sqexp
    m_a_theta = rep(2, nPars), 
    m_b_theta = rep(2, nPars), #changed to gamma(2,2) from gammma (2,1 to reflect wider range of correlations)
    ##########################
    
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget_a = 2, m_nugget_b = 2, #uncomment with nugget
    
    c_beta_m = rep(0, (nPars+1)), 
    c_beta_s = var.beta[1:(nPars+1),1:(nPars+1)],
    
    ##########################
    #for sqexp
    c_a_theta = rep(2, nPars), 
    c_b_theta = rep(2, nPars), #changed to gamma(2,2) from gammma (2,1 to reflect wider range of correlations)
    ##########################
    
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 2, #uncomment with nugget
    
    m_rho = 1, 
    s_rho = 1/3
    
  )
  
  
  temp <- list()
  
  ## fit multilevel GP
  
  find.mode <- function(x){
    #squared exponential
    #rstan::optimizing(homML, data = ml.data, verbose = F, as_vector = F)
    #exponential
    #rstan::optimizing(homMLExp, data = ml.data, verbose = F, as_vector = F)
    #matern
    #rstan::optimizing(homMLMat3_2nug, data = ml.data, verbose = F, as_vector = F)
    #exponential with nugget
    #rstan::optimizing(homMLExpnug, data = ml.data, verbose = F, as_vector = F)
    #squared exponential with nugget
    rstan::optimizing(homMLSqExpnug, data = ml.data, verbose = F, as_vector = F)
  }
  st1 <- Sys.time()
  temp <- parallel::mclapply(1:3, find.mode, mc.cores = 3)
  en1 <- Sys.time()
  time= en1-st1
  #beepr::beep(4)
  
  best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))
  c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value)
  
  ml.fit <-  temp[[best.emulator]]
  pars <- ml.fit$par
  
  #squared exponential
  #save(pars,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homML.PC",k,".RData",sep=""))
  #save(time,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/time.homML.PC",k,".RData",sep=""))
  
  #exponential
  #save(pars,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLExp.PC",k,".RData",sep=""))
  #save(time,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/time.homMLExp.PC",k,".RData",sep=""))
  
  #matern
  #save(pars,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLMat3_2nug.PC",k,".RData",sep=""))
  #save(time,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/time.homMLMat3_2nug.PC",k,".RData",sep=""))
  
  #exponential with nugget
  #save(pars,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLExpnug.PC",k,".RData",sep=""))
  #save(time,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/time.homMLExpnug.PC",k,".RData",sep=""))
  
  #squared exponential with nugget
  save(pars,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLSqExpnug.PC",k,".RData",sep=""))
  save(time,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/time.homMLSqExpnug.PC",k,".RData",sep=""))
  
}

stopCluster(cl)

#nPCs=30

#times_MR= rep(NA,nPCs)
#for(k in 1:nPCs){
#  print(k)
#  load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLExp.PC",k,".RData",sep=""))
#  #load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/downscale/",nPCs,"PCs/time.homMLExp.PC",k,".RData",sep=""))
#  print(pars)
#  #times_MR[k]<- time
#}

#mean(times_MR); sd(times_MR)