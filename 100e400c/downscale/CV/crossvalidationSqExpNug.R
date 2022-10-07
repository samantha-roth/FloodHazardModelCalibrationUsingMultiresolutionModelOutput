#disaggregate then get PCs with all model output at 10m resolution
rm(list=ls())

#library(raster)
library(boot)
library(ggplot2)
library(viridis)
library(stats)

library(rstan)
library(boot)
library("R.matlab")
library(foreach)
library(doParallel)
library(emulator)
library(mgcv)
source("/storage/work/svr5482/FloodingModelCalibrationProject/multires/code/myFunctions.R")

##################################with nugget###################################
#homMLExpnug<- stan_model("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main/hom-MLExp-nugget.stan")
homMLSqExpnug<- stan_model("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main/hom-ML-nugget.stan")
#homMLMat3_2nug<- stan_model("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main/hom-MLMat3_2-nugget.stan")
#homMLMat5_2nug<- stan_model("/storage/work/svr5482/FloodingModelCalibrationProject/sml-athena-main/hom-MLMat5_2-nugget.stan")
################################################################################

nRuns10m<- 100
nRuns50m<- 400
nRuns= nRuns10m+nRuns50m
res.e<-10
res.c<-50
nLoc10m<- 126791
nLoc10min50m<- 125050
nLoc50m<- 5146
nPars=2
nCV= nExtra=10

J_y= nPCs=18
load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))

#load parameters
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/100e400c/nChRWE_lhs_samples_allU.RData")

parVals10m<- data.frame("run"=1:nRuns10m, "n_ch"= samp.E[,1], "rwe"= samp.E[,2])
parVals50m<- data.frame("run"=1:nRuns50m, "n_ch"= samp.C[,1], "rwe"= samp.C[,2])

#load predictions
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/all10min50mVals.RData")

#load downscaled 50m values
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/Runs50mDownscaled.RData")

##concatenate the matrix of 10m preds with the matrix of 50m preds downscaled via bilinear interpolation to 10m
preds10m50mat10m<- rbind(all10min50mVals,Runs50mDownscaled)

#load locations of zero observations
load("/storage/work/svr5482/FloodingModelCalibrationProject/multires/outputData/nCh_RWE/100e400c/whichAllZero.RData")

preds10m50mat10m.NZ<- preds10m50mat10m[,-whichAllZero]

all10mVals.NZ<- all10min50mVals[,-whichAllZero]


#divide into sets 
set.seed(79)
randomInds<- sample(1:nrow(all10mVals.NZ), nrow(all10mVals.NZ), replace = FALSE)
CVsize=10
nExtra= CVsize
allInds<- 1:nRuns10m

RMSE_MR_CV<- rep(NA, nRuns10m)

#setup parallel backend to use many processors
#cores=detectCores()
#cl <- parallel::makeCluster(cores[1]-1) # -1 not to overload system
#registerDoParallel(cl)

v=CVsize-(nPars+1)

USE_MR_CV<- matrix(NA, nrow= nRuns10m, ncol= nPCs)
eta.store= matrix(NA, nrow= nRuns10m, ncol= nPCs)

for(i in 1:nCV){
  #foreach(i = 1:nCV)%dopar%{
  testInds<- randomInds[((CVsize*(i-1))+1):(i*CVsize)]
  trainInds<- randomInds[-testInds]
  
  Y_Rtrain= rbind(Y_R[-c(testInds,nRuns10m+testInds),],Y_R[nRuns10m+testInds,])
  holdOutPCs= Y_R[testInds,]
  
  x.c <- rbind(parVals50m[-testInds,-1], parVals50m[testInds,-1])
  x.e <- parVals10m[-testInds,-1]
  x.test<- parVals10m[testInds,-1]
  
  ##all10mVals.NZ.train<- all10mVals.NZ[trainInds,]
  #preds10m50mat10m.NZ.train<-preds10m50mat10m.NZ[-c(testInds,nRuns10m+testInds),]
  all10mVals.NZ.test<- all10mVals.NZ[testInds,]
  
  
  #preds10m50mat10m.NZ.C<- scale(preds10m50mat10m.NZ.train,scale=FALSE) #center the data
  
  ##X= UDV'= ULA'
  ##the columns of V are the right singular vectors aka the vectors of PC loadings
  
  #mat.Y1=t(preds10m50mat10m.NZ.C) #transpose the centered observations
  
  #svd.Y=svd(mat.Y1) #apply SVD 
  
  ##Determine how many eigenvectors to use by proportion of variation explained
  #explained.vars=cumsum(svd.Y$d)/sum(svd.Y$d) #the proportion of explained variance
  
  ##Calculate proportion of explained variance
  
  #prop.var.g95<- which(explained.vars>.95)
  #nPCs=prop.var.g95[1]
  #J_y= nPCs
  
  #prop.var.g50<- which(explained.vars>.50)
  #prop.var.g63<- which(explained.vars>.63)
  #prop.var.g75<- which(explained.vars>.75)
  #prop.var.g90<- which(explained.vars>.90)
  
  #atleast50<- prop.var.g50[1]; atleast63<- prop.var.g63[1]
  #atleast75<- prop.var.g75[1]; atleast90<- prop.var.g90[1]
  #save(atleast50,atleast63,atleast75,atleast90,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/emulatorPCsExplainHomMLSqExpnugCV",i,".RData",sep=""))
  
  ################################################################################
  
  
  #K_y=-svd.Y$u%*%diag(svd.Y$d)[,1:J_y]/sqrt(p) #orthonormalized Principal component
  
  ##get the matrix Y_R (uncentered since I have a mean function)
  
  #Y_R<- t(solve(t(K_y)%*%K_y)%*%t(K_y)%*%t(preds10m50mat10m.NZ.train))
  
  #holdOutPCs<- t(solve(t(K_y)%*%K_y)%*%t(K_y)%*%t(all10mVals.NZ.test))
  
  allPCpreds<- matrix(NA, nrow=nExtra,ncol=nPCs)
  
  #setup parallel backend to use many processors
  #cores=detectCores()
  #cl <- parallel::makeCluster(cores[1]-1) # -1 not to overload system
  #registerDoParallel(cl)
  
  for(k in 1:nPCs){
    print(k)
    library(emulator)
    library(stats)
    
    Y_R1E<- Y_Rtrain[1:(nRuns10m-CVsize),k]
    Y_R1C<- Y_Rtrain[(nRuns10m-CVsize+1):nrow(Y_Rtrain),k] #Y_R1C<- Y_Rtrain[(nRuns10m-CVsize+1):(nRuns-2*CVsize),k]
    
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
      m_a_theta = rep(2, nPars), 
      m_b_theta = rep(2, nPars),
      m_a_sigma = 2, m_b_sigma = 2,
      m_nugget_a = 2, m_nugget_b = 2, #uncomment with nugget
      
      c_beta_m = rep(0, (nPars+1)), 
      c_beta_s = var.beta[1:(nPars+1),1:(nPars+1)],
      c_a_theta = rep(2, nPars), 
      c_b_theta = rep(2, nPars),
      c_a_sigma = 2, c_b_sigma = 2,
      c_nugget_a = 2, c_nugget_b = 2, #uncomment with nugget
      
      m_rho = 1, 
      s_rho = 1/3
      
    )
    
    temp <- list()
    
    ## fit multilevel GP
    
    find.mode <- function(x){
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
    
    designmat <- matrix(0, ncol=ncol(m.h), nrow=length(c(y.c, y.e)))
    designmat[1:length(y.c),1:(nPars+1)] <- cbind(1, x.c.rev)
    designmat[-(1:length(y.c)), (1:(nPars+1))] <- pars$rho*cbind(1, x.e.rev)
    designmat[-(1:length(y.c)), -(1:(nPars+1))] <- cbind(1, x.e.rev)
    
    #exponential
    dat.covmat <- dataCovmatHomNug(x.c.std.rev, x.e.std.rev, pars, var.beta, designmat, covfunc="sqexp")
    
    precmat.ml <- chol2inv(chol(dat.covmat))
    rm(dat.covmat)
    
    H.c <- cbind(1, x.c.rev)
    H.e <- cbind(1, as.matrix(x.c[n.e:1,]))
    H <- cbind(1, as.matrix(x.test))
    
    var.bc <- var.beta[1:(nPars+1),1:(nPars+1)]
    var.be <- var.beta[(nPars+2):(ncol(m.h)), (nPars+2):(ncol(m.h))]
    
    #exponential
    crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be, covfunc="sqexp")
    
    #assume a priori E(y|x) = 0 for any x
    mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))
    
    #exponential
    var.ml1 <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,covfunc="sqexp") + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,covfunc="sqexp")
    
    var.ml <- var.ml1 + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
    
    #no nugget
    #var.ml.full <- var.ml
    
    #nugget
    var.ml.full <- var.ml+ diag(pars$m_nugget,nrow(var.ml))
    
    allPCpreds[,k]<- mean.ml
    
    if(dir.exists(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs",sep=""))==F){
      dir.create(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs",sep=""))}
    
    save(mean.ml,var.ml.full,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/emulatorPredsHomMLSqExpnug.PC",k,"CV",i,".RData",sep=""))
    #load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/emulatorPredsHomMLSqExpnug.PC",k,"CV",i,".RData",sep=""))
    
    #USE_MR_CV[testInds,k]<- (mean.ml- holdOutPCs[,k])/sqrt(diag(var.ml.full))
    
    
    if(k<5){
      #compute uncorrelated standard errors using pivoted cholesky decomposition
      
      #For PC k: 
      #Sigma11.2 = covariance matrix for prediction at held out parameter settings
      #eta = emulator prediction at held out parameter settings
      #mat.Y.red[leave.out,k]= observation at held out parameter
      #q number held out
      
      ######################################################
      # Compute standardized uncorrelated residuals        #
      ######################################################
      #leave.out=c(25,75,125,175,225) #...for these parameter settings
      
      #cov.red=zeta.mat[i]*diag(1,p)+kappa.mat[i]*exp(-(C1/phi.mat[i,1])^2-(C2/phi.mat[i,2])^2-(C3/phi.mat[i,3])^2)
      #Sigma11=cov.red[leave.out,leave.out]
      #Sigma22=cov.red[-leave.out,-leave.out] #Sigma_22 matrix
      #Sigma12=matrix(cov.red[leave.out,-leave.out],q,250-q) #Sigma matrix
      #eta=Sigma12%*%solve(Sigma22,mat.Y.red[-leave.out,i])
      #eta.store=c(eta.store,eta)
      #Sigma11.2=Sigma11-Sigma12%*%solve(Sigma22,t(Sigma12))
      
      var.ml.full[lower.tri(var.ml.full)]=t(var.ml.full)[lower.tri(t(var.ml.full))]
      #what does this mean?
      USE_MR_CV[testInds,k]=mgcv::mroot(solve(var.ml.full),method="svd",rank=CVsize)%*%(holdOutPCs[,k]-mean.ml)
      
      ######################################################
      # Creat diagnostics plots                            #
      ######################################################
      #fitted value vs residuals
      jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/plotMeanVSUSE_HomMLSqExpnug.PC",k,"CV",i,".jpeg",sep=""),
           width = 800, height = 600)
      #par(mar=c(5,5,1,1))
      plot(mean.ml,USE_MR_CV[testInds,k],ylab="Standardized Residual",xlab="Fitted Value",cex.lab=1.6,cex.axis=1.6,cex=3, main= paste("CV ",i,"PC ",k,sep=""))
      dev.off()
      
      #q-q plots for df=v
      library(car)
      #png("diag2.png")
      jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/plotqqtUSE_HomMLSqExpnug.PC",k,"CV",i,".jpeg",sep=""),
           width = 800, height = 600)
      #par(mar=c(5,5,1,1))
      #qqt(resid,df=5,cex.lab=1.6,cex.axis=1.6,main="",cex=3)
      qqPlot(USE_MR_CV[testInds,k], dist="t", df=v, ylab= "Standardized Residual", main= paste("CV ",i,"PC ",k,sep=""))
      #qqplot(resid,distribution="t",df=q)
      #qqline(resid)
      dev.off()
      
      #Channel roughness vs residuals
      jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/nChUSE_HomMLSqExpnug.PC",k,"CV",i,".jpeg",sep=""),
           width = 800, height = 600)
      #par(mar=c(5,5,1,1))
      plot(x.test[,1],USE_MR_CV[testInds,k],ylab="Standardized Residual",xlab="Channel Roughness",cex.lab=1.6,cex.axis=1.6,cex=3, main= paste("CV ",i,"PC ",k,sep=""))
      #par(mar=c(5,5,1,1))
      dev.off()
      
      #RWE vs residuals
      jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/RWEUSE_HomMLSqExpnug.PC",k,"CV",i,".jpeg",sep=""),
           width = 800, height = 600)
      #par(mar=c(5,5,1,1))
      plot(x.test[,2],USE_MR_CV[testInds,k],ylab="Standardized Residual",xlab="River width error",cex.lab=1.6,cex.axis=1.6,cex=3, main= paste("CV ",i,"PC ",k,sep=""))
      #par(mar=c(5,5,1,1))
      dev.off()
    }
    
  }
  
  #stopCluster(cl)
  
  centeredPredsAtScale<- allPCpreds%*%t(K_y)  #nPars by nLoc = nPars by nPCs %*% nPCs by nLoc
  
  #predsAtScale<- apply(centeredPredsAtScale,1,function(x) x+ mean.Y )

  for(j in 1:nExtra){
    RMSE_MR_CV[testInds[j]]<- sqrt(mean((centeredPredsAtScale[j,] + mean.Y - all10mVals.NZ.test[j,])^2))
  }
  
  if(i==10){
    MRpredsAtScale<- matrix(0,nrow=nrow(centeredPredsAtScale),ncol=ncol(centeredPredsAtScale))
    for(i1 in 1:nExtra){
      MRpredsAtScale[i1,]<- centeredPredsAtScale[i1,] + mean.Y
    }
    
    save(MRpredsAtScale,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/predsAtScaleCV10SqExpnugCV.RData",sep=""))
  }
  
}



#exponential
#save(RMSE_downscale,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/RMSE_ExpExtra.RData",sep=""))

#exponential with nugget
save(RMSE_MR_CV,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/RMSE_SqExpnugCV.RData",sep=""))

save(USE_MR_CV,file=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/USE_SqExpnugCV.RData",sep=""))

#stopCluster(cl)

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/RMSE_SqExpnugCV.RData",sep=""))

#load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/USE_SqExpnugCV.RData",sep=""))


################################################################################

rm(list=ls())

library(car)

nRuns10m<- 100
nRuns50m<- 400
nRuns= nRuns10m+nRuns50m
res.e<-10
res.c<-50
nLoc10m<- 126791
nLoc10min50m<- 125050
nLoc50m<- 5146
nPars=2
nCV= nExtra=10

J_y= nPCs=18
load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/Y_R_",nPCs,"PCs.RData",sep=""))

#load parameters
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/100e400c/nChRWE_lhs_samples_allU.RData")

parVals10m<- data.frame("run"=1:nRuns10m, "n_ch"= samp.E[,1], "rwe"= samp.E[,2])
parVals50m<- data.frame("run"=1:nRuns50m, "n_ch"= samp.C[,1], "rwe"= samp.C[,2])



#divide into sets 
set.seed(79)
randomInds<- sample(1:nrow(parVals10m), nrow(parVals10m), replace = FALSE)
CVsize=10
nExtra= CVsize
allInds<- 1:nRuns10m


#v=CVsize-(nPars+1)
v= CVsize-1

allPCpreds<- matrix(NA,nrow=nrow(parVals10m),ncol=nPCs)
for(i in 1:nCV){
  for(k in 1:nPCs){
    load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/emulatorPredsHomMLSqExpnug.PC",k,"CV",i,".RData",sep=""))
    testInds<- randomInds[((CVsize*(i-1))+1):(i*CVsize)]
    allPCpreds[testInds,k]<- mean.ml
  }
}

load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/USE_SqExpnugCV.RData",sep=""))

allUSE_MR_CV<- as.vector(USE_MR_CV[,1:4])

plot(rep(parVals10m$n_ch,times=4),allUSE_MR_CV)

plot(rep(parVals10m$rwe,times=4),allUSE_MR_CV)

plot(as.vector(allPCpreds[,1:4]),allUSE_MR_CV)

qqPlot(allUSE_MR_CV, dist="t", df=v, ylab= "Standardized Residual", main= paste("PC ",k,sep=""))

length(which(abs(allUSE_MR_CV)>2.37))/length(abs(allUSE_MR_CV)) #.05 level two sided alpha for df=7

summary(allUSE_MR_CV)


allUSE_MR_CV<- as.vector(USE_MR_CV[,1:4])

#channel roughness
jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/nChUSE_HomMLSqExpnug.jpeg",sep=""),
     width = 800, height = 600)
plot(rep(parVals10m$n_ch,times=4),allUSE_MR_CV, xlab= "Channel roughness",ylab= "Standardized Residual", main= "MR 100e400c")
dev.off()

#RWE vs residuals
jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/RWEUSE_HomMLSqExpnug.jpeg",sep=""),
     width = 800, height = 600)
plot(rep(parVals10m$rwe,times=4),allUSE_MR_CV, xlab= "River width error",ylab= "Standardized Residual", main= "MR 100e400c")
dev.off()

jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/plotMeanVSUSE_HomMLSqExpnug.PC.jpeg",sep=""),
     width = 800, height = 600)
plot(as.vector(allPCpreds[,1:4]),allUSE_MR_CV,  xlab= "Predicted value", ylab= "Standardized Residual", main= "MR 100e400c")
dev.off()

jpeg(filename=paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/plotqqtUSE_HomMLSqExpnug.PC.jpeg",sep=""),
     width = 800, height = 600)
qqPlot(allUSE_MR_CV, dist="t", df=v, ylab= "Standardized Residual", main= "MR 100e400c")
dev.off()
