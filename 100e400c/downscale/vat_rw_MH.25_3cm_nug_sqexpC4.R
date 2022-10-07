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

#setup parallel backend to use many processors
#cores=detectCores()
#cl <- makeCluster(cores[1]-1) # -1 not to overload system
#registerDoParallel(cl)


#Code the Gibbs sampler.
mh.alg<- function(init, n.sample) {
  x.t <- init[1:nPars] #ch,rwe
  s2.t <- init[length(init)] #sigma^2
  
  x.out <- matrix(NA,nrow= n.sample+1,ncol=nPars)
  s2.out <- rep(NA,n.sample+1)
  
  x.out[1,]= x.t
  s2.out[1]= s2.t
  
  x.test<- matrix(x.out[1,], nrow=1, ncol=nPars)
  
  predMeanCov<- pred_mean_cov(x.test)
  
  predMean<- predMeanCov[[1]]
  predCov<- predMeanCov[[2]]
  
  proposals<- matrix(NA, nrow= n.sample, ncol= nPars+4)
  
  for (i in 1 : n.sample) {
    #if(i%%1000==0){print(paste("Iteration",i))}
    
    ############################################################################
    
    #n_ch
    u<- runif(1)
    
    propDistMean<- UtoR(as.numeric(x.out[i,1]),a=.02,b=.1)
    
    propX<- rnorm(1, mean = propDistMean, sd = .25) #was .5
    
    prop<- RtoU(propX,a=.02,b=.1)
    
    #record proposals
    proposals[i,1]<- prop
    
    prop.x.test<- matrix(c(prop,x.out[i,2]),nrow=1,ncol= nPars)
    
    propPredMeanCov<- pred_mean_cov(prop.x.test)
    
    propPredMean<- propPredMeanCov[[1]]
    propPredCov<- propPredMeanCov[[2]]
    
    logMetRatio1<- ll_func(mean=propPredMean, cov=propPredCov, sigma2= s2.out[i])
    
    logMetRatio2<- ll_func(mean=predMean, cov=predCov, sigma2= s2.out[i])
    
    
    logMetRatio<- logMetRatio1-logMetRatio2+
      logDetJ(as.numeric(x.out[i,1]),a= .02,b=.1)-
      logDetJ(prop,a=.02,b=.1)
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio){
      x.out[i+1,1]<- prop
      
      x.test<- matrix(c(x.out[i+1,1],x.out[i,2]),nrow=1,ncol=nPars)
      
      #pt<-proc.time() # Start Time
      predMeanCov<- pred_mean_cov(x.test)
      #ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3]
      
      predMean<- predMeanCov$pred_mean
      predCov<- predMeanCov$pred_cov
      
      #record acceptance
      proposals[i,nPars+2]<- 1
    }
    if(log(u)>logMetRatio){
      x.out[i+1,1]<- x.out[i,1]
      #record rejection
      proposals[i,nPars+2]<- 0
    } 
    
    ############################################################################
    
    #rwe
    u<- runif(1)
    
    propDistMean<- UtoR(as.numeric(x.out[i,2]),a=.95,b=1.05)
    
    propX<- rnorm(1, mean = propDistMean, sd = .25) #was .5
    
    prop<- RtoU(propX,a=.95,b=1.05)
    
    #record proposals
    proposals[i,2]<- prop
    
    prop.x.test<- matrix(c(x.out[i+1,1],prop),nrow=1,ncol= nPars)
    
    propPredMeanCov<- pred_mean_cov(prop.x.test)
    
    propPredMean<- propPredMeanCov[[1]]
    propPredCov<- propPredMeanCov[[2]]
    
    logMetRatio1<- ll_func(mean=propPredMean, cov=propPredCov, sigma2= s2.out[i])
    
    logMetRatio2<- ll_func(mean=predMean, cov=predCov, sigma2= s2.out[i])
    
    
    logMetRatio<- logMetRatio1-logMetRatio2+
      logDetJ(as.numeric(x.out[i,2]),a=.95, b=1.05)-
      logDetJ(prop,a=.95, b=1.05)
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio){
      x.out[i+1,2]<- prop
      
      x.test<- matrix(x.out[i+1,],nrow=1,ncol=nPars)
      
      #pt<-proc.time() # Start Time
      predMeanCov<- pred_mean_cov(x.test)
      #ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3]
      
      predMean<- predMeanCov$pred_mean
      predCov<- predMeanCov$pred_cov
      
      #record acceptance
      proposals[i,nPars+3]<- 1
    }
    if(log(u)>logMetRatio){
      x.out[i+1,2]<- x.out[i,2]
      #record rejection
      proposals[i,nPars+3]<- 0
    } 
    
    
    ############################################################################
    #MH algorithm for sigma^2
    u<- runif(1)
    propX<- rnorm(1,log(s2.out[i]),sd= 1) #was 3
    prop<- exp(propX); 
    
    #record proposal
    proposals[i,nPars+1]<- prop
    
    #prop<- rnorm(1,s2.out[i],sd= 5) #symmetric proposal
    
    logMetRatio1<- ll_func(mean=predMean, cov=predCov, sigma2= prop) - 
      (.2+1)*log(prop) - (.2/prop)
    logMetRatio2<- ll_func(mean=predMean, cov=predCov, sigma2= s2.out[i]) - 
      (.2+1)*log(s2.out[i]) - (.2/s2.out[i])
    
    logMetRatio<- logMetRatio1-logMetRatio2+propX-log(s2.out[i])
    
    #metropolis ratio in the metropolis hastings algorithm 
    if(log(u)<=logMetRatio){
      s2.out[i+1]<- prop
      #record acceptance
      proposals[i,nPars+4]=1
    } 
    
    if(log(u)>logMetRatio){
      s2.out[i+1]<- s2.out[i]
      
      #record rejection
      proposals[i,nPars+4]=0
    } 
    
  }
  out <- cbind(x.out, s2.out)
  list(out=out,proposals=proposals)
}


if(dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale")==F)
{dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale")}

if(dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal")==F)
{dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal")}

if(dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25")==F)
{dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25")}

if(dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25/3cm")==F)
{dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25/3cm")}

if(dir.exists(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs",sep=""))==F)
{dir.create(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs",sep=""))}

if(dir.exists(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25",sep=""))==F)
{dir.create(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25",sep=""))}

if(dir.exists(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))==F)
{dir.create(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm",sep=""))}

if(dir.exists(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))==F)
{dir.create(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))}

#set.seed(8)
#init_vals<- c(runif(1,.02,.1),runif(1,.95,1.05),.2/(.2+1))

set.seed(2)
iv2<- c(runif(1,.02,.1),runif(1,.95,1.05),.2/(.2+1))
set.seed(3)
iv3<- c(runif(1,.02,.1),runif(1,.95,1.05),.2/(.2+1))
set.seed(4)
iv4<- c(runif(1,.02,.1),runif(1,.95,1.05),.2/(.2+1))
init_vals2<- c(iv2[1],iv4[2],.2/(.2+1))
init_vals3<- iv3
init_vals4<- c(iv4[1],iv2[2],.2/(.2+1))

init_vals<- init_vals4


#second 100k with nugget
#setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))
#load("vat_rw_MH_homMRnug.RData")
#init_vals<- c(res[nrow(res),]) #use last step of first Markov chain
#rm(res)

set.seed(51) #for first 100k steps
#set.seed(52) #for last 200k steps
pt<-proc.time() # Start Time


#output <- mh.alg(init = init_vals, n.sample = 201000)
#second 100000
output <- mh.alg(init = init_vals, n.sample = 100000)
#output <- mh.alg(init = init_vals, n.sample = 1000)


proposals<- output$proposals
res<- output$out

ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3] # End Time to be used in Effective Samples per Second Calculation

res<- res[-1,]

if(nug==FALSE){
  #setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25/3cm")
  setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))
  #save(res,file="vat_rw_MH_homMR.RData")
  #save(proposals,file="vatprops_rw_MH_homMR.RData")
  #save(ptFinal,file="vat_time_rw_MH_homMR.RData")
  #save(res,file="vat_rw_MH_homMR2.RData")
  #save(proposals,file="vatprops_rw_MH_homMR2.RData")
  #save(ptFinal,file="vat_time_rw_MH_homMR2.RData")
  
  #setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25/3cm")
  #setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))
  #load("vat_rw_MH_homMR.RData")
  #load("vatprops_rw_MH_homMR.RData")
  #load("vat_time_rw_MH_homMR.RData")
}

##################################with nugget###################################
if(nug==TRUE){
  #setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25/3cm")
  setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))
  save(res,file="vat_rw_MH_homMR_sqexp_nug.RData")
  save(proposals,file="vatprops_rw_MH_homMR_sqexp_nug.RData")
  save(ptFinal,file="vat_time_rw_MH_homMR_sqexp_nug.RData")
  #2nd 100k
  #save(res,file="vat_rw_MH_homMRnug2.RData")
  #save(proposals,file="vatprops_rw_MH_homMRnug2.RData")
  #save(ptFinal,file="vat_time_rw_MH_homMRnug2.RData")
  
  nPCs=18
  #setwd("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/.25/3cm")
  setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))
  load("vat_rw_MH_homMR_sqexp_nug.RData")
  load("vatprops_rw_MH_homMR_sqexp_nug.RData")
  load("vat_time_rw_MH_homMR_sqexp_nug.RData")
  
  #res1<- res
  #ptFinal1<- ptFinal
  
  #setwd(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/cal/",nPCs,"PCs/.25/3cm/chain4",sep=""))
  #load("vat_rw_MH_homMRnug2.RData")
  #load("vatprops_rw_MH_homMRnug2.RData")
  #load("vat_time_rw_MH_homMRnug2.RData")
  
  #res2<- res
  #ptFinal2<- ptFinal
  
  #res<- rbind(res1,res2)
}

################################################################################


plot(1:nrow(res),res[,1],main="Channel roughness",type="l")
plot(1:nrow(res),res[,2],main="River width error",type="l")
plot(1:nrow(res),res[,3],main="sigma^2",type="l")

length(unique(res[1:nrow(res),1]))/nrow(res)
length(unique(res[1:nrow(res),2]))/nrow(res)

sum(proposals[,4])/nrow(proposals); sum(proposals[,5])/nrow(proposals)

hist(res[,1],breaks=50,xlim=c(.02,.1))
abline(v=.0305,col="red")
hist(res[,2],breaks=50,xlim=c(.95,1.05))
abline(v=1,col="red")
hist(res[,3],breaks=50)

library(batchmeans)

for(i in 1:ncol(res)){
  print(i)
  print(paste("MR ESS is ",ess(res[,i],imse=FALSE),sep=""))
  print(bm(res[,i]))
}

#(ptFinal/(ess(res[,2],imse=FALSE)/4000))/(60^2)