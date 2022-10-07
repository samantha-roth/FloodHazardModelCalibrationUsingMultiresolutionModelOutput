#calculated the predictive mean and covariance for Z_R
#from HomMR approach

pred_mean_cov<- function(x.test){
  pred_mean<- rep(NA,nPCs)
  pred_cov<- matrix(0,nrow=nPCs,ncol=nPCs)
  
  ##scale test parameters
  x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
  
  #pt<-proc.time() # Start Time
  for(k in 1:nPCs){
    #pt<-proc.time() # Start Time
    y.c2 <- y.c2mat[,k]
    y.e2 <- y.e2mat[,k]
    
    #get emulator parameter estimates
    #load(paste("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/100e400c/downscale/",nPCs,"PCs/pars.homMLExp.PC",k,".RData",sep=""))
    load(paste(emulatorParsLoc,k,".RData",sep=""))
    
    precmat.ml<- precmat[,,k]
    
    H <- cbind(1, as.matrix(x.test))
    
    #squared exponential
    #crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be, power= 2)
    #exponential (with or without nugget)
    #crosscov.ml <- ml.crosscov(x.v2, x.c.std.rev, x.e.std.rev, pars, H.c, H.e, H, var.bc, var.be, covfunc=cf)
    crosscov.ml <- ml.crosscovFAST(x.v2, x.c.std.rev, x.e.std.rev, pars, H, covfunc=cf, k=k)
    
    #assume a priori E(y|x) = 0 for any x
    #mean.ml <- crosscov.ml %*% (precmat.ml %*%(c(y.c2, y.e2)))
    mean.ml <- crossprod(t(crosscov.ml), crossprod(t(precmat.ml),c(y.c2, y.e2)))
    
    #squared exponential
    #var.ml <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,power=2) + cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta,power=2)
    
    #exponential
    var.ml1 <- cov.x1.x2(x.v2, x.v2, pars$m_sigma, pars$m_theta,covfunc=cf) + 
      cov.x1.x2(x.v2, x.v2, pars$c_sigma*pars$rho, pars$c_theta, covfunc=cf)
    
    #var.ml <- var.ml1 + H %*% (pars$rho^2*var.bc + var.be) %*% t(H) - 
    #  crosscov.ml %*% precmat.ml %*% t(crosscov.ml)
    
    var.ml<- var.ml1 + quad.tform(pars$rho^2*var.bc + var.be, H) -
      quad.tform(precmat.ml, crosscov.ml)
    
    #no nugget
    if(nug==FALSE) var.ml.full <- var.ml
    
    #nugget
    if(nug==TRUE) var.ml.full <- var.ml + pars$m_nugget
    
    
    pred_mean[k]<- as.numeric(mean.ml)
    pred_cov[k,k]<- as.numeric(var.ml.full)
    #ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3]
  }
  #ptFinal<-proc.time()-pt ; ptFinal<-ptFinal[3]
  list(pred_mean=pred_mean,pred_cov=pred_cov)
}


ll_func<- function(mean, cov, sigma2){
  #-.5*log(det(cov+ sigma2*invtKK))-.5*t(Z_R-mean)%*%solve(cov+ sigma2*invtKK)%*%(Z_R-mean)
  
  -.5*log(det(cov+ sigma2*invtKK))-.5*quad.form(chol2inv(chol(cov+ sigma2*invtKK)),Z_R-mean)
}

logDetJ<- function(p,a=c(.02,.02,.95,-5),b=c(.1,.4,1.05,5)){
  -sum(log(p-a)+log(b-p))
}

#logit<- function(p){log(p/(1-p))}

UtoR<- function(p,a=c(.02,.02,.95,-5),b=c(.1,.4,1.05,5)){
  log((p-a)/(b-a))-log(1-((p-a)/(b-a)))
}

RtoU<- function(x,a=c(.02,.02,.95,-5),b=c(.1,.4,1.05,5)){
  (b-a)*(1/(1+exp(-x)))+a
}
