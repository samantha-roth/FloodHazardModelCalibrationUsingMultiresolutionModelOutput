#plot where we have expensive parameter samples from LHS

rm(list=ls())

nRuns10m<- 50
nRuns50m<- 200

#load parameters
load("/storage/work/svr5482/FloodingModelCalibrationProject/parameterSamples/50e200c/nChRWE_lhs_samples_allU.RData")

parVals10m<- data.frame("run"=1:nRuns10m, "n_ch"= samp.E[,1], "rwe"= samp.E[,2])
parVals50m<- data.frame("run"=1:nRuns50m, "n_ch"= samp.C[,1], "rwe"= samp.C[,2])

plot(parVals10m$n_ch,parVals10m$rwe)

plot(parVals50m$n_ch,parVals50m$rwe)

indsChless.028<- which(parVals10m$n_ch<.028)
indsRWEless.96<- which(parVals10m$rwe<.96)
indsChgreat.092<- which(parVals10m$n_ch>.092)
indsRWEgreat1.04<- which(parVals10m$rwe>1.04)

union(indsChless.028,indsRWEless.96)
union(indsChless.028,indsRWEgreat1.04)
union(indsChgreat.092,indsRWEless.96)
union(indsChgreat.092,indsRWEgreat1.04)

lolo<- union(indsChless.028,indsRWEless.96)
lohi<- union(indsChless.028,indsRWEgreat1.04)
hilo<- union(indsChgreat.092,indsRWEless.96)
hihi<- union(indsChgreat.092,indsRWEgreat1.04)

parVals50m[lolo,]==parVals10m[lolo,] 
parVals50m[lohi,]==parVals10m[lohi,]
parVals50m[hilo,]==parVals10m[hilo,]
parVals50m[hihi,]==parVals10m[hihi,]


indsChless.028_50m<- which(parVals50m$n_ch<.028)
indsRWEless.96_50m<- which(parVals50m$rwe<.96)
indsChgreat.092_50m<- which(parVals50m$n_ch>.092)
indsRWEgreat1.04_50m<- which(parVals50m$rwe>1.04)

union(indsChless.028_50m,indsRWEless.96_50m)
union(indsChless.028_50m,indsRWEgreat1.04_50m)
union(indsChgreat.092_50m,indsRWEless.96_50m)
union(indsChgreat.092_50m,indsRWEgreat1.04_50m)

lolo_50m<- union(indsChless.028_50m,indsRWEless.96_50m)
lohi_50m<- union(indsChless.028_50m,indsRWEgreat1.04_50m)
hilo_50m<- union(indsChgreat.092_50m,indsRWEless.96_50m)
hihi_50m<- union(indsChgreat.092_50m,indsRWEgreat1.04_50m)

#look for good lohi values
plot(parVals50m$n_ch[lohi_50m],parVals50m$rwe[lohi_50m])
plot(parVals50m$n_ch[lohi],parVals50m$rwe[lohi])

#also use 69, 92


lohi_50mNot10m<- lohi_50m[which(lohi_50m>100)]

plot(parVals50m$n_ch[lohi_50mNot10m],parVals50m$rwe[lohi_50mNot10m])
plot(parVals50m$n_ch[c(99,lohi_50mNot10m)],parVals50m$rwe[c(99,lohi_50mNot10m)])

parVals10m[lohi,]
parVals50m[lohi,] #use run 10
parVals50m[lohi_50m,] #use run 10

hilo_50mNot10m<- hilo_50m[which(hilo_50m>100)]

#look for good hilo values
plot(parVals50m$n_ch[hilo_50m],parVals50m$rwe[hilo_50m])
plot(parVals50m$n_ch[hilo],parVals50m$rwe[hilo])

plot(parVals50m$n_ch[hilo_50mNot10m],parVals50m$rwe[hilo_50mNot10m])
plot(parVals50m$n_ch[c(2,hilo_50mNot10m)],parVals50m$rwe[c(2,hilo_50mNot10m)])

parVals50m[hilo,] #use run 2

parVals50m[hilo_50m,]


if(dir.exists("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c")==F){
  dir.create("/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c")
}
save(lolo,lohi,hilo,hihi,lolo_50m,lohi_50m,hilo_50m,hihi_50m,
     file="/storage/work/svr5482/FloodingModelCalibrationProject/multires/spatial/changPC2014/10m50m/outputData/nCh_RWE/lohi/50e200c/5thExtremes.RData")


indsChless.024<- which(parVals10m$n_ch<.024)
indsRWEless.955<- which(parVals10m$rwe<.955)
indsChgreat.096<- which(parVals10m$n_ch>.096)
indsRWEgreat1.045<- which(parVals10m$rwe>1.045)

union(indsChless.024,indsRWEless.955)
union(indsChless.024,indsRWEgreat1.045)
union(indsChgreat.096,indsRWEless.955)
union(indsChgreat.096,indsRWEgreat1.045)

length(union(union(indsChless.024,indsRWEless.955),union(indsChgreat.096,indsRWEgreat1.045)))


indsChless.024_50m<- which(parVals50m$n_ch<.024)
indsRWEless.955_50m<- which(parVals50m$rwe<.955)
indsChgreat.096_50m<- which(parVals50m$n_ch>.096)
indsRWEgreat1.045_50m<- which(parVals50m$rwe>1.045)


union(indsChless.024_50m,indsRWEless.955_50m)
union(indsChless.024_50m,indsRWEgreat1.045_50m)
union(indsChgreat.096_50m,indsRWEless.955_50m)
union(indsChgreat.096_50m,indsRWEgreat1.045_50m)

length(union(union(indsChless.024_50m,indsRWEless.955_50m),union(indsChgreat.096_50m,indsRWEgreat1.045_50m)))