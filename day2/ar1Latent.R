
rm(list=ls())

#Estimate model
library(TMB)
compile("ar1Latent.cpp")
dyn.load("ar1Latent")

load("cpue.RData")
data = list(y = y)
par = list(logSigma = -2,
           phiTrans = 1,
           gamma = rep(0,length(y)))

obj = MakeADFun(data,par,random = "gamma",DLL = "ar1Latent")
opt = nlminb(obj$par,obj$fn,obj$gr,control = list(trace = 1))
rep = sdreport(obj,getJointPrecision = TRUE)


library(SparseM)
nameIndex = which(colnames(rep$jointPrecision)=="gamma")
Q = rep$jointPrecision[nameIndex,nameIndex]
Q[Q!=0] = 1  # We are only interested if they are 0 or not
image(Q, main = "Sparsness structure of AR(1)")

