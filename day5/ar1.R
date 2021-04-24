load("cpue.RData")

library(TMB)
compile("ar1.cpp")
dyn.load(dynlib("ar1"))

data = list(y = y)
par = list(logSigma = -2, phiTrans = 1, gamma = rep(0,length(y)))

obj <- MakeADFun(data,par,random = "gamma",DLL = "ar1")
fit <- nlminb(obj$par,obj$fn,obj$gr)
