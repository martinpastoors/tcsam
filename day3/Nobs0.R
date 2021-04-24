load("Nobs.RData")

matplot(log(Nobs$Nobs), main="logN")
library(TMB)
compile("Nobs0.cpp")
dyn.load(dynlib("Nobs0"))

par <- list()
par$logsdR <- 0
par$logsdS <- 0
par$logsd <- 0
par$logN <- matrix(0, nrow=nrow(Nobs$Nobs), ncol=ncol(Nobs$Nobs))

obj <- MakeADFun(Nobs, par, random="logN", DLL="Nobs0")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
matplot(as.list(sdr,"Est")$logN, type="l", add=TRUE)
