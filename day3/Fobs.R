load("Fobs.RData")

matplot(Fobs$year, log(Fobs$Fobs), xlab="Year", ylab="logF", pch=colnames(Fobs$Fobs))

library(TMB)
compile("Fobs.cpp")
dyn.load(dynlib("Fobs"))

Fobs$cormode <- 0

par <- list()
par$logsdF <- rep(0,ncol(Fobs$Fobs))
par$transPsi <- if(Fobs$cormode==0){numeric(0)}else{0.1}
par$logsd <- 0
par$logF <- matrix(0, nrow=nrow(Fobs$Fobs), ncol=ncol(Fobs$Fobs))

map=list(logsdF=factor(rep(1,ncol(Fobs$Fobs))))
obj <- MakeADFun(Fobs, par, random="logF", DLL="Fobs", map=map)
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
matplot(Fobs$year, as.list(sdr,"Est")$logF, type="l", add=TRUE)

