load("Robs.RData")

library(TMB)
compile("Robs0.cpp")
dyn.load(dynlib("Robs0"))

plot(Robs$year, log(Robs$Robs))

# RW

Robs$mode <- 0
par <- list()
par$logsdo <- 0
par$logsdp <- 0
par$logR <- rep(0,length(Robs$Robs))

obj <- MakeADFun(Robs, par, random="logR", DLL="Robs0")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
logR<-as.list(sdr, "Est")$logR
sdlogR<-as.list(sdr, "Std")$logR
lines(Robs$year, logR, col="black")
lines(Robs$year, logR-2*sdlogR, col="black", lty="dashed")
lines(Robs$year, logR+2*sdlogR, col="black", lty="dashed")
