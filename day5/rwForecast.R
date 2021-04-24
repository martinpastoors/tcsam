library(TMB)
compile("rwForecast.cpp")
dyn.load(dynlib("rwForecast"))
load("F.RData")

parameters <- list(
  logSdRw=0,
  logSdObs=0,
  lam=rep(0,length(data$y))
  )

obj <- MakeADFun(data,parameters,random="lam",DLL="rwForecast")

opt<-nlminb(obj$par,obj$fn,obj$gr)

sdr<-sdreport(obj)
pl <- as.list(sdr,"Est")
plsd <- as.list(sdr,"Std")

plot(data$y, xlim=c(0,50))
lines(pl$lam)
lines(pl$lam-2*plsd$lam, lty="dotted")
lines(pl$lam+2*plsd$lam, lty="dotted")

nsim <- 1000
nYears = 6 #Forecast years added base year
sim <- matrix(NA, ncol=nsim, nrow=nYears)
set.seed(123)
sim[1,] <- rnorm(nsim, pl$lam[length(data$y)], ...) #Simulate in base year
for(y in 2:nYears){
  sim[y,] <- rnorm(nsim, ... , ...)#Simulate the process in future years
}

#Plot future realisations
plot(data$y, xlim=c(0,50),main = "Forecast",ylab = "F",cex.main = 2,
     cex.lab = 1.6)
lines(pl$lam)
lines(pl$lam-2*plsd$lam, lty="dotted")
lines(pl$lam+2*plsd$lam, lty="dotted")
matplot((length(data$y)):(length(data$y) + nYears-1), sim, type="l", add=TRUE, col=gray(.8))
matplot((length(data$y)):(length(data$y) + nYears-1), t(apply(sim, 1, quantile, c(0.025, 0.5, 0.975))), type="l", col="red", add=TRUE, lwd=3,lty = c(3,1,3))
points((length(data$y) +1):(length(data$y) + nYears-1), data$futureObs)

