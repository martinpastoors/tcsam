load("allfleetsblock.RData")

library(TMB)
compile("allfleetsblock.cpp")
dyn.load(dynlib("allfleetsblock"))

allfleetsblock$keyQ <- rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,NA),
                        c(NA, 0, 1, 2, 3, 4, 5, 6,NA),
                        c( 7, 8, 9,10,11,12,NA,NA,NA))

allfleetsblock$keySd <- rbind(c( 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         c(NA, 1, 1, 1, 1, 1, 1, 1,NA),
                         c( 2, 2, 2, 2, 2, 2,NA,NA,NA))

par <- list()
par$logQ <- numeric(max(allfleetsblock$keyQ, na.rm=TRUE)+1)
par$logsd <- numeric(max(allfleetsblock$keySd, na.rm=TRUE)+1)
par$missing <- numeric(sum(is.na(allfleetsblock$obs)))

obj <- MakeADFun(allfleetsblock, par, random="missing", DLL="allfleetsblock")
fit <- nlminb(obj$par, obj$fn, obj$gr)
est <- obj$report()$logPred

par(mfrow=c(1,3))
for(f in 1:3){
  idx<-which(allfleetsblock$aux[,2]==f)
  matplot(xtabs(log(allfleetsblock$obs[idx])~allfleetsblock$aux[idx,1]+allfleetsblock$aux[idx,3]), ylab="Log Obs")
  matplot(xtabs(est[idx]~allfleetsblock$aux[idx,1]+allfleetsblock$aux[idx,3]), type="l", add=TRUE)
}
