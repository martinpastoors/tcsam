load("Cobs.RData")
library(TMB)

compile("linkEx.cpp")
dyn.load(dynlib("link"))

par <- list()
par$logAlpha <- rep(0, length(unique(Cobs$aux[,3])))
par$logBeta <- rep(log(2), length(unique(Cobs$aux[,3])))

map=list(logAlpha=factor(rep(1,length(par$logAlpha))),
         logBeta=factor(rep(1,length(par$logBeta))))

obj <- MakeADFun(Cobs, par, DLL="link", map = map)
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)


matplot(rownames(Cobs$N), xtabs(log(Cobs$Cobs)~Cobs$aux[,1]+Cobs$aux[,3]), ylab="Log C", xlab="Year")
matplot(rownames(Cobs$N), xtabs(obj$report()$logPred~Cobs$aux[,1]+Cobs$aux[,3]), type="l", add=TRUE)

