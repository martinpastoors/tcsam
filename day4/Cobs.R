load("Cobs.RData")
library(TMB)

compile("Cobs.cpp")
dyn.load(dynlib("Cobs"))

par <- list()
par$logsd <- rep(0, length(unique(Cobs$aux[,3])))

map = list()
map$logsd = as.factor(rep(1,length(par$logsd)))
obj <- MakeADFun(Cobs, par, DLL="Cobs",map = map)
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)

matplot(rownames(Cobs$N), xtabs(log(Cobs$Cobs)~Cobs$aux[,1]+Cobs$aux[,3]), ylab="Log C", xlab="Year")
matplot(rownames(Cobs$N), xtabs(obj$report()$logPred~Cobs$aux[,1]+Cobs$aux[,3]), type="l", add=TRUE)


