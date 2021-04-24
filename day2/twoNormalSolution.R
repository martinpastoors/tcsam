library(TMB)
compile("lda.cpp")
dyn.load(dynlib("lda"))
load("lda.RData")

param <- list()
param$muRed <- c(0,0)
param$logSigmaRed <- c(0,0)
param$logitRhoRed <- 1
param$muBlue <- c(0,0)
param$logSigmaBlue <- c(0,0)
param$logitRhoBlue <- 1

obj <- MakeADFun(lda, param, DLL="lda")
opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj))


plot(rbind(t(lda$red),t(lda$blue)), type="n", xlab="x", ylab="y")
points(t(lda$red), col="red")
points(t(lda$blue), col="blue")
res<-obj$report()$res
points(t(lda$black), col=ifelse(res==0, "red", "blue"), pch=20)
