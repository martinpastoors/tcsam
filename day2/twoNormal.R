library(TMB)

setwd(file.path(getwd(),"tcsam", "day2"))
Sys.setenv(PATH=paste("C:/Rtools/bin", Sys.getenv("PATH"), sep =";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

compile("twoNormal.cpp")
dyn.load(dynlib("twoNormal"))
load("twoNormal.RData")

param <- list()
param$muRed <- c(0,0)
param$logSigmaRed <- c(0,0)
param$logitRhoRed <- 1
param$muBlue <- c(0,0)
param$logSigmaBlue <- c(0,0)
param$logitRhoBlue <- 1

obj <- MakeADFun(lda, param, DLL="twoNormal")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)


plot(rbind(t(lda$red),t(lda$blue)), type="n", xlab="x", ylab="y")
points(t(lda$red), col="red")
points(t(lda$blue), col="blue")
res<- obj$report()$res
points(t(lda$black), col=ifelse(res==0, "red", "blue"), pch=20)
