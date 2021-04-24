library(TMB)
compile("ar1Solution.cpp")
dyn.load(dynlib("ar1Solution"))
x<-scan("ar1.dat", quiet=TRUE)

data <- list(x=x, code=0)
parameters <- list(
  logSigma=0,
  itPhi=1
  )

data$code=0
obj <- MakeADFun(data,parameters,DLL="ar1Solution")
opt<-nlminb(obj$par,obj$fn,obj$gr)
c(opt$obj, opt$par)

data$code=1
obj2 <- MakeADFun(data,parameters,DLL="ar1Solution")
opt2<-nlminb(obj2$par,obj2$fn,obj2$gr)
c(opt2$obj, opt2$par)

data$code=2
obj3 <- MakeADFun(data,parameters,DLL="ar1Solution")
opt3<-nlminb(obj3$par,obj3$fn,obj3$gr)
c(opt3$obj, opt3$par)

data$code=3
obj4 <- MakeADFun(data,parameters,DLL="ar1Solution")
opt4<-nlminb(obj4$par,obj4$fn,obj4$gr)
c(opt4$obj, opt4$par)



