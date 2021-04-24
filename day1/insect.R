library(TMB)
compile("insect.cpp")
dyn.load(dynlib("insect"))

#For data we use the built-in InsectSprays
par <- list()
par$logAlpha=rep(0,nlevels(InsectSprays$spray))
obj <- MakeADFun(InsectSprays, par, DLL="insect")
opt <- nlminb(obj$par, obj$fn, obj$gr)



map1 <- ...
obj1 <- ...
opt1 <- ...

#test the hypotesis
1-pchisq(2*(opt1$obj-opt$obj),2)

map2 <- ...
...
obj2 <- ...
opt2 <- ...

#test the hypotesis
1-pchisq(2*(opt2$obj-opt1$obj),1)