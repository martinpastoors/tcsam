load("allfleets.RData")

library(TMB)
compile("allfleets0.cpp")
dyn.load(dynlib("allfleets0"))

allfleets$keyQ <- rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,NA),
                        c(NA, 0, 1, 2, 3, 4, 5, 6,NA),
                        c( 7, 8, 9,10,11,12,NA,NA,NA))

allfleets$keySd <- rbind(c( 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         c(NA, 1, 1, 1, 1, 1, 1, 1,NA),
                         c( 2, 2, 2, 2, 2, 2,NA,NA,NA))

par <- list()
par$logQ <- numeric(max(allfleets$keyQ, na.rm=TRUE)+1)
par$logsd <- numeric(max(allfleets$keySd, na.rm=TRUE)+1)
par$missing <- numeric(sum(is.na(allfleets$obs)))

obj <- MakeADFun(allfleets, par, random="missing", DLL="allfleets0")
fit <- nlminb(obj$par, obj$fn, obj$gr)
est <- obj$report()$logPred

