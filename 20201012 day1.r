# install.packages("TMB")
# install.packages("pkgbuild")
# pkgbuild::has_rtools()
# pkgbuild::rtools_path()

# download.file("https://github.com/kaskr/adcomp/archive/master.zip", "adcomp-master.zip")
# unzip("adcomp-master.zip")

# writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# Sys.which("make")   # check if it works; it does

# devtools::install_github('fishfollower/SAM/stockassessment')
library(stockassessment)
library(TMB)

setwd("D:/GIT/tcsam/day1")

# need to set the system setenv link to Rtools
Sys.setenv(PATH=paste("C:/Rtools/bin", Sys.getenv("PATH"), sep =";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")


example(sam.fit)
fit


# C++ example of optimization

compile("a1.cpp")              # compile c code
dyn.load(dynlib("a1"))      # link c code to R

dat <- list()
param <- list(x=20)

obj <- MakeADFun(dat, param, DLL="a1")
obj$fn()     # function output
obj$gr(42)   # gradient

fit <- nlminb(obj$par, obj$fn)           # without gradient
fit <- nlminb(obj$par, obj$fn, obj$gr)   # with gradient

# Exercise 2
data(nscodData)
data(nscodConf)
data(nscodParameters)
fit <- sam.fit(nscodData, nscodConf, nscodParameters)
fc<-forecast(fit, fscale=c(1,1,1,1))
ssbplot(fc)
attributes(fit)

for (i in names(fit)) {
  print(i)
  print(summary(fit[[i]]))
}

sdreport(fit$obj)
as.list(sdreport(fit$obj),"Est",report= TRUE)
as.list(sdreport(fit$obj),"Std",report= TRUE)

# Beverton-Holt example

dat<-read.table("bh.dat", header=TRUE)
compile("bh.cpp")
dyn.load(dynlib("bh"))

data <- list(SSB=dat$SSB,logR=dat$logR)
parameters <- list(
  logA=0,
  logB=0,
  logSigma=0
)

obj <- MakeADFun(data,parameters,DLL="bh")
opt <- nlminb(obj$par,obj$fn,obj$gr)
rep <- sdreport(obj)

plot(dat$SSB,dat$logR, main = "Beverton-Holt", cex.main = 2,cex.lab = 1.5,
     xlab = "SSB",ylab = "logR")
lines(dat$SSB, rep$value, lwd=3)
lines(dat$SSB, rep$value+2*rep$sd, lty="dotted")
lines(dat$SSB, rep$value-2*rep$sd, lty="dotted")

obj$fn()
obj$gr()
obj$gr(opt$par)  # the real minimum; shows derivates to logA, logB and logSigma


rl  <- as.list(sdreport(obj), "Est", report=TRUE)
rl$pred
rsd <- as.list(sdreport(obj), "Std", report=TRUE)
rsd$pred

# Insect spray example

compile("insect.cpp")
dyn.load(dynlib("insect"))

#For data we use the built-in InsectSprays
par <- list()
par$logAlpha=rep(0,nlevels(InsectSprays$spray))
obj <- MakeADFun(InsectSprays, par, DLL="insect")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par

# InsectSprays$spray

map1 <- list(logAlpha=factor(c(1,1,3,4,5,1)))
obj1 <- MakeADFun(InsectSprays, par, DLL="insect", map=map1)
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

#test the hypotesis
1-pchisq(2*(opt1$obj-opt$obj),2)

map2 <- list(logAlpha=factor(c(NA,NA,1,2,3,NA)))
par2 <- par
par2$logAlpha[c(1,2,6)] <- 15
obj2 <- MakeADFun(InsectSprays, par2, DLL="insect", map=map2)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)

#test the hypotesis
1-pchisq(2*(opt2$obj-opt1$obj),1)

# Lecture 3 - FSA Fish Stock Assessment ===================

load("fsa.RData") # gets "dat"
# min(dat$age)
# max(dat$year)
# cbind(dat$year, dat$fleet, dat$age, dat$obs)

compile("fsa.cpp")
dyn.load(dynlib("fsa"))

parameters <- list(
  logN1Y=rep(0,nrow(dat$M)),
  logN1A=rep(0,ncol(dat$M)-1),
  logFY=rep(0,ncol(dat$M)),
  logFA=rep(0,nrow(dat$M)),
  logVarLogCatch=0, 
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logVarLogSurvey=0
)

# FA initialized at 0; F5-7 are at NA. 
# obj <- MakeADFun(dat,parameters,DLL="fsa2", map=list(logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE)
obj <- MakeADFun(dat,parameters,DLL="fsa", map=list(logFA=factor(c(1,1,2,2,2,NA,NA))), silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=1000,eval.max=1000))
rep <- sdreport(obj)
ssb <- rep$value[names(rep$value)=="ssb"]
ssb.sd <- rep$sd[names(rep$value)=="ssb"]


plot(ssb, type="l", lwd=5, col="red", ylim=c(0,550000))
lines(ssb-2*ssb.sd, type="l", lwd=1, col="red")
lines(ssb+2*ssb.sd, type="l", lwd=1, col="red")

# plot F
pl <- as.list(rep, "Est", report=TRUE)
matplot(t(pl$F), type="l")


# 2 plus group fsa2.cpp

load("fsa.RData") # gets "dat"

compile("fsa2.cpp")
dyn.load(dynlib("fsa2"))

parameters2 <- list(
  logN1Y=rep(0,nrow(dat$M)),
  logN1A=rep(0,ncol(dat$M)-1),
  logFY=rep(0,ncol(dat$M)),
  logFA=rep(0,nrow(dat$M)),
  logVarLogCatch=0, 
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logVarLogSurvey=0
)

# FA initialized at 0; F5-7 are at NA. 
obj2 <- MakeADFun(dat,parameters2,DLL="fsa2", map=list(logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE)

opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, control=list(iter.max=1000,eval.max=1000))
rep2 <- sdreport(obj2)
ssb2 <- rep2$value[names(rep2$value)=="ssb"]
ssb2.sd <- rep2$sd[names(rep2$value)=="ssb"]

plot(ssb2, type="l", lwd=5, col="red", ylim=c(0,550000))
lines(ssb2-2*ssb2.sd, type="l", lwd=1, col="red")
lines(ssb2+2*ssb2.sd, type="l", lwd=1, col="red")

# plot F
pl2 <- as.list(rep2, "Est", report=TRUE)
matplot(t(pl2$F), type="l")

# 3 variances for catch fsa3.cpp
compile("fsa3.cpp")
dyn.load(dynlib("fsa3"))

parameters3 <- list(
  logN1Y=rep(0,nrow(dat$M)),
  logN1A=rep(0,ncol(dat$M)-1),
  logFY=rep(0,ncol(dat$M)),
  logFA=rep(0,nrow(dat$M)),
  logVarLogCatch=rep(0,nrow(dat$M)), 
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logVarLogSurvey=0
)

# FA initialized at 0; F5-7 are at NA. 
obj3 <- MakeADFun(dat,parameters3,DLL="fsa3", map=list(logFA=factor(c(1:4,NA,NA,NA)),
                                            logVarLogCatch=factor(c(1,2,2,2,2,2,2))), silent=TRUE)

opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr, control=list(iter.max=1000,eval.max=1000))
rep3 <- sdreport(obj3)
ssb3 <- rep3$value[names(rep3$value)=="ssb"]
ssb3.sd <- rep3$sd[names(rep3$value)=="ssb"]

plot(ssb3, type="l", lwd=5, col="red", ylim=c(0,550000))
lines(ssb3-2*ssb3.sd, type="l", lwd=1, col="red")
lines(ssb3+2*ssb3.sd, type="l", lwd=1, col="red")

# plot F
pl3 <- as.list(rep, "Est", report=TRUE)
matplot(t(pl3$F), type="l")

# 4 survey catchabilities fsa4.cpp

compile("fsa4.cpp")
dyn.load(dynlib("fsa4"))

parameters4 <- list(
  logN1Y=rep(0,nrow(dat$M)),
  logN1A=rep(0,ncol(dat$M)-1),
  logFY=rep(0,ncol(dat$M)),
  logFA=rep(0,nrow(dat$M)),
  logVarLogCatch=rep(0,nrow(dat$M)), 
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logVarLogSurvey=rep(0,ncol(dat$M))
)

fy <- min(dat$year)
ly <- max(dat$year)

# FA initialized at 0; F5-7 are at NA. 
obj4 <- MakeADFun(dat,parameters4,DLL="fsa4", 
                 map=list(logFA=factor(c(1:4,NA,NA,NA)),
                          logVarLogCatch=factor(c(1,2,2,2,2,2,2)),
                          logVarLogSurvey=factor(c(rep(1,2000-fy+1), 
                                                   rep(2,ly-2000)))), 
                 silent=TRUE)

opt4 <- nlminb(obj4$par, obj4$fn, obj4$gr, control=list(iter.max=1000,eval.max=1000))
#opt$par

rep4 <- sdreport(obj4)
ssb4 <- rep$value[names(rep4$value)=="ssb"]
ssb4.sd <- rep$sd[names(rep4$value)=="ssb"]

plot(ssb4, type="l", lwd=5, col="red", ylim=c(0,550000))
lines(ssb4-2*ssb4.sd, type="l", lwd=1, col="red")
lines(ssb4+2*ssb4.sd, type="l", lwd=1, col="red")

# plot F
pl4 <- as.list(rep, "Est", report=TRUE)
matplot(t(pl4$F), type="l")


plot(ssb, type="l", lwd=5, col="red", ylim=c(0,550000))
lines(ssb2, type="l", lwd=1, col="blue")
lines(ssb3, type="l", lwd=1, col="green")
lines(ssb4, type="l", lwd=1, col="gray")


L0 <- opt$objective
L2 <- opt2$objective
L3 <- opt3$objective
L4 <- opt4$objective

