# install.packages("TMB")
# install.packages("pkgbuild")
# pkgbuild::has_rtools()
# pkgbuild::rtools_path()

# download.file("https://github.com/kaskr/adcomp/archive/master.zip", "adcomp-master.zip")
# unzip("adcomp-master.zip")

# devtools::install_github('fishfollower/SAM/stockassessment')
library(stockassessment)
library(TMB)

setwd("D:/GIT/tcsam/day3")

Sys.setenv(PATH=paste("C:/Rtools/bin", Sys.getenv("PATH"), sep =";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

cn <- read . ices (" testdata /cn. dat ")
cw <- read . ices (" testdata /cw. dat ")
dw <- read . ices (" testdata /dw. dat ")
lf <- read . ices (" testdata /lf. dat ")
lw <- read . ices (" testdata /lw. dat ")
mo <- read . ices (" testdata /mo. dat ")
nm <- read . ices (" testdata /nm. dat ")
pf <- read . ices (" testdata /pf. dat ")
pm <- read . ices (" testdata /pm. dat ")
sw <- read . ices (" testdata /sw. dat ")
surveys <- read . ices (" testdata / survey .dat ")
dat <-setup .sam. data ( surveys = surveys , residual . fleet =cn , prop . mature =mo , stock . mean . weight =sw ,
                         catch . mean . weight =cw , dis . mean . weight =dw , land . mean . weight =lw , prop .f=pf ,
                         prop .m=pm , natural . mortality =nm , land . frac =lf)
conf <- stockassessment::defcon (dat)  # define configuration
par <- stockassessment::defpar (dat , conf ) # setup parameter file
fit <- stockassessment::sam.fit(dat ,conf , par ) # fit the model

stockassessment::saveConf(conf , file ="model.cfg")
stockassessment::saveConf(conf , file ="model.cfg", overwrite=TRUE)

fitweb <- stockassessment::fitfromweb("WGWIDE2020_MAC")
fit    <- get(load("D:/iWGWIDE/2020/06. Data/mac.27.nea/output/model fit.RData"))

stockassessment::modeltable(c(fit,fitweb))

fit$obj$par
names(fitweb$obj$par)

# N process exercises

# -----------------------------------------------------------
# Recruitment modelling
# -----------------------------------------------------------

rm(list=ls())

compile("Robs0.cpp")
dyn.load(dynlib("Robs0"))

load("Robs.RData")
plot(Robs$year, log(Robs$Robs))

Robs$mode <- 0  # RW

par <- list()
par$logsdo <- 0 # log sd observations
par$logsdp <- 0  # log sd process error
par$logR <- rep(0,length(Robs$Robs))
par$rickerpar <- if(Robs$mode==1){numeric(2)}else{numeric(0)}
par$bhpar <- if(Robs$mode==2){numeric(2)}else{numeric(0)}

obj <- MakeADFun(Robs, par, random="logR", DLL="Robs0")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
logR<-as.list(sdr, "Est")$logR
sdlogR<-as.list(sdr, "Std")$logR

lines(Robs$year, logR, col="black")
lines(Robs$year, logR-2*sdlogR, col="black", lty="dashed")
lines(Robs$year, logR+2*sdlogR, col="black", lty="dashed")

# Ricker

Robs$mode <- 1  # Ricker

par <- list()
par$logsdo <- 0 # log sd observations
par$logsdp <- 0  # log sd process error
par$logR <- rep(0,length(Robs$Robs))
par$rickerpar <- if(Robs$mode==1){numeric(2)}else{numeric(0)}
par$bhpar <- if(Robs$mode==2){numeric(2)}else{numeric(0)}

obj <- MakeADFun(Robs, par, random="logR", DLL="Robs0")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
logR<-as.list(sdr, "Est")$logR
sdlogR<-as.list(sdr, "Std")$logR
lines(Robs$year, logR, col="red")

# BH

Robs$mode <- 2  # BH

par <- list()
par$logsdo <- 0 # log sd observations
par$logsdp <- 0  # log sd process error
par$logR <- rep(0,length(Robs$Robs))
par$rickerpar <- if(Robs$mode==1){numeric(2)}else{numeric(0)}
par$bhpar <- if(Robs$mode==2){numeric(2)}else{numeric(0)}

obj <- MakeADFun(Robs, par, random="logR", DLL="Robs0")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
logR<-as.list(sdr, "Est")$logR
sdlogR<-as.list(sdr, "Std")$logR

lines(Robs$year, logR, col="blue")

# =========================================================================================

# plot of random walk
plot(cumsum(rnorm(100)), type="l", lwd=50)

# plot of white noise
plot(rnorm(100), type="l", lwd=50)

# =========================================================================================
# Survival exercise
# =========================================================================================

rm(list=ls())

load("Nobs.RData")
matplot(log(Nobs$Nobs), main="logN")
# str(Nobs)

compile("Nobs0.cpp")
dyn.load(dynlib("Nobs0"))

par <- list()
par$logsdR <- 0
par$logsdS <- 0
par$logsd <- 0
par$logN <- matrix(0, nrow=nrow(Nobs$Nobs), ncol=ncol(Nobs$Nobs))

obj <- MakeADFun(Nobs, par, random="logN", DLL="Nobs0")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
matplot(as.list(sdr,"Est")$logN, type="l", add=TRUE)

# =========================================================================================
# Fishing mortality exercise
# =========================================================================================

rm(list=ls())

load("Fobs.RData")

matplot(Fobs$year, log(Fobs$Fobs), xlab="Year", ylab="logF", pch=colnames(Fobs$Fobs))

compile("Fobs.cpp")
dyn.load(dynlib("Fobs"))

Fobs$cormode <- 2

par <- list()
par$logsdF <- rep(0,ncol(Fobs$Fobs))
par$transPsi <- if(Fobs$cormode==0){numeric(0)}else{0.1}
par$logsd <- 0
par$logF <- matrix(0, nrow=nrow(Fobs$Fobs), ncol=ncol(Fobs$Fobs))

map=list(logsdF=factor(rep(1,ncol(Fobs$Fobs))))
obj <- MakeADFun(Fobs, par, random="logF", DLL="Fobs", map=map)
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)

matplot(Fobs$year, as.list(sdr,"Est")$logF, type="l", add=TRUE)
matplot(Fobs$year, as.list(sdr,"Est")$logF, type="l", col="blue", add=TRUE)
matplot(Fobs$year, as.list(sdr,"Est")$logF, type="l", col="red", add=TRUE)


Fobs$cormode <- 0
par <- list()
par$logsdF <- rep(0,ncol(Fobs$Fobs))
par$transPsi <- if(Fobs$cormode==0){numeric(0)}else{0.1}
par$logsd <- 0
par$logF <- matrix(0, nrow=nrow(Fobs$Fobs), ncol=ncol(Fobs$Fobs))
map=list(logsdF=factor(rep(1,ncol(Fobs$Fobs))))
obj0 <- MakeADFun(Fobs, par, random="logF", DLL="Fobs", map=map)
fit0 <- nlminb(obj0$par, obj0$fn, obj0$gr)
sdr0<-sdreport(obj0)

Fobs$cormode <- 1
par <- list()
par$logsdF <- rep(0,ncol(Fobs$Fobs))
par$transPsi <- if(Fobs$cormode==0){numeric(0)}else{0.1}
par$logsd <- 0
par$logF <- matrix(0, nrow=nrow(Fobs$Fobs), ncol=ncol(Fobs$Fobs))
map=list(logsdF=factor(rep(1,ncol(Fobs$Fobs))))
obj1 <- MakeADFun(Fobs, par, random="logF", DLL="Fobs", map=map)
fit1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
sdr1<-sdreport(obj1)

AIC0 <- 2*fit0$obj + 2*length(fit0$par)
AIC1 <- 2*fit1$obj + 2*length(fit1$par)


# Practical exercise

setwd("D:/GIT/mptools/tcsam/day3/PracticalExercise")

cn<-read.ices("cn.dat")
cw<-read.ices("cw.dat")
dw<-read.ices("dw.dat")
lf<-read.ices("lf.dat")
lw<-read.ices("lw.dat")
mo<-read.ices("mo.dat")
nm<-read.ices("nm.dat")
pf<-read.ices("pf.dat")
pm<-read.ices("pm.dat")
sw<-read.ices("sw.dat")
surveys<-read.ices("survey.dat")

dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

conf0<-loadConf(dat,"model.cfg")
par0<-defpar(dat,conf0)
fit0<-sam.fit(dat,conf0,par0)

conf1<-loadConf(dat,"model.cfg") 
conf1$corFlag <- 1
par1<-defpar(dat,conf1)
fit1<-sam.fit(dat,conf1,par1)

conf2<-loadConf(dat,"model.cfg") 
conf2$corFlag <- 2
# saveConf(conf2,"model1,cfg")
par2<-defpar(dat,conf2)
fit2<-sam.fit(dat,conf2,par2)

conf3<- conf2 
par3<-defpar(dat,conf3)
par3$itrans_rho <- 7
nages <- conf3$maxAge - conf3$minAge
map <- list()
map$itrans_rho <- as.factor(NA)
fit3<-sam.fit(dat,conf3,par3, map=map)

modeltable(c(fit0, fit1,fit2, fit3))

ssbplot(c(fit0, fit2, fit3))
fbarplot(c(fit0,fit2,fit3))

fselectivityplot(fit2)
fselectivityplot(fit3)

matplot(faytable(fit3), type="b")
