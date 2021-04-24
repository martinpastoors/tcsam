# install.packages("TMB")
# install.packages("pkgbuild")
# pkgbuild::has_rtools()
# pkgbuild::rtools_path()

# download.file("https://github.com/kaskr/adcomp/archive/master.zip", "adcomp-master.zip")
# unzip("adcomp-master.zip")

# devtools::install_github('fishfollower/SAM/stockassessment')
library(stockassessment)
library(TMB)


setwd("D:/GIT/tcsam/day4")

Sys.setenv(PATH=paste("C:/Rtools/bin", Sys.getenv("PATH"), sep =";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# Catch observations

rm(list=ls())

load("Cobs.RData")

compile("Cobs.cpp")
dyn.load(dynlib("Cobs"))

par <- list()
par$logsd <- rep(0, length(unique(Cobs$aux[,3])))

map = list()
map$logsd = as.factor(rep(1,length(par$logsd)))
obj <- MakeADFun(Cobs, par, DLL="Cobs",map = map)
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)

# first plot the observations
matplot(rownames(Cobs$N), xtabs(log(Cobs$Cobs)~Cobs$aux[,1]+Cobs$aux[,3]), ylab="Log C", xlab="Year")
# then add the fitted values
matplot(rownames(Cobs$N), xtabs(obj$report()$logPred~Cobs$aux[,1]+Cobs$aux[,3]), type="l", add=TRUE)

aic1 <- 2*fit$obj + 2*length(fit$par)

map = list()
map$logsd = as.factor(1:length(par$logsd))
obj2 <- MakeADFun(Cobs, par, DLL="Cobs",map = map)
fit2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
sdr<-sdreport(obj2)

aic2 <- 2*fit$obj2 + 2*length(fit2$par)
aic1

# =========================================================================
# variance difference by year
# =========================================================================

# linkEx.R

rm(list=ls())

load("Cobs.RData")

compile("linkEx.cpp")
dyn.load(dynlib("linkEx"))

par <- list()
par$logAlpha <- rep(0, length(unique(Cobs$aux[,3])))
par$logBeta <- rep(log(2), length(unique(Cobs$aux[,3])))

map=list(logAlpha=factor(rep(1,length(par$logAlpha))),
         logBeta=factor(rep(1,length(par$logBeta))))

obj <- MakeADFun(Cobs, par, DLL="linkEx", map = map)
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)


matplot(rownames(Cobs$N), xtabs(log(Cobs$Cobs)~Cobs$aux[,1]+Cobs$aux[,3]), ylab="Log C", xlab="Year")
matplot(rownames(Cobs$N), xtabs(obj$report()$logPred~Cobs$aux[,1]+Cobs$aux[,3]), type="l", add=TRUE)

fit$objective
# -129.2638
rl <- as.list(sdreport(obj), "Est", report=TRUE)

# =========================================================================
# Practical with variance covariance matrix
# =========================================================================

rm(list=ls())

setwd("D:/GIT/mptools/tcsam/day4/practical obsCatch")

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

conf<-loadConf(dat,"model.cfg", patch=TRUE)
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

ages    <- conf$minAge:conf$maxAge
fleets  <- 1:nrow(conf$keyLogFsta)
nr      <- nrow(conf$keyLogFsta)
nc      <- ncol(conf$keyLogFsta)

t <- 
  as.data.frame(conf$keyLogFpar) %>% 
  rownames_to_column(var="fleet") %>% 
  pivot_longer(names_to="test", values_to="id", cols=1:nc+1) %>% 
  bind_cols(age=rep(ages, nr)) %>% 
  filter(id != -1) %>% 
  group_by(id, fleet) %>% 
  summarise(ages = paste(age, collapse="_")) %>% 
  ungroup() %>% 
  distinct(id, fleet, ages) %>% 
  mutate(id = paste0("keyLogFpar_",id)) 








    

#Load covariance matrices
load("covCatch.RData") #catch covariance matrices are now stored in the variable "covCatch"

#include covariance
attr(cn, "cov-weight") <- covCatch

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

conf<-loadConf(dat,"model.cfg", patch=TRUE)
par<-defpar(dat,conf)
fit2<-sam.fit(dat,conf,par)

AIC(fit,fit2)
ssbplot(c(fit,fit2))

setwd("D:/GIT/mptools/tcsam/day4")

# =========================================================================
# Using partly total catches; partly catches by age
# =========================================================================

# =========================================================================
# Index observations
# =========================================================================

rm(list=ls())

load("allfleets.RData")
# sum(is.na(allfleets$obs))

compile("allfleets0.cpp")
dyn.load(dynlib("allfleets0"))

allfleets$keyQ <- rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,NA),
                        c(NA, 0, 1, 2, 3, 4, 5, 6,NA),
                        c( 7, 8, 9,10,11,12,NA,NA,NA))

allfleets$keySd <- rbind(c( 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         c(NA, 1, 1, 1, 1, 1, 1, 1,NA),
                         c( 2, 2, 2, 2, 2, 2,NA,NA,NA))

# sum(!is.na(unique(as.numeric(allfleets$keyQ))))

par <- list()
par$logQ <- numeric(max(allfleets$keyQ, na.rm=TRUE)+1)
par$logsd <- numeric(max(allfleets$keySd, na.rm=TRUE)+1)
par$missing <- numeric(sum(is.na(allfleets$obs)))

obj <- MakeADFun(allfleets, par, random="missing", DLL="allfleets0")
fit <- nlminb(obj$par, obj$fn, obj$gr)
est <- obj$report()$logPred

fit$objective

library(tidyverse)
cbind(allfleets$aux, obs=log(allfleets$obs), est) %>% 
  data.frame() %>% 
  ggplot(aes(x=year, y=obs)) +
  theme_bw() +
  geom_point(aes(colour=factor(age))) +
  geom_line(aes(y=est, colour=factor(age))) +
  facet_wrap(~fleet)




# =========================================================================
# Index observations in blocks (preparing for being to include correlations)
# =========================================================================

rm(list=ls())

load("allfleetsblock.RData")

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




# =========================================================================
# 
# =========================================================================

