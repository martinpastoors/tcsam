# install.packages("TMB")
# install.packages("pkgbuild")
# pkgbuild::has_rtools()
# pkgbuild::rtools_path()

# download.file("https://github.com/kaskr/adcomp/archive/master.zip", "adcomp-master.zip")
# unzip("adcomp-master.zip")

# devtools::install_github('fishfollower/SAM/stockassessment')
library(stockassessment)
library(TMB)

setwd("D:/GIT/tcsam/day5")
Sys.setenv(PATH=paste("C:/Rtools/bin", Sys.getenv("PATH"), sep =";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# One step ahead residuals

rm(list=ls())

load("cpue.RData")
compile("ar1.cpp")
dyn.load(dynlib("ar1"))

data = list(y = y)
par = list(logSigma = -2, phiTrans = 1, gamma = rep(0,length(y)))

obj <- MakeADFun(data,par,random = "gamma",DLL = "ar1")
fit <- nlminb(obj$par,obj$fn,obj$gr)
res <- oneStepPredict(obj, observation.name = "y", 
                      data.term.indicator = "keep",
                      discrete = TRUE,
                      discreteSupport = 0:20,
                      method = "oneStepGeneric")

hist(res$residual)
qqnorm(res$residual)
abline(0,1)

# wrong/traditional residuals
sdr <- sdreport(obj)
pl <- as.list(sdr, "Est")
sdpl <- as.list(sdr, "Std")

hist((data$y-pl$gamma)/sdpl$gamma)

# example to calculate 
estX <- summary (sdr ,"random")
C <- solve(obj$env$spHess(obj$env$last.par.best, random = TRUE )) # solve is inverting a matrix; in standard R
gamma.star <- MASS::mvrnorm (1, estX [,1] ,C)

rep <- obj$report()

# remove the first point and the last point
res <- (gamma.star[-1]-rep$phi*gamma.star[-length(gamma.star)])/rep$sd

qqnorm(res); abline(0,1)

fit.jit <- jit(fit, nojit=10)


# Explore example diagnostics

fit <- (fitfromweb("WHOM_2020"))
fit

res  <- residuals(fit) # one observation ahead residuals
plot(res)

res2 <- procres(fit)   # process residuals (single joint sample)
plot(res2)

ret  <- retro(fit, year=5, ncores=2) # retrospective
ssbplot(ret)

lo  <- leaveout(fit) # leaveout
ssbplot(lo)

jt  <- jit(fit, nojit=10) # leaveout
jt  <- jit(fit, nojit=10, sd=0.1) # leaveout
ssbplot(jt)
jt

sim  <- simstudy(fit, nsim=10) # simulation of assessment

ssbplot(sim)

# ===============================================================================
# cross validation
# ===============================================================================

# function for cross - validation
xval <- function(fit, year=NULL, fleet=NULL, age=NULL, ...){
  data <- fit$data
  nam <- c("year", "fleet", "age")[c(length(year)>0,length(fleet)>0,length(age)>0)]
  if((length(year)==0) & (length(fleet)==0) & (length(age)==0)){
    idx <- rep(TRUE,nrow(data$aux))
  }else{
    idx <- !do.call(paste, as.data.frame(data$aux[,nam,drop=FALSE])) %in% do.call(paste, as.data.frame(cbind(year=year, fleet=fleet, age=age)))
  }
  idx <- !idx
  data$logobs[idx] <- NA
  idx2 <- which(is.na(data$logobs))
  conf <- fit$conf
  par <- defpar(data,conf)
  thisfit <- sam.fit(data, conf, par, rm.unidentified = TRUE, newtonsteps=0, silent=TRUE,...)
  ret <- as.data.frame(cbind(data$aux[idx2,], obs=fit$data$logobs[idx2], pred=thisfit$pl$missing, predSd=thisfit$plsd$missing))
  ret <- ret[complete.cases(ret),]
  attr(ret, "fit") <- thisfit
  return(ret)
}

#predicting the observations in the last 5 years
pred <- xval (fit , year =c (2015:2020) )
plot(pred$pred, pred$obs)
abline(0,1)

sum((pred$obs - pred$pred)^2)

p <- as.data.frame(pred)
p %>% 
  ggplot(aes(x=obs, y=pred)) +
  theme_bw() +
  # theme(legend.position = "none") +
  geom_point(aes(colour=(year))) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fleet)

p %>% 
  filter(fleet==1) %>% 
  ggplot(aes(x=obs, y=pred)) +
  theme_bw() +
  # theme(legend.position = "none") +
  geom_point(aes(colour=(year))) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~age)


# cod prediction of observation
fit <- sam.fit(nscodData, nscodConf, nscodParameters)
pred <- xval (fit , year =c (2011:2025) )
p <- as.data.frame(pred)
p %>% 
  ggplot(aes(x=obs, y=pred)) +
  theme_bw() +
  # theme(legend.position = "none") +
  geom_point(aes(colour=(year))) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~fleet)

p %>% 
  filter(fleet==1) %>% 
  ggplot(aes(x=obs, y=pred)) +
  theme_bw() +
  # theme(legend.position = "none") +
  geom_point(aes(colour=(year))) +
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~age)


# ===============================================================================
# Forecasting methods
# ===============================================================================

compile("rwForecast.cpp")
dyn.load(dynlib("rwForecast"))
load("F.RData")

parameters <- list(
  logSdRw=0,
  logSdObs=0,
  lam=rep(0,length(data$y))
)

obj  <- MakeADFun(data,parameters,random="lam",DLL="rwForecast")
opt  <-nlminb(obj$par,obj$fn,obj$gr)
sdr  <-sdreport(obj)
pl   <- as.list(sdr,"Est")   # parameters
plsd <- as.list(sdr,"Std")   # parameter uncertainties


plot(data$y, xlim=c(0,50))
lines(pl$lam)
lines(pl$lam-2*plsd$lam, lty="dotted")
lines(pl$lam+2*plsd$lam, lty="dotted")

nsim <- 1000
nYears = 6 #Forecast years added base year
sim <- matrix(NA, ncol=nsim, nrow=nYears)
set.seed(123)
sim[1,] <- rnorm(nsim, pl$lam[length(data$y)], sd=plsd$lam[length(data$y)]) #Simulate in base year
for(y in 2:nYears){
  sim[y,] <- rnorm(n=nsim, mean=sim[y-1,], sd=exp(pl$logSdRw) )#Simulate the process in future years
  # sim[y,] <- rnorm(n=nsim, mean=data$futureObs[y-1], sd=exp(pl$logSdRw) )#Simulate the process in future years
}

#Plot future realisations
plot(data$y, xlim=c(0,50),main = "Forecast",ylab = "F",cex.main = 2,
     cex.lab = 1.6)
lines(pl$lam)
lines(pl$lam-2*plsd$lam, lty="dotted")
lines(pl$lam+2*plsd$lam, lty="dotted")
matplot((length(data$y)):(length(data$y) + nYears-1), sim, type="l", add=TRUE, col=gray(.8))
matplot((length(data$y)):(length(data$y) + nYears-1), t(apply(sim, 1, quantile, c(0.025, 0.5, 0.975))), type="l", col="red", add=TRUE, lwd=3,lty = c(3,1,3))
points((length(data$y) +1):(length(data$y) + nYears-1), data$futureObs)


# ===============================================================================
# Forecasting of North Sea code
# ===============================================================================

setwd("D:/GIT/mptools/tcsam/day5/practicalForecast")

#Fit SAM with currently (2019) used configurations
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
conf<-loadConf(dat,"model.cfg")
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

#Forecast
forecastSAM = forecast(fit, 
                       catchval=c(800000,NA,NA), 
                       fval=c(NA, 0.503,NA),
                       fscale=c(NA,NA,1))

#Forecast exercise
set.seed(123)
forecastSAM = forecast(fit, 
                       catchval=c(697000,NA,NA,NA), 
                       fval=c(NA, 0.503,NA,NA),
                       fscale=c(NA,NA,1,1))

t <-
  as.data.frame(attr(forecastSAM, "tab")) %>% 
  rownames_to_column(var="year") %>% 
  pivot_longer(names_to="variable", values_to="value", c("fbar:median":"catch:high")) %>% 
  separate(variable, into=c("var","metric"), sep=":")   

t %>% 
  filter(year > 2019) %>% 
  group_by(var, metric) %>% 
  summarise(value = mean(value, na.rm=TRUE))

# calculate from the tabl object
mean(attributes(forecastSAM)$tab[c("2022","2020","2021"),"catch:median"])

ffc <- 
  function(x) forecast(fit, 
                       fscale= c(1,NA,NA,NA,NA), 
                       fval  = c(NA,x,x,x,x), 
                       processNoiseF = FALSE)

ffc(0.00000001)

x <- 0.0
df <- data.frame(stringsAsFactors = FALSE)
for (x in c(10^-6, 0.2, 0.4, 0.503, 0.6, 0.8)) {
   df <- 
     as.data.frame(attr(ffc(x), "tab")) %>% 
     rownames_to_column(var="year") %>% 
     pivot_longer(names_to="variable", values_to="value", c("fbar:median":"catch:high")) %>% 
     separate(variable, into=c("var","metric"), sep=":")  %>% 
     mutate(f = x) %>%
     mutate(year = as.integer(year)) %>% 
     pivot_wider(names_from="metric", values_from="value") %>% 
     bind_rows(df, .)
}

df %>% 
  ggplot(aes(x=year, y=median, group=factor(f))) +
  theme_bw() +
  geom_line(aes(colour=factor(f))) +
  geom_ribbon(aes(ymin=low, ymax=high, fill=factor(f)), alpha=0.5) +
  facet_wrap(~var, scales="free_y")

df %>% 
  filter(var=="catch") %>% 
  
  ggplot(aes(x=year, y=median, group=factor(f))) +
  theme_bw() +
  geom_line(aes(colour=factor(f))) +
  geom_ribbon(aes(ymin=low, ymax=high, fill=factor(f)), alpha=0.5) +
  labs(x="", y="catch") +
  facet_wrap(~f)


setwd("D:/GIT/mptools/tcsam/day5")









# Catch observations

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
names(surveys)

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


library(stockassessment)
library(tidyverse)
fitfromwebtemp<-function (stockname, character.only = FALSE) 
{
  if (!character.only) 
    stockname <- as.character(substitute(stockname))
  fit <- NA
  load(url(sub("SN", stockname, "https://www.stockassessment.org/datadisk/stockassessment/userdirs/user3/SN/run/model.RData")))
  fit
}

fit <- fitfromwebtemp("BW_2020")
set.seed(123)
fc  <- forecast(fit, 
                catchval=c(1400000,NA,NA,NA,NA,NA,NA), 
                fval=c(NA, 0.32,NA,NA,NA,NA,NA),
                fscale=c(NA,NA,1,1,1,1,1))

t <-
  as.data.frame(attr(fc, "tab")) %>% 
  rownames_to_column(var="year") %>% 
  pivot_longer(names_to="variable", values_to="value", c("fbar:median":"catch:high")) %>% 
  separate(variable, into=c("var","metric"), sep=":")   

t %>% 
  filter(year > 2019) %>% 
  group_by(var, metric) %>% 
  summarise(value = mean(value, na.rm=TRUE))

myvar <- "catch"

t %>% 
  filter(var== myvar) %>% 
  pivot_wider(names_from = "metric", values_from="value") %>% 
  mutate(year=as.integer(year)) %>% 
  ungroup() %>% 
  
  ggplot(aes(x=year, y=median)) +
  theme_bw() +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.5) +
  labs(x="", y=myvar) +
  expand_limits(y=0)
