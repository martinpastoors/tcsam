# --------------------------------------------------------------------------------------
# SAM assessment for Western Horse Mackerel
#
# Original: Vanessa Trijoulet
# 02/09/2019 Adapted by Martin Pastoors for use in IBPWHOM
# 19/10/2020 Adapted after TCSAM 2020 by Martin Pastoors
#
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] TMB_1.7.18            stockassessment_0.9.0
# 
# loaded via a namespace (and not attached):
# [1] lattice_0.20-41  crayon_1.3.4     dplyr_1.0.2      MASS_7.3-51.6    grid_4.0.2       R6_2.4.1        
# [7] lifecycle_0.2.0  magrittr_1.5     ellipse_0.4.2    pillar_1.4.6     rlang_0.4.7      rstudioapi_0.11 
# [13] Matrix_1.2-18    ellipsis_0.3.1   vctrs_0.3.2      generics_0.0.2   tools_4.0.2      glue_1.4.1      
# [19] purrr_0.3.4      tinytex_0.25     xfun_0.16        parallel_4.0.2   compiler_4.0.2   pkgconfig_2.0.3 
# [25] tidyselect_1.1.0 tibble_3.0.3    
# --------------------------------------------------------------------------------------


library(TMB)
library(stockassessment)   # devtools::install_github("fishfollower/SAM/stockassessment", subdir = "SAM")
library(tidyverse)


# get HOM data on stockassessment.org  
fit <- (fitfromweb("WHOM_2020"))
fit


rm(list=ls())

setwd("D:/TEMP/WHOM_2020")

cn<-read.ices("data/cn.dat")
cw<-read.ices("data/cw.dat")
dw<-read.ices("data/dw.dat")
lf<-read.ices("data/lf.dat")
lw<-read.ices("data/lw.dat")
mo<-read.ices("data/mo.dat")
nm<-read.ices("data/nm.dat")
pf<-read.ices("data/pf.dat")
pm<-read.ices("data/pm.dat")
sw<-read.ices("data/sw.dat")
surveys<-read.ices("data/survey.dat")

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

conf<-loadConf(dat,"conf/model.cfg", patch=TRUE)
par<-defpar(dat,conf)

fit<-sam.fit(dat,conf,par)
res<- residuals(fit)
procres <- procres(fit)
aic<- AIC(fit)
retro<- retro(fit, year=5)
loo <-leaveout(fit)
runs <- simstudy(fit,nsim=1000)

fit.df <- 
  as.data.frame(summary(fit)) %>% 
  rownames_to_column() %>% 
  setNames(c("year","recr","lowrecr","highrecr","ssb","lowssb","highssb","f","lowf","highf")) %>% 
  mutate(year = as.integer(year)) %>% 
  gather(key="variable", value="data", recr:highf) %>% 
  mutate(est = ifelse(grepl("low", variable), "low", NA), 
         est = ifelse(grepl("high", variable), "high", est),
         est = ifelse(is.na(est), "est", est),
         variable = gsub("low|high","", variable),
         assess="sam") %>% 
  spread(key=est, value=data)

save(fit, file="run/model.RData")
save(res, file="run/residuals.RData")
save(procres, file="run/procresiduals.RData")
save(retro, file="run/retro.RData")
save(loo, file="run/leaveout.RData")
save(fit.df, file="run/fit.df.RData")

plot(fit)
plot(retro)
plot(loo)


#### Default config ####
conf0    <- defcon(fit$data)
par0     <- defpar(fit$data, conf0)
fit0     <- sam.fit(fit$data,conf0,par0) 
res0     <- residuals(fit0)
aic0     <- AIC(fit0)

# retro0   <- retro(fit.def, year=5)
# plot(fit.def)
plot(res)
# plot(retro)
# big residuals for ages 0-2 and opposite residuals for ages 14-15


#### Change obs variance catch first three age groups ####
conf2                   <- conf0
conf2$keyVarObs[1,2:3]  <- c(1,2)   # set separate variances for catch age  0, 1 and 2
conf2$keyVarObs[1,4:16] <- 3        # same variances from age 3 onwards
conf2$keyVarObs[2:4,1]  <- 4:6      # set 1 variance for each of the surveys
par2                    <- defpar(fit$data,conf2)
fit2                    <- sam.fit(fit$data,conf2,par2) 
res2                    <- residuals(fit2)
aic2                    <- AIC(fit2)
plot(res2)
# year effect in residuals catch at beginning time series

#### Change obs variance catch first three age groups and the last age group ####
conf2b                   <- conf0
conf2b$keyVarObs[1,2:3]  <- c(1,2)
conf2b$keyVarObs[1,4:15] <- 3
conf2b$keyVarObs[1,16]   <- 4
conf2b$keyVarObs[2:4,1]  <- 5:7
par2b                      <- defpar(fit$data,conf2b)

fit2b                    <- sam.fit(fit$data,conf2b,par2b) 
#res2b                    <- residuals(fit2b)
aic2b                    <- AIC(fit2b)

# #### Remove 1st year of data ####
fit2c <- runwithout(fit2b, year=1982)
# res2c <- residuals(fit2c)
# modeltable(c(fit0,fit2, fit2b, fit2c))

#### Test with variance in F changing at age ####
conf3                 <- conf2b
conf3$keyVarF[1,2:16] <- c(1:14,14)
par3                  <- defpar(fit$data,conf3)
fit3                  <- sam.fit(fit$data,conf3,par3) 
aic3                  <- AIC(fit3)
# res3                  <- residuals(fit3)
# fit3$sdrep

#### Test with grouped variance in F for certain ages given values in fit3$sdrep ####
conf4                  <- conf2b
conf4$keyVarF[1,5:6]   <- 1:2
conf4$keyVarF[1,7:9]   <- 3
conf4$keyVarF[1,10:12] <- 4
conf4$keyVarF[1,13:16] <- 5
par4                   <- defpar(fit$data,conf4)

fit4                   <- sam.fit(fit$data,conf4,par4) 
aic4                   <- AIC(fit4)
# res4                   <- residuals(fit4)
# retro4                 <- retro(fit4, year=5)
# plot(res4)
# empirobscorrplot(res4)


#### Test with correlation in catch ages (AR) ####
conf5                   <- conf4
conf5$obsCorStruct[1]   <- "AR"
conf5$keyCorObs[1,1:6]  <- 0
conf5$keyCorObs[1,7:15] <- 1
par5                    <- defpar(fit$data,conf5)

fit5                    <- sam.fit(fit$data,conf5,par5) 
aic5                    <- AIC(fit5)
res5                    <- residuals(fit5)
plot(res5)
procres5                <- procres(fit5)
plot(procres5)
# retro5                  <- retro(fit5, year=5)
# plot(retro5)
matplot(t(faytable(fit5)), type="l")

# #### Decouple last ages for F ####
# conf6 <- conf5
# conf6$keyLogFsta[1,16] <- 15
# par <- defpar(fit$data,conf6)
# fit6 <- sam.fit(fit$data,conf6,par) 
# AIC(fit6)
# res6 <- residuals(fit6)
# plot(res6)
# retro6 <- retro(fit6, year=5)
# plot(retro6)


#### Constant F selectivity to reduce retro ####
conf6                   <- conf5
conf6$keyCorObs[1,7:15] <- 0
par6                    <- defpar(fit$data,conf6)
rho                     <- 0.999999
par6$itrans_rho         <- -log(2/(1+rho)-1)/2

fit6                    <- sam.fit(fit$data,conf6,par6, 
                                   map = list(itrans_rho = as.factor(NA))) 
aic6                    <- AIC(fit6)
res6                    <- residuals(fit6)
retro6                  <- retro(fit6, year=5)
plot(res6)
plot(retro6)
fit6$sdrep


#### Constant F selectivity but change variance in F process (looking at sdrep of fit6, similar values for logSdLogFsta) ####
#### Best fit with compromise AIC and retro
def.conf                <- defcon(fit$data)

conf2                   <- def.conf
conf2$keyVarObs[1,2:3]  <- c(1,2)
conf2$keyVarObs[1,4:16] <- 3
conf2$keyVarObs[2:4,1]  <- 4:6

conf2b                   <- def.conf
conf2b$keyVarObs[1,2:3]  <- c(1,2)
conf2b$keyVarObs[1,4:15] <- 3
conf2b$keyVarObs[1,16]   <- 4
conf2b$keyVarObs[2:4,1]  <- 5:7

conf4                  <- conf2b
conf4$keyVarF[1,5:6]   <- 1:2
conf4$keyVarF[1,7:9]   <- 3
conf4$keyVarF[1,10:12] <- 4
conf4$keyVarF[1,13:16] <- 5

conf5                   <- conf4
conf5$obsCorStruct[1]   <- "AR"
conf5$keyCorObs[1,1:6]  <- 0
conf5$keyCorObs[1,7:15] <- 1

conf6                   <- conf5
conf6$keyCorObs[1,7:15] <- 0

conf7                  <- conf6
conf7$keyVarF[1,7:16]  <- 2
par                    <- defpar(fit$data,conf7)
rho                     <- 0.999999
par$itrans_rho          <- -log(2/(1+rho)-1)/2

fit7                    <- sam.fit(fit$data,conf7,par, map = list(itrans_rho = as.factor(NA))) 
res7                    <- residuals(fit7)
procres7                <- procres(fit7)
retro7                  <- retro(fit7, year=5)
fit7.df                 <- as.data.frame(summary(fit7)) %>% 
  rownames_to_column() %>% 
  setNames(c("year","recr","lowrecr","highrecr","ssb","lowssb","highssb","f","lowf","highf")) %>% 
  mutate(year = as.integer(year)) %>% 
  gather(key="variable", value="data", recr:highf) %>% 
  mutate(est = ifelse(grepl("low", variable), "low", NA), 
         est = ifelse(grepl("high", variable), "high", est),
         est = ifelse(is.na(est), "est", est),
         variable = gsub("low|high","", variable),
         assess="sam") %>% 
  spread(key=est, value=data)






AIC(fit7)
fit7$sdrep
plot(fit7)
fitplot(fit7, fleets=1) 
fitplot(fit7, fleets=2:4) 
ssbplot(fit7)
recplot(fit7)
matplot(t(faytable(fit7)), type="l")

# residuals
plot(res7)
plot(res7, type="summary")

# retro analysis
plot(retro7)

# progressive residuals
plot(procres7)

# include the SS assessment of 2018
ss <- 
  get(load("D:/Dropbox/iAdvice/RData/iAssess.RData")) %>% 
  filter(stockkeylabelold == "hom-west", assessmentyear == 2018) %>% 
  dplyr::select(year, 
                recr=recruitment, highrecr=highrecruitment, lowrecr=lowrecruitment,
                ssb = stocksize, highssb=highstocksize, lowssb=lowstocksize,
                f   = fishingpressure, highf = highfishingpressure, lowf=lowfishingpressure) %>% 
  gather(key="variable", value="data", recr:lowf) %>% 
  mutate(est = ifelse(grepl("low", variable), "low", NA), 
         est = ifelse(grepl("high", variable), "high", est),
         est = ifelse(is.na(est), "est", est),
         variable = gsub("low|high","", variable),
         assess="ss") %>% 
  spread(key=est, value=data)

fit.df <- 
  as.data.frame(summary(fitfromweb("WHOM_2019"))) %>% 
  rownames_to_column() %>% 
  setNames(c("year","recr","lowrecr","highrecr","ssb","lowssb","highssb","f","lowf","highf")) %>% 
  mutate(year = as.integer(year)) %>% 
  gather(key="variable", value="data", recr:highf) %>% 
  mutate(est = ifelse(grepl("low", variable), "low", NA), 
         est = ifelse(grepl("high", variable), "high", est),
         est = ifelse(is.na(est), "est", est),
         variable = gsub("low|high","", variable),
         assess="sam") %>% 
  spread(key=est, value=data)


# Compare SS assessment with SAM assessment
bind_rows(fit.df, ss) %>% 
  filter(year >= 1985) %>% 
  ggplot(aes(x=year,y=est)) +
  theme_bw() +
  geom_line(aes(colour=assess)) +
  geom_ribbon(aes(x=year, ymin=low, ymax=high, fill=assess), alpha=0.5) +
  expand_limits(y=0) +
  facet_wrap(~variable, scales = "free_y")

# Scaled to the mean
bind_rows(fit.df, ss) %>% 
  filter(year >= 1985) %>% 
  group_by(variable, assess) %>% 
  mutate(z = (est - mean(est, na.rm=TRUE))/sd(est, na.rm=TRUE)) %>% 
  ggplot(aes(x=year,y=z)) +
  theme_bw() +
  geom_line(aes(colour=assess)) +
  # geom_ribbon(aes(x=year, ymin=low, ymax=high, fill=assess), alpha=0.5) +
  expand_limits(y=0) +
  facet_wrap(~variable, scales = "free_y")

# table with averages
bind_rows(fit.df, ss) %>% 
  filter(year >= 1985) %>% 
  group_by(variable, assess) %>% 
  summarize(mean = mean(est, na.rm=TRUE),
            sd   = sd(est, na.rm=TRUE)) 
  


# #### Shorter time series ####
fit7b   <- runwithout(fit7, year=1982:1997)
res7b   <- residuals(fit7b)
retro7b <- retro(fit7b, year=5)
plot(fit7b)
catchplot(fit7b)
matplot(t(faytable(fit7b)), type="l")
plot(res7b)
plot(retro7b)
# very bad retrospective




##### Exploration of 2015 - 2016 assessments
fit2015only        <- runwithout(fit7, year=2016:2018)
fit2015withsurvey  <- sam.fit(fitfromweb("WHOMb_VT")$data, conf7,par, map = list(itrans_rho = as.factor(NA))) 
fit2016            <- runwithout(fit7, year=2017:2018)

ssbplot(fit2015only)
ssbplot(fit2015withsurvey)
ssbplot(fit2016)

fitplot(fit2015only, fleets=c(1)) 
fitplot(fit2015only, fleets=c(2,3)) 
fitplot(fit2015withsurvey, fleets=c(1)) 
fitplot(fit2015withsurvey, fleets=c(2,3)) 
fitplot(fit2016, fleets=c(2,3)) 

t1 <- as.data.frame(ssbtable(fit2015only)) %>% rownames_to_column() %>%  mutate(assess="2015only")
t2 <- as.data.frame(ssbtable(fit2015withsurvey)) %>% rownames_to_column() %>%  mutate(assess="2015withsurvey")
t3 <- as.data.frame(ssbtable(fit2016)) %>% rownames_to_column() %>%  mutate(assess="2016")

bind_rows(t1,t2,t3) %>% 
  mutate(year = as.integer(rowname)) %>%
  filter(year > 2000, year < 2017) %>% 
  ggplot(aes(x=year, y=Estimate)) +
  geom_line(aes(colour=assess)) +
  # geom_ribbon(aes(ymin=Low, ymax=High, fill=assess), alpha=0.5) +
  geom_point(aes(colour=assess)) +
  expand_limits(y=0) +
  facet_wrap(~assess)

fitfromweb("WHOM3b")$data

stockassessment::modeltable(fit2015only)

save(fit2015only, fit2015withsurvey, fit2016, file="samfits.RData")
stockassessment::data.plot(fit2015only)
dataplot(fit2015withsurvey)


