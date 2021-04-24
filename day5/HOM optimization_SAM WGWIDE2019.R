# --------------------------------------------------------------------------------------
# SAM assessment for Western Horse Mackerel
#
# Original: Vanessa Trijoulet
# 02/09/2019 Adapted by Martin Pastoors for use in IBPWHOM
# 
# Session info: 
# R version 3.5.0 (2018-04-23)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] forcats_0.3.0         stringr_1.3.1         dplyr_0.8.0.1         purrr_0.2.5           readr_1.1.1          
# [6] tidyr_0.8.2           tibble_2.1.1          ggplot2_3.1.0         tidyverse_1.2.1       stockassessment_0.8.1
# [11] TMB_1.7.15           
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.2       cellranger_1.1.0 pillar_1.3.1     compiler_3.5.0   plyr_1.8.4       tools_3.5.0      lubridate_1.7.4 
# [8] jsonlite_1.6     nlme_3.1-137     gtable_0.2.0     lattice_0.20-35  pkgconfig_2.0.2  rlang_0.3.4      Matrix_1.2-15   
# [15] cli_1.0.1        rstudioapi_0.8   yaml_2.2.0       parallel_3.5.0   haven_1.1.2      withr_2.1.2      xml2_1.2.0      
# [22] httr_1.3.1       generics_0.0.2   hms_0.4.2        grid_3.5.0       tidyselect_0.2.5 glue_1.3.0       ellipse_0.4.1   
# [29] R6_2.3.0         readxl_1.1.0     modelr_0.1.2     magrittr_1.5     backports_1.1.2  scales_1.0.0     rvest_0.3.2     
# [36] assertthat_0.2.0 colorspace_1.3-2 stringi_1.2.4    lazyeval_0.2.1   munsell_0.5.0    broom_0.5.2      crayon_1.3.4  
# --------------------------------------------------------------------------------------

library(TMB)
library(stockassessment)   # devtools::install_github("fishfollower/SAM/stockassessment", subdir = "SAM")
library(tidyverse)


# get HOM data on stockassessment.org  
fit=fitfromweb("WHOM2")
# load("D:/TEMP/WHOM3/run/model.RData")


#### Default config ####
def.conf <- defcon(fit$data)
par      <- defpar(fit$data,def.conf)

fit.def  <- sam.fit(fit$data,def.conf,par) 
res      <- residuals(fit.def)
retro    <- retro(fit.def, year=5)
AIC(fit.def)
plot(fit.def)
plot(res)
plot(retro)
# same that fit
# big residuals for ages 0-2 and opposite residuals for ages 14-15

def.conf <- defcon(fit$data)
par      <- defpar(fit$data,def.conf)

#### Change obs variance catch first age groups ####
conf2                   <- def.conf
conf2$keyVarObs[1,2:3]  <- c(1,2)
conf2$keyVarObs[1,4:16] <- 3
conf2$keyVarObs[2:4,1]  <- 4:6
par                     <- defpar(fit$data,conf2)

fit2                    <- sam.fit(fit$data,conf2,par) 
res2                    <- residuals(fit2)
AIC(fit2)
# year effect in residuals catch at beginning time series

#### Change obs variance catch first and last age groups ####
conf2b                   <- def.conf
conf2b$keyVarObs[1,2:3]  <- c(1,2)
conf2b$keyVarObs[1,4:15] <- 3
conf2b$keyVarObs[1,16]   <- 4
conf2b$keyVarObs[2:4,1]  <- 5:7
par                      <- defpar(fit$data,conf2b)

fit2b                    <- sam.fit(fit$data,conf2b,par) 
res2b                    <- residuals(fit2b)
retro2b                  <- retro(fit2b)
AIC(fit2b)

# #### Remove 1st year of data ####
# fit3 <- runwithout(fit2b, year=1982)
# res3 <- residuals(fit3)
# 

#### Test with variance in F changing at age ####
conf3                 <- conf2b
conf3$keyVarF[1,2:16] <- c(1:14,14)
par                   <- defpar(fit$data,conf3)

fit3                  <- sam.fit(fit$data,conf3,par) 
res3                  <- residuals(fit3)
AIC(fit3)
fit3$sdrep

#### Test with grouped variance in F for certain ages given values in fit3$sdrep ####
conf4                  <- conf2b
conf4$keyVarF[1,5:6]   <- 1:2
conf4$keyVarF[1,7:9]   <- 3
conf4$keyVarF[1,10:12] <- 4
conf4$keyVarF[1,13:16] <- 5
par                    <- defpar(fit$data,conf4)

fit4                   <- sam.fit(fit$data,conf4,par) 
res4                   <- residuals(fit4)
retro4                 <- retro(fit4, year=5)
AIC(fit4)
plot(res4)
empirobscorrplot(res4)


#### Test with correlation in catch ages (AR) ####
conf5                   <- conf4
conf5$obsCorStruct[1]   <- "AR"
conf5$keyCorObs[1,1:6]  <- 0
conf5$keyCorObs[1,7:15] <- 1
par                     <- defpar(fit$data,conf5)

fit5                    <- sam.fit(fit$data,conf5,par) 
res5                    <- residuals(fit5)
retro5                  <- retro(fit5, year=5)
AIC(fit5)
plot(res5)
plot(retro5)
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
par                     <- defpar(fit$data,conf6)
rho                     <- 0.999999
par$itrans_rho          <- -log(2/(1+rho)-1)/2

fit6                    <- sam.fit(fit$data,conf6,par, map = list(itrans_rho = as.factor(NA))) 
res6                    <- residuals(fit6)
retro6                  <- retro(fit6, year=5)
AIC(fit6)
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


