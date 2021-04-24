# install.packages("TMB")
# install.packages("pkgbuild")
# pkgbuild::has_rtools()
# pkgbuild::rtools_path()

# download.file("https://github.com/kaskr/adcomp/archive/master.zip", "adcomp-master.zip")
# unzip("adcomp-master.zip")

# devtools::install_github('fishfollower/SAM/stockassessment')
library(stockassessment)
library(TMB)

setwd("D:/GIT/tcsam/day2")

Sys.setenv(PATH=paste("C:/Rtools/bin", Sys.getenv("PATH"), sep =";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

rm(list=ls())

