#Hyperparameters
phi = 0.90
sigma = 0.4
n = 100

#Sample latent effect
x = rep(0,n)
set.seed(123)
x = rnorm(1,0,sd = sqrt(sigma^2/(1-phi^2)))
for(i in 2:n){
  x[i] = phi*x[i-1] +rnorm(1,mean = 0,sd = sigma) 
}

#Sample observations
y = rpois(n,exp(x))
save(y,file = "cpue.RData")

library(TMB)
compile("ar1LatentSolution.cpp")
dyn.load("ar1LatentSolution")

data = list(y = y)
par = list(logSigma = -2,
           phiTrans = 1,
           gamma = rep(0,length(y)))

data$code=0

obj = MakeADFun(data,par,random = "gamma",DLL = "ar1LatentSolution")
opt = nlminb(obj$par,obj$fn,obj$gr,control = list(trace = 1))
rep = sdreport(obj,getJointPrecision = TRUE)

pl = as.list(rep,"Est")
plSd = as.list(rep,"Std")

estPhi =  2/(1 + exp(-2*pl$phiTrans))-1
estPhiU =  2/(1 + exp(-2*(pl$phiTrans + 2*plSd$phiTrans)))-1
estPhiL =  2/(1 + exp(-2*(pl$phiTrans - 2*plSd$phiTrans)))-1
c(estPhi,estPhiL,estPhiU)

sigma =  exp(pl$logSigma)
sigmaU =  exp(pl$logSigma + 2*plSd$logSigma)
sigmaL =  exp(pl$logSigma - 2*plSd$logSigma)
c(sigma,sigmaL,sigmaU)


library(SparseM)
nameIndex = which(colnames(rep$jointPrecision)=="gamma")
Q = rep$jointPrecision[nameIndex,nameIndex]
Q[Q!=0] = 1  # We are only interested if they are 0 or not

#pdf("../presentation/sparsness.pdf")
image(Q, main = "Sparsness structure of AR(1)")
#dev.off()



#Plot toy data
library(maps)
library(maptools)
newmap <- map("world", c("Norway","Denmark","UK","Faroe Islands","Sweden", "Finland","Russia"),fill = TRUE,plot = FALSE, col = "transparent")
IDs <- sapply(strsplit(newmap$names, ":"), function(x) x[1])
map.sp <- map2SpatialPolygons(
  newmap, IDs = IDs,
  proj4string = CRS("+proj=longlat +datum=WGS84"))


lineX = c(6,0)
lineY = c(56,62)
dist = sqrt((lineX[1]-lineX[2])^2 + (lineY[1] - lineY[2])^2)
dd = matrix(0,n,2)
dd[1,] = c(lineX[1],lineY[1])
for(i in 2:n){
  dd[i,] = c(lineX[1],lineY[1]) + i/n *(c(lineX[2],lineY[2])  -c(lineX[1],lineY[1]) )
}

#pdf("../presentation/mapObs.pdf")
plot(lineX-99,lineY-99,xlim = c(-3,12),ylim = c(56,63),xlab = "Lon", ylab = "Lat", main = "Catch per unit effort",cex.main = 2)
points(dd,cex =log(y + 1.1))
plot(map.sp,add = TRUE, col = gray(0.9))
#dev.off()

