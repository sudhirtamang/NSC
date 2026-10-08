library(mvtnorm)
library(glasso)
source("getnorm.R")

n<-300
p<-100
Theta<-matrix(0,p,p)
Theta[1:3,1:3]<-0.5
diag(Theta)<-1


Sigma<-solve(Theta)
Lambda<-diag(Sigma)
Sigma<-diag(Lambda^(-1/2))%*%Sigma%*%diag(Lambda^(-1/2))
Theta<-diag(Lambda^(1/2))%*%Theta%*%diag(Lambda^(1/2))

rho<-(1:100)/100
nruns<-500

TP<-matrix(0,nruns,length(rho))
TPR<-matrix(0,nruns,length(rho))

FP<-matrix(0,nruns,length(rho))
FPR<-matrix(0,nruns,length(rho))

Theta.nd<-Theta
diag(Theta.nd)<-0

nz.total<-sum(c(Theta.nd)!=0)
zero.total<-sum(c(Theta)==0)

for(iruns in 1:nruns){
  set.seed(20151207+iruns)
  x<-rmvnorm(n,rep(0,p),Sigma)
    
  
  x.norm<-matrix(0,n,p)
  
  for(i in 1:p){
    obj.norm<-getnorm(x[,i],type="NS")
    x.norm[,i]<-obj.norm$x.norm
  }
  
  Sigma.hat<-cov(x.norm)
for(i in 1:length(rho)){  
obj<-glasso(Sigma.hat,rho=rho[i],penalize.diagonal=FALSE)
Theta.hat<-obj$wi
Theta.hat.nd<-Theta.hat
diag(Theta.hat.nd)<-0
TPR[iruns,i]<-sum(c(Theta.hat.nd)!=0&c(Theta.nd)!=0)/nz.total
FPR[iruns,i]<-sum(c(Theta.hat.nd)!=0&c(Theta)==0)/zero.total
}
 # print(paste('i =',iruns,Sys.time()))#
}

TPR.m<-apply(TPR,2,mean)
FPR.m<-apply(FPR,2,mean)

TPR.m
FPR.m

plot(FPR.m,TPR.m,type="l") #

