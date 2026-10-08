# ---------------------------------------------------------
# Project: NSC
# Description: Used the original code for the paper "Coordinatewise Gaussianization: Theories and Applications" by
#              Qing Mai, Di He & Hui Zou and added the part where the simulation is repeated the desired number of times(RUNs)
# Author: Sudhir Tamang
# Date: September 30, 2026
# ---------------------------------------------------------

rm(list = ls())
library(mvtnorm)
library(pamr)
source("functions.R")


n.test=10000

pi <- 0.5

delta=(2*qnorm(0.10))^2

n=200

p <- 2000

x <- matrix(0,p,nrow=n)

sigma <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    sigma[i,j] <- 0^abs(i-j)
  }
}

mu1 <- c(rep(1,10),rep(0,p-10))
a <- drop(sqrt(delta/(mu1%*%solve(sigma)%*%mu1)))
mu1 <- (a*mu1)

mu2 <- rep(0,p)

RUNs <- 500
error <- numeric(RUNs)
C <- numeric(RUNs)
IC <- numeric(RUNs)

for(h in 1:RUNs){

  x.test <- matrix(0,p,nrow=n.test)

  seed <- 200901*h
  set.seed(seed, kind = NULL, normal.kind = NULL)

  y.test <- rep(0,n.test)
  z <- runif(n.test,0,1)
  id <- which(z>0.5)
  n2.test <- length(id)
  n1.test <- n.test-n2.test
  y.test[id] <- 1
  x.test <- matrix(rnorm(n.test*p),n.test,p)
  x.test[id,] <- sweep(x.test[id,],2,mu1)

  ## nruns <- 500
  ## error <- rep(0,nruns)
  ## C <- rep(0,nruns)
  ## IC <- rep(0,nruns)


  # h=1

  ##=== We demonstrate the following procedure with only one run. For a
  #simulation study with more replicates, readers can wrap the following
  #code with a loop, modify the objects "error, C, IC" into vectors, and
  #report the median of "error, C, IC".

  seed <- 30901+h*7
  set.seed(seed, kind = NULL, normal.kind = NULL)

  y <- rep(0,n)
  z <- runif(n,0,1)
  id <- which(z>0.5)
  n2 <- length(id)
  n1 <- n-n2
  y[id] <- 1
  x <- matrix(rnorm(n*p),n,p)
  x[id,] <- sweep(x[id,],2,mu1)

  x.uncontam <- x[]
  p.contam <- 0.05
  for(i in 1:p){
    idx.contam <- which(rbinom(200, 1, p.contam) == 1)
    x[idx.contam, i] <- rnorm(length(idx.contam), 0, sqrt(3))
  }

  x.norm <- matrix(0,n,p)
  x.test.norm <- matrix(0,n.test,p)
  for(i in 1:p){
    obj.norm <- InvNorm(x[,i],c=0)
    x.norm[,i] <- obj.norm$x.norm
    x.test.norm[,i] <- obj.norm$transform(x.test[,i])
  }

  dat.train <- list(x=t(x.norm),y=y,geneid=c(1:p))
  dat.test <- list(x=t(x.test.norm),y=y.test,geneid=c(1:p))
  junk <- capture.output(obj.pam <- pamr.train(dat.train))
  set.seed(1256*h)
  junk <- capture.output(obj.pamr.cv <- pamr.cv(obj.pam,dat.train))
  t <- obj.pamr.cv$threshold[which.min(obj.pamr.cv$error)]
  fit <- pamr.predict(obj.pam,t(x.test.norm),type="class",threshold=t)
  fit <- as.numeric(fit)-1
  error[[h]] <- mean(abs(fit-y.test))#

  pam.sel <- rep(0,p)
  junk <- capture.output(pam.sel[as.numeric(pamr.listgenes(obj.pam,dat.train,threshold=t)[,1])] <- 1)

  Iv.true <- rep(0,p)
  Iv.true[mu1!=0] <- 1

  Iv.hat <- rep(0,p)
  Iv.hat[pam.sel>0] <- 1

  C[[h]] <- sum((Iv.true+Iv.hat)==2)#
  IC[[h]] <- sum((Iv.true-Iv.hat)== -1)#
}



#== results ===
cat("Error Median: ", median(error), "\n")
cat("FP Median: ", median(10-C), "\n")
cat("FN Median: ", median(IC), "\n")
cat("Standard Error Results!!!\n\n")
cat("Error Median: ", sd(error*100)/sqrt(RUNs), "\n")
cat("FP Median: ", sd(10-C)/sqrt(RUNs), "\n")
cat("FN Median: ", sd(IC)/sqrt(RUNs), "\n")




