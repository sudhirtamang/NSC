# On Oct 06, 2026 added a function do.simulation to do as many RUNs of simulation with given percentage of contamination(p.contam)
# Number of obserbations on test set(n.test), sample size(n), number of predictors(p)

InvNorm<-function(x,c=NULL){
  f<-ecdf(x)
  n<-length(x)
  
  transform<-function(t){
    t.norm<-(f(t)-c/n)*n/(n-2*c+1)
    t.norm[t.norm<1e-6]<-(1-c)/(n-2*c+1)
    t.norm<-qnorm(t.norm)
    t.norm
  }
  x.norm<-transform(x)
  list(x.norm=x.norm,transform=transform)
}

getnorm.w<-function(x,W.parameter=NULL){
  f<-ecdf(x)
  n<-length(x)
  
    if(missing(W.parameter))W.parameter=1/n
    transform<-function(t){
      t.norm<-f(t)
      t.norm[t.norm<W.parameter]<-W.parameter
      t.norm[t.norm>1-W.parameter]<-1-W.parameter
      t.norm<-qnorm(t.norm)
      t.norm
    }

  x.norm<-transform(x)
  list(x.norm=x.norm,transform=transform)
}

ModelN4 <- function(n, p, d = 10, n.test, rand.seed) {
  set.seed(rand.seed)
  
  # Ensures length of MU2 is exactly p
  MU2 <- 0.5 * c(rep(1, d), numeric(p - d))
  
  # --- Training Data ---
  y <- numeric(n)
  z <- runif(n, 0, 1)
  idx <- which(z > 0.5)
  y[idx] <- 1
  
  x <- matrix(rweibull(n * p, shape = 1, scale = 1), nrow = n, ncol = p)
  x[idx, ] <- sweep(x[idx, ], 2, MU2, FUN = "+")
  
  # --- Test Data ---
  y.test <- numeric(n.test)
  z.test <- runif(n.test, 0, 1)
  id <- which(z.test > 0.5)
  y.test[id] <- 1  
  
  x.test <- matrix(rweibull(n.test * p, shape = 1, scale = 1), nrow = n.test, ncol = p)
  x.test[id, ] <- sweep(x.test[id, ], 2, MU2, FUN = "+")
  
  list(x = x, y = y, x.test = x.test, y.test = y.test, MU2 = MU2)
}

ModelN1 <- function(n, p, d = 10, n.test, rand.seed) {
  set.seed(rand.seed)
  
  delta <- (2 * qnorm(0.10))^2
  sigma <- diag(p) 
  
  MU2 <- c(rep(1, d), rep(0, p - d))
  # Scaling factor calculation
  # For sigma = I, MU2 %*% solve(sigma) %*% MU2 simplifies to sum(MU2^2) = d
  quad_form <- as.numeric(t(MU2) %*% solve(sigma, MU2))
  a <- sqrt(delta / quad_form)
  MU2 <- a * MU2
  
  # --- Training Data ---
  y <- numeric(n)
  z <- runif(n, 0, 1)
  idx <- which(z > 0.5)
  y[idx] <- 1
  
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  x[idx, ] <- sweep(x[idx, ], 2, MU2, FUN = "+")
  
  # --- Test Data ---
  y.test <- numeric(n.test)
  z.test <- runif(n.test, 0, 1)
  id <- which(z.test > 0.5)
  y.test[id] <- 1  
  
  x.test <- matrix(rnorm(n.test * p), nrow = n.test, ncol = p)
  x.test[id, ] <- sweep(x.test[id, ], 2, MU2, FUN = "+")
  
  list(x = x, y = y, x.test = x.test, y.test = y.test, MU2 = MU2)
}

ModelN5 <- function(n, p, d = 10, n.test, rand.seed) {
  set.seed(rand.seed)
  
  # Exact constant 'a' for 10% Bayes error under independent Exp(1) / Weibull(1,1)
  a <- log(5) / d
  MU2 <- c(rep(a, d), rep(0, p - d))
  
  # --- Training Data ---
  y <- numeric(n)
  z <- runif(n, 0, 1)
  idx <- which(z > 0.5)
  y[idx] <- 1
  
  x <- matrix(rweibull(n * p, shape = 1, scale = 1), nrow = n, ncol = p)
  x[idx, ] <- sweep(x[idx, ], 2, MU2, FUN = "+")
  
  # --- Test Data ---
  y.test <- numeric(n.test)
  z.test <- runif(n.test, 0, 1)
  id <- which(z.test > 0.5)
  y.test[id] <- 1  
  
  x.test <- matrix(rweibull(n.test * p, shape = 1, scale = 1), nrow = n.test, ncol = p)
  x.test[id, ] <- sweep(x.test[id, ], 2, MU2, FUN = "+")
  
  list(x = x, y = y, x.test = x.test, y.test = y.test, MU2 = MU2)
}
do.simulation <- function(n, p, d, n.test, p.contam, RUNs, isNSTransformed, isOracle, Model, err.sigma){
  pi <- 0.5
  delta <- (2 * qnorm(0.10))^2
  
  error <- numeric(RUNs)
  TP <- numeric(RUNs)
  FP <- numeric(RUNs)
  FN <- numeric(RUNs)
  
  for(h in 1:RUNs){
    # Data Generation
    all.data <- Model(n, p, d, n.test, 1245*h)
    x <- all.data[["x"]]
    y <- all.data[["y"]]
    x.test <- all.data[["x.test"]]
    y.test <- all.data[["y.test"]]
    MU2 <- all.data[["MU2"]]
    
    # Contamination Loop
    x.uncontam <- x[]
    for(i in 1:n){
      idx.contam <- which(rbinom(p, 1, p.contam) == 1)
      if (length(idx.contam) > 0) {
        x[i, idx.contam] <- rnorm(length(idx.contam), 0, err.sigma)
      }
    }
    
    if(isNSTransformed){
      x.norm <- matrix(0, n, p)
      x.test.norm <- matrix(0, n.test, p)
      if(isOracle){
        for(i in 1:p){
          x.norm[, i] <- 0.5 * (pnorm(x.uncontam[, i]) + pnorm(x.uncontam[, i], mean=MU2[i], sd=1))
          x.test.norm[, i] <- 0.5 * (pnorm(x.test[, i]) + pnorm(x.test[, i], mean=MU2[i], sd=1))
        }
      }else{
        for(i in 1:p){
          obj.norm <- InvNorm(x[, i], c = 0)
          x.norm[, i] <- obj.norm$x.norm
          x.test.norm[, i] <- obj.norm$transform(x.test[, i])
        }
      }
    } else {
      x.norm <- x
      x.test.norm <- x.test
    }
    
    dat.train <- list(x = t(x.norm), y = y, geneid = 1:p)
    dat.test <- list(x = t(x.test.norm), y = y.test, geneid = 1:p)
    
    junk <- capture.output(obj.pam <- pamr.train(dat.train))
    set.seed(1256 * h)
    junk <- capture.output(obj.pamr.cv <- pamr.cv(obj.pam, dat.train))
    
    t_thresh <- obj.pamr.cv$threshold[which.min(obj.pamr.cv$error)]
    fit <- pamr.predict(obj.pam, t(x.test.norm), type = "class", threshold = t_thresh)
    fit <- as.numeric(fit) - 1
    error[h] <- mean(abs(fit - y.test))
    
    pam.sel <- rep(0,p)
    junk <- capture.output(pam.sel[as.numeric(pamr.listgenes(obj.pam,dat.train,threshold=t_thresh)[,1])] <- 1)
    
    Iv.true <- as.numeric(MU2 != 0)
    Iv.hat  <- as.numeric(pam.sel > 0)
    
    TP[h] <- sum(Iv.true == 1 & Iv.hat == 1)
    FP[h] <- sum(Iv.true == 0 & Iv.hat == 1)
    FN[h] <- sum(Iv.true == 1 & Iv.hat == 0)
  }
  
  cat("----------------------------------------------------------------------------------------------\n")
  cat("Error :", median(error), "(", sd(error * 100) / sqrt(RUNs), ")", "\n")
  cat("FN    :", median(FN), "(", sd(FN) / sqrt(RUNs), ")","\n")
  cat("FP    :", median(FP), "(", sd(FP) / sqrt(RUNs), ")","\n")
  cat("----------------------------------------------------------------------------------------------\n")
}


do.simulation2 <- function(n, p, d, n.test, p.contam, RUNs, isNSTransformed, isOracle, Model, err.shape){
  pi <- 0.5
  delta <- (2 * qnorm(0.10))^2
  
  error <- numeric(RUNs)
  TP <- numeric(RUNs)
  FP <- numeric(RUNs)
  FN <- numeric(RUNs)
  
  for(h in 1:RUNs){
    # Data Generation
    all.data <- Model(n, p, d, n.test, 1245*h)
    x <- all.data[["x"]]
    y <- all.data[["y"]]
    x.test <- all.data[["x.test"]]
    y.test <- all.data[["y.test"]]
    MU2 <- all.data[["MU2"]]
    
    # Contamination Loop
    x.uncontam <- x[]
    for(i in 1:n){
      idx.contam <- which(rbinom(p, 1, p.contam) == 1)
      if (length(idx.contam) > 0) {
        x[i, idx.contam] <- rweibull(length(idx.contam), shape = err.shape, scale=1)
      }
    }
    
    if(isNSTransformed){
      x.norm <- matrix(0, n, p)
      x.test.norm <- matrix(0, n.test, p)
      if(isOracle){
        for(i in 1:p){
          x.norm[, i] <- 0.5 * (pweibull(x.uncontam[, i], shape = 1, scale = 1) + pweibull(x.uncontam[, i] - MU2[i], shape = 1, scale = 1))
          x.test.norm[, i] <- 0.5 * (pweibull(x.test[, i], shape = 1, scale = 1) + pweibull(x.test[, i] - MU2[i], shape = 1, scale = 1))
        }
      }else{
        for(i in 1:p){
          obj.norm <- InvNorm(x[, i], c = 0)
          x.norm[, i] <- obj.norm$x.norm
          x.test.norm[, i] <- obj.norm$transform(x.test[, i])
        }
      }
    } else {
      x.norm <- x
      x.test.norm <- x.test
    }
    
    dat.train <- list(x = t(x.norm), y = y, geneid = 1:p)
    dat.test <- list(x = t(x.test.norm), y = y.test, geneid = 1:p)
    
    junk <- capture.output(obj.pam <- pamr.train(dat.train))
    set.seed(1256 * h)
    junk <- capture.output(obj.pamr.cv <- pamr.cv(obj.pam, dat.train))
    
    t_thresh <- obj.pamr.cv$threshold[which.min(obj.pamr.cv$error)]
    fit <- pamr.predict(obj.pam, t(x.test.norm), type = "class", threshold = t_thresh)
    fit <- as.numeric(fit) - 1
    error[h] <- mean(abs(fit - y.test))
    
    pam.sel <- rep(0,p)
    junk <- capture.output(pam.sel[as.numeric(pamr.listgenes(obj.pam,dat.train,threshold=t_thresh)[,1])] <- 1)
    
    Iv.true <- as.numeric(MU2 != 0)
    Iv.hat  <- as.numeric(pam.sel > 0)
    
    TP[h] <- sum(Iv.true == 1 & Iv.hat == 1)
    FP[h] <- sum(Iv.true == 0 & Iv.hat == 1)
    FN[h] <- sum(Iv.true == 1 & Iv.hat == 0)
  }
  
  cat("----------------------------------------------------------------------------------------------\n")
  cat("Error :", median(error), "(", sd(error * 100) / sqrt(RUNs), ")", "\n")
  cat("FN    :", median(FN), "(", sd(FN) / sqrt(RUNs), ")","\n")
  cat("FP    :", median(FP), "(", sd(FP) / sqrt(RUNs), ")","\n")
  cat("----------------------------------------------------------------------------------------------\n")
}
