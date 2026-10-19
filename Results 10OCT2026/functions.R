# On Oct 06, 2026 added a function do.simulation to do as many RUNs of simulation with given percentage of contamination(p.contam)
# Number of obserbations on test set(n.test), sample size(n), number of predictors(p)

# This function.R is for folder:: Results 10OCT2026

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



do.simulationN1 <- function(n, p, d, n.test, p.contam, RUNs, isNSTransformed, isOracle, err.shape){
  pi <- 0.5
  delta <- (2 * qnorm(0.10))^2
  sigma <- diag(p) 
  MU2 <- c(rep(1, d), rep(0, p - d))
  quad_form <- as.numeric(t(MU2) %*% solve(sigma, MU2))
  a <- sqrt(delta / quad_form)
  MU2 <- a * MU2
  error <- numeric(RUNs)
  TP <- numeric(RUNs)
  TN <- numeric(RUNs)
  FP <- numeric(RUNs)
  FN <- numeric(RUNs)
  
  for(h in 1:RUNs){
    set.seed(1234*h)
    z <- runif(n, 0, 1)
    idx <- which(z > 0.5)
    
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    x <- rbind(sweep(x[-idx, , drop = FALSE], 2, MU2, FUN = "+"), x[idx, ])
    y <- c(rep(1, n - length(idx)), rep(0, length(idx)))

    
    # --- Test Data ---
    z.test <- runif(n.test, 0, 1)
    idx.test <- which(z.test > 0.5)
    
    x.test <- matrix(rnorm(n.test * p), nrow = n.test, ncol = p)
    x.test <- rbind(sweep(x.test[-idx.test, , drop = FALSE], 2, MU2, FUN = "+"), x.test[idx.test, ])
    y.test <- c(rep(1, n.test - length(idx.test)), rep(0, length(idx.test)))
    
    # Contamination Loop
    x.uncontam <- x[]
    for(i in 1:n){
      idx.contam <- which(rbinom(p, 1, p.contam) == 1)
      if (length(idx.contam) > 0) {
        x[i, idx.contam] <- rweibull(length(idx.contam), shape=err.shape, scale = 1)
      }
    }
    
    if(isNSTransformed){
      x.norm <- matrix(0, n, p)
      x.test.norm <- matrix(0, n.test, p)
      if(isOracle){
        for(i in 1:p){
          tmp1 <- pnorm(x.uncontam[, i], mean=0, sd=1) + pnorm(x.uncontam[, i], mean=MU2[i], sd=1)
          x.norm[, i] <- qnorm(0.5 * tmp1 * (n/(n+1)))
          tmp1 <- pnorm(x.test[, i], mean=0, sd=1) + pnorm(x.test[, i], mean=MU2[i], sd=1)
          x.test.norm[, i] <- qnorm(0.5 * tmp1 * (n/(n+1)))
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
    gene_list <- tryCatch({
      junk <- capture.output(res <- pamr.listgenes(obj.pam, dat.train, threshold = t_thresh))
      res
    }, error = function(e) {
      NULL 
    })
    if (!is.null(gene_list) && is.matrix(gene_list) && nrow(gene_list) > 0) {
      gene_indices <- as.numeric(gene_list[, 1])
      pam.sel[gene_indices] <- 1
    } else {
      pam.sel[] <- 0 
    }
    
    Iv.true <- as.numeric(MU2 != 0)
    Iv.hat  <- as.numeric(pam.sel > 0)
    
    TP[h] <- sum(Iv.true == 1 & Iv.hat == 1)
    TN[h] <- sum(Iv.true == 0 & Iv.hat == 0)
    FP[h] <- sum(Iv.true == 0 & Iv.hat == 1)
    FN[h] <- sum(Iv.true == 1 & Iv.hat == 0)
  }
  
  cat("----------------------------------------------------------------------------------------------\n")
  cat("Error :", median(error * 100), "(", sd(error * 100) / sqrt(RUNs), ")", "\n")
  cat("FN    :", median(FN), "(", sd(FN) / sqrt(RUNs), ")","\n")
  cat("FP    :", median(FP), "(", sd(FP) / sqrt(RUNs), ")","\n")
  cat("TN    :", median(TN), "(", sd(TN) / sqrt(RUNs), ")","\n")
  cat("TP    :", median(TP), "(", sd(TP) / sqrt(RUNs), ")","\n")
  cat("----------------------------------------------------------------------------------------------\n")
}



do.simulationN5 <- function(n, p, d, n.test, p.contam, RUNs, isNSTransformed, isOracle, err.shape){
  MU2 <- c(rep(1, d), rep(0, p - d))
  sigma <- 0.5 ^ abs(outer(1:p, 1:p, "-"))
  R <- chol(sigma)
  
  tmp1 <- crossprod(MU2, solve(sigma, MU2))
  a <- -2 * drop(qnorm(0.1) / sqrt(tmp1))
  MU2 <- a * MU2
  error <- numeric(RUNs)
  TP <- numeric(RUNs)
  TN <- numeric(RUNs)
  FP <- numeric(RUNs)
  FN <- numeric(RUNs)
  
  for(h in 1:RUNs){
    set.seed(1234*h)
    z <- runif(n, 0, 1)
    idx <- which(z > 0.5)
    x <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% R
    x <- rbind(sweep(x[-idx, , drop = FALSE], 2, MU2, FUN = "+"), x[idx, ])
    y <- c(rep(1, n - length(idx)), rep(0, length(idx)))
    
    # --- Test Data ---
    z.test <- runif(n.test, 0, 1)
    idx.test <- which(z.test > 0.5)
    x.test <- matrix(rnorm(n.test * p), nrow = n.test, ncol = p) %*% R
    x.test <- rbind(sweep(x.test[-idx.test, , drop = FALSE], 2, MU2, FUN = "+"), x.test[idx.test, ])
    y.test <- c(rep(1, n.test - length(idx.test)), rep(0, length(idx.test)))
    
    # Contamination Loop
    x.uncontam <- x[]
    for(i in 1:n){
      idx.contam <- which(rbinom(p, 1, p.contam) == 1)
      if (length(idx.contam) > 0) {
        x[i, idx.contam] <- rweibull(length(idx.contam), shape=err.shape, scale = 1)
      }
    }
    
    if(isNSTransformed){
      x.norm <- matrix(0, n, p)
      x.test.norm <- matrix(0, n.test, p)
      if(isOracle){
        for(i in 1:p){
          tmp1 <- pnorm(x.uncontam[, i], mean=0, sd=1) + pnorm(x.uncontam[, i], mean=MU2[i], sd=1)
          x.norm[, i] <- qnorm(0.5 * tmp1 * (n/(n+1)))
          tmp1 <- pnorm(x.test[, i], mean=0, sd=1) + pnorm(x.test[, i], mean=MU2[i], sd=1)
          x.test.norm[, i] <- qnorm(0.5 * tmp1 * (n/(n+1)))
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
    gene_list <- tryCatch({
      junk <- capture.output(res <- pamr.listgenes(obj.pam, dat.train, threshold = t_thresh))
      res
    }, error = function(e) {
      NULL 
    })
    if (!is.null(gene_list) && is.matrix(gene_list) && nrow(gene_list) > 0) {
      gene_indices <- as.numeric(gene_list[, 1])
      pam.sel[gene_indices] <- 1
    } else {
      pam.sel[] <- 0 
    }
    
    Iv.true <- as.numeric(MU2 != 0)
    Iv.hat  <- as.numeric(pam.sel > 0)
    
    TP[h] <- sum(Iv.true == 1 & Iv.hat == 1)
    TN[h] <- sum(Iv.true == 0 & Iv.hat == 0)
    FP[h] <- sum(Iv.true == 0 & Iv.hat == 1)
    FN[h] <- sum(Iv.true == 1 & Iv.hat == 0)
  }
  
  cat("----------------------------------------------------------------------------------------------\n")
  cat("Error :", median(error * 100), "(", sd(error * 100) / sqrt(RUNs), ")", "\n")
  cat("FN    :", median(FN), "(", sd(FN) / sqrt(RUNs), ")","\n")
  cat("FP    :", median(FP), "(", sd(FP) / sqrt(RUNs), ")","\n")
  cat("TN    :", median(TN), "(", sd(TN) / sqrt(RUNs), ")","\n")
  cat("TP    :", median(TP), "(", sd(TP) / sqrt(RUNs), ")","\n")
  cat("----------------------------------------------------------------------------------------------\n")
}


do.simulationN4 <- function(n, p, d, n.test, p.contam, RUNs, isNSTransformed, isOracle, err.shape){
  # Ensures length of MU2 is exactly p
  MU2 <- 0.5 * c(rep(1, d), numeric(p - d))

  error <- numeric(RUNs)
  TP <- numeric(RUNs)
  FP <- numeric(RUNs)
  FN <- numeric(RUNs)
  TN <- numeric(RUNs)
  
  for(h in 1:RUNs){
    set.seed(1234*h)
    z <- runif(n, 0, 1)
    idx <- which(z > 0.5)
    
    x <- matrix(rweibull(n * p, shape = 1, scale = 1), nrow = n, ncol = p)
    x <- rbind(sweep(x[-idx, , drop = FALSE], 2, MU2, FUN = "+"), x[idx, ])
    y <- c(rep(1, n - length(idx)), rep(0, length(idx)))
    
    # --- Test Data ---
    z.test <- runif(n.test, 0, 1)
    idx.test <- which(z.test > 0.5)
    
    x.test <- matrix(rweibull(n.test * p, shape = 1, scale = 1), nrow = n.test, ncol = p)
    x.test <- rbind(sweep(x.test[-idx.test, , drop = FALSE], 2, MU2, FUN = "+"), x.test[idx.test, ])
    y.test <- c(rep(1, n.test - length(idx.test)), rep(0, length(idx.test)))
    
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
          tmp1 <- pweibull(x.uncontam[, i], shape=1, scale=1) + pweibull(x.uncontam[, i] - MU2[i], shape=1, scale=1)
          x.norm[, i] <- qnorm(0.5 * tmp1 * (n/(n+1)))
          tmp1 <- pweibull(x.test[, i], shape=1, scale=1) + pweibull(x.test[, i] - MU2[i], shape=1, scale=1)
          x.test.norm[, i] <- qnorm(0.5 * tmp1 * (n/(n+1)))
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
    gene_list <- tryCatch({
      junk <- capture.output(res <- pamr.listgenes(obj.pam, dat.train, threshold = t_thresh))
      res
    }, error = function(e) {
      NULL 
    })
    if (!is.null(gene_list) && is.matrix(gene_list) && nrow(gene_list) > 0) {
      gene_indices <- as.numeric(gene_list[, 1])
      pam.sel[gene_indices] <- 1
    } else {
      pam.sel[] <- 0 
    }
    
    Iv.true <- as.numeric(MU2 != 0)
    Iv.hat  <- as.numeric(pam.sel > 0)
    
    TP[h] <- sum(Iv.true == 1 & Iv.hat == 1)
    FP[h] <- sum(Iv.true == 0 & Iv.hat == 1)
    FN[h] <- sum(Iv.true == 1 & Iv.hat == 0)
  }
  
  cat("----------------------------------------------------------------------------------------------\n")
  cat("Error :", median(error * 100), "(", sd(error * 100) / sqrt(RUNs), ")", "\n")
  cat("FN    :", median(FN), "(", sd(FN) / sqrt(RUNs), ")","\n")
  cat("FP    :", median(FP), "(", sd(FP) / sqrt(RUNs), ")","\n")
  cat("TN    :", median(TN), "(", sd(TN) / sqrt(RUNs), ")","\n")
  cat("TP    :", median(TP), "(", sd(TP) / sqrt(RUNs), ")","\n")
  cat("----------------------------------------------------------------------------------------------\n")
}
