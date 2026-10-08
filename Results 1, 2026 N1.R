# ---------------------------------------------------------
# Project: NSC
# Description: We want to see how well Nearest Shrunken Centroid do in the presence of contamination, the data being Normal Score
#              transformed.
# Author: Sudhir Tamang
# Date: September 22, 2026
# ---------------------------------------------------------
rm(list = ls())
library(Tlasso)
library(tensr)
library(glasso)
library(expm)
library(rTensor)
library(doParallel)
library(furrr)
library(purrr)
library(caret)
library(mvnfast)



source("C:/Users/sudhi/Desktop/Coordinatewise Gaussianization/replication/List_functions_for_cellwise_corruption.R")
source("C:/Users/sudhi/Desktop/PhD projects/NS-Estimator/Lfunctions.R")


n <- 200
n.test <- 10000
p <- 2000
RUNs <- 10
d <- 10
MU1 <- numeric(p)
MU2 <- c(rep(1, 10), numeric(p-d))


Error <- numeric(RUNs)
FP <- numeric(RUNs)
FN <- numeric(RUNs)
for(i in 1:RUNs){
  
set.seed(1234*i)
  train.Epsilons <- rmvnorm_precision_n(n/2, p, diag(p))
  test.Epsilons <- rmvnorm_precision_n(n.test/2, p, diag(p))
  bayes.scale <- 0.811
  contam.sigma <- sqrt(2)
  
  
  train.data.U1 <- train.Epsilons
  test.data.U1 <- test.Epsilons
  train.data.U2 <- sweep(train.Epsilons, 2, c(bayes.scale*rep(1, d), numeric(p-d)), `+`)
  test.data.U2 <- sweep(test.Epsilons, 2, c(bayes.scale*rep(1, d), numeric(p-d)), `+`)
  
  
  
  NS.transform <- function(data) {
    for (i in 1:ncol(data)) {
      data[, i] <- qnorm(rank(data[, i]) / (nrow(data) + 1))
    }
    data
  }
  
  train.data.U1.T <- NS.transform(train.data.U1)
  test.data.U1.T <- NS.transform(test.data.U1)
  train.data.U2.T <- NS.transform(train.data.U2)
  test.data.U2.T <- NS.transform(test.data.U2)
  train.data.U1.T <- train.data.U1
  test.data.U1.T <- test.data.U1
  train.data.U2.T <- train.data.U2
  test.data.U2.T <- test.data.U2
  
  # Now, contaminating the data
  # contam.cellwise <- function(data, pct.contam, contam.sigma){
  #   for(i in 1:nrow(data)){
  #     ind <- rbinom(ncol(data), 1, prob = pct.contam/100) == 1
  #     data[i, ind] <- rnorm(sum(ind), 10, contam.sigma)
  #   }
  #   data
  # }
  # 
  # train.data.X1 <- contam.cellwise(train.data.U1, 10, contam.sigma)
  # train.data.X2 <- contam.cellwise(train.data.U2, 10, contam.sigma)
  # 
  # train.data.X1.T <- NS.transform(train.data.X1)
  # train.data.X2.T <- NS.transform(train.data.X2)
  
  
  train.data <- rbind(train.data.U1.T, train.data.U2.T)
  labels.train <- factor(c(rep(0, n/2), rep(1, n/2)))
  colnames(train.data) <- 1:p
  # Set up cross-validation
  ctrl <- trainControl(method = "cv", number = 10)
  
  # Train Nearest Shrunken Centroid ('pam')
  fit_caret <- train(
    x = train.data, 
    y = labels.train,
    method = "pam",
    trControl = ctrl,
    tuneGrid = expand.grid(threshold = seq(0, 4, by = 0.5))
  )
  
  # # View optimal shrinkage threshold and accuracy
  # print(fit_caret)
  
  # Predict on new data directly (rows = new samples)
  test.data <- rbind(test.data.U1.T, test.data.U2.T)
  colnames(test.data) <- 1:p
  test.labels <- factor(c(rep(0, n.test/2), rep(1, n.test/2)))
  predictions <- predict(fit_caret, newdata = test.data)
  pam.obj <- fit_caret$finalModel
  t <- fit_caret$bestTune$threshold
  
  selected.idx <- as.numeric(
    pamr.listgenes(
      pam.obj,
      list(x = t(train.data), y = labels.train, geneid = 1:p),
      threshold = t
    )[,1]
  )
  
  pam.sel <- rep(0, p)
  pam.sel[selected.idx] <- 1
  
  false_positive <- function(actual, predicted) {
    sum(actual == 0 & predicted == 1)
  }
  
  false_negative <- function(actual, predicted) {
    sum(actual == 1 & predicted == 0)
  }
  
  FP[[i]] <- false_positive(MU2, pam.sel)
  FN[[i]] <- false_negative(MU2, pam.sel)
  Error[[i]] <- sum(test.labels != predictions)/n
}

cat("Printing Results!!!\n\n")
cat("Error Median: ", median(Error), "\n")
cat("FP Median: ", median(FP), "\n")
cat("FN Median: ", median(FN), "\n")
