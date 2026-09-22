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
library(mvnfast)



source("C:/Users/sudhi/Desktop/Coordinatewise Gaussianization/replication/List_functions_for_cellwise_corruption.R")
source("C:/Users/sudhi/Desktop/PhD projects/NS-Estimator/Lfunctions.R")


n <- 200
p <- 2000
RUNS <- 100


Epsilons <- rmvnorm_precision_n(n, p, diag(p))
bayes.scale <- 0.811


data.U1 <- Epsilons
data.U2 <- sweep(Epsilons, 2, c(bayes.scale*rep(1, d), numeric(p-d)), `+`)

