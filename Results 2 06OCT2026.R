# ---------------------------------------------------------
# Project: NSC
# Description: Used the original code from the paper "Coordinatewise Gaussianization: Theories and Applications" by
#              Qing Mai, Di He & Hui Zou and added the part where the simulation is repeated the desired number of times(RUNs)
#              with some outliers being added
# Author: Sudhir Tamang
# Date: 06OCT2026
# ---------------------------------------------------------


rm(list = ls())
library(mvtnorm)
library(MASS)
library(pamr)
source("functions.R")




n <- 200
p <- 2000
d <- 10
n.test <- 10000
RUNs <- 500

# sink("Results 2 06OCT2026 Model N1.txt")
# cat("===========================================================================================\n")
# cat("-----------------------------START::: MODEL N1---------------------------------------\n")
# cat("===========================================================================================\n")
# # --------------------------------------------------------------------------------------------------
# # p.contam <- 0; isNSTransformed <- FALSE; isOracle <- FALSE
# # cat("Simulation for Model N1, Raw ", p.contam*100,  "%",  " Contamination", "\n")
# # do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN1, sqrt(5))
# # # --------------------------------------------------------------------------------------------------
# # # do.simulation <- function(n, p, d, n.test, p.contam, RUNs, MU2, isNSTransformed)
# # 
# # # --------------------------------------------------------------------------------------------------
# # p.contam <- 0; isNSTransformed <- TRUE; isOracle <- FALSE
# # cat("Simulation for Model N1, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
# # do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN1, sqrt(5))
# # # --------------------------------------------------------------------------------------------------
# # 
# # # --------------------------------------------------------------------------------------------------
# # p.contam <- 0; isNSTransformed <- TRUE; isOracle <- TRUE
# # cat("Simulation for Model N1, Oracle ", p.contam*100, "%",   " Contamination", "\n")
# # do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN1, sqrt(5))
# # --------------------------------------------------------------------------------------------------
# 
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0.05; isNSTransformed <- FALSE; isOracle <- FALSE
# cat("Simulation for Model N1, Raw ", p.contam*100,  "%",  " Contamination", "\n")
# do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN1, sqrt(5))
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0.05; isNSTransformed <- TRUE; isOracle <- FALSE
# cat("Simulation for Model N1, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
# do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN1, sqrt(5))
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# # p.contam <- 0.05; isNSTransformed <- TRUE; isOracle <- TRUE
# # cat("Simulation for Model N1, Oracle ", p.contam*100, "%",   " Contamination", "\n")
# # do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, MU2, isNSTransformed, isOracle)
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0.1; isNSTransformed <- FALSE; isOracle <- FALSE
# cat("Simulation for Model N1, Raw ", p.contam*100,  "%",  " Contamination", "\n")
# do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN1, sqrt(5))
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0.1; isNSTransformed <- TRUE; isOracle <- FALSE
# cat("Simulation for Model N1, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
# do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN1, sqrt(5))
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# # p.contam <- 0.1; isNSTransformed <- TRUE; isOracle <- TRUE
# # cat("Simulation for Model N1, Oracle ", p.contam*100, "%",   " Contamination", "\n")
# # do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, MU2, isNSTransformed, isOracle)
# # --------------------------------------------------------------------------------------------------
# cat("===========================================================================================\n")
# cat("-----------------------------END::: MODEL N1---------------------------------------\n")
# cat("===========================================================================================\n")
# sink()
# 
# sink("Results 2 06OCT2026 Model N4.txt")
# cat("===========================================================================================\n")
# cat("-----------------------------START::: MODEL N4---------------------------------------\n")
# cat("===========================================================================================\n")
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0; isNSTransformed <- FALSE; isOracle <- FALSE
# cat("Simulation for Model N4, Raw ", p.contam*100,  "%",  " Contamination", "\n")
# do.simulation2(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN4, err.shape=0.75)
# # --------------------------------------------------------------------------------------------------
# # do.simulation <- function(n, p, d, n.test, p.contam, RUNs, MU2, isNSTransformed)
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0; isNSTransformed <- TRUE; isOracle <- FALSE
# cat("Simulation for Model N4, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
# do.simulation2(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN4, err.shape=0.75)
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0; isNSTransformed <- TRUE; isOracle <- TRUE
# cat("Simulation for Model N4, Oracle ", p.contam*100, "%",   " Contamination", "\n")
# do.simulation2(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN4, err.shape=0.75)
# # --------------------------------------------------------------------------------------------------
# 
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0.05; isNSTransformed <- FALSE; isOracle <- FALSE
# cat("Simulation for Model N4, Raw ", p.contam*100,  "%",  " Contamination", "\n")
# do.simulation2(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN4, err.shape=0.75)
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0.05; isNSTransformed <- TRUE; isOracle <- FALSE
# cat("Simulation for Model N4, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
# do.simulation2(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN4, err.shape=0.75)
# # --------------------------------------------------------------------------------------------------
# 
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0.1; isNSTransformed <- FALSE; isOracle <- FALSE
# cat("Simulation for Model N4, Raw ", p.contam*100,  "%",  " Contamination", "\n")
# do.simulation2(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN4, err.shape=0.75)
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0.1; isNSTransformed <- TRUE; isOracle <- FALSE
# cat("Simulation for Model N4, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
# do.simulation2(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN4, err.shape=0.75)
# # --------------------------------------------------------------------------------------------------
# 
# cat("===========================================================================================\n")
# cat("-----------------------------END::: MODEL N4---------------------------------------\n")
# cat("===========================================================================================\n")
# sink()




sink("Results 2 06OCT2026 Model N5.txt")
cat("===========================================================================================\n")
cat("-----------------------------START::: MODEL N5---------------------------------------\n")
cat("===========================================================================================\n")
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0; isNSTransformed <- FALSE; isOracle <- FALSE
# cat("Simulation for Model N5, Raw ", p.contam*100,  "%",  " Contamination", "\n")
# do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN5, sqrt(5))
# # --------------------------------------------------------------------------------------------------
# # do.simulation <- function(n, p, d, n.test, p.contam, RUNs, MU2, isNSTransformed)
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0; isNSTransformed <- TRUE; isOracle <- FALSE
# cat("Simulation for Model N5, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
# do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN5, sqrt(5))
# # --------------------------------------------------------------------------------------------------
# 
# # --------------------------------------------------------------------------------------------------
# p.contam <- 0; isNSTransformed <- TRUE; isOracle <- TRUE
# cat("Simulation for Model N5, Oracle ", p.contam*100, "%",   " Contamination", "\n")
# do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN5, sqrt(5))
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
p.contam <- 0.05; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N5, Raw ", p.contam*100,  "%",  " Contamination", "\n")
do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN5, sqrt(5))
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0.05; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N5, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN5, sqrt(5))
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
p.contam <- 0.1; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N5, Raw ", p.contam*100,  "%",  " Contamination", "\n")
do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN5, sqrt(5))
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0.1; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N5, NS-NSC ", p.contam*100, "%",   " Contamination", "\n")
do.simulation(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, ModelN5, sqrt(5))
# --------------------------------------------------------------------------------------------------

cat("===========================================================================================\n")
cat("-----------------------------END::: MODEL N5---------------------------------------\n")
cat("===========================================================================================\n")
sink()
