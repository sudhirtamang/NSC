# ---------------------------------------------------------
# Project: NSC
# Description: Used the original code from the paper "Coordinatewise Gaussianization: Theories and Applications" by
#              Qing Mai, Di He & Hui Zou and added the part where the simulation is repeated the desired number of times(RUNs)
#              with some outliers being added
# Author: Sudhir Tamang
# Date: 18OCT2026
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

sink("Results 1 18OCT2026 Model N1.txt")
cat("===========================================================================================\n")
cat("-----------------------------START::: MODEL N1---------------------------------------\n")
cat("===========================================================================================\n")
# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N1, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs ", "\n")
do.simulationN1(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 124)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N1, NS-NSC ", p.contam*100, "%",   " Contamination", " With ", RUNs, " Runs ", "\n")
do.simulationN1(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 453)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- TRUE; isOracle <- TRUE
cat("Simulation for Model N1, Oracle ", p.contam*100, "%",   " Contamination", " With ", RUNs, " Runs ", "\n")
do.simulationN1(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 4132)
# --------------------------------------------------------------------------------------------------


cat("===========================================================================================\n")
cat("-----------------------------END::: MODEL N1---------------------------------------\n")
cat("===========================================================================================\n")
sink()

sink("Results 1 18OCT2026 Model N4.txt")
cat("===========================================================================================\n")
cat("-----------------------------START::: MODEL N4---------------------------------------\n")
cat("===========================================================================================\n")
# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N4, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs ", "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 13)
# --------------------------------------------------------------------------------------------------
# do.simulation <- function(n, p, d, n.test, p.contam, RUNs, MU2, isNSTransformed)

# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N4, NS-NSC ", p.contam*100, "%",   " Contamination", " With ", RUNs, " Runs ", "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 4564)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- TRUE; isOracle <- TRUE
cat("Simulation for Model N4, Oracle ", p.contam*100, "%",   " Contamination", " With ", RUNs, " Runs ", "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 132)
# --------------------------------------------------------------------------------------------------


cat("===========================================================================================\n")
cat("-----------------------------END::: MODEL N4---------------------------------------\n")
cat("===========================================================================================\n")
sink()




sink("Results 1 18OCT2026 Model N5.txt")
cat("===========================================================================================\n")
cat("-----------------------------START::: MODEL N5---------------------------------------\n")
cat("===========================================================================================\n")
# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N5, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs ", "\n")
do.simulationN5(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 13)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N5, NS-NSC ", p.contam*100, "%",   " Contamination", " With ", RUNs, " Runs ", "\n")
do.simulationN5(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 453)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- TRUE; isOracle <- TRUE
cat("Simulation for Model N5, Oracle ", p.contam*100, "%",   " Contamination", "\n")
do.simulationN5(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, 844)
# --------------------------------------------------------------------------------------------------

cat("===========================================================================================\n")
cat("-----------------------------END::: MODEL N5---------------------------------------\n")
cat("===========================================================================================\n")
sink()
