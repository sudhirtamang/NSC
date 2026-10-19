# ---------------------------------------------------------
# Project: NSC
# Description: Used the original code from the paper "Coordinatewise Gaussianization: Theories and Applications" by
#              Qing Mai, Di He & Hui Zou and added the part where the simulation is repeated the desired number of times(RUNs)
#              with some outliers being added. The Normal contamination has variance 2 in this simulation.
# Author: Sudhir Tamang
# Date: 12OCT2026
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

sink("Results 3 06OCT2026 Model N1.txt")
cat("===========================================================================================\n")
cat("-----------------------------START::: MODEL N1---------------------------------------\n")
cat("===========================================================================================\n")



# --------------------------------------------------------------------------------------------------
p.contam <- 0.05; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N1, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs", "\n")
do.simulationN1(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, sqrt(2))
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0.05; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N1, NS-NSC ", p.contam*100, "%",   " Contamination"," With ", RUNs, " Runs",  "\n")
do.simulationN1(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, sqrt(2))
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
p.contam <- 0.1; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N1, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs", "\n")
do.simulationN1(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, sqrt(2))
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0.1; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N1, NS-NSC ", p.contam*100, "%",   " Contamination"," With ", RUNs, " Runs",  "\n")
do.simulationN1(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, sqrt(2))
# --------------------------------------------------------------------------------------------------

cat("===========================================================================================\n")
cat("-----------------------------END::: MODEL N1---------------------------------------\n")
cat("===========================================================================================\n")
sink()

sink("Results 3 06OCT2026 Model N4.txt")
cat("===========================================================================================\n")
cat("-----------------------------START::: MODEL N4---------------------------------------\n")
cat("===========================================================================================\n")
# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N4, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs", "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, err.shape=0.95)
# --------------------------------------------------------------------------------------------------
# do.simulation <- function(n, p, d, n.test, p.contam, RUNs, MU2, isNSTransformed)

# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N4, NS-NSC ", p.contam*100, "%",   " Contamination"," With ", RUNs, " Runs",  "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, err.shape=0.95)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0; isNSTransformed <- TRUE; isOracle <- TRUE
cat("Simulation for Model N4, Oracle ", p.contam*100, "%",   " Contamination"," With ", RUNs, " Runs",  "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, err.shape=0.95)
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
p.contam <- 0.05; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N4, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs", "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, err.shape=0.95)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0.05; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N4, NS-NSC ", p.contam*100, "%",   " Contamination"," With ", RUNs, " Runs",  "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, err.shape=0.95)
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
p.contam <- 0.1; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N4, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs", "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, err.shape=0.95)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0.1; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N4, NS-NSC ", p.contam*100, "%",   " Contamination"," With ", RUNs, " Runs",  "\n")
do.simulationN4(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, err.shape=0.95)
# --------------------------------------------------------------------------------------------------

cat("===========================================================================================\n")
cat("-----------------------------END::: MODEL N4---------------------------------------\n")
cat("===========================================================================================\n")
sink()




sink("Results 3 06OCT2026 Model N5.txt")
cat("===========================================================================================\n")
cat("-----------------------------START::: MODEL N5---------------------------------------\n")
cat("===========================================================================================\n")



# --------------------------------------------------------------------------------------------------
p.contam <- 0.05; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N5, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs", "\n")
do.simulationN5(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, sqrt(2))
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0.05; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N5, NS-NSC ", p.contam*100, "%",   " Contamination"," With ", RUNs, " Runs",  "\n")
do.simulationN5(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, sqrt(2))
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
p.contam <- 0.1; isNSTransformed <- FALSE; isOracle <- FALSE
cat("Simulation for Model N5, Raw ", p.contam*100,  "%",  " Contamination", " With ", RUNs, " Runs", "\n")
do.simulationN5(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, sqrt(2))
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
p.contam <- 0.1; isNSTransformed <- TRUE; isOracle <- FALSE
cat("Simulation for Model N5, NS-NSC ", p.contam*100, "%",   " Contamination"," With ", RUNs, " Runs",  "\n")
do.simulationN5(n=n, p=p, d=d, n.test=n.test, p.contam=p.contam, RUNs=RUNs, isNSTransformed, isOracle, sqrt(2))
# --------------------------------------------------------------------------------------------------

cat("===========================================================================================\n")
cat("-----------------------------END::: MODEL N5---------------------------------------\n")
cat("===========================================================================================\n")
sink()
