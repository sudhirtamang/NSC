
#===================================== DATA GENERATION ==============================
#====================================================================================

require(mvtnorm)
#require(mvnfast)
#source("Methods & Evaluation.R")
#====================== DC-SIS Example 1 ==============================

n = 200
p = 2000

rho.seq = c(0.5,0.8)

Sigma1 = outer(1:p,1:p,function(a,b) rho.seq[1]^abs(a-b))
Sigma2 = outer(1:p,1:p,function(a,b) rho.seq[2]^abs(a-b))



#set.seed(12345678)
#beta0 = (-1)^rbinom(4,1,0.4)*( 4*log(n)/sqrt(n) + abs(rnorm(4)) )


data.X <- function(n,case=c(1,2)){
    X = rmvnorm(n, sigma = switch(case,Sigma1,Sigma2), method="chol" )
    return(X)
}



data.y <- function(X, model = c(1,2,3,4)){
    
    n = dim(X)[1]
    p = dim(X)[2]
    const = c(2,0.5,3,2)
    
    signal = switch(model,
        cbind(X[,1], X[,2], as.numeric(X[,12]<0), X[,22]) %*% (const*beta0),
        cbind(X[,1]*X[,2], as.numeric(X[,12]<0),  X[,22]) %*% (const[c(1,3,4)]*beta0[1:3]),
        cbind(X[,1]*X[,2], as.numeric(X[,12]<0)* X[,22]) %*% (const[c(1,3)]*beta0[1:2]),
        cbind(X[,1], X[,2], as.numeric(X[,12]<0)) %*% (const[1:3]*beta0[1:3])           )
    
    epsilon = if (model != 4){
        rnorm(n)
    } else {
         exp(const[4]*abs(X[,22]))*rnorm(n)
    }
    
    return( signal + epsilon )
}



