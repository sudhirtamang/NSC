
source("data.generation.R")
source("Methods & Evaluation.R")

#====================== DC-SIS Example 1 ==============================

n = 200 
p = 2000

set.seed(12345678)
beta0 = (-1)^rbinom(4,1,0.4)*( 4*log(n)/sqrt(n) + abs(rnorm(4)) )
beta0 = round(beta0,1)


##== We demonstrate the following screening procedure with only one
#run and report the minimum model size assuming we know the true
#predictors. In practice, we select the active set with the top $d$
#predictors (here we choose d=2*n/log(n)), since we do not know the
#true predictors. The comparing methods following are DCS,
#UDCS, NS-DCS and Wins-DCS in the paper, respectively. For a
#simulation study with more replicates, readers can wrap the following
#code with a loop and a slight modification, and report the 5%, 25%,
#50%, 75%, and 95% quantiles of the minimum model size.

X = data.X(n,1)
y = data.y(X, model = 1)


wk.DC = DC.wk(X,y)
min.model.size( wk.DC, true.pred = c(1,2,12,22) )
Active.set(wk.DC,2)

wk.DC.unif = DC.unif.wk(X,y)
min.model.size( wk.DC.unif, true.pred = c(1,2,12,22) )
Active.set(wk.DC.unif,2)

wk.DC.mns = DC.normal.wk(X,y, method="mns")
min.model.size( wk.DC.mns, true.pred = c(1,2,12,22) )
Active.set(wk.DC.mns,2)

wk.DC.wins = DC.normal.wk(X,y, method="wins")
min.model.size( wk.DC.wins, true.pred = c(1,2,12,22) )
Active.set(wk.DC.wins,2)















