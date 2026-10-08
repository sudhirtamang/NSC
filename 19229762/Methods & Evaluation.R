

#==============================================================================
#=================================   METHODS  =================================

#==================================== DC-SIS ==================================

require(energy)

DC.wk <- function(X,y){
    
    p = dim(X)[2]
    
    wk = numeric(p)
    for (k in 1:p){
        wk[k] = dcor(X[,k],y)
    }
    
    return(wk)
}


unif.tran <- function(X,max.catgr=10){
    X = as.matrix(X)
    n = dim(X)[1]
    p = dim(X)[2]

    catgr.id = NULL
    for (i in 1:p){
        if (length(unique(X[,i])) <= max.catgr){# max # of categorical variables
            catgr.id = c(catgr.id,i)
        }
    }

    if (is.null(catgr.id)){
        for (k in 1:p){
            temp = outer(X[,k],X[,k], function(a,b) a<=b )
            X[,k] = colSums(temp)/n
        }
    }else if(length(catgr.id) != p){
        for (k in setdiff(1:p,catgr.id)){
            temp = outer(X[,k],X[,k], function(a,b) a<=b )
            X[,k] = colSums(temp)/n
                }
    }
    
    return(X)
}


DC.unif.wk <- function(X,y){
    p = dim(X)[2]
    
    X = unif.tran(X)
    y = unif.tran(y)    
    
    wk = numeric(p)
    for (k in 1:p){
        wk[k] = dcor(X[,k],y)
    }
    
    return(wk)
}


mod.norm.scr <- function(X,max.catgr=10){
    X = as.matrix(X)
    n = dim(X)[1]
    p = dim(X)[2]

    X = unif.tran(X)*n/(n+1)

    catgr.id = NULL
    for (i in 1:p){
        if (length(unique(X[,i])) <= max.catgr){# max # of categorical variables
            catgr.id = c(catgr.id,i)
        }
    }

    if (is.null(catgr.id)){
        for (k in 1:p){
            X[,k] = qnorm(X[,k])
        }
    }else if(length(catgr.id) != p){
        for (k in setdiff(1:p,catgr.id)){
            X[,k] = qnorm(X[,k])
        }
    }
    
    return(X)
}


Winsorized.norm <- function(X,max.catgr=10){# skip transformation for categorical variables
    X = as.matrix(X)
    n = dim(X)[1]
    p = dim(X)[2]

    catgr.id = NULL
    for (i in 1:p){
        if (length(unique(X[,i])) < max.catgr){# max # of categorical variables
            catgr.id = c(catgr.id,i)
        }
    }

    if (is.null(catgr.id)){
        X = unif.tran(X)
        X[X==1] = 1 - 1/n
        for (k in 1:p){
            X[,k] = qnorm(X[,k])
        }
    }else if(length(catgr.id) != p){
        X[,-catgr.id] = unif.tran(X[,-catgr.id])
        X[,-catgr.id][X[,-catgr.id]==1] = 1 - 1/n
        for (k in setdiff(1:p,catgr.id)){
            X[,k] = qnorm(X[,k])
        }
    }
    
    return(X)
}


DC.normal.wk <- function(X,y,method=c("mns","wins")){
    p = dim(X)[2]

    if (method == "mns"){
        X = mod.norm.scr(X)
        y = mod.norm.scr(y)

        wk = numeric(p)
        for (k in 1:p){
            wk[k] = dcor(X[,k],y)
        }
        
    }else if (method == "wins"){
        X = Winsorized.norm(X)
        y = Winsorized.norm(y)

        wk = numeric(p)
        for (k in 1:p){
            wk[k] = dcor(X[,k],y)
        }
    }
    
    return(wk)
}

#==============================================================================
#================================   EVALUATION   ==============================
#==============================================================================

Active.set <- function(wk,di = c(1,2,3)){
    nid = round(di* floor(n/log(n))) # hard-thresholding
 #   nid = round(di*n) 
    if (length(unique(wk)) == 1){
        id = NULL
    }else{
        id = order(wk,decreasing = T)[1:nid]
    }
    return(id)
}


min.model.size <- function(wk,true.pred = c(1,2,12,22) )
{ # Maximum index of true predictor in wk
    wk_sort = sort(wk,decreasing = T)

    true.pred = sort(true.pred)
    n.true.pred = length(true.pred)
    id = numeric( n.true.pred )

    for(i in 1:n.true.pred){
        temp = which(wk_sort == wk[true.pred[i]])
        id[i] = temp[length(temp)]  # get last index when ties occur
    }
    return( max(id) )
}






















