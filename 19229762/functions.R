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
