
parametric.homothetic.func<-function(x){
  a<-matrix(rep(0,nrow(x)),nrow(x),1)
  for(i in 1:nrow(x)){
    a[i]<-((x[i,1])*(x[i,2]))^(0.5)
    a[i]<-15/(1+exp(-5*log(a[i])))
  }
  return(a)
}

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nonparametric.homothetic.func<-function(x){
  si<-1.51
  be<-0.45
  a<-matrix(rep(0,nrow(x)),nrow(x),1)
  for(i in 1:nrow(x)){      
    a[i]<-((be*((x[i,1])^((si-1)/si))) + ((1-be)*((x[i,2])^((si-1)/si))))^(si/(si-1))
    a[i]<-15/(1+exp(-5*log(a[i])))
  }
return(a)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


nonparametric.nonhomothetic.func<-function(x){
  si<-1.51
  bet<-function(y){
    return((0.25*(y/15))+0.3)
  }
  comp_func3<-function(be,yy){
    ac<-((be*((yy[1])^((si-1)/si))) + ((1-be)*((yy[2])^((si-1)/si))))^(si/(si-1))
    return(15/(1+exp(-5*log(ac))))
  }
  a<-matrix(rep(0,nrow(x)),nrow(x),1)
   storeyc<-seq(0,15,by=0.01)
    minvec<-seq(0,15,by=0.01)
  for(i in 1:nrow(x)){
    print(i)
    for (ix in 1:length(storeyc)){
      minvec[ix]<-(storeyc[ix]-comp_func3(bet(storeyc[ix]),x[i,]))^2
    }
    ystar<-storeyc[which.min(minvec)]
    a[i]<-comp_func3(bet(ystar),c(x[i,1],x[i,2]))
  }
  return(a)
}


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nonparametric.nonhomothetic<-function(sample.size,error.variance)
{
  si<-1.51
  bet<-function(y)
  {
    return((0.25*(y/15))+0.3)
  }
  comp_func3<-function(be,x)
  {
    a<-((be*((x[1])^((si-1)/si))) + ((1-be)*((x[2])^((si-1)/si))))^(si/(si-1))
    return(15/(1+exp(-5*log(a))))
  }
  X<-matrix(0,sample.size,2)
  y<-rep(0,sample.size)
  z<-rep(0,sample.size)
  for(i in 1:(sample.size))
  {
    psi<-(2.5)*runif(1)
    eta<-(0.05)+((pi/2)-0.1)*runif(1)
    X[i,1]<-psi*cos(eta)
    X[i,2]<-psi*sin(eta)
    minvec<-seq(0,15,by=0.01)
    storeyc<-seq(0,15,by=0.01)
    for (ix in 1:length(storeyc))
    {
      minvec[ix]<-(storeyc[ix]-comp_func3(bet(storeyc[ix]),X[i,]))^2
    }
    ystar<-storeyc[which.min(minvec)]
    # if(i==1)
    # {
    #   print(minvec)
    # }
    #print(c(i,ystar-comp_func3(bet(ystar),X[i,])))
    z[i] <- comp_func3(bet(ystar),c(X[i,1],X[i,2]))
    y[i]<- z[i]+rnorm(1,mean=0,sd=sqrt(error.variance))
  }
 ret<- cbind(as.matrix(X),as.matrix(y)) #, as.matrix(z))
 return(ret)
}

## Comments please
nonparametric.homothetic<-function(sample.size,error.variance)
{
  si<-1.51
  be<-0.45
  comp_func2<-function(x)
  {
    a<-((be*((x[1])^((si-1)/si))) + ((1-be)*((x[2])^((si-1)/si))))^(si/(si-1))
    return(15/(1+exp(-5*log(a))))
  }
  X<-matrix(0,sample.size,2)
  y<-rep(0,sample.size)
  for(i in 1:(sample.size))
  {
    psi<-(2.5)*runif(1)
    eta<-(0.05)+((pi/2)-0.1)*runif(1)
    X[i,1]<-psi*cos(eta)
    X[i,2]<-psi*sin(eta)
    y[i]<-comp_func2(c(X[i,1],X[i,2]))+rnorm(1,mean=0,sd=sqrt(error.variance))
  }
 ret<- cbind(as.matrix(X),as.matrix(y))
 return(ret)
}

#In Andrew's paper, this function was run 9 times, with sample.size=(100,500,1000) and error.variance=(1,4,9).
parametric.homothetic<-function(sample.size,error.variance)
{
  comp_func<-function(x)
  {
    a<-((x[1])*(x[2]))^(0.5)
    return(15/(1+exp(-5*log(a))))
  }
  X<-matrix(0,sample.size,2)
  y<-rep(0,sample.size)
  for(i in 1:(sample.size))
  {
    psi<-(2.5)*runif(1)
    eta<-(0.05)+((pi/2)-0.1)*runif(1)
    X[i,1]<-psi*cos(eta)
    X[i,2]<-psi*sin(eta)
    y[i]<-comp_func(c(X[i,1],X[i,2]))+rnorm(1,mean=0,sd=sqrt(error.variance))
  }
 ret<- cbind(as.matrix(X),as.matrix(y))
 return(ret)
}

