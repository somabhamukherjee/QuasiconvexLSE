#Additional functions
#Do not edit. 
interpol.func.Dec<-function(X,cplex.fhat,nondatapoint)
{
  X<-as.matrix(X)
  z<-as.matrix(cplex.fhat)
  mat<-cbind(X,z)
  mat<-mat[order(mat[,ncol(mat)]),]
  out<-mat[nrow(mat),ncol(mat)]
  Y<-mat
  counter<-1
  if(nrow(Y)>1){
    counter<-chullupcheck(Y[,1:(ncol(Y)-1)],nondatapoint)
  }
  else{
    counter<-chullupcheck(t(as.matrix(Y[,1:(ncol(Y)-1)])),nondatapoint)
  }
  while(counter){
    out<-Y[nrow(Y),ncol(Y)]
    if(nrow(Y)>2) {
      Y<-Y[1:(nrow(Y)-1),]
      counter<-chullupcheck(Y[,1:(ncol(Y)-1)],nondatapoint)
    }
    else if(nrow(Y)==2){
      Y<-t(as.matrix(Y[1:(nrow(Y)-1),]))
      counter<-chullupcheck(t(as.matrix(Y[,1:(ncol(Y)-1)])),nondatapoint)
    }
    else{
      break
    }
  }
  return(out)
}  





interpol.func.Inc<-function(X,cplex.fhat,nondatapoint)
{
  X<-as.matrix(X)
  z<-as.matrix(cplex.fhat)
  mat<-cbind(X,z)
  mat<-mat[order(mat[,ncol(mat)]),]
  out<-mat[nrow(mat),ncol(mat)]
  Y<-mat
  counter<-1
  if(nrow(Y)>1){
    counter<-chulldowncheck(Y[,1:(ncol(Y)-1)],nondatapoint)
  }
  else{
    counter<-chulldowncheck(t(as.matrix(Y[,1:(ncol(Y)-1)])),nondatapoint)
  }
  while(counter){
    out<-Y[nrow(Y),ncol(Y)]
    if(nrow(Y)>2) {
      Y<-Y[1:(nrow(Y)-1),]
      counter<-chulldowncheck(Y[,1:(ncol(Y)-1)],nondatapoint)
    }
    else if(nrow(Y)==2){
      Y<-t(as.matrix(Y[1:(nrow(Y)-1),]))
      counter<-chulldowncheck(t(as.matrix(Y[,1:(ncol(Y)-1)])),nondatapoint)
    }
    else{
      break
    }
  }
  return(out)
}  





interpol.func<-function(X,cplex.fhat,nondatapoint)
{
  X<-as.matrix(X)
  z<-as.matrix(cplex.fhat)
  mat<-cbind(X,z)
  mat<-mat[order(mat[,ncol(mat)]),]
  out<-mat[nrow(mat),ncol(mat)]
  Y<-mat
  counter<-1
  if(nrow(Y)>1){
    counter<-chullcheck(Y[,1:(ncol(Y)-1)],nondatapoint)
  }
  else{
    counter<-chullcheck(t(as.matrix(Y[,1:(ncol(Y)-1)])),nondatapoint)
  }
  while(counter){
    out<-Y[nrow(Y),ncol(Y)]
    if(nrow(Y)>2) {
      Y<-Y[1:(nrow(Y)-1),]
      counter<-chullcheck(Y[,1:(ncol(Y)-1)],nondatapoint)
    }
    else if(nrow(Y)==2){
      Y<-t(as.matrix(Y[1:(nrow(Y)-1),]))
      counter<-chullcheck(t(as.matrix(Y[,1:(ncol(Y)-1)])),nondatapoint)
    }
    else{
      break
    }
  }
  return(out)
} 





QuasiConv.only.control <- function(X,y, Max.b , MM , ep) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.vector(y)) y <- as.vector(y)
  nn <- nrow(X)
  dd <- ncol(X)
  cvec<- c( -2*y, rep(0,  nn*dd+ nn^2))
  Qmat <- matrix(0, nrow= nn+ nn*dd+ nn^2, ncol =nn+ nn*dd+ nn^2)
  Qmat[1: nn, 1:nn] <- 2* diag(nn)
  Amat <- matrix(0,(2*nn*(nn-1)),(nn+(nn*dd)+nn^2))
  negid<- -diag(nn-1)
  jdmat<-cbind(matrix(1,nn-1,1),negid)
  for(i in 2:(nn-1))
  {
    jd<-cbind(negid[,1:(i-1)],matrix(1,nn-1,1),negid[,i:(nn-1)])
    jdmat<-rbind(jdmat,jd)
  }
  jd<-cbind(negid,matrix(1,nn-1,1))
  jdmat<-rbind(jdmat,jd)  
  Amat[1:(nn*(nn-1)),1:nn]<-jdmat
  for(alp in 1:nn)
  {
    for(bet in setdiff(c(1:nn),alp))
    {
      if(bet<alp)
      {
        Amat[(nn*(nn-1))+((nn-1)*(alp-1))+bet,(nn+(alp-1)*dd+1):(nn+alp*dd)]<-X[alp,]-X[bet,]
      }
      else
      {
        Amat[(nn*(nn-1))+((nn-1)*(alp-1))+bet-1,(nn+(alp-1)*dd+1):(nn+alp*dd)]<-X[alp,]-X[bet,]
      }
    }
  }
  Amat[1:(nn-1),(nn+nn*dd+1):(nn+nn*dd+nn)]<-cbind(matrix(0,nn-1,1),MM*negid)
  for(i in 2:(nn-1))
  {
    Amat[((i-1)*(nn-1)+1):(i*(nn-1)),((nn*dd)+(i*nn)+1):(nn*dd + i*nn + nn)]<-cbind(MM*negid[,1:(i-1)],matrix(0,nn-1,1),MM*negid[,i:(nn-1)])
  }
  i<-nn
  Amat[((i-1)*(nn-1)+1):(i*(nn-1)),((nn*dd)+(i*nn)+1):(nn*dd + i*nn + nn)]<- cbind(MM*negid,matrix(0,nn-1,1))
  Amat[(nn*(nn-1)+1):(2*nn*(nn-1)),(nn+nn*dd+1):(nn+nn*dd+nn^2)]<- - Amat[1:(nn*(nn-1)),(nn+nn*dd+1):(nn+nn*dd+nn^2)]
  bvec <-c(rep(0, nn*(nn-1)), rep(MM-ep, nn*(nn-1)))
  lb <- c(rep(-Max.b, nn), rep(-Max.b, nn*dd), rep(-1, nn^2))
  ub <-  c(rep(Max.b, nn+ nn*dd), rep(2, nn^2))
  vtype <- c(rep("C", nn+ nn*dd), rep("B", nn^2)) 
  list(cvec = cvec,
       Amat = Amat,
       lb = lb,
       ub = ub,
       bvec = bvec,
       vtype = vtype,
       Qmat = Qmat,
       dd = dd,
       varlength = length(cvec) )
}





QuasiConv.Dec.control <- function(X,y, Max.b , MM, ep ) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.vector(y)) y <- as.vector(y)
  nn <- nrow(X)
  dd <- ncol(X)
  cvec<- c( -2*y, rep(0,  nn*dd+ nn^2))
  Qmat <- matrix(0, nrow= nn+ nn*dd+ nn^2, ncol =nn+ nn*dd+ nn^2)
  Qmat[1: nn, 1:nn] <- 2* diag(nn)
  Amat <- matrix(0,(2*nn*(nn-1)),(nn+(nn*dd)+nn^2))
  negid<- -diag(nn-1)
  jdmat<-cbind(matrix(1,nn-1,1),negid)
  for(i in 2:(nn-1))
  {
    jd<-cbind(negid[,1:(i-1)],matrix(1,nn-1,1),negid[,i:(nn-1)])
    jdmat<-rbind(jdmat,jd)
  }
  jd<-cbind(negid,matrix(1,nn-1,1))
  jdmat<-rbind(jdmat,jd)  
  Amat[1:(nn*(nn-1)),1:nn]<-jdmat
  for(alp in 1:nn)
  {
    for(bet in setdiff(c(1:nn),alp))
    {
      if(bet<alp)
      {
        Amat[(nn*(nn-1))+((nn-1)*(alp-1))+bet,(nn+(alp-1)*dd+1):(nn+alp*dd)]<-X[alp,]-X[bet,]
      }
      else
      {
        Amat[(nn*(nn-1))+((nn-1)*(alp-1))+bet-1,(nn+(alp-1)*dd+1):(nn+alp*dd)]<-X[alp,]-X[bet,]
      }
    }
  }
  Amat[1:(nn-1),(nn+nn*dd+1):(nn+nn*dd+nn)]<-cbind(matrix(0,nn-1,1),MM*negid)
  for(i in 2:(nn-1))
  {
    Amat[((i-1)*(nn-1)+1):(i*(nn-1)),((nn*dd)+(i*nn)+1):(nn*dd + i*nn + nn)]<-cbind(MM*negid[,1:(i-1)],matrix(0,nn-1,1),MM*negid[,i:(nn-1)])
  }
  i<-nn
  Amat[((i-1)*(nn-1)+1):(i*(nn-1)),((nn*dd)+(i*nn)+1):(nn*dd + i*nn + nn)]<- cbind(MM*negid,matrix(0,nn-1,1))
  Amat[(nn*(nn-1)+1):(2*nn*(nn-1)),(nn+nn*dd+1):(nn+nn*dd+nn^2)]<- - Amat[1:(nn*(nn-1)),(nn+nn*dd+1):(nn+nn*dd+nn^2)]  
  bvec <-c(rep(0, nn*(nn-1)), rep(MM-ep, nn*(nn-1)))
  lb <- c(rep(-Max.b, nn), rep(0, nn*dd), rep(-1, nn^2))
  ub <-  c(rep(Max.b, nn+ nn*dd), rep(2, nn^2))
  vtype <- c(rep("C", nn+ nn*dd), rep("B", nn^2))
  list(cvec = cvec,
       Amat = Amat,
       lb = lb,
       ub = ub,
       bvec = bvec,
       vtype = vtype,
       Qmat = Qmat,
       dd = dd,
       varlength = length(cvec) )
}





QuasiConv.Inc.control <- function(X,y, Max.b , MM, ep ) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.vector(y)) y <- as.vector(y)
  nn <- nrow(X)
  dd <- ncol(X)
  cvec<- c( -2*y, rep(0,  nn*dd+ nn^2))
  Qmat <- matrix(0, nrow= nn+ nn*dd+ nn^2, ncol =nn+ nn*dd+ nn^2)
  Qmat[1: nn, 1:nn] <- 2* diag(nn)
  Amat <- matrix(0,(2*nn*(nn-1)),(nn+(nn*dd)+nn^2))
  negid<- -diag(nn-1)
  jdmat<-cbind(matrix(1,nn-1,1),negid)
  for(i in 2:(nn-1))
  {
    jd<-cbind(negid[,1:(i-1)],matrix(1,nn-1,1),negid[,i:(nn-1)])
    jdmat<-rbind(jdmat,jd)
  }
  jd<-cbind(negid,matrix(1,nn-1,1))
  jdmat<-rbind(jdmat,jd)  
  Amat[1:(nn*(nn-1)),1:nn]<-jdmat
  for(alp in 1:nn)
  {
    for(bet in setdiff(c(1:nn),alp))
    {
      if(bet<alp)
      {
        Amat[(nn*(nn-1))+((nn-1)*(alp-1))+bet,(nn+(alp-1)*dd+1):(nn+alp*dd)]<-X[alp,]-X[bet,]
      }
      else
      {
        Amat[(nn*(nn-1))+((nn-1)*(alp-1))+bet-1,(nn+(alp-1)*dd+1):(nn+alp*dd)]<-X[alp,]-X[bet,]
      }
    }
  }
  Amat[1:(nn-1),(nn+nn*dd+1):(nn+nn*dd+nn)]<-cbind(matrix(0,nn-1,1),MM*negid)
  for(i in 2:(nn-1))
  {
    Amat[((i-1)*(nn-1)+1):(i*(nn-1)),((nn*dd)+(i*nn)+1):(nn*dd + i*nn + nn)]<-cbind(MM*negid[,1:(i-1)],matrix(0,nn-1,1),MM*negid[,i:(nn-1)])
  }
  i<-nn
  Amat[((i-1)*(nn-1)+1):(i*(nn-1)),((nn*dd)+(i*nn)+1):(nn*dd + i*nn + nn)]<- cbind(MM*negid,matrix(0,nn-1,1))
  Amat[(nn*(nn-1)+1):(2*nn*(nn-1)),(nn+nn*dd+1):(nn+nn*dd+nn^2)]<- - Amat[1:(nn*(nn-1)),(nn+nn*dd+1):(nn+nn*dd+nn^2)]  
  bvec <-c(rep(0, nn*(nn-1)), rep(MM-ep, nn*(nn-1)))
  lb <- c(rep(-Max.b, nn+ nn*dd), rep(-1, nn^2))
  ub <-  c(rep(Max.b, nn), rep(0, nn*dd), rep(2, nn^2))
  vtype <- c(rep("C", nn+ nn*dd), rep("B", nn^2))
  list(cvec = cvec,
       Amat = Amat,
       lb = lb,
       ub = ub,
       bvec = bvec,
       vtype = vtype,
       Qmat = Qmat,
       dd = dd,
       varlength = length(cvec) )
}
