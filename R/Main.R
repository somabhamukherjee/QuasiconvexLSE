#' Computes quasiconvex and decreasing least squares estimate using CPLEX .
#' This function will take in the multivariate regressors and the response, and
#' return the least squares estimates at the regressors, based on the response.
#' X,y and Monotone must be supplied, though all other parameters have pre-set values
#' the user can proceed with unless they wish to change the prior specification.
#'
#' @param X              An n by d matrix of regressors.
#' @param y              An n by 1 vector of responses.
#' @param Monotone       A categorical variable indicating the type of regression. 
#'                       The user must input one of the following three types:
#'                       "non.inc": means nonincreasing, quasiconvex regression,
#'                       "non.dec": means nondecreasing, quasiconvex regression,
#'                       "non": means quasiconvex regression.
#' @param Max.b          Bound on the absolute values of the variables.
#' @param MM             The 'M' parameter in the big-M method applied by us.
#' @param ep             A small positive quantity to convert open constraints to close constraints.
#' @param time.limit     Time limit in seconds of call to optimizer. Can be any nonnegative number.   
#' @param tol            Relative MIP optimality gap tolerance. Can be any nonnegative number. Default 1e-4.
#' @param output.flag    Turn CPLEX output on (1) or off(0).           
#' 
#' @return               A list of objects that contain among other things,the least squares estimates at
#'                       the regressors (within f.hat) and negative of the optimal objective function
#'                       (within least.squares).  
#'                       
#' @examples
#' n=20
#' d=8
#' library(MASS)
#' ## We demonstrate an example, where the regressor matrix X is Gaussian, and the response vector y
#' ## is a noisy version of a real-valued quasiconvex, decreasing function applied to the rows of X.
#' X = mvrnorm(n,rep(0,d), diag(d), tol=1e-06,empirical=FALSE)
#' y = exp(-rowSums(X)) + rnorm(n)
#'
#' ret = Quasiconvex.cplex(X, y,   Monotone = "non.inc", Max.b = 10^4, MM =10^4, ep =0.001,  time.limit= 200, tol =  1e-06,  output.flag = 1)
Quasiconvex.cplex <- function(X, y,   Monotone = c("non.inc", "non.dec","none"), Max.b = 10^4, MM =10^4, ep =0.001,  time.limit= 200, tol =  1e-06,  output.flag = 1){
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.vector(y)) y <- as.vector(y)
  dd <- ncol(X)
  if( Monotone == "non.inc"){
    Cons <-  QuasiConv.Dec.control(X, y, Max.b = Max.b, MM = MM, ep = ep )
  }
  else if(Monotone=="none") {
    Cons <-  QuasiConv.only.control(X, y, Max.b = Max.b, MM = MM, ep = ep )
  }
  else if(Monotone=="non.dec"){
    Cons <-  QuasiConv.Inc.control(X, y, Max.b = Max.b, MM = MM, ep = ep )
  }
  
  t <- Sys.time()
  out <- Rcplex(cvec = Cons$cvec, Amat = Cons$Amat, bvec = Cons$bvec, Qmat = Cons$Qmat, lb = Cons$lb, ub = Cons$ub, objsense = "min",sense = "L", vtype = Cons$vtype, control = list(tilim = time.limit ,epgap = tol, trace = output.flag)) 
  opt.fhat <- out$xopt[(1:nrow(X))] 
  opt.zeta <- out$xopt[(nrow(X)+1) : (nrow(X)+ nrow(X)*dd)] 
  
  t1 <- Sys.time()
  ret <- list( out = out,
               f.hat  =  opt.fhat,
               zeta.hat = opt.zeta,
               least.squares  =  -out$obj,
               time  =  t1-t,
               tol =  tol,
               method = "cplex",
               monotone = Monotone,
               status = out$status);
  ret$call <- match.call()
  class(ret) <- "Quasi.convex"
  return(ret)
}






#'Checks whether the points in the vector y are actually the outputs of a single quasiconvex
#'function applied to the corresponding rows of the matrix X. 
#'This function will take the matrix X and the vector y as inputs, and return a string, 
#'indicating whether the vector y is or not a quasiconvex realization of X.
#'
#' @param X              An n by d matrix.
#' @param y              An n by 1 vector.
#' @param tol            A small positive quantity.
#' 
#' @return               A string, indicating whether the vector y is, or not 
#'                       a quasiconvex realization of X.
#' @examples
#' library(MASS)
#' ## We demonstrate an example, where the regressor matrix X is Gaussian, and the response vector y
#' ## is a real-valued quasiconvex function applied to the rows of X.
#' X = mvrnorm(20,rep(0,8), diag(8), tol=1e-06,empirical=FALSE)
#' y = exp(-rowSums(X))
#'
#' out = qconvcheck(X,y, tol = 1e-3)
qconvcheck<-function(X,y, tol = 1e-3){
  X<-as.matrix(X)
  y<-as.matrix(y)
  n<-nrow(X)
  d<-ncol(X)
  if(length(y)!=n)
  {
    print("ERROR in qconvcheck: Dimension of the testing point does not match the dimension of the rows of the base matrix!")
  }
  else
  {
    output.qconv<-"Quasiconvex"
    for(i in 1:n)
    {
      if(sum(y<(y[i]-tol))>0)
      {
        Xchull<-X[(y<(y[i]-tol)),]
        if(sum(y<(y[i]-tol))==1)
        {
          Xchull<-t(as.matrix(Xchull))
        }
        flag<-chullcheck(Xchull,X[i,])
        if(flag==1)
        {
          output.qconv<-"Not Quasiconvex"
          break
        }
      }
    }
    return(output.qconv)
  }
}





#'Checks whether the point test is in the convex hull of the rows of X. 
#'This function will take the matrix X and the vector test as inputs,
#'and return a binary variable, which is 1 if the point test is in the
#'convex hull of the rows of X, and 0 otherwise.
#'
#' @param X              An n by d matrix.
#' @param test           An n by 1 vector. 
#' 
#' @return               A binary variable, which is 1 if the point test is in the
#'                       convex hull of the rows of X, and 0 otherwise.    
#' @examples
#' 
#' ## We demonstrate an example, where the point test is in the convex hull of the rows of X.
#' X<-matrix(c(0,0,1,1,0,1,0,1),4,2)
#' test<-c(0.3,0.4)
#'
#' out = chullcheck(X,test)
chullcheck<-function(X,test){
  X<-as.matrix(X)
  test<-as.matrix(test)
  n<-nrow(X)
  d<-ncol(X)
  if(length(test)!=d)
  {
    print("ERROR in chullcheck: Dimension of the testing point does not match the dimension of the rows of the base matrix!")
  }
  else
  {
    objective.in<-rep(0,n)
    const.mat<-rbind(t(X),rep(1,n),diag(n))
    const.dir<-c(rep("=",(d+1)),rep(">=",n))
    const.rhs<-rbind(test,1,as.matrix(rep(0,n)))
    out<-lp(direction="min", objective.in,const.mat,const.dir,const.rhs)
    if(out$status==2)
    {
      return(0)
    }
    else
    {
      return(1)
    }
  }
}





#'Checks whether the point test is in the upper orthant of the convex 
#'hull of the rows of X. This function will take the matrix X and the
#'vector test as inputs, and return a binary variable, which is 1 if 
#'the point test is in the upper orthant of the convex hull of the 
#'rows of X, and 0 otherwise.
#'
#' @param X              An n by d matrix.
#' @param test           An n by 1 vector.  
#' 
#' @return               A binary variable, which is 1 if the point test 
#'                       is in the upper orthant of the convex hull of the 
#'                       rows of X, and 0 otherwise.
#' @examples
#' 
#' ## We demonstrate an example, where the point test is in the upper orthant of the convex hull 
#' ## of the rows of X.
#' X<-matrix(c(0,0,1,1,0,1,0,1),4,2)
#' test<-c(2,7)
#'
#' out = chullupcheck(X,test)
chullupcheck<-function(X,test){
  X<-as.matrix(X)
  test<-as.matrix(test)
  n<-nrow(X)
  d<-ncol(X)
  if(length(test)!=d){
    print("ERROR in chullcheck: Dimension of the testing point does not match the dimension of the rows of the base matrix!")
  } else{
    objective.in<-rep(0,(n+d))
    const.mat1<-cbind(t(X),diag(d))
    const.mat2<-c(rep(1,n),rep(0,d))
    const.mat<-rbind(const.mat1,const.mat2,diag(n+d))
    const.dir<-c(rep("=",(d+1)),rep(">=",(n+d)))
    const.rhs<-rbind(test,1,as.matrix(rep(0,(n+d))))
    out<-lp(direction="min", objective.in,const.mat,const.dir,const.rhs)
    if(out$status==2){
      return(0)
    } else{
      return(1)
    }
  }
}





#'Checks whether the point test is in the lower orthant of the convex 
#'hull of the rows of X. This function will take the matrix X and the
#'vector test as inputs, and return a binary variable, which is 1 if 
#'the point test is in the lower orthant of the convex hull of the 
#'rows of X, and 0 otherwise.
#'
#' @param X              An n by d matrix.
#' @param test           An n by 1 vector. 
#' 
#' @return               A binary variable, which is 1 if the point test 
#'                       is in the lower orthant of the convex hull of the 
#'                       rows of X, and 0 otherwise.
#' @examples
#' 
#' ## We demonstrate an example, where the point test is in the lower orthant of the convex hull 
#' ## of the rows of X.
#' X<-matrix(c(0,0,1,1,0,1,0,1),4,2)
#' test<-c(-1,-2)
#'
#' out = chulldowncheck(X,test)

chulldowncheck<-function(X,test){
  X<-as.matrix(X)
  test<-as.matrix(test)
  n<-nrow(X)
  d<-ncol(X)
  if(length(test)!=d){
    print("ERROR in chullcheck: Dimension of the testing point does not match the dimension of the rows of the base matrix!")
  } else{
    objective.in<-rep(0,(n+d))
    const.mat1<-cbind(t(X),-diag(d))
    const.mat2<-c(rep(1,n),rep(0,d))
    const.mat<-rbind(const.mat1,const.mat2,diag(n+d))
    const.dir<-c(rep("=",(d+1)),rep(">=",(n+d)))
    const.rhs<-rbind(test,1,as.matrix(rep(0,(n+d))))
    out<-lp(direction="min", objective.in,const.mat,const.dir,const.rhs)
    if(out$status==2){
      return(0)
    } else{
      return(1)
    }
  }
}





#'Assumes that cplex.fhat is a decreasing/increasing/simple quasiconvex,
#'fit obtained from the matrix X of regressors and some response vector.
#'This function will take as inputs, the matrix X, the fitted vector cplex.fhat, 
#'a matrix nondatapoint_matrix whose rows may not be among those of X, and a 
#'categorical variable Monotone which denotes the type of regression in terms of 
#'monotonicity. It will return a quasiconvex and monotone (type of monotonicity
#'determined by the user entered option for the variable Monotone) function
#'evaluated at the rows of the matrix nondatapoint_matrix, which is an 
#'extrapolation of the fitted values in the vector cplex.fhat.
#'
#' @param X                         An n by d matrix of regressors.
#' @param cplex.fhat                An n by 1 vector of fitted values.
#' @param nondatapoint_matrix       An N by d matrix whose rows may not be among the rows of X.
#' @param Monotone                  A categorical variable denoting the type of regression in terms of monotonicity.
#'                                  The user must input one of the following three types:
#'                                  "non.inc": means nonincreasing, quasiconvex regression,
#'                                  "non.dec": means nondecreasing, quasiconvex regression,
#'                                  "non": means quasiconvex regression.
#'                                  
#' @return                          A vector, representing a quasiconvex and monotone (type of monotonicity
#'                                  determined by the user entered option for the input variable Monotone) 
#'                                  function evaluated at the rows of the matrix nondatapoint_matrix, 
#'                                  which is an extrapolation of the fitted values in the vector cplex.fhat.          
#' @examples
#' library(MASS)
#' X = mvrnorm(20,rep(0,8), diag(8), tol=1e-06,empirical=FALSE)
#' y = exp(-rowSums(X)) + rnorm(20)
#'
#' ret = Quasiconvex.cplex(X, y,   Monotone = "non.inc", Max.b = 10^4, MM =10^4, ep =0.001,  time.limit= 200, tol =  1e-06,  output.flag = 1)
#' Xnondata = mvrnorm(60,rep(0,8), diag(8), tol=1e-06,empirical=FALSE)
#' X = rbind(X,Xnondata)
#' out=interpolation.matrix.func(X,ret$f.hat,Xnondata,Monotone="non.inc")

interpolation.matrix.func<-function(X,cplex.fhat,nondatapoint_matrix,Monotone = c("non.inc", "non.dec","none"))
{
  rr<-nrow(nondatapoint_matrix)
  out<-rep(0,rr)
  if( Monotone == "non.inc")
  {
    for(i in 1:rr){
      out[i]<-interpol.func.Dec(X,cplex.fhat,nondatapoint_matrix[i,])
    }
  }
  else if( Monotone == "non.dec")
  {
    for(i in 1:rr){
      out[i]<-interpol.func.Inc(X,cplex.fhat,nondatapoint_matrix[i,])
    }
  }
  else if( Monotone == "none")
  {
    for(i in 1:rr){
      out[i]<-interpol.func(X,cplex.fhat,nondatapoint_matrix[i,])
    }
  }
  return(out)
}