# see Help for Replication.txt for instructions
sessionInfo()
library(Rcplex)
library(lpSolve)
library(Rearrangement)
library(sp)
library(rgeos)
library(np)
library(QuasiconvexLSE)
# source("../../Functions.R")

source("SimulationFunction.R")
source("monotonize.R")
source("qconvexify_2d.R")
nn <- 100
rep.no <- 200
grid.size.s <- 800
cases <- c("PH","NPH","NPNH")
varainceS <- c(1,2,3,4)


# Setting up caseno and ind
ind <- in_ind%%rep.no
caseno <- (in_ind-ind)/rep.no+1

for(vv  in 1:length(varainceS)){
  print(paste("At variance = ", varainceS[vv],"############"))
  case <- cases[caseno]
  Save.file.name<-paste(getwd(),'/data/Oper_nn_', nn,"case", case, "var_", varainceS[vv], "_ind_",ind,".RData",sep='')
  if(!file.exists(Save.file.name)){
    if( caseno == 1){
      data.mat<-parametric.homothetic(nn, varainceS[vv])
    } else if( caseno == 2){
      data.mat<-nonparametric.homothetic(nn, varainceS[vv])
    } else if( caseno == 3){
      data.mat<-nonparametric.nonhomothetic(nn, varainceS[vv])
    }
    
    x.cord<- seq(min(data.mat[,1]), max(data.mat[,1]), length =floor(sqrt(grid.size.s)))
    y.cord<- seq(min(data.mat[,2]), max(data.mat[,2]), length =floor(sqrt(grid.size.s)))
    X.grid <- expand.grid(x=x.cord, y=y.cord)

    LSE <- Quasiconvex.cplex(data.mat[,1:2], -data.mat[,3], MM=10^9, ep=.1, time.limit = 20200, tol =  1e-04)
    if( caseno == 1){
      f0.grid<-parametric.homothetic.func(X.grid)
      f0.data<-parametric.homothetic.func(data.mat[,1:2])
    } else if( caseno == 2){
      f0.grid<-nonparametric.homothetic.func(X.grid)
      f0.data<-nonparametric.homothetic.func(data.mat[,1:2])
    } else if( caseno == 3){
      f0.grid<-nonparametric.nonhomothetic.func(X.grid)
      f0.data<-nonparametric.nonhomothetic.func(data.mat[,1:2])
    }
    #Interpolation of the LSE at the grid points
    Hat.Y.grid<- - interpolation.matrix.func(data.mat[,1:2], LSE$f.hat, as.matrix(X.grid))

    h.vec.RoT<- sqrt(diag(var(data.mat)))[-ncol(data.mat)]* length(data.mat[,1])^(-1/(4+ncol(data.mat))) #Rule of Thumb bandwidth for kernel estimation

    np.fit.RoT<- fitted(npreg(bws= h.vec.RoT, txdat = data.mat[,1:2], tydat = data.mat[,3], bwtype ="fixed",ckertype ="gaussian" ))
    #Non parametric fit at the grid points
    np.fit.grid.RoT<- fitted(npreg(bws= h.vec.RoT, txdat = data.mat[,1:2], tydat = data.mat[,3],exdat=X.grid,  bwtype ="fixed",ckertype ="gaussian" ))
    # Shape restricted "projection"
    Shape.fit.grid.RoT<- -qconvexify_2d(X.grid, -monotonize(X.grid, np.fit.grid.RoT))
    # Shape restricted "projection" at the data points

    Shape.fit.dat.RoT <- - interpolation.matrix.func(X.grid, -Shape.fit.grid.RoT, data.mat[,1:2])

    #All estimates evaluated at the data points
    Estimates.data <- as.data.frame(list(f0= f0.data, LSE = -LSE$f.hat, npfit.RoT = np.fit.RoT, Sfit.Rot= Shape.fit.dat.RoT))
    #All estimates evaluated at the grid points.
    Estimates.grid <- as.data.frame(list(f0= f0.grid, LSE = Hat.Y.grid, npfit.RoT = np.fit.grid.RoT, Sfit.Rot = Shape.fit.grid.RoT))

    save(data.mat,LSE, bw.fit,X.grid,bw.fit,Estimates.data, Estimates.grid,  file=Save.file.name)
  }

}
