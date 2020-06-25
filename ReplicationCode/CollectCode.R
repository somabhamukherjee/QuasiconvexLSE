 rm(list=ls())

rep.no <- 200
cases <- c("PH","NPH","NPNH")
varainceS <- c(1,2,3,4)
col.mat<- NULL
col.mat$Example <- NULL
col.mat$Variance <- NULL
col.mat$Replication <- NULL
col.mat$grid.loss <- NULL
col.mat$L2Pn.loss <- NULL
col.mat$Estim.name <- NULL
# nn <- 50
# grid.size <- 800

nn <- 100
grid.size <- 800
ind<- 1
caseno <- 1
vv <- 1
count<-0
count0<- 0
for (ind in 1:rep.no){
  if(ind%%10==1) print(ind)
  for (caseno in 1:length(cases)){
    for(vv  in 1:length(varainceS)){
      case <- cases[caseno]
      Save.file.name<-paste(getwd(),'/data/Oper_nn_', nn,"case", case, "var_", varainceS[vv], "_ind_",ind,".RData",sep='')
        # Save.file.name<-paste(getwd(),'/data/Oper_nn_', nn,"grid", grid.size, "case", case, "var_", varainceS[vv], "_ind_",ind,".RData",sep='')
      if(file.exists(Save.file.name)){
        data.mat <- NULL
        LSE <- NULL
        bw.fit <- NULL
        X.grid <- NULL
        bw.fit <- NULL
        Estimates.data <- NULL
        Estimates.grid <- NULL
        load(Save.file.name)
        if(nrow(Estimates.grid)==(floor(sqrt(grid.size)))^2){
          for(jj in 1:ncol(Estimates.grid)){
            col.mat$Example <- c(col.mat$Example,case)
            col.mat$Variance <- c(col.mat$Variance, varainceS[vv])
            col.mat$Replication <- c(col.mat$Replication,ind)
            col.mat$grid.loss <- c(col.mat$grid.loss, sqrt(mean((Estimates.grid[,jj]- Estimates.grid[,1])^2))) 
            col.mat$L2Pn.loss <- c(col.mat$L2Pn.loss, sqrt(mean((Estimates.data[,jj]- Estimates.data[,1])^2)))
            col.mat$Estim.name <- c(col.mat$Estim, colnames(Estimates.grid)[jj])
          }
        } else{
          print(paste("New grids need to be done", case, vv, ind))
          count<- count+1
        }
      } else{
        print(paste("File doesnt exist", case, vv, ind))
        count0<- count0+1
      }
    }
  }
}

fff<- as.data.frame(col.mat)
varainceSs <- c(1,2,3,4)
for (cases in c("PH","NPH","NPNH")){
  for(vv in varainceSs){
    print(paste(cases, "vv", vv, "count",sum((fff$Example ==cases)*( fff$Estim.name=="f0")*(fff$Variance ==vv))))
  }
}

save(col.mat, file=paste0("Oper_",nn,"Grid",grid.size,"Collected_data".RData"))