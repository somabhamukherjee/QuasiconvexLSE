rm(list=ls())
#Use the file created by "CollectCode.R". File name can be found at line 68 of CollectCode.R
load("Opern100Grid800Collected_data.RData")
library(reshape2)
library(ggplot2)
ls()
fff<- as.data.frame(col.mat)
colnames(fff)[6] <- "Estimators"
head(fff)
levels(fff[,6])[c(3,4)]<- c("Kernel", "ChenEtAl")
#seeing counts  
varainceSs <- c(1,2,3,4)
for (cases in c("PH","NPH","NPNH")){
  for(vv in varainceSs){
    print(paste(cases, "vv", vv, "count",sum((fff$Example ==cases)*( fff$Estimators=="f0")*(fff$Variance ==vv))))
  }
}
nn<- 100
grid.size<- 800

fff<- fff[fff$Estimators!="f0",]
str(fff)
levels(fff[,1]) <- c("Setting II", "Setting III", "Setting I")
head(fff)
tail(fff)

fff$Example_f = factor(fff$Example, levels=c("Setting I", "Setting II", "Setting III"))


pdf(paste0("L2pn_", nn,"grid", grid.size,Sys.Date(),".pdf"), width =10, height=5)
 ggplot(data= fff, mapping= aes(x=Estimators, y =L2Pn.loss, group = Estimators, color= Estimators)) +
   geom_boxplot() +
   facet_grid(Example_f~as.factor(Variance)) +
   theme(legend.position = "none",strip.text.y = element_text(size=12, face="bold"),axis.text.x = element_text(face="bold",size=10))+ labs(y=expression(paste("In sample ",L[2] ,"-loss")), x = "")
 dev.off()

pdf(paste0("Gridloss_", nn,"grid", grid.size,Sys.Date(),".pdf"), width =10, height=7)
 ggplot(data= fff, mapping= aes(x=Estimators, y =grid.loss, group = Estimators, color= Estimators))+  geom_boxplot() + facet_grid(Example~as.factor(Variance)) +theme(strip.text.y = element_text(size=12, face="bold"))
dev.off()
