library(EMMIXcskew)
library(EMMIXskew)
library(ggplot2)
library(mnormt)
library(locfdr)
library(qvalue)
library(cepp)
library(mixtools)
source("fit.R")
source("util.R")
source("density.R")


data("hedenfalk")
hedenfalk$z<-qnorm(1-hedenfalk$p)
data("Colon")

for(i in 1:2000){
  Colon$t[i]<-t.test( Colon$X[Colon$Y==1,i] , Colon$X[Colon$Y==2,i] )$statistic
}
Colon$p<-pt(Colon$t,60)
Colon$z<-qnorm(1-Colon$p)

data("hivdata")# already pval

x<-fmcfust((hedenfalk$z),zeta=1, emp=T)
y<-fmcfust((hivdata),zeta=3,emp=T)
z<-fmcfust(Colon$z,zeta=1, emp=T)

x2<-fmcfust((hedenfalk$z),1)
y2<-fmcfust((hivdata),3)
z2<-fmcfust(Colon$z,1)

hiv_fit_norm <- normalmixEM(x=hivdata, lambda =pi_0(hivdata,3), k=2)
hed_fit_norm <- normalmixEM(x=hedenfalk$z,k=2, lambda = pi_0(hedenfalk$z,0), mean.constr = c(0, NA), sd.constr=c(1, NA))
col_fit_norm <- normalmixEM(x=Colon$z, lambda =pi_0(Colon$z,0), k=2, mean.constr = c(0, NA), sd.constr=c(1, NA))

plot_normal_densities(hiv_fit_norm,hivdata, 1000, 50, "Hiv Emp. Normal")
plot_normal_densities(col_fit_norm,Colon$z, 1000, 50, "Colon Normal")
plot_normal_densities(hed_fit_norm,hedenfalk$z, 1000, 50, "Hedenfalk Normal")

plot_skew_densities(z, Colon$z, 1000, 50, "Colon Skew  Emp. Null")
plot_skew_densities(y, hivdata, 1000, 50, "HIV Skew Emp. Null")
plot_skew_densities(x, hedenfalk$z, 1000, 50, "Hedenfalk Skew  Emp. Null")

plot_skew_densities(z2, Colon$z, 1000, 50, "Colon Skew")
plot_skew_densities(y2, hivdata, 1000, 50, "HIV Skew")
plot_skew_densities(x2, hedenfalk$z, 1000, 50, "Hedenfalk Skew")


heden<-namer(x)
FDR_etc(heden, 0.1)
FDR_etc(y, 0.1)
FDR_etc(z, 0.1)






