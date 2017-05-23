library(EMMIXcskew)
library(EMMIXskew)

library(mnormt)
library(locfdr)
library(qvalue)
library(cepp)
source("fmcfust.R")
source("fmmst.R")


data("hedenfalk")
hedenfalk$z<-qnorm(1-hedenfalk$p)
data("Colon")

for(i in 1:2000){
  Colon$t[i]<-t.test( Colon$X[Colon$Y==1,i] , Colon$X[Colon$Y==2,i] )$statistic
}
Colon$p<-pt(Colon$t,60)
Colon$z<-qnorm(1-Colon$p)

data("hivdata")# already pval


x<-fmcfust2((hedenfalk$z),1, verbose=T, emp=F)
x2<-fmcfust2((hedenfalk$z),1, verbose=T, emp=T)

y<-fmcfust2((hivdata),3, verbose=T, emp=T)
z<-fmcfust2(Colon$z,1, verbose=T, emp=T, itmax=20)

heden<-namer(x)
log(length(hed_fit_norm$x))*3-2*hed_fit_norm$loglik
FDR_etc(heden, 0.1)

plot_skew_densities(z, Colon$z, 1000, 50, "Colon Skew")
plot_skew_densities(y, hivdata, 1000, 50, "HIV Skew")
plot_skew_densities(x, hedenfalk$z, 1000, 50, "Hedenfalk Skew")


