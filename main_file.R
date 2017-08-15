library(EMMIXcskew)
library(EMMIXskew)
library(mnormt)
library(locfdr)
library(ggplot2)
# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")
library(qvalue)
library(cepp)
library(mixtools)

source("fmcfust.R")
source("fmmst.R")
source("util.R")

###############################################################
data("hedenfalk")
hedenfalk$z<-qnorm(1-hedenfalk$p)
data("Colon")

for(i in 1:2000){
  Colon$t[i]<-t.test( Colon$X[Colon$Y==1,i] , Colon$X[Colon$Y==2,i] )$statistic
}
Colon$p<-pt(Colon$t,60)
Colon$z<-qnorm(1-Colon$p)

data("hivdata")# already pval

###############################################################


#normal fits
hiv_fit_norm <-normalmixEM(x=hivdata, lambda =pi_0(hivdata,3), k=2)
hed_fit_norm <- normalmixEM(x=hedenfalk$z,k=2, lambda = pi_0(hedenfalk$z,0), mean.constr = c(0, NA), sd.constr=c(1, NA))
col_fit_norm <- normalmixEM(x=Colon$z, lambda =pi_0(Colon$z,0), k=2, mean.constr = c(0, NA), sd.constr=c(1, NA))

plot_normal_densities(hiv_fit_norm,hivdata, 1000, 40, "HIV Normal")
plot_normal_densities(hed_fit_norm,hedenfalk$z, 1000, 40, "Hedenfalk Normal")
plot_normal_densities(col_fit_norm, Colon$z, 1000, 50, "Colon Normal")

FDR_etc(namer(hiv_fit_norm),0.1)
FDR_etc(namer(hed_fit_norm),0.1)
FDR_etc(namer(col_fit_norm),0.1)


#free fits
hiv_fit <- fmcfust(2, hivdata, initial =list(pro=pi_0(hivdata,3)))
hed_fit <- fmcfust(2, hedenfalk$z,initial =list(pro=pi_0(hedenfalk$z,0)))
col_fit <- fmcfust(2, Colon$z,initial =list(pro=pi_0(Colon$z, 0)))

plot_skew_densities(col_fit , Colon$z, 1000, 50, "Colon Skew")
plot_skew_densities(hiv_fit , hivdata, 1000, 50, "HIV Skew")
plot_skew_densities(hed_fit , hedenfalk$z, 1000, 50, "Hedenfalk Skew")

FDR_etc(namer(hiv_fit),0.1)
FDR_etc(namer(hed_fit),0.1)
FDR_etc(namer(col_fit),0.1)

#partly free fits
#emprical null - bad actually

x<-c(rnorm(1000), rnorm(200, mean = 3))
x1<-fmcfust2(x,1.5, verbose=T, emp=T)
plot_skew_norm_densities(x1 , x, 1000, 50, "Test")

# hed_fit_new2<-fmcfust2((hedenfalk$z),1, verbose=T, emp=T)
# plot_skew_densities(hed_fit_new2 , hedenfalk$z, 1000, 50, "Hedenfalk Skew/Norm")
# FDR_etc(namer(hed_fit_new2),0.1)
#         

hiv_fit_new2<-fmcfust2((hivdata),3, verbose=T, emp=T)
plot_skew_norm_densities(hiv_fit_new2 , hivdata, 1000, 50, "HIV Skew/Norm")
FDR_etc(namer(hiv_fit_new2),0.1)

# col_fit_new2<-fmcfust2(Colon$z,1, verbose=T, emp=T)
# plot_skew_densities(col_fit_new2 , Colon$z, 1000, 50, "Colon Skew/Norm")
# FDR_etc(namer(col_fit_new2),0.1)

#partly free fits
#theoretocal null
###############################################################################
hed_fit_new<-fmcfust2(hedenfalk$z,1, verbose=T, emp=F)
plot_skew_norm_densities(hed_fit_new , hedenfalk$z, 1000, 50, "Hedenfalk Skew/Norm")
FDR_etc(namer(hed_fit_new),0.1)

# hiv_fit_new<-fmcfust2(hivdata,3, verbose=T, emp=F, g=3)
# plot_skew_densities3(hiv_fit_new , hivdata, 1000, 50, "HIV Skew/Norm")
# FDR_etc(namer(hiv_fit_new),0.1)

col_fit_new<-fmcfust2(Colon$z,1, verbose=T, emp=F)
plot_skew_norm_densities(col_fit_new , Colon$z, 1000, 50, "Colon Skew/Norm")
FDR_etc(namer(col_fit_new),0.1)




