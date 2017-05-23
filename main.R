source("util.R")

library(qvalue)
library(locfdr)
library(EMMIXcskew)
library(EMMIXskew)
library(ggplot2)
library(cepp)
library('mixtools')

data("hedenfalk") #$p
data("hivdata")# already pval
data("Colon")

bins=40
k=1000

for(i in 1:2000){
Colon$t[i]<-t.test( Colon$X[Colon$Y==1,i] , Colon$X[Colon$Y==2,i] )$statistic
}
Colon$p<-pt(Colon$t,60)
Colon$z<-qnorm(1-Colon$p)

hedenfalk$z<-qnorm(1-hedenfalk$p)#isnt exactly same



#free fits

hiv_fit <- fmcfust(2, hivdata, initial =list(pro=pi_0(hivdata,3)), itmax=10)
hed_fit <- fmcfust(2, hedenfalk$z,initial =list(pro=pi_0(hedenfalk$z,0)))
col_fit <- fmcfust(2, Colon$z,initial =list(pro=pi_0(Colon$z, 0)))

#normal fits

hiv_fit_norm <-normalmixEM(x=hivdata, lambda =pi_0(hivdata,3), k=2)
hed_fit_norm <- normalmixEM(x=hedenfalk$z,k=2, lambda = pi_0(hedenfalk$z,0), mean.constr = c(0, NA), sd.constr=c(1, NA))
col_fit_norm <- normalmixEM(x=Colon$z, lambda =pi_0(Colon$z,0), k=2, mean.constr = c(0, NA), sd.constr=c(1, NA))

#skew and normal combined fits
hiv_fit_new <- fmcfust(1, hivdata[1:100])
hed_fit_new <- fmcfust(2, hedenfalk$z)
col_fit_new <- fmcfust(2, Colon$z)

hed_fit_norm <- normalmixEM(x=hedenfalk$z,k=2, lambda = pi_0(hedenfalk$z,0), mean.constr = c(0, NA), sd.constr=c(1, NA))
plot_normal_densities(hed_fit_norm,hedenfalk$z, 1000, 40, "Hedenfalk")

hiv_fit_norm <-normalmixEM(x=hivdata, lambda =pi_0(hivdata,3), k=2)
plot_normal_densities(hiv_fit_norm,hivdata, 1000, 40, "HIV")

col_fit_norm <- normalmixEM(x=Colon$z, lambda =pi_0(Colon$z,0), k=2, mean.constr = c(0, NA), sd.constr=c(1, NA))
plot_normal_densities(col_fit_norm, Colon$z, 1000, 50, "Colon")
FDR_etc(col_fit_norm,.1)
FDR_etc(hed_fit,.1)



plot_normal_densities(hed_fit_norm, 1000, 40)
plot_normal_densities(hiv_fit_norm, 1000, 40)
plot_normal_densities(col_fit_norm, 1000, 40)

plot_skew_densities(hed_fit, hedenfalk$z, 1000,50, "Hedenfalk Skew")
plot_skew_densities(hiv_fit, hivdata, 1000, 50, "HIV Skew")
plot_skew_densities(col_fit, Colon$z, 1000, 50, "Colon Skew")

plot_skew_densities(hed_fit, hedenfalk$z, 1000,40)
plot_skew_densities(hiv_fit, hivdata, 1000, 40)
plot_skew_densities(col_fit, Colon$z, 1000, 40)


#tidy names on other vars too
fit$isNull=factor(fit$clusters)
levels(fit$isNull)=c("NonNull","Null")
#########################
#now we are ready to plot#

data_plot<-data.frame(values=data,isNull= )
p=ggplot(data=data_plot)+geom_histogram(aes(y=..density.., x=values))+geom_density(aes(x=values,colour=isNull, fill=isNull), alpha=0.4)
p+scale_colour_discrete(guide=FALSE)
print(p)

known<-list(list(matrix(0,1,1),matrix(1,1,1)))


