# source("https://bioconductor.org/biocLite.R")
# biocLite("lemma")
pi_0<-function(data, zeta){
  p<-sum(data<zeta)/(length(data)*pnorm(zeta))
  return(c(p,1-p))
}

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

for(i in 1:2000){
Colon$t[i]<-t.test( Colon$X[Colon$Y==1,i] , Colon$X[Colon$Y==2,i] )$statistic
}
Colon$p<-pt(Colon$t,60)
Colon$z<-qnorm(1-Colon$p)

hedenfalk$z<-qnorm(1-hedenfalk$p)#isnt exactly same



#free fits
<<<<<<< HEAD
hiv_fit <- fmcfust(2, hivdata, initial =list(pro=pi_0(hivdata,3)), itmax=10)
=======
hiv_fit <- fmcfust(2, hivdata, initial =list(pro=pi_0(hivdata,3)))
>>>>>>> b5e8d8445499a913be0882a7bb593a4a4dddcf31
hed_fit <- fmcfust(2, hedenfalk$z,initial =list(pro=pi_0(hedenfalk$z,0)))
col_fit <- fmcfust(2, Colon$z,initial =list(pro=pi_0(Colon$z, 0)))

#normal fits

hiv_fit_norm <-normalmixEM(x=hivdata, lambda =pi_0(hivdata,3), k=2)
<<<<<<< HEAD
hed_fit_norm <- normalmixEM(x=hedenfalk$z,k=2, lambda = pi_0(hedenfalk$z,0), mean.constr = c(0, NA), sd.constr=c(1, NA))
=======
hed_fit_norm <- normalmixEM(x=hedenfalk$z,k=2, lambda = pi_0(hedenfalk$z,0) , maxit=2000, maxrestarts=100, mean.constr = c(0, NA), sd.constr=c(1, NA))
>>>>>>> b5e8d8445499a913be0882a7bb593a4a4dddcf31
col_fit_norm <- normalmixEM(x=Colon$z, lambda =pi_0(Colon$z,0), k=2, mean.constr = c(0, NA), sd.constr=c(1, NA))



#skew and normal combined fits
hiv_fit_new <- fmcfust(1, hivdata[1:100])
hed_fit_new <- fmcfust(2, hedenfalk$z)
col_fit_new <- fmcfust(2, Colon$z)

bins=40
k=1000

normal_densities<-function(fit, k){
  fit$G<-length(fit$lambda)
  y<-array(NA, c(k,fit$G+1))
  xseq<-seq(min(fit$x), max(fit$x), length.out=k)
  y[,1]<-xseq
  for(i in 2:(fit$G+1)){
    y[,i]<-fit$lambda[i-1]*dnorm(xseq, mean=fit$mu[i-1], sd=fit$sigma[i-1])
   
  }
  
  colnames(y)<-c("x", paste0('y', seq(1:fit$G)))
  return(y)
}
<<<<<<< HEAD

plot_normal_densities<-function(fit,data, k, bins,name){
  data1<-as.data.frame(normal_densities(fit,k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}

hed_fit_norm <- normalmixEM(x=hedenfalk$z,k=2, lambda = pi_0(hedenfalk$z,0), mean.constr = c(0, NA), sd.constr=c(1, NA))
plot_normal_densities(hed_fit_norm,hedenfalk$z, 1000, 40, "Hedenfalk")

hiv_fit_norm <-normalmixEM(x=hivdata, lambda =pi_0(hivdata,3), k=2)
plot_normal_densities(hiv_fit_norm,hivdata, 1000, 40, "HIV")

col_fit_norm <- normalmixEM(x=Colon$z, lambda =pi_0(Colon$z,0), k=2, mean.constr = c(0, NA), sd.constr=c(1, NA))
plot_normal_densities(col_fit_norm, Colon$z, 1000, 50, "Colon")
FDR_etc(col_fit_norm,.1)
FDR_etc(hed_fit,.1)
=======
plot_normal_densities<-function(fit, k, bins){
  data1<-as.data.frame(normal_densities(fit,k))
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")
  print(p)
}

plot_normal_densities(hed_fit_norm, 1000, 40)
plot_normal_densities(hiv_fit_norm, 1000, 40)
plot_normal_densities(col_fit_norm, 1000, 40)


>>>>>>> b5e8d8445499a913be0882a7bb593a4a4dddcf31

skew_densities<-function(fit,data, k){
  G<-length(fit$pro)
  y<-array(NA, c(k,G+1))
  xseq<-seq(min(data), max(data), length.out=k)
  y[,1]<-xseq
  
  for(i in 2:(G+1)){
<<<<<<< HEAD
    y[,i]<-fit$pro[i-1]*dcfust(dat=as.matrix(xseq,nrow=1), mu=as.vector(fit$mu[[i-1]]), sigma=fit$sigma[[i-1]], delta=fit$delta[[i-1]], dof=fit$dof[i-1])
=======
    y[,i]<-fit$pro[i-1]*dcfust(dat=as.matrix(xseq,nrow=1), mu=as.vector(hed_fit$mu[[i-1]]), sigma=hed_fit$sigma[[i-1]], delta=hed_fit$delta[[i-1]], dof=hed_fit$dof[i-1])
>>>>>>> b5e8d8445499a913be0882a7bb593a4a4dddcf31

    
  }
  
  colnames(y)<-c("x", paste0('y', seq(1:G)))
  return(y)
}

<<<<<<< HEAD
plot_skew_densities<-function(fit, data, k, bins, name){
  data1<-as.data.frame(skew_densities(fit,data, k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
=======
plot_skew_densities<-function(fit, data, k, bins){
  data1<-as.data.frame(skew_densities(fit,data, k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")
>>>>>>> b5e8d8445499a913be0882a7bb593a4a4dddcf31
  print(p)
}


<<<<<<< HEAD
plot_skew_densities(hed_fit, hedenfalk$z, 1000,50, "Hedenfalk Skew")
plot_skew_densities(hiv_fit, hivdata, 1000, 50, "HIV Skew")
plot_skew_densities(col_fit, Colon$z, 1000, 50, "Colon Skew")
=======
plot_skew_densities(hed_fit, hedenfalk$z, 1000,40)
plot_skew_densities(hiv_fit, hivdata, 1000, 40)
plot_skew_densities(col_fit, Colon$z, 1000, 40)
>>>>>>> b5e8d8445499a913be0882a7bb593a4a4dddcf31


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

pi_0<-function(data, zeta){
  p<-sum(data<zeta)/(length(data)*pnorm(zeta))
  return(c(1-p,p))
}
