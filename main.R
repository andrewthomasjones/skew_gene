# source("https://bioconductor.org/biocLite.R")
# biocLite("lemma")

library(qvalue)
library(locfdr)
library(EMMIXcskew)
library(EMMIXskew)
library(ggplot2)
library(cepp)
library(mclust)

data("hedenfalk") #$p
data("hivdata")# already pval
data("Colon")

for(i in 1:2000){
Colon$t[i]<-t.test( Colon$X[Colon$Y==1,i] , Colon$X[Colon$Y==2,i] )$statistic
}
Colon$p<-pt(Colon$t,60)
Colon$z<-qnorm(1-Colon$p)

hedenfalk$z<-qnorm(1-hedenfalk$p)



#free fits
hiv_fit <- fmcfust(2, hivdata)
hed_fit <- fmcfust(2, hedenfalk$z)
col_fit <- fmcfust(2, Colon$z)

#normal fits
hiv_fit_norm <- Mclust(hivdata, G = 2, modelNames = 'V' )
hed_fit_norm <- Mclust( hedenfalk$p, G = 2, modelNames = 'V' )
col_fit_norm <- Mclust(Colon$z, G = 2, modelNames = 'V' )

#skew and normal combined fits
hiv_fit_new <- fmcfust(2, hivdata)
hed_fit_new <- fmcfust(2, hedenfalk$z)
col_fit_new <- fmcfust(2, Colon$z)

