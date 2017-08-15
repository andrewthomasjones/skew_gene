pi_0<-function(data, zeta){
  p<-sum(data<zeta)/(length(data)*pnorm(zeta))
  return(c(p,1-p))
}


fmcfust2 <- function(Y, zeta=1, emp,  initial=NULL, known=NULL, clust=NULL, itmax=100, eps=1e-6, nkmeans=20, verbose=T, method=c("moments","transformation","EMMIXskew","EMMIXuskew"), convergence=c("Aitken","likelihood","parameters")) {
  g=2
  q=1
  Y<-as.matrix(Y)
  p <- ncol(Y); n <- nrow(Y); fulldebug=F; 
  q<-1
  if (itmax > 1000) itmax <- 1000   #do not allow more than 1000 iterations (too much)
  if(n > 5000 && p>=3) {
    #        cat("  NOTE: For large datasets, please consider using EmSkew for faster computation.\n")
    fulldebug = T
  }
  if(n <= 20*g) stop("sample size is too small!")      
  
  if(verbose) {
    cat("Finite Mixture of Multivariate CFUST Distributions\n")
    if(g==1) cat("with 1 component\n")
    else cat("with ", g, "components\n")
    cat("  ----------------------------------------------------\n\n")
  }  
  initial_1 <-fmcfust(2, Y, initial =list(pro=pi_0(Y,zeta)), itmax=10, verbose=F)
  initial_1$pro<-pi_0(Y,zeta)
  ndelta<-4
  res <- fmfustt2( Y, initial_1, ndelta, emp, itmax=itmax)
  return(res)
}

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

plot_normal_densities<-function(fit,data, k, bins,name){
  data1<-as.data.frame(normal_densities(fit,k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}





skew_densities<-function(fit,data, k){
  G<-length(fit$pro)
  y<-array(NA, c(k,G+2))
  xseq<-seq(min(data), max(data), length.out=k)
  y[,1]<-xseq
  y[,2]<-dnorm(y[,1])*fit$pro[1]
  for(i in 3:(G+2)){
    y[,i]<-fit$pro[i-2]*dcfust(dat=as.matrix(xseq,nrow=1), mu=as.vector(fit$mu[[i-2]]), sigma=fit$sigma[[i-2]], delta=fit$delta[[i-2]], dof=fit$dof[i-2])
    
    
  }
  
  colnames(y)<-c("x", "norm", paste0('y', seq(1:G)))
  return(y)
}

plot_skew_densities<-function(fit, data, k, bins, name){
  data1<-as.data.frame(skew_densities(fit,data, k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}



namer<-function(fit){
  #tidy names on other vars too
  fit$isNull=factor(fit$clusters)
  levels(fit$isNull)=c("NonNull","Null")
  return(fit)
}

FDR_etc<-function(fit, c0){
  
  if(!is.null(fit$tau)){
    tau0<-(fit)$tau[1,]
    tau1<-(fit)$tau[2,]
  }
  if(!is.null(fit$posterior)){
    tau0<-(fit)$posterior[,1]
    tau1<-(fit)$posterior[,2]
  }
  
  N_r<-sum(tau0<c0)
  N<-length(tau1)
  FDR<-sum(tau0*(tau0<c0))/N_r
  FNDR<-sum(tau1*(tau0>c0))/(N-N_r)
  FPR<-sum(tau0*(tau0<c0))/sum(tau0)
  FNR<-sum(tau1*(tau0>c0))/sum(tau1)
  
  return(list(N_r=N_r, FDR=FDR, FNDR=FNDR, FNR=FNR,FPR=FPR))
}
