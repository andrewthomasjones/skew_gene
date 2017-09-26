
fmcfust2 <- function(Y, zeta=1, emp,g=2,  initial=NULL, known=NULL, clust=NULL, itmax=100, eps=1e-6, nkmeans=20, verbose=T, method=c("moments","transformation","EMMIXskew","EMMIXuskew"), convergence=c("Aitken","likelihood","parameters")) {
  
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
  
  if(g==2){
    initial_1 <-fmcfust(g, Y, initial =list(pro=pi_0(Y,zeta)), itmax=10, verbose=F)
    initial_1$pro<-pi_0(Y,zeta)
  }else{
    initial_1 <-fmcfust(g, Y, itmax=10, verbose=F)
    
  }
  
  ndelta<-4
  res <- fmfustt2((Y), initial_1, ndelta, emp, itmax=itmax)
  return(res)
}




pi_0<-function(data, zeta){
  p<-sum(data<zeta)/(length(data)*pnorm(zeta))
  return(c(p,1-p))
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

plot_normal_densities<-function(fit, data, k, bins,name){
  data1<-as.data.frame(normal_densities(fit,k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}


skew_densities<-function(fit,data, k){
  G<-length(fit$pro)
  y<-array(NA, c(k,G+1))
  xseq<-seq(min(data), max(data), length.out=k)
  y[,1]<-xseq
  for(i in 2:(G+1)){
    y[,i]<-fit$pro[i-1]*dcfust(dat=as.matrix(xseq,nrow=1), mu=as.vector(fit$mu[[i-1]]), sigma=fit$sigma[[i-1]], delta=fit$delta[[i-1]], dof=fit$dof[i-1])
  }
  
  colnames(y)<-c("x", paste0('y', seq(1:G)))
  return(y)
}


plot_skew_densities<-function(fit, data, k, bins, name){
  data1<-as.data.frame(skew_densities(fit,data, k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(data=data1, aes(x=x,  y=y1+y2))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}

skew_norm_densities<-function(fit,data, k){
  G<-length(fit$pro)
  y<-array(NA, c(k,G+1))
  xseq<-seq(min(data), max(data), length.out=k)
  y[,1]<-xseq
  y[,2]<-dnorm(y[,1], mean=fit$mu[[1]], sd=sqrt(fit$sigma[[1]]))*fit$pro[1]
  for(i in 3:(G+1)){
    y[,i]<-fit$pro[i-1]*dcfust(dat=as.matrix(xseq,nrow=1), mu=as.vector(fit$mu[[i-1]]), sigma=fit$sigma[[i-1]], delta=fit$delta[[i-1]], dof=fit$dof[i-1])
  }
  
  colnames(y)<-c("x", paste0('y', seq(1:G)))
  return(y)
}

plot_skew_norm_densities<-function(fit, data, k, bins, name){
  data1<-as.data.frame(skew_norm_densities(fit,data, k))
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

FDR_table<-function(fit, cList){
  n<-length(cList)
  table1<-data.frame(c0=array(NA, n), N_r=array(NA, n), FDR=array(NA, n), FNDR=array(NA, n), FNR=array(NA, n),FPR=array(NA, n))
  
  for(i in 1:n){
    table1[i,]<-round(c(cList[i],unlist(FDR_etc(fit, cList[i]))),3)
  }
  
  return(table1)
}

#
#  EM algorithm for Mixture of Multivariate Canonical Fundamental Skew t-distributioins
#  Package: EMMIX-cskew
#  Version: 0.9-1
#
#  Code by S. X. Lee
#  Updated on 07 Sep 2015
#
# Lee S.X. and Mclachlan, G.J. (2015) Finite mixtures of canonical fundamental 
# skew t-distributions: the unification of the restricted and unrestricted 
# skew t-mixture models. Statistics and Computing. doi:10.1007/s11222-015-9545-x
#
# The following code are adapted from other packages
#
################################################################################
#  misc.r
#                           Miscellaneous Tools
#
################################################################################

#skewness
skewness <- function (x, na.rm = FALSE){
  if (is.matrix(x))
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm) x <- x[!is.na(x)]
    n <- length(x)
    (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x))
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

#mahalanobis
mahalanobis. <- function(x, center, cov, inverted=FALSE, ...) {
  x <- if(is.vector(x)) matrix(x, ncol=length(x)) else as.matrix(x)
  x <- t(sweep(x, 2, center))# = (x - center)
  retval <- colSums(x * if(inverted) cov%*%x else solve(cov, x, ...))
  names(retval) <- rownames(x)
  retval
}

mahalanobis <- function(x, center, cov, inverted=FALSE, ...){
  x <- if(is.vector(x)) matrix(x, ncol=length(x)) else as.matrix(x)
  x <- sweep(x, 2, center)# = (x - center)
  if(!inverted)
    cov <- solve(cov, ...)
  retval <- rowSums((x%*%cov) * x)
  names(retval) <- rownames(x)
  retval
}

#
isPosDef <- function(M) { 
  if ( all(M == t(M) ) ) { 
    if (  all(eigen(M)$values > 0) ) {TRUE}
    else {FALSE} 
  }else {FALSE}  
} 

#
is.whole <- function(a) { 
  (is.numeric(a) && floor(a)==a) ||
    (is.complex(a) && floor(Re(a)) == Re(a) && floor(Im(a)) == Im(a))
}


#multivariate normal density
dmn <- function(x, mu=rep(0,d), Sigma, log=FALSE) {
  d  <- if(is.matrix(Sigma)) ncol(Sigma) else 1
  if(d>1 & is.vector(x)) x <- matrix(x, 1, d)
  n  <- if(d==1)  length(x) else nrow(x) 
  X  <- t(matrix(x, nrow=n, ncol=d)) - mu
  Q  <- apply((solve(Sigma)%*% X)* X, 2, sum) 
  logDet <- sum(logb(abs(diag(qr(Sigma)[[1]]))))
  logPDF <- as.vector(Q + d*logb(2*pi)+logDet)/(-2)
  if(log) logPDF else exp(logPDF)
}

#multivariate t density
dmt <- function(dat, mu, sigma, dof = Inf, log = FALSE) {
  if (dof == Inf)  return(dmn(dat, mu, sigma, log = log))
  dof <- round(dof)
  d  <- if(is.matrix(sigma)) ncol(sigma) else 1
  n  <- if(d==1) length(dat) else nrow(dat)
  X <- t(matrix(dat, nrow = n, ncol = d)) - mu
  Q <- apply((solve(sigma) %*% X) * X, 2, sum)
  logDet <- sum(logb(abs(diag(qr(sigma)$qr))))
  logPDF <- (lgamma((dof + d)/2) - 0.5 * (d * logb(pi * dof) + logDet)
             - lgamma(dof/2) - 0.5 * (dof + d) * logb(1 + Q/dof))
  if(log) logPDF else exp(logPDF)
}


#error rate
error.rate <- function(clust1,clust2){
  if((n=length(clust1))!=length(clust2))  stop("error: length not equal")
  if( (g=length(table(clust1)))!=length(table(clust2))) stop("the number of clusters are not equal")
  permute<-function(a){
    n<-length(a)
    if(n==1) f<-a
    else{
      nm<-gamma(n)
      f<-array(0,c(n,n*nm))
      j<-1
      
      for(i in a){
        f[1, (j-1)*nm+1:nm]<-i
        f[-1,(j-1)*nm+1:nm]<-permute(setdiff(a,i))
        j<-j+1
      }
    }    
    f
  }
  
  id<-1:n
  cmb<-permute(1:g)
  nperm<-ncol(cmb)
  rate<-rep(0,nperm)
  
  for(i in 1:nperm){
    tmp<-rep(0,g)
    tc<-rep(0,n)
    for(j in 1:g)
      tc[clust2==j]=cmb[j,i]
    
    for(j in 1:g){  
      tmp1<-0 
      for(k in (1:g)[-j])
        tmp1<-tmp1+length(intersect(id[clust1==j],id[tc==k]))
      tmp[j]<-tmp1
    }
    rate[i]<-sum(tmp)/n
  }
  min(rate)
}




plot_skew_densities3<-function(fit, data, k, bins, name){
  data1<-as.data.frame(skew_densities(fit,data, k))
  data2<-data.frame(data=data)
  p=ggplot(data=data1)+geom_histogram(data=data2, aes(y=..density..,x=data), fill='white', colour='grey', bins=bins)+geom_line(aes(x=x,  y=y1))+geom_line(aes(x=x,  y=y2))+geom_line(aes(x=x,  y=y3))+geom_line(data=data1, aes(x=x,  y=y1+y2+y3))
  p=p+scale_colour_discrete(guide=FALSE)+theme_bw()+xlab("z-scores")+ylab("Density")+ggtitle(name)
  print(p)
}




