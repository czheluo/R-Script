#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
    'input','i',1,'character',
	'output','o',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
contact:meng.luo@majorbio.com or czheluo@gmail.com
Usage example:
Usage:
    --input input file name
	--output output file name
	--help		usage
\n")
	q(status=1);
}
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec)}
times<-Sys.time()
library(R.matlab)
library(glmnet)
library(rrBLUP)
library(BGLR)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

mat<-readMat('PHS_all.mat')
Y<-read.csv('phs_sim0.5.csv',header = TRUE)

mn<-21
acc1<-matrix(0,10,mn)
acc2<-matrix(0,10,mn)
ytest<-matrix(0,185,mn)
nIter=3000; burnIn=1000
nfolds<-10
time1 <- Sys.time()
for (ij in 2:22) {
  luoy<-Y[,c(1,ij)]
  
  #for (j in 1:mn) {
    
    #luoy<-luoy[sample(nrow(luoy)),]
    folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    for (i in 1:nfolds) {
      #X<-luoGD[,-1]
      M<-as.matrix(mat$x)
      #M<-as.matrix(X)
      testIndexes <- which(folds==i)
      yte <- luoy[testIndexes,2]
      ytr <- luoy[-testIndexes,2]
      #Y<-yte[,2]
      T<-as.matrix(M[-testIndexes,])
      t<-as.matrix(M[testIndexes,])
      fmBB=BGLR(y=ytr,ETA=list( list(X=T,model='BayesB')), 
                nIter=nIter,burnIn=burnIn,saveAt='bb')
      acc1[i,ij-1]<-cor(t%*%fmBB$ETA[[1]]$b+fmBB$mu,yte)
      acc2[i,ij-1]<-cor(T%*%fmBB$ETA[[1]]$b+fmBB$mu,ytr)
      ytest[testIndexes,ij-1]<-t%*%fmBB$ETA[[1]]$b+fmBB$mu
    }
 # }
  
}
writeMat(paste('BB/phssim/y_sim05_PHS_BB.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
t1 = f(time1)
t1


mn<-21
acc1<-matrix(0,10,mn)
acc2<-matrix(0,10,mn)
ytest<-matrix(0,185,mn)
nIter=3000; burnIn=1000
nfolds<-10
time1 <- Sys.time()
for (ij in 2:22) {
  luoy<-Y[,c(1,ij)]
  
  #for (j in 1:mn) {
  
  #luoy<-luoy[sample(nrow(luoy)),]
  folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
  #ans1<-NULL
  for (i in 1:nfolds) {
    #X<-luoGD[,-1]
    M<-as.matrix(mat$x)
    #M<-as.matrix(X)
    testIndexes <- which(folds==i)
    yte <- luoy[testIndexes,2]
    ytr <- luoy[-testIndexes,2]
    #Y<-yte[,2]
    T<-as.matrix(M[-testIndexes,])
    t<-as.matrix(M[testIndexes,])
    fmBA=BGLR(y=ytr,ETA=list( list(X=T,model='BayesA')), 
              nIter=nIter,burnIn=burnIn,saveAt='ba')
    acc1[i,ij-1]<-cor(t%*%fmBA$ETA[[1]]$b+fmBA$mu,yte)
    acc2[i,ij-1]<-cor(T%*%fmBA$ETA[[1]]$b+fmBA$mu,ytr)
    ytest[testIndexes,ij-1]<-t%*%fmBA$ETA[[1]]$b+fmBA$mu
  }
  # }
  
}
writeMat(paste('BA/phssim/y_sim05_PHS_BA.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
t1 = f(time1)
t1


mn<-21
acc1<-matrix(0,10,mn)
acc2<-matrix(0,10,mn)
ytest<-matrix(0,185,mn)
nIter=3000; burnIn=1000
nfolds<-10
time1 <- Sys.time()
for (ij in 2:22) {
  luoy<-Y[,c(1,ij)]
  
  #for (j in 1:mn) {
  
  #luoy<-luoy[sample(nrow(luoy)),]
  folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
  #ans1<-NULL
  for (i in 1:nfolds) {
    #X<-luoGD[,-1]
    M<-as.matrix(mat$x)
    #M<-as.matrix(X)
    testIndexes <- which(folds==i)
    yte <- luoy[testIndexes,2]
    ytr <- luoy[-testIndexes,2]
    #Y<-yte[,2]
    T<-as.matrix(M[-testIndexes,])
    t<-as.matrix(M[testIndexes,])
    fmBL=BGLR(y=ytr,ETA=list( list(X=T,model='BL')), 
              nIter=nIter,burnIn=burnIn,saveAt='bl')
    acc1[i,ij-1]<-cor(t%*%fmBL$ETA[[1]]$b+fmBL$mu,yte)
    acc2[i,ij-1]<-cor(T%*%fmBL$ETA[[1]]$b+fmBL$mu,ytr)
    ytest[testIndexes,ij-1]<-t%*%fmBL$ETA[[1]]$b+fmBL$mu
  }
  # }
  
}
writeMat(paste('BL/phssim/y_sim05_PHS_BL.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
t1 = f(time1)
t1

mn<-21
acc1<-matrix(0,10,mn)
acc2<-matrix(0,10,mn)
ytest<-matrix(0,185,mn)
nIter=3000; burnIn=1000
nfolds<-10
time1 <- Sys.time()
for (ij in 2:22) {
  luoy<-Y[,c(1,ij)]
  
  #for (j in 1:mn) {
  
  #luoy<-luoy[sample(nrow(luoy)),]
  folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
  #ans1<-NULL
  for (i in 1:nfolds) {
    #X<-luoGD[,-1]
    M<-as.matrix(mat$x)
    #M<-as.matrix(X)
    testIndexes <- which(folds==i)
    yte <- luoy[testIndexes,2]
    ytr <- luoy[-testIndexes,2]
    #Y<-yte[,2]
    T<-as.matrix(M[-testIndexes,])
    t<-as.matrix(M[testIndexes,])
    fmBC=BGLR(y=ytr,ETA=list( list(X=T,model='BayesC')), 
              nIter=nIter,burnIn=burnIn,saveAt='bc')
    acc1[i,ij-1]<-cor(t%*%fmBC$ETA[[1]]$b+fmBC$mu,yte)
    acc2[i,ij-1]<-cor(T%*%fmBC$ETA[[1]]$b+fmBC$mu,ytr)
    ytest[testIndexes,ij-1]<-t%*%fmBC$ETA[[1]]$b+fmBC$mu
  }
  # }
  
}
writeMat(paste('BC/phssim/y_sim05_PHS_BC.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
t1 = f(time1)
t1

library(R.matlab)
library(glmnet)
library(rrBLUP)
library(BGLR)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

mat<-readMat('PHS_all.mat')
Y<-read.csv('phs_sim0.5.csv',header = TRUE)

mn<-20
acc1<-matrix(0,10,mn)
acc2<-matrix(0,10,mn)
ytest<-matrix(0,185,mn)
nIter=3000; burnIn=1000
nfolds<-10
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]
  
  #for (j in 1:mn) {
    #luoy<-luoy[sample(nrow(luoy)),]
    folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    for (i in 1:nfolds) {
      #X<-luoGD[,-1]
      M<-as.matrix(mat$x)
      #M<-as.matrix(X)
      testIndexes <- which(folds==i)
      yte <- luoy[testIndexes,2]
      ytr <- luoy[-testIndexes,2]
      #Y<-yte[,2]
      T<-as.matrix(M[-testIndexes,])
      t<-as.matrix(M[testIndexes,])
      RR<- mixed.solve(y=ytr,Z=T) 
      acc1[i,ij-1]<-cor(t%*%RR$u,yte)
      acc2[i,ij-1]<-cor(T%*%RR$u,ytr)
      ytest[testIndexes,ij-1]<-t%*%RR$u
    }
  #}
    writeMat(paste('RR/phssim/y_sim05_PHS_RR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)

}
t1 = f(time1)
t1
escaptime<-Sys.time()-times;
print("Done!");
print(escaptime)
