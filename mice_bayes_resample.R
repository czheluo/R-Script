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

setwd('D:\\MATLAB\\isreg\\rice2017\\mouse')

mat<-readMat('cfw_R.mat')
mx=as.matrix(mat$x);rm(mat)
Y<-read.table('CFW_ALL_MLMM.txt',header = T)
mn<-20;
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:14) {
  luoy<-Y[,c(1,ij)]
  name.of.trait=noquote(names(luoy)[2])
  nn<-which(is.na(luoy[,2]));luoy<-luoy[-nn,]
  M=mx[-nn,];n=1161-length(nn);
  testing <- sample(n,round(n/5),replace=F)
  ytest<-matrix(0,length(testing),mn);test<-matrix(0,length(testing),mn)
  for (i in 1:20) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    #M<-as.matrix(mat$x)
    #n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBA=BGLR(y=ytr,ETA=list( list(X=T,model='BayesA')), 
              nIter=nIter,burnIn=burnIn,saveAt='ba')
    acc1[i]<-cor(t%*%fmBA$ETA[[1]]$b+fmBA$mu,yte)
    acc2[i]<-cor(T%*%fmBA$ETA[[1]]$b+fmBA$mu,ytr)
    ytest[,i]<-t%*%fmBA$ETA[[1]]$b+fmBA$mu;test[,i]=testing
    
    writeMat(sprintf('BA/mice/%s%d_mmice_BA.mat',name.of.trait,i),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
  }
  writeMat(sprintf('BA/mice/%s_mmice_BA.mat',name.of.trait),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
t1 = f(time1)
t1


mat<-readMat('cfw_R.mat')
mx=as.matrix(mat$x);rm(mat)
Y<-read.table('CFW_ALL_MLMM.txt',header = T)
mn<-20;
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:14) {
  luoy<-Y[,c(1,ij)]
  name.of.trait=noquote(names(luoy)[2])
  nn<-which(is.na(luoy[,2]));luoy<-luoy[-nn,]
  M=mx[-nn,];n=1161-length(nn);
  testing <- sample(n,round(n/5),replace=F)
  ytest<-matrix(0,length(testing),mn);test<-matrix(0,length(testing),mn)
  for (i in 1:20) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    #M<-as.matrix(mat$x)
    #n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBB=BGLR(y=ytr,ETA=list( list(X=T,model='BayesB')), 
              nIter=nIter,burnIn=burnIn,saveAt='bb')
    acc1[i]<-cor(t%*%fmBB$ETA[[1]]$b+fmBB$mu,yte)
    acc2[i]<-cor(T%*%fmBB$ETA[[1]]$b+fmBB$mu,ytr)
    ytest[,i]<-t%*%fmBB$ETA[[1]]$b+fmBB$mu;test[,i]=testing
    
    writeMat(sprintf('BB/mice/%s%d_mmice_BB.mat',name.of.trait,i),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
  }
  writeMat(sprintf('BB/mice/%s_mmice_BB.mat',name.of.trait),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
t1 = f(time1)
t1


mat<-readMat('cfw_R.mat')
mx=as.matrix(mat$x);rm(mat)
Y<-read.table('CFW_ALL_MLMM.txt',header = T)
mn<-20;
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:14) {
  luoy<-Y[,c(1,ij)]
  name.of.trait=noquote(names(luoy)[2])
  nn<-which(is.na(luoy[,2]));luoy<-luoy[-nn,]
  M=mx[-nn,];n=1161-length(nn);
  testing <- sample(n,round(n/5),replace=F)
  ytest<-matrix(0,length(testing),mn);test<-matrix(0,length(testing),mn)
  for (i in 1:20) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    #M<-as.matrix(mat$x)
    #n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBC=BGLR(y=ytr,ETA=list( list(X=T,model='BayesC')), 
              nIter=nIter,burnIn=burnIn,saveAt='bc')
    acc1[i]<-cor(t%*%fmBC$ETA[[1]]$b+fmBC$mu,yte)
    acc2[i]<-cor(T%*%fmBC$ETA[[1]]$b+fmBC$mu,ytr)
    ytest[,i]<-t%*%fmBC$ETA[[1]]$b+fmBC$mu;test[,i]=testing
    
    writeMat(sprintf('BC/mice/%s%d_mmice_BC.mat',name.of.trait,i),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
  }
  writeMat(sprintf('BC/mice/%s_mmice_BC.mat',name.of.trait),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
t1 = f(time1)
t1


mat<-readMat('cfw_R.mat')
mx=as.matrix(mat$x);rm(mat)
Y<-read.table('CFW_ALL_MLMM.txt',header = T)
mn<-20;
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:14) {
  luoy<-Y[,c(1,ij)]
  name.of.trait=noquote(names(luoy)[2])
  nn<-which(is.na(luoy[,2]));luoy<-luoy[-nn,]
  M=mx[-nn,];n=1161-length(nn);
  testing <- sample(n,round(n/5),replace=F)
  ytest<-matrix(0,length(testing),mn);test<-matrix(0,length(testing),mn)
  for (i in 1:20) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    #M<-as.matrix(mat$x)
    #n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-M[trainning,]
    t<-M[testing,]
    fmBL=BGLR(y=ytr,ETA=list( list(X=T,model='BL')), 
              nIter=nIter,burnIn=burnIn,saveAt='bl')
    acc1[i]<-cor(t%*%fmBL$ETA[[1]]$b+fmBL$mu,yte)
    acc2[i]<-cor(T%*%fmBL$ETA[[1]]$b+fmBL$mu,ytr)
    ytest[,i]<-t%*%fmBL$ETA[[1]]$b+fmBL$mu;test[,i]=testing
    
    writeMat(sprintf('BL/mice/%s%d_mmice_BL.mat',name.of.trait,i),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
  }
  writeMat(sprintf('BL/mice/%s_mmice_BL.mat',name.of.trait),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
t1 = f(time1)
t1




mat<-readMat('cfw_R.mat')
mx=as.matrix(mat$x);rm(mat)
Y<-read.table('CFW_ALL_MLMM.txt',header = T)
mn<-20;
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
time1 <- Sys.time()
for (ij in 2:2) {
  luoy<-Y[,c(1,ij)]
  name.of.trait=noquote(names(luoy)[2])
  nn<-which(is.na(luoy[,2]));luoy<-luoy[-nn,]
  M=mx[-nn,];n=1161-length(nn);
  testing <- sample(n,round(n/5),replace=F)
  ytest<-matrix(0,length(testing),mn);test<-matrix(0,length(testing),mn)
  for (i in 1:10) {
    #M<-as.matrix(mat$x)
    #n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-M[trainning,]
    t<-M[testing,]
    RR<- mixed.solve(y=ytr,Z=T) 
    acc1[i]<-cor(t%*%RR$u,yte)
    acc2[i]<-cor(T%*%RR$u,ytr)
    ytest[,i]<-t%*%RR$u;test[,i]<-testing
    writeMat(sprintf('RR/mice/%s%d_mmice_RR1.mat',name.of.trait,i),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
  }
  writeMat(sprintf('RR/mice/%s_mmice_RR2.mat',name.of.trait),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
t1 = f(time1)
t1




