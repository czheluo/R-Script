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

acc1<-matrix(0,20,1)
#ans2<-NULL
acc2<-matrix(0,20,1)
ytest<-matrix(0,37,20)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:2) {
  luoy<-Y[,c(1,ij)]
  
  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(mat$x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBB=BGLR(y=ytr,ETA=list( list(X=T,model='BayesB')), 
              nIter=nIter,burnIn=burnIn,saveAt='bb')
    acc1[i,]<-cor(t%*%fmBB$ETA[[1]]$b+fmBB$mu,yte)
    acc2[i,]<-cor(T%*%fmBB$ETA[[1]]$b+fmBB$mu,ytr)
    ytest[,i]<-t%*%fmBB$ETA[[1]]$b+fmBB$mu
  }
  #writeMat(paste('BB/y_',ij-1,'_cattle_BB.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
t1 = f(time1)
t1
