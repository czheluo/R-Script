
library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

Y<-read.csv('PHS_pheno.csv',header = TRUE)
geno.file<-'PHS_wheat.txt.gz'
pheno.file<-'PHS_wheat.pheno.txt'
out.file<-'PHS_wheat'
mn<-20
acc1<-matrix(0,10,mn)
acc2<-matrix(0,10,mn)
ytest<-matrix(0,185,mn)
nfolds<-10
time1 <- Sys.time()
for (ij in 11:11) {
  luoy<-Y[,c(1,ij)]
  
  for (j in 1:mn) {
    
    #luoy<-luoy[sample(nrow(luoy)),]
    folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    nfolds<-10
    
    time1 <- Sys.time()
    for (i in 1:nfolds) {
      y<-matrix(NA,185,1)
      testing<- which(folds==i)
      y[-testing]<-luoy[,2][-testing]
      write.table(y,file = 'PHS_wheat.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
      system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                     geno.file,pheno.file,out.file),ignore.stdout = FALSE)
      system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                     geno.file,pheno.file,out.file),ignore.stdout = FALSE)
      system(sprintf("./gemma -g %s -p %s -epm output/PHS_wheat.param.txt -emu output/PHS_wheat.log.txt -ebv output/PHS_wheat.bv.txt -k output/PHS_wheat.cXX.txt -predict 1 -o %s",
                     geno.file,pheno.file,out.file),ignore.stdout = FALSE)
      pred<-read.table(file='output/PHS_wheat.prdt.txt',header = FALSE)
      acc1[i,j]<-cor(pred$V1[testing],luoy[,2][testing])
      acc2[i,j]<-cor(pred$V1[-testing],luoy[,2][-testing])
      ytest[testing,j]<-pred$V1[testing]
    }
  }
  writeMat(paste('BSLMM/PHS/y_',ij-1,'_PHS_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
t1 = f(time1)
t1


geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim08qtn500.pheno.txt'
out.file<-'BS_cattle_5024sim08qtn500'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.8_qtn500.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 16:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim08qtn500.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim08qtn500.param.txt -emu output/BS_cattle_5024sim08qtn500.log.txt -ebv output/BS_cattle_5024sim08qtn500.bv.txt -k output/BS_cattle_5024sim08qtn500.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim08qtn500.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim08/sim08qtn500',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim08/sim08qtn500_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1

