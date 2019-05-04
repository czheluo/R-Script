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
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

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



geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim05qtn500.pheno.txt'
out.file<-'BS_cattle_5024sim05qtn500'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.5_qtn500.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 15:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim05qtn500.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim05qtn500.param.txt -emu output/BS_cattle_5024sim05qtn500.log.txt -ebv output/BS_cattle_5024sim05qtn500.bv.txt -k output/BS_cattle_5024sim05qtn500.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim05qtn500.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim05/sim05qtn500',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim05/sim05qtn500_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1


geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim02qtn500.pheno.txt'
out.file<-'BS_cattle_5024sim02qtn500'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.2_qtn500.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 14:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim02qtn500.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim02qtn500.param.txt -emu output/BS_cattle_5024sim02qtn500.log.txt -ebv output/BS_cattle_5024sim02qtn500.bv.txt -k output/BS_cattle_5024sim02qtn500.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim02qtn500.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim02/sim02qtn500',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim02/sim02qtn500_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1


library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim08qtn150.pheno.txt'
out.file<-'BS_cattle_5024sim08qtn150'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.8_qtn150.csv',header = TRUE)
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
    write.table(y,file = 'BS_cattle_5024sim08qtn150.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim08qtn150.param.txt -emu output/BS_cattle_5024sim08qtn150.log.txt -ebv output/BS_cattle_5024sim08qtn150.bv.txt -k output/BS_cattle_5024sim08qtn150.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim08qtn150.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim08/sim08qtn150',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim08/sim08qtn150_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1


library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim05qtn150.pheno.txt'
out.file<-'BS_cattle_5024sim05qtn150'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.5_qtn150.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 15:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim05qtn150.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim05qtn150.param.txt -emu output/BS_cattle_5024sim05qtn150.log.txt -ebv output/BS_cattle_5024sim05qtn150.bv.txt -k output/BS_cattle_5024sim05qtn150.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim05qtn150.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim05/sim05qtn150',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim05/sim05qtn150_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1



library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim02qtn150.pheno.txt'
out.file<-'BS_cattle_5024sim02qtn150'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.2_qtn150.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 15:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim02qtn150.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim02qtn150.param.txt -emu output/BS_cattle_5024sim02qtn150.log.txt -ebv output/BS_cattle_5024sim02qtn150.bv.txt -k output/BS_cattle_5024sim02qtn150.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim02qtn150.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim02/sim02qtn150',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim02/sim02qtn150_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1

library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024.pheno.txt'
out.file<-'BS_cattle_5024'
Y<-read.table('FileS2.txt',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 2:2) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:20) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024.param.txt -emu output/BS_cattle_5024.log.txt -ebv output/BS_cattle_5024.bv.txt -k output/BS_cattle_5024.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024.prdt.txt',header = FALSE)
    acc1[j]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[j]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,j]<-pred$V1[testing];test[,j]=testing
  }
  writeMat(paste('BSLMM/y',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
#writeMat(paste('BSLMM/_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1


geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024.pheno.txt'
out.file<-'BS_cattle_5024'
Y<-read.table('FileS2.txt',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 3:3) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:20) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024.param.txt -emu output/BS_cattle_5024.log.txt -ebv output/BS_cattle_5024.bv.txt -k output/BS_cattle_5024.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024.prdt.txt',header = FALSE)
    acc1[j]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[j]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,j]<-pred$V1[testing];test[,j]=testing
  }
  writeMat(paste('BSLMM/y',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
#writeMat(paste('BSLMM/_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1


geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024.pheno.txt'
out.file<-'BS_cattle_5024'
Y<-read.table('FileS2.txt',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 4:4) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:20) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024.param.txt -emu output/BS_cattle_5024.log.txt -ebv output/BS_cattle_5024.bv.txt -k output/BS_cattle_5024.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024.prdt.txt',header = FALSE)
    acc1[j]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[j]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,j]<-pred$V1[testing];test[,j]=testing
  }
  writeMat(paste('BSLMM/y',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
#writeMat(paste('BSLMM/_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1

geno.file<-'HDRA-G6-4-RDP1-RDP2-NIAS-filter-pheo-maf0.1'
#pheno.file<-'cattle_5024sim05qtn500.pheno.txt'
out.file<-'rice1000_rootlength_BS'
map<-read.table('rice_pheo-maf0.1.fam',header = F)
luoy<-data.frame(map$V1,map[,6]);n<-dim(luoy)[1]
#mn<-20;
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,226,20);test<-matrix(0,226,20)
time1 <- Sys.time()
for (j in 1:20) {
  #luoy<-Y[,c(1,ij)]
  #for (j in 1:1) {
  y<-matrix(NA,nrow(luoy),1)
  testing <- sample(n,round(n/5),replace=F)
  y[-testing]<-luoy[,2][-testing];map[,6]<-y
  write.table(map,file = 'HDRA-G6-4-RDP1-RDP2-NIAS-filter-pheo-maf0.1.fam',quote = FALSE,col.names = FALSE,row.names = FALSE)
  #y<-matrix(NA,nrow(luoy),1)
  #testing <- sample(n,round(n/5),replace=F)
  #y[-testing]<-luoy[,2][-testing]
  #write.table(y,file = 'BS_cattle_5024.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
  system(sprintf("./gemma -bfile %s -n 1 -gk 1 -maf 0.00001 -o %s",
                 geno.file,out.file),ignore.stdout = FALSE)
  system(sprintf("./gemma -bfile %s -n 1 -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                 geno.file,out.file),ignore.stdout = FALSE)
  system(sprintf("./gemma -bfile %s -n 1 -epm output/rice1000_rootlength_BS.param.txt -emu output/rice1000_rootlength_BS.log.txt -ebv output/rice1000_rootlength_BS.bv.txt -k output/rice1000_rootlength_BS.cXX.txt -predict 1 -o %s",
                 geno.file,out.file),ignore.stdout = FALSE)
  pred<-read.table(file='output/rice1000_rootlength_BS.prdt.txt',header = FALSE)
  acc1[j]<-cor(pred$V1[testing],luoy[,2][testing])
  acc2[j]<-cor(pred$V1[-testing],luoy[,2][-testing])
  ytest[,j]<-pred$V1[testing];test[,j]=testing
  
    #}
  writeMat(paste('rice/rice1000_rootlength_BS',j,'gemma.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('rice/rice1000_rootlength_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1
timess<-Sys.time()-times
print(timess)
