library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}


####simulation#####

geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'cattle_sim08_qtn500_va.pheno.txt'
out.file<-'cattle_sim08_qtn500_va'
Y<-read.csv('cattle_sim0.8_qtn500_va.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'cattle_sim08_qtn500_va.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./DPR -g %s -p %s -dpr 2 -maf 0.00000001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./DPR -g %s -p %s -epm output/cattle_sim08_qtn500_va.param.txt -emu output/cattle_sim08_qtn500_va.log.txt -predict -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/cattle_sim08_qtn500_va.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing]
  }
  writeMat(paste('cattle/sim08/sim08_qtn500_va',ij-1,'_cattle_DPR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
writeMat(paste('cattle/sim08/cattle_sim08_qtn500_va_DPR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
t1 = f(time1)
t1



####simulation#####

geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'cattle_sim02_qtn500_va.pheno.txt'
out.file<-'cattle_sim02_qtn500_va'
Y<-read.csv('cattle_sim0.2_qtn500_va.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'cattle_sim02_qtn500_va.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./DPR -g %s -p %s -dpr 2 -maf 0.00000001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./DPR -g %s -p %s -epm output/cattle_sim02_qtn500_va.param.txt -emu output/cattle_sim02_qtn500_va.log.txt -predict -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/cattle_sim02_qtn500_va.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing]
  }
  writeMat(paste('cattle/sim02/sim02_qtn500_va',ij-1,'_cattle_DPR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
writeMat(paste('cattle/sim02/cattle_sim02_qtn500_va_DPR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
t1 = f(time1)
t1


####simulation#####

geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'cattle_sim05_qtn500_va.pheno.txt'
out.file<-'cattle_sim05_qtn500_va'
Y<-read.csv('cattle_sim0.5_qtn500_va.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'cattle_sim05_qtn500_va.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./DPR -g %s -p %s -dpr 2 -maf 0.00000001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./DPR -g %s -p %s -epm output/cattle_sim05_qtn500_va.param.txt -emu output/cattle_sim05_qtn500_va.log.txt -predict -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/cattle_sim05_qtn500_va.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing]
  }
  writeMat(paste('cattle/sim05/sim05_qtn500_va',ij-1,'_cattle_DPR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
writeMat(paste('cattle/sim05/cattle_sim05_qtn500_va_DPR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
t1 = f(time1)
t1
