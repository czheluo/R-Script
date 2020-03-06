
#####human1000#####


setwd('J:\\human1000\\sim2\\')
library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

sim<-readMat('sim_human1000_60qtl_0.6_0.5_1rnd_1002_100_R.mat')
setwd('/home/luomeng/plink/human1000/out_plink/sim1')
time1 <- Sys.time()
plmr1<-matrix(0,502503,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,502503,100)
for (i in 1:100) {
  #setwd('/home/luomeng/plink/human1000/out_plink/sim2/1')
  plink.result1<-read.table(sprintf("1/20qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/20qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1002,i]<-plink.result2$PROP;plmr1[1003:502503,i]<-plink.result1$P;
  effr1[1003:502503,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim1_20qtl_0.6_0.5_1.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1


time1 <- Sys.time()
plmr1<-matrix(0,502503,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,502503,100)
for (i in 1:100) {
  #setwd('/home/luomeng/plink/human1000/out_plink/sim2/1')
  plink.result1<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1002,i]<-plink.result2$PROP;plmr1[1003:502503,i]<-plink.result1$P;
  effr1[1003:502503,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim2_60qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human1000/out_plink/sim3')
time1 <- Sys.time()
plmr1<-matrix(0,502503,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,502503,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1002,i]<-plink.result2$PROP;plmr1[1003:502503,i]<-plink.result1$P;
  effr1[1003:502503,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim3_100qtl_0.6_0.5.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human1000/out_plink/sim3')
time1 <- Sys.time()
plmr1<-matrix(0,502503,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,502503,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1002,i]<-plink.result2$PROP;plmr1[1003:502503,i]<-plink.result1$P;
  effr1[1003:502503,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim2_100qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human1000/out_plink/sim4')
time1 <- Sys.time()
plmr1<-matrix(0,502503,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,502503,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1002,i]<-plink.result2$PROP;plmr1[1003:502503,i]<-plink.result1$P;
  effr1[1003:502503,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim4_60qtl_0.6_0.5.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human1000/out_plink/sim4')
time1 <- Sys.time()
plmr1<-matrix(0,502503,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,502503,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1002,i]<-plink.result2$PROP;plmr1[1003:502503,i]<-plink.result1$P;
  effr1[1003:502503,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim4_60qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1


setwd('/home/luomeng/plink/human1000/out_plink/sim5')
time1 <- Sys.time()
plmr1<-matrix(0,502503,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,502503,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1002,i]<-plink.result2$PROP;plmr1[1003:502503,i]<-plink.result1$P;
  effr1[1003:502503,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim5_100qtl_0.6_0.5.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human1000/out_plink/sim5')
time1 <- Sys.time()
plmr1<-matrix(0,502503,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,502503,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1002,i]<-plink.result2$PROP;plmr1[1003:502503,i]<-plink.result1$P;
  effr1[1003:502503,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim5_100qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1


#####huamn2000####
setwd('/home/luomeng/plink/human2000/out_plink/sim1')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/20qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/20qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim1_20qtl_0.6_0.5_1.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human2000/out_plink/sim1')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/20qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/20qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim1_20qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1


setwd('/home/luomeng/plink/human2000/out_plink/sim2')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim2_60qtl_0.6_0.5_1.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human2000/out_plink/sim2')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim2_60qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1


setwd('/home/luomeng/plink/human2000/out_plink/sim3')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim3_100qtl_0.6_0.5_1.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human2000/out_plink/sim3')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim3_100qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1




setwd('/home/luomeng/plink/human2000/out_plink/sim4')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim4_60qtl_0.6_0.5_1.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human2000/out_plink/sim4')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim4_60qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1



setwd('/home/luomeng/plink/human2000/out_plink/sim5')
time1 <- Sys.time()
plmr1<-matrix(0,503506,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,503506,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1003,i]<-plink.result2$PROP;plmr1[1004:503506,i]<-plink.result1$P;
  effr1[1004:503506,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim5_100qtl_0.6_0.5_1.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/human2000/out_plink/sim5')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim5_100qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1


setwd('/home/luomeng/plink/rice1000/out_plink/sim1')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/20qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/20qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim1_20qtl_0.6_0.5.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/rice1000/out_plink/sim1')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/20qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/20qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim1_20qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/rice1000/out_plink/sim2')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim2_60qtl_0.6_0.5.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/rice1000/out_plink/sim2')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim2_60qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1



setwd('/home/luomeng/plink/rice1000/out_plink/sim3')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim3_100qtl_0.6_0.5.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/rice1000/out_plink/sim3')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim3_100qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1



setwd('/home/luomeng/plink/rice1000/out_plink/sim4')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/60qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim4_60qtl_0.6_0.5.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/rice1000/out_plink/sim4')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/60qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim4_60qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1



setwd('/home/luomeng/plink/rice1000/out_plink/sim5')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("1/100qtl_0.6_0.5_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim5_100qtl_0.6_0.5.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

setwd('/home/luomeng/plink/rice1000/out_plink/sim5')
time1 <- Sys.time()
plmr1<-matrix(0,508536,100)
#pver1<-matrix(0,502503,100)
effr1<-matrix(0,508536,100)
for (i in 1:100) {
  #
  plink.result1<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt",i),header = T)
  plink.result2<-read.table(sprintf("2/100qtl_0.6_0.8_%d.epi.qt.summary",i),header = T)
  plmr1[1:1008,i]<-plink.result2$PROP;plmr1[1009:508536,i]<-plink.result1$P;
  effr1[1009:508536,i]<-plink.result1$BETA_INT
}
#setwd('/home/luomeng/plink/human1000/out_plink/sim2')
writeMat(paste('sim5_100qtl_0.6_0.8.mat'),plmr1=plmr1,effr1=effr1)
t1 = f(time1)
t1

