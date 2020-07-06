# Plotting normal curve over histogram using ggplot2 and normal distribution test

library(ggpubr)
library(ggplot2)
binwidth <- 4 ### the number of bar
Trait <- read.csv("Trit.csv",header=T)
resulttable <- matrix(0,dim(Trait)[2],6)
n <- dim(Trait)[1]
mean<-NULL;max<-NULL;sd <- NULL;nt1<-NULL;nt2<-NULL;min<-NULL;
for (i in 2:dim(Trait)[2]) {
  Y <- Trait[, c(1, i)]
  name.of.trait <- noquote(names(Y)[2])
  mean[i-1]<-mean(Trait[,i])
  sd[i-1]<-sd(Trait[,i])
  min[i-1] <- min(Trait[,i])
  max[i-1] <- max(Trait[,i])
  nt1[i-1]<- shapiro.test(Trait[,i])[1]
  nt2[i-1]<- shapiro.test(Trait[,i])[2]
  #png( paste(name.of.trait,'.png'))
  mn<- min[i-1];mx<-max[i-1];
  qplot(Trait[i], geom = "histogram", breaks = seq(min[i-1], max[i-1], binwidth),
        colour = I("black"), fill = I("white"),
        xlab = name.of.trait, ylab = "Count") +
    stat_function(
      fun = function(x, mean, sd, n, bw){
        dnorm(x = x, mean = mean[i-1], sd = sd[i-1]) * n * bw
      },
      args = c(mean = mean[i-1], sd = sd[i-1], n = n, bw = binwidth))
  #dev.off()
  #pdf(paste(name.of.trait,'.pdf'))
  ggsave(paste(name.of.trait,'.png'))
  ggsave(paste(name.of.trait,'.png'))

}
W_VALUE<-NULL;P_VALUE<-NULL;
for (i in 1:4) {
  W_VALUE[i]<-nt1[[i]]
  P_VALUE[i]<-nt2[[i]]
  }
cb<-cbind(min,max,mean,sd,W_VALUE,P_VALUE)
rownames(cb)<-colnames(Trait[,-1])
write.csv(cb,file = 'Trait_summary.csv',quote = F,row.names = T)

ph<-read.csv("traits.csv",header = T,na.strings = F)
#Big pvalue and the normer,

ph<-read.csv("traits.csv",header = T)
nor<-NULL
for (i in 2:16) {
  test<-shapiro.test(ph[,i])
  nor[i-1]<-test$p.value
}

norm.test(ph[,i])
library(rMVP)
MVP.Hist(phe=ph, file.type="jpg", breakNum=18, dpi=300)


