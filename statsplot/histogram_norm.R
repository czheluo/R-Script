# Plotting normal curve over histogram using ggplot2 
# and normal distribution test
#!/usr/bin/env Rscript
##library(ggplot2)
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
  'input','i',1,'character',
  'output','o',0,'character',
  'help','h',0,'logical'
), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
      Usage example: 
      Rscript emap.R --input --output --HB --LB
      Usage:
      --input *vcf.table
      --output	output dir
      --help		usage
      \n")
  q(status=1);
}

times<-Sys.time()
Trait <- read.csv(opt$input,header=T)
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
  png( paste(name.of.trait,'.png'))
  mn<- min[i-1];mx<-max[i-1];
  x <- Trait[,i] 
  pdf(paste(name.of.trait,'.pdf'))
  h<-hist(x, breaks=20, col="blue", xlab= name.of.trait, 
          main="Histogram with Normal Curve") 
  xfit<-seq(min(x),max(x),length=100) 
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
  yfit <- yfit*diff(h$mids[1:2])*length(x) 
  lines(xfit, yfit, col="black", lwd=2)
  dev.off()
  png( paste(name.of.trait,'.png'))
  h<-hist(x, breaks=20, col="blue", xlab="name.of.trait", 
          main="Histogram with Normal Curve") 
  xfit<-seq(min(x),max(x),length=100) 
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
  yfit <- yfit*diff(h$mids[1:2])*length(x) 
  lines(xfit, yfit, col="black", lwd=2)
  dev.off()
}

W_VALUE<-NULL;P_VALUE<-NULL;
for (i in 1:4) {
  W_VALUE[i]<-nt1[[i]]
  P_VALUE[i]<-nt2[[i]]
}
cb<-cbind(min,max,mean,sd,W_VALUE,P_VALUE)
rownames(cb)<-colnames(Trait[,-1])
write.csv(cb,file = 'Trait_summary.csv',quote = F,row.names = T)

escaptime<-Sys.time()-times;
print("Done!");
print(escaptime)
