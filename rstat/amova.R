#!/usr/bin/env Rscript
##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 Majorbio
#  contact: meng.luo@majorbio.com
##########################################

library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',1,'character',
	'output','o',1,'character',
	'gro','g',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:
	Rscript Ammova.R --input pop.recode.vcf --output outputdir --gro groupdir
Usage:
	--input the pop.recode.vcf file dir
	--output	output dir
	--gro 	the group list file dir
	--help		usage
\n")
	q(status=1);
}
if (is.null(opt$input)) { print_usage(spec)}
if (is.null(opt$output)){ print_usage(spec) }
if (is.null(opt$gro)){ print_usage(spec) }

times<-Sys.time()

setwd(opt$output)

library(pegas)
library(adegenet)
library(poppr)
#info <- VCFloci("/mnt/ilustre/users/tong.wang/Project/MJ20180320018-yaozhangxiu/evo20180814/step01.vcf-filter/pop.recode.vcf")
info <-VCFloci(paste(opt$input,"pop.recode.vcf",sep=""))
popdata<-read.table(paste(opt$gro,"group.list",sep=""))
colnames(popdata)=c("sampleID","groupID")
SNP <- is.snp(info)
x <- length(which(SNP))
x<-read.vcf(paste(opt$input,"pop.recode.vcf",sep=""), from = 1, to = x)
g<-loci2genind(x)
g<-as.genclone(g)
ploidy(g)<-2
strata(g)<-data.frame(popdata)
amovacc<-poppr.amova(g,~groupID)
signif   <- randtest(amovacc, nrepet = 999)
pdf("Amova.signif.pdf")
plot(signif)
dev.off()
png("Amova.signif.png")
plot(signif)
dev.off()
d<-cbind(signif$names,signif$alter,signif$expvar,signif$pvalue)
colnames(d)<-c("Test","Alter","Obs","Expectation","Variance","Pvalue")
write.table(file="Amova.signif",d)
write.table(cbind(amovacc$results,amovacc$componentsofcovariance), file = "AMOVA.csv")
write.table(amovacc$statphi, file = "AMOVA.phi")

escaptime=Sys.time()-times;
print("Done!");
print(escaptime)

