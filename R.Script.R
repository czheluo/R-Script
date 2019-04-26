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
setwd(opt$input)
na<-list.files(path = getwd())
for(i in 1:length(na)){
	setwd(opt$output)
	map<-read.csv(na[i],header = T)
	nam<-strsplit(na[i],split=".",fixed=TRUE)[[1]][1]
	map1<-map[which(map[,5] %in% 1),c(3,6)]
	write.table(map1,file=(paste(nam,".female.map",sep="")),quote=F,sep="\t",row.names=F,col.names=F)
	map2<-map[which(map[,5] %in% 2),c(3,6)]
	write.table(map2,file=(paste(nam,".sexAver.map",sep="")),quote=F,sep="\t",row.names=F,col.names=F)
	map3<-map[which(map[,5] %in% 3),c(3,6)]
	write.table(map3,file=(paste(nam,".male.map",sep="")),quote=F,sep="\t",row.names=F,col.names=F)
	mat<-matrix(i,length(dim(map1)[1]))
	map4<-cbind(map1,mat)
	write.table(map4,file=paste(na[i],".1.mat",sep = ""),quote=F,sep = "\t",row.names = F,col.names = F)
	mat<-matrix(i,length(dim(map2)[1]))
	map5<-cbind(map2,mat)
	write.table(map5,file=paste(na[i],".2.mat",sep = ""),quote=F,sep = "\t",row.names = F,col.names = F)
	mat<-matrix(i,length(dim(map3)[1]))
	map6<-cbind(map3,mat)
	write.table(map6,file=paste(na[i],".3.mat",sep = ""),quote=F,sep = "\t",row.names = F,col.names = F)
}
nnSNP<-list()
chr<-unique(ntable$CHROM)
for (k in 1:length(chr)){
	nSNPs<-NULL
	int<-round((max(ntable[which(ntable[,1] %in% chr[k] ),2])-min(ntable[which(ntable[,1] %in% chr[k] ),2]))/1000000)+1
	pos2<-NULL
	pos2[1]<-ntable[which(ntable[,1] %in% chr[k] )[1],2]+1000000
	for (i in 2:int){
		pos2[i]<-pos2[i-1]+1000000
	}
	pos1<-matrix(0,int,1)
	pos1[1]=ntable[which(ntable[,1] %in% chr[k] )[1],2]
	pos1[c(2:int)]=pos2[c(1:int-1)]
	for (j in 1:int){
		nsnp<-length(which(ntable[which(ntable[,1] %in% chr[k] ),2]>=pos1[j] & ntable[which(ntable[,1] %in% chr[k] ),2]<pos2[j]))
		nSNPs[which(ntable[which(ntable[,1] %in% chr[k] ),2]>=pos1[j] & ntable[which(ntable[,1] %in% chr[k] ),2]<pos2[j])]<-nsnp
	}

	nnSNP[[k]]<-as.data.frame(nSNPs)
}
ntable$nSNPs<-do.call(rbind,unname(nnSNP))$nSNPs


png(paste("nSNPs",".png",sep=""),width=1000, height=800)
CMplot(data, plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,
       amplify=FALSE,main=expression("The Number of SNPs IN 1M windowsize"),chr.labels=NULL,
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)

dev.off()

pdf(paste("nSNPs",".pdf",sep=""),height=9,width=16)
CMplot(data, plot.type="m", LOG10=FALSE, line=2.5,cex.lab=2,
       amplify=FALSE,main=expression("The Number of SNPs IN 1M windowsize"),chr.labels=NULL,
       bin.size=1e6,type="l",multracks=FALSE,file.output=FALSE)
dev.off()

#library(data.table)
#ntable$nSNPs<-rbindlist(nnSNP)
#nSNPs<-do.call(rbind,unname(nnSNP))
escaptime<-Sys.time()-times;
print("Done!");
print(escaptime)



