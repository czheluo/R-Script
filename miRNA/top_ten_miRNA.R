args<-commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
stop("Rscript top_ten_miRNA.r --args <novel/known_miR_norm.xls>")
}
me<-read.table(args[2],sep="\t",head=T,row.names=1,check.names=F)
N<-ncol(me)
widthP<-ifelse(N>=7,7+1*(N-6),7)
heightP<-ifelse(N>=7,7+0.5*(N-6),7)
milist<-list()
library(ggplot2)
theme_set(theme_bw())
binddata<-data.frame()
for(i in 1:N){
milist[[i]]<-me[order(-me[,i]),][1:10,i,drop=F]
milist[[i]]$percent_weight<-((milist[[i]][,1])/sum(milist[[i]][,1]))*100
milist[[i]]$sample<-colnames(me)[i]
milist[[i]]$microRNA<-rownames(milist[[i]])
colnames(milist[[i]])[1]<-"counts"
rownames(milist[[i]])<-NULL
binddata<-rbind(binddata,milist[[i]])
}
if(length(args)==3){
    data<-binddata[which(binddata$sample %in% args[3]),]
    binddata<-data
}
library(reshape2)
mr<-dcast(binddata,microRNA~sample,value.var="percent_weight")
for(i in 1:nrow(mr)){
mr$sum[i]<-mean(as.numeric(as.character(mr[i,2:ncol(mr)])),na.rm=T)
}
mr<-mr[order(-mr$sum),]
binddata$microRNA<-factor(binddata$microRNA,levels=mr$microRNA)
length_col<-length(unique(binddata$microRNA))
widthP<-ifelse(length_col>15,widthP,widthP+3)
mycol <- c(34, 51, 142, 26, 31, 371, 36, 7, 12, 30, 84, 88, 116, 121, 77, 56, 386, 373, 423, 435, 438, 471, 512, 130, 52, 47, 6, 11, 43, 54, 367, 382, 422, 4
, 8, 375, 124, 448, 419, 614, 401, 403, 613, 583, 652, 628, 633, 496, 638, 655, 132, 503, 24)
colb<-colors()[mycol[1:length_col]]
pdfname<-paste(args[2],"top10.pdf",sep=".")
pdf(pdfname,width=widthP,height=heightP)
ggplot(binddata,aes(x=sample,y=percent_weight,fill=microRNA))+geom_bar(stat="identity")+scale_fill_manual(values=colb)+theme(axis.text.x  = element_text(size=5,angle=90, hjust=1))+
theme(axis.title.x=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(2,1,1,2),"cm"))+ylab("percent_weight(%)")+guides(fill=guide_legend(reverse=TRUE))+
ggtitle("Top 10 expressed miRNAs in each sample")
dev.off()
