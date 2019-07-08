##########################################
#  https://github.com/czheluo/
#  Copyright (C) 2019  by Meng Luo
#  contact: czheluo@gmail.com
##########################################

library(RColorBrewer)
gro<-read.table("group.list",header=F)
g<-unique(gro[,1])
for (i in 1:length(g)){
	print(g[i])
	loc1<-read.table(paste(g[i],"_1_rfam_stat.xls",sep=""),sep="\t",header=T)
	loc2<-read.table(paste(g[i],"_2_rfam_stat.xls",sep=""),sep="\t",header=T)
	loc3<-read.table(paste(g[i],"_3_rfam_stat.xls",sep=""),sep="\t",header=T)
	loc1[,c(3,5)]<-data.frame(lapply(loc1[,c(3,5)], function(x) as.numeric(sub("%", "", x))) )
	loc2[,c(3,5)]<-data.frame(lapply(loc2[,c(3,5)], function(x) as.numeric(sub("%", "", x))) )
	loc3[,c(3,5)]<-data.frame(lapply(loc3[,c(3,5)], function(x) as.numeric(sub("%", "", x))) )
	pal1<-round((loc1[,c(3,5)]+loc2[,c(3,5)]+loc3[,c(3,5)])/3,2)
	namety<-as.character(loc1[,1])
	palp1<-paste(pal1[,1],"%",sep="")
	palp2<-paste(pal1[,2],"%",sep="")
	pal2<-round((loc1[,c(2,4)]+loc2[,c(2,4)]+loc3[,c(2,4)])/3)
	pall<-cbind(namety,pal2[,1],palp1,pal2[,2],palp2)
	p <- pal2[,2]
	#rownames(pal)<-rownames(loc1)
	colnames(pall)<-colnames(loc1)
	print(pall)
	p1<- loc1[,4]
	p2<- loc2[,4]
	p3<- loc3[,4]
	#print(p)
	write.table(pall,file=paste(g[i],"_rfam.xls",sep=""),sep="\t",quote=F,col.name=T,row.name=F)
	pdf(paste(g[i],"_rfam.pdf",sep=""),width=10)
	#slices<-c(pall[11,4],pall[14,4],pall[13,4],pall[1,4]-pall[2,4]-pall[11,4]-pall[13,4]-pall[14,4],pall[1,4]-pall[2,4])
	slices<-c(p[11],p[14],p[13],p[2]-p[11]-p[14]-p[13],p[1]-p[2])
	print(slices)
	lbls<-c("rRNA","miRNA","tRNA","Other","Unmatched")
	pct<-round(slices/sum(slices)*100)
	lbls<-paste(lbls,pct)
	lbls<-paste(lbls,"%",sep="")
	ColorList<-brewer.pal(length(lbls),"Paired")
	pie(slices,labels=lbls,col=ColorList,main=paste("Pie Chart of ",g[i],sep=""),border="white")
	dev.off()
}
