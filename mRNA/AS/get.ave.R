

gro<-read.table("group.list",header=F)
g<-unique(gro[,1])
for (i in 1:length(g)){
	print(g[i])
	loc1<-read.table(paste(g[i],"_1/",g[i],"_1_AS_statistics.xls",sep=""),sep="\t",header=T)
	loc2<-read.table(paste(g[i],"_2/",g[i],"_2_AS_statistics.xls",sep=""),sep="\t",header=T)
	loc3<-read.table(paste(g[i],"_3/",g[i],"_3_AS_statistics.xls",sep=""),sep="\t",header=T)
	p <- round((loc1[,2]+loc2[,2]+loc3[,2])/3)
	colnames(p)<-colnames(loc1[,2])
	p2<-cbind(as.character(loc1[,1]),p)
	colnames(p2)<-colnames(loc1)
	dir.create(file.path(getwd(),g[i]), showWarnings = FALSE)
	write.table(p2,paste(g[i],"/",g[i],"_AS_statistics.xls",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

}
