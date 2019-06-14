

gro<-read.table("group.list",header=F)
g<-unique(gro[,1])
for (i in 1:length(g)){
	print(g[i])
	loc1<-read.table(paste(g[i],"_1_loc_bias_per.xls",sep=""),sep="\t",header=T)
	loc2<-read.table(paste(g[i],"_2_loc_bias_per.xls",sep=""),sep="\t",header=T)
	loc3<-read.table(paste(g[i],"_3_loc_bias_per.xls",sep=""),sep="\t",header=T)
	p1<- t(loc1)[-1,]
	colnames(p1) <- t(loc1)[1,]
	p2<- t(loc2)[-1,]
	colnames(p2) <- t(loc2)[1,]
	p3<- t(loc3)[-1,]
	colnames(p3) <- t(loc3)[1,]
	p <- (p1+p2+p3)/3
	write.table(t(p),paste(g[i],"_loc_bias_per.xls",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
	pdf(paste(g[i],"_loc_bias_per.pdf",sep=""),10,6)
	par(pin=c(5,1.8),fig=c(0,0.95,0,1),xpd=T,cex.axis=0.8)
	barplot(p,col=c("blue","yellow","red","purple"),yaxt="n",ann=FALSE)
	title(main="miRNA nucleotide bias at each position")
	axis(2,at=seq(0,1,0.25),labels=seq(0,100,25))
	legend(32+6,1,c("A","G","C","U"),col=c("blue","yellow","red","purple"),pch=c(15),bty="n",cex=0.7)
	mtext("Percent(%)",side=2,las=0,cex.lab=0.7,line=2)
	mtext("Position",side=1,las=0,cex.lab=0.7,line=2)
	segments(-1.5,0,32+5,0)
	dev.off()
	loc1<-read.table(paste(g[i],"_1_first_bias_per.xls",sep=""),sep="\t",header=T)
	loc2<-read.table(paste(g[i],"_2_first_bias_per.xls",sep=""),sep="\t",header=T)
	loc3<-read.table(paste(g[i],"_3_first_bias_per.xls",sep=""),sep="\t",header=T)
	p1<- t(loc1)[-1,]
	colnames(p1) <- t(loc1)[1,]
	p2<- t(loc2)[-1,]
	colnames(p2) <- t(loc2)[1,]
	p3<- t(loc3)[-1,]
	colnames(p3) <- t(loc3)[1,]
	p <- (p1+p2+p3)/3
	write.table(t(p),paste(g[i],"_first_bias_per.xls",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
	pdf(paste(g[i],"_first_bias_per.pdf",sep=""),10,6)
	par(pin=c(5,1.8),fig=c(0,0.95,0,1),xpd=T,cex.axis=0.8)
	barplot(p,col=c("blue","yellow","red","purple"),yaxt="n",ann=FALSE)
	title(main="miRNA nucleotide bias at each position")
	axis(2,at=seq(0,1,0.25),labels=seq(0,100,25))
	legend(32+6,1,c("A","G","C","U"),col=c("blue","yellow","red","purple"),pch=c(15),bty="n",cex=0.7)
	mtext("Percent(%)",side=2,las=0,cex.lab=0.7,line=2)
	mtext("Position",side=1,las=0,cex.lab=0.7,line=2)
	segments(-1.5,0,32+5,0)
	dev.off()
}
