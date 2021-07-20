  
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
data<-read.table("express.txt",header=T,stringsAsFactors=F,sep="\t")  
##volcano plot here##
#anno<-read.table("anno.txt",header=T,stringsAsFactors=F,sep="\t")  
data$LogPval <- -log10(data[,"padjust"])
data$Log2FC <- data[,grep("^log2",colnames(data))]
data$type <- factor(data$regulate,levels=c("nosig","up","down"))

Log2FC <- data[,grep("^log2",colnames(data))]
LogPval <- -log10(data[,"padjust"])
LogPval <- -log10(data[,"pvalue"])
res<-as.data.frame(cbind(Log2FC=Log2FC,LogPval=LogPval,anno=data[,"anno"],gene=data[,1]))
res$type <- factor(data$regulate,levels=c("nosig","up","down"))
anno<-cbind(anno=data[,"anno"],gene=data[,1])
#rownames(res)<-data[,1]
#res[which(res$LogPval>5),2]=5
volcano<-ggplot(data,aes(x=Log2FC,y=LogPval))+
	geom_point(aes(color=type,shape=type))+
	scale_color_manual(values=c("up"="#FF0000","down"="#009933","nosig"="grey"))+
	scale_shape_manual(values=c("up"=18,"down"=15,"nosig"=19))+
	labs(x=paste0("Log2FC"),y=paste0("-Log10(Padjust)"))+#,title=paste0("control_vs_sema4D",".volcano\n"))+
	#scale_y_continuous(expand=c(0,0),limits = c(0,ceiling(max(res$LogPval))+1))+
	scale_y_continuous(expand=c(0,0))+
	#scale_x_continuous(expand=c(0,0),limits = c(floor(min(res$Log2FC)),ceiling(max(res$Log2FC))))+
	scale_x_continuous(expand=c(0,0))+
	guides(shape=guide_legend(override.aes=list(size=3)))+
	geom_text_repel(data = subset(data, anno == 1),
			aes(label = seq_id),size = 5,	
			box.padding = unit(1.4, "lines"),
			point.padding = unit(1.3, "lines")
			)+
	theme(panel.border = element_blank(),
		panel.grid = element_blank(),
		axis.line = element_line(),
		panel.background =  element_blank(),
		plot.background = element_rect(fill = "transparent", colour = NA)) +
		theme(plot.margin = margin(2.5, 2.5, 2, 2.5, "cm"))+
		theme(legend.position=c(1,1),legend.title=element_blank(),
		legend.key = element_blank(),legend.text = element_text(size = 14))+
	theme(axis.text.x=element_text(size = 12),axis.title.x = element_text(size = 14, color = "black"),
		axis.text.y=element_text(size = 12),axis.title.y = element_text(size = 14, color = "black"),
		plot.title = element_text(hjust = 0.5,size = 14, color = "black"))
ggsave(paste0("test1",".diffexp.volcano.png"), height=9, width=8,plot = volcano)
ggsave(paste0("test1",".diffexp.volcano.pdf"), height=9, width=8,plot = volcano)
## MA plot ##
res$express <- apply(res[,2:(nrow(colData)+1)], 1, mean)
MA<-ggplot(res,aes(x=log10(express),y=Log2FC))+
	geom_point(aes(color=type,shape=type))+
	scale_color_manual(values=c("up"="#FF0000","down"="#009933","nosig"="grey"))+
	scale_shape_manual(values=c("up"=18,"down"=15,"nosig"=19))+
	labs(x=paste0("Log10(",exp_type,")"),y=paste0("Log2(",fc,")"),title=paste0(cmp,".MA\n"))+ 
	scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
	guides(shape=guide_legend(override.aes=list(size=3)))+
	theme(panel.border = element_blank(),
		panel.grid = element_blank(),
		axis.line = element_line(),
		panel.background =  element_blank(),
		plot.background = element_rect(fill = "transparent", colour = NA)) +
	theme(plot.margin = margin(2.5, 2.5, 2, 2.5, "cm"))+
	theme(legend.position=c(1,1),legend.title=element_blank(),
		legend.key = element_blank(),legend.text = element_text(size = 14))+
	theme(axis.text.x=element_text(size = 12),axis.title.x = element_text(size = 14, color = "black"),
		axis.text.y=element_text(size = 12),axis.title.y = element_text(size = 14, color = "black"),
		plot.title = element_text(hjust = 0.5,size = 14, color = "black"))
ggsave(paste0("test",".diffexp.scatter.pdf"), height=9, width=8,plot = MA)

