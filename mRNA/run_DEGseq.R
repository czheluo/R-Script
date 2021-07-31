#!/usr/bin/env Rscript
# Author meng.luo
# update 2021/07/22

suppressPackageStartupMessages(library(optparse))
#suppressPackageStartupMessages(library(DESeq2))
#suppressPackageStartupMessages(library(doParallel))
#suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DEGseq))


########parameter check
timeStart<-Sys.time()
Usage="Rscript %prog -r counts.xls -g group.txt -c compare.txt [-e tpm.xls -t TPM -a annotation.xls -q 0.05 -f 2]"
option_list = list(
    make_option(c("-r", "--readcounts"), type="character", default=NULL,
          help="read counts matrix file name", metavar="character"),
    make_option(c("-g", "--group"), type="character", default=NULL, 
          help="sample group information file name with header: sample<tab>group<tab>batch", metavar="character"),
    make_option(c("-c", "--compare"), type="character", default=NULL, 
          help="group comparison information file name with header: control<tab>case", metavar="character"),
    make_option(c("-f", "--foldchange"),type="double", default=2.0,
          help="[Optional]foldchange threshold [default %default]", metavar="double"),
    make_option(c("-p", "--pvalue"), type="double",default= 0.05,
          help="[Optional]pvalue threshold [default %default]", metavar="double"),
    make_option(c("-q", "--padjust"), type="double", default=NULL,
          help="[Optional]padjuast threshold.Filtering can be performed using any one of (-p), (-q) at a time", metavar="double"),
    make_option(c("-e", "--expression"), type="character", default=NULL, 
          help="[Optional]expression matrix file name", metavar="character"),
    make_option(c("-t", "--exp_type"), type="character", default="Express", 
          help="[Optional]express type: FPKM/TPM/RPM/SRPBM.[default %default]", metavar="character"),
    make_option(c("-a", "--annotation"),type="character", default=NULL,
          help="[Optional]annotation file name", metavar="character"),
    make_option(c("-o", "--outputdir"),type="character", default="./",
          help="[Optional]output directory for results [default %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list,usage=Usage);
opt = parse_args(opt_parser);
if (is.null(opt$readcounts) |is.null(opt$compare)){
    print_help(opt_parser)
    stop("read counts matrix file, group and comparison information file must be supplied", call. = FALSE)
}

count<-read.table(normalizePath(opt$readcounts),header=T)
compare<-read.table(normalizePath(opt$compare),header=F,stringsAsFactors = FALSE)

#count<- round(read.delim(normalizePath(opt$readcounts), row.names = 1,strip.white = TRUE))
#compare <- read.delim(normalizePath(opt$compare),comment.char = "",stringsAsFactors = FALSE)
if (! is.null(opt$annotation)){annotation <- read.delim(normalizePath(opt$annotation), quote = "")}else{annotation <- NULL}
if (! is.null(opt$expression)){expression <- read.delim(normalizePath(opt$expression), quote = "")}else{expression <- NULL}
if (! dir.exists(opt$outputdir)){dir.create(opt$outputdir)}
group <- read.delim(normalizePath(opt$group),comment.char = "",stringsAsFactors=FALSE, quote = "",col.names=c("sample","group"))

for(i in 1:dim(compare)[1]){
	ctrl <- compare[i,1]
	case <- compare[i,2]
	rownames(group)<-group[,1]
	case.group <- group[which(group$group == case),]
	ctrl.group <- group[which(group$group == ctrl),]
	colData <- rbind(case.group[order(rownames(case.group)),],ctrl.group[order(rownames(ctrl.group)),])
	countData <- count[,rownames(colData)]
	colData <- apply(colData, 2, factor)
	DEGexp(geneExpMatrix1=count, expCol1=compare[i,1], groupLabel1=compare[i,1], 
		  geneExpMatrix2=count, expCol2=compare[i,2], groupLabel2=compare[i,2],
		  method="LRT",outputDir=paste(compare[i,1],".vs.",compare[i,2],sep=""),
		  thresholdKind=3,qValue=0.05,normalMethod="loess")#"none", "loess", "median" #fdr
	result<-read.table(file=file.path(paste(compare[i,1],".vs.",compare[i,2],sep=""),"output_score.txt"),header = T,sep = "\t")
	colnames(result)<-c("GeneNames", compare[i,1], compare[i,2], "log2FoldChange", "log2Foldchange_normalized","Pvalue","Padjust","qval2003","Signature")
	rownames(result)<-result[,1]
	res <- as.data.frame(result[count[,1],])
	res[, "FoldChange"] <- 2 ^ res[, "log2FoldChange"]
	res <- res[, c("GeneNames", "FoldChange", "log2FoldChange", "Pvalue","Padjust")]
	#res[is.na(res)]<-1.0000
	fc <- paste0("FoldChange(",compare[i,2],"/",compare[i,1],")")
	colnames(res) <- c("seq_id",fc, paste0("log2(",fc,")"), "Pvalue", "Padjust")
	Pval <- ifelse(! is.null(opt$padjust),"Padjust","Pvalue")
	res$Significant <- "no"
	res[which(res[,Pval] < opt$padjust & abs(res[,3]) > log2(opt$foldchange)),"Significant"] <- "yes"
	#if (! is.null(opt$padjust)){res[which(res[,Pval] < opt$padjust & abs(res[,3]) > log2(opt$foldchange)),"Significant"] <- "yes"}
	res$Regulate <-ifelse(res[,2]==1,"no change", ifelse(res[,2] > 1,"up","down"))
	
	if (! is.null(expression)){
		exp.tmp <- expression[,c(colnames(expression)[1],rownames(colData))]
		colnames(exp.tmp)[2:ncol(exp.tmp)] <- paste(colnames(exp.tmp)[2:ncol(exp.tmp)],opt$exp_type,sep = "_")
		res <- merge(exp.tmp,res,by.y=colnames(res)[1],by.x=colnames(exp.tmp)[1],all.y=TRUE)
	}
	if (! is.null(annotation)) {res <- merge(res,annotation,by.x=colnames(res)[1],by.y=colnames(annotation)[1],all.x=TRUE)}
	write.table(res,paste0(opt$outputdir,"/",compare[i,1],".vs.",compare[i,2],".DEGseq.annot.xls"),sep ="\t", col.names=TRUE,row.names=FALSE,quote=FALSE,na = "")
	resSig <- res[which(res$Significant=="yes"),]
	resSig_out <- paste0(opt$outputdir,"/",compare[i,1],".vs.",compare[i,2],".DEGseq.pval", opt$padjust, "-fc", opt$foldchange,".xls")
	if (! is.null(opt$padjust)){resSig_out <- paste0(opt$outputdir,"/",compare[i,1],".vs.",compare[i,2],".DEGseq.padj-", opt$padjust, "-fc-", opt$foldchange, ".xls")}
	write.table(resSig,resSig_out, quote = F, row.names = F, sep = "\t", na = "")
	write.table(resSig[,1],file=paste0(opt$outputdir,"/",compare[i,1],".vs.",compare[i,2],'.',opt$padjust,".DE.list",sep=''), quote = F, row.names = F,col.names=F, sep = "\t", na = "")
	##volcano plot here##
	res$type <-ifelse(res$Significant=="no","nosig", ifelse(res$Regulate=="down","down","up"))
	res$type <- factor(res$type,levels=c("nosig","up","down"))
	res$Log2FC <- res[,grep("^log2",colnames(res))]
	res$LogPval <- -log10(res[,Pval])
	volcano<-ggplot(res,aes(x=Log2FC,y=LogPval))+
	geom_point(aes(color=type,shape=type))+
	scale_color_manual(values=c("up"="#FF0000","down"="#009933","nosig"="grey"))+
	scale_shape_manual(values=c("up"=18,"down"=15,"nosig"=19))+
	labs(x=paste0("Log2(",fc,")"),y=paste0("-Log10(",Pval,")"),title=paste0(compare[i,1],".vs.",compare[i,2],".volcano\n"))+
	#scale_y_continuous(expand=c(0,0),limits = c(0,ceiling(max(res$LogPval))+1))+
	scale_y_continuous(expand=c(0,0))+
	#scale_x_continuous(expand=c(0,0),limits = c(floor(min(res$Log2FC)),ceiling(max(res$Log2FC))))+
	scale_x_continuous(expand=c(0,0))+
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
	ggsave(paste0(opt$outputdir,"/",compare[i,1],".vs.",compare[i,2],".diffexp.volcano.pdf"), height=9, width=8,plot = volcano)
	#save.image("result.RData")
	## MA plot ##
	if (! is.null(expression)){
	res$express <- apply(res[,2:(nrow(colData)+1)], 1, mean)
	MA<-ggplot(res,aes(x=log10(express),y=Log2FC))+
	geom_point(aes(color=type,shape=type))+
	scale_color_manual(values=c("up"="#FF0000","down"="#009933","nosig"="grey"))+
	scale_shape_manual(values=c("up"=18,"down"=15,"nosig"=19))+
	labs(x=paste0("Log10(",opt$exp_type,")"),y=paste0("Log2(",fc,")"),title=paste0(compare[i,1],".vs.",compare[i,2],".MA\n"))+ 
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
		ggsave(paste0(opt$outputdir,"/",compare[i,1],".vs.",compare[i,2],".diffexp.scatter.pdf"), height=9, width=8,plot = MA)
    }
}

print(Sys.time()-timeStart)
