
library(ggplot2)
library(dplyr)
library(forcats)

files<-list.files()
kegg<-files[grep("*kegg_enrichment.xls$",files)]
go<-files[grep("*enrichment.detail.xls$",files)]
#top20_KEGG_enrichment_terms_by_p-value
for (i in 1:length(kegg))
{
read<-read.table(kegg[i],sep="\t",fill=T)
read$input_number<-as.numeric(read[,5])
read$background<-as.numeric(read[,6])
read$padjudt<-as.numeric(read[,8])
read$rich_factor<-read[,5]/read[,6]
read<-read[order(read$padjudt,na.last=TRUE),]
read10<-read[1:10,]
pp<-read10 %>% mutate(name = fct_reorder(V2, padjudt)) %>% ggplot(aes(rich_factor,name))+
	geom_point(aes(size=input_number,color=padjudt))+
	scale_color_gradient(low="red",high="green")+
	labs(x="Rich factor",y="Top 10 terms")
name<-gsub("xls","padjudt.png",kegg[i])
ggsave(pp, file = name, width = 10, height = 6.75)
name1<-gsub("xls","",kegg[i])
svg(file = paste(name1,"padjudt.svg",sep=""), width = 10, height = 6.75)
pp
dev.off()
#ggsave(pp, file = paste(name1,".svg",sep=""), width = 10, height = 6.75,device="svg")
ggsave(pp, file = paste(name1,"padjudt.tiff",sep=""), width = 10, height = 6.75,device="tiff")
ggsave(pp, file = paste(name1,"padjudt.pdf",sep=""), width = 10, height = 6.75,device="pdf")
}

#top20_GO_enrichment_terms_by_p-value
#for (i in 1:length(go))
#{
#read<-read.table(go[i],sep="\t",header=T)
#ratio_instudy<-as.character(read$ratio_in_study)
#ratio_inpop<-as.character(read$ratio_in_pop)
#read$input<-sapply(ratio_instudy,function(x)unlist(strsplit(x,"/"))[1])
#read$background<-sapply(ratio_inpop,function(x)unlist(strsplit(x,"/"))[2])
#read$input_number<-as.numeric(read$input)
#read$background_number<-as.numeric(read$background)
#read$pvalue<-read[,6]
#read$rich_factor<-read$input_number/read$background_number
#read<-read[order(read$pvalue,na.last=TRUE),]
#read20<-read[1:20,]
#pp<-ggplot(read20,aes(rich_factor,description))+geom_point(aes(size=input_number,color=-log10(pvalue)))+scale_color_gradient(low="green",high="red")+labs(x="Rich factor",y="Top 20 terms")
#name<-gsub("xls","png",go[i])
#ggsave(pp, file = name, width = 10, height = 6.75)
#}
