
options(warn=-100)
library("VennDiagram")
fillColor<-c("dodgerblue", "goldenrod1")
files<-list.files("list/")
#files<-unlist(strsplit("PF_known_exp.list,PFI_known_exp.list",",",fix=T))
#Lables<-unlist(strsplit("0,0",",",fix=T))
setwd("list/")
InputList<-list()
for(i in 1:(length(files)-1)){
    genes<-scan(file=files[i],what=character())
    InputList[[1]]<-genes
    #InputList[[1]]<-read.table(files[i],header = F)
  for (j in 2:length(files)) {
    genes<-scan(file=files[j],what=character())
    InputList[[2]]<-genes
    #InputList[[2]]<-read.table(files[j],header = F)
    name1 <- strsplit(files[i],split = "[known]")[[1]][1]
    name2<- strsplit(files[j],split = "[known]")[[1]][1]
    name <- c(name1,name2)
    outname<-paste(name1,name2,sep="vs_")
    names(InputList)<-name
    pdf(file = paste(outname,"known","pdf",sep="."),width=10,height=10)
    venn.plot<-venn.diagram(InputList,filename = NULL,col = "black",fill = fillColor,alpha = 0.50,cat.cex = 1,cat.fontface = "bold",margin = 0.15,cex=2,scale=TRUE)
    grid.draw(venn.plot)
    dev.off()
  }
}
