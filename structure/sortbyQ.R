

setwd("H:\\PROJECT\\liqiang\\Result\\ALL_sample\\Structure")

###READ ALL THE BEST RESULT FROM STRUCTURE

best<-read.table(paste("pop.",15,".xls",sep=""), header = F)


lm<-list()


for(i in 2:15){
  if (i==2){
    
    lm[[i-1]]<-best[order(best[,i], decreasing = F),]
  
  }else{
    lm[[i-1]]<-lm[[i-2]][order(lm[[i-2]][,i], decreasing = F),]
    }
    
}


best<-read.csv("bestpop15.csv",header = F)

lm<-best[with(order(best[,c(2:16)],decreasing=TRUE)),]

write.csv(lm,file="lm.pop15.csv",quote = F,row.names = F)

id<-lm[,1]

ldata<-lm[[14]]

write.csv(ldata,file = "ldata.csv",quote = F,row.names = T)

barplot(as.matrix(t(best[,-1])),axisnames=F,width = 2,space = 0,
        col=rcolor[1:col-1],
        border=NA,las=3)

library(RColorBrewer)

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
rcolor<-color[sample(1:length(color),length(color))]

png(paste("tobest.structure",".png",sep=""),width=1600, height=900)
pdf(paste("tobest.structure",".pdf",sep=""),width=16, height=9)

par(mfrow = c(15,1),
    #oma = c(1,4,1,0) + 0.1,
    mar = c(0,4,0.2,0.2) + 0.1)


for(i in 2:15) {
  if(i<15){
    
    data<-read.table(paste("pop.",i,".xls",sep=""), header = F)
    
    col<-length(colnames(data))
    
    lmdata<-data[which(id %in% data[,1]),-1]
    
    barplot(t(as.matrix(lmdata)),axisnames=F,width = 2,space = 0,
            col=rcolor[1:col-1],
            border=NA,las=3)
    
    title(ylab = paste("k=",i,sep = ""),line=2.4, cex.lab=1.2)
    
  } else {
    
    #data<-read.table(paste("pop.",i,".xls",sep=""), header = F)
    #col<-length(colnames(data))
    #lmdata<-data[match(id,data[,1]),]
    
    barplot(t(as.matrix(best[,-1])),names.arg = id,axisnames=T,width = 2,space = 0,
            col=rcolor[1:col-1],
            border=NA,las=3)
    
    title(ylab = paste("k=",i,sep = ""),line=2.4, cex.lab=1.2)
    
  }
}


dev.off()



