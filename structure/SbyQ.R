
setwd("H:\\PROJECT\\liqiang\\Result\\ALL_sample\\Structure")

library(RColorBrewer)

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]


rcolor<-color[sample(1:length(color),length(color))]


png(paste("tobest.structure",".png",sep=""),width=1600, height=900)


pdf(paste("tobest.structure",".pdf",sep=""),width=16, height=9)


par(mfrow = c(15,1),
    #oma = c(1,4,1,0) + 0.1,
    mar = c(0,4,0.2,0.2) + 0.1)


for(i in 2:16){
  if(i<15){

    data<-read.table(paste("pop.",i,".xls",sep=""), header = F)
    id<-data$V1;
    col<-length(colnames(data))
    #coldata<-rainbow(20)
    data<-data[,2:col]

    #pick max cluster, match max to cluster

    maxval <- apply(data,1,max)
    matchval <- vector(length=nrow(data))

    for(j in 1:nrow(data)) {
      matchval[j] <- match(maxval[j],data[j,])
    }

    #add max and match to df
    lmdata <- data
    lmdata$maxval <- maxval
    lmdata$matchval <- matchval
    #order dataframe ascending match and decending max
    lmdata <- lmdata[with(lmdata, order(matchval,-maxval)), ]
    write.csv(lmdata,file=paste("popnew.",i,".csv"),quote = F,row.names=T)
    lmdata$maxval <- NULL
    lmdata$matchval <- NULL

    #lmdata<-lmdata[lmdataid,]

    barplot(t(as.matrix(lmdata)),axisnames=F,width = 2,space = 0,
            col=rcolor[1:col-1],
            border=NA,las=3)


    title(ylab = paste("k=",i,sep = ""),line=2.4, cex.lab=1.2)

  } else {
    data<-read.table(paste("pop.",i,".xls",sep=""), header = F)
    id<-data$V1;
    col<-length(colnames(data))
    #coldata<-rainbow(20)
    data<-data[,2:col]

    #pick max cluster, match max to cluster

    maxval <- apply(data,1,max)
    matchval <- vector(length=nrow(data))

    for(j in 1:nrow(data)) {
      matchval[j] <- match(maxval[j],data[j,])
    }

    #add max and match to df
    lmdata <- data
    lmdata$maxval <- maxval
    lmdata$matchval <- matchval
    #order dataframe ascending match and decending max
    lmdata <- lmdata[with(lmdata, order(matchval,-maxval)), ]
    write.csv(lmdata,file=paste("popnew.",i,".csv"),quote = F,row.names=T)
    lmdata$maxval <- NULL
    lmdata$matchval <- NULL
    lmdataid<-as.numeric(rownames(lmdata))
    id<-id[lmdataid]

    barplot(t(as.matrix(lmda[,c(2:16)])),names.arg = id,cex.names=0.8,
            width = 2,space = 0,
            col=rcolor[1:col-1],
            border=NA,las=1)
    title(ylab = paste("k=",i,sep = ""),line=2.4, cex.lab=1.4)


  }
}

dev.off()

library(RColorBrewer)

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

rcolor<-color[sample(1:length(color),length(color))]


##chose the best group

g<-read.table("group.list")
id<-as.numeric(rownames(g[order(g[,2]),]))

#png(paste("tobest.structure",".png",sep=""),width=1600, height=900)

pdf(paste("tobest.structure",".pdf",sep=""),width=16, height=9)

par(mfrow = c(16,1),
    #oma = c(1,4,1,0) + 0.1,
    mar = c(0,4,0.2,0.2) + 0.1)

for(i in 2:16) {
  if(i<17){

    data<-read.table(paste("pop.",i,".xls",sep=""), header = F)

    col<-length(colnames(data))

    lmdata<-data[id,]

    barplot(t(as.matrix(lmdata[,-1])),axisnames=F,width = 2,space = 0,
            col=rcolor[1:col-1],
            border=NA,las=3)

    title(ylab = paste("k=",i,sep = ""),line=2.4, cex.lab=1.2)

  } else if (i==15) {

    data<-read.table(paste("pop.",i,".xls",sep=""), header = F)
    #id<-data$V1;
    col<-length(colnames(data))
    #coldata<-rainbow(20)
    data<-data[,2:col]
    #pick max cluster, match max to cluster

    maxval <- apply(data,1,max)
    matchval <- vector(length=nrow(data))

    for(j in 1:nrow(data)) {
      matchval[j] <- match(maxval[j],data[j,])
    }

    #add max and match to df
    lmdata <- data
    lmdata$maxval <- maxval
    lmdata$matchval <- matchval
    #order dataframe ascending match and decending max
    lmdata <- lmdata[with(lmdata, order(matchval,-maxval)), ]
    #write.csv(lmdata,file=paste("popnew.",i,".csv"),quote = F,row.names=T)
    lmdata$maxval <- NULL
    lmdata$matchval <- NULL
    #id<-as.numeric(rownames(lmdata))
    #id<-id[lmdataid]

    barplot(t(as.matrix(lmdata)),axisnames=F,#names.arg = id,cex.names=0.8,
            width = 2,space = 0,
            col=rcolor[1:col-1],
            border=NA,las=3)

    title(ylab = paste("k=",i,sep = ""),line=2.4, cex.lab=1.4)


  } else {
    data<-read.table(paste("pop.",i,".xls",sep=""), header = F)

    col<-length(colnames(data))

    lmdata<-data[id,]

    barplot(t(as.matrix(lmdata[,-1])),axisnames=F,width = 2,space = 0,#,names.arg = id
            col=rcolor[1:col-1],
            border=NA,las=3)

    title(ylab = paste("k=",i,sep = ""),line=2.4, cex.lab=1.2)
  }
}

dev.off()






