

setwd("H:\\PROJECT\\liqiang\\Result\\ALL_sample\\Structure")

library(RColorBrewer)

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]


rcolor<-color[sample(1:length(color),length(color))]


png(paste("tobest.structure",".png",sep=""),width=1600, height=900)


pdf(paste("tobest.structure",".pdf",sep=""),width=16, height=9)

par(mfrow = c(15,1),
    #oma = c(1,4,1,0) + 0.1,
    mar = c(0,4,0.2,0.2) + 0.1)

#

for(i in 2:15){
  if(i<16){
    
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


setwd("H:\\PROJECT\\liqiang")
data<-read.table("pop.15.xls",header = F)
lm<-order(example[,1], decreasing = F)



lm<-list()

for(i in 2:15){
  if (i==2){
    lm[[i-1]]<-data[order(data[,i], decreasing = T),]
  }else
  lm[[i-1]]<-lm[[i-2]][order(lm[[i-2]][,i], decreasing = F),]
}

lmda<-lm[[14]]




dev.off()

tbl<-read.table(paste("pop.",4,".xls",sep=""), header = F)
ord = tbl[order(tbl$V2,tbl$V3,tbl$V4,tbl$V5),]
bp = barplot(t(as.matrix(ord[-1])), 
             space = c(0.2),
             col=rainbow(4),
             xlab="Individual #", 
             ylab="Ancestry",
             border=NA)

barplot(t(as.matrix(df_q[,c(1:5)])), width = 2,space = 0,
        col=rcolor[1:col-1],
        border=NA,las=2)




dev.off() 


for(i in 2:19){
  
  da<-data[order(data[,1]),]
  
}


for(i in 1:length(data[,1])){
  da[i,]<-data[i,order(data[i,])]
}

df=data

#pick max cluster, match max to cluster
maxval <- apply(df,1,max)
matchval <- vector(length=nrow(df))
for(j in 1:nrow(df)) matchval[j] <- match(maxval[j],df[j,])

#add max and match to df
df_q <- df
df_q$maxval <- maxval
df_q$matchval <- matchval

#order dataframe ascending match and decending max
df_q <- df_q[with(df_q, order(matchval,-maxval)), ]

#remove max and match
df_q$maxval <- NULL
df_q$matchval <- NULL










pdf(paste("ALL CLUSTER4",".pdf",sep=""),width=16, height=9)

par(mfrow = c(18,1),
    #oma = c(1,4,1,0) + 0.1,
    mar = c(0,4,0,0.2) + 0.1)


for(i in 2:19){
  if(i<19){
    data<-read.table(paste("pop.",i,".xls",sep=""), header = F)
    id<-data$V1;
    col<-length(colnames(data))
    #coldata<-rainbow(20)
    data<-data[,2:col]
    barplot(t(as.matrix(data)), width = 2,space = 0,
            col=rcolor[1:col-1],
            border=NA,las=2)
    title(ylab = paste("k=",i,sep = ""),line=2, cex.lab=1.2)
  } else {
    data<-read.table(paste("pop.",i,".xls",sep=""), header = F)
    id<-data$V1;
    col<-length(colnames(data))
    #coldata<-rainbow(20)
    data<-data[,2:col]
    barplot(t(as.matrix(data)), names.arg=id,width = 2,space = 0,
            col=rcolor[1:col-1],
            border=NA,las=2)
    title(ylab = paste("k=",i,sep = ""),line=2.3, cex.lab=1.2)
    
  }
}

dev.off() 

