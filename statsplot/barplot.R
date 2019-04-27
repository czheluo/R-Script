
library(RColorBrewer)

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
rcolor<-color[sample(1:length(color),length(color))]


pq_before_convert <- parental_quantities(dosage_matrix = all_genotype,
                                         parent1 = P1,parent2 = P2, las=2)


png(paste("ALL segregation summary",".png",sep=""),width=1600, height=900)

par(mfrow=c(1,1))
col<-length(pq_before_convert)
cor<-rcolor[1:col]
bp<-barplot(pq_before_convert,#,space = 0,
        col=cor,border = NULL,ylab = "Nr.Markers",
        main = "Marker segregation summary",cex.names =1.2,cex.axis = 1, cex.lab=1.2)
mtext(paste("Total number markers:", sum(pq_before_convert)),side=c(1,1,1,1), font = 3)
text(x=bp,y=pq_before_convert,labels = round(pq_before_convert,0),pos=3,xpd=NA)
dev.off()



pdf(paste("ALL segregation summary",".pdf",sep=""),width=16, height=9)

par(mfrow=c(1,1))
col<-length(pq_before_convert)
cor<-rcolor[1:col]
bp<-barplot(pq_before_convert,#,space = 0,
        col=cor,border = NULL,ylab = "Nr.Markers",
        main = "Marker segregation summary",cex.names =1.2,cex.axis = 1, cex.lab=1.2)
mtext(paste("Total number markers:", sum(pq_before_convert)),side=c(1,1,1,1), font = 3)

text(x=bp,y=pq_before_convert,labels = round(pq_before_convert,0),pos=3,xpd=NA)

dev.off()



pq_before_convert <- parental_quantities(
  dosage_matrix = ALL_dosages_a,parent1 = P1,parent2 = P2,
  las = 2)




png(paste("segregation summary",".png",sep=""),width=1600, height=900)

par(mfrow=c(1,1))
col<-length(pq_before_convert)
cor<-rcolor[1:col]
bp<-barplot(pq_before_convert,#,space = 0,
        col=cor,border = NULL,ylab = "Nr.Markers",
        main = "Marker segregation summary",cex.names =1.2,cex.axis = 1, cex.lab=1.2)
mtext(paste("Total number markers:", sum(pq_before_convert)),side=c(1,1,1,1), font = 3)
text(x=bp,y=pq_before_convert,labels = round(pq_before_convert,0),pos=3,xpd=NA)
dev.off()



pdf(paste("segregation summary",".pdf",sep=""),width=16, height=9)

par(mfrow=c(1,1))
col<-length(pq_before_convert)
cor<-rcolor[1:col]
bp<-barplot(pq_before_convert,#,space = 0,
        col=cor,border = NULL,ylab = "Nr.Markers",
        main = "Marker segregation summary",cex.names =1.2,cex.axis = 1, cex.lab=1.2)
mtext(paste("Total number markers:", sum(pq_before_convert)),side=c(1,1,1,1), font = 3)

text(x=bp,y=pq_before_convert,labels = round(pq_before_convert,0),pos=3,xpd=NA)

dev.off()



for(i in 1:4){
  text(i-0.2,pq_before_convert[i]+200,t(pq_before_convert)[i])
}


pq_before_convert <- parental_quantities(
  dosage_matrix = ALL_dosages_a,parent1 = P1,parent2 = P2,
  las = 2)

tmp = c(mean(1.0000001:100),mean(100.0000001:200), mean(200.0000001:300), mean(300.0000001:400))
bp = barplot(tmp, names=c("site 1", "site 2", "site 3", "site 4") )
# numbers above bars
text(x=bp, y=tmp, labels=round(tmp,0), pos=3, xpd=NA)
# numbers within bars
text(x=bp, y=tmp, labels=round(tmp,0), pos=1)




