setwd("H:\\PROJECT\\RNA\\stats")

library("ggplot2")
df<-read.csv(file = "C_6w_VD.circ.stat.csv",header = T)
df$chromosome<-factor(df$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))
png(paste("C_6w_VD.chr.distribution", ".png", sep = ""), width = 1000, height = 800)
ggplot(data=df, aes(x=chromosome, y=num,fill=chromosome)) +
  geom_bar(stat="identity", width=0.5)+xlab("chromosome")+ylab("circRNA_number")

dev.off()  

pdf("C_6w_VD.chr.distribution.pdf", width = 10, height = 10)
df$chromosome<-factor(df$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))
#png(paste("C_6w_VD", ".png", sep = ""), width = 1000, height = 800)
ggplot(data=df, aes(x=chromosome, y=num,fill=chromosome)) +
  geom_bar(stat="identity", width=0.5)+xlab("chromosome")+ylab("circRNA_number")
dev.off()


df<-read.csv(file = "C_6w.circ.stat.csv",header = T)
df$chromosome<-factor(df$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))
png(paste("C_6w.chr.distribution", ".png", sep = ""), width = 1000, height = 800)
ggplot(data=df, aes(x=chromosome, y=num,fill=chromosome)) +
  geom_bar(stat="identity", width=0.5)+xlab("chromosome")+ylab("circRNA_number")

dev.off()  

pdf("C_6w.chr.distribution.pdf", width = 10, height = 10)
#df$chromosome<-factor(df$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))
#png(paste("C_6w_VD", ".png", sep = ""), width = 1000, height = 800)
ggplot(data=df, aes(x=chromosome, y=num,fill=chromosome)) +
  geom_bar(stat="identity", width=0.5)+xlab("chromosome")+ylab("circRNA_number")
dev.off()



df<-read.csv(file = "D_6w.circ.stat.csv",header = T)
df$chromosome<-factor(df$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))
png(paste("D_6w.chr.distribution", ".png", sep = ""), width = 1000, height = 800)
ggplot(data=df, aes(x=chromosome, y=num,fill=chromosome)) +
  geom_bar(stat="identity", width=0.5)+xlab("chromosome")+ylab("circRNA_number")

dev.off()  

pdf("D_6w.chr.distribution.pdf", width = 10, height = 10)
#df$chromosome<-factor(df$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))
#png(paste("C_6w_VD", ".png", sep = ""), width = 1000, height = 800)
ggplot(data=df, aes(x=chromosome, y=num,fill=chromosome)) +
  geom_bar(stat="identity", width=0.5)+xlab("chromosome")+ylab("circRNA_number")
dev.off()


df<-read.csv(file = "D_6w_VD.circ.stat.csv",header = T)
df$chromosome<-factor(df$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))
png(paste("D_6w_VD.chr.distribution", ".png", sep = ""), width = 1000, height = 800)
ggplot(data=df, aes(x=chromosome, y=num,fill=chromosome)) +
  geom_bar(stat="identity", width=0.5)+xlab("chromosome")+ylab("circRNA_number")

dev.off()  

pdf("D_6w_VD.chr.distribution.pdf", width = 10, height = 10)
#df$chromosome<-factor(df$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY"))
#png(paste("C_6w_VD", ".png", sep = ""), width = 1000, height = 800)
ggplot(data=df, aes(x=chromosome, y=num,fill=chromosome)) +
  geom_bar(stat="identity", width=0.5)+xlab("chromosome")+ylab("circRNA_number")
dev.off()



