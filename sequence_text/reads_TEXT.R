
library(ggplot2)
theme<-theme(axis.line=element_blank(),
             axis.text = element_text(size = 16),
             #axis.text.x=element_text(angle = 90),
             axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             #axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             legend.position="none",
             panel.background=element_blank(),
             panel.border=element_blank(),
             #panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank(),
             plot.margin = unit(c(1, 1, 1, 1), "line")
)
png(paste("A", ".png", sep = ""), width = 1000, height = 800)
ggplot(data = read, aes(x = V1, y = V4)) +       
  geom_line(aes(group = V2),size=1) +theme+xlab("")+ylab("")+ 
  #geom_hline(aes(yintercept=2),size=1)+
  #geom_text(aes(label = V5),size = 3.5)+
  #scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(name ="B87_segment_A (nt)",limits=c(0,3300),
                     breaks=seq(0,3300,500))
#scale_x_discrete(name ="B87_segment_A",breaks=c(50, 502, 1014,1510,2199,2745,3260),
#labels=c("50", "502", "1014","1510","2199","2745","3260"))
dev.off()
setwd("H:\\PROJECT\\OUYANGWEI/all/3/")
read<-read.csv("22.AA.csv",header = F,stringsAsFactors = F)
example <- read.csv("example.csv",header = F,stringsAsFactors = F)
ggplot(data = example, aes(x = V1, y = V5)) +       
  #geom_line(aes(group = V3),size=1,color="red") +
  theme+xlab("")+ylab("")+ 
  #geom_hline(aes(yintercept=2),size=1)+
  geom_text(aes(label = V2), fontface="bold",size = 2.5,vjust=-0.006)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(name ="B87_segment_A (nt)",limits=c(0,3300),
                     breaks=seq(0,3300,500))
pdf("A.pdf", width = 10, height = 2)
png(paste("A", ".png", sep = ""), width = 1000, height = 100)
pdf("Atext.pdf", width = 10, height = 3)
png(paste("Atext", ".png", sep = ""), width = 1000, height = 300)

ggplot(data = read, aes(x = V1, y = V5)) +       
  geom_line(aes(group = V3),size=1,color="red") +theme+xlab("")+ylab("")+ 
  #geom_hline(aes(yintercept=2),size=1)+
  #geom_text(aes(label = V3),size = 1,vjust=-0.006)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(name ="B87_segment_A (nt)",limits=c(-400,3300),
                     breaks=seq(0,3300,500))

ggplot(data = read, aes(x = V1, y = V5)) +       
  #geom_line(aes(group = V4),size=1,color="red") +
  theme+xlab("")+ylab("")+ 
  #geom_hline(aes(yintercept=2),size=1)+
  geom_text(aes(label = V2),size = 1,vjust=-0.006)+
  scale_y_continuous(limits = c(0.0000001,0.0000001574),expand = c(0, 0))+
  scale_x_continuous(name ="B87_segment_A (nt)",limits=c(0,3300),
                     breaks=seq(0,3300,500))
dev.off()
theme<-theme(axis.line=element_blank(),
             axis.text = element_text(size = 16),
             #axis.text.x=element_blank(),
             #axis.ticks.x=element_blank(),
             axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             #axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             legend.position="none",
             panel.background=element_blank(),
             panel.border=element_blank(),
             #panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank(),
             plot.margin = unit(c(1, 1, 1, 1), "line")
)
read<-read.csv("B.csv",header = F,stringsAsFactors = F)
pdf("B.pdf", width = 10, height = 2)
png(paste("B", ".png", sep = ""), width = 1000, height = 200)
pdf("Btext.pdf", width = 10, height = 4)
png(paste("Btext", ".png", sep = ""), width = 1000, height = 400)

ggplot(data = read, aes(x = V1, y = V5)) +       
  geom_line(aes(group = V2),size=1,color="white") +theme+xlab("")+ylab("")+ 
  geom_hline(aes(yintercept=2),size=1)+
  geom_text(aes(label = V3),size = 3.5,vjust=-0.4)+
  scale_y_continuous(limits=c(1.8,4.1),expand = c(0, 0))+
  scale_x_continuous(name ="B87_segment_B (nt)",limits=c(-400,3300),
                     breaks=seq(0,3300,500))
ggplot(data = read, aes(x = V1, y = V5,colour=V4)) +       
  geom_line(aes(group = V2),size=1) +theme+xlab("")+ylab("")+ 
  geom_hline(aes(yintercept=2),size=1)+
  #geom_text(aes(label = V3),size = 3.5,vjust=-0.4)+
  scale_y_continuous(limits=c(1.8,4.1),expand = c(0, 0))+
  scale_x_continuous(name ="B87_segment_B (nt)",limits=c(0,3300),
                     breaks=seq(0,3300,500))
dev.off()
theme<-theme(axis.line=element_blank(),
             axis.text = element_text(size = 16),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             #axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             legend.position="none",
             panel.background=element_blank(),
             panel.border=element_blank(),
             #panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank(),
             plot.margin = unit(c(1, 1, 1, 1), "line")
)
