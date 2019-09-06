
library(ggplot2)
setwd("H:\\PROJECT\\OUYANGWEI/new/3/")
read<-read.csv("B.csv",header = F,stringsAsFactors = T)
fill<- c("#00AFBB", "#E7B800", "#FC4E07")
themenew<-theme(axis.line=element_blank(),
             axis.text = element_text(size = 16),
             #axis.text.x=element_text(angle = 90),
             axis.line.x = element_blank(),#element_line(size = 1, linetype = "solid", colour = "black"),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             axis.title.x=element_text(hjust=0.09,size=16),
             #axis.title.y=element_text(angle = 180,size=16),
             legend.position="left",
             panel.background=element_blank(),
             panel.border=element_blank(),
             #panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank(),
             plot.margin = unit(c(1, 1, 1, 1), "line"),
             legend.background = element_rect(fill = "white")

)
B <- ggplot(data = read, aes(x = V1, y = V6,group=V5,colour=factor(V5))) +       
  #geom_line(aes(group = V2),size=1) +
  themenew+
  theme(legend.key = element_rect(fill = "white", color = NA))+
  #theme(legend.text=element_text(colour=factor(V5),size=12))+
  geom_hline(aes(yintercept=2),size=1)+
  geom_text(aes(label = V3),show.legend = NA,size = 3.5,vjust=-0.4)+
  #scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(limits=c(-100,500))+
  #scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"), 
                    name="",
                    labels=c("21 nt", "22 nt", "23 nt"))+
  labs(title="Positive-strand",x ="Negative-strand",y="")

pdf("B text.pdf", width = 10, height = 7)
print(B)
dev.off()
png(paste("Btext", ".png", sep = ""), width = 1000, height = 500)
print(B)
dev.off()
read<-read.csv("A.csv",header = F,stringsAsFactors = T)
fill<- c("#00AFBB", "#E7B800", "#FC4E07")
themenew<-theme(axis.line=element_blank(),
                axis.text = element_text(size = 16),
                #axis.text.x=element_text(angle = 90),
                axis.line.x = element_blank(),#element_line(size = 1, linetype = "solid", colour = "black"),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.title.x=element_text(hjust=0.09,size=16),
                #axis.title.y=element_text(angle = 180,size=16),
                legend.position="left",
                panel.background=element_blank(),
                panel.border=element_blank(),
                #panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                plot.margin = unit(c(1, 1, 1, 1), "line"),
                legend.background = element_rect(fill = "white")
                
)

A <- ggplot(data = read, aes(x = V1, y = V6,group=V5,colour=factor(V5))) +       
  #geom_line(aes(group = V2),size=1) +
  themenew+
  theme(legend.key = element_rect(fill = "white", color = NA))+
  #theme(legend.text=element_text(colour=factor(V5),size=12))+
  geom_hline(aes(yintercept=2),size=1)+
  geom_text(aes(label = V3),show.legend = NA,size = 3.5,vjust=-0.4)+
  #scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(limits=c(-60,560))+
  #scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"), 
                     name="",
                     labels=c("21 nt", "22 nt", "23 nt"))+
  labs(title="Positive-strand",x ="Negative-strand",y="")

A
pdf("A text.pdf", width = 12, height = 8)
print(A)
dev.off()
png(paste("A text", ".png", sep = ""), width = 1000, height = 500)
print(A)
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
