
##########################################
#  https://github.com/mengl1993
#  Copyright (C) MengLuo
#  contact: czheluo@gmail.com
##########################################

#requre packages
library(ggplot2)

#color
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

# prepare dataset
prim<-read.csv("primary.result.csv",header = T)
type<-unique(data$type)
snptype<- NULL
indeltype<-NULL
for (i in 1:length(type)) {
  snptype[i]<-sum(data[which(data$type %in% type[i]),3])
  indeltype[i]<-sum(data[which(data$type %in% type[i]),4])
}
result<-cbind(snptype,indeltype)
rownames(result)<-type
write.csv(result,file = "primary.csv",quote = F,row.names = T)
# read dataset
prim<-read.csv('primary.csv',header = T)

#Function to transform data to y positions
trans <- function(x){pmin(x,40) + 0.05*pmax(x-40,0)}
#yticks <- c(0, 20, 40, 500, 1000, 1500, 2000)

yticks <- c(0, 10, 30, 100,500, 700)

#Transform the data onto the display scale
prim$mean_t <- trans(prim$num)
prim$type <- factor(prim$type, levels = c("SNP","Indel"))

png(paste("breaks", ".png", sep = ""), width = 1000, height = 800)

p <- ggplot(prim, aes(x=name, y=mean_t,group=type,fill=type)) +
  geom_bar(position="dodge", stat="identity")+
  geom_col(position="dodge") +ylab("Number of markers")+xlab("segregation patterns")+
  geom_rect(aes(xmin=0, xmax=6, ymin=30, ymax=35), fill="white") +
  scale_y_continuous(limits=c(0,NA), breaks=trans(yticks), 
                     labels=c("0","10000","25000","300000","500000","700000","750000"))+
  scale_fill_manual("Marker Type",values = rcolor[c(1:2)]) +
  geom_text(aes(x=name, y=mean_t,label=txt,group=type),
            position = position_dodge(width = 1),
            vjust = -0.5, size = 4)+
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.background = element_rect(color = "black"),
    
    #panel.border = element_rect(fill = NA),
    #panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    #strip.background = element_blank(),
    #axis.line = element_line(color="black", size = 1),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "line")
  )

print(p)
dev.off()
pdf("breaks.pdf", width = 10, height = 10)

print(p)

dev.off()
# SET THEM
theme <- theme(
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  panel.background = element_blank(),
  #panel.border = element_rect(fill = NA),
  #panel.grid.major = element_blank(),
 # panel.grid.minor = element_blank(),
  #strip.background = element_blank(),
  #axis.line = element_line(color="black", size = 1),
  axis.text.x = element_text(colour = "black"),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)
