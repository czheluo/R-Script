# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("tidyverse")

#import packages
require(ggplot2)
require(ggrepel)
require(tidyverse)
require(ggsci)
#input file 
options = commandArgs(trailingOnly = T)
file_name = options[1]
out_name = options[2]

cirosclass<-read.table(file_name,header=T)
cirosclass$type <- factor(cirosclass$type, levels = cirosclass$type)

p<-ggplot(cirosclass, aes(x = "" , y = num, fill = type)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set2") +
  #geom_label_repel(data = cirosclass,
   #                aes(y = num-5, label = paste0(num)),
   #                size = 4.5, nudge_x = 0.7, show.legend = FALSE) +
  geom_text_repel(data = cirosclass,
                   aes(y = num-5, label = paste0(num)),
    nudge_x = .7,
    box.padding = 0.5,
    nudge_y = 1,
    point.padding = 0, # additional padding around each point
    min.segment.length = 0, # draw all line segments
    #segment.curvature = -0.1,
    segment.ncp = 1,
    #segment.angle = 20,
    #min.segment.length = -0.1,
   # arrow = arrow(length = unit(0.02, "npc"))
    )+
  guides(fill = guide_legend(title = " ")) +
  theme_void()+
  labs(title = 'cirRNA_class',x = '', y = '')+
  theme(plot.title = element_text(size = 12, color = 'black', hjust = 0.5),
	text=element_text(family="Times", size=12,color = 'black'))+
  scale_fill_npg()

png(paste('cirRNAclass.png',sep=""), width = 1000,height=800)
p
dev.off()

pdf(paste('cirRNAclass.pdf',sep=""), width = 5,height=5)
p
dev.off()
