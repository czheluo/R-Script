rm(list = ls())
#import packages
library(ggplot2)
library(ggsci)
library(extrafont)

#input file 
options = commandArgs(trailingOnly = T)
file_name = options[1]
data = read.table(file_name, sep = '\t', header = T, quote = "",stringsAsFactors=F)

#font_import(paths = "~/.local/share/fonts/",prompt = F)

#data<-read.table("kegg_statistics.xls",header=T,sep="\t",stringsAsFactors=F)
data$sum=data[,3]+data[,4]

data$second_category <- factor(data$second_category, levels = data$second_category)
data$first_category <- factor(data$first_category, levels = unique(data$first_category))
pdf('Histogram_of_Kegg.pdf', width = 10)
ggplot(data = data, aes(x=second_category, y=sum, fill=first_category))+
  geom_bar( stat="identity",position = position_dodge(1.3), width = 0.8)+
  labs(title = 'Histogram of KEGG',x = 'KEGG pathways', y = 'Number of genes') +
  scale_y_continuous(expand = c(0, 0))+
  #scale_fill_brewer(palette="Set3")+
  #geom_text(aes(label = percent, vjust = -0.3, hjust = 0.6))+
  #scale_y_continuous(limits = c(0, floor(max(plot_data$percent)*1.2)), expand = c(0,0))+
  theme_classic()+
  theme(plot.margin=unit(rep(3,4),'lines'),
	#panel.grid.minor.y = element_blank(),
        #panel.grid.major.y = element_line(size = 0.1, color = 'black'),
        #panel.grid.major.x = element_blank(),
        #panel.background = element_rect(fill = 'white', color = 'white'),
        #plot.margin=unit(rep(3,4),'lines'),
        plot.title = element_text(size = 12, color = 'black', hjust = 0.5),
        axis.text.x = element_text(angle=70,hjust=1,vjust=1,color = 'black'),
        axis.title.y = element_text(size = 12, color = 'black', vjust = 5),
        axis.title.x = element_text(size = 12, color = 'black'),
	text=element_text(family="Times", size=12,color = 'black'),
        #axis.title.y = element_text(size = 12, color = 'black'),
	#axis.title.x = element_blank(),
	#panel.grid = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	panel.border = element_blank(),
	panel.background = element_blank(),
	axis.line = element_line(),
	legend.title = element_blank(),
        legend.position = c(0.9,0.9),
        legend.key.size = unit(0.5,'cm'))+
	scale_fill_npg()
dev.off()
        #axis.ticks = element_blank())


