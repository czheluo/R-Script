#import packages
library(ggplot2)
library(ggsci)


data = read.table("go.txt",sep = '\t',  header = T)

data$second_category <- factor(data$second_category, levels = data$second_category)
data$first_category <- factor(data$first_category, levels = unique(data$first_category))




p<-ggplot(data = data, aes(x=second_category))+
  geom_bar(aes(y=num,fill=first_category), color="black", stat="identity", size=.1, 
           alpha=.4)+
  #geom_bar(aes(y=percent / 10,fill=first_category, color=first_category), stat="identity", size=.1, 
  #alpha=.4)+
  #geom_bar( stat="identity",position = position_dodge(1.3), width = 0.8)+
  geom_label(aes(x=second_category,y=num,label = round(percent,2)*100, vjust = -0.3, hjust = 0.6))+
  scale_y_continuous(
    limits = c(0, 70),
    # Features of the first axis
    name = "Number of genes",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./data$sum*100, name="percent"),
    expand = c(0,0)
  )+
  labs(title = 'Histogram of KEGG',x = 'KEGG pathways', y = 'Number of genes') +
  #scale_fill_brewer(palette="Set3")+
  
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
