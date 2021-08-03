args<-commandArgs(T)

library(ggplot2)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]


rpm = read.table(args[1], header=TRUE, check.names=0, row.names=1 )
logrpm=log(rpm)/log(2)

dens=c()
for(i in 1:dim(logrpm)[2]){
	dens=c(dens,logrpm[,i])
}
distribs = data.frame(
  values = dens,
  type = gl(n = dim(logrpm)[2], k = dim(logrpm)[1]))
 
# plot
pdf(paste(args[1],".density.pdf",sep=""),width=10,height=10)
ggplot(data = distribs, aes(x = values, fill= type)) +
  geom_density(color = NA, alpha = 0.8)+
  xlim(-10, 20)+
  scale_fill_manual(values = rcolor[1:dim(logrpm)[2]],labels = colnames(logrpm))+
  guides(fill = guide_legend(title = ""))+
  scale_y_continuous(expand = c(-0.001,0))+
  theme_classic()+
  labs(title = 'Expression density distribution',x = 'log2RPM', y = 'density')+
  theme(plot.margin=unit(rep(3,4),'lines'),
	#panel.grid.minor.y = element_blank(),
        #panel.grid.major.y = element_line(size = 0.1, color = 'black'),
        #panel.grid.major.x = element_blank(),
        #panel.background = element_rect(fill = 'white', color = 'white'),
        #plot.margin=unit(rep(3,4),'lines'),
        plot.title = element_text(size = 12, color = 'black', hjust = 0.5),
        #axis.text.x = element_text(angle=70,hjust=1,vjust=1,color = 'black'),
        #axis.title.y = element_text(size = 12, color = 'black', vjust = 5),
        #axis.title.x = element_text(size = 12, color = 'black'),
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
        legend.key.size = unit(0.5,'cm'))

dev.off()
