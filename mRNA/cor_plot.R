rm(list=ls())
#import package
library(Hmisc)
library(ggplot2)
#library(ggprism)
library(optparse)
#input table file
option_list <-  list(
  make_option(c("-c","--data1"),type="character",help="file with log2fc"),
  make_option(c("-s","--data2"),type="character",help="file with log2fc"),
  make_option(c("-o","--outpwd"),type="character",help="out pwd")
)
opt <- parse_args(OptionParser(option_list=option_list))

fname1 = opt$data1
fname2 = opt$data2
outfile = opt$outpwd
#################################################
df1 = read.table(fname1,sep = '\t', header = T,quote = "" )
df2 = read.table(fname2,sep = '\t', header = T, quote = "" )
#extract columns
FC_data = cbind(df1[,c("log2fc")], df2[,c("log2fc")])
#columns rename
colnames(FC_data) = c('df1', 'df2')

#calculate correlation and p-value
CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor) 
  data.frame(row = rownames(cor)[row(cor)[ut]] ,
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] )
}
res <- rcorr(as.matrix(FC_data))
cor_value = round(CorMatrix (res$r, res$P) [,'cor'], 2)
p_value = CorMatrix (res$r, res$P) [,'p']
if (p_value < 1e-200) {
  title_p = 'p < 1.0e-200'
}else{title_p = format(p_value, digits = 2, scientific = T)
  }
#draw plot
data = as.data.frame(FC_data)


#data$change <- as.factor(ifelse(abs(data$df1) > 1 & abs(data$df2) < 1,
#                                ifelse(abs(data$df1) < 1 &  abs(data$df2) > 1, 'UP','DOWN'),'NOT'))

data$change <- as.factor(ifelse(abs(data$df1) > 1 & abs(data$df2) < 1,'UP',
                                ifelse(abs(data$df1) < 1 &  abs(data$df2) < 1, 'DOWN','NOT')))



p = ggplot(data = data,aes(x = df1,y = df2, color = change))+
  geom_point(size = 2)+
  labs(title = paste0('Cor=', cor_value, ',' , title_p))+
      #x =  paste(strsplit(fname1, '.', fixed = T)[[1]][1], '[log2FC]'),
      #y = paste(strsplit(fname2, '.', fixed = T)[[1]][1], '[log2FC]'))+
  scale_color_manual(values = c("blue", "grey", "black"), limits = c("UP", "DOWN", "NOT"))  +
  geom_hline(yintercept = c(-1,1), linetype='dashed')+
  geom_vline(xintercept = c(-1,1), linetype='dashed')+
  theme(panel.background = element_rect(fill = 'white',color = 'black',size = 1),
        plot.title =element_text(colour = "black", size = 15, hjust = 0.5),
        axis.text.x = element_text(size = 10, color='black'),
        axis.text.y = element_text(size = 10, color='black'),
        axis.title.x = element_text(size = 15, color = 'black'),
        axis.title.y = element_text(size = 15, color = 'black'),
        plot.margin=unit(rep(3,4),'lines'),
        legend.position = "none",
	axis.ticks.length = unit(-1, "mm"))
	#annotation_ticks(outside = F))
ggsave(p,file=file.path(outfile,"Cor_plot.pdf"),dpi=600)
        
       
