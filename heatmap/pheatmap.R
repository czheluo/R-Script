args<- commandArgs(T)
library(pheatmap)
data <- read.table(args[1],sep="\t",header=T,row.names=1,check.name=F)
df_new <- log2(data+0.01)
filename=sub("$","heatmap.pdf",args[1])
cairo_pdf(filename=filename,width = 6,height = 8)
pheatmap(df_new,cluster_cols = FALSE,show_rownames = F, fontsize_col = 14,angle_col=90,
         scale="row", color =colorRampPalette(c("blue", "white", "red"))(256))
dev.off()
