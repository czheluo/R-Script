library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

KFex<-fkpm[,gro$seq]
KFexrow<-rowMeans(KFex)
KFexN<-KFex[which(KFexrow > 0),]

human.pca <- prcomp(as.matrix(t(KFexN)), scale = T)
human.pca.out <- as.data.frame(human.pca$x)
human.pca.out$group <- gro[, 3]


p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size = 6) + theme
p

# add label text
p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group, label = row.names(human.pca.out)))
p <- p + geom_point(size = 6) + geom_text(size = 7) + theme
p
# explation

percentage <- round(human.pca$sdev / sum(human.pca$sdev) * 100, 2)
percentage <- paste(
  colnames(human.pca.out),
  "(", paste(as.character(percentage), "%", ")", sep = "")
)

p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size = 6) + theme + xlab(percentage[1]) + ylab(percentage[2])
p

rcolor <- color[sample(1:length(color), length(color))]
pdf("KF.pca.pdf",width=15,height=10)
p <- ggplot(human.pca.out, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size = 6) +
  xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_manual(values = rcolor[c(1:242)]) + theme

p
dev.off()

data = as.matrix(KFexN)
centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
pdf("KF.cluster.pdf",width=80,height=30)
plot(hc_samples,cex=0.8)
dev.off()




YLex<-fkpm[,ylgro$seq]
YLexrow<-rowMeans(YLex)
YLexN<-YLex[which(YLexrow > 0),]

human.pcayl <- prcomp(as.matrix(t(YLexN)), scale = T)
human.pca.outyl <- as.data.frame(human.pcayl$x)
human.pca.outyl$group <- ylgro[, 3]




theme <- theme(
  axis.text = element_text(size = 16),
  axis.title = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.x = element_text(colour = "black"),
  axis.text.y = element_text(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1, 1, 1, 1), "line")
)

percentage <- round(human.pcayl$sdev / sum(human.pcayl$sdev) * 100, 2)
percentage <- paste(
  colnames(human.pca.outyl),
  "(", paste(as.character(percentage), "%", ")", sep = "")


pdf("YL.pca.pdf",width=15,height=10)
p <- ggplot(human.pca.outyl, aes(x = PC1, y = PC2, color = group))
p <- p + geom_point(size = 6) +
  xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_manual(values = rcolor[c(1:242)]) + theme

p
dev.off()

datayl = as.matrix(YLexN)
centered_datayl = t(scale(t(datayl), scale=F)) # center rows, mean substracted
#hc_genes = agnes(centered_datayl, diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_datayl, method="spearman")), method="complete") # cluster conditions


pdf("YL.cluster.pdf",width=80,height=30)
plot(hc_samples,cex=0.8)
dev.off()


