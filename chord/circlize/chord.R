
setwd("C:\\Users\\meng.luo\\Desktop")
library(chorddiag)
library(circlize)
setwd("C:\\Users/meng.luo/Desktop/")

circ<-read.csv("geneid.csv", header = T)
go<-unique(circ$GO.term)
gene<-unique(circ$gene.ID)
na <-matrix(0,length(gene),length(go))
for (i in 1:length(go)) {
  for (j in 1:length(gene)) {
    if (gene[j] %in% circ[which(circ[,1] %in% go[i]),2] ){
      na[j,i]=5
    }else{
      na[j,i]=0
    }
    
  }
  #chord[which(ko[which(ko[,1] %in% term[i]),2] %in% gene),i]=5
}
mat<-t(na)
rownames(na)<-gene
colnames(na)<-go
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]
sum(dim(na))
col<-rcolor[c(1:sum(dim(mat)))]
#col1<-rcolor[c(1:sum(dim(mat)))]

#col2<-rcolor[c(1:sum(dim(mat)))]

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)

png("orf.GO.png", width = 1000, height = 1000)

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()
pdf("orf.GO.pdf", width = 10, height = 10)

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()
