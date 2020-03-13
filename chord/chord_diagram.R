# Libraries
# install.packages("devtools")
devtools::install_github("thomasp85/patchwork")
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")
library(htmltools)
library(RColorBrewer)

# Load dataset from github
data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/13_AdjacencyDirectedWeighted.csv", header=TRUE)

# short names
colnames(data) <- c("Africa", "East Asia", "Europe", "Latin Ame.",   "North Ame.",   "Oceania", "South Asia", "South East Asia", "Soviet Union", "West.Asia")
rownames(data) <- colnames(data)

# I need a long format
data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname)

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor <- viridis(10, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:10)]

# Base plot
chordDiagram(
  x = data_long, 
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      minor.ticks = 1, 
      major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  }
)


library(chorddiag)

students = data.frame(Math = c(50, 25, 10, 5,15),
                      Art = c(0, 25, 10, 5,0),
                      Science = c(0,25,0, 5,15),
                      PE = c(0,0,25,5,15))
# the size represent the size of vulue color was samll to big changged

students = as.matrix(students)
row.names(students) = c("Section A", "Section B", "Section C", "Section E","Section D")

chorddiag(students, type = "bipartite", showTicks = F, 
          palette = "Set3", palette2 = "Blues",
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)

chordd(students, type = "bipartite", showTicks = F, 
          #palette = "ygopp", palette2 = "rainbow",
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)

chorddiag(students, type = "bipartite", showTicks = F, 
       palette = "Set3", palette2 = "Blues",
       groupnameFontsize = 14, groupnamePadding = 10, margin = 90)
library(RColorBrewer)
palette(brewer.pal(n = 8, name = "Set2"))




library(chorddiag)

students = data.frame(Math = c(50, 25, 5, 12),
                      Art = c(10, 55, 5, 20),
                      Science = c(45,12,29, 20),
                      PE = c(24,67,27,15))

students = as.matrix(students)
row.names(students) = c("Section A", "Section B", "Section C", "Section D")

chorddiag(students, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)

setwd("I:\\MAJORBIO\\Pipeline\\R")

pro<-read.csv("protein.csv",header = T)
pro<-read.csv("kegg_down.csv",header = T)
PO<-pro[,-1]
rownames(PO)<-pro[,1]
head(PO)
po<-as.matrix(PO)
chorddiag(po, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)
chordd(po, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)



# Load package
# devtools::install_github("mattflor/chorddiag")
library(chorddiag)

# Create dummy data
m <- matrix(c(11975,  5871, 8916, 2868,
              1951, 10048, 2060, 6171,
              8010, 16145, 8090, 8045,
              1013,   990,  940, 6907),
            byrow = TRUE,
            nrow = 4, ncol = 4)

# A vector of 4 colors for 4 groups
haircolors <- c("black", "blonde", "brown", "red")
dimnames(m) <- list(have = haircolors,
                    prefer = haircolors)
groupColors <- c("#000000", "#FFDD89", "#957244", "#F26223")

# Build the chord diagram:
p <- chorddiag(m, groupColors = groupColors, groupnamePadding = 20)
p



setwd("I:\\MAJORBIO\\Pipeline\\R")
library(chorddiag)

ko<-read.csv("case2.kegg.down.csv",header = T)
ko<-read.csv("case1414.csv",header = T)

term<-unique(ko$KO)
term<-unique(ko$term)
gene<-unique(ko$gene)
chord<-matrix(0,length(term),length(gene))
na <-matrix(0,length(gene),length(term))
for (i in 1:length(term)) {
  for (j in 1:length(gene)) {
    if (ko[which(ko[,1] %in% term[i]),2] %in% gene[j]){
      na[j,i]=5
    }else{
      na[j,i]=0
    }
    
  }
  #chord[which(ko[which(ko[,1] %in% term[i]),2] %in% gene),i]=5
}
rownames(na)<-gene
colnames(na)<-term
chorddiag(na, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)

na<-as.matrix(ko[,-1])
rownames(na)<-ko$id
chorddiag(na, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)

source("I:\\MAJORBIO\\Pipeline\\R\\chordd.R")

chordd(na, type = "bipartite", showTicks = F, 
       groupnameFontsize = 14, groupnamePadding = 10, margin = 90)

for (i in 1:length(gene)) {
  #for (j in 1:length(gene)) {
  #if (ko[which(ko[,1] %in% term[i]),2] %in% gene[j]){
  #na[j,i]=5
  #}else{
  #na[j,i]=0
  #}
  
  #}
  chord[which(ko[which(ko[,2] %in% gene[i]),1] %in% term),i]=5
}

rownames(chord)<-gene
colnames(chord)<-term
chorddiag(chord, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)
chordd(chord, type = "bipartite", showTicks = F, 
       groupnameFontsize = 14, groupnamePadding = 10, margin = 90)



setwd("I:\\MAJORBIO\\Pipeline\\R\\case14/KO/")
setwd("F:\\MAJORBIO\\Pipeline\\R\\case14/KO/")
library(chorddiag)

ko<-read.csv("case2.kegg.down.csv",header = T)
ko<-read.csv("orf.GO.csv",header = T)

ko<-read.csv("euk.GO.csv",header = T)

ko<-read.csv("orf.KO.csv",header = T)

ko<-read.csv("euk.KO.csv",header = T)


term<-unique(ko$KO)
term<-unique(ko$term)
gene<-unique(ko$gene)

na<-as.matrix(t(ko[,-1]))
rownames(na)<-ko$id
colnames(na)<-ko$id
chorddiag(na, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding =2, margin = 250)


ko<-read.csv("orf.KO.csv",header = T)

ko<-read.csv("euk.KO.csv",header = T)

na<-as.matrix(t(ko[,-1]))
#rownames(na)<-ko$id
colnames(na)<-ko$id
library(circlize)
chordDiagram(na)
mat<-na

chordDiagram(mat, annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.001),
             annotationTrackHeight = convert_height(c(1, 1), "mm"))#,
             #preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(na))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df
chordDiagram(df)
setwd("F:\\MAJORBIO\\Pipeline\\R\\case14/KO/")
library(chorddiag)
library(circlize)
ko<-read.csv("orf.KO.csv",header = T)
ko<-read.csv("euk.KO.csv",header = T)

mat<-as.matrix(ko[,-1])
rownames(mat)<-ko$id
#mat<-as.matrix(t(ko[,-1]))
#colnames(mat)<-ko$id
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]
#sum(dim(mat))
col<-rcolor[c(1:sum(dim(mat)))]
col1<-rcolor[c(1:sum(dim(mat)))]

col2<-rcolor[c(1:sum(dim(mat)))]

chordDiagram(mat, grid.col =col1 ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)

png("euk.KO.png", width = 1000, height = 1000)

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()
pdf("euk.KO.pdf", width = 10, height = 10)

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()


png("orf.KO.png", width = 1000, height = 1000)

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()
pdf("orf.KO.pdf", width = 10, height = 10)

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()

setwd("F:\\MAJORBIO\\Pipeline\\R\\case14/GO/")
library(chorddiag)
library(circlize)
ko<-read.csv("orf.GO.csv",header = T)
ko<-read.csv("euk.GO.csv",header = T)

mat<-as.matrix(ko[,-1])
rownames(mat)<-ko$id
#mat<-as.matrix(t(ko[,-1]))
#colnames(mat)<-ko$id
#library(RColorBrewer)
#color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
#rcolor <- color[sample(1:length(color), length(color))]
#sum(dim(mat))
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


png("euk.GO.png", width = 1000, height = 1000)

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_METAylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()
pdf("uek.GO.pdf", width = 10, height = 10)

chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()


# save image

pdf("orf.KO.pdf", width = 10, height = 10)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]
#sum(dim(mat))
col<-rcolor[c(1:sum(dim(mat)))]
chordDiagram(mat, grid.col =col ,annotationTrack = c("grid"),#annotationTrack = c("name", "grid")
             annotationTrackHeight = c(0.03, 0.01),preAllocateTracks=2)
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()

ko<-read.csv("orf.GO.csv",header = T)
ko<-read.csv("euk.GO.csv",header = T)
