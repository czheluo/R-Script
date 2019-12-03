devtools::install_github("mattflor/chorddiag")
library(chorddiag)
library("RColorBrewer")
display.brewer.pal(Dark)

setwd("H:\\PROJECT\\RNA\\chord_diagram\\case6\\baotong\\GO")

ko<-read.table("GO.UP",header = T)
ko<-read.table("GO.DOWN",header = T)
ko<-read.csv("go.down.csv",header = T)
term<-unique(ko$term)
gene<-unique(ko$gene)
#chord<-matrix(0,length(term),length(gene))
na <-matrix(0,length(gene),length(term))

for (i in 1:length(term)) {
  ann<-table( ko[which(ko[,1] %in% unique(ko[,1])[i]),2])
  na[,i]<-as.matrix(ann)
  #chord[which(ko[which(ko[,2] %in% gene[i]),1] %in% term),i]=5
}
chorddiag(na, type = "bipartite", showTicks = F,
          palette = "Set3", palette2 = "BrBG", 
          showGroupnames = TRUE, 
          groupnameFontsize = 20, groupnamePadding = 10, margin = 200)
chordd(na, type = "bipartite", showTicks = F,palette = "BrBG", palette2 = "BrBG", showGroupnames = TRUE, 
       groupNames = NULL, groupColors = NULL, 
       groupnameFontsize = 14, groupnamePadding = 10, margin = 200)



setwd("H:\\PROJECT\\RNA\\chord_diagram\\case6\\baotong\\GO/")

ko<-read.table("KEGG.DOWN",header = T,sep = "\t")
ko<-read.table("KEGG.UP",header = T,sep = "\t")
ko<-read.csv("go.down.csv",header = T)
term<-unique(ko$term)
gene<-unique(ko$gene)
#chord<-matrix(0,length(term),length(gene))
na <-matrix(0,length(gene),length(term))

for (i in 1:length(term)) {
  ann<-table( ko[which(ko[,1] %in% unique(ko[,1])[i]),2])
  na[,i]<-as.matrix(ann)
  #chord[which(ko[which(ko[,2] %in% gene[i]),1] %in% term),i]=5
}

colnames(na)<-term
rownames(na)<-gene
chorddiag(na, type = "bipartite", showTicks = F,
          palette = "Set3", palette2 = "BrBG", 
          showGroupnames = TRUE, 
          groupnameFontsize = 20, groupnamePadding = 10, margin = 200)
chordd(na, type = "bipartite", showTicks = F,palette = "BrBG", palette2 = "BrBG", showGroupnames = TRUE, 
       groupNames = NULL, groupColors = NULL, 
       groupnameFontsize = 14, groupnamePadding = 10, margin = 200)

chorddiag(na, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding = 10, margin = 90)

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
          palette = "Set3", palette2 = "BrBG", 
          showGroupnames = TRUE, 
          groupnameFontsize = 20, groupnamePadding = 10, margin = 200)
chordd(na, type = "bipartite", showTicks = F,palette = "BrBG", palette2 = "BrBG", showGroupnames = TRUE, 
       groupNames = NULL, groupColors = NULL, 
       groupnameFontsize = 14, groupnamePadding = 10, margin = 150)

for (i in 1:length(gene)) {
  ann<-table( ko[which(ko[,1] %in% unique(ko[,1])[3]),2])
  chord[which(ko[which(ko[,2] %in% gene[i]),1] %in% term),i]=5
}

rownames(chord)<-gene
colnames(chord)<-term
chorddiag(chord, type = "bipartite", showTicks = F, 
          groupnameFontsize = 14, groupnamePadding = 10, margin = 100)
chordd(chord, type = "bipartite", showTicks = F, 
       groupnameFontsize = 14, groupnamePadding = 10, margin = 90)
