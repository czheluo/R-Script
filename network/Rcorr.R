#Formatting the correlation matrix in 4 column tables
#The following custom function is used:
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}
data<-read.table('1.txt',header=T )
tma <- t(data[,-1])
colnames(tma)<-data[,1]
library(Hmisc)
gene_cor<-rcorr(tma,type="spearman")
mat <- flattenCorrMatrix(gene_cor$r,gene_cor$P)
write.csv(mat,'1.csv')
# Viz
#library(PerformanceAnalytics)
#chart.Correlation(mydata, histogram=TRUE, pch=19)


## Network plot from expression data in R using igraph ##

library(igraph)

# Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph.adjacency(
  as.matrix(as.dist(cor(t(data[,-1]), method="pearson"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
#g <- graph.adjacency(
#  as.matrix(dist(estrogenMainEffects, method="euclidean")),
#  mode="undirected",
#  weighted=TRUE,
#  diag=FALSE
#)

# Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"

# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"

# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)

# Change arrow size
# For directed graphs only
#E(g)$arrow.size <- 1.0

# Remove edges below absolute Pearson correlation 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.5)])

# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name

# Change shape of graph vertices
V(g)$shape <- "sphere"

# Change colour of graph vertices
V(g)$color <- "skyblue"

# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"

# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(data[,-1], 1, mean)) + 1.0) * 10

# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")

# Plot the tree object
plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph"
)


## step3

mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

par(mfrow=c(1,2))
plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph"
)

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph"
)

## network art

### https://www.data-to-viz.com/graph/edge_bundling.html


