
options(warn=-100)
input_matrix<-"fat_diffcirc_intron.txt"
logNorm<-10
width<-12
height<-16
clsize<-1
rlsize<-0
mlow<-12
mright<-12
cltype<-"row" #### both row column none

library(cluster)
library(gplots)
library(Biobase)

data = read.table(input_matrix, header=T, check.names=F, sep="\t",stringsAsFactors=F)
rownames(data) = data[,1] # set rownames to gene identifiers
data = data[,2:length(data[1,])] # remove the gene column since its now the rowname value
data = as.matrix(data) # convert to matrix
cdata<-colnames(data)
myheatcol = redgreen(75)[75:1]
#colnames(data)
if(logNorm!=0){
data = log(data+1,base=logNorm)
centered_data = t(scale(t(data), scale=T)) # center rows, mean substracted
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
final_data<-centered_data
}
#colnames(final_data)
if(logNorm==0){
hc_genes = agnes(data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(data, method="spearman")), method="complete") # cluster conditions
final_data<-data
}
if(cltype=="both"){Rowv=as.dendrogram(hc_genes);Colv=as.dendrogram(hc_samples)}
if(cltype=="row"){Rowv=as.dendrogram(hc_genes);Colv=NA}
if(cltype=="column"){Rowv=NV;Colv=as.dendrogram(hc_samples)}
if(cltype=="none"){Rowv=NA;Colv=NA}

#gene_partition_assignments <- cutree(as.hclust(hc_genes), k=subclustNum);
#partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
#gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")

### cexRow cexCol
if(rlsize==0){ 
NR<-nrow(final_data)
if(NR<=100){rlsize = 1.2/log10(NR);}
if(NR>100 && NR<=300){rlsize = 1/log10(NR)}
if(NR>300 && NR<=500){rlsize = 0.5/log10(NR)}
if(NR>500 && NR<=1000){rlsize = 0.3/log10(NR)}
if(NR>1000){rlsize = 0.3/log10(NR);height=24}
}

####### output the odered matrix after clustered
if(cltype=="both"){
	order_mat<-data[hc_genes$order[nrow(data):1],hc_samples$order]
	write.table(order_mat,paste("heatmap",".pdf.orderedmat",sep=""),sep="	",col.names=NA,row.names=T,quote=F)
	}
if(cltype=="row"){
	order_mat<-data[hc_genes$order[nrow(data):1],]
	write.table(order_mat,paste("heatmap",".pdf.orderedmat",sep=""),sep="	",col.names=NA,row.names=T,quote=F)
	 }

colnames(final_data)
if(1==1){subclustNum<-ifelse(nrow(order_mat)>200, 10, 5)}
####### heatmap-plot
heatmap_filename<-paste("heatmap",".pdf",sep="")
pdf(file=heatmap_filename, width=width,height=height, paper="special")

if(1==0){
	heatmap.2(final_data, dendrogram=cltype,Rowv=Rowv,Colv=Colv,col=myheatcol, scale="none", density.info="none", trace="none",cexCol=clsize, cexRow=rlsize,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(mlow,mright),labRow = F)
}
if(1==1){
	gene_partition_assignments <- cutree(as.hclust(hc_genes), k=subclustNum);
	partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
	gene_colors = partition_colors[gene_partition_assignments]
	heatmap.2(final_data, dendrogram=cltype,Rowv=Rowv,Colv=Colv,col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=clsize, cexRow=rlsize,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(mlow,mright),labRow = F)
}
dev.off()

####### plot the subclusters trendlines
subcluster_out = paste("subclusters_",subclustNum,sep="")
dir.create(subcluster_out)
gene_names = rownames(final_data)
num_cols = length(final_data[1,])
for (i in 1:subclustNum) {
	partition_i = (gene_partition_assignments == i)
	partition_centered_data = final_data[partition_i,]
	# if the partition involves only one row, then it returns a vector instead of a table
	if (sum(partition_i) == 1){
	dim(partition_centered_data) = c(1,num_cols)
	colnames(partition_centered_data) = colnames(final_data)
	rownames(partition_centered_data) = gene_names[partition_i]
	}
	outfile = paste(subcluster_out, "/subcluster_", i, sep='')
	write.table(partition_centered_data, file=outfile, quote=F, col.names=NA, row.names=T,sep="	")
}
tmp = as.data.frame(list.files(subcluster_out))
colnames(tmp)=c("file")
tmp$num = as.numeric(gsub("subcluster_","",tmp$file))
tmp = tmp[order(tmp$num),]
files <- as.character(tmp$file)
ncols<-ceiling(length(files)/2)
pdf(file=paste("heatmap","_trendlines_for_",subclustNum,"_subclusters.pdf",sep=""))
par(mfrow=c(2,2))
par(cex=0.6)
par(mar=c(7,4,4,2))
for (i in 1:length(files)) {
    data = read.delim(paste(subcluster_out,files[i],sep="/"), header=T, row.names=1)
    ymin = min(data); ymax = max(data);
    plot_label = paste(files[i], ', ', length(data[,1]), " genes", sep='')
    plot(as.numeric(data[1,]), type='l', ylim=c(ymin,ymax), main=plot_label, col='lightgray', xaxt='n', xlab='', ylab='centered log2(fpkm+1)')
    axis(side=1, at=1:length(data[1,]), labels=as.character(cdata), las=2)
    for(r in 2:length(data[,1])) {
        points(as.numeric(data[r,]), type='l', col='lightgray')
    }
    points(as.numeric(colMeans(data)), type='o', col='blue')
}
dev.off()

