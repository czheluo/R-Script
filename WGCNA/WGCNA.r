args<-commandArgs(trailingOnly=TRUE)
if(length(args) < 6)
{
	stop("
Rscript WGCNA.r --args <exp.matrix> <outputdir> <annot.txt> <type> <num> <Trait>
	
options:
	exp.matrix : tpm/fpkm matrix
	outputdir  : directory name for output
	annot.txt  : gene annotation file(no header)
	type       : signed or unsigned is support
	num        : number of genes to draw TOM heatmap, less than 2000 is recommend, the bigger the slower

data: 2020-03-27
contact: meng.luo@majorbio.com
Notice: Trait file is optional , if NULL, a Sample-module correlation will be done
###annot.txt format---2 column###################
gene1	description1(could be names ...)
gene2	description2
gene3   description3
...
#################################################
###Optional : Trait format#######################
S/T      Trait1     Trait2     Trait3    ...
Sample1  0          3          0         ...
Sample2  1          0          5         ...
Sample3  2          4          1         ...
...
#################################################
")
}
###step0 basic setup
library(WGCNA)
library(flashClust)
options(warn=-1)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=20)
outdir <- args[3]
dir.create(outdir)
NetType <- args[5]
myData <- read.table(args[2], sep="\t", header=TRUE)
myData <- myData[rowMeans(myData[,-1]) > 1,]
rawExpr <- as.data.frame(t(myData[,-1]))
names(rawExpr) = colnames(myData[,-1])
colnames(rawExpr) = myData[,1]
nit <- dim(rawExpr)[2]
cv <- NULL
for(i in 1:nit){
  cv[i] <- sd(rawExpr[,i])/mean(rawExpr[,i])
}
rawExpr <- rawExpr[,cv > 0.1]

###step1 sample cluster
cat("########## cluster raw sample tree ##########")
cat("\n")
sampleTree = flashClust(dist(rawExpr), method = "average")
pdf(file = paste(outdir, "/RawSampleClustering.pdf", sep = ""), width = 12, height = 9)
par(mar = c(5,5,5,5))
plot(sampleTree, main = "Sample Clustering", sub="", xlab="")
dev.off()

###step2 filter
cat("########## remove low quality genes or samples ##########")
cat("\n")
gsg = goodSamplesGenes(rawExpr, verbose = 3)
if (!gsg$allOK)
{
if (sum(!gsg$goodGenes)>0)
	printFlush(paste("Removing genes:", paste(names(rawExpr)[!gsg$goodGenes], collapse = ", ")))
if (sum(!gsg$goodSamples)>0)
	printFlush(paste("Removing samples:", paste(rownames(rawExpr)[!gsg$goodSamples], collapse = ", ")))
	datExpr = rawExpr[gsg$goodSamples, gsg$goodGenes]
}else{
	datExpr = rawExpr
}
write.table(names(rawExpr)[!gsg$goodGenes], file=paste(outdir, "/removeGene.xls", sep = ""), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(names(rawExpr)[!gsg$goodSamples], file=paste(outdir, "/removeSample.xls", sep = ""), row.names=FALSE, col.names=FALSE, quote=FALSE)

###step3 sample cluster after filter
cat("########## cluster filtered tree ##########")
cat("\n")
sampleTree = flashClust(dist(datExpr), method = "average")
pdf(file = paste(outdir, "/SampleClustering.pdf", sep = ""), width = 12, height = 9)
par(mar = c(5,5,5,5))
plot(sampleTree, main = "Sample Clustering", sub="", xlab="")
dev.off()

###step4 read in Trait data, Trait/sample module cor
cat("########## creat Trait infor and plot trait_sample ##########")
cat("\n")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
if(length(args) == 7)
{
	Trait<-read.table(args[7],sep="\t",header=T,row.names=1)
	Trait<-Trait[rownames(datExpr),]
	Trait_colors <- numbers2colors(Trait,signed = FALSE)
}else{
	Trait<-as.data.frame(diag(x=1,nrow=nSamples))
	rownames(Trait)<-rownames(datExpr)
	colnames(Trait)<-rownames(datExpr)
	Trait_colors <- numbers2colors(Trait, colors = c("white","forestgreen"),signed = FALSE)
}

#plot sample tree + trait heatmap
pdf(file = paste(outdir, "/Sample_Trait.pdf", sep = ""), width = 15, height = 15)
par(mar = c(1,4,3,1),cex=0.6)
plotDendroAndColors(sampleTree, Trait_colors, groupLabels = colnames(Trait), cex.dendroLabels = 0.8, marAll = c(1, 4, 3, 1), cex.rowText = 0.01, main = "Sample dendrogram and trait heatmap")
dev.off()

#plot pca result
cat("########## PCA analysis ##########")
cat("\n")
pca = prcomp(datExpr,na.action=na.omit)
sampletype<-rownames(datExpr)
pdf(file = paste(outdir, "/PCA.pdf", sep = ""), width = 12, height = 12)
par(mar = c(5,5,5,8))
plot(pca$x[,c(1,2)],pch=16,col=rep(rainbow(nSamples),each=1),cex=1.5,main = "PCA map")
text(pca$x[,c(1,2)],row.names(pca$x),col="black",pos=3,cex=0.8)
legend("right",legend=sampletype,ncol = 1,xpd=T,inset = -0.15, pch=16,cex=1,col=rainbow(length(sampletype)),bty="n")
dev.off()

#plot pca_3d result
library(scatterplot3d)
pdf(file = paste(outdir, "/PCA_3D.pdf", sep = ""), width = 15, height = 15)
par(mar = c(4,4,4,4))
scatterplot3d(pca$x[,1:3], highlight.3d=F, col.axis="black",color = rep(rainbow(nSamples),each=1),cex.symbols=1.5,cex.lab=1,cex.axis=1, col.grid="lightblue", main="PCA map", pch=16)
legend("topleft",legend = row.names(pca$x) ,pch=16,cex=1,col=rainbow(nSamples), ncol = 2,bty="n")
dev.off()


###step5 pick beta value
cat("########## pick beta value ##########")
cat("\n")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = NetType) ### define the network type: signed or unsigned
pdf(file = paste(outdir, "/RawSoftPower.pdf", sep = ""), width = 9, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()
cat("which softPower: ")
softPower = 16

#print beta with ref line
pdf(file = paste(outdir, "/SoftPower.pdf", sep = ""), width = 9, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.9,col="red")
#height_line = (-sign(sft$fitIndices[,3])*sft$fitIndices[,2])[which (sft$fitIndices$Power==softPower)]
abline(h=(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])[which (sft$fitIndices$Power==softPower)],col="red",lty=2)
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
abline(h=sft$fitIndices[,5][which (sft$fitIndices$Power==softPower)],col="red",lty=2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()

###step6 scale_free plot
cat("########## plot scale free model ###########")
cat("\n")
k.datExpr=softConnectivity(datExpr,power=softPower,type = NetType)-1
pdf(file = paste(outdir, "/scale-free.pdf", sep = ""), width = 9, height = 5)
scaleFreePlot(k.datExpr, main=paste("data set, power=", softPower, sep=""), truncated=F)
dev.off()

###step7 calculate adjacency
cat("########## calculate adjacency matrix ##########")
cat("\n")
adjacency = adjacency(datExpr, power = softPower,type = NetType)
write.table(adjacency,paste(outdir,"/adjacency.xls",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)

###step8 calculate TOM
cat("########## calculate topological overlap measure ##########")
cat("\n")
#mycol = c("brown2" , "chartreuse4" , "gold" , "blue" , "blueviolet" , "hotpink4" , "brown4" , "antiquewhite4" , "aquamarine4" , "blue4" , "darkmagenta" , "darkolivegreen3" , "deeppink" , "deepskyblue" , "darkgoldenrod2" , "chocolate4" , "khaki4" , "indianred1" , "lightpink4" , "lightslateblue" , "lightsteelblue" , "mediumpurple4" , "orchid4" , "dodgerblue2" , "chocolate" , "chartreuse" , "antiquewhite3" , "aquamarine3" , "cadetblue1" , "chocolate2" , "hotpink" , "khaki" , "lightpink3" , "antiquewhite1" , "aquamarine" , "indianred3" , "deepskyblue3" , "limegreen" , "lightpink" , "springgreen4" , "lightblue2" , "lightblue4" , "springgreen3" , "seashell4" , "yellow" , "thistle3" , "tomato3" , "olivedrab3" , "turquoise3" , "yellow3" , "dodgerblue4" , "orangered" , "brown1" , "chartreuse3" , "blanchedalmond" , "hotpink3" , "brown3" , "blue3" , "darkkhaki" , "darkolivegreen2" , "darkviolet" , "deeppink4" , "darkgoldenrod1" , "chocolate3" , "khaki3" , "indianred" , "lightskyblue4" , "lightslategrey" , "mediumpurple3" , "orchid3" , "dodgerblue1" , "cadetblue4" , "antiquewhite2" , "aquamarine2" , "cadetblue" , "chocolate1" , "honeydew4" , "lightpink2" , "antiquewhite" , "indianred2" , "deepskyblue2" , "lightyellow4" , "lightgrey" , "lightblue1" , "lightblue3" , "springgreen2" , "seashell3" , "whitesmoke" , "thistle2" , "tomato2" , "olivedrab2" , "turquoise2" , "yellow2" , "dodgerblue3" , "orange4" , "bisque4" , "brown" , "chartreuse2" , "gainsboro" , "hotpink2" , "blue2" , "darkgrey" , "darkolivegreen1" , "darkturquoise" , "deeppink3" , "darkgoldenrod" , "khaki2" , "lightskyblue3" , "lightslategray" , "mediumpurple2" , "orchid2" , "dodgerblue" , "cadetblue3" , "aquamarine1" , "burlywood4" , "honeydew3" , "lightpink1" , "aliceblue" , "deepskyblue1" , "lightyellow3" , "lightgreen" , "lightblue" , "springgreen1" , "seashell2" , "wheat4" , "thistle1" , "tomato1" , "olivedrab1" , "turquoise1" , "yellow1" , "orange3" , "bisque3")

#remove bad colors
mycol = standardColors()[-c(7,17,24,27,43)]
TOM = TOMsimilarity(adjacency,TOMType = NetType)
dissTOM = 1-TOM
geneTree = flashClust(as.dist(dissTOM), method = "average")
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods,colorSeq = mycol)
MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower = softPower)
MEs = MEList$eigengenes

###step9 identify inital modules
cat("########## identify dissimilarity TOM ##########")
cat("\n")
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method = "average")
pdf(file = paste(outdir, "/RawEigengeneClustering.pdf", sep = ""), width = 7, height = 9)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#input disSimilarity
cat("choose a disSimilarity(0.25): ")
MEDissThres = 0.4

#draw ref line
pdf(file = paste(outdir, "/EigengeneClustering.pdf", sep = ""), width = 7, height = 9)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres,col="red",lty=2)
dev.off()

###step10 merge modules
cat("########## merge module ##########")
cat("\n")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file = paste(outdir, "/ModuleTree.pdf", sep = ""), width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs
write.table(paste(colnames(datExpr), moduleColors, sep = "\t"), file = paste(outdir, "/netcolor2gene.xls", sep = ""), row.names=FALSE, quote=FALSE)


###step11 plot TOM heatmap 
cat("########## network heatmap ##########")
cat("\n")
nSelect = as.integer(args[6])
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = flashClust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
plotDiss = selectTOM^7
diag(plotDiss) = NA
pdf(file = paste(outdir, "/TomHeatmap.pdf", sep = ""), width = 15, height = 15)
TOMplot(plotDiss, selectTree, selectColors, main = "Tom heatmap plot")
dev.off()


###step12 plot module-trait heatmap
cat("########## module_trait_analysis ##########")
cat("\n")
#MEs = moduleEigengenes(datExpr, moduleColors, softPower=softPower)$eigengenes
MET = orderMEs(MEs)
moduleTraitCor = cor(MEs, Trait, use = "p", method="spearman")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

pdf(file = paste(outdir, "/ModuleTrait.pdf", sep = ""), width = 16, height = 9)
par(mar = c(6, 10, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(Trait),yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,cex.lab.x = 0.8,cex.lab.y = 0.85,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()
write.table(moduleTraitCor, paste(outdir,"/module_trait_cor.xls",sep=""), sep="\t", quote=F, row.names=T, col.names=NA)
write.table(moduleTraitPvalue, paste(outdir,"/module_trait_pvalue.xls",sep=""), sep="\t", quote=F, row.names=T, col.names=NA)

###step13 plot module-eigengene exp trend ###
Pearson<-cor(datExpr,use = 'p')
Pearson_p<-corPvalueStudent(Pearson, nGenes)
write.table(Pearson,paste(outdir,"/person_corr.xls",sep=""),sep="\t",quote=F,row.names=T, col.names=NA)
write.table(Pearson_p,paste(outdir,"/person_corr_pvalue.xls",sep=""),sep="\t",quote=F,row.names=T, col.names=NA)
Corr<-TOM
Colors<-mergedColors
colnames(Corr)<-colnames(datExpr)
rownames(Corr)<-colnames(datExpr)
names(Colors)<-colnames(datExpr)

colnames(Pearson)<-colnames(datExpr)
rownames(Pearson)<-colnames(datExpr)

umc = unique(mergedColors)
lumc = length(umc)

for (i in c(1:lumc)){
    if(umc[i]== "grey"){
        next
    }
    ME=MEs[, paste("ME",umc[i], sep="")]
	pdf(file = paste(outdir, "/Module_",umc[i],"_eigengene.pdf", sep = ""),width = 12, height = 16)
    par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#	Mat = t(scale(datExpr[,Colors==umc[i]]))
#	plotMat(Mat,nrgcols=30,rlabels=rownames(Mat),clabels=colnames(Mat),rcols=umc[i], main=umc[i], cex.main=2)
	plotMat(t(scale(datExpr[,Colors==umc[i]])),nrgcols=30,rlabels=F,rcols=umc[i], main=umc[i], cex.main=2)
    par(mar=c(5, 4.2, 0, 0.7))
    barplot(ME, col=umc[i], main="", cex.main=2,ylab="Eigengene expression",xlab="Samples")
	dev.off()
}

###step14 plot eigengene heatmap #########
cat("########## eigengene heatmap ##########")
cat("\n")
pdf(file = paste(outdir, "/EigenGeneHeatmap.pdf", sep = ""), width = 6, height = 6)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,3,3,3), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

###step15 plot MM - GS scatterplot #######
cat("########## MM - GS scatterplot ##########")
cat("\n")
kME = signedKME(datExpr, mergedMEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
kME_p = as.data.frame(corPvalueStudent(as.matrix(kME), nSamples))
write.table(abs(kME), paste(outdir,"/kME.xls",sep=""), sep="\t", quote=F, row.names=T, col.names=NA)
write.table(kME_p, paste(outdir,"/kME_pvalue.xls",sep=""), sep="\t", quote=F, row.names=T, col.names=NA)
geneTraitSignificance = as.data.frame(cor(datExpr, Trait, use = "p", method="spearman"))
names(geneTraitSignificance) = paste("GS.", colnames(Trait), sep="")
#GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
#names(GSPvalue) = paste("GS.", colnames(Trait), sep="")

#MM_GS scatter plot
modNames = substring(names(MET), 3)
for (module in modNames){
	if(module == "grey"){
		next
	}
	column = match(module, modNames)  # col number of interesting modules
	moduleGenes = Colors==module
	pdf(file = paste(outdir, "/MM_GS_",module,".pdf", sep = ""),width = 12, height = 12)
	par(mfrow = c(1,1))
	verboseScatterplot(abs(kME[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]), xlab = paste("Module Membership in", module, "module"), ylab = "Gene significance", main = paste("Module membership vs. gene significance"),cex.main = 1.2, cex.lab = 1.2, pch=19,cex.axis = 1.2, col = module)
	dev.off()
}

#MS barplot
GS = abs(geneTraitSignificance)
names(GS)=colnames(Trait)

for (nGS in 1:(ncol(GS))){
	nameGS = names(GS)[nGS]
	pdf (file = paste(outdir,"/MS_",nameGS,".pdf",sep=""), width = 12, height = 12)
	plotModuleSignificance(GS[,nGS],mergedColors,ylim=c(0,0.5), main= paste(nameGS," module significance"))
	dev.off()
	MS = tapply (GS[,nGS], mergedColors, mean, na.rm=T)
	write.table(as.data.frame(MS), paste(outdir,"/MS_",nameGS,".xls",sep=""), row.names=T, col.names=NA, quote=F, sep= "\t")
}


###step16 plot gene tree and trait corralation
cat("########## gene tree and trait corralation ##########")
cat("\n")
geneTraitColor=as.data.frame(numbers2colors(geneTraitSignificance,signed=TRUE,colors = colorRampPalette(c("blue","white","red"))(100)))
names(geneTraitColor)= colnames(Trait)
pdf (file = paste(outdir, "/GeneTree_TraitCorr.pdf", sep = ""), width = 12, height = 16)
par(mar = c(3.5, 7, 2, 1))
plotDendroAndColors(geneTree, cbind(mergedColors, geneTraitColor), c("Module",colnames(Trait)),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

###step17 output whole network
cat("########## output whole network ##########")
cat("\n")
cyt = exportNetworkToCytoscape(TOM,edgeFile = paste(getwd(),"/",outdir,"/all_edges.txt",sep=""),nodeFile =paste(getwd(),"/",outdir,"/all_nodes.txt",sep=""),weighted = TRUE, threshold = 0, nodeNames = names(datExpr),nodeAttr = moduleColors)

###step18 output to cytoscape
cat("########## network visualization ##########")
cat("\n")
annot = read.table(file = args[4], sep = "\t", header = TRUE)
#TOM = 1 - dissTOM
thConnect = 0
for (i in c(1:lumc)){
	if(umc[i] == "grey"){
		next
	}
	module = umc[i]
	probes = names(datExpr)
	inModule = is.finite(match(moduleColors, module))
	modProbes = probes[inModule]
	modGenes = annot[[2]][match(modProbes, annot[[1]])]
	modTOM = TOM[inModule, inModule]
	dimnames(modTOM) = list(modProbes, modProbes)
	cyt = exportNetworkToCytoscape(modTOM,
		edgeFile = paste(outdir, "/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
		nodeFile = paste(outdir, "/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
		weighted = TRUE, threshold = thConnect, nodeNames = modProbes, altNodeNames = modGenes,
		nodeAttr = moduleColors[inModule])
}
