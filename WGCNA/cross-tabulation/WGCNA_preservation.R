
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


dat1 = read.csv("co.xls",sep="\t", header=T)
dim(dat1)
names(dat1)



indexB=c(2:37)
indexnts=c(39:74)
dat=cbind(dat1[,indexB],dat1[,indexnts])
datExpr=t(dat)
indexB=c(1:36)
indexnts=c(36:72)

# Number of data sets that we work with
nSets = 2;
# Object that will contain the expression data
multiExpr = list();
multiExpr[[1]] = list(data = data.frame(datExpr[indexB, ]));
multiExpr[[2]] = list(data = data.frame(datExpr[indexnts, ]));
# Names for the two sets
setLabels = c("B", "nts");
# Important: components of multiExpr must carry identificating names
names(multiExpr) = setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExpr, lapply, dim)


colo=read.csv("modualcolor.csv",header=T)
B=as.data.frame(colo$B)
nts=colo$nts
#x = load("HumanChimp-OldhamAnalysis-colorHuman-colorChimp-inNetwork.RData")
# Create an object (list) holding the module labels for each set:
colorList = list(B, nts)
# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels;



system.time( {
	mp = modulePreservation(multiExpr, colorList,
				referenceNetworks = c(1:2),
				loadPermutedStatistics = FALSE,
				nPermutations = 10,
				verbose = 3)
} )

# Save the results
save(mp, file = "B_nts.RData");



library(impute)
# Impute missing data and calculate eigengenes
impExpr = list();
for (set in 1:nSets){
	impExpr[[set]] = list(data = t(impute.knn(t(multiExpr[[set]]$data))$data));
}

eigengenes = list();
for (set in 1:nSets){
	eigengenes[[set]] = multiSetMEs(impExpr, universalColors = colorList[[set]], excludeGrey = TRUE);
	for (ss in 1:nSets){
		rownames(eigengenes[[set]][[ss]]$data) = rownames(multiExpr[[ss]]$data);
	}
}



ref = 1;
dendrograms = list();
for (set in 2:nSets)
{
adj = abs(cor(multiExpr[[set]]$data[, inNetwork], use = "p"))^9;
dtom = TOMdist(adj);
dendrograms[[set]] = flashClust(as.dist(dtom), method = "a");
}
# Get eigengenes
mes = list()
for (set in 1:nSets)
{
mes[[set]] = moduleEigengenes(multiExpr[[set]]$data, colorList[[ref]])$eigengenes
}

r_sample_logical(32977, name = "Var")
table(r_sample_logical(1000))
c("B", "W")[r_sample_logical(10)]

library(wakefield)
inNetwork=r_sample_logical(32976, name = "Var")

# Calculate the contingency table and p-values
overlap = overlapTable(B[inNetwork], nts[inNetwork]);
# The numMat will encode color. We use -log of the p value.
numMat = -log10(overlap$pTable);
numMat[numMat >50] = 50;
# Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of
# counts and corresponding p-values.
textMat = paste(overlap$countTable, "\n", signif(overlap$pTable, 2));
dim(textMat) = dim(numMat)
# Additional information for the plot. These will be used shortly.
xLabels = paste("M", sort(unique(nts)));
yLabels = paste("M", sort(unique(B)));
xSymbols = paste(sort(unique(nts)), ": ", table(nts[inNetwork]), sep = "")
ySymbols = paste(sort(unique(B)), ": ", table(B[inNetwork]), sep = "")



# Open a graphical window. If plotting into a file, skip the next line.
sizeGrWindow(10,10); fp = FALSE;
pdf("b_NTS-motivationFigure-dendrosAndTable.pdf", w = 10, h = 10);
layout(matrix(c(1,2,5, 3,4,5), 3, 2),
heights = c(3, 1, 5.5), widths = c(1, 1));
#layout.show(5);
par(mgp = c(3, 1, 0));
plotDendroAndColors(dendrograms[[1]],
		cbind(B[inNetwork], nts[inNetwork]),
		c("B modules", "nts modules"),
		setLayout = FALSE,
		marAll = c(1, 6, 2.7, 0.2),
		addGuide = FALSE,
		main = "B gene dendrogram\nand module colors", cex.main = 1.2,
		dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7)#abHeight = 0.95
par(mgp = c(3, 1, 0));
plotDendroAndColors(dendrograms[[2]],
			cbind(nts[inNetwork],B[inNetwork]),
			c("nts modules","B modules"),
			setLayout = FALSE,
			marAll = c(1, 6, 2.7, 0.2),
			addGuide = FALSE,
			main = "nts gene dendrogram\nand module colors", cex.main = 1.2,
			dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7 )#abHeight = 0.95
# Plot the overlap table
#sizeGrWindow(20,15); fp = FALSE;
#pdf("B_nts-motivationFigure-dendrosAndTable.pdf", w = 10, h = 10)
fcex = 1.00;
pcex = 1.0
fcexl = 1.00;
pcexl = 1.00;
par(mar = c(11, 11, 2, 1.0));
labeledHeatmap(Matrix = numMat,
		xLabels = xLabels, yLabels = yLabels,
		xSymbols = xSymbols, ySymbols = ySymbols,
		colorLabels = TRUE,
		colors = greenWhiteRed(100)[50:100],
		textMatrix = textMat, cex.text = if (fp) fcex else pcex, setStdMargins = FALSE,
		cex.lab = if (fp) fcexl else pcexl,
		#xColorWidth = 0.08,
		#xColorWidth = 2 * 8,#strheight("M")
		#yColorWidth = 2 * 8,
		main = "B modules (rows) vs. nts modules (columns)", cex.main = 1.2)#,xlab="B modules",ylab="nts modules")#,xLabels=,yLabels="nts modules")
mtext(side=1, text="nts modules", line=9)
mtext(side=2, text="B modules", line=9)

dev.off()


























