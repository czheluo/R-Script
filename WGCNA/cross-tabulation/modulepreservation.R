
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


file = bzfile("Dataset1.csv.bz2");
dat1 = read.csv(file, header=T)
dim(dat1)
names(dat1)

datExpr=data.frame(t(dat1[dat1$Brain_variant_H>0,2:39]))
indexHuman=c(19:36)
indexChimp=c(1:18)

# Number of data sets that we work with
nSets = 2;
# Object that will contain the expression data
multiExpr = list();
multiExpr[[1]] = list(data = datExpr[indexHuman, ]);
multiExpr[[2]] = list(data = datExpr[indexChimp, ]);
# Names for the two sets<u></u>
setLabels = c("Human", "Chimp");
# Important: components of multiExpr must carry identificating names
names(multiExpr) = setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExpr, lapply, dim)



x = load("HumanChimp-OldhamAnalysis-colorHuman-colorChimp-inNetwork.RData")
# Create an object (list) holding the module labels for each set:
colorList = list(colorHuman, colorChimp);
# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels;



system.time( {
	mp = modulePreservation(multiExpr, colorList,
				referenceNetworks = c(1:2),
				loadPermutedStatistics = FALSE,
				nPermutations = 200,
				verbose = 3)
} )

# Save the results
save(mp, file = "HumanChimp-HumanSpecific-modulePreservation.RData");



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




# Calculate the contingency table and p-values
overlap = overlapTable(colorHuman[inNetwork], colorChimp[inNetwork]);
# The numMat will encode color. We use -log of the p value.
numMat = -log10(overlap$pTable);
numMat[numMat >50] = 50;
# Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of
# counts and corresponding p-values.
textMat = paste(overlap$countTable, "\n", signif(overlap$pTable, 2));
dim(textMat) = dim(numMat)
# Additional information for the plot. These will be used shortly.
xLabels = paste("M", sort(unique(colorChimp)));
yLabels = paste("M", sort(unique(colorHuman)));
xSymbols = paste(sort(unique(colorChimp)), ": ", table(colorChimp[inNetwork]), sep = "")
ySymbols = paste(sort(unique(colorHuman)), ": ", table(colorHuman[inNetwork]), sep = "")



# Open a graphical window. If plotting into a file, skip the next line.
sizeGrWindow(7, 7); fp = FALSE;
pdf(fi = spaste("humanChimp-motivationFigure-dendrosAndTable.pdf"), w = 7, h = 7.0); fp = TRUE
layout(matrix(c(1,2,5, 3,4,5), 3, 2),
heights = c(3, 1, 5.5), widths = c(1, 1));
#layout.show(5);
par(mgp = c(3, 1, 0));
plotDendroAndColors(dendrograms[[1]],
		cbind(colorHuman[inNetwork], colorChimp[inNetwork]),
		c("Human modules", "Chimp modules"),
		setLayout = FALSE,
		marAll = c(1, 6, 2.7, 0.2),
		addGuide = FALSE,
		main = "A. Human gene dendrogram\nand module colors", cex.main = 1.2,
		dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.95)
par(mgp = c(3, 1, 0));
plotDendroAndColors(dendrograms[[2]],
			cbind(colorHuman[inNetwork], colorChimp[inNetwork]),
			c("Human modules", "Chimp modules"),
			setLayout = FALSE,
			marAll = c(1, 6, 2.7, 0.2),
			addGuide = FALSE,
			main = "B. Chimp gene dendrogram\nand module colors", cex.main = 1.2,
			dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.95)
# Plot the overlap table
fcex = 1.00;
pcex = 1.0
fcexl = 1.00;
pcexl = 1.00;
par(mar = c(6, 7, 2, 1.0));
labeledHeatmap(Matrix = numMat,
		xLabels = xLabels, xSymbols = xSymbols,
		yLabels = yLabels, ySymbols = ySymbols,
		colorLabels = TRUE,
		colors = greenWhiteRed(100)[50:100],
		textMatrix = textMat, cex.text = if (fp) fcex else pcex, setStdMargins = FALSE,
		cex.lab = if (fp) fcexl else pcexl,
		xColorWidth = 0.08,
		main = "C. Human modules (rows) vs. Chimp modules (columns)", cex.main = 1.2)
dev.off()




