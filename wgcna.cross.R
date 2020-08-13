# modulePreservation simulation study in which none of the modules are preserved.

# Generate two data sets with the same module structure, then in one
# of them permute the columns. Thus, none of the modules should be preserved.

source("../../CommonFunctions/networkFunctions-extras-04.R");
source("../../RLibs/WGCNA/R/Functions.R");
source("../../RLibs/WGCNA/R/modulePreservation.R");
# Set up simulation parameters

nSamples = c(100, 100);
nMods = 10;
nSets = length(nSamples);
nAllSamples = sum(nSamples);

# Simulate eigengenes

set.seed(1)

allEigengenes = matrix(rnorm(nAllSamples * nMods), nAllSamples, nMods);

modSizes = c(1000, 500, 250, 100, 50, 1000, 500, 250, 100, 50, 300);
nGenes = as.integer(sum(modSizes)*1.2);
modProps = modSizes/nGenes;

eigengenes = list();
eigengenes[[1]] = list(data = allEigengenes[1:nSamples[1], ]);
eigengenes[[2]] = list(data = allEigengenes[(nSamples[1]+1):nAllSamples, ]);

# Simulate data

leaveOut = matrix(FALSE, nMods, nSets);

data = simulateMultiExpr(eigengenes, nGenes, modProps, leaveOut = leaveOut,
                         minCor = 0.3, maxCor = 0.9, backgroundNoise = 0.2);

multiExpr = data$multiExpr;

colnames(multiExpr[[1]]$data) = spaste("Gene", c(1:nGenes));
colnames(multiExpr[[2]]$data) = spaste("Gene", c(1:nGenes));

# permute columns in test data set

set.seed(3)
perm = sample(c(1:nGenes));
multiExpr[[2]]$data = multiExpr[[2]]$data[, perm];
colnames(multiExpr[[2]]$data) = spaste("Gene", c(1:nGenes));

# Identify modules

mods = list();
for (set in 1:nSets)
   mods[[set]] = blockwiseModules(multiExpr[[set]]$data, numericLabels = TRUE, verbose = 4)

table(mods[[1]]$colors, data$allLabels[, 1])
table(mods[[2]]$colors, data$allLabels[, 2])

colorList = list();

colorList[[1]] = mods[[1]]$colors;
colorList[[2]] = mods[[2]]$colors;

expr0 = multiExpr;
color0 = colorList;

# Calculate module preservation

doMP = TRUE;
if (doMP)
{
  print( system.time( {
  mp = modulePreservation(multiExpr, colorList, referenceNetworks=1,
                          nPermutations = 100,
                          networkType = "unsigned",
                          randomSeed = 2345,
                          permutedStatisticsFile = "preservedButPermuted-noReg-permStats.RData",
                          verbose = 4, indent = 0)
             } ));
  # Save the results
  save(mp, file= "simulation-preservedButPermuted.RData");
}

# Calculate IGP using clusterRepro

doClusterRepro = TRUE;
if (doClusterRepro)
{
  eigengenes = as.matrix(moduleEigengenes(multiExpr[[2]]$data, colorList[[1]], excludeGrey = TRUE)$eigengenes);
  library(clusterRepro);

  rownames(eigengenes) = rownames(multiExpr[[2]]$data);
  set.seed(20)
  print(system.time( {
    cr = clusterRepro(Centroids = eigengenes, New.data = multiExpr[[2]]$data, Number.of.permutations = 10000);
                      } ));
  save(cr, file = "simulation-preservedButPermuted-cr.RData");
}

stop("All done.");

load(file= "simulation-preservedButPermuted.RData");
load(file= "simulation-preservedButPermuted-cr.RData");

# Prepare statistics for plotting

stats = cbind(mp$quality$Z[[1]][[2]][, -1],
              mp$referenceSeparability$Z[[1]][[2]][, -1, drop = FALSE],
              mp$preservation$Z[[1]][[2]][, -1],
              mp$accuracy$Z[[1]][[2]][, -1],
              mp$testSeparability$Z[[1]][[2]][, -1, drop = FALSE]);

order = order(as.numeric(rownames(stats)))
stats = stats[order, ]
labelsX = as.numeric(rownames(stats))
labelsX[labelsX==0.1] = 25
colors = labels2colors(labelsX);

moduleSizes = as.numeric(table(colorList[[1]]));

letter = labelsX;
presColor = rep("black", length(letter));
useStats = c(1:ncol(stats))[-19];

useStats = list(c(1:8), c(9:ncol(stats)));
sectioning = list(c(2,4), c(4,5));
dims = list(sectioning[[1]] * 2.0, sectioning[[2]]*2.0);
figNames = c("quality", "preservation");

numLabels2 = letter;
numLabels2[c(3, 5, 7, 9, 11)] = seq(from = 5, to=1, by = -1);
numLabels2[c(4, 6, 8, 10, 12)] = seq(from = 10, to=6, by = -1);
numLabels2 = as.numeric(numLabels2)

# Plot results into two pdf files, one for quality and one for preservation statistics

#sizeGrWindow(10,7);
for (f in 1:2)
{
  pdf(file=spaste("simulation-preservedButPermuted-",figNames[f],".pdf"), w=dims[[f]][2], h=dims[[f]][1])
  par(mfrow = sectioning[[f]])
  par(mar = c(3.2, 3.2, 2, 0.5))
  par(mgp = c(2.0, 0.6, 0))

  for (s in useStats[[f]])
  {
    # Shift module sizes on modules 6--10 to make them distinct from 1-5 in the plot.
    ms = moduleSizes[-1];
    ms[numLabels2[-c(1:2)] > 5] = ms[numLabels2[-c(1:2)] > 5] + 20;

    min = min(stats[-c(1:1), s], na.rm = TRUE);
    max = max(stats[-c(1:1), s], na.rm = TRUE);
    if (min > -2.5) min = -2.5
    if (max <  11  ) max =  11
    if (min > -max/5) min = -max/5;
    plot(ms, stats[-c(1:2), s], col = colors, pch = 20,
         main = colnames(stats)[s],
         cex = 2,
         ylab = colnames(stats)[s], type = "n", xlab = "Module size",
         cex.main = 1, ylim = c(min, max))
    text(ms, stats[-c(1:2), s], labels = numLabels2[-c(1:2)], col = presColor[-c(1:2)]);
    box = par("usr");
    abline(h=0)
    abline(h=2, col = "blue", lty = 2);
    abline(h=10, col = "darkgreen", lty = 2);
  }
  dev.off();
}

msx = ms
msx[numLabels2[-c(1:2)] > 5] = ms[numLabels2[-c(1:2)] > 5] + 40;

# Plot results of clusterRepro

sizeGrWindow(4,8);
#pdf(file="Plots/simulation-preservedButPermuted-clusterRepro.pdf", wi=4, he=8);
par(mfrow = c(2,1))
par(mar = c(3.2, 3.2, 2, 0.5))
par(mgp = c(2.0, 0.6, 0))
plot(ms, cr$Actual.IGP[-1], col = colors, pch = 20,
     main = "B. Preserved but permuted\nObserved IGP",
     cex = 2,
     ylab = "Observed IGP", type = "n", xlab = "Module size",
     cex.main = 1)
text(ms, cr$Actual.IGP[-1], labels = numLabels2[-c(1:2)], col = presColor[-c(1:2)]);
addGrid();

plot(msx, cr$p.value[-1], col = colors, pch = 20,
     main = "E. Preserved but permuted\nIGP permutation p-value",
     cex = 2,
     ylab = "p-value", type = "n", xlab = "Module size",
     cex.main = 1)
text(msx, cr$p.value[-1], labels = numLabels2[-c(1:2)], col = presColor[-c(1:2)]);
addGrid();

dev.off();




# Plot cross-tabulation


sizeGrWindow(4,4);
pdf(file="simulation-preservedButPermuted-crossTabulation.pdf", wi=4, he=4);
par(mar = c(3.2, 3.2, 2, 0.5))
par(mgp = c(2.0, 0.6, 0))
plot(ms, mp$accuracy$observed[[1]][[2]][-1, "minusLogFisherP"], col = colors, pch = 20,
     main = "B. Preserved but permuted\n-log(Fisher test p-value)",
     cex = 2,
     ylab = "-log(Fisher test p-value)", type = "n", xlab = "Module size",
     cex.main = 1)
text(ms, mp$accuracy$observed[[1]][[2]][-1, "minusLogFisherP"],
     labels = numLabels2[-c(1:2)], col = presColor[-c(1:2)]);
addGrid();

dev.off();

sizeGrWindow(4,4);
#pdf(file="Plots/simulation-preservedButPermuted-coClustering.pdf", wi=4, he=4);
par(mar = c(3.2, 3.2, 2, 0.5))
par(mgp = c(2.0, 0.6, 0))
plot(ms, mp$accuracy$observed[[1]][[2]][-1, "coClustering"], col = colors, pch = 20,
     main = "E. Preserved but permuted\nCo-clustering",
     cex = 2,
     ylab = "Co-clustering", type = "n", xlab = "Module size",
     cex.main = 1)
text(ms, mp$accuracy$observed[[1]][[2]][-1, "coClustering"],
     labels = numLabels2[-c(1:2)], col = presColor[-c(1:2)]);
addGrid();

dev.off();






#Is My Network Module Preserved and Reproducible






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
# Names for the two sets
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

# Here comes the IGP calculation
library(clusterRepro)
cr = list();
set.seed(20);
for (ref in 1:nSets){
	cr[[ref]] = list();
	for (test in 1:nSets){
		printFlush(system.time({
		cr[[ref]][[test]] = clusterRepro(Centroids = as.matrix(eigengenes[[ref]][[test]]$data),
		New.data = as.matrix(impExpr[[test]]$data),
		Number.of.permutations = 100); }));#Number.of.permutations = 10000
		collectGarbage();
	}
}

# Save the results
save(cr, file = "HumanChimp-HumanSpecific-clusterRepro.RData");



# Load the module preservation statistics
load(file = "HumanChimp-HumanSpecific-clusterRepro.RData");
load(file = "HumanChimp-HumanSpecific-modulePreservation.RData")
ref = 1 # Select the human data as reference
test = 2 # Select the chimp data as test
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);


print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2));
# Compare preservation to quality:
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))


ref = 1;
test = 2;
ind = 1;
stats= mp$preservation$observed[[ref]][[test]];
labelsX = rownames(stats)
labelsX[labelsX=="gold"] = "orange"
modColors = labelsX;
plotMods = !(modColors %in% c("grey", "orange"));
moduleSizes = stats[plotMods, 1];
textLabels = match(modColors, standardColors(20))[plotMods];
colorLabels = labelsX[plotMods];


nModules = sum(plotMods);
nPlots = 6
plotData = list();
# Fill up the plotData
plotData[[1]] = plotData[[2]] = matrix(0, nModules, nPlots);
plotData[[1]][, c(1:4)] = moduleSizes;
plotData[[2]][, 1] = mp$preservation$Z[[ref]][[test]]$Zsummary.pres[plotMods];
plotData[[2]][, 2] = mp$preservation$observed[[ref]][[test]]$medianRank.pres[plotMods];
# Match the modulePreservation ordering of modules to that of clusterRepro
crLabels = sort(unique(colorLabels));
mp2cr = match(colorLabels, crLabels);
# Scatterplots of IGP and p-value vs. module size
plotData[[2]][, 3] = cr[[ref]][[test]]$Actual.IGP
plotData[[2]][, 4] = -log10(cr[[ref]][[test]]$p.value + 1e-4);
# Scatterplot of observed IGP vs. Zsummary and medianRank
plotData[[1]][, c(5,6)] = plotData[[2]][, c(1:2)];
plotData[[2]][, c(5,6)] = plotData[[2]][, 3];
# Plot annotation
xLabs = c(rep("Module size", 4), "Zsummary", "Median rank");
yLabs = c("Zsummary", "Median rank", "Observed IGP", "-log10(IGP perm p)", "Observed IGP", "Observed IGP");
mains = spaste(LETTERS[1:nPlots], ". ", #rep("Ref: Human, Test: Chimp\n", nPlots),
c(yLabs[1:4], paste(yLabs[5:6], "vs.", xLabs[5:6])),
c("", "", "", "", "\n", "\n"));
# Scatterplot options
verbose = c(rep(FALSE, 4), rep(TRUE, 2));
ablines = list(c(0, 2, 10), NA, NA, c(-log10(0.05), -log10(0.05/nModules)), NA, NA);
abColors = list(c("black", "blue", "darkgreen"), NA, NA, c("blue", "red"), NA, NA);
logs = c("x", "x", "x", "x", "", "");
invertY = c(FALSE, TRUE, rep(FALSE, 4));
verSP = function(...) { verboseScatterplot(..., abline = TRUE) }


cexLabels = 1.4
sizeGrWindow(6,9);
#pdf(file = "Plots/HumanSpecific-NetworkAndIGPStatistics.pdf", w=6, h=9, onefile = FALSE);
par(mfrow = c(3,2));
par(mar = c(3.3, 3.3, 3.2, 0.5));
par(mgp = c(2, 0.6, 0))
for (p in 1:nPlots)
{
x = plotData[[1]][, p];
y = plotData[[2]][, p]
miny = min(y, ablines[[p]], na.rm = TRUE);
maxy = max(y, ablines[[p]], na.rm = TRUE);
miny = miny - (maxy-miny)*0.1;
maxy = maxy + (maxy-miny)*0.1;
(if (verbose[p]) verSP else plot ) (plotData[[1]][, p], plotData[[2]][, p],
main = mains[p],
xlab = xLabs[p],
ylab = yLabs[p],
cex.main = cexLabels, cex.lab = cexLabels, cex.axis = cexLabels,
bg = colorLabels,
col = colorLabels, cex = 2.2,
ylim = if (invertY[p]) c(maxy, miny) else c(miny, maxy),
pch = 21,
log = logs[p]);
labelPoints(plotData[[1]][, p], plotData[[2]][, p], textLabels, cex = cexLabels, offs = 0.06);
if (!is.na(ablines[[p]][[1]]))
for (al in 1:length(ablines[[p]]))
abline(h = ablines[[p]][[al]], col = abColors[[p]][[al]], lty = 2);
}
# If plotting into a pdf file, close the file. An un-closed pdf file is not readable.
dev.off();




figNames = c("quality", "preservation");
useStats = list(c(1:8)+1, c(9:27)+1);
sectioning = list(c(2,4), c(4,5));
dims = list(sectioning[[1]] * 1.8, sectioning[[2]]*1.8);
for (f in 1:2)
{
pdf(file=spaste("Plots/HumanChimp-HumanSpecific-modulePreservation-",figNames[f],"-%02d.pdf"),
w=dims[[f]][2], h=dims[[f]][1], onefile = FALSE)
for (ref in 1:nSets) for (test in 1:nSets) if (ref!=test)
{
stats = cbind(mp$quality$Z[[ref]][[test]],
mp$referenceSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$Z[[ref]][[test]][, -1],
mp$accuracy$Z[[ref]][[test]][, -1],
mp$testSeparability$Z[[ref]][[test]][, -1, drop = FALSE]);
labelsX = rownames(stats)
labelsX[labelsX=="gold"] = "orange"
modColors = labelsX;
moduleSizes = stats[, 1];
plotMods = !(modColors %in% c("grey", "orange"));
plotStats = useStats[[f]]
textLabels = match(modColors, standardColors(20))[plotMods];
par(mfrow = sectioning[[f]]);
par(mar = c(3.0, 3.0, 4, 0.4));
par(mgp = c(1.7, 0.6, 0))
for (s in plotStats)
{
min = min(stats[plotMods, s], na.rm = TRUE);
max = max(stats[plotMods, s], na.rm = TRUE);
minMS = min(moduleSizes[plotMods]);
maxMS = max(moduleSizes[plotMods]);
nms = colnames(stats);
if (max < 10) max = 10;
if (min > -max/5) min = -max/5
plot(moduleSizes[plotMods], stats[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
main = paste("Ref.:", setLabels[ref], "\nTest:", setLabels[test], "\n", nms[s]),
cex = 1.2,
cex.main = 1.0,
ylab = nms[s], xlab = "Module size", log = "x",
ylim = c(min, max + 0.1 * (max-min)),
xlim = c(minMS/(maxMS/minMS)^0.1, maxMS*(maxMS/minMS)^0.1)
)
labelPoints(moduleSizes[plotMods], stats[plotMods, s], textLabels, cex = 0.90, offs = 0.06);
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
}




# This variable will contain the summary table
summaryTable = NULL
# Loop over all combinations of reference and tests sets
for (ref in 1:nSets) for (test in 1:nSets) if (ref!=test)
{
modules = rownames(mp$preservation$Z[[ref]][[test]]);
nMods = length(modules);
sizes = mp$preservation$Z[[ref]][[test]][, 1];
acc = matrix(NA, nMods, 3);
if (test!=4)
{
acc[match(rownames(mp$accuracy$observed[[ref]][[test]]), modules), ] =
mp$accuracy$observed[[ref]][[test]][, -1, drop = FALSE];
colnames(acc) = colnames(mp$accuracy$observed[[ref]][[test]])[-1];
accZ = mp$accuracy$Z[[ref]][[test]][, -1, drop = FALSE];
acc.log.p = mp$accuracy$log.p[[ref]][[test]][, -1, drop = FALSE];
acc.log.pBonf = mp$accuracy$log.pBonf[[ref]][[test]][, -1, drop = FALSE];
} else {
accZ = matrix(NA, nMods, 3);
acc.log.p = matrix(NA, nMods, 3);
acc.log.pBonf = matrix(NA, nMods, 3);
colnames(acc) = colnames(mp$accuracy$observed[[1]][[2]])[-1];
colnames(accZ) = colnames(mp$accuracy$Z[[1]][[2]])[-1];
colnames(acc.log.p) = colnames(mp$accuracy$log.p[[1]][[2]])[-1];
colnames(acc.log.pBonf) = colnames(mp$accuracy$log.pBonf[[1]][[2]])[-1];
}
# Table of results for this reference-test combination
tab = cbind(referenceSet = rep(setLabels[ref], nMods),
testSet = rep(setLabels[test], nMods),
moduleLabel = modules,
moduleSize = sizes,
mp$quality$observed[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$observed[[ref]][[test]][, -1, drop = FALSE],
acc,
mp$referenceSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
mp$testSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
mp$quality$Z[[ref]][[test]][, -1, drop = FALSE],
mp$quality$log.p[[ref]][[test]][, -1, drop = FALSE],
mp$quality$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$Z[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$log.p[[ref]][[test]][, -1, drop = FALSE],
mp$preservation$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
accZ,
acc.log.p,
acc.log.pBonf,
mp$referenceSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
mp$referenceSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
mp$referenceSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
mp$testSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
mp$testSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
mp$testSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE]
)
# Add the table to the main table.
if (is.null(summaryTable)) summaryTable = tab else summaryTable = rbind(summaryTable, tab);
}
# Save the table in csv format.
write.table(summaryTable, file = "HumanChimp-HumanSpecific-completeResults.csv", row.names = FALSE,
sep = ",", quote = FALSE);

dendrograms = list();
for (set in 1:nSets)
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

r_sample_logical(100, name = "Var")
table(r_sample_logical(1000))
c("B", "W")[r_sample_logical(10)]

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


# Use human set as reference, chimp as test
ref = 1;
test = 2
# The area from which each sample waas taken is encoded in the sample names.
sampleNames = rownames(multiExpr[[ref]]$data)
cortical = substring(sampleNames, 4) %in% c("Brocas.", "prv", "prf", "acc", "Brocas")
cn = substring(sampleNames, 4) %in% c("CN");
cerebellum = substring(sampleNames, 4) %in% c("VC");
# Cortical samples will be coded lightblue
sampleIndicator = c("lightblue", "white")[2-as.numeric(cortical) ];
# CN samples are coded orange
sampleIndicator[cn] = "orange";
# Cerebellum samples are code turquoise
sampleIndicator[cerebellum] = "magenta";



# We will plot these two modules. Modify to get plots for other modules.
plotMods = c("yellow", "blue");
# Open a graphical window. If plotting into a file, skip the next line and instead uncomment and paste the
# line starting with #pdf.
sizeGrWindow(10,5.5);
#pdf(fi = spaste("Plots/humanChimp-motivationFigure-moduleSpecific.pdf"), w = 10, 5.5); fp = TRUE
# Set the plot layout (sectioning). Type help(layout) to see what layout does.
layout(matrix(c(1,2,6,7, 3,3,8,8, 4,4,9,9, 5,5,10,10), 4, 4),
heights = rep(c(1.0, 1.0), 2));
# Optional: show the sectioning. Skip if plotting into a file.
layout.show(10);
# Keep track of number of plotted panels.
ind = 1;
# Loop over modules to be displayed.
for (mod in plotMods)
{
# The next big chunk of the code creates a module heatmap and an eigengene plot beneath. Samples are
# ordered by the area from which they were taken. We promise to one day create a simple function that
# will create this plot in a single call.
# Here we prepare various variables for the plot
order = order(sampleIndicator)
refX = 1;
modx = match(mod, substring(colnames(mes[[refX]]), 3));
if (refX==1) par(mar = c(0.2,2.5,3,0.2)) else par(mar = c(0.2,0.2,3,0.5))
nModGenes = sum(colorList[[ref]]==mod);
modExpr = scale(multiExpr[[ref]]$data[, colorList[[ref]]==mod])
geneOrder = hclust(as.dist(1-cor(modExpr, use="p")), method = "average")$order
# sampleOrder = hclust(dist(modExpr), method = "average")$order;
sampleOrder = order;
MaxExpr = max(modExpr, na.rm = TRUE);
MinExpr = min(modExpr, na.rm = TRUE);
MaxZ = max(MaxExpr, abs(MinExpr));
zlim = c(-MaxZ, MaxZ);
palette = c("grey", greenWhiteRed(50));
exprColIndex = (modExpr + MaxZ)/2/MaxZ * 49 + 2;
exprColIndex[is.na(exprColIndex)] = 1;



eigengene = mes[[refX]][, modx]
# Set up the plotting coordinates by plotting an empty barplot
par(mgp = c(0, 0, 0));
if (refX==1)
{
bp = barplot(as.vector(eigengene[sampleOrder]), col = "white", border = "white", axisnames = FALSE,
main = spaste(LETTERS[ind], ". Module ", mod, " in human"),
axes = FALSE, ylab = "Module expression\nheatmap", cex.lab = 0.9);
} else {
bp = barplot(as.vector(eigengene[sampleOrder]), col = "white", border = "white", axisnames = FALSE,
main = spaste("male CTX95"),
axes = FALSE, ylab = "", cex.lab = 0.9);
}
ind = ind + 1;
# Get the coordinates of the plot rectangle
plotbox = par("usr");
xstep = bp[2]-bp[1]; xLeft = bp - xstep/2; xRight = bp + xstep/2;
nrows = ncol(modExpr);
yRange = plotbox[4]-plotbox[3];
yBot = plotbox[3] + c(0:(nrows-1)) * yRange/nrows;
yTop = yBot + yRange/nrows;
# Plot the actual heatmap
for (sample in 1:nrow(modExpr))
{
colorInd = as.integer(exprColIndex[sampleOrder[sample], geneOrder ]);
rect(xleft = rep(xLeft[sample], nrows), xright = rep(xRight[sample], nrows),
ybottom = yBot, ytop = yTop, col = palette[colorInd], border = palette[colorInd]);
}
# -----------------------------------------------------------------------------------
# The next chunk of code plots the color-coded eigengene barplot.
if (refX==1) par(mar = c(2.4, 2.5, 0.2, 0.2)) else par(mar = c(2,0.2,0.2,0.5))
par(mgp = c(0.5, 0, 0));
MaxE = max(abs(eigengene), na.rm = TRUE);
eigengene[is.na(eigengene)] = 0;
index = as.integer(as.vector(eigengene + MaxE)/(2*MaxE) * 50 + 1);
Colors = greenWhiteRed(50)[index];
barplot(as.vector(eigengene[sampleOrder]), col = Colors[sampleOrder],
xlab = "", ylab = if (refX==1) "Eigengene expression" else "", cex.lab = 0.9, yaxt = "none");
# -----------------------------------------------------------------------------------
# This part plots the sample indicator as a color-coded bar.
plotbox = par("usr"); yRange = plotbox[4] - plotbox[3]
rect(xleft = xLeft, xright = xRight, ytop = plotbox[3] - yRange/10, ybot = plotbox[3] - 1.5 *yRange/10,
col = sampleIndicator[sampleOrder],
border = sampleIndicator[sampleOrder], xpd = TRUE);
text(mean(c(xLeft, xRight)), plotbox[3] - 2.2*yRange/10, "Sample brain area indicator", adj = c(0.5,1),
xpd = TRUE, cex = 0.9);
# This concludes the heatmap and eigengene barplots. The scatterplots are much easier to generate.
# Restrict to module genes:
modGenes = colorList[[ref]]==mod;
# -----------------------------------------------------------------------------------
# Calculate intramodular correlations and adjacencies
corRef = vectorizeMatrix(cor(multiExpr[[ref]]$data[, modGenes], use = "p"));
corTest= vectorizeMatrix(cor(multiExpr[[test]]$data[, modGenes], use = "p"));
adjRef = adjacency(multiExpr[[ref]]$data[, modGenes], power = 6, type = "signed hybrid");



adjTest= adjacency(multiExpr[[test]]$data[, modGenes], power = 6, type = "signed hybrid");
# Select a sample to plot (using all gene-gene pairs would result in a needlessly large plot)
set.seed(10)
select = sample(1:length(corRef), 3000);
# Set margin parameters
mar = c(3.8, 3.8, 3.8, 0.7)
mgp = c(2.1, 0.7, 0);
par(mar = mar);
par(mgp = mgp);
# Create the scatterplot. Note that the correlation and p-value are calculated from all data, the sample
# only used for plotting.
verboseScatterplot(corRef, corTest, sample = select,
xlab = spaste("Correlation in human data"),
ylab = spaste("Correlation in chimp data"),
abline = TRUE,
main = spaste(LETTERS[ind], ". Intramodular correlation\n",
"in human ", mod, " module\n"),
corLabel = "cor.cor",
cex.main = 1.2, cex.lab = 1, cex.axis = 1, cex = 0.6)
ind = ind + 1;
# -----------------------------------------------------------------------------------
# We now plot a scatterplot of intramodular connectivities. We first calculate kIM:
diag(adjRef) = NA;
diag(adjTest) = NA;
kIMRef = apply(adjRef, 1, sum, na.rm = TRUE);
kIMTest = apply(adjTest, 1, sum, na.rm = TRUE);
# Plot the scatterplot:
color = 1;
verboseScatterplot(kIMRef, kIMTest,
xlab = spaste("kIM in human data"),
ylab = spaste("kIM in chimp data"),
abline = TRUE,
main = spaste(LETTERS[ind], ". Intramodular connectivity\n",
"in human ", mod, " module\n"), corLabel = "cor.kIM",
cex.main = 1.2, cex.lab = 1, cex.axis = 1, cex = 0.6, col = color)
ind = ind + 1;
# -----------------------------------------------------------------------------------
# Lastly, we plot a scatterplot of module eigengene-based connectivities:
# Calculate module eigengene-based connectivities
kMERef = cor(multiExpr[[ref]]$data[, modGenes], mes[[ref]][, modx], use = "p")
kMETest= cor(multiExpr[[test]]$data[, modGenes], mes[[test]][, modx], use = "p")
# Plot the scatterplot
verboseScatterplot(kMERef, kMETest,
xlab = spaste("kME in human data"),
ylab = spaste("kME in chimp data"),
abline = TRUE,
main = spaste(LETTERS[ind], ". Eigengene-based connect.\n",
"in human ", mod, " module\n"), corLabel = "cor.kME",
cex.main = 1.2, cex.lab = 1, cex.axis = 1, cex = 0.6, col = color)
ind = ind + 1;
}
# If plotting into a pdf file, close the file. An opened file cannot be used.
dev.off();




