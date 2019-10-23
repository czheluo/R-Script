
setwd("I:\\MAJORBIO\\liuhaifang\\Result\\case6/time1/")


library(NetGenerator)
nod<-read.table("nodename.txt",header=T)
nodeNames <- as.character(colnames(nod))
lfc.pathogen <- read.table("case6.path.txt",sep="\t",header=TRUE,row.names=1)
lfc.host <- read.table("case6.host.txt",sep="\t",header=TRUE,row.names=1)
sd.pathogen <- read.table("case6.path.txt",sep="\t",header=TRUE,row.names=1)
sd.host <- read.table("case6.host.txt",sep="\t",header=TRUE,row.names=1)

experData <- list(list(y=NULL, ty=NULL, u=NULL, tu=NULL, yUncertainties=NULL))
experData[[1]]$y <- as.matrix(cbind(lfc.pathogen,lfc.host)) ## output data (lfc)
experData[[1]]$ty <- c(0, 5, 10, 15, 20) ## output time vector
experData[[1]]$u <- matrix(rep(1, 5), 5, 1) ## input (stimulus)
experData[[1]]$tu <- experData[[1]]$ty ## input time vector
experData[[1]]$yUncertainties <- as.matrix(cbind(sd.pathogen,sd.host))
experNames <- "experiment"
result <- netGenGenerateVariables(experData,
                                  experimentNames = experNames,
                                  nodeNames = nodeNames,
                                  allowedError = 0.003,
                                  maxNoForwardConnect = 3,
                                  maxNoFeedbacks = 0,
                                  weightingStruct = 0.08,
                                  weightingInterpol = 0.1,
                                  interpolSteps = 50,
                                  maxDynamicOrder = 1,
                                  flagLinearModel = TRUE,
                                  flagZeros=TRUE,
                                  
)


geneNet <- result$geneNet
dataSet <- result$dataSet
geneNetOpt <- netGenTrain(geneNet, dataSet, verbose = TRUE)


netGenPlot(geneNetOpt,dataSet,layout="matrix",plotName="case6.h1.timeCourse")
netGenGraph(geneNetOpt, graphName="p1_h1.graph",outputFileFormat="pdf")

save.image("case6.time1.RData")
