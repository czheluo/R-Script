
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

a=experData[[1]]$y
z <- array(0, dim=c(dim(experData[[1]]$y),20))
for(i in 1:dim(z)[1]){
  for(j in 1:dim(z)[2]){
    z[i,j,]=rnorm(20,mean=a[i,j],sd=abs(a[i,j]*0.5))
  }
}
test_result=array(0,c(dim(geneNetOpt$model$parametersW),20))


for(t in 1:20){
  b=z[,,t]
  rownames(b)=rownames(a)
  colnames(b)=colnames(a)
  experData[[1]]$y=b
  result <- netGenGenerateVariables(experData, ## experimental data
                                    experimentNames = experNames, ## experiment name(s)
                                    nodeNames = nodeNames, ## gene/node names
                                    allowedError = 0.003, ## parameter weighting and thresholding error between measurement and simulation
                                    maxNoForwardConnect = 3, ## maximal number of forward edges
                                    maxNoFeedbacks = 0, ## maximal number of (global) feedback edges
                                    weightingStruct = 0.08, ## parameter weighting the overall influence of structural information (prior knowledge)
                                    weightingInterpol = 0.1, ## parameter weighting interpolated values
                                    interpolSteps = 50, ## number of interpolation steps
                                    maxDynamicOrder = 1, ## maximal dynamic order of each submodel
                                    #connectionsW.flex = priorflex, ## prior knowledge for flexible integration
                                    #scoresW.flex = score, ## scores of prior knowledge
                                    flagLinearModel = TRUE, ## generate linear model
                                    #flagConsiderUncertainties = TRUE ## consider uncertainties
  )
  geneNet <- result$geneNet
  dataSet <- result$dataSet
  geneNetOptTest <- netGenTrain(geneNet, dataSet, verbose = TRUE)
  test_result[,,t]=geneNetOptTest$model$parametersW
}


dataRuns=test_result
dataRuns[dataRuns<0]=-1
dataRuns[dataRuns>0]=1
netGenBubbleMap(dataRuns,rep(0.05,20),names = geneNetOpt$model$nodeNames,fileName = "bubbleMap")

save.image("case6.time1.RData")
