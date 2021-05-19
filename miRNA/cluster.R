require(ppclust)
require(factoextra)
require(dplyr)
require(cluster)
require(fclust)

data(iris)
x=iris[,-5]
x[1:5,]

c1<-read.table("cluster_1",header=T,row.names=1,sep="\t")
#x<-as.matrix(t(c1))
res.fcm <- fcm(c1, centers=5)
membership<-res.fcm$u[,1]
colorbar.plot( 2, 75, y, horizontal=FALSE, col=rainbow(105))
#plot(c(1:5),c1[1,],type="l") 
pdf("membership.pdf",width=10,height=10)
par(mar=c(4.5, 4.5, 1, 6))
plot(c(1:5),c1[1,c(1:5)],xlim=c(0.5,5), ylim=c(-2,1),type="l",col=rainbow(105)[1],axes=FALSE, ann=FALSE)
for (i in 2:105){# Draw first line
	lines(c(1:5),c1[i,c(1:5)], type = "l", col = rainbow(105)[i])                   # Add second line
}  # Add third line

axis(1, at=1:5, lab=c("W8","W4","W2","W12","D1"))
# Make x axis using Mon-Fri labels

# Make y axis with horizontal labels that display ticks at 
# every 4 marks. 4*0:g_range[2] is equivalent to c(0,4,8,12).
axis(2, las=1)

# Create box around plot
#box()
# Label the x and y axes with dark green text
title(xlab="Stages", col.lab=rgb(0,0.5,0))
title(ylab="expression level", col.lab=rgb(0,0.5,0))
box()

#library(s2dverification)
#lims <- seq(-1, 1, 0.2)
#ColorBar(c1[,6], cols=rainbow(105))
zr<- range(c1[,6])
image.plot( legend.only=TRUE, zlim= zr,col=rainbow(105)) 

dev.off()

text(0,1,"Christmas strip", cex=2)
library(fields)
colorbar.plot(6, 0, c1[,6], horizontal=FALSE, col=rainbow(105),legend.only=TRUE)

text(0,1,"Christmas strip", cex=2)
# NOT RUN {
# set up a plot but don't plot points  and no "box"
plot( 1:10, (1:10)*10, type="n", bty="n") 
# of course this could be anything 

y<- cbind( 1:15, (1:15)+25)

colorbar.plot( 2.5, 30, y)
points( 2.5,30, pch="+", cex=2, adj=.5)
