library(MASS)
mat<-as.data.frame(matrix(0,5,2))
pmat<-as.data.frame(matrix(0,6606,2))
for(j in 1:6606){
	c<-c(1,2)
	for(i in 1:5){
		mat[i,]=ta[j,c]
		c<-c+2
	}
	chi<- chisq.test(mat)
	pmat[j,]<-c(chi$statistic,chi$p.value)
}

colnames(pmat)<-c("X2","Pvalue")

write.csv(pmat,file="chisq.csv",quote=F,row.names=F)
