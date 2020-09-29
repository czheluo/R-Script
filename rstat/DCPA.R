## prepare data
library('vcfR')
vcf <- read.vcfR("pop.vcf.gz")
x <- vcfR2genlight(vcf)

# DAPC requires the adegenet package. Let's load this package:
library("adegenet")
group<-read.table("3.group.list",header=F)
pop(x)<-group[,2]

dapc <- dapc(x, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(x) - 1)

pdf("Dapc.1.pdf",width=10,height=10)

scatter(dapc, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

dev.off()

#system.time(dapc <- dapc(x, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(x) - 1))

pdf("Dapc.predic.pdf",width=10,height=10)

system.time(pramx <- xvalDapc(tab(x, NA.method = "mean"), pop(x),
                             n.pca = 10:40, n.rep = 1000,
                             parallel = "multicore", ncpus = 4L))
dev.off()


pdf("Dapc.2.pdf",width=10,height=10)

scatter(pramx$DAPC,  cex = 2, legend = TRUE,#col = other(x)$comparePal,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

dev.off()
save.image("result.RData")

