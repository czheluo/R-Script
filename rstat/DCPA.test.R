## prepare data
library('vcfR')
vcf <- read.vcfR("pop.vcf.gz")
x <- vcfR2genlight(vcf)

# DAPC requires the adegenet package. Let's load this package:
library("adegenet")
dapc <- dapc(x, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(x) - 1)
scatter(dapc, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)



system.time(dapc <- dapc(x, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(x) - 1))

system.time(pramx <- xvalDapc(tab(x, NA.method = "mean"), pop(x),
                             n.pca = 5:20, n.rep = 1000,
                             parallel = "multicore", ncpus = 4L))



pdf("Dapc.34.pdf",width=10,height=10)

scatter(pramx$DAPC,  cex = 2, legend = TRUE,#col = other(x)$comparePal,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

dev.off()



library("poppr")
data("Pram", package = "poppr")
Pram


set.seed(999)
system.time(pramx <- xvalDapc(tab(Pram, NA.method = "mean"), pop(Pram),
                             n.pca = 5:20, n.rep = 1000,
                             parallel = "multicore", ncpus = 4L))


system.time(pramx <- xvalDapc(tab(H3N2, NA.method = "mean"), pop(H3N2),
                             n.pca = 5:20, n.rep = 1000,
                             parallel = "multicore", ncpus = 4L))



