
## installation

install.packages("rMVP")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("XiaoleiLiuBio/rMVP")

library(rMVP)
setwd("J:\\MAJORBIO\\project\\GWAS\\juanzhu")
setwd("G:\\gwas")
# Full-featured function (Recommended)

MVP.Data(fileNum="numeric.txt",
         filePhe="phen.txt",
         fileMap="map.txt",
         sep.num="\t",
         sep.map="\t", 
         sep.phe="\t",
         fileKin=T,
         filePC=FALSE,
         out="mvp.num",
         #priority="memory"£¬
         #maxLine=10000
)

# Only convert genotypes
#MVP.Data.Numeric2MVP("numeric.txt", out='mvp', maxLine=1e4, priority='speed', auto_transpose=T)

genotype <- attach.big.matrix("mvp.num.geno.desc")
phenotype <- read.table("mvp.num.phe",head=TRUE)
map <- read.table("mvp.num.geno.map" , head = TRUE)
#MVP.Hist(phe=phenotype[,c(1,i)], file="png", breakNum=30, dpi=300)
for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    nPC.MLM=3,
    nPC.FarmCPU=3,
    #priority="speed",
    ncpus=1,
    vc.method="EMMA",
    maxLoop=5,
    method.bin="FaST-LMM",
    threshold=0.05,
    method=c( "MLM", "FarmCPU")
  )
}
