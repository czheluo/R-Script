#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
  'hapmap','h',2,'character',
  'method','m',2,'character',
  'vcf','v',2,'character',
  'maf','a',2,'numeric',
  'minDP','d',2,'numeric',
  'missing','s',2,'numeric',
  'trait','t',1,'character',
  'output','o',1,'character',
  'permutation','p',2,'numeric',
  'help','c',0,'logical'
), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
Usage example:
	Rscript indel_len.r --i  --o

Usage:
	--hapmap	the input hapmap file
	--vcf	the input vcf file
	--method	1:GLM,2:MLM,3:FarmCPU,4:THREE (default: 2)
	--trait	the trait file
	--output	the output dir
	--permutation the number of permutaion
	--maf
	--minDP
	--missing
	--help		usage
\n")
  q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$trait)) { print_usage(spec)}
if ( is.null(opt$hapmap) && is.null(opt$vcf)){ print_usage(spec) }
if ( is.null(opt$method)){ opt$method=2 }
if ( is.null(opt$output)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
if ( is.null(opt$maf)){ opt$maf=0.05 }
if ( is.null(opt$minDP)){ opt$minDP=5 }
if ( is.null(opt$missing)){ opt$missing=0.7 }

library(rMVP)
setwd(opt$output);


if (!is.null(opt$vcf)){
       vcf.file=opt$vcf
       out.file="pop.filter"
       maf=opt$maf
       minDP=opt$minDP
       missing=opt$missing
       system(sprintf("vcftools --vcf %s --maf %s --minDP %s --max-missing %s --remove-indels --out %s --recode",
                      vcf.file,maf,minDP,missing,out.file),ignore.stdout = FALSE)
       # Full-featured function (Recommended)
       MVP.Data(fileVCF="pop.filter.recode.vcf",
                filePhe=opt$trait,
                sep.hmp="\t",
                sep.phe="\t",
                fileKin=FALSE,
                filePC=FALSE,
                out="mvp.vcf"
       )
       genotype <- attach.big.matrix("mvp.vcf.geno.desc")
       phenotype <- read.table("mvp.vcf.phe",head=TRUE)
       map <- read.table("mvp.vcf.geno.map" , head = TRUE)
}else if(!is.null(opt$hapmap)){
       MVP.Data(fileHMP=opt$hapmap,
                filePhe=opt$trait,
                sep.hmp="\t",
                sep.phe="\t",
                SNP.effect="Add",
                fileKin=FALSE,
                filePC=FALSE,
                out="mvp.hmp",
                priority="memory",
       )
       genotype <- attach.big.matrix("mvp.hmp.geno.desc")
       phenotype <- read.table("mvp.hmp.phe",head=TRUE)
       map <- read.table("mvp.hmp.map" , head = TRUE)
     }

#MVP.Hist(phe=phenotype[,c(1,i)], file="png", breakNum=30, dpi=300)
if (!is.null(opt$permutation)){
       imMVP <- MVP(
         phe=phenotype,
         geno=genotype,
         map=map,
         nPC.GLM=5,
         nPC.MLM=3,
         nPC.FarmCPU=3,
         perc=1,
         priority="speed",
         ncpus=10,
         vc.method="EMMA",
         maxLoop=10,
         method.bin="FaST-LMM",
         permutation.threshold=TRUE,
         permutation.rep=opt$permutation,
         file="pdf",
         method="MLM")

}else if (opt$method==1){
       imMVP <- MVP(
         phe=phenotype,
         geno=genotype,
         map=map,
         #K=Kinship,
         #CV.GLM=Covariates,  ##if you have additional covariates, please keep there open.
         #CV.MLM=Covariates,
         #CV.FarmCPU=Covariates,
         nPC.GLM=5,   ##if you have added PCs into covariates, please keep there closed.
         #nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
         #nPC.FarmCPU=3,
         priority="speed",   ##for Kinship construction
         #ncpus=10,
         vc.method="BRENT",  ##only works for MLM
         maxLoop=10,
         method.bin="static",   ## "FaST-LMM", "static" (#only works for FarmCPU)
         #permutation.threshold=TRUE,
         #permutation.rep=100,
         threshold=0.05,
         method=c("GLM")
       )
}else if(opt$method==2){
       imMVP <- MVP(
         phe=phenotype,
         geno=genotype,
         map=map,
         #K=Kinship,
         #CV.GLM=Covariates,  ##if you have additional covariates, please keep there open.
         #CV.MLM=Covariates,
         #CV.FarmCPU=Covariates,
         #nPC.GLM=5,   ##if you have added PCs into covariates, please keep there closed.
         nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
         #nPC.FarmCPU=3,
         priority="speed",   ##for Kinship construction
         #ncpus=10,
         vc.method="BRENT",  ##only works for MLM
         maxLoop=10,
         #method.bin="static",   ## "FaST-LMM", "static" (#only works for FarmCPU)
         #permutation.threshold=TRUE,
         #permutation.rep=100,
         threshold=0.05,
         method=c("MLM")
       )
}else if(opt$method==3){
       imMVP <- MVP(
         phe=phenotype,
         geno=genotype,
         map=map,
         #K=Kinship,
         #CV.GLM=Covariates,  ##if you have additional covariates, please keep there open.
         #CV.MLM=Covariates,
         #CV.FarmCPU=Covariates,
         #nPC.GLM=5,   ##if you have added PCs into covariates, please keep there closed.
         #nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
         nPC.FarmCPU=3,
         priority="speed",   ##for Kinship construction
         #ncpus=10,
         #vc.method="BRENT",  ##only works for MLM
         maxLoop=10,
         method.bin="FaST-LMM",   ## "FaST-LMM", "static" (#only works for FarmCPU)
         #permutation.threshold=TRUE,
         #permutation.rep=100,
         threshold=0.05,
         method=c("FarmCPU")
       )

}else if(opt$method==4){
       imMVP <- MVP(
         phe=phenotype,
         geno=genotype,
         map=map,
         #K=Kinship,
         #CV.GLM=Covariates,  ##if you have additional covariates, please keep there open.
         #CV.MLM=Covariates,
         #CV.FarmCPU=Covariates,
         nPC.GLM=5,   ##if you have added PCs into covariates, please keep there closed.
         nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
         nPC.FarmCPU=3,
         priority="speed",   ##for Kinship construction
         #ncpus=10,
         vc.method="BRENT",  ##only works for MLM
         maxLoop=10,
         method.bin="static",   ## "FaST-LMM", "static" (#only works for FarmCPU)
         #permutation.threshold=TRUE,
         #permutation.rep=100,
         threshold=0.05,
         method=c("GLM", "MLM", "FarmCPU")
       )
}

     #save.image(paste(colnames(phenotype)[2],".result.RData",sep=""))
     #MVP.Report(imMVP, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),col=c("grey60","grey30"), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,chr.den.col=c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3"),bin.size=1e6,signal.col=c("red","green"),signal.cex=c(1,1),signal.pch=c(19,19),file="jpg",memo="",dpi=300)

