#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'vcf','v',1,'character',
	'maf','a',2,'numeric',
	'minDP','d',2,'numeric',
	'missing','s',2,'numeric',
	'output','o',1,'character',
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
	--vcf	the input vcf file
	--output	the output dir
	--maf
	--minDP
	--missing
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$vcf)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
if ( is.null(opt$maf)){ opt$maf=0.05 }
if ( is.null(opt$minDP)){ opt$minDP=5 }
if ( is.null(opt$missing){ opt$missing=0.7 }

library(rMVP)
library(CMplot)
setwd(opt$output);
vcf.file=opt$vcf
out.file="pop.filter"
maf=opt$maf
minDP=opt$minDP
missing=opt$missing
system(sprintf("vcftools --vcf %s --maf %s --minDP %s --max-missing %s --remove-indels --out %s --recode",
		vcf.file,maf,minDP,missing,out.file),ignore.stdout = FALSE)

MVP.Data(fileVCF="pop.filter.recode.vcf",
	fileKin=FALSE,
	filePC=FALSE,
	out="mvp.vcf"
	)
map <- read.table("mvp.vcf.geno.map" , head = TRUE)

CMplot(map,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
      main="illumilla_60K",file.output=TRUE,verbose=TRUE,width=9,height=6)
