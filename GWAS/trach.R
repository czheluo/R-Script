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


MVP.Data(fileHMP="pop.hapmap",
		filePhe="T1.trt",
		sep.hmp="\t",
		sep.phe="\t",
		SNP.effect="Add",
		fileKin=FALSE,
		filePC=FALSE,
		out="mvp.hmp",
		priority="memory",
	)


vcf.file="pop.recode.vcf"
out.file="pop.filter"
maf=0.05
minDP=5
missing=0.7
system(sprintf("vcftools --vcf %s --maf %s --minDP %s --max-missing %s --remove-indels --out %s --recode",
                   vcf.file,maf,minDP,missing,out.file),ignore.stdout = FALSE)

CMplot(map,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
    main="illumilla_60K",file.output=TRUE,verbose=TRUE,width=9,height=6)