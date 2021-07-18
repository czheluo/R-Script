#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## PCA, HEATMAP AND SCATTERPLOTS FOR SAMPLES IN COUNTS FILE
## - SAMPLE NAMES HAVE TO END IN e.g. "_R1" REPRESENTING REPLICATE ID. LAST 3 CHARACTERS OF SAMPLE NAME WILL BE TRIMMED TO OBTAIN GROUP ID FOR DESEQ2 COMPARISONS.
## - PACKAGES BELOW NEED TO BE AVAILABLE TO LOAD WHEN RUNNING R

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-f", "--count_col"     ), type="integer"  , default=2       , metavar="integer", help="First column containing sample count data."                                             ),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , metavar="integer", help="Column containing identifiers to be used."                                              ),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''      , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='deseq2', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$count_file)){
    print_help(opt_parser)
    stop("Please provide a counts file.", call.=FALSE)
}
times<-Sys.time()

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################




escaptime<-Sys.time()-times;
print("Done!");
print(escaptime)


