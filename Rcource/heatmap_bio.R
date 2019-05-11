setwd("G:\\MAJORBIO\\R Course")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
.libPaths('F:\\R\\win-library\\3.4')
install.packages("ComplexHeatmap")
install.packages("circllize")
library(ComplexHeatmap)
library(circlize)
#install.packages('ComplexHeatmap')
#install.packages('circlize')
#source("http://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")

#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
#library('ComplexHeatmap')



expr = readRDS(paste0(system.file(package = "ComplexHeatmap"), "/extdata/gene_expression.rds"))
mat = as.matrix(expr[, grep("cell", colnames(expr))])
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))

type = gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(df = data.frame(type = type))

Heatmap(mat_scaled, name = "expression", km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        top_annotation = ha, top_annotation_height = unit(4, "mm"), 
        show_row_names = FALSE, show_column_names = FALSE) +
  Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, width = unit(5, "mm")) +
  Heatmap(expr$length, name = "length", col = colorRamp2(c(0, 1000000), c("white", "orange")),
          heatmap_legend_param = list(at = c(0, 200000, 400000, 60000, 800000, 1000000), 
                                      labels = c("0kb", "200kb", "400kb", "600kb", "800kb", "1mb")),
          width = unit(5, "mm")) +
  Heatmap(expr$type, name = "type", width = unit(5, "mm"))







