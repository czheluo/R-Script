
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]


tmtl<-as.matrix(read.csv("TMTL.csv",header = T,row.names = 1))

col = rev(rainbow(10))
cn = "log2(FC(UpTL/UpTM))"
pdf(file="TMTL.pdf",width = 6,height = 8)
Heatmap(tmtl,show_column_names = FALSE,cluster_rows = FALSE,
        heatmap_legend_param = list(title=" ",grid_height = unit(1.5, "cm"),grid_width = unit(0.7, "cm")),
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(cn, rot = 45, offset = unit(1, "npc"), just = "right"),
          annotation_height = max_text_width(cn)
        )
)
dev.off()

png(file="TMTL.png",width = 600,height = 800)
Heatmap(tmtl,show_column_names = FALSE,cluster_rows = FALSE,
        heatmap_legend_param = list(title=" ",grid_height = unit(1.5, "cm"),grid_width = unit(0.7, "cm")),
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(cn, rot = 45, offset = unit(1, "npc"), just = "right"),
          annotation_height = max_text_width(cn)
        )
)
dev.off()

tmtl<-as.matrix(read.csv("TMTH.csv",header = T,row.names = 1))
cn = "log2(FC(UpTH/UpTM))"
pdf(file="TMTH.pdf",width = 6,height = 8)
Heatmap(tmtl,show_column_names = FALSE,cluster_rows = FALSE,
        heatmap_legend_param = list(title=" ",grid_height = unit(1.5, "cm"),grid_width = unit(0.7, "cm")),
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(cn, rot = 45, offset = unit(1, "npc"), just = "right"),
          annotation_height = max_text_width(cn)
        )
)
dev.off()

png(file="TMTH.png",width = 600,height = 800)
Heatmap(tmtl,show_column_names = FALSE,cluster_rows = FALSE,
        heatmap_legend_param = list(title=" ",grid_height = unit(1.5, "cm"),grid_width = unit(0.7, "cm")),
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(cn, rot = 45, offset = unit(1, "npc"), just = "right"),
          annotation_height = max_text_width(cn)
        )
)
dev.off()
