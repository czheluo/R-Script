
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

rm(list=ls())

# load plotting functions

source("plotting_functions.R")

# directories

results.dir <- "./"
figures.dir <- "./"

# files

design.file <- paste(results.dir, "group.txt", sep="")
taxonomy.file <- paste(results.dir, "tax.txt", sep="")
otu_table.file <- paste(results.dir, "asv_taxon.xls", sep="")

# load data

design <- read.table(design.file, header=F, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

# re-order data matrices

idx <- design$V1%in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$V1, colnames(otu_table))
otu_table <- otu_table[, idx]

library(edgeR)

# parameters

scale <- 1000               # scale for the ternary plots
alpha <- 0.05               # significance threshold
p.adj.method <- "fdr"       # FDR p-value adjustment method

# subset samples

#idx <- design$group %in% c("root", "soil", "rhizosphere") & 
       #design$genotype %in% c("gifu", "soil")
idx<- design$V2 %in% c("springCK", "springP")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset$V2
#groups <- droplevels(groups)

d <- DGEList(counts=otu_table_subset, group=groups)
d <- calcNormFactors(d)

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)
d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTagwiseDisp(d2, design.mat)

fit <- glmFit(d2, design.mat)

lrt_rhizo_root <- glmLRT(fit, contrast=c(0, 1))
de_rhizo_root <- decideTestsDGE(lrt_rhizo_root, adjust.method=p.adj.method, p.value=alpha)
rhizo_otus <- rownames(otu_table_subset)[abs(de_rhizo_root) ==1]

lrt_rhizo_root$table$padj<-p.adjust(lrt_rhizo_root$table$PValue, method = p.adj.method)

##enriched

# root_pvals <- lrt_root_rhizo$table
root_pvals <- lrt_rhizo_root$table

### get mean R.A. in gifu root samples

otu_table_norm_all <- apply(otu_table, 2, function(x) x / sum(x))
 
idx<- design$V2 %in% c("springCK")

root_gifu_samples <- design$V1[idx]
otu_table_norm_root_gifu <- otu_table_norm_all[, colnames(otu_table_norm_all) %in% root_gifu_samples]

### get mean R.A. in gifu rhizosphere samples
idx<- design$V2 %in% c("springP")

rhizo_gifu_samples <- design$V1[idx]
otu_table_norm_rhizo_gifu <- otu_table_norm_all[, colnames(otu_table_norm_all) %in% rhizo_gifu_samples]

### gifu root-enriched OTUs

### generate data frame for plotting

root_pvals$otu <- rownames(root_pvals)

# P values

root_pvals$neglogp <- -log10(root_pvals$PValue)

# enrichment status

root_pvals$enrichment <- root_pvals$otu %in% rhizo_otus

# log fold change

#idx <- root_pvals$logFC < 0
#root_pvals$neglogp[idx] <- 0

# order OTUs according to taxonomy

taxonomy <- taxonomy[order(taxonomy[, 3], taxonomy[, 4], taxonomy[, 5]), ]
idx <- rownames(taxonomy) %in% root_pvals$otu
taxonomy <- taxonomy[idx, ]

idx <- match(rownames(taxonomy), root_pvals$otu)
root_pvals <- root_pvals[idx, ]

root_pvals$tax <- taxonomy[, 5]
root_pvals$asv <- taxonomy[, 1]
root_pvals$g <- taxonomy[, 8]
#~ root_pvals$tax <- as.character(root_pvals$tax)
#~ idx <- sort(root_pvals$tax, index.return=T)$ix
#~ root_pvals <- root_pvals[idx, ]

# root_pvals$tax <- as.character(root_pvals$tax)
# idx <- match(taxonomy[, 1], root_pvals$otu)
# root_pvals <- root_pvals[idx, ]

root_pvals$otu <- factor(root_pvals$otu, levels=root_pvals$otu)

# relative abundances


ra <- apply(otu_table_norm_root_gifu, 1, mean)
names(ra)<-paste(c(1:dim(otu_table_norm_root_gifu)[1]),sep="")
root_pvals$ra <- ra[match(root_pvals$otu, names(ra))]
root_pvals$ra <- ra[root_pvals$otu]


# generate vector of colors
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

utax<-unique(root_pvals$tax)
colors<-as.data.frame(cbind(taxon=utax,colors=rcolor[1:length(utax)]))

#colors <- read.table("taxon_colors.txt", sep="\t", header=T)
taxon <- sort(unique(root_pvals$tax))
colors <- data.frame(taxon=taxon, colors=colors$color[match(taxon, colors$taxon)])
colors$colors <- as.character(colors$colors)
colors$colors[is.na(colors$colors)] <- "grey"
enriched_taxa <- unique(taxonomy[rownames(taxonomy) %in% c(rhizo_otus), 5])
colors$colors[!colors$taxon %in% enriched_taxa] <- "grey"

colors[which(!(colors[,1] %in% enriched_taxa)),2] <- "grey"



# generate vector of colors
library(RColorBrewer)
color <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
rcolor <- color[sample(1:length(color), length(color))]

utax<-unique(root_pvals$g)
colors<-as.data.frame(cbind(taxon=utax,colors=rcolor[1:length(utax)]))

#colors <- read.table("taxon_colors.txt", sep="\t", header=T)
taxon <- sort(unique(root_pvals$g))
colors <- data.frame(taxon=taxon, colors=colors$color[match(taxon, colors$taxon)])
colors$colors <- as.character(colors$colors)
colors$colors[is.na(colors$colors)] <- "grey"

enriched_taxa <- unique(taxonomy[rownames(taxonomy) %in% c(rhizo_otus), 8])

colors$colors[!colors$taxon %in% g$g] <- "grey"

colors[which(!(colors[,1] %in% g$g)),2] <- "grey"


# multiple testing correction thresholds

BF <- -log10(0.05 / dim(root_pvals)[1])
FDR <- min(root_pvals$neglogp[root_pvals$enrichment==TRUE])
 root_pvals$enrich<-"NA"
root_pvals[which(root_pvals$enrichment==TRUE),12]<-"significantly enriched"
root_pvals[which(!(root_pvals$enrichment==TRUE)),12]<-"n.s."


p1 <- ggplot(root_pvals, aes(x=otu, y=neglogp, color=g, size=ra, shape=enrich)) +
             geom_point(alpha=.8) +
	    geom_hline(yintercept=BF, linetype=2, color="lightgrey") +
             scale_color_manual(values=colors[,2]) +
             scale_shape_manual(values=c(21, 19)) +
             scale_size(breaks=c(0, 0.04, 0.08, 0.12, 0.16)) +
             #ylim(c(0, 1)) +
             labs(x="OTU", y="-log10(P)") +
             main_theme +
             theme(axis.ticks.x=element_blank(),
                   axis.text.x=element_blank(),
                   legend.position="bottom")+guides(color = "none")+labs(shape = "",size="relative abundance (%)")

#ggsave(paste("manhattan_spring1.pdf", sep=""), p1, width=10, height=4, useDingbats=F)


ggsave(paste("manhattan_springg.pdf", sep=""), p1, width=10, height=4, useDingbats=F)
write.csv(root_pvals,"spring.result.csv",quote=F,row.names=F)
save.image("spring.RData")



datas<-read.csv("test.csv",header=F)
p1 <- ggplot() +	    geom_rect(data=datas,aes(xmin = V1,
                xmax = V2,
               ymin = - Inf,
              ymax = Inf,fill =V3), alpha = 0.5)+
             geom_point(data=root_pvals, aes(x=otu, y=neglogp, color=tax, size=ra, shape=enrich),alpha=.8) +
             #geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
             #scale_color_manual(values=colors[,2]) +
            #scale_shape_manual(values=c(21, 19)) +
             #scale_size(breaks=c(0, 0.04, 0.08, 0.12, 0.16)) +
             #ylim(c(0, 1)) +
             labs(x="OTU", y="-log10(P)") +
             main_theme +
             theme(axis.ticks.x=element_blank(),
                   axis.text.x=element_blank(),
                   legend.position="bottom")+guides(color = "none")+labs(shape = "",size="relative abundance (%)")


ggsave(paste("manhattan_winter.pdf", sep=""), p1, width=10, height=4, useDingbats=F)



pdf("test.pdf",width=15,height=10)

p1

dev.off()























lrt_root_rhizo <- glmLRT(fit, contrast=c(-1, 1, 0))
lrt_root_soil <- glmLRT(fit, contrast=c(0, 1, -1))

de_root_rhizo <- decideTestsDGE(lrt_root_rhizo, adjust.method=p.adj.method, p.value=alpha)
de_root_soil <- decideTestsDGE(lrt_root_soil, adjust.method=p.adj.method, p.value=alpha)

root_otus <- rownames(otu_table_subset)[de_root_rhizo==1 & de_root_soil==1]
root_otus_gifu <- root_otus

lrt_soil_rhizo <- glmLRT(fit, contrast=c(-1, 0, 1))
lrt_soil_root <- glmLRT(fit, contrast=c(0, -1, 1))

de_soil_rhizo <- decideTestsDGE(lrt_soil_rhizo, adjust.method=p.adj.method, p.value=alpha)
de_soil_root <- decideTestsDGE(lrt_soil_root, adjust.method=p.adj.method, p.value=alpha)

soil_otus <- rownames(otu_table_subset)[de_soil_rhizo==1 & de_soil_root==1]

enriched_otus <- data.frame(root_enriched=rownames(otu_table) %in% root_otus,
                            rhizo_enriched=rownames(otu_table) %in% rhizo_otus,
                            soil_enriched=rownames(otu_table) %in% soil_otus)
row.names(enriched_otus) <- rownames(otu_table)

write.table(enriched_otus, paste(results.dir, "gifu_root_rhizo_soil_OTUs.txt", sep=""),
            quote=F, sep="\t", col.names=T, row.names=T)

### ternary plots

# normalize subsetted OTU table and apply log transform

otu_table_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))
otu_table_norm_log <- log2(otu_table_norm + 1)

# create vectors of mean reltive abundances

idx <- design_subset$compartment=="rhizosphere"
rhizo_means <- apply(otu_table_norm[, idx], 1, mean)

idx <- design_subset$compartment=="root"
root_means <- apply(otu_table_norm[, idx], 1, mean)

idx <- design_subset$compartment=="soil"
soil_means <- apply(otu_table_norm[, idx], 1, mean)

# create matrix of average r.a. per group

df <- data.frame(root=root_means, rhizosphere=rhizo_means, soil=soil_means)
df <- df[rowSums(df)!=0, ]
df <- log2(df * scale + 1)

# sort the rows by decreasing abundance (looks better)

idx <- sort(rowSums(df), decreasing=F, index.return=T)$ix
df <- df[idx, ]

# create vector of colors according to enrichment

colors <- rep(c_grey, dim(df)[1])
colors[rownames(df) %in% root_otus] <- c_very_dark_green
colors[rownames(df) %in% soil_otus] <- c_dark_brown
colors[rownames(df) %in% rhizo_otus] <- c_dark_red

idx <- sort(colors==c_grey, decreasing=T, index.return=T)$ix
df <- df[idx, ]
colors <- colors[idx]

# plot colored by enrichment

pdf(file=paste(figures.dir, "gifu_root_rhizo_soil_enrichment.pdf", sep=""))

tern_e(df, prop_size=T, col=colors, grid_color="grey",
       labels_color="transparent", pch=19, main="gifu")

dev.off()

### boxplots of aggregated relative abundances

idx <- design_subset$compartment=="rhizosphere"

rhizo_rhizo <- colSums(otu_table_norm[rownames(otu_table_norm) %in% rhizo_otus, idx])
root_rhizo <- colSums(otu_table_norm[rownames(otu_table_norm) %in% root_otus, idx])
soil_rhizo <- colSums(otu_table_norm[rownames(otu_table_norm) %in% soil_otus, idx])

df_rhizo <- rbind(data.frame(otus="rhizo", compartment="rhizo", ra=rhizo_rhizo),
                  data.frame(otus="root", compartment="rhizo", ra=root_rhizo),
                  data.frame(otus="soil", compartment="rhizo", ra=soil_rhizo))

idx <- design_subset$compartment=="root"

rhizo_root <- colSums(otu_table_norm[rownames(otu_table_norm) %in% rhizo_otus, idx])
root_root <- colSums(otu_table_norm[rownames(otu_table_norm) %in% root_otus, idx])
soil_root <- colSums(otu_table_norm[rownames(otu_table_norm) %in% soil_otus, idx])

df_root <- rbind(data.frame(otus="rhizo", compartment="root", ra=rhizo_root),
                 data.frame(otus="root", compartment="root", ra=root_root),
                 data.frame(otus="soil", compartment="root", ra=soil_root))

idx <- design_subset$compartment=="soil"

rhizo_soil <- colSums(otu_table_norm[rownames(otu_table_norm) %in% rhizo_otus, idx])
root_soil <- colSums(otu_table_norm[rownames(otu_table_norm) %in% root_otus, idx])
soil_soil <- colSums(otu_table_norm[rownames(otu_table_norm) %in% soil_otus, idx])

df_soil <- rbind(data.frame(otus="rhizo", compartment="soil", ra=rhizo_soil),
                 data.frame(otus="root", compartment="soil", ra=root_soil),
                 data.frame(otus="soil", compartment="soil", ra=soil_soil))

df <- rbind(df_rhizo, df_root, df_soil)

df$compartment <- factor(df$compartment, levels=c("soil", "rhizo", "root"))
df$otus <- factor(df$otus, levels=c("soil", "rhizo", "root"))

p1 <- ggplot(df, aes(x=compartment, y=ra, color=otus)) +
             geom_boxplot(alpha=1, outlier.size=0, size=0.6, width=0.8,
                          position=position_dodge(width=0.7)) +
             scale_color_manual(values=c(c_dark_brown, c_dark_red, c_very_dark_green)) +
             scale_y_continuous(labels=percent, limits=c(0, 1)) +
             labs(x="", y="relative abundance") +
             main_theme +
             theme(legend.position="none")

ggsave(paste(figures.dir, "boxplots_gifu_root_rhizo_soil.pdf", sep=""), p1, width=5, height=3)

