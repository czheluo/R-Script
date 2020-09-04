fpkm<-read.table("unigene.TMM.fpkm.matrix",header=T,sep="\t")
fpkm[1,]
fpkm<-read.table("unigene.TMM.fpkm.matrix",header=T,sep="\t",rownames=1)
fpkm<-read.table("unigene.TMM.fpkm.matrix",header=T,sep="\t",row.names=1)
fpkm[1,]
fpkm$sum<-with(fpkm,mh_1+mh_2+mh_3+mkh_1+mkh_2+mkh_3)
fpkm[1,]
fpkm_order<-fpkm[order(-sum),]
fpkm_order<-fpkm[order(-fpkm$sum),]
fpkm_order[1,]
write.table(fpkm_order[1:5000,],"fpkm_order.xls",quote=F,colnames=T,rownames=T)
write.table(fpkm_order[1:5000,],"fpkm_order.xls",quote=F,col.names=T,row.names=T)
savehistory()
savehistory("fpkm_order.txt")
quit()
j<-read.table("all.xls",sep="\t",header=T,row.names=1)
j[1,]
j$sum<-with(j,mh_1+mh_2+mh_3+mkh_1+mkh_2+mkh_3)
a<-j[order(-j$sum),]
a[1,]
write.table(a,"all.a.xls",sep="\t",col.names=T,row.names=T,quote=F)
savehistory("a.all.txt")
quit()
al<-read.table("all.xls",sep="\t",head=T,row.names=1)
al$sum<-with(al,sum(colnames(al)))
colnames(al)
(al)
length(colnames(al))
al$sum<-with(al,sum(al[,1:length(colnames(al))]))
al[1,]
al[2,]
al$sum<-colSums(al)
al[2,]
al[1,]
al$sum<-rowSums(al)
al[1,]
al<-al[,-7]
al[1,]
al$sum<-rowSums(al)
al[1,]
new<-al[order(-al$sum),]
new[1:20,]
new<-al[order(-al$sum),][1:5000,]
write.table(new,"new.xls",sep="\t",head=T,quote=F,colnames=T,rownames=F)
write.table(new,"new.xls",sep="\t",head=T,quote=F,col.names=T,row.names=F)
write.table(new,"new.xls",sep="\t",quote=F,col.names=T,row.names=F)
savehistory("get.r")
quit()
library(DESeq)
count<-read.table('gene_count.table.xls',sep="\t",row.names=1)
head (count)
count<-read.table('gene_count.table.xls',sep="\t",row.names=1)
head (count)
count<-read.table('gene_count.table.xls',row.names=1,sep="\t")
head (count)
count[1,]
count<-read.table('gene_count.table.xls',row.names=2,sep="\t")
count<-read.table('gene_count.table.xls',row.names=0,sep="\t")
count<-read.table('gene_count.table.xls',sep="\t")
head (count)
count<-read.table('gene_count.table.xls',row.names=1,sep="\t")
head (count)
count<-read.table('gene_count.table.xls',row.names=1,col.names=1,sep="\t")
count<-read.table('gene_count.table.xls',row.names=1,col.names=2,sep="\t")
count<-read.table('gene_count.table.xls',row.names=1,sep="\t")
head (count)
a<-as.data.frame(count)
head (a)
rm(a)
count<-read.table('gene_count.table.xls',row.names=-1,sep="\t")
count<-read.table('gene_count.table.xls',row.names=1,sep="\t")
head (count)
count<-read.table('gene_count.table.xls',sep="\t")
head (count)
count<-read.table('gene_count.table.xls',row.names=1,sep="\t")
head (count)
count<-read.table('gene_count.table.xls',row.names=1,sep="\t",header=T)
head (count)
head (counts(count))
counts(count)
condition<-factor(c("A", "A", "A","B","B", "B"))
condition<-factor(c("S", "S", "S","D","D", "D"))
condition
cds<-newCountDataSet(count, condition)
cds
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
res<-nbinomTest(cds,"A","B") 
res<-nbinomTest(cds,"D","S") 
res
head (res)
condition
head (counts(cds))
head (counts(cds,normalized=TRUE))
nor_count (counts(cds,normalized=TRUE))
nor_count<- (counts(cds,normalized=TRUE))
head (nor_count)
write.table(nor_count,file="count.txt",sep="\t")
write.table(nor_count,file="count.txt",sep="\t",row.names=1)
write.table(nor_count,file="count.txt",sep="\t",header=T)
write.table(nor_count,file="count.txt")
write.table(nor_count,file="count.txt",sep="\t")
q()
head (res)
write.table(res,file=ref)
write.table(res,file="res.txt",sep="\t")
q()
dir()
GO<-read.delim("gene.GO.list",header=F,sep="\t")
GO[1:4,]
Gene1<-as.character(GO[[1]])
table(Gene1)[1:4]
q()
dir()
t2g<-read.delim("isoform_gene_IdMapping",sep="\t",header=F)
t2g[1:4,]
t2g<-read.delim("isoform_gene_IdMapping",sep="\t",header=T,check.names=F)
t2g[1:4,]
allGO<-read.delim("GO.list",header=F,sep="\t")
names(allGO)
names(allGO)<-c("tracking_id","GOs")
merge(allGO,t2g,by="tracking_id")->t2GO2g
t2GO2g[1:4,]
g2GOs<-unique(t2GO2g[,c("gene_id","GOs")])
dim(g2GOs)
g2GOs[1:10,]
allGenes<-as.character(g2GOs[[1]])
table(allGenes)[1:4]
table(allGenes)[1:40]
table_result<-table(allGenes)
table_result[1:10]
cbind(names(table_result),table_result)->ll
ll[1:40,]
dim(ll)
ll[1:4,]
table_result<-sort(table(allGenes))
table_result[1:10]
?sort
table_result<-sort(table(allGenes),decreasing=T)
table_result[1:10]
table_result<-table(allGenes)
table_result[1:4,]
table_result[1:4]
g2GOs[1:4,]
which(table_result==1)
which(table_result==1)->uniq_lines
uniq_lines[1:10]
table_result[1:10]
which(table_result>1)->multi_lines
multi_lines[1:10]
g2GOs[multi_lines,][1:4,]
uniq_lines_g2GOs<-g2GOs[uniq_lines,]
uniq_lines_g2GOs[1:4,]
multi_lines_g2GOs<-g2GOs[multi_lines,]
multi_lines_g2GOs[1:10,]
dim(multi_lines_g2GOs)
ls()
table_result[1:4]
Genes_uniq<-names(table_result)[uniq_lines]
Genes_uniq[1:4,]
Genes_uniq[1:4]
Genes_multi<-names(table_result)[multi_lines]
Genes_multi[1:4]
Genes_uniq[1:4,]
Genes_uniq[1:4]
Genes_uniq[1:40]
Genes_uniq[100:150]
Genes_uniq_1<-as.data.frame(Genes_uniq)
names(Genes_uniq_1)
names(Genes_uniq_1)<-"gene_id"
Genes_multi_1<-as.data.frame(Genes_multi)
Genes_multi_1[1:4,]
names(Genes_multi_1)<-"gene_id"
g2GOs[1:4,]
merge(Genes_uniq_1,g2GOs,by="gene_id")->Genes_uniq_1_GOs
merge(Genes_multi_1,g2GOs,by="gene_id")->Genes_multi_1_GOs
dim(Genes_uniq_1_GOs)
dim(Genes_multi_1_GOs)
Genes_multi_1_GOs[1:4,]
Genelist<-unique(as.character(Genes_multi_1_GOs[[1]]))
length(Genelist)
Genes_multi_1_GOs_new<-data.frame()
for(i in 1:length(Genelist)){
gos<-Genes_multi_1_GOs[which(Genes_multi_1_GOs[["gene_id"]]==Genelist[i]),"GOs"]
gos_1<-paste(gos,collapse=";")
gos_2<-unique(unlist(strsplit(gos_1,";",fix=T)))
Genes_multi_1_GOs_new[i,1]<-Genelist[i]
Genes_multi_1_GOs_new[i,2]<-gos_2
}
Genes_multi_1_GOs[which(Genes_multi_1_GOs[["gene_id"]]==Genelist[1]),"GOs"]
as.character(Genes_multi_1_GOs[which(Genes_multi_1_GOs[["gene_id"]]==Genelist[1]),"GOs"])
Genes_multi_1_GOs_new<-data.frame()
for(i in 1:length(Genelist)){
gos<-as.character(Genes_multi_1_GOs[which(Genes_multi_1_GOs[["gene_id"]]==Genelist[i]),"GOs"])
gos_1<-paste(gos,collapse=";")
gos_2<-unique(unlist(strsplit(gos_1,";",fix=T)))
Genes_multi_1_GOs_new[i,1]<-Genelist[i]
Genes_multi_1_GOs_new[i,2]<-gos_2
}
gos<-Genes_multi_1_GOs[which(Genes_multi_1_GOs[["gene_id"]]==Genelist[1]),"GOs"]
gos
gos<-as.character(Genes_multi_1_GOs[which(Genes_multi_1_GOs[["gene_id"]]==Genelist[1]),"GOs"])
gos
paste(gos,collapse=";")
gos_1<-paste(gos,collapse=";")
gos_1
gos_2<-unique(unlist(strsplit(gos_1,";",fix=T)))
gos_2
Genes_multi_1_GOs_new<-data.frame()
for(i in 1:length(Genelist)){
gos<-as.character(Genes_multi_1_GOs[which(Genes_multi_1_GOs[["gene_id"]]==Genelist[i]),"GOs"])
gos_1<-paste(gos,collapse=";")
gos_2<-unique(unlist(strsplit(gos_1,";",fix=T)))
Genes_multi_1_GOs_new[i,1]<-Genelist[i]
Genes_multi_1_GOs_new[i,2]<-paste(gos_2,collapse=";")
}
Genes_multi_1_GOs_new[1:4,]
dim(Genes_multi_1_GOs_new)
dim(Genes_uniq_1_GOs)
names(Genes_uniq_1_GOs)
names(Genes_multi_1_GOs_new)
names(Genes_multi_1_GOs_new)<-c("gene_id","GOs")
rbind(Genes_uniq_1_GOs,Genes_multi_1_GOs_new)->g2GOs_uniq
dim(g2GOs_uniq)
q()
write.table(g2GOs_uniq,"gene.GO.list.uniq",sep="\t",col.names=F,row.names=F,quote=F)
q()
dir()
pathway<-read.delim("pathway.txt.new",sep="\t",header=T,check.names=F)
pathway[1:4,]
pathway[1:4,1]
pathway<-read.delim("pathway.txt.new",sep="\t",header=T,check.names=F)
pathway[1:4,]
dir()
t2g<-read.delim("isoform_gene_IdMapping",sep="\t",header=T,check.names=F)
t2g[1:4,]
names(t2g)[1]<-"isoform_id"
merge(pathway,t2g,by="isoform_id")->t2g2pathway
t2g2pathway[1:4,]
g2KOs<-unique(t2g2pathway[,c("gene_id","KO")])
g2KOs[1:4,]
g2KOs[1:40,]
length(as.character(g2KOs[[1]]))
length(unique(as.character(g2KOs[[1]])))
g2KOs[1:4,]
write.table(g2KOs,"gene.pathway",sep="\t",col.names=F,row.names=F,quote=F)
q()
ls()
dim(g2KOs)
g2KOs[1:4,]
dir()
Genes_multi_pathway<-read.delim("Genes_multi_pathway",sep="\t",header=F)
Genes_multi_pathway[1:4,]
Genes_multi_pathway<-read.delim("Genes_multi_pathway",sep=" ",header=F)
Genes_multi_pathway[1:4,]
Genes_multi_pathway<-read.delim("Genes_multi_pathway.new",sep=" ",header=F)
Genes_multi_pathway[1:4,]
Genes_multi_pathway<-read.delim("Genes_multi_pathway.new",sep="\t",header=F)
Genes_multi_pathway[1:4,]
Nonelines<-vector()
for(i in nrow(Genes_multi_pathway)){
}
ls()
Genes_multi_pathway<-read.delim("Genes_multi_pathway.new.new",sep="\t",header=F)
Genes_multi_pathway
for(i in 1:nrow(g2KOs)){
}
for(i in 1:nrow(Genes_multi_pathway)){
lines_1<-which(g2KOs[[1]]==Genes_multi_pathway[i,2])
which(g2KOs[[1]]==Genes_multi_pathway[1,2]
)
which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2])
Genes_multi_pathway[which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2]),2]
g2KOs[which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2]),2]
g2KOs[which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2]),]
g2KOs[1:4,]
NoneLines<-which(g2KOs[[2]]=="None")
NoneLines[1:10]
g2KOs[1:10,]
g2KOs[intersect(which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2]),NoneLines)]
which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2])
NoneLines[1:10]
intersect((which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2])),NoneLines)
g2KOs[intersect((which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2])),NoneLines)]
g2KOs[intersect((which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2])),NoneLines),]
for(i in 1:nrow(Genes_multi_pathway)){
Nonelines[i]<-intersect((which(as.character(g2KOs[[1]])==Genes_multi_pathway[1,2])),NoneLines)
}
Nonelines
Nonelines<-vector()
for(i in 1:nrow(Genes_multi_pathway)){
Nonelines[i]<-intersect((which(as.character(g2KOs[[1]])==Genes_multi_pathway[i,2])),NoneLines)
}
i
Genes_multi_pathway[i,2]
intersect((which(as.character(g2KOs[[1]])==Genes_multi_pathway[i,2])),NoneLines)
history()
Nonelines<-vector()
for(i in 1:nrow(Genes_multi_pathway)){
line1<-intersect((which(as.character(g2KOs[[1]])==Genes_multi_pathway[i,2])),NoneLines)
if(length(line1)!=0){Nonelines[i]<-line1}
if(length(line1)==0){next}
}
Nonelines
Nonelines<-vector()
flag<-1
for(i in 1:nrow(Genes_multi_pathway)){
line1<-intersect((which(as.character(g2KOs[[1]])==Genes_multi_pathway[i,2])),NoneLines)
if(length(line1)!=0){Nonelines[flag]<-line1;flag<-flag+1}
if(length(line1)==0){next}
}
Nonelines
g2KOs[Nonelines,]
g2KOs_new<-g2KOs[-Nonelines,]
g2KOs_new[1:40,]
write.table(g2KOs_new,"gene2KOs",sep="\t",col.names=F,row.names=F,quote=F)
q()
R
library(cummeRbund)
cuff<-readCufflinks()
cuff<-readCufflinks()
cuff
q()
library(cummeRbund)
cuff<-readCufflinks()
cuff<-readCufflinks()
cuff
dens<-csDensity(genes(cuff))
dens
savePlot(filename = "FPKM_distribution_gene.png",type = "png",device = dev.cur())
b<-csBoxplot(genes(cuff))
b
savePlot(filename = "Boxplot_gene.png",type = "png",device = dev.cur())
#以下需两两样本执行samples(cuff)
samples(cuff)
s12<-csScatter(genes(cuff),x="X3960E_SAM_DEG",y="YB_SAM_DEG",smooth=T)
#样品名如果是以数字开头，软件会自动在样品名前加
Xs12
savePlot(filename = "X3960E_vs_YB.png",type = "png",device = dev.cur())
S12<-csScatter(genes(cuff),x="X3960E_SAM_DEG",y="YB_SAM_DEG",smooth=T)
#样品名如果是以数字开头，软件会自动在样品名前加
XS12
savePlot(filename = "X3960E_vs_YB.png",type = "png",device = dev.cur())
ss<-csScatter(genes(cuff),x="X3960E_SAM_DEG",y="YB_SAM_DEG",smooth=T)
ss
savePlot(filename = "X3960E_vs_YB.png",type = "png",device = dev.cur())
ss<-csScatter(genes(cuff),x="X3960E_SAM_DEG",y="YB_SAM_DEG",smooth=T)
ss
savePlot(filename = "X3960E_vs_YB_Compare.png",type = "png",device = dev.cur())
v12<-csVolcano(genes(cuff),"X3960E_SAM_DEG","YB_SAM_DEG")
v12
savePlot(filename = "X3960E_vs_YB_Volcano.png",type = "png",device = dev.cur())
dens_isoform<-csDensity(isoforms(cuff))
dens_isoform
savePlot(filename = "FPKM_distribution_isoform.png",type = "png",device = dev.cur())
b_isoform<-csBoxplot(isoforms(cuff))
b_isoform
savePlot(filename = "Boxplot_isoform.png",type = "png",device = dev.cur())
ss<-csScatter(isoforms(cuff),x="X3960E_SAM_DEG",y="YB_SAM_DEG",smooth=T)
ss
savePlot(filename = "X3960E_vs_YB_Compare_isoform.png",type = "png",device = dev.cur())
v12<-csVolcano(isoforms(cuff),"X3960E_SAM_DEG","YB_SAM_DEG")
v12
savePlot(filename = "X3960E_vs_YB_Volcano_isoform.png",type = "png",device = dev.cur())
gene.diff<-read.delim("gene_exp.diff",header=TRUE,sep="\t")
gene.diff.sig<-gene.diff[which(gene.diff[["significant"]]=="yes"),]
print(dim(gene.diff.sig)[1])
gene.diff.sig[1:4,]
gene.diff.sig[1:100,]
DEgene<-as.character(gene.diff.sig[,"gene"])
write.table(gene.diff.sig,"3960E_vs_YB_DEGinfo_gene.xls",sep="\t",row.names=F)
write.table(DEgene,"3960E_vs_YB_DEGs_gene.txt",sep="\n",col.names=F,row.names=F,quote=F)
isoform.diff<-read.delim("isoform_exp.diff",header=TRUE,sep="\t")
isoform.diff.sig<-isoform.diff[which(isoform.diff[["significant"]]=="yes"),]
print(dim(isoform.diff.sig)[1])
isoform.diff.sig[1:20,]
DEisoform<-as.character(isoform.diff.sig[,"test_id"])
write.table(isoform.diff.sig,"3960E_vs_YB_DEGinfo_isoform.xls",sep="\t",row.names=F)
write.table(DEisoform,"3960E_vs_YB_DEGs_isoform.txt",sep="\n",col.names=F,row.names=F,quote=F)
length(DEisoform)
length(DEgene)
dim(isoform.diff)
dim(gene.diff)
print(dim(isoform.diff.sig)[1])
print(dim(gene.diff.sig)[1])
q()
miRNA_Exp_Matrix_nonZero<-read.delim("miRNA_3960E_YB_nonZero",sep="\t",header=F)
gene_Exp_Matrix_nonZero<-read.delim("geneID_3960E_YB_nonZero",sep="\t",header=T)
miRNA_target_map<-read.delim("miranda_miRNA_target_map",sep="\t",header=T)
miRNA_target_map<-read.delim("miranda_miRNA_target_map",sep="\t",header=F)
names(miRNA_Exp_Matrix_nonZero)<-c("miRNA","3960E_Normal","YB_Normal")
names(gene_Exp_Matrix_nonZero)<-c("geneName","3960E_FPKM","YB_FPKM")
names(miRNA_target_map)<-c("miRNA","geneName")
miR_Tar_Mmir<-merge(miRNA_target_map,miRNA_Exp_Matrix_nonZero,by="miRNA")
miR_Tar_Mmir[1:2,]
miR_Tar_MM<-merge(miR_Tar_Mmir,gene_Exp_Matrix_nonZero,by="geneName")
miR_Tar_MM[1:4,]
miR_Tar_Mmir
ls()
names(miR_Tar_Mmir)
dim(miR_Tar_Mmir)
names(gene_Exp_Matrix_nonZero)
dim(gene_Exp_Matrix_nonZero)
class(gene_Exp_Matrix_nonZero)
class(miR_Tar_Mmir)
miR_Tar_MM<-merge(miR_Tar_Mmir,gene_Exp_Matrix_nonZero,by="geneName")
dim(miR_Tar_MM)
miR_Tar_Mmir[1:10,]
gene_Exp_Matrix_nonZero[1:4,]
gene_Exp_Matrix_nonZero[100:105,]
isoform_gene_IdMapping<-read.delim("isoform_gene_IdMapping",sep="\t",header=T)
names(miRNA_Exp_Matrix_nonZero)<-c("miRNA","3960E_Normal","YB_Normal")
names(gene_Exp_Matrix_nonZero)<-c("geneName","3960E_FPKM","YB_FPKM")
names(miRNA_target_map)<-c("miRNA","geneName")
names(isoform_gene_IdMapping)<-c("isoformID","geneID")
names(isoform_gene_IdMapping)<-c("isoformID","geneName")
miRNA_target_map_new<-merge(miRNA_target_map,isoform_gene_IdMapping,by="geneName")
miRNA_target_map_new[1:4,]
names(miRNA_target_map)<-c("miRNA","isoform")
names(isoform_gene_IdMapping)<-c("isoform","geneName")
miRNA_target_map_new<-merge(miRNA_target_map,isoform_gene_IdMapping,by="isoform")
miRNA_target_map_new[1:10,]
miRNA_target_map_new[1000:1050,]
miRNA_target_map_new[5000,]
miRNA_target_map_new[8000,]
miRNA_target_map_new[9000,]
miRNA_target_map_new[4000,]
dim(miRNA_target_map_new)
miRNA_target_map_new[10000,]
miRNA_target_map_new[50000,]
miR_Tar_Mmir_new<-miR_Tar_Mmir[,c("miRNA","geneName")]
miR_Tar_MM<-merge(miR_Tar_Mmir_new,gene_Exp_Matrix_nonZero,by="geneName")
miR_Tar_MM[1:4,]
miR_Tar_Mmir_new[1:4,]
dim(miR_Tar_Mmir_new)
miR_Tar_Mmir_new[40000,]
miR_Tar_Mmir_new[50000,]
miR_Tar_Mmir[40000:40010,]
miRNA_target_map_new_new<-miRNA_target_map_new[,c("miRNA","geneName")]
miRNA_target_map_new_new[40000,]
miRNA_target_map_new_new[40000:40010,]
miR_Tar_Mmir<-merge(miRNA_target_map_new_new,miRNA_Exp_Matrix_nonZero,by="miRNA")
miR_Tar_Mmir[1:4,]
miR_Tar_MM<-merge(miR_Tar_Mmir,gene_Exp_Matrix_nonZero,by="geneName")
dim(miR_Tar_MM)
miR_Tar_MM[1:4,]
miR_Tar_MM[1:4,]
names(miR_Tar_MM)
miR_Tar_MM_new<-miR_Tar_MM[,c("miRNA","geneName","3960E_Normal","YB_Normal","3960E_FPKM","YB_FPKM")]
miR_Tar_MM_new[1:20,]
miR_Tar_MM_new_new<-miR_Tar_MM_new[order(miR_Tar_MM_new[["miRNA"]],decreasing=FALSE),]
miR_Tar_MM_new_new[1:10,]
miR_Tar_MM_new<-unique(miR_Tar_MM[,c("miRNA","geneName","3960E_Normal","YB_Normal","3960E_FPKM","YB_FPKM")])
miR_Tar_MM_new_new<-miR_Tar_MM_new[order(miR_Tar_MM_new[["miRNA"]],decreasing=FALSE),]
miR_Tar_MM_new_new[1:40,]
write.table(miR_Tar_MM_new_new,"miR_Tar_MM.txt",sep="\t",col.names=TRUE,row.names=F,quote=F)
miRNAs<-as.character(miR_Tar_MM_new_new[["miRNA"]])
aa<-table(miRNAs)
aa[1:10,]
aa[1:10]
miRNA_tarNum<-as.data.frame(cbind(names(aa)),aa)
miRNA_tarNum[1:10,]
miRNA_tarNum<-as.data.frame(cbind(names(aa),as.numeric(aa)))
miRNA_tarNum[1:10,]
miRNA_tarNum_new<-miRNA_tarNum[order(as.numeric(miRNA_tarNum[[2]]),decreasing=TRUE),]
miRNA_tarNum_new[1:10,]
as.numeric(miRNA_tarNum[[2]]
)
DIM(miRNA_tarNum)
dim(miRNA_tarNum)
miRNA_tarNum
miRNA_tarNum_new<-miRNA_tarNum[order(as.numeric(miRNA_tarNum[,2]),decreasing=TRUE),]
miRNA_tarNum_new[1:10,]
names(miRNA_tarNum)<-c("miRNA","targetNumber")
write.table(miRNA_tarNum,"miRNA_tarNum.txt",sep="\t",col.names=TRUE,row.names=F,quote=F)
q()
ComputeCorrForMiM<-function(MM_matrix,method,OutPutFileName){
   CCbox<-vector()
   Pbox<-vector()
   N1<-dim(MM_matrix)[1]
   N2<-dim(MM_matrix)[2]
   miRNA_col<-3:(((N2-2)/2)+2)
   mRNA_col<-(((N2-2)/2)+3):N2
   for(i in 1:N1){
vector_miRNA<-as.numeric(MM_matrix[i,miRNA_col])
vector_mRNA<-as.numeric(MM_matrix[i,mRNA_col])
if(method=="pearson"){
    if(length(vector_miRNA)<3){print("No enough values for caculating the pvalues,please input at least 3 pairs of values!");break}
else{CorrCoeffcient<-cor.test(vector_miRNA,vector_mRNA,method = c("pearson"))}
}
        if(method=="spearman"){CorrCoeffcient<-cor.test(vector_miRNA,vector_mRNA,method = c("spearman"))}
if(method=="kendall"){CorrCoeffcient<-cor.test(vector_miRNA,vector_mRNA,method = c("kendall"))}
CCbox[i]<-CorrCoeffcient$estimate
Pbox[i]<-CorrCoeffcient$p.value
print(i)
    }
FDRbox<-p.adjust(as.numeric(Pbox),"BH") ## "BH":Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
ResultBox<-cbind(MM_matrix,CCbox,Pbox,FDRbox)
    colnames(ResultBox)<-c(names(MM_matrix),"CorrC","P_value","FDR")
ResultBox<-ResultBox[order(ResultBox[["FDR"]],decreasing=FALSE),]
    write.table(ResultBox,OutPutFileName,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
    return(ResultBox) 
}
ls()
Result_CC<-ComputeCorrForMiM(miR_Tar_MM_new_new,method="pearson",OutPutFileName="test")
dim(miR_Tar_MM_new_new)
Result_CC<-ComputeCorrForMiM(miR_Tar_MM_new_new,method="spearman",OutPutFileName="test")
ls -rtl
ls()
test_matrix<-miR_Tar_MM_new_new[1:1000,]
Result_CC<-ComputeCorrForMiM(test_matrix,method="spearman",OutPutFileName="test")
ls()
Result_CC[1:10,]
Result_CC[950:1000,]
q()
miRNA_target_map<-read.delim("target_predict",sep="\t",header=F)
miRNA_target_map_new<-data.frame()
flag<-1
for(i in 1:dim(miRNA_target_map)[1]){
    print(i)
    miRNA<-as.character(miRNA_target_map[i,1])
    targetList<-unique(unlist(strsplit(as.character(miRNA_target_map[i,3]),";",fixed=F)))
    N<-length(targetList)
    miRNA_target_map_new[flag:(flag+N-1),1]<-rep(miRNA,N)
    miRNA_target_map_new[flag:(flag+N-1),2]<-targetList
    flag<-flag+N
}
names(miRNA_target_map_new)<-c("miRNA","isoform")
dim(miRNA_target_map_new)
miRNA_target_map_new[1:4,]
miRNA_target_map_new_new<-unique(miRNA_target_map_new)
dim(miRNA_target_map_new_new)
write.table(miRNA_target_map_new_new,"miRNA_target_map.txt",sep="\t",col.names=TRUE,row.names=F,quote=F)
miRNA_Exp_Matrix_nonZero<-read.delim("known_and_novo_miR_expMatrix_nonZero",sep="\t",header=F)
gene_Exp_Matrix_nonZero<-read.delim("geneID_3960E_YB_nonZero",sep="\t",header=T)
miRNA_target_map<-read.delim("miRNA_target_map.txt",sep="\t",header=T)
isoform_gene_IdMapping<-read.delim("isoform_gene_IdMapping",sep="\t",header=T)
names(miRNA_Exp_Matrix_nonZero)<-c("miRNA","3960E_Normal","YB_Normal")
names(gene_Exp_Matrix_nonZero)<-c("geneName","3960E_FPKM","YB_FPKM")
names(miRNA_target_map)<-c("miRNA","isoform")
names(isoform_gene_IdMapping)<-c("isoform","geneName")
miRNA_target_map_new<-merge(miRNA_target_map,isoform_gene_IdMapping,by="isoform")
miRNA_target_map_new_new<-miRNA_target_map_new[,c("miRNA","geneName")]
miR_Tar_Mmir<-merge(miRNA_target_map_new_new,miRNA_Exp_Matrix_nonZero,by="miRNA")
miR_Tar_MM<-merge(miR_Tar_Mmir,gene_Exp_Matrix_nonZero,by="geneName")
miR_Tar_MM_new<-unique(miR_Tar_MM[,c("miRNA","geneName","3960E_Normal","YB_Normal","3960E_FPKM","YB_FPKM")])
miR_Tar_MM_new_new<-miR_Tar_MM_new[order(miR_Tar_MM_new[["miRNA"]],decreasing=FALSE),]
write.table(miR_Tar_MM_new_new,"miR_Tar_MM.txt",sep="\t",col.names=TRUE,row.names=F,quote=F)
miR_Tar_MM_new_new[1:4,]
q()
dir()
gf_1<-read.delim("3960E_vs_YB_diff_exp_result_gene",header=T,sep="\t",check.names=F)
gf_2<-read.delim("MIM_vs_YB_diff_exp_result_gene",header=T,sep="\t",check.names=F)
gf_12<-merge(gf_1,gf_2,by="gene_id")
dim(gf_1)
dim(gf_2)
dim(gf_12)
names(gf_12)
gf_12_new<-gf_12[,-14]
names(gf_12_new)
names(gf_12_new)[2]
names(gf_12_new)[2]<-"locus"
names(gf_12_new)
gf_12_new<-gf_12[,-c(14,9,11,21,23)]
names(gf_12_new)
names(gf_12_new)[2]<-"locus"
names(gf_12_new)
gf_12_new[1:10,]
names(gf_12_new)
names(gf_12)
gf_12_new<-gf_12[,-c(14,9,11,22,24)]
gf_12_new[1:10,]
write.table(gf_12_new,"TwoTimes_expdiff",sep="\t",col.names=T,row.names=F,quote=F)
write.table(gf_12_new,"TwoTimes_expdiff_noTitle",sep="\t",col.names=F,row.names=F,quote=F)
names(gf_12_new)
all_F<-as.character(read.delim("all_1",header=F)[[1]])
length(all_F)
all_T<-as.character(read.delim("all.gene.list",header=F)[[1]])
length(all_T)
intersect(all_T,all_F)->inter_gene
length(inter_gene)
setdiff(all_T,all_F)->diff_gene
diff_gene
diff_gene_new<-as.data.frame(diff_gene)
names(diff_gene_new)<-"gene_id"
merge(diff_gene_new,gf_12_new,by="gene_id")->diff_gene_new_gf_12
diff_gene_new_gf_12[1:4,]
write.table(diff_gene_new_gf_12,"diff_gene_new_gf_12",sep="\t",col.names=T,row.names=F,quote=F)
history()
history(200)
q()
dir()
fpkmg<-read.delim("genes.fpkm.matrix",header=T,sep="\t",check.names=F)
fpkmg[1:10,]
fpkmg<-read.delim("genes.fpkm.matrix",header=T,sep="\t",check.names=F,rowname=1)
?read.delim
fpkmg<-read.delim("genes.fpkm.matrix",header=T,sep="\t",check.names=F,row.name=1)
fpkmg[1:10,]
s1<-as.numeric(fpkmg[[1]])
class(s1)
min(s1)
max(s1)
s1[1:20]
sort(s1)[1:10]
sort(s1)[length(s1)-100:length(s1)]
sort(s1)[(length(s1)-100):length(s1)]
which(s1>=500)
length(which(s1>=500))
length(which(s1>=100))
i=1
  s1<-as.numeric(fpkmg[[i]])
    
s1
hist(s1)
hist(s1,breaks=50)
hist(s1,breaks=500)
 s1<-as.numeric(fpkmg[[i]])
    p1<-length(intersect(which(s1>=1),which(s1<=5)))
    p2<-length(intersect(which(s1>=6),which(s1<=10)))
    p3<-length(intersect(which(s1>=11),which(s1<=20)))
    p4<-length(intersect(which(s1>=21),which(s1<=50)))
    p5<-length(intersect(which(s1>=51),which(s1<=100)))
    p6<-length(which(s1>100))
    SUM<-nrow(fpkmg)
    lb1<-paste("[1,5](",p1,",",round(p1/SUM,2),"%",sep="")
    lb
    lb1
matrix(c(1,1,2,2,3,1,1,2,2,3,1,1,2,2,3),nrow=3,ncol=5,byrow=T)
matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T)
    s1<-as.numeric(fpkmg[[i]])
    p1<-length(intersect(which(s1>=1),which(s1<=5)))
    p2<-length(intersect(which(s1>=6),which(s1<=10)))
    p3<-length(intersect(which(s1>=11),which(s1<=20)))
    p4<-length(intersect(which(s1>=21),which(s1<=50)))
    p5<-length(intersect(which(s1>=51),which(s1<=100)))
    p6<-length(which(s1>100))
    SUM<-nrow(fpkmg)
    lb1<-paste("[1,5](",p1,",",round(p1/SUM,2),"%)",sep="")
    lb2<-paste("[6,10](",p1,",",round(p2/SUM,2),"%)",sep="")
    lb3<-paste("[11,20](",p1,",",round(p3/SUM,2),"%)",sep="")
    lb4<-paste("[21,50](",p1,",",round(p4/SUM,2),"%)",sep="")
    lb5<-paste("> 100(",p1,",",round(p5/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5)
    slices<-c(p1,p2,p3,p4,p5)
 pie(slices)
splices
slices
order(slices)
  ORDER<-order(slices,decreasing=T)
    slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    
order(slices_sort)
 pie(slices_sort)
names(fpkmg)
names(fpkmg)[i]
s1<-as.numeric(fpkmg[[i]])
    p1<-length(intersect(which(s1>=1),which(s1<=5)))
    p2<-length(intersect(which(s1>=6),which(s1<=10)))
    p3<-length(intersect(which(s1>=11),which(s1<=20)))
    p4<-length(intersect(which(s1>=21),which(s1<=50)))
    p5<-length(intersect(which(s1>=51),which(s1<=100)))
    p6<-length(which(s1>100))
    SUM<-nrow(fpkmg)
    lb1<-paste("[1,5](",p1,",",round(p1/SUM,2),"%)",sep="")
    lb2<-paste("[6,10](",p1,",",round(p2/SUM,2),"%)",sep="")
    lb3<-paste("[11,20](",p1,",",round(p3/SUM,2),"%)",sep="")
    lb4<-paste("[21,50](",p1,",",round(p4/SUM,2),"%)",sep="")
    lb5<-paste("> 100(",p1,",",round(p5/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5)
    slices<-c(p1,p2,p3,p4,p5)
    ORDER<-order(slices,decreasing=T)
    slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    col=rainbow(length(lbls))
 layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1/SUM,2),"%",sep=""),paste(round(p2/SUM,2),"%",sep=""),paste(round(p3/SUM,2),"%",sep=""),paste(round(p4/SUM,2),"%",sep=""))
    pie(slices_sort,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i] Gene Expression))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls_sort,fill=col,cex=1.5,bty="n")
s1<-as.numeric(fpkmg[[i]])
    p1<-length(intersect(which(s1>=1),which(s1<=5)))
    p2<-length(intersect(which(s1>=6),which(s1<=10)))
    p3<-length(intersect(which(s1>=11),which(s1<=20)))
    p4<-length(intersect(which(s1>=21),which(s1<=50)))
    p5<-length(intersect(which(s1>=51),which(s1<=100)))
    p6<-length(which(s1>100))
    SUM<-nrow(fpkmg)
    lb1<-paste("[1, 5](",p1,",",round(p1/SUM,2),"%)",sep="")
    lb2<-paste("[6, 10](",p1,",",round(p2/SUM,2),"%)",sep="")
    lb3<-paste("[11, 20](",p1,",",round(p3/SUM,2),"%)",sep="")
    lb4<-paste("[21, 50](",p1,",",round(p4/SUM,2),"%)",sep="")
    lb5<-paste(">  100(",p1,",",round(p5/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5)
    slices<-c(p1,p2,p3,p4,p5)
    ORDER<-order(slices,decreasing=T)
    slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    col=rainbow(length(lbls))
 layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1/SUM,2),"%",sep=""),paste(round(p2/SUM,2),"%",sep=""),paste(round(p3/SUM,2),"%",sep=""),paste(round(p4/SUM,2),"%",sep=""))
    pie(slices_sort,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls_sort,fill=col,cex=1.5,bty="n")
p1<-length(intersect(which(s1>=1),which(s1<=5)))
    p2<-length(intersect(which(s1>=6),which(s1<=10)))
    p3<-length(intersect(which(s1>=11),which(s1<=20)))
    p4<-length(intersect(which(s1>=21),which(s1<=50)))
    p5<-length(intersect(which(s1>=51),which(s1<=100)))
    p6<-length(which(s1>100))
p1
p2
p3
p4
p5
p6
 s1<-as.numeric(fpkmg[[i]])
    p1<-length(intersect(which(s1>=1),which(s1<=5)))
    p2<-length(intersect(which(s1>=6),which(s1<=10)))
    p3<-length(intersect(which(s1>=11),which(s1<=20)))
    p4<-length(intersect(which(s1>=21),which(s1<=50)))
    p5<-length(intersect(which(s1>=51),which(s1<=100)))
    p6<-length(which(s1>100))
    SUM<-nrow(fpkmg)
    lb1<-paste("[1, 5](",p1,",",round(p1/SUM,2),"%)",sep="")
    lb2<-paste("[6, 10](",p2,",",round(p2/SUM,2),"%)",sep="")
    lb3<-paste("[11, 20](",p3,",",round(p3/SUM,2),"%)",sep="")
    lb4<-paste("[21, 50](",p4,",",round(p4/SUM,2),"%)",sep="")
    lb5<-paste("[51, 100](",p5,",",round(p5/SUM,2),"%)",sep="")
    lb6<-paste(">  100(",p6,",",round(p6/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5,lb6)
    slices<-c(p1,p2,p3,p4,p5,p6)
    ORDER<-order(slices,decreasing=T)
    slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    col=rainbow(length(lbls))
 layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1/SUM,2),"%",sep=""),paste(round(p2/SUM,2),"%",sep=""),paste(round(p3/SUM,2),"%",sep=""),paste(round(p4/SUM,2),"%",sep=""),paste(round(p5/SUM,2),"%",sep=""),paste(round(p6/SUM,2),"%",sep=""))
    pie(slices_sort,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls_sort,fill=col,cex=1.5,bty="n")
ss<-as.numeric(fpkmg[[i]])
    s1<-ss[which(ss>0)];SUM<-length(s1)   
length(ss)
length(s1)
  s1<-as.numeric(fpkmg[[i]])
    #s1<-ss[which(ss>0)];
    SUM<-length(s1)   
    p1<-length(intersect(which(s1>=0),which(s1<=5)))
    p2<-length(intersect(which(s1>=6),which(s1<=10)))
    p3<-length(intersect(which(s1>=11),which(s1<=20)))
    p4<-length(intersect(which(s1>=21),which(s1<=50)))
    p5<-length(intersect(which(s1>=51),which(s1<=100)))
    p6<-length(which(s1>100))
     
    lb1<-paste("[0, 5](",p1,",",round(p1/SUM,2),"%)",sep="")
    lb2<-paste("[6, 10](",p2,",",round(p2/SUM,2),"%)",sep="")
    lb3<-paste("[11, 20](",p3,",",round(p3/SUM,2),"%)",sep="")
    lb4<-paste("[21, 50](",p4,",",round(p4/SUM,2),"%)",sep="")
    lb5<-paste("[51, 100](",p5,",",round(p5/SUM,2),"%)",sep="")
    lb6<-paste(">  100(",p6,",",round(p6/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5,lb6)
    slices<-c(p1,p2,p3,p4,p5,p6)
    ORDER<-order(slices,decreasing=T)
    slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    col=rainbow(length(lbls))
layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1/SUM,2),"%",sep=""),paste(round(p2/SUM,2),"%",sep=""),paste(round(p3/SUM,2),"%",sep=""),paste(round(p4/SUM,2),"%",sep=""),paste(round(p5/SUM,2),"%",sep=""),paste(round(p6/SUM,2),"%",sep=""))
    pie(slices_sort,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls_sort,fill=col,cex=1.5,bty="n")
  
    lb1<-paste("[0, 5](",p1,",",round(p1*100/SUM,2),"%)",sep="")
    lb2<-paste("[6, 10](",p2,",",round(p2*100/SUM,2),"%)",sep="")
    lb3<-paste("[11, 20](",p3,",",round(p3*100/SUM,2),"%)",sep="")
    lb4<-paste("[21, 50](",p4,",",round(p4*100/SUM,2),"%)",sep="")
    lb5<-paste("[51, 100](",p5,",",round(p5*100/SUM,2),"%)",sep="")
    lb6<-paste(">  100(",p6,",",round(p6*100/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5,lb6)
    slices<-c(p1,p2,p3,p4,p5,p6)
    ORDER<-order(slices,decreasing=T)
    slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    col=rainbow(length(lbls))
layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1/SUM,2),"%",sep=""),paste(round(p2/SUM,2),"%",sep=""),paste(round(p3/SUM,2),"%",sep=""),paste(round(p4/SUM,2),"%",sep=""),paste(round(p5/SUM,2),"%",sep=""),paste(round(p6/SUM,2),"%",sep=""))
    pie(slices_sort,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls_sort,fill=col,cex=1.5,bty="n")
 layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1/SUM,2),"%",sep=""),paste(round(p2/SUM,2),"%",sep=""),paste(round(p3/SUM,2),"%",sep=""),paste(round(p4/SUM,2),"%",sep=""),paste(round(p5/SUM,2),"%",sep=""),paste(round(p6/SUM,2),"%",sep=""))
    pie(slices,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls,fill=col,cex=1.5,bty="n")
  
    lb1<-paste("[0, 5](",p1,",",round(p1*100/SUM,2),"%)",sep="")
    lb2<-paste("[6, 10](",p2,",",round(p2*100/SUM,2),"%)",sep="")
    lb3<-paste("[11, 20](",p3,",",round(p3*100/SUM,2),"%)",sep="")
    lb4<-paste("[21, 50](",p4,",",round(p4*100/SUM,2),"%)",sep="")
    lb5<-paste("[51, 100](",p5,",",round(p5*100/SUM,2),"%)",sep="")
    lb6<-paste(">  100(",p6,",",round(p6*100/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5,lb6)
    slices<-c(p1,p2,p3,p4,p5,p6)
    #ORDER<-order(slices,decreasing=T)
    #slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    col=rainbow(length(lbls))
  layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1/SUM,2),"%",sep=""),paste(round(p2/SUM,2),"%",sep=""),paste(round(p3/SUM,2),"%",sep=""),paste(round(p4/SUM,2),"%",sep=""),paste(round(p5/SUM,2),"%",sep=""),paste(round(p6/SUM,2),"%",sep=""))
    pie(slices,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls,fill=col,cex=1.5,bty="n")
 s1<-as.numeric(fpkmg[[i]])
    #s1<-ss[which(ss>0)];
    SUM<-length(s1)   
    p1<-length(intersect(which(s1>=0),which(s1<=20)))
    p2<-length(intersect(which(s1>=21),which(s1<=40)))
    p3<-length(intersect(which(s1>=41),which(s1<=60)))
    p4<-length(intersect(which(s1>=61),which(s1<=80)))
    p5<-length(intersect(which(s1>=81),which(s1<=100)))
    p6<-length(which(s1>100))
     
    lb1<-paste("[0, 20](",p1,",",round(p1*100/SUM,2),"%)",sep="")
    lb2<-paste("[21, 40](",p2,",",round(p2*100/SUM,2),"%)",sep="")
    lb3<-paste("[41, 60](",p3,",",round(p3*100/SUM,2),"%)",sep="")
    lb4<-paste("[61, 80](",p4,",",round(p4*100/SUM,2),"%)",sep="")
    lb5<-paste("[81, 100](",p5,",",round(p5*100/SUM,2),"%)",sep="")
    lb6<-paste(">  100(",p6,",",round(p6*100/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5,lb6)
    slices<-c(p1,p2,p3,p4,p5,p6)
    #ORDER<-order(slices,decreasing=T)
    #slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    col=rainbow(length(lbls))
layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1/SUM,2),"%",sep=""),paste(round(p2/SUM,2),"%",sep=""),paste(round(p3/SUM,2),"%",sep=""),paste(round(p4/SUM,2),"%",sep=""),paste(round(p5/SUM,2),"%",sep=""),paste(round(p6/SUM,2),"%",sep=""))
    pie(slices,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls,fill=col,cex=1.5,bty="n")
layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1*100/SUM,2),"%",sep=""),paste(round(p2*100/SUM,2),"%",sep=""),paste(round(p3*100/SUM,2),"%",sep=""),paste(round(p4*100/SUM,2),"%",sep=""),paste(round(p5*100/SUM,2),"%",sep=""),paste(round(p6*100/SUM,2),"%",sep=""))
    pie(slices,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls,fill=col,cex=1.5,bty="n")
for(i in 1:ncol(fpkmg)){
    s1<-as.numeric(fpkmg[[i]])
    #s1<-ss[which(ss>0)];
    SUM<-length(s1)   
    p1<-length(intersect(which(s1>=0),which(s1<=20)))
    p2<-length(intersect(which(s1>=21),which(s1<=40)))
    p3<-length(intersect(which(s1>=41),which(s1<=60)))
    p4<-length(intersect(which(s1>=61),which(s1<=80)))
    p5<-length(intersect(which(s1>=81),which(s1<=100)))
    p6<-length(which(s1>100))
     
    lb1<-paste("[0, 20](",p1,",",round(p1*100/SUM,2),"%)",sep="")
    lb2<-paste("[21, 40](",p2,",",round(p2*100/SUM,2),"%)",sep="")
    lb3<-paste("[41, 60](",p3,",",round(p3*100/SUM,2),"%)",sep="")
    lb4<-paste("[61, 80](",p4,",",round(p4*100/SUM,2),"%)",sep="")
    lb5<-paste("[81, 100](",p5,",",round(p5*100/SUM,2),"%)",sep="")
    lb6<-paste(">  100(",p6,",",round(p6*100/SUM,2),"%)",sep="")    
    lbls<-c(lb1,lb2,lb3,lb4,lb5,lb6)
    slices<-c(p1,p2,p3,p4,p5,p6)
    #ORDER<-order(slices,decreasing=T)
    #slices_sort<-slices[ORDER];lbls_sort<-lbls[ORDER]
    col=rainbow(length(lbls))
    pdf(file=paste("GeneExpDistribution_",names(fpkmg)[i],".pdf",sep=""),height=5,width=10)   
    layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(0.1,0.1,0.1,0.1))
    pie_lbs<-c(paste(round(p1*100/SUM,2),"%",sep=""),paste(round(p2*100/SUM,2),"%",sep=""),paste(round(p3*100/SUM,2),"%",sep=""),paste(round(p4*100/SUM,2),"%",sep=""),paste(round(p5*100/SUM,2),"%",sep=""),paste(round(p6*100/SUM,2),"%",sep=""))
    pie(slices,labels=pie_lbs,col=col,main=paste("Distribution of ",names(fpkmg)[i],"Gene Expression",sep=""))
    par(mar=c(0.1,0.1,0.1,0.05))
    plot(1,1,col="white",axes=F)
    legend("left",legend=lbls,fill=col,cex=1.5,bty="n")
    dev.off()
 
}
q()
read.delim("map_ratios.new",header=T,sep="\t")
read.delim("map_ratios.new",header=T,sep="\t",check.names=F)
aa<-read.delim("map_ratios.new",header=T,sep="\t",check.names=F)
aa
q()
ls()
fpkmg[1:4,]
history(200)
history(1000)
s1<-as.numeric(fpkmg[[1]])
s1
ls()
q()
ls()
mapratio<-read.delim("map_ratios.new",header=T,sep="\t",check.names=F)
mapratio
q()
mapratio<-read.delim("map_ratios.new1",header=T,sep="\t",check.names=F)
mapratio
mapratio<-read.delim("map_ratios.new1",header=T,sep="\t",check.names=F,rownames=1)
mapratio<-read.delim("map_ratios.new1",header=T,sep="\t",check.names=F,row.names=1)
mapratio
q()
mapratio<-read.delim("map_ratios.new1",header=T,sep="\t",check.names=F,row.names=1)
mapratio
q()
mapratio<-read.delim("map_ratios.new1",header=T,sep="\t",check.names=F,row.names=1)
mapratio
q()
mapratio
mapratio
i
i=1
mapratio[i,]
mapratio[,i]
mapratio[,i+1]
as.character(mapratio[,i+1])
colnames(mapratio)
colnames(mapratio)[i]
gsub("_num","",colnames(mapratio)[i])
 lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    par(mar=c(4,4,4,0.1))
    pie(slices,labels=lbs,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
    par(mar=c(4,0.1,4,4))
 lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    par(mar=c(4,4,4,0.1))
    pie(slices,labels=lbls,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2)
)
matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T)
    lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(4,4,4,0.1))
    pie(slices,labels=lbls,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
    par(mar=c(4,0.1,4,4))
    plot(1,1,col="white",axes=F)
    legend("left",legend=rownames(mapratio),fill=col,cex=1.5,bty="n")
    lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(4,1,4,0.1))
    pie(slices,labels=lbls,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
    par(mar=c(4,0.1,4,1))
    plot(1,1,col="white",axes=F)
    legend("left",legend=rownames(mapratio),fill=col,cex=1,bty="n")
    
    
gsub("_num","",colnames(mapratio)[i])
lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    pdf(file=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),".pdf",sep=""),w=10,h=8)
    par(mar=c(4,1,4,0.1))
    pie(slices,labels=lbls,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
    par(mar=c(4,0.1,4,1))
    plot(1,1,col="white",axes=F)
    legend("left",legend=rownames(mapratio),fill=col,cex=1,bty="n")
    dev.off()
    
    dev.off()
lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    pdf(file=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),".pdf",sep=""),w=10,h=8)
    par(mar=c(4,1,4,0.1))
    pie(slices,labels=lbls,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
    par(mar=c(4,0.1,4,1))
    plot(1,1,col="white",axes=F)
    legend("left",legend=rownames(mapratio),fill=col,cex=1,bty="n")
    dev.off()
    
q()
for(i in 1:4){
    lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    pdf(file=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),".pdf",sep=""),w=10,h=8)
    layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(4,1,4,0.1))
    pie(slices,labels=lbls,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
    par(mar=c(4,0.1,4,1))
    plot(1,1,col="white",axes=F)
    legend("left",legend=rownames(mapratio),fill=col,cex=1,bty="n")
    dev.off()  
}
i
i
    lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    pdf(file=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),".pdf",sep=""),w=10,h=8)
    layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(4,1,4,0.1))
    pie(slices,labels=lbls,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
    par(mar=c(4,0.1,4,1))
    plot(1,1,col="white",axes=F)
    legend("left",legend=rownames(mapratio),fill=col,cex=1,bty="n")
    dev.off()  
}
paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep="")
colnames(mapratio)
for(i in c(1,3,5,7)){
    lbls<-as.character(mapratio[,i+1])
    slices<-mapratio[,i]
    col=rainbow(length(lbls)) 
    pdf(file=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),".pdf",sep=""),w=10,h=8)
    layout(matrix(c(1,1,1,2,2,1,1,1,2,2,1,1,1,2,2),nrow=3,ncol=5,byrow=T))
    par(mar=c(4,1,4,0.1))
    pie(slices,labels=lbls,col=col,main=paste("Mapping of ",gsub("_num","",colnames(mapratio)[i]),sep=""))
    par(mar=c(4,0.1,4,1))
    plot(1,1,col="white",axes=F)
    legend("left",legend=rownames(mapratio),fill=col,cex=1,bty="n")
    dev.off()  
}
q()
matrix0_path<-"IFN-0_vs_FIFN-0_all_GO.list.txt"
matrix1_path<-"IFN-0_vs_FIFN-0_up_GO.list.txt"
matrix2_path<-"IFN-0_vs_FIFN-0_down_GO.list.txt"
scan0<-scan(matrix0_path,what=character(),nlines=1)
scan1<-scan(matrix1_path,what=character(),nlines=1)
scan2<-scan(matrix2_path,what=character(),nlines=1)
matrix0_path
Num0<-as.numeric(scan0[length(scan1)])
Num1<-as.numeric(scan1[length(scan1)])
Num2<-as.numeric(scan2[length(scan2)])
Matrix0<-read.delim(matrix0_path,header=T,check.names=F,skip=1)
Matrix1<-read.delim(matrix1_path,header=T,check.names=F,skip=1)
Matrix2<-read.delim(matrix2_path,header=T,check.names=F,skip=1)
Matrix0[1:4,]
names(Matrix0)<-c("term_type_1","term_1","number_1","percent_1","GO")
names(Matrix1)<-c("term_type_1","term_1","number_1","percent_1","GO")
names(Matrix2)<-c("term_type_2","term_2","number_2","percent_2","GO")
MergeMatrix0<-unique(merge(Matrix0,Matrix1,by="GO",all=T))
MergeMatrix<-unique(merge(MergeMatrix0,Matrix2,by="GO",all=T))
MergeMatrix[1:4,]
names(MergeMatrix)
names(Matrix0)<-c("term_type","term","number_0","percent_0","GO")
names(Matrix1)<-c("term_type","term","number_1","percent_1","GO")
names(Matrix2)<-c("term_type","term","number_2","percent_2","GO")
MergeMatrix0<-unique(merge(Matrix0,Matrix1,by="GO",all=T))
MergeMatrix<-unique(merge(MergeMatrix0,Matrix2,by="GO",all=T))
MergeMatrix[1:4,]
names(MergeMatrix)
MergeMatrix[1:4,]
MergeMatrix
     
MergeMatrix
names(Matrix0)<-c("term_type_0","term_0","number_0","percent_0","GO")
names(Matrix1)<-c("term_type_1","term_1","number_1","percent_1","GO")
names(Matrix2)<-c("term_type_2","term_2","number_2","percent_2","GO")
MergeMatrix0<-unique(merge(Matrix0,Matrix1,by="GO",all=T))
MergeMatrix<-unique(merge(MergeMatrix0,Matrix2,by="GO",all=T))
MergeMatrix[1:4,]
MergeMatrix[1:40,]
MergeMatrix[1:20,]
NA0_lines<-which(is.na(MergeMatrix[["term_0"]])==TRUE)
NA1_lines<-which(is.na(MergeMatrix[["term_1"]])==TRUE)
NA2_lines<-which(is.na(MergeMatrix[["term_2"]])==TRUE)
NA0_lines
NA1_lines
NA2_lines
names(MergeMatrix)
Num0
Num0<-as.numeric(scan0[length(scan0)])
Num0
q()
ls()
matrix0_path<-"IFN-0_vs_FIFN-0_all_GO.list.txt"
matrix1_path<-"IFN-0_vs_FIFN-0_up_GO.list.txt"
matrix2_path<-"IFN-0_vs_FIFN-0_down_GO.list.txt"
scan0<-scan(matrix0_path,what=character(),nlines=1)
scan1<-scan(matrix1_path,what=character(),nlines=1)
scan2<-scan(matrix2_path,what=character(),nlines=1)
Num0<-as.numeric(scan0[length(scan0)])
Num1<-as.numeric(scan1[length(scan1)])
Num2<-as.numeric(scan2[length(scan2)])
Matrix0<-read.delim(matrix0_path,header=T,check.names=F,skip=1)
Matrix1<-read.delim(matrix1_path,header=T,check.names=F,skip=1)
Matrix2<-read.delim(matrix2_path,header=T,check.names=F,skip=1)
### 
All_GO<-unique(rbind(Matrix0[,c("GO","term","term_type")],Matrix1[,c("GO","term","term_type")],Matrix2[,c("GO","term","term_type")]))
for(i in 1:nrow(All_GO)){
    GOID<-as.character(All_GO[i,"GO"])
    number_all<-Matrix0[which(Matrix0[["GO"]]==GOID),"number"]
    number_up<-Matrix0[which(Matrix1[["GO"]]==GOID),"number"]
    number_down<-Matrix0[which(Matrix2[["GO"]]==GOID),"number"]
    if(length(number_all)==0){number_all=0};if(length(number_up)==0){number_up=0};if(length(number_down)==0){number_down=0}
    All_GO[i,4]<-paste(number_up,"(",round(100*(number_up/Num0),1),"%)",sep="")
    All_GO[i,5]<-paste(number_down,"(",round(100*(number_down/Num0),1),"%)",sep="")
    All_GO[i,6]<-paste(number_all,"(",round(100*(number_all/Num0),1),"%)",sep="")
}
names(All_GO)[4:6]<-c("Number(Proportion)_up","Number(Proportion)_up","Number(Proportion)_up")
All_GO_new<-as.matrix(All_GO)
All_GO_new
All_GO_new<-as.matrix(All_GO)
All_GO_new[i+1,1]<-"Total"
All_GO_new[i+1,2]<-"Total"
All_GO_new[i+1,3]<-"Total"
All_GO_new[i+1,4]<-paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep="")
All_GO_new[i+1,5]<-paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep="")
All_GO_new[i+1,6]<-paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep="")
All_GO_new_new<-as.data.frame(All_GO_new)
i
All_GO_new[1:10,]
All_GO_new<-as.matrix(All_GO)
All_GO_new[i+1,1]<-"Total"
All_GO_new[i+1,2]<-"Total"
All_GO_new[i+1,3]<-"Total"
names(All_GO)[4:6]<-c("Number(Proportion)_up","Number(Proportion)_up","Number(Proportion)_up")
All_GO_new<-as.matrix(All_GO)
dim(All_GO_new)
All_GO_new<-as.data.frame(as.matrix(All_GO))
All_GO_new[i+1,1]<-"Total"
All_GO_new[i+1,2]<-"Total"
All_GO_new[i+1,3]<-"Total"
All_GO_new[i+1,4]<-paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep="")
All_GO_new[i+1,5]<-paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep="")
All_GO_new[i+1,6]<-paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep="")
All_GO_new_new<-as.data.frame(All_GO_new)
All_GO_new
q()
ls()
All_GO
i
All_GO_new[i+1,1]<-"Total"
All_GO_new<-as.matrix(All_GO)
All_GO_new[i+1,1]<-"Total"
dim(All_GO_new)
i + 1
All_GO_new<-as.data.frame(as.matrix(All_GO))
All_GO_new[i+1,1]<-"Total"
All_GO_new[i+1,2]<-"Total"
All_GO_new[i+1,3]<-"Total"
ll
names(All_GO)<-NULL
tailline<-as.data.frame(rep("Total",3),paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep=""))
All_GO_new<-rbind(All_GO,tailline)
names(All_GO_new)<-c("GO","term","term_type","Number(Proportion)_up","Number(Proportion)_up","Number(Proportion)_up")
tailline<-as.data.frame(c(rep("Total",3),paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep="")))
tailline
All_GO_new<-rbind(All_GO,tailline)
names(All_GO_new)<-c("GO","term","term_type","Number(Proportion)_up","Number(Proportion)_up","Number(Proportion)_up")
?as.data.frame
?data.frame
names(All_GO)<-NULL
tailline<-data.frame(c(rep("Total",3),paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep="")),6,1,byrow=T)
All_GO_new<-rbind(All_GO,tailline)
All_GO
tailline
tailline<-data.frame(c(rep("Total",3),paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep="")),1,6,byrow=T)
All_GO_new<-rbind(All_GO,tailline)
tailline
tailline<-data.frame(c(rep("Total",3),paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep="")),1,6)
tailline
dim(tailline)
c(rep("Total",3),paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep=""))
?data.frame
tailline<-c(rep("Total",3),paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep=""))
All_GO[i+1,]<-tailline
All_GO
names(All_GO_new)<-c("GO","term","term_type","Number(Proportion)_up","Number(Proportion)_up","Number(Proportion)_up")
All_GO
names(All_GO)<-c("GO","term","term_type","Number(Proportion)_up","Number(Proportion)_up","Number(Proportion)_up")
All_GO
q()
group_name<-scan("group_name",what=character())
group_name
p=1
paste(getwd(),"/",group_name[p],"/",group_name[p],"_all_GO.list.txt",sep="")
group_name<-scan("group_name",what=character())
for(p in 1:length(group_name)){
matrix0_path<-paste(getwd(),"/",group_name[p],"/",group_name[p],"_all_GO.list.txt",sep="")
matrix1_path<-paste(getwd(),"/",group_name[p],"/",group_name[p],"_up_GO.list.txt",sep="")
matrix2_path<-paste(getwd(),"/",group_name[p],"/",group_name[p],"_down_GO.list.txt",sep="")
scan0<-scan(matrix0_path,what=character(),nlines=1)
scan1<-scan(matrix1_path,what=character(),nlines=1)
scan2<-scan(matrix2_path,what=character(),nlines=1)
Num0<-as.numeric(scan0[length(scan0)])
Num1<-as.numeric(scan1[length(scan1)])
Num2<-as.numeric(scan2[length(scan2)])
Matrix0<-read.delim(matrix0_path,header=T,check.names=F,skip=1)
Matrix1<-read.delim(matrix1_path,header=T,check.names=F,skip=1)
Matrix2<-read.delim(matrix2_path,header=T,check.names=F,skip=1)
### 
All_GO<-unique(rbind(Matrix0[,c("GO","term","term_type")],Matrix1[,c("GO","term","term_type")],Matrix2[,c("GO","term","term_type")]))
for(i in 1:nrow(All_GO)){
    GOID<-as.character(All_GO[i,"GO"])
    number_all<-Matrix0[which(Matrix0[["GO"]]==GOID),"number"]
    number_up<-Matrix0[which(Matrix1[["GO"]]==GOID),"number"]
    number_down<-Matrix0[which(Matrix2[["GO"]]==GOID),"number"]
    if(length(number_all)==0){number_all=0};if(length(number_up)==0){number_up=0};if(length(number_down)==0){number_down=0}
    All_GO[i,4]<-paste(number_up,"(",round(100*(number_up/Num0),1),"%)",sep="")
    All_GO[i,5]<-paste(number_down,"(",round(100*(number_down/Num0),1),"%)",sep="")
    All_GO[i,6]<-paste(number_all,"(",round(100*(number_all/Num0),1),"%)",sep="")
}
names(All_GO)<-NULL
tailline<-c(rep("Total",3),paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep=""))
All_GO[i+1,]<-tailline
names(All_GO)<-c("GO","term","term_type","Number(Proportion)_up","Number(Proportion)_up","Number(Proportion)_up")
write.table(All_GO,paste(getwd(),"/",group_name[p],"/",group_name[p],"_GO_num_distribution.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
print(group_name[p])
}
q()
ls()
dir()
KEGG_path2ENST<-read.delim("KEGG_path2ENST",sep="\t",header=T,check.names=F)
KEGG_path2ENST[1:10,]
KEGG_path2ENST[1:,]
KEGG_path2ENST[1:2,]
names(KEGG_path2ENST)
q()
transcript2gene<-read.delim("transcript2gene",header=F,sep="\t")
ls()
KEGG_path2ENST[1:4,]
ls()
transcript2gene[1:4,]
names(KEGG_path2ENST)
names(transcript2gene)
ls()
for(i in 1:nrow(KEGG_path2ENST)){
}
path2gene<-data.frame()
for(i in 1:nrow(KEGG_path2ENST)){
i=1
unique(unlist(strsplit(as.character(KEGG_path2ENST[i,2]),";",fix=T))
)
sapply(unique(unlist(strsplit(as.character(KEGG_path2ENST[i,2]),";",fix=T))),function(x) unlist(strsplit(x,"(",fix=T))[1])
transcript_list<-sapply(unique(unlist(strsplit(as.character(KEGG_path2ENST[i,2]),";",fix=T))),function(x) unlist(strsplit(x,"(",fix=T))[1])
    
transcript_list
 names(transcript_list)<-NULL
    
    
transcript_list
    transcript_list<-as.data.frame(sapply(unique(unlist(strsplit(as.character(KEGG_path2ENST[i,2]),";",fix=T))),function(x) unlist(strsplit(x,"(",fix=T))[1]))
    names(transcript_list)<-NULL
names(transcript_list)
transcript_list<-as.data.frame(sapply(unique(unlist(strsplit(as.character(KEGG_path2ENST[i,2]),";",fix=T))),function(x) unlist(strsplit(x,"(",fix=T))[1]))
transcript_list
row.names(transcript_list)<-NULL
    col.names(transcript_list)<-"transcript"
    
    
    colnames(transcript_list)<-"transcript"
 rownames(transcript_list)<-NULL
    colnames(transcript_list)<-"transcript"
    gene_list<-unique(as.character(merge(transcript_list,transcript2gene,by="transcript",all.x=T)[["gene"]]))
    
    names(transcript2gene)<-c("transcript","gene")
  gene_list<-unique(as.character(merge(transcript_list,transcript2gene,by="transcript",all.x=T)[["gene"]]))
    
    
    gene_list
as.character(KEGG_path2ENST[i,1])
paste(gene_list,sep=";")
paste(gene_list,collapse=";")
path2gene<-data.frame()
flag<-1
for(i in 1:nrow(KEGG_path2ENST)){
    transcript_list<-as.data.frame(sapply(unique(unlist(strsplit(as.character(KEGG_path2ENST[i,2]),";",fix=T))),function(x) unlist(strsplit(x,"(",fix=T))[1]))
    rownames(transcript_list)<-NULL
    colnames(transcript_list)<-"transcript"
    gene_list<-unique(as.character(merge(transcript_list,transcript2gene,by="transcript",all.x=T)[["gene"]]))
    N<-length(gene_list)
    path2gene[i,1]<-as.character(KEGG_path2ENST[i,1])
    path2gene[i,2]<-paste(gene_list,collapse=";")
}
path2gene[1:4,]
path2gene[1:4,1]
ls()
group_name
group_name<-scan("group_name",what=character())
group_name
p=1
matrix0_path<-paste(getwd(),"/",group_name[p],"/","all.list",sep="")
matrix1_path<-paste(getwd(),"/",group_name[p],"/","up.list",sep="")
matrix2_path<-paste(getwd(),"/",group_name[p],"/","down.list",sep="")
scan0<-scan(matrix0_path,what=character())
scan1<-scan(matrix1_path,what=character())
scan2<-scan(matrix2_path,what=character())
allSDEGinpathway<-unique(unlist(strsplit(paste(as.character(path2gene[[2]]),collapse=";"),";",fix=T)))
allSDEGinpathway
length(allSDEGinpathway)
Num0<-length(intersect(scan0,allGeneinpathway))
Num1<-length(intersect(scan1,allGeneinpathway))
Num2<-length(intersect(scan2,allGeneinpathway))
allGeneinpathway<-unique(unlist(strsplit(paste(as.character(path2gene[[2]]),collapse=";"),";",fix=T)))
Num0<-length(intersect(scan0,allGeneinpathway))
Num1<-length(intersect(scan1,allGeneinpathway))
Num2<-length(intersect(scan2,allGeneinpathway))
Num0
Num1
Num2
length(scan0)
p
matrix0_path<-paste(getwd(),"/",group_name[p],"/","all.list",sep="")
matrix1_path<-paste(getwd(),"/",group_name[p],"/","up.list",sep="")
matrix2_path<-paste(getwd(),"/",group_name[p],"/","down.list",sep="")
scan0<-scan(matrix0_path,what=character())
scan1<-scan(matrix1_path,what=character())
scan2<-scan(matrix2_path,what=character())
all_SDEGinpathway<-intersect(scan0,allGeneinpathway)
up_SDEGinpathway<-intersect(scan1,allGeneinpathway)
down_SDEGinpathway<-intersect(scan2,allGeneinpathway)
Num0<-length(all_SDEGinpathway)
Num1<-length(up_SDEGinpathway)
Num2<-length(down_SDEGinpathway)
path_up_down_all<-data.frame()
for(i in 1:nrow(path2gene)){
    
    path_up_down_all[i,1]<-as.character(path2gene[i,1])
    gene_list<-unlist(strsplit(as.character(path2gene[i,2]),";",fix=T))    
    number_all<-length(intersect(gene_list,all_SDEGinpathway))
    number_up<-length(intersect(gene_list,up_SDEGinpathway))
    number_down<-length(intersect(gene_list,down_SDEGinpathway))
    if(length(number_all)==0){number_all=0};if(length(number_up)==0){number_up=0};if(length(number_down)==0){number_down=0}
    path_up_down_all[i,2]<-paste(number_up,"(",round(100*(number_up/Num0),1),"%)",sep="")
    path_up_down_all[i,3]<-paste(number_down,"(",round(100*(number_down/Num0),1),"%)",sep="")
    path_up_down_all[i,4]<-paste(number_all,"(",round(100*(number_all/Num0),1),"%)",sep="")
}
path_up_down_all
path_up_down_all<-data.frame()
flag<-1
for(i in 1:nrow(path2gene)){
    
    path_up_down_all[flag,1]<-as.character(path2gene[i,1])
    gene_list<-unlist(strsplit(as.character(path2gene[i,2]),";",fix=T))    
    number_all<-length(intersect(gene_list,all_SDEGinpathway))
    number_up<-length(intersect(gene_list,up_SDEGinpathway))
    number_down<-length(intersect(gene_list,down_SDEGinpathway))
    if(length(number_all)==0){number_all=0};if(length(number_up)==0){number_up=0};if(length(number_down)==0){number_down=0}
    if(number_all==0 && number_up==0 && number_down==0){next}
    path_up_down_all[flag,2]<-paste(number_up,"(",round(100*(number_up/Num0),1),"%)",sep="")
    path_up_down_all[flag,3]<-paste(number_down,"(",round(100*(number_down/Num0),1),"%)",sep="")
    path_up_down_all[flag,4]<-paste(number_all,"(",round(100*(number_all/Num0),1),"%)",sep="")
    }
path_up_down_all
p
path_up_down_all<-data.frame()
flag<-1
for(i in 1:nrow(path2gene)){
    
    path_up_down_all[flag,1]<-as.character(path2gene[i,1])
    gene_list<-unlist(strsplit(as.character(path2gene[i,2]),";",fix=T))    
    number_all<-length(intersect(gene_list,all_SDEGinpathway))
    number_up<-length(intersect(gene_list,up_SDEGinpathway))
    number_down<-length(intersect(gene_list,down_SDEGinpathway))
    if(length(number_all)==0){number_all=0};if(length(number_up)==0){number_up=0};if(length(number_down)==0){number_down=0}
    if(number_all==0 && number_up==0 && number_down==0){next}
    path_up_down_all[flag,2]<-paste(number_up,"(",round(100*(number_up/Num0),1),"%)",sep="")
    path_up_down_all[flag,3]<-paste(number_down,"(",round(100*(number_down/Num0),1),"%)",sep="")
    path_up_down_all[flag,4]<-paste(number_all,"(",round(100*(number_all/Num0),1),"%)",sep="")
    flag<-flag+1
    }
path_up_down_all
dim(path_up_down_all)
flag
names(path_up_down_all)<-NULL
tailline<-c("Total",paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep=""))
path_up_down_all[flag,]<-tailline
path_up_down_all
names(path_up_down_all)<-NULL
tailline<-c("Total",paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep=""))
path_up_down_all[flag,]<-tailline
names(path_up_down_all)<-c("KEGG Pathway","Number(Proportion)_up","Number(Proportion)_down","Number(Proportion)_all")
write.table(path_up_down_all,paste(getwd(),"/",group_name[p],"/",group_name[p],"_KEGGPathway_num_distribution.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
print(group_name[p])
### all genes with pathway
allGeneinpathway<-unique(unlist(strsplit(paste(as.character(path2gene[[2]]),collapse=";"),";",fix=T)))
for(p in 1:length(group_name)){
matrix0_path<-paste(getwd(),"/",group_name[p],"/","all.list",sep="")
matrix1_path<-paste(getwd(),"/",group_name[p],"/","up.list",sep="")
matrix2_path<-paste(getwd(),"/",group_name[p],"/","down.list",sep="")
scan0<-scan(matrix0_path,what=character())
scan1<-scan(matrix1_path,what=character())
scan2<-scan(matrix2_path,what=character())
all_SDEGinpathway<-intersect(scan0,allGeneinpathway)
up_SDEGinpathway<-intersect(scan1,allGeneinpathway)
down_SDEGinpathway<-intersect(scan2,allGeneinpathway)
Num0<-length(all_SDEGinpathway)
Num1<-length(up_SDEGinpathway)
Num2<-length(down_SDEGinpathway)
path_up_down_all<-data.frame()
flag<-1
for(i in 1:nrow(path2gene)){    
    path_up_down_all[flag,1]<-as.character(path2gene[i,1])
    gene_list<-unlist(strsplit(as.character(path2gene[i,2]),";",fix=T))    
    number_all<-length(intersect(gene_list,all_SDEGinpathway))
    number_up<-length(intersect(gene_list,up_SDEGinpathway))
    number_down<-length(intersect(gene_list,down_SDEGinpathway))
    if(length(number_all)==0){number_all=0};if(length(number_up)==0){number_up=0};if(length(number_down)==0){number_down=0}
    if(number_all==0 && number_up==0 && number_down==0){next}
    path_up_down_all[flag,2]<-paste(number_up,"(",round(100*(number_up/Num0),1),"%)",sep="")
    path_up_down_all[flag,3]<-paste(number_down,"(",round(100*(number_down/Num0),1),"%)",sep="")
    path_up_down_all[flag,4]<-paste(number_all,"(",round(100*(number_all/Num0),1),"%)",sep="")
    flag<-flag+1
}
names(path_up_down_all)<-NULL
tailline<-c("Total",paste(Num1,"(",round(100*(Num1/Num0),1),"%)",sep=""),paste(Num2,"(",round(100*(Num2/Num0),1),"%)",sep=""),paste(Num0,"(",round(100*(Num0/Num0),1),"%)",sep=""))
path_up_down_all[flag,]<-tailline
names(path_up_down_all)<-c("KEGG Pathway","Number(Proportion)_up","Number(Proportion)_down","Number(Proportion)_all")
write.table(path_up_down_all,paste(getwd(),"/",group_name[p],"/",group_name[p],"_KEGGPathway_num_distribution.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
print(group_name[p])
}
q()
a<-c("a","b")
 write.table(a,"test",sep="\t")
q()
dir()
tmp<-read.delim("genes.diff_exp_ann_result",header=T,sep="\t",check.names=F)
tmp[1:4,]
names(tmp)
names(tmp)
tmp[1:4,]
tmp[1:4,]
 
tmp[1:4,]
?which
data1<-tmp[,c("gene_id","IFN-0_FPKM","IFN-5_FPKM","significant","Direction"),]
data1[1:10,]
tmp[1:40,"Direction"]
data1[1:10,]
tmp[1:60,"Direction"]
data1[1:60,]
data1[1:10,]
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]]))
log10IFN0[1:10]
log10(0)
log10IFN0[1:20]
data2<-data1[-(which(data1[[5]]=="_")),]
data2[1:20,]
which(data1[[5]]=="_")
class(data1[10,5])
data2<-data1[-(which(as.character(data1[[5]]=="-"))),]
data2<-data1[-(which(as.character(data1[[5]])=="-")),]
data2[1:10,]
data2[1:20,]
q()
ls()
history()
history()
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]]))
log10IFN5<-log10(as.numeric(data1[["IFN-5_FPKM"]]))
data2[1:10,]
which(data2[["Direction"]]=="-")
aa<-which(data2[["Direction"]]=="_")
aa
data2[1:20,]
data2[1:40,]
lines0<-unique(c(which(as.numeric(data2[[2]])==0),which(as.numeric(data2[[3]])==0)))
lines0
data2[1lines0][1:100,]
data2[1lines0,][1:100,]
data2[1lines0[1:20],]
data2[lines0[1:20],]
data2[lines0[1:40],]
log10(0)
?smoothScatter
?smoothScatter
x1  <- matrix(rnorm(n), ncol = 2)
     x2  <- matrix(rnorm(n, mean = 3, sd = 1.5), ncol = 2)
     x   <- rbind(x1, x2)
     
     oldpar <- par(mfrow = c(2, 2))
     smoothScatter(x, nrpoints = 0)
     smoothScatter(x)
?smoothScatter
 n <- 10000
     x1  <- matrix(rnorm(n), ncol = 2)
     x2  <- matrix(rnorm(n, mean = 3, sd = 1.5), ncol = 2)
     x   <- rbind(x1, x2)
     
     oldpar <- par(mfrow = c(2, 2))
     smoothScatter(x, nrpoints = 0)
     smoothScatter(x)
     smoothScatter(x, nrpoints = 0)
?smoothScatter
?colorRampPalette
illed.contour(volcano,
                    color.palette =
                        colorRampPalette(c("red", "white", "blue")),
                    asp = 1)
     filled.contour(volcano,
                    color.palette =
                        colorRampPalette(c("red", "white", "blue"),
                                         space = "Lab"),
                    asp = 1)
     
     ## Interpolating a 'sequential' ColorBrewer palette
     YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
     filled.contour(volcano,
                    color.palette = colorRampPalette(YlOrBr, space = "Lab"),
                    asp = 1)
     filled.contour(volcano,
                    color.palette = colorRampPalette(YlOrBr, space = "Lab",
                                                     bias = 0.5),
                    asp = 1)
     
     ## 'jet.colors' is "as in Matlab"
     ## (and hurting the eyes by over-saturation)
     jet.colors <-
       colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
     filled.contour(volcano, color = jet.colors, asp = 1)
     
     ## space="Lab" helps when colors don't form a natural sequence
     m <- outer(1:20,1:20,function(x,y) sin(sqrt(x*y)/3))
     rgb.palette <- colorRampPalette(c("red", "orange", "blue"),
                                     space = "rgb")
     Lab.palette <- colorRampPalette(c("red", "orange", "blue"),
                                     space = "Lab")
     filled.contour(m, col = rgb.palette(20))
     filled.contour(m, col = Lab.palette(20))
     filled.contour(volcano, color = jet.colors, asp = 1)
ls()
history(200)
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]])+1)
log10IFN5<-log10(as.numeric(data1[["IFN-5_FPKM"]])+1)
log10IFN0[1:100]
history(200)
names(data2)
data2[1:10,]
data2[1:50,]
uplines<-intersect(which(data2[["significant"]]=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(data2[["significant"]]=="yes"),which(data2[["Direction"]]=="down"))
length(uplines)
length(downlines)
downlines
uplines
data2[1:20,]
data2<-data1[-(which(as.character(data1[[5]])=="-")),]
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]])+1)
log10IFN5<-log10(as.numeric(data1[["IFN-5_FPKM"]])+1)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
uplines
which(as.character(data2[["significant"]])=="yes")
uplines<-intersect(which(data2[,"significant"])=="yes"),which(data2[,"Direction"]=="up"))
uplines<-intersect(which(data2[,"significant"]=="yes"),which(data2[,"Direction"]=="up"))
uplines
data2[15,]
data2<-data1[-(which(as.character(data1[[5]])=="-")),]
rownames(data2)<-NULL
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]])+1)
log10IFN5<-log10(as.numeric(data1[["IFN-5_FPKM"]])+1)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
data2[1:20,]
uplines
downlines
data2[downlines,]
plot(log10IFN5,log10IFN0)
data1[1:220,]
data1[1:20,]
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]])+0.001)
log10IFN5<-log10(as.numeric(data1[["IFN-5_FPKM"]])+0.001)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
plot(log10IFN5,log10IFN0)
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]])+0.01)
log10IFN5<-log10(as.numeric(data1[["IFN-5_FPKM"]])+0.01)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
plot(log10IFN5,log10IFN0)
plot(log10IFN5,log10IFN0,xlim=c(-2,max(c(log10IFN0,log10IFN5))),ylim=c(-2,max(c(log10IFN0,log10IFN5))))
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]])+1)
log10IFN5<-log10(as.numeric(data1[["IFN-5_FPKM"]])+1)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
plot(log10IFN5,log10IFN0,xlim=c(-2,max(c(log10IFN0,log10IFN5))),ylim=c(-2,max(c(log10IFN0,log10IFN5))))
log10IFN0<-log10(as.numeric(data1[["IFN-0_FPKM"]])+0.01)
log10IFN5<-log10(as.numeric(data1[["IFN-5_FPKM"]])+0.01)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
plot(log10IFN5,log10IFN0,xlim=c(-2,max(c(log10IFN0,log10IFN5))),ylim=c(-2,max(c(log10IFN0,log10IFN5))))
getwd()
data2[1:50,]
data2<-data1[-(which(as.character(data1[[5]])=="-")),]
rownames(data2)<-NULL
log10IFN0<-log10(as.numeric(data2[["IFN-0_FPKM"]])+0.01)
log10IFN5<-log10(as.numeric(data2[["IFN-5_FPKM"]])+0.01)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
plot(log10IFN5,log10IFN0)
data2<-data1[-(which(as.character(data1[[5]])=="-")),]
rownames(data2)<-NULL
log10IFN0<-log10(as.numeric(data2[["IFN-0_FPKM"]])+0.01)
log10IFN5<-log10(as.numeric(data2[["IFN-5_FPKM"]])+0.01)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
plot(log10IFN5,log10IFN0)
write.table(data2,"data2",sep="\t",col.names=T,row.names=F,quote=F)
q()
ls()
history()
history()
log10IFN0<-log10(as.numeric(data2[["IFN-0_FPKM"]])+0.01)
log10IFN5<-log10(as.numeric(data2[["IFN-5_FPKM"]])+0.01)
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
plot(log10IFN5,log10IFN0)
    
uplines
pchlist<-rep(1,length(log10IFN5))
pchlist[c(uplines,downlines)]<-rep(16,length(c(uplines,downlines)))
collist<-rep("black",length(log10IFN5))
collist[uplines]<-rep("red",length(uplines))
collist[downlines]<-rep("blue",length(uplines))
plot(log10IFN5,log10IFN0,pch=pchlist,col=collist)
collist<-rep("black",length(log10IFN5))
collist[uplines]<-rep("#EE4F2B",length(uplines))
collist[downlines]<-rep("#2C60AD",length(uplines))
plot(log10IFN5,log10IFN0,pch=pchlist,col=collist)
  
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist)
   
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist,alpha=0.5)
collist[uplines]<-rep(rgb(238,79,43,alpha=0.2),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,alpha=0.2),length(uplines))
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist)
?rgb
collist[uplines]<-rep(rgb(238,79,43,max=255,alpha=0.2),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,max=255,alpha=0.2),length(uplines))
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist)
collist
rgb(238,79,43,max=255,alpha=0.2)
?rgb
?rgb
collist[uplines]<-rep(rgb(238,79,43,max=255,alpha=250),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,max=255,alpha=250),length(uplines))
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist)
collist<-rep(rgb(34,30,31,max=255,alpha=250),length(log10IFN5))
collist[uplines]<-rep(rgb(238,79,43,max=255,alpha=250),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,max=255,alpha=250),length(uplines))
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist)
alpha<-200
collist<-rep(rgb(34,30,31,max=255,alpha=alpha),length(log10IFN5))
collist[uplines]<-rep(rgb(238,79,43,max=255,alpha=alpha),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,max=255,alpha=alpha),length(uplines))
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist)
?legend
nochange<-length(log10IFN5)-length(uplines)-length(downlines)
pchlist<-rep(1,length(log10IFN5))
pchlist[c(uplines,downlines)]<-rep(16,length(c(uplines,downlines)))
alpha<-200
collist<-rep(rgb(34,30,31,max=255,alpha=alpha),length(log10IFN5))
collist[uplines]<-rep(rgb(238,79,43,max=255,alpha=alpha),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,max=255,alpha=alpha),length(uplines))
legend<-c(paste("up-regulated genes (",length(uplines),")",sep=""),paste("down-regulated genes (",length(downlines),")",sep=""),paste("not differential expressed (",length(nochange),")")
)
legend
length(log10IFN5)
length(log10IFN5)-length(uplines)-length(downlines)
legend<-c(paste("up-regulated genes (",length(uplines),")",sep=""),paste("down-regulated genes (",length(downlines),")",sep=""),paste("not differential expressed (",nochange,")"))
legend
nochange<-length(log10IFN5)-length(uplines)-length(downlines)
pchlist<-rep(1,length(log10IFN5))
pchlist[c(uplines,downlines)]<-rep(16,length(c(uplines,downlines)))
alpha<-200
sigs<-"FDR<=0.05"
collist<-rep(rgb(34,30,31,max=255,alpha=alpha),length(log10IFN5))
collist[uplines]<-rep(rgb(238,79,43,max=255,alpha=alpha),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,max=255,alpha=alpha),length(uplines))
legend<-c(sigs,paste("up-regulated genes (",length(uplines),")",sep=""),paste("down-regulated genes (",length(downlines),")",sep=""),paste("not differential expressed (",nochange,")",sep=""))
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist,xlab=,ylab=)
legend("topleft",pch=15,legend=legend)
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist,xlab=,ylab=)
legend("topleft",pch=15,legend=legend,inset=0.1,col=c("white",rgb(238,79,43,max=255,alpha=alpha),rgb(44,96,173,max=255,alpha=alpha),rgb(34,30,31,max=255,alpha=alpha)))
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist,xlab=,ylab=)
legend("topleft",pch=15,legend=legend,inset=0.015,col=c("white",rgb(238,79,43,max=255,alpha=alpha),rgb(44,96,173,max=255,alpha=alpha),rgb(34,30,31,max=255,alpha=alpha)),cex=3)
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist,xlab=,ylab=)
legend("topleft",pch=15,legend=legend,inset=0.015,col=c("white",rgb(238,79,43,max=255,alpha=alpha),rgb(44,96,173,max=255,alpha=alpha),rgb(34,30,31,max=255,alpha=alpha)),cex=1)
legend("topleft",pch=15,legend=legend,inset=0.015,col=c("white",rgb(238,79,43,max=255,alpha=alpha),rgb(44,96,173,max=255,alpha=alpha),rgb(34,30,31,max=255,alpha=alpha)),cex=1.2)
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist,xlab="IFN-0 expression (log10)",ylab="IFN-5 expression (log10)")
legend("topleft",pch=16,legend=legend,inset=0.015,col=c("white",rgb(238,79,43,max=255,alpha=alpha),rgb(44,96,173,max=255,alpha=alpha),rgb(34,30,31,max=255,alpha=alpha)),cex=1)
plot(log10IFN0,log10IFN5,pch=pchlist,col=collist,xlab="IFN-0 expression (log10)",ylab="IFN-5 expression (log10)")
legend("topleft",pch=16,legend=legend,inset=0.015,col=c("white",rgb(238,79,43,max=255),rgb(44,96,173,max=255),rgb(34,30,31,max=255)),cex=1)
tmp[1:4,]
names(data1)[2:3]
gsub("_FPKM","",names(data1)[2:3])
tmp<-read.delim("genes.diff_exp_ann_result",header=T,sep="\t",check.names=F)
#data1<-tmp[,c("gene_id","IFN-0_FPKM","IFN-5_FPKM","significant","Direction"),]
data1<-tmp[,c(1,7,8,12,13)]
data2<-data1[-(which(as.character(data1[[5]])=="-")),]
rownames(data2)<-NULL
sample1_log10<-log10(as.numeric(data2[[7]])+0.01)
sample2_log10<-log10(as.numeric(data2[[8]])+0.01)
samplenames<-gsub("_FPKM","",names(data1)[2:3])
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
nochange<-length(sample1_log10)-length(uplines)-length(downlines)
pchlist<-rep(1,length(sample1_log10))
pchlist[c(uplines,downlines)]<-rep(16,length(c(uplines,downlines)))
alpha<-200
sigs<-"FDR<=0.05"
collist<-rep(rgb(34,30,31,max=255,alpha=alpha),length(sample1_log10))
collist[uplines]<-rep(rgb(238,79,43,max=255,alpha=alpha),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,max=255,alpha=alpha),length(uplines))
legend<-c(sigs,paste("up-regulated genes (",length(uplines),")",sep=""),paste("down-regulated genes (",length(downlines),")",sep=""),paste("not differential expressed (",nochange,")",sep=""))
plot(sample1_log10,sample2_log10,pch=pchlist,col=collist,xlab=paste(samplenames[1]," expression (log10)",sep=""),ylab=paste(samplenames[1]," expression (log10)",sep=""))
legend("topleft",pch=16,legend=legend,inset=0.015,col=c("white",rgb(238,79,43,max=255),rgb(44,96,173,max=255),rgb(34,30,31,max=255)),cex=1)
rownames(data2)<-NULL
sample1_log10<-log10(as.numeric(data2[[2]])+0.01)
sample2_log10<-log10(as.numeric(data2[[3]])+0.01)
sample1_log10
samplenames<-gsub("_FPKM","",names(data1)[2:3])
uplines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="up"))
downlines<-intersect(which(as.character(data2[["significant"]])=="yes"),which(data2[["Direction"]]=="down"))
nochange<-length(sample1_log10)-length(uplines)-length(downlines)
pchlist<-rep(1,length(sample1_log10))
pchlist[c(uplines,downlines)]<-rep(16,length(c(uplines,downlines)))
alpha<-200
sigs<-"FDR<=0.05"
collist<-rep(rgb(34,30,31,max=255,alpha=alpha),length(sample1_log10))
collist[uplines]<-rep(rgb(238,79,43,max=255,alpha=alpha),length(uplines))
collist[downlines]<-rep(rgb(44,96,173,max=255,alpha=alpha),length(uplines))
legend<-c(sigs,paste("up-regulated genes (",length(uplines),")",sep=""),paste("down-regulated genes (",length(downlines),")",sep=""),paste("not differential expressed (",nochange,")",sep=""))
plot(sample1_log10,sample2_log10,pch=pchlist,col=collist,xlab=paste(samplenames[1]," expression (log10)",sep=""),ylab=paste(samplenames[1]," expression (log10)",sep=""))
legend("topleft",pch=16,legend=legend,inset=0.015,col=c("white",rgb(238,79,43,max=255),rgb(44,96,173,max=255),rgb(34,30,31,max=255)),cex=1)
plot(sample1_log10,sample2_log10,pch=pchlist,col=collist,xlab=paste(samplenames[1]," expression (log10)",sep=""),ylab=paste(samplenames[2]," expression (log10)",sep=""))
legend("topleft",pch=16,legend=legend,inset=0.015,col=c("white",rgb(238,79,43,max=255),rgb(44,96,173,max=255),rgb(34,30,31,max=255)),cex=1)
q()
annotation_matrix<-read.csv("mart_export.txt",header=T)
annotation_matrix[1:4,]
annotation_matrix[1:40,]
write.table(annotation_matrix,"annotation.matrix",sep="\t",col.names=T,row.names=F,quote=F)
q()
ls
ls()
sortSigDEmiRs<-function(groupName,path_in,path_out,flag){
   setwd(path_in)
   ##groupName<-list.files()
  ## all_sig_miRs<-vector()
   for(i in 1:length(groupName)){
       print(groupName[i])
       DE_Matrix<-read.delim(paste(groupName[i],"//output_score.txt",sep=""),sep="\t",header=TRUE)
       sig_miRs<-DE_Matrix[which(as.character(DE_Matrix[["Signature.p.value...0.001."]])=="TRUE"),"GeneNames"]
       sig_miRs<-unique(as.character(sig_miRs))
       write.table(sig_miRs,paste(path_out,"//",groupName[i],"sigMir_",flag,sep=""),sep="\n",col.names=F,row.names=F,quote=F)
       ## all_sig_miRs<-c(all_sig_miRs,sig_miRs)
       ## all_sig_miRs<-unique(all_sig_miRs)
    }
       ## print(length(all_sig_miRs))
       ## write.table(all_sig_miRs,paste(path_out,"//","all_sig_miRs","_",flag,sep=""),sep="\n",col.names=F,row.names=F,quote=F)
}
path1<-"/data/users/wangyan/project/fangmeixia_SR/all/known_miRNA"
path2<-"/data/users/wangyan/project/fangmeixia_SR/all/novo_miRNA"
path_out<-"/data/users/wangyan/project/fangmeixia_SR/all/sigDEmiRs"
groupName<-as.character(read.delim("groupName",header=F)[[1]])
groupName
sortSigDEmiRs(groupName[4],path1,path_out,"known")
sortSigDEmiRs(groupName[4],path2,path_out,"novo")
q()
groupName<-list.files()
groupName
q()
groupName
all_sig_miRs<-vector()
for(i in 1:length(groupName)){
   print(groupName[i])
   DE_Matrix<-read.delim(paste(groupName[i],"//output_score.txt",sep=""),sep="\t",header=TRUE)
   sig_miRs<-DE_Matrix[which(as.character(DE_Matrix[["Signature(p-value < 0.001)"]])=="TRUE"),"GeneNames"]
   sig_miRs<-unique(as.character(sig_miRs))
   all_sig_miRs<-c(all_sig_miRs,sig_miRs)
}
all_sig_miRs<-unique(all_sig_miRs)
length(all_sig_miRs)
write.table(all_sig_miRs,"all_sig_miRs",sep="\n",col.names=F,row.names=F,quote=F)
i
i=1
 DE_Matrix<-read.delim(paste(groupName[i],"//output_score.txt",sep=""),sep="\t",header=TRUE)
   sig_miRs<-DE_Matrix[which(as.character(DE_Matrix[["Signature(p-value < 0.001)"]])=="TRUE"),"GeneNames"]
sig_miRs
DE_Matrix[1:10,]
dim(DE_Matrix)
 DE_Matrix<-read.delim(paste(groupName[i],"//output_score.txt",sep=""),sep="\t",header=TRUE)
dim(DE_Matrix)
DE_Matrix[1000:1283,]
names(DE_Matrix)
 sig_miRs<-DE_Matrix[which(as.character(DE_Matrix[["Signature.p.value...0.001."]])=="TRUE"),"GeneNames"]
   sig_miRs<-unique(as.character(sig_miRs))
   sig_miRs
all_sig_miRs<-vector()
for(i in 1:length(groupName)){
   print(groupName[i])
   DE_Matrix<-read.delim(paste(groupName[i],"//output_score.txt",sep=""),sep="\t",header=TRUE)
   sig_miRs<-DE_Matrix[which(as.character(DE_Matrix[["Signature.p.value...0.001."]])=="TRUE"),"GeneNames"]
   sig_miRs<-unique(as.character(sig_miRs))
   all_sig_miRs<-c(all_sig_miRs,sig_miRs)
}
all_sig_miRs<-unique(all_sig_miRs)
length(all_sig_miRs)
write.table(all_sig_miRs,"all_sig_miRs",sep="\n",col.names=F,row.names=F,quote=F)
length(all_sig_miRs)
all_sig_miRs[1:10]
q()
abline(a = 50, b = 180, col = "red")
barplot(seq(1:5), horiz=TRUE)
barplot(seq(1:5), horiz=TRUE, space=0)
barplot(seq(1:5), horiz=TRUE, space=0, col = "yellow")
library(doBy)
PuroA.mean = summaryBy(rate ~ conc, data = PuroA,
FUN = mean)
plot(rate ~ conc, data = PuroA, pch = 16, col = 4,
cex = 1.5)
points(mean.rate ~ conc, data = PuroA.mean, col = "cyan",
lwd = 10, pch = "x")
plot(rate ~ conc, data = PuroA)
smooth1 <- with(PuroA, lowess(rate ~ conc, f = 0.9))
smooth2 <- with(PuroA, lowess(rate ~ conc, f = 0.3))
lines(smooth1, col = "red")
lines(smooth2, col = "blue")
m1 <- lm(rate ~ conc, data = PuroA)
m2 <- lm(rate ~ conc + I(conc^2), data = PuroA)
m3 <- lm(rate ~ conc + I(conc^2) + I(conc^3),
data = PuroA)
lines(fitted(m1) ~ conc, data = PuroA, col = "red")
lines(fitted(m2) ~ conc, data = PuroA, col = "blue")
lines(fitted(m3) ~ conc, data = PuroA, col = "cyan")
plot(rate ~ conc, data = PuroA)
abline(lm(rate ~ conc, data = PuroA))
abline(a = 100, b = 105, col = "blue")
abline(h = 200, col = "red")
abline(v = 0.6, col = "green")
mysymb <- c(1, 2)[Puromycin$state]
plot(rate ~ conc, data = Puromycin, col = mysymb,
pch = mysymb)
PuroB <- subset(Puromycin, state == "untreated")
rm(list=ls(0))
rm(list=ls())
x <- pmin(3, pmax(-3, stats::rnorm(50)))
x
pmin(1:5)
x
pmin(1:5)
y <- pmin(3, pmax(-3, stats::rnorm(50)))
xhist <- hist(x, breaks=seq(-3,3,0.5), plot=FALSE)
fix(xhist)
yhist <- hist(y, breaks=seq(-3,3,0.5), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(-3,3); yrange <- c(-3,3)
layout(matrix(c(2,0,1,3), 2, 2, byrow=TRUE),
c(3,1), c(1,3), TRUE)
layout.show(3)
layout.show(n = 2)
layout.show(n = 1)
layout.show(n = 4)
layout.show(n = 3)
layout.show(n = 4)
layout.show(n = 1)
layout.show(n = 3)
plot(x, y, xlim=xrange, ylim=yrange, xlab="", ylab="")
par(mar=c(0,3,1,1))
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
par(mar=c(3,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top),
space=0, horiz=TRUE)
x <- seq(from=1, to= 10, by= 0.01)
y <- sin(x)
rm(list=ls())
pos = seq(from=1, to= 1000, by = 0.1)
pos
clear
pos
cov_plus = sin(x)
cov_plus = sin(pos)
cov_plus = cos(pos)
cov_plus = sin(pos)
cov_minus = cos(pos)
plot(pos, cov_plus)
plot(pos, cov_plus, lty= "l")
plot(pos, cov_plus, lty= "1")
plot(pos, cov_plus)
pos = seq(from=1, to= 100, by = 0.1)
cov_plus = sin(pos)
plot(pos, cov_plus)
plot(pos, cov_plus, lty = 1)
plot(pos, cov_plus, lty = 2)
?lty
?par
plot(pos, cov_plus, lty = 1)
plot(pos, cov_plus, lty = c(1))
plot(pos, cov_plus, lty = "solid")
plot(pos, cov_plus, lty = 1, pch = 2)
plot(pos, cov_plus, lty = 1, pch = 1)
line(pos, cov_plus)
plot(0)
line(pos, cov_plus)
line(pos, cov_plus, lty = 1)
line(pos, cov_plus)
plot(0)
line(pos, cov_plus)
lines(pos, cov_plus)
plot(0,type="n")
lines(pos, cov_plus, type= 's')
plot(0,type="n", xlim = c(min(pos), max(pos)))
lines(pos, cov_plus, type= 's')
layout(mat=matrix(1:4, 4, 1), heights= c(3,1,1,1))
layout.show
layout.show()
layout.show()
layout.show(2)
layout.show(3)
layout.show(4)
layout.show(1)
lines(pos, cov_plus, type= 's')
plot(0,type="n", xlim=c(min(pos), max(pos)))
lines(pos, cov_plus, type= 's')
plot(0,type="n")
plot(0,type="n", xlim=c(min(pos), max(pos)))
layout.show()
layout.show()
rm(list=ls())
layout(matrix(c(2,0,1,3), 2, 2, byrow=TRUE),
c(3,1), c(1,3), TRUE)
layout.show(3)
layout.show(1)
layout.show(2)
layout.show(3)
layout.show(0)
layout.show()
layout.show()
layout.show(3)
x <- 1:10
y <- sin(x)
plot(x,y)
layout.show(3)
plot(x,y)
plot(x,2y)
plot(x,2*y)
plot(x, 1)
rep(1,10)
plot(x, rep(1,10))
rm(list=ls())
layout(mat=matrix(c(2,1,4,3), 2, 2, byrow=TRUE))
layout.show(n=)
layout.show()
layout.show(2)
layout.show(4)
layout.show(3)
layout(mat=matrix(c(1,2,3,4), 4, 1))
layout.show(1)
layout.show(2)
layout.show(3)
layout.show(4)
layout(mat=matrix(c(1,2,3,4), 4, 1), heights=c(1,2,3,4))
layout.show(1)
layout.show(2)
layout.show(4)
layout(mat=matrix(c(1,2,3,4), 4, 1), c(4,2,2,2))
layout.show(1)
layout.show(2)
layout.show(3)
layout.show(4)
layout(mat=matrix(c(1,2,3,4), 4, 1), heights=c(1,1,1,1))
layout.show(1)
layout(mat=matrix(c(1,2,3,4), 4, 1), heights=c(2,1,1,1))
layout.show(1)
layout.show(2)
layout.show(3)
layout.show(4)
layout(mat=matrix(c(1,2,3,4), 4, 1), heights=c(4,1,1,1))
layout.show(1)
layout.show(4)
layout(mat=matrix(c(1,2,3,4), 4, 1), heights=c(1,2,3,4))
layout.show(4)
layout.show(1)
layout.show(2)
layout.show(3)
layout.show(4)
x = 1:10:0.01
x
x = 1:10
x
y = 1* x
plot(x,y)
layout.show(4)
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
layout.show(4)
par(mar = rep(2, 4))
plot(x,y)
x
y
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
par(mar = rep(10, 4))
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
par(mar = rep(20, 4))
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
par(mar = rep(0, 4))
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
plot(x,2*y, ylim = c(min(2*y), max(2*y)), xlim = c(min(x), max(x)))
plot(x,3*y, ylim = c(min(3*y), max(3*y)), xlim = c(min(x), max(x)))
plot(x,4*y, ylim = c(min(4*y), max(4*y)), xlim = c(min(x), max(x)))
layout.show(4)
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
par(mar = rep(0, 4))
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
title("test")
title("test")
title("test2")
plot(x,2*y, ylim = c(min(2*y), max(2*y)), xlim = c(min(x), max(x)))
layout(mat=matrix(c(1,2,3,4), 4, 1), heights=c(1,2,3,4),respect=TRUE)
par(mar = rep(0, 4))
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
plot(x,2*y, ylim = c(min(2*y), max(2*y)), xlim = c(min(x), max(x)))
layout(mat=matrix(c(1,2,3,4), 4, 1), heights=c(1,2,3,4),respect=TRUE)
layout.show()
layout.show(2)
layout.show(3)
layout.show(4)
par(mar = rep(0, 4))
plot(x,y, ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
plot(x,2*y, ylim = c(min(2*y), max(2*y)), xlim = c(min(x), max(x)))
title("picture2")
par(mar = rep(2, 4))
plot(x,3*y, ylim = c(min(3*y), max(3*y)), xlim = c(min(x), max(x)))
par(mar = rep(1, 4))
plot(x,3*y, ylim = c(min(3*y), max(3*y)), xlim = c(min(x), max(x)))
title("picture3")
par(mar = rep(2, 4))
plot(x,4*y, ylim = c(min(4*y), max(4*y)), xlim = c(min(x), max(x)))
par(mar = rep(1, 4))
plot(x,4*y, ylim = c(min(4*y), max(4*y)), xlim = c(min(x), max(x)))
title("picture4")
layout.show()
layout.show(1)
layout.show(2)
layout.show(3)
layout.show(4)
par()
?par
?pdf
par()
par(no.readonly=TRUE)
par(no.readonly=TRUE)
op = par(no.readonly=TRUE)
fix(op)
x
x = seq(1:10)
y = sin(x)
plot(x,y)
title("sin(x)")
par(mar=c(5,5,5,5))
plot(x,y)
title("sin(x)")
par(mar=c(10,10,10,10))
plot(x,y)
title("sin(x)")
par(mar=rep(100,4))
plot(x,y)
title("sin(x)")
par(mar=rep(0,4))
plot(x,y)
title("sin(x)")
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)
?par
x <- c(1:10)
y <- x
z <- 10/x
opar <- par(no.readonly=TRUE)
par(mar=c(5, 4, 4, 8) + 0.1)
plot(x, y, type="b",
pch=21, col="red",
yaxt="n", lty=3, ann=FALSE)
lines(x, z, type="b", pch=22, col="blue", lty=2)
axis(2, at=x, labels=x, col.axis="red", las=2)
axis(4, at=z, labels=round(z, digits=2),
col.axis="blue", las=2, cex.axis=0.7, tck=-.01)
mtext("y=1/x", side=4, line=3, cex.lab=1, las=2, col="blue")
title("An Example of Creative Axes",
xlab="X values",
ylab="Y=X")
par(opar)
?legend
library(Hmisc)
install.packages("Hmisc")
library(Hmisc)
opar <- par(no.readonly=TRUE)
par(lwd=2, cex=1.5, font.lab=2)
plot(dose, drugA, type="b",
pch=15, lty=1, col="red", ylim=c(0, 60),
main="Drug A vs. Drug B",
xlab="Drug Dosage", ylab="Drug Response")
lines(dose, drugB, type="b",
pch=17, lty=2, col="blue")
abline(h=c(30), lwd=1.5, lty=2, col="gray")
library(Hmisc)
minor.tick(nx=3, ny=3, tick.ratio=0.5)
legend("topleft", inset=.05, title="Drug Type", c("A","B"),
lty=c(1, 2), pch=c(15, 17), col=c("red", "blue"))
par(opar)
attach(mtcars)
plot(wt, mpg,
main="Mileage vs. Car Weight",
xlab="Weight", ylab="Mileage",
pch=18, col="blue")
text(wt, mpg,
row.names(mtcars),
cex=0.6, pos=4, col="red")
detach(mtcars)
row.names(mtcars)
mtcars
col.names(mtcars)
names(mtcars)
opar <- par(no.readonly=TRUE)
par(cex=1.5)
plot(1:7,1:7,type="n")
text(3,3,"Example of default text")
text(4,4,family="mono","Example of mono-spaced text")
text(5,5,family="serif","Example of serif text")
par(opar)
attach(mtcars)
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,2))
plot(wt,mpg, main="Scatterplot of wt vs. mpg")
plot(wt,disp, main="Scatterplot of wt vs disp")
hist(wt, main="Histogram of wt")
boxplot(wt, main="Boxplot of wt")
par(opar)
detach(mtcars)
attach(mtcars)
opar <- par(no.readonly=TRUE)
par(mfrow=c(3,1))
hist(wt)
hist(mpg)
hist(disp)
par(opar)
detach(mtcars)
demo(plotmath)
plotmathTranslate("hello")
plotmathTranslate("x^y")
plotmathTranslate(x^y)
help(plotmath)
plot(10,10,type="n")
text(5,5,labels=expression(x^y+2))
text(10,10,labels=expression(x^y+2))
plot(10,10,type="n")
text(10,10,labels=expression(x%+-%y))
attach(mtcars)
opar <- par(no.readonly=TRUE)
par(mfrow=c(3,1))
hist(wt)
hist(mpg)
hist(disp)
par(opar)
detach(mtcars)
attach(mtcars)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
hist(wt)
hist(mpg)
hist(disp)
detach(mtcars)
attach(mtcars)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
widths=c(3, 1), heights=c(1, 2))
hist(wt)
hist(mpg)
hist(disp)
detach(mtcars)
attach(mtcars)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
widths=c(3, 1), heights=c(1, 2))
hist(wt)
hist(mpg)
hist(disp)
detach(mtcars)
attach(mtcars)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
widths=c(3, 1), heights=c(1, 2))
hist(wt)
hist(mpg)
hist(disp)
detach(mtcars)
opar = par(no.readonly=TRUE)
par(mar=c(0,0,0,0))
attach(mtcars)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
widths=c(3, 1), heights=c(1, 2))
hist(wt)
hist(mpg)
hist(disp)
detach(mtcars)
par(mar=c(1,1,1,1))
attach(mtcars)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
widths=c(3, 1), heights=c(1, 2))
hist(wt)
hist(mpg)
hist(disp)
detach(mtcars)
par(mar=c(2,2,2,2))
attach(mtcars)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
widths=c(3, 1), heights=c(1, 2))
hist(wt)
hist(mpg)
hist(disp)
detach(mtcars)
par(mar=c(4,4,4,4))
attach(mtcars)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
widths=c(3, 1), heights=c(1, 2))
hist(wt)
hist(mpg)
hist(disp)
detach(mtcars)
library(ggplot2)
install.packages("ggplot2")
library(ggplot2)
transmission <- factor(mtcars$am, levels=c(0, 1),
labels=c("Automatic", "Manual"))
qplot(wt,mpg, data=mtcars,
color=transmission, shape=transmission,
geom=c("point", "smooth"),
method="lm", formula=y~x,
xlab="Weight", ylab="Miles Per Gallon",
main="Regression Example")
library(ggplot2)
mtcars$cylinder <- as.factor(mtcars$cyl)
qplot(cylinder, mpg, data=mtcars, geom=c("boxplot", "jitter"),
fill=cylinder,
main="Box plots with superimposed data points",
xlab= "Number of Cylinders",
ylab="Miles per Gallon")
library(playwith)
library(lattice)
playwith(
xyplot(mpg~wt|factor(cyl)*factor(am),
data=mtcars, subscripts=TRUE,
type=c("r", "p"))
)
install.packages("playwith")
library(playwith)
library(lattice)
playwith(
xyplot(mpg~wt|factor(cyl)*factor(am),
data=mtcars, subscripts=TRUE,
type=c("r", "p"))
)
library(iplots)
attach(mtcars)
cylinders <- factor(cyl)
gears <- factor(gear)
transmission <- factor(am)
ihist(mpg)
ibar(gears)
iplot(mpg, wt)
ibox(mtcars[c("mpg", "wt", "qsec", "disp", "hp")])
ipcp(mtcars[c("mpg", "wt", "qsec", "disp", "hp")])
imosaic(transmission, cylinders)
detach(mtcars)
install.packages("iplots")
library(iplots)
attach(mtcars)
cylinders <- factor(cyl)
gears <- factor(gear)
transmission <- factor(am)
ihist(mpg)
ibar(gears)
iplot(mpg, wt)
ibox(mtcars[c("mpg", "wt", "qsec", "disp", "hp")])
ipcp(mtcars[c("mpg", "wt", "qsec", "disp", "hp")])
imosaic(transmission, cylinders)
detach(mtcars)
?plot
?pie
pie
setwd("C:/Users/gang/Desktop/yuanwencheng")
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
dev.off()
dev.off()
dev.off()
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
pie
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
diff(level2Values)
dx
dx = diff(level2Values)
dx
nx <- length(dx)
nx
col
col = c("white", "lightblue", "mistyrose", "lightcyan",
"lavender", "cornsilk")
col <- rep(col, length.out = nx)
col
border <- rep(border, length.out = nx)
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
source('C:/Users/gang/Desktop/yuanwencheng/go_pie.R')
q()
pdf(file=pdfname,width=20,height=8)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(2,2501,5001,7500,9970)],as.character(c(0,0.25,0.5,0.75,1)),cex=1,adj=0,xpd=T)
text(rep(0,length(continue_cols)+100),"Pvalue",adj=0,xpd=T)
dev.off()
length(continue_cols)+100
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=20,height=8)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(2,2501,5001,7500,9970)],as.character(c(0,0.25,0.5,0.75,1)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+100,"Pvalue",adj=0,xpd=T)
dev.off()
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=20,height=8)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(2,2501,5001,7500,9970)],as.character(c(0,0.25,0.5,0.75,1)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+1000,"Pvalue",adj=0,xpd=T)
dev.off()
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=20,height=8)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(2,2501,5001,7500,9900)],as.character(c(0,0.25,0.5,0.75,1)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=20,height=8)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(2,2501,5001,7500,9900)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=20,height=8)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(2,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
typeII<-c("Pathway Class:\n","EIP: Environmental Information Processing","GIP: Genetic Information Processing","CP : Cellular Processes","OS : Organismal Systems","DD : Drug Development","HD : Human Diseases","M  : Metabolism")
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=20,height=8)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
typeII<-c("Pathway Class:","EIP: Environmental Information Processing","GIP: Genetic Information Processing","CP : Cellular Processes","OS : Organismal Systems","DD : Drug Development","HD : Human Diseases","M  : Metabolism")
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=20,height=8)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.9,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.89,1,0,0.9),mar = c(10,0.1,6,2)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.95,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.94,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.90,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.3,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.90,1,0,0.9),mar = c(10,0.1,6,6)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.4,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(10,6,6,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c(10,0.1,6,4)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(10,6,4,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c(10,0.1,4,4)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
q()
history(2000)
history(20000)
history(200000)
history(2000000)
history(20000000)
history(200000000)
history(100000)
enrichBox<-read.delim("1_vs_3_enrich_short_type2",sep="\t",header=T,check.names=F)
names(enrichBox)
ramp_cc<-colorRamp(as.vector(c("#008B45","white")))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10000),5)),max=255),"E5",sep="")
### read enrich table
enrichBox_sort<-enrichBox[order(enrichBox[,7],enrichBox[,6],-(enrichBox[[2]]/enrichBox[[3]]),-enrichBox[,4],decreasing=F),]
RichFactor<-round(enrichBox_sort[[2]]/enrichBox_sort[[3]],2)
enrichBox[1:4,]
 enrichBox<-read.delim("1_vs_3_enrich_short_type2",sep="\t",header=T,check.names=F)
ramp_cc<-colorRamp(as.vector(c("#008B45","white")))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10000),5)),max=255),"E5",sep="")
### read enrich table
enrichBox_sort<-enrichBox[order(enrichBox[,7],enrichBox[,6],-(enrichBox[[2]]/enrichBox[[3]]),-enrichBox[,4],decreasing=F),]
RichFactor<-round(enrichBox_sort[[2]]/enrichBox_sort[[3]],2)
enrichBox[1:10,]
 enrichBox<-read.delim("1_vs_3_enrich_short_type2",sep="\t",header=T,check.names=F)
ramp_cc<-colorRamp(as.vector(c("#008B45","white")))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10000),5)),max=255),"E5",sep="")
### read enrich table
enrichBox_sort<-enrichBox[order(enrichBox[,7],enrichBox[,6],-(enrichBox[[2]]/enrichBox[[3]]),-enrichBox[,4],decreasing=F),]
RichFactor<-round(enrichBox_sort[[2]]/enrichBox_sort[[3]],2)
### significant stars:
Pvalue<-as.numeric(enrichBox_sort[[4]])
n_3tars<-which(Pvalue<=0.001)
n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
barcols<-getMyColor(targets=Pvalue)
### pathway class:
name2type<-enrichBox_sort[,c(1,7)]
shortName<-gsub("Organismal Systems","OS",gsub("Human Diseases","HD",gsub("Genetic Information Processing","GIP",gsub("Environmental Information Processing","EIP",gsub("Drug Development","DD",gsub("Cellular Processes","CP",gsub("Metabolism","M",name2type[[2]])))))))
xlabs<-paste(enrichBox_sort[[1]],shortName,sep=": ")
typeII<-c("Pathway Class:","EIP: Environmental Information Processing","GIP: Genetic Information Processing","CP : Cellular Processes","OS : Organismal Systems","DD : Drug Development","HD : Human Diseases","M  : Metabolism")
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(10,6,4,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c(10,0.1,4,4)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
enrichBox<-read.delim("2_vs_3_enrich_short_type2",sep="\t",header=T,check.names=F)
##### define 100,000 continue colors
##### springgreen4: #008B45; steelblue: #4682B4 (#4782B9)
getMyColor<-function(col1="#008B45",col2="white",cutn=10000,targets){
targets=round(targets,4)
ramp_cc<-colorRamp(as.vector(c(col1,col2)))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=cutn),4)),max=255),"E5",sep="")
num2col<-as.data.frame(cbind(round(seq(0,1,length=cutn),4),continue_cols))
names(num2col)<-c("number","color")
targetcolor<-vector()
for(i in 1:length(targets)){
targetcolor[i]<-as.character(num2col[which(num2col[[1]]==targets[i]),2])
}
return(targetcolor)
}
### color card:
ramp_cc<-colorRamp(as.vector(c("#008B45","white")))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10000),5)),max=255),"E5",sep="")
### read enrich table
enrichBox_sort<-enrichBox[order(enrichBox[,7],enrichBox[,6],-(enrichBox[[2]]/enrichBox[[3]]),-enrichBox[,4],decreasing=F),]
RichFactor<-round(enrichBox_sort[[2]]/enrichBox_sort[[3]],2)
### significant stars:
Pvalue<-as.numeric(enrichBox_sort[[4]])
n_3tars<-which(Pvalue<=0.001)
n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
barcols<-getMyColor(targets=Pvalue)
### pathway class:
name2type<-enrichBox_sort[,c(1,7)]
shortName<-gsub("Organismal Systems","OS",gsub("Human Diseases","HD",gsub("Genetic Information Processing","GIP",gsub("Environmental Information Processing","EIP",gsub("Drug Development","DD",gsub("Cellular Processes","CP",gsub("Metabolism","M",name2type[[2]])))))))
xlabs<-paste(enrichBox_sort[[1]],shortName,sep=": ")
typeII<-c("Pathway Class:","EIP: Environmental Information Processing","GIP: Genetic Information Processing","CP : Cellular Processes","OS : Organismal Systems","DD : Drug Development","HD : Human Diseases","M  : Metabolism")
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(10,6,4,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c(10,0.1,4,4)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
q()
history(20000)
enrichBox<-read.delim("1_vs_2_enrich_short_type2",sep="\t",header=T,check.names=F)
##### define 100,000 continue colors
##### springgreen4: #008B45; steelblue: #4682B4 (#4782B9)
getMyColor<-function(col1="#008B45",col2="white",cutn=10000,targets){
targets=round(targets,4)
ramp_cc<-colorRamp(as.vector(c(col1,col2)))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=cutn),4)),max=255),"E5",sep="")
num2col<-as.data.frame(cbind(round(seq(0,1,length=cutn),4),continue_cols))
names(num2col)<-c("number","color")
targetcolor<-vector()
for(i in 1:length(targets)){
targetcolor[i]<-as.character(num2col[which(num2col[[1]]==targets[i]),2])
}
return(targetcolor)
}
enrichBox[1:4,]
enrichBox<-read.delim("1_vs_2_enrich_short_type2",sep="\t",header=T,check.names=F)
##### define 100,000 continue colors
##### springgreen4: #008B45; steelblue: #4682B4 (#4782B9)
getMyColor<-function(col1="#008B45",col2="white",cutn=10000,targets){
targets=round(targets,4)
ramp_cc<-colorRamp(as.vector(c(col1,col2)))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=cutn),4)),max=255),"E5",sep="")
num2col<-as.data.frame(cbind(round(seq(0,1,length=cutn),4),continue_cols))
names(num2col)<-c("number","color")
targetcolor<-vector()
for(i in 1:length(targets)){
targetcolor[i]<-as.character(num2col[which(num2col[[1]]==targets[i]),2])
}
return(targetcolor)
}
### color card:
ramp_cc<-colorRamp(as.vector(c("#008B45","white")))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10000),5)),max=255),"E5",sep="")
### read enrich table
enrichBox_sort<-enrichBox[order(enrichBox[,7],enrichBox[,6],-(enrichBox[[2]]/enrichBox[[3]]),-enrichBox[,4],decreasing=F),]
RichFactor<-round(enrichBox_sort[[2]]/enrichBox_sort[[3]],2)
### significant stars:
Pvalue<-as.numeric(enrichBox_sort[[4]])
n_3tars<-which(Pvalue<=0.001)
n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
barcols<-getMyColor(targets=Pvalue)
### pathway class:
name2type<-enrichBox_sort[,c(1,7)]
shortName<-gsub("Organismal Systems","OS",gsub("Human Diseases","HD",gsub("Genetic Information Processing","GIP",gsub("Environmental Information Processing","EIP",gsub("Drug Development","DD",gsub("Cellular Processes","CP",gsub("Metabolism","M",name2type[[2]])))))))
xlabs<-paste(enrichBox_sort[[1]],shortName,sep=": ")
typeII<-c("Pathway Class:","EIP: Environmental Information Processing","GIP: Genetic Information Processing","CP : Cellular Processes","OS : Organismal Systems","DD : Drug Development","HD : Human Diseases","M  : Metabolism")
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(10,6,4,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c(10,0.1,4,4)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
RichFactor
history(2000)
enrichBox[1:4,]
names(enrichBox)
enrichBox<-read.delim("1_vs_2_enrich_short_type2",sep="\t",header=T,check.names=F)[,-1]
##### define 100,000 continue colors
##### springgreen4: #008B45; steelblue: #4682B4 (#4782B9)
getMyColor<-function(col1="#008B45",col2="white",cutn=10000,targets){
targets=round(targets,4)
ramp_cc<-colorRamp(as.vector(c(col1,col2)))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=cutn),4)),max=255),"E5",sep="")
num2col<-as.data.frame(cbind(round(seq(0,1,length=cutn),4),continue_cols))
names(num2col)<-c("number","color")
targetcolor<-vector()
for(i in 1:length(targets)){
targetcolor[i]<-as.character(num2col[which(num2col[[1]]==targets[i]),2])
}
return(targetcolor)
}
### color card:
ramp_cc<-colorRamp(as.vector(c("#008B45","white")))
continue_cols<-paste(rgb(ramp_cc(round(seq(0,1,length=10000),5)),max=255),"E5",sep="")
### read enrich table
enrichBox_sort<-enrichBox[order(enrichBox[,7],enrichBox[,6],-(enrichBox[[2]]/enrichBox[[3]]),-enrichBox[,4],decreasing=F),]
RichFactor<-round(enrichBox_sort[[2]]/enrichBox_sort[[3]],2)
RichFactor
### significant stars:
Pvalue<-as.numeric(enrichBox_sort[[4]])
n_3tars<-which(Pvalue<=0.001)
n_2tars<-setdiff(which(Pvalue<=0.01),n_3tars)
n_1tars<-setdiff(which(Pvalue<=0.05),c(n_3tars,n_2tars))
barcols<-getMyColor(targets=Pvalue)
### pathway class:
name2type<-enrichBox_sort[,c(1,7)]
shortName<-gsub("Organismal Systems","OS",gsub("Human Diseases","HD",gsub("Genetic Information Processing","GIP",gsub("Environmental Information Processing","EIP",gsub("Drug Development","DD",gsub("Cellular Processes","CP",gsub("Metabolism","M",name2type[[2]])))))))
xlabs<-paste(enrichBox_sort[[1]],shortName,sep=": ")
typeII<-c("Pathway Class:","EIP: Environmental Information Processing","GIP: Genetic Information Processing","CP : Cellular Processes","OS : Organismal Systems","DD : Drug Development","HD : Human Diseases","M  : Metabolism")
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(10,6,4,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c(10,0.1,4,4)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
q()
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(12,8,4,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c(10,0.1,4,4)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
q()
### ploting
pdfname <-"pathway_enrichment.pdf"
pdf(file=pdfname,width=16,height=7)
skip<-max(RichFactor)/40
par(fig=c(0,0.92,0,1),mar = c(12,10,4,0.1)+0.1,new=F)
bar<-barplot(RichFactor,col=barcols,horiz=FALSE,ylim=c(0,max(RichFactor)+max(RichFactor)/3),axisnames=FALSE,ylab="EnrichmentRatio: Sample_number/Background_number",cex.lab=1,font.lab=1,col.lab="black",border=F)
text(bar,rep(-max(RichFactor)/20,length(bar)),labels=xlabs,srt=30,adj=1,xpd=T,cex=0.5,font=1)
text(bar[n_3tars],RichFactor[n_3tars]+skip,"***",cex=0.6,col="black")
text(bar[n_2tars],RichFactor[n_2tars]+skip,"**",cex=0.6,col="black")
text(bar[n_1tars],RichFactor[n_1tars]+skip,"**",cex=0.6,col="black")
legend("topleft",typeII,cex=1,bty="n",inset=0.01)
### legend for bar colors:
par(fig=c(0.91,1,0,0.9),mar = c(10,0.1,4,4)+0.1,new=T)
Bars<-barplot(rep(0.5,length(continue_cols)),horiz=T,col=continue_cols[10000:1],border=F,space=0,axes=F,xlim=c(0,1),cex.main=1,font.main=1,ylim=c(0,length(continue_cols)+100))
text(rep(0.6,5),Bars[c(200,2501,5001,7500,9800)],as.character(c(1,0.75,0.5,0.25,0)),cex=1,adj=0,xpd=T)
text(0,length(continue_cols)+500,"Pvalue",adj=0,xpd=T)
dev.off()
q()
pathways<-read.delim("enrich_pathway.xls",header=T,sep="\t",check.names=F)
plantKO<-read.delim("Kegg_plant_kopath.list",header=F)
names(pathways)
names(plantKO)
names(plantKO)<-"Id"
merge(pathways,plantKO,by="Id")->plant_pathways
dim(plant_pathways)
dim(pathways)
write.table(plant_pathways,"enrich_pathway_plant.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
history()
pathways<-read.delim("vsRice.list.1_kegg_enrichment.xls.pathway_addclass.xls",header=T,sep="\t",check.names=F)
names(pathway)
names(pathways)
history
history()
merge(pathways,plantKO,by="Id")->plant_pathways
dim(plant_pathways)
plant_pathways[1:4,]
history()
write.table(plant_pathways,"enrich_pathway_plant.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
ls()
plant_pathways[1:4,]
plant_pathways_short<-plant_pathways[,c("typeI","typeII","Id","#Term","Sample number","Genes")]
plant_pathways_short[1:4,]
plant_pathways_short_sort<-plant_pathways_short[order(as.character(plant_pathways_short[[1]]),as.character(plant_pathways_short[[2]]),as.character(plant_pathways_short[[4]])),]
plant_pathways_short_sort[1:10,]
plant_pathways_short_sort[1:10,1:5]
write.table(plant_pathways_short_sort,"plant_pathways_short_sort.xls",sep="\t",col.names=T,row.names=F,quote=F)
plant_pathways_short_sort<-plant_pathways_short[order(as.character(plant_pathways_short[[1]]),as.character(plant_pathways_short[[2]]),as.character(plant_pathways_short[[4]]),as.character(plant_pathways_short[[5]]),descreasing=T),]
?order
plant_pathways_short_sort<-plant_pathways_short[order(as.character(plant_pathways_short[[1]]),as.character(plant_pathways_short[[2]]),as.character(plant_pathways_short[[4]]),as.numeric(plant_pathways_short[[5]]),decreasing=T),]
plant_pathways_short_sort[1:10,1:5]
plant_pathways_short_sort<-plant_pathways_short[order(as.character(plant_pathways_short[[1]]),as.character(plant_pathways_short[[2]]),as.numeric(plant_pathways_short[[5]]),decreasing=T),]
plant_pathways_short_sort[1:10,1:5]
history()
write.table(plant_pathways_short_sort,"plant_pathways_short_sort.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
HDL_read.delim
HDL_read.delim("HDL.short",header=T,sep="\t")
HDL<-read.delim("HDL.short",header=T,sep="\t")
Rice<-read.delim("Rice.short",header=T,sep="\t")
ATH<-read.delim("ATH.short",header=T,sep="\t")
merge(HDL,Rice,by="Id",all=T)->11
merge(HDL,Rice,by="Id",all=T)->m1
dim(m1)
merge(HDL,Rice,by="Id",all.x=T,all.y=T)->m1
dim(m1)
dim(HDL)
dim(Rice)
merge(m1,ATH,by="Id",all.x=T,all.y=T)->m2
dim(m2)
q()
history()
ATH<-read.delim("ATH.short",header=T,sep="\t")
history()
merge(m1,ATH,by="Id",all.x=T,all.y=T)->m2
dim(m2)
m2
write.table(m2,"3orgs_enrich_path.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
histoey()
history()
HDL<-read.delim("HDL.short",header=T,sep="\t")
Rice<-read.delim("Rice.short",header=T,sep="\t")
ATH<-read.delim("ATH.short",header=T,sep="\t")
HDL<-read.delim("HDL",header=T,sep="\t")
Rice<-read.delim("Rice",header=T,sep="\t")
ATH<-read.delim("ATH",header=T,sep="\t")
ATH[1:4,]
q()
ATH<-read.delim("ATH",header=T,sep="\t")
ATH[1:4,]
history()
history(1000)
merge(HDL,Rice,by="Id",all.x=T,all.y=T)->m1
dim(m1)
dim(HDL)
dim(Rice)
merge(m1,ATH,by="Id",all.x=T,all.y=T)->m2
m2[1:4,]
names(m2)
Three_org_pathway<-m2[,c(2,3,4,1,5,6,10,11,15,16)]
Three_org_pathway[1:4,]
names(Three_org_pathway)
write.table(Three_org_pathway,"Three_org_pathway.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
history(1000)
()
Q()
q()
Three_org_pathway_new<-read.delim("Three_org_pathway.txt",sep="\t",header=T,check.names=F)
names(Three_org_pathway_new)
q()
ls()
AT[1:4,]
ATH[1:4,]
ATH_name2number<-ATH[,c(3,5)]
ATH_name2number
ATH_name2number_sort<-ATH_name2number[order(as.numeric(ATH_name2number[[2]]),decreasing=T),]
ATH_name2number_sort
numbers<-as.numeric(ATH_name2number_sort[[2]])[1:20]
lables<-as.character(ATH_name2number_sort[[1]])[1:20]
jump<-max(numbers)/10
pdfname <-"ATH_top20.pdf"
numbers
lables
numbers<-as.numeric(ATH_name2number_sort[[2]])[1:30]
lables<-as.character(ATH_name2number_sort[[1]])[1:30]
jump<-max(numbers)/10
pdfname <-"ATH_top30.pdf"
par(mar = c(10,5,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="lightblue",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Unigenes ",sep=""),cex.lab=1,font.lab=1)
text(bar,rep(-(jump/3),length(bar)),labels=labels,srt=45,adj=1,xpd=T,cex=1,font=1)
text(bar,(numbers+jump/3),labels=as.character(numbers),cex=1,font=1)
par(mar = c(10,5,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="lightblue",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Unigenes ",sep=""),cex.lab=1,font.lab=1)
text(bar,rep(-(jump/3),length(bar)),labels=labels,srt=45,adj=1,xpd=T,cex=1,font=1)
text(bar,(numbers+jump/3),labels=lables,cex=1,font=1)
par(mar = c(10,5,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="lightblue",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Unigenes ",sep=""),cex.lab=1,font.lab=1)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=1,font=1)
text(bar,(numbers+jump/3),labels=numbers,cex=1,font=1)
par(mar = c(15,6,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.6,font=1)
text(bar,(numbers+jump/3),labels=numbers,cex=0.6,font=1)
par(mar = c(10,6,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.6,font=1)
text(bar,(numbers+jump/3),labels=numbers,cex=0.6,font=1)
pdf(file=pdfname,width=12,height=8)
par(mar = c(10,6,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.6,font=1)
text(bar,(numbers+jump/3),labels=numbers,cex=0.6,font=1)
dev.off()
pdf(file=pdfname,width=12,height=8)
par(mar = c(12,6,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.8,font=1)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=1)
dev.off()
pdf(file=pdfname,width=12,height=8)
par(mar = c(13,6,2,2)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.8,font=3)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=3)
dev.off()
par(mar = c(13,6,2,2)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.8,font=2)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=2)
pdf(file=pdfname,width=12,height=8)
par(mar = c(13,6,2,2)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.8,font=2)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=2)
dev.off()
pdf(file=pdfname,width=12,height=8)
par(mar = c(15,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.8,font=2)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=2)
dev.off()
q()
history(1000)
ATH_name2number<-ATH[,c(3,5)]
ATH_name2number
ATH_name2number_sort<-ATH_name2number[order(as.numeric(ATH_name2number[[2]]),decreasing=T),]
Rice_name2number<-Rice[,c(3,5)]
Rice_name2number
Rice_name2number_sort<-Rice_name2number[order(as.numeric(Rice_name2number[[2]]),decreasing=T),]
Rice_name2number_sort
history()
names(Rice)
Rice_name2number<-Rice[,c(4,5)]
Rice_name2number_sort<-Rice_name2number[order(as.numeric(Rice_name2number[[2]]),decreasing=T),]
Rice_name2number_sort
numbers<-as.numeric(Rice_name2number_sort[[2]])[1:30]
lables<-as.character(Rice_name2number_sort[[1]])[1:30]
jump<-max(numbers)/10
pdfname <-"Rice_top30.pdf"
par(mar = c(15,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=45,adj=1,xpd=T,cex=0.8,font=2)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=2)
par(mar = c(15,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=60,adj=1,xpd=T,cex=0.8,font=2)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=2)
par(mar = c(15,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8,font=2)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=2)
par(mar = c(15,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8,font=1)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=1)
pdfname <-"Rice_top30.pdf"
pdf(file=pdfname,width=12,height=8)
par(mar = c(15,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8,font=1)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=1)
dev.off()
numbers<-as.numeric(Rice_name2number_sort[[2]])[1:30]
lables<-as.character(Rice_name2number_sort[[1]])[1:30]
jump<-max(numbers)/10
pdfname <-"Rice_top30.pdf"
pdf(file=pdfname,width=12,height=8)
par(mar = c(15,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA,axis.font=3)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8,font=3)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=3)
dev.off()
?barplot
?barplot
par(mar = c(15,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA,font=3)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8,font=3)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8,font=3)
pdfname <-"Rice_top30.pdf"
pdf(file=pdfname,width=12,height=8)
par(mar = c(12,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8)
dev.off()
numbers<-as.numeric(ATH_name2number_sort[[2]])[1:30]
lables<-as.character(ATH_name2number_sort[[1]])[1:30]
jump<-max(numbers)/10
pdfname <-"ATH_top30.pdf"
pdf(file=pdfname,width=12,height=8)
par(mar = c(12,6,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8)
dev.off()
pdf(file=pdfname,width=10,height=8)
par(mar = c(11,8,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8)
dev.off()
pdf(file=pdfname,width=10,height=6)
par(mar = c(11,8,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8)
dev.off()
pdf(file=pdfname,width=12,height=8)
par(mar = c(11,8,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8)
dev.off()
history(1000)
ATH_name2number<-ATH[,c(3,5)]
ATH_name2number
ATH_name2number_sort<-ATH_name2number[order(as.numeric(ATH_name2number[[2]]),decreasing=T),]
HDL_name2number<-HDL[,c(3,5)]
names(HDL)
HDL_name2number<-HDL[,c(4,5)]
HDL_name2number
HDL_name2number_sort<-HDL_name2number[order(as.numeric(HDL_name2number[[2]]),decreasing=T),]
HDL_name2number_sort
numbers<-as.numeric(HDL_name2number_sort[[2]])[1:30]
lables<-as.character(HDL_name2number_sort[[1]])[1:30]
jump<-max(numbers)/10
pdfname <-"HDL_top30.pdf"
par(mar = c(11,8,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8)
numbers<-as.numeric(HDL_name2number_sort[[2]])[1:30]
lables<-as.character(HDL_name2number_sort[[1]])[1:30]
jump<-max(numbers)/10
pdfname <-"HDL_top30.pdf"
pdf(file=pdfname,width=12,height=8)
par(mar = c(10,10,2,1)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="red",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Genes ",sep=""),cex.lab=1,font.lab=1,border=NA)
text(bar,rep(-(jump/3),length(bar)),labels=lables,srt=30,adj=1,xpd=T,cex=0.8)
text(bar,(numbers+jump/3),labels=numbers,cex=0.8)
dev.off()
q()
enrich_pathway<-read.delim("vsATH.list.1_kegg_enrichment.xls.pathway_addclass.xls",header=T,sep="\t",check.names=F)
q()
enrich_pathway<-read.delim("vsATH.list.1_kegg_enrichment.xls.pathway_addclass.xls",header=T,sep="\t",check.names=F)
plant_pathway<-read.delim("Kegg_plant_kopath.list",header=F)
names(enrich_pathway)
names(plant_pathway)
names(plant_pathway)<-"Id"
merge(enrich_pathway,plant_pathway,by="Id")->enrich_pathway_plant
le enrich_pathway[1:10,]
enrich_pathway[1:10,]
names(enrich_pathway)
names(enrich_pathway_plant)
enrich_pathway_plant_short<-enrich_pathway_plant[,c("typeI","typeII","#Term","Id","Sample number","Genes")]
enrich_pathway_plant_short[1:4,]
enrich_pathway_plant[1:4,]
enrich_pathway[1:4,]
q()
history()
plant_pathway<-read.delim("Kegg_plant_kopath.list.ath",header=F)
history()
names(plant_pathway)<-"Id"
merge(enrich_pathway,plant_pathway,by="Id")->enrich_pathway_plant
le enrich_pathway[1:10,]
enrich_pathway[1:10,]
names(enrich_pathway)
names(enrich_pathway_plant)
enrich_pathway_plant_short<-enrich_pathway_plant[,c("typeI","typeII","#Term","Id","Sample number","Genes")]
enrich_pathway_plant_short[1:4,]
names(enrich_pathway_plant_short)
enrich_pathway_plant_short_sort<-enrich_pathway_plant_short[order(as.character(enrich_pathway_plant_short[[1]]),as.character(enrich_pathway_plant_short[[1]]),as.character(enrich_pathway_plant_short[[1]]))]
enrich_pathway_plant_short_sort<-enrich_pathway_plant_short[order(as.character(enrich_pathway_plant_short[[1]]),as.character(enrich_pathway_plant_short[[2]]),as.numeric(enrich_pathway_plant_short[[5]]),decreasing=T),]
enrich_pathway_plant_short_sort[1:10,1:4]
enrich_pathway_plant_short_sort[1:10,1:5]
enrich_pathway_plant_short_sort[1:20,1:5]
write.table(enrich_pathway_plant_short_sort,"enrich_pathway_plant_short_sort.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
pathways<-read.delim("pathway.txt.2",sep="\t",header=F)
ORF<-read.delim("ORF.list")
names(pathways)
pathways[1:4,]
merge(ORF,pathways,by="V1")->ORF_pathways
names(ORF)
ORF<-read.delim("ORF.list",header=F)
merge(ORF,pathways,by="V1")->ORF_pathways
ORF_pathways[1:4,]
write.table(ORF_pathways
write.table(ORF_pathways,"HDL_ORF_pathway.xls",sep="\t",col.names=F,row.names=F,quote=F)
q()
enrich_pathway<-read.delim("ORF.list_kegg_enrichment.xls.pathway_addclass.xls",sep="\t",header=T,check.names=F)
plant_pathway<-read.delim("Kegg_plant_kopath.list",header=F)
names(plant_pathway)<-"Id"
merge(enrich_pathway,plant_pathway,by="Id")->enrich_pathway_plant
enrich_pathway_plant[1:10,]
names(enrich_pathway_plant0
names(enrich_pathway_plant)
enrich_pathway_plant_short<-enrich_pathway_plant[,c("typeI","typeII","Id","#Term","Sample number","Genes")]
enrich_pathway_plant_short[1:4,]
enrich_pathway_plant_short[1:4,1:4,]
enrich_pathway_plant_short_sort<-enrich_pathway_plant_short[order(as.character(enrich_pathway_plant_short[[1]]),as.character(enrich_pathway_plant_short[[4]]),as.character(enrich_pathway_plant_short[[5]]),decreasing=T),]
enrich_pathway_plant_short_sort[1:10,1:5]
enrich_pathway_plant_short_sort<-enrich_pathway_plant_short[order(as.character(enrich_pathway_plant_short[[1]]),as.character(enrich_pathway_plant_short[[2]]),as.character(enrich_pathway_plant_short[[5]]),decreasing=T),]
enrich_pathway_plant_short_sort[1:10,1:5]
write.table(enrich_pathway_plant_short_sort,"enrich_pathway_plan.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
List_1<-scan("macthed_assemble.list",what=character())
List_2<-scan("All_transcript.list",what=character())
ls()
List_1[1:4]
List_2[1:4]
List_co<-intersect(List_1,List_2)
List_specific<-setdiff(List_2,List_1)
length(List_specific)
length(List_specific)
 write.table(List_specific,"assemble_specific.list",sep="\n",col.names=F,row.names=F,quote=F)
List_specific[1:4]
List_specific_box<-as.data.frame(List_specific)
annBox<-read.delim("All.ann",header=T,sep="\t",check.names=F)
annBox[1:4,]
names(List_specific)
names(List_specific)<-"transcript_id"
names(List_specific_box)
names(List_specific_box)<-"transcript_id"
merge(List_specific_box,annBox,by="transcript_id")->List_specific_box_ann
merge(List_specific_box,annBox,by="transcript_id",all.x=T)->List_specific_box_ann
 write.table(List_specific_box_ann,"assemble_specific.list.ann",sep="\t",col.names=T,row.names=F,quote=F)
q()
common_list<-scan("macthed_assemble.list",what=character())
common_list<-as.data.frame(scan("macthed_assemble.list",what=character()))
names(common_list)
names(common_list)<-"transcript_id"
merge(common_list,annBox,by="transcript_id",all.x=T)->common_list_box_ann
common_list_box_ann[1:4,]
 write.table(common_list_box_ann,"common.list.ann",sep="\t",col.names=T,row.names=F,quote=F)
q()
plant_kpist<-read.delim("Kegg_plant_kopath.list",header=F)
1_vs_2<-read.delim("1_vs_2_kegg_enrich.xls",header=T,sep="\t")
p_1_vs_2<-read.delim("1_vs_2_kegg_enrich.xls",header=T,sep="\t")
p_1_vs_3<-read.delim("1_vs_3_kegg_enrich.xls",header=T,sep="\t")
p_2_vs_3<-read.delim("2_vs_3_kegg_enrich.xls",header=T,sep="\t")
p_1_vs_2[1:4,1:5]
names(plant_kpist)
names(plant_kpist)<-"Id"
p_1_vs_2_m<-merge(p_1_vs_2,plant_kpist,by="Id")
dim(p_1_vs_2_m)
dim(p_1_vs_2)
p_1_vs_3_m<-merge(p_1_vs_3,plant_kpist,by="Id")
p_2_vs_3_m<-merge(p_2_vs_3,plant_kpist,by="Id")
write.table(p_1_vs_2_m,"1_vs_2_kegg_enrichment_plant.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(p_1_vs_3_m,"1_vs_3_kegg_enrichment_plant.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(p_2_vs_3_m,"2_vs_3_kegg_enrichment_plant.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
aa<-read.delim("1_vs_2.sdegs",header=T)
bb<-read.delim("Rice_blast_function.xls.xls",sep="\t",header=T,check.names=D)
bb<-read.delim("Rice_blast_function.xls.xls",sep="\t",header=T,check.names=F)
merge(aa,bb,by="transcript_id")->cc
cc[1:4,]
dim(cc)
write.table(cc,"1_vs_2.sdegs.Rice",col.names=T,row.names=F,quote=F)
write.table(cc,"1_vs_2.sdegs.Rice",col.names=T,row.names=F,quote=F,sep="\t")
aa<-read.delim("1_vs_3.sdegs",header=T)
 merge(aa,bb,by="transcript_id")->cc
dim(cc)
write.table(cc,"1_vs_3.sdegs.Rice",col.names=T,row.names=F,quote=F,sep="\t")
ls()
aa[1:4,]
bb[1:4,]
Rice<-bb
Ara<-read.delim("Ara_blast_function.xls.xls",sep="\t",header=T,check.names=F)
 merge(Rice,Ara,by="transcript_id")->Rice_Ara
 merge(Rice,Ara,by="transcript_id",all=T)->Rice_Ara
 merge(Rice,Ara,by="transcript_id",all.x=T,all.y=T)->Rice_Ara
 merge(Rice,Ara,by="transcript_id")->Rice_Ara
names(Rice)
names(Ara)
 merge(Rice,Ara,by="transcript_id",all.x=T)->Rice_Ara
 merge(Rice,Ara,by="transcript_id",all.x=T,all.y=T)->Rice_Ara
?merge
 merge(Rice,Ara,by="transcript_id",all=T)->Rice_Ara
names(Rice)
Rice1<-Rice[,-9]
 merge(Rice1,Ara,by="transcript_id",all=T)->Rice_Ara
Rice_Ara[1:4,]
all_isigene<-read.delim("26565_isogenes.list",header=T)
 merge(all_isigene,Rice_Ara,by="transcript_id",all=T)->Rice_Ara_isogene
Rice_Ara_isogene[1:4,]
Rice_Ara_isogene[1:40,]
write.table(Rice_Ara_isogene,"Rice_Ara_isogene",sep="\t",col.names=T,row.names=F,quote=F)
q()
names(Rice_Ara_isogene)
Rice_Ara_isogene_sort<-Rice_Ara_isogene[order(Rice_Ara_isogene[,2]),]
Rice_Ara_isogene_sort[1:4,]
Rice_Ara_isogene_sort[1:40,]
write.table(Rice_Ara_isogene_sort,"Rice_Ara_isogene",sep="\t",col.names=T,row.names=F,quote=F)
names(Rice_Ara_isogene_sort)
Rice_Ara_isogene_sort_short<-Rice_Ara_isogene_sort[,c(1,2,3,9,16)]
Rice_Ara_isogene_sort_short[1:50,]
Rice_Ara_isogene_sort_short
unique(Rice_Ara_isogene_sort_short[[5]])
write.table(Rice_Ara_isogene_sort_short,"Rice_Ara_isogene_short",sep="\t",col.names=T,row.names=F,quote=F)
q()
hisory
history(100)
names(Rice_Ara_isogene_sort)
Rice[1:10,]
Rice_Ara_isogene_sort_short<-Rice_Ara_isogene_sort[,c(1,2,4,9,16)]
Rice_Ara_isogene_sort_short[1:4,]
write.table(Rice_Ara_isogene_sort_short,"Rice_Ara_isogene_short",sep="\t",col.names=T,row.names=F,quote=F)
q()
history(2000)
mm<-read.delim("2_vs_3.sdegs",header=F)
mm
history(2000)
ls()
names(bb)
names(mm)
names(mm)<-"transcript_id"
history(2000)
merge(mm,bb,by="transcript_id")->nn
nn[1:4,]
history(2000)
write.table(nn,"2_vs_3.sdegs.Rice",col.names=T,row.names=F,quote=F,sep="\t")
q()
history()
history(1000)
ls()
names(aa)
names(Ar)
names(Ara)
merge(mm,Ara,by="transcript_id")->nn
nn[1:10,]
write.table(nn,"2_vs_3.sdegs.Ara",col.names=T,row.names=F,quote=F,sep="\t")
q()
Pfam<-read.delim("Trinotatepfam.out.new.tab.0.1.title",sep="\t",header=T,check.names=F)
q()
Pfam<-read.delim("Trinotatepfam.out.new.tab.0.1.title",sep="\t",header=T,check.names=F)
q()
Pfam<-read.delim("Trinotatepfam.out.new.tab.0.1.title",sep="\t",header=T,check.names=F)
Pfam[1:4,]
q
()
q()
ls()
Pfam[1:4,]
ls()
Pfam_short<-Pfam[,c(1:2,4,7,8)]
Pfam_short[1:4,]
Pfam_short<-Pfam[,c(1:2,4,7,8,23)]
Pfam_short[1:4,]
Pfam_short_sort<-Pfam_short[order(Pfam_short[,3],Pfam_short[,4],Pfam_short[,5]),]
Pfam_short_sort[1:4,]
Pfam_short_sort[1:10,]
Pfam_short_sort<-Pfam_short[order(Pfam_short[,3],-Pfam_short[,4],Pfam_short[,5],decreasing=F),]
Pfam_short_sort[1:10,]
Pfam_short_sort<-Pfam_short[order(Pfam_short[,3],-Pfam_short[,4],Pfam_short[,5],decreasing=T),]
Pfam_short_sort[1:20,]
Gene<-as.character(Pfam_short_sort[[3]])
Gene<-unique(as.character(Pfam_short_sort[[3]]))
length(Gene)
Hit_line<-vector()
for(i in 1:length(Gene)){
Hit_line[i]<-which(Pfam_short_sort[[3]]==Gene[i])[1]
}
Hit_line
ls()
Pfam_short_sort_best<-Pfam_short_sort[Hit_line,]
Pfam_short_sort_best[1:4,]
write.table(Pfam_short_sort_best,"Pfam_best_evalue0.1",sep="\t",col.names=T,row.names=F)
q()
write.table(Pfam_short_sort_best,"Pfam_best_evalue0.1",sep="\t",col.names=T,row.names=F,quote=F)
q()
Box<-read.delim("Pfam_type2gene",sep="\t",header=F)
Box[1:4,]
Pfams<-as.character(Box[[1]])
Pfams
Box[1:4,]
Box_new，-data.frame()
Box_new<-data.frame()
for(i in 1:length(Pfams)){
tmp<-which(Box[,1]==Pfams[i])
Box_new[i,1]<-Box[tmp[1],1]
Box_new[i,2]<-Box[tmp[1],2]
Box_new[i,3]<-Box[tmp[1],3]
Box_new[i,4]<-paste(as.character(Box[tmp,4]),collapse=";")
}
Box_new[1:10,]
length(Pfams)
Pfams[1:10]
history()
Pfams<-unique(as.character(Box[[1]]))
history()
Box_new<-data.frame()
for(i in 1:length(Pfams)){
tmp<-which(Box[,1]==Pfams[i])
Box_new[i,1]<-Box[tmp[1],1]
Box_new[i,2]<-Box[tmp[1],2]
Box_new[i,3]<-Box[tmp[1],3]
Box_new[i,4]<-paste(as.character(Box[tmp,4]),collapse=";")
}
Box_new[1:10,]
history()
Box_new<-data.frame()
for(i in 1:length(Pfams)){
tmp<-which(Box[,1]==Pfams[i])
Box_new[i,1]<-Box[tmp[1],1]
Box_new[i,2]<-Box[tmp[1],2]
Box_new[i,3]<-Box[tmp[1],3]
Box_new[i,4]<-length(tmp)
 Box_new[i,5]<-paste(as.character(Box[tmp,4]),collapse=";")
}
Box_new[1:4,]
Box_new_sort<-Box_new[order(Box_new[[4]],decreasing=T),]
Box_new_sort[1:4,]
Box_new_sort[1:4,1:4]
Box_new_sort[1:40,1:4]
Box_new_sort[1:10,1:4]
write.table(Box_new_sort,"Pfam_distribution.xls",sep="\t",col.names=T,row.names=F,quote=F)
dim(Box_new_sort)
BarMat<-Box_new_sort[,3:4]
BarMat[1:4,]
alls<-sum(BarMat[[2]])
alls
BarMat_n<-as.numeric(BarMat[[2]])/alls
BarMat_n
BarMat_n[1:20]
BarMat_n<-as.numeric(BarMat[[2]])*100/alls
BarMat_n
BarMat_n<-(as.numeric(BarMat[[2]])*100/alls)[1:30]
BarMat_n
BarMat[1:30,]
BarMat_new<-BarMat[1:30,]
BarMat_new[31,1]<-"Others"
class(BarMat)
BarMat_new
BarMat_new<-BarMat[1:30,
BarMat_new
names(BarMat_new)
rownames(BarMat_new)
rownames(BarMat_new)<-as.character(BarMat_new[[1]])
rownames(BarMat_new)<-as.character(BarMat_new[,1])
history()
BarMat_new
BarMat
BarMat_name<-unique(as.character(BarMat[[1]]))
BarMat_name_num<-as.data.frame()
BarMat_name_num<-data.frame()
for(i in 1:length(BarMat_name)){
BarMat_name_num[i,1]<-BarMat_name[i]
BarMat_name_num[i,2]<-sum(BarMat[which(BarMat[[1]]==BarMat_name[i]),2])
}
BarMat_name_num
dim(BarMat_name_num)
dim(BarMat)
BarMat_name_num[1:10,]
plot_mat<-BarMat_name_num[1:30,]
plot_mat
plot_mat[31,1]<-"Others"
plot_mat[31,2]<-sum(BarMat_name_num[31:nrow(BarMat_name_num),2])
plot_mat
BarMat_name_num
BarMat_name_num[1:100,]
BarMat_name_num_sort<-BarMat_name_num[order(BarMat_name_num[[2]],decreasing=T),]
BarMat_name_num_sort[1:20,]
plot_mat<-BarMat_name_num[1:60,]
plot_mat[61,2]<-sum(BarMat_name_num_sort[61:nrow(BarMat_name_num_sort),2])
plot_mat
plot_mat[61,1]<-"Others"
plot_mat
BarMat_name_num_sort<-BarMat_name_num[order(as.numeric(BarMat_name_num[[2]]),decreasing=T),]
BarMat_name_num_sort[1:60,]
plot_mat<-BarMat_name_num_sort[1:60,]
plot_mat[61,1]<-"Others"
plot_mat[61,2]<-sum(BarMat_name_num_sort[61:nrow(BarMat_name_num_sort),2])
plot_mat
barplot(plot_mat[[2]])
plot_mat<-BarMat_name_num_sort[1:60,]
barplot(plot_mat[[2]])
numbers<-as.numeric(plot_mat[[2]])
labels<-as.character(plot_mat[[1]])
jump<-max(numbers)/10
par(mar = c(12,6,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="darkred",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Unigenes ",sep=""),cex.lab=1,font.lab=1)
text(bar,rep(-(jump/3),length(bar)),labels=labels,srt=45,adj=1,xpd=T,cex=1,font=1)
text(bar,(numbers+jump/3),labels=as.character(numbers),cex=1,font=1)
par(mar = c(12,6,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="darkred",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Unigenes ",sep=""),cex.lab=1,font.lab=1)
text(bar,rep(-(jump/3),length(bar)),labels=labels,srt=30,adj=1,xpd=T,cex=0.6,font=1)
text(bar,(numbers+jump/3),labels=as.character(numbers),cex=0.6,font=1)
par(mar = c(12,10,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="darkred",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Unigenes ",sep=""),cex.lab=1,font.lab=1)
text(bar,rep(-(jump/3),length(bar)),labels=labels,srt=30,adj=1,xpd=T,cex=0.6,font=1)
text(bar,(numbers+jump/3),labels=as.character(numbers),cex=0.6,font=1)
pdfname <-"Pfam_top60_distribution.pdf"
pdf(file=pdfname,width=12,height=8)
par(mar = c(12,10,3,3)+0.1)
ylims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="darkred",horiz=FALSE,ylim=ylims,axisnames=FALSE,ylab=paste(" Number of Unigenes ",sep=""),cex.lab=1,font.lab=1)
text(bar,rep(-(jump/3),length(bar)),labels=labels,srt=30,adj=1,xpd=T,cex=0.6,font=1)
text(bar,(numbers+jump/3),labels=as.character(numbers),cex=0.6,font=1)
dev.off()
ls()
BarMat_name_num_sort[1:4,]
BarMat[1:4,]
Pfam_short_sort_best[1:4,]
Box_new_sort[1:4,]
Box_new_sort[1:4,1:4]
Box_new_sort[1:40,1:4]
write.table(Box_new_sort,"Pfam_ID_name_num_list.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
ls()
Pfam_12_Pfam2Num2Gene<-read.delim("")
Pfam_12_Pfam2Num2Gene<-read.delim("Pfam_12_Pfam2Num2Gene.xls",header=F,sep="\t")
Pfam_12_Pfam2Num2Gene[1:4,]
Pfam_12_Pfam2Num2Gene[1:4,1:4]
Pfam_12_Pfam2Num2Gene[1:4,1:2]
Pfam_12_Pfam2Num2Gene<-read.delim("Pfam_12_Pfam2Num2Gene.xls",header=T,sep="\t")
Pfam_12_Pfam2Num2Gene[1:4,1:2]
Pfam_12_Pfam2Num2Gene<-read.delim("Pfam_12_Pfam2Num2Gene.xls",header=T,sep="\t")
Pfam_12_Pfam2Num<-Pfam_12_Pfam2Num2Gene[,1:2]
Pfam_12_Pfam2Num
Pfam_12_Pfam2Num[1:20,]
Pfam_12_Pfam2Num<-Pfam_12_Pfam2Num2Gene[-1,1:2]
Pfam_12_Pfam2Num[1:20,]
barplot(Pfam_12_Pfam2Num[[2]])
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30])
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30])
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="lightblue",border=NA)
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30])
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA)
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+0.2,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=0.5)
par(mar=c(10,10,3,3))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5))
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+0.2,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=0.5)
par(mar=c(10,10,3,3))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5))
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+0.2,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=0.2)
par(mar=c(10,10,3,3))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5))
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.6)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+0.2,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=0.5)
pdf("1_vs_2_sdeg_top30_Pfam_distribution.pdf",w=10,h=8)
par(mar=c(10,5,3,3))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5))
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.6)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=0.5)
dev.off()
pdf("1_vs_2_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(10,5,3,1))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=1)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=1)
dev.off()
?title
par(mar=c(10,5,3,1))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=1)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for griup 1_vs_2")
pdf("1_vs_2_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(10,5,3,1))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=1)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for group 1_vs_2")
dev.off()
Pfam_12_Pfam2Num2Gene<-read.delim("Pfam_12_Pfam2Num2Gene.xls",header=T,sep="\t")
Pfam_12_Pfam2Num<-Pfam_12_Pfam2Num2Gene[-1,1:2]
pdf("1_vs_2_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(10,5,3,1))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.8)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for group 1_vs_2")
dev.off()
?title
pdf("1_vs_2_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(12,6,3,1))
bars<-barplot(Pfam_12_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_12_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_12_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.8)
text(bars,Pfam_12_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_12_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for group 1_vs_2",cex.main=2)
dev.off()
Pfam_13_Pfam2Num2Gene<-read.delim("Pfam_13_Pfam2Num2Gene.xls",header=T,sep="\t")
Pfam_13_Pfam2Num<-Pfam_13_Pfam2Num2Gene[-1,1:2]
pdf("1_vs_3_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(12,6,3,1))
bars<-barplot(Pfam_13_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_13_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_13_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.8)
text(bars,Pfam_13_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_13_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for group 1_vs_3",cex.main=2)
dev.off()
Pfam_13_Pfam2Num2Gene<-read.delim("Pfam_13_Pfam2Num2Gene.xls",header=T,sep="\t")
Pfam_13_Pfam2Num<-Pfam_13_Pfam2Num2Gene[-1,1:2]
pdf("1_vs_3_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(12,10,3,1))
bars<-barplot(Pfam_13_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_13_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_13_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.8)
text(bars,Pfam_13_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_13_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for group 1_vs_3",cex.main=2)
dev.off()
Pfam_13_Pfam2Num2Gene<-read.delim("Pfam_13_Pfam2Num2Gene.xls",header=T,sep="\t")
Pfam_13_Pfam2Num<-Pfam_13_Pfam2Num2Gene[-1,1:2]
pdf("1_vs_3_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(12,10,3,5))
bars<-barplot(Pfam_13_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_13_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_13_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.6)
text(bars,Pfam_13_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_13_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for group 1_vs_3",cex.main=2)
dev.off()
Pfam_23_Pfam2Num2Gene<-read.delim("Pfam_23_Pfam2Num2Gene.xls",header=T,sep="\t")
Pfam_23_Pfam2Num<-Pfam_23_Pfam2Num2Gene[-1,1:2]
pdf("2_vs_3_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(12,10,3,5))
bars<-barplot(Pfam_23_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_23_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_23_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.6)
text(bars,Pfam_23_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_23_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for group 2_vs_3",cex.main=2)
dev.off()
Pfam_23_Pfam2Num2Gene<-read.delim("Pfam_23_Pfam2Num2Gene.xls",header=T,sep="\t")
Pfam_23_Pfam2Num<-Pfam_23_Pfam2Num2Gene[-1,1:2]
pdf("2_vs_3_sdeg_top30_Pfam_distribution.pdf",w=12,h=8)
par(mar=c(12,10,3,5))
bars<-barplot(Pfam_23_Pfam2Num[[2]][1:30],col="#6BA0D2",border=NA,ylim=c(0,max(Pfam_23_Pfam2Num[[2]][1:30])+5),ylab="Number of unigenes")
text(bars,rep(-0.5,30),Pfam_23_Pfam2Num[[1]][1:30],srt=30,adj=1,xpd=T,cex=0.8)
text(bars,Pfam_23_Pfam2Num[[2]][1:30]+1,as.numeric(Pfam_23_Pfam2Num[[2]][1:30]),cex=1)
title(main="Pfam distribution of SDEGs for group 2_vs_3",cex.main=2)
dev.off()
q()
plant_list<-read.delim("Kegg_plant_kopath.list.new",header=T)
annotation<-read.delim("pathway_table.xls",header=T,sep="\t")
names(plant_list)
plant_pathways<-merge(annotation,plant_list,by="PathWay")
dim(plant_pathways)
dim(annotation)
write.table(plant_pathways,"pathway_table.xls.plant.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
Ara_blast<-read.delim("Ara_blast.list",header=T,sep="\t")
Ara_function<-read.delim("Ara_function.list",header=T,sep="\t")
names(Ara_function)
names(Ara_function)[1]<-names(Ara_blast)[2]
names(Ara_function)
merge(Ara_blast,Ara_function,by="Ara_transcript_id",all.x=T)->Ara_blast_function
Ara_blast_function[1:4,]
write.table(Ara_blast_function,"Ara_blast_function.list",sep="\t",col.names=T,row.names=F,quote=F)
q()
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0)
text(2,7.5,"diphosphate synthase activity",adj=0)
text(2,14.5,"oxidoreductase activity,",adj=0)
text(2,13.5,"acting on CH or CH2 groups,",adj=0)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=30)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=30)
par(fig=c(0.85,1,0,1),new=TRUE,mar=c(4,1,4,4))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA)
box()
par(fig=c(0,0.8,0,1),new=F,mar=c(4,4,4,1))
plot(c(1,35),c(1,25),type="n",xlab="",ylab="",ylim=c(0,25),xlim=c(0,35),axes=T)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0)
text(2,7.5,"diphosphate synthase activity",adj=0)
text(2,14.5,"oxidoreductase activity,",adj=0)
text(2,13.5,"acting on CH or CH2 groups,",adj=0)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=30)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=30)
par(fig=c(0.85,1,0,1),new=TRUE,mar=c(4,1,4,4))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA)
box()
par(fig=c(0,0.8,0,1),new=F,mar=c(4,4,4,1))
plot(c(1,35),c(1,25),type="n",xlab="",ylab="",ylim=c(0,25),xlim=c(0,35),axes=T)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0)
text(2,7.5,"diphosphate synthase activity",adj=0)
text(2,14.5,"oxidoreductase activity,",adj=0)
text(2,13.5,"acting on CH or CH2 groups,",adj=0)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=30)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=30)
par(fig=c(0.8,1,0,1),new=TRUE,mar=c(4,1,4,4))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA)
par(fig=c(0,0.8,0,1),new=F,mar=c(4,4,4,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=T)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0)
text(2,7.5,"diphosphate synthase activity",adj=0)
text(2,14.5,"oxidoreductase activity,",adj=0)
text(2,13.5,"acting on CH or CH2 groups,",adj=0)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40)
par(fig=c(0.8,1,0,1),new=TRUE,mar=c(4,1,4,4))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA)
as.matrix(rep(1,length(Colors)))
as.matrix(rep(1,length(Colors)*2))
as.matrix(rep(1,length(Colors)*2),17,2,byrow=T)
as.matrix(rep(1,length(Colors)*2),17,2)
rep(1,length(Colors)*2)
as.matrix(rep(1,length(Colors)*2),nrow=17,ncol=2)
matrix(rep(1,length(Colors)*2),nrow=17,ncol=2)
barplot(matrix(rep(1,length(Colors)*2),nrow=17,ncol=2))
barplot(matrix(rep(1,length(Colors)*2),nrow=17,ncol=2),col=c(rep("red",17),rep("blue",17)))
box()
barplot(matrix(rep(1,length(Colors)*2),nrow=17,ncol=2),col=c(rep("red",17),rep("blue",17)),xlim=c(1,5))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,4))
box()
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,4))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=T,border=NA,xlim=c(1,4))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,4))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=T,border=NA,xlim=c(1,4))
Sort_values
text(rep(1,17),Sort_values)
text(rep(2,17),Sort_values)
text(rep(2,17),Sort_values,adj=0)
text(rep(2,17),1:17,as.character(Sort_values),adj=0)
text(rep(2,17),seq(0.5,16.5,by=0.5),as.character(Sort_values),adj=0)
Sort_values
seq(0.5,16.5,by=0.5)
seq(0.5,16.5,by=1)
 barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=T,border=NA,xlim=c(1,4))
text(rep(2,17),seq(0.5,16.5,by=1),as.character(Sort_values),adj=0)
text(rep(2,17),seq(0.5,16.5,by=1),as.character(Sort_values),adj=0,cex=0.8)
 barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=T,border=NA,xlim=c(1,4))
text(rep(2,17),seq(0.5,16.5,by=1),as.character(Sort_values),adj=0,cex=0.8)
 barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=T,border=NA,xlim=c(1,6))
text(rep(2,17),seq(0.5,16.5,by=1),as.character(Sort_values),adj=0,cex=0.8)
par(fig=c(0,0.8,0,1),new=F,mar=c(4,4,4,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0,cex=0.8)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=0.8)
text(2,7.5,"diphosphate synthase activity",adj=0,cex=0.8)
text(2,14.5,"oxidoreductase activity,",adj=0,cex=0.8)
text(2,13.5,"acting on CH or CH2 groups,",adj=0,cex=0.8)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=0.8)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40)
par(fig=c(0.8,1,0,1),new=TRUE,mar=c(4,1,4,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
text(rep(2,17),seq(0.5,16.5,by=1),as.character(Sort_values),adj=0,cex=0.8)
par(fig=c(0,0.8,0,1),new=F,mar=c(4,4,4,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0,cex=0.8)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=0.8)
text(2,7.5,"diphosphate synthase activity",adj=0,cex=0.8)
text(2,14.5,"oxidoreductase activity,",adj=0,cex=0.8)
text(2,13.5,"acting on CH or CH2 groups,",adj=0,cex=0.8)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=0.8)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40)
par(fig=c(0.8,1,0.2,0.8),new=TRUE,mar=c(4,1,4,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
text(rep(2,17),seq(0.5,16.5,by=1),as.character(Sort_values),adj=0,cex=0.8)
par(fig=c(0,0.8,0,1),new=F,mar=c(4,4,4,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0,cex=0.8)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=0.8)
text(2,7.5,"diphosphate synthase activity",adj=0,cex=0.8)
text(2,14.5,"oxidoreductase activity,",adj=0,cex=0.8)
text(2,13.5,"acting on CH or CH2 groups,",adj=0,cex=0.8)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=0.8)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40,cex=0.8)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40,cex=0.8)
par(fig=c(0.8,1,0.2,0.8),new=TRUE,mar=c(4,1,4,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
text(rep(2,17),seq(0.5,16.5,by=1),as.character(Sort_values),adj=0,cex=0.8)
par(fig=c(0,0.8,0,1),new=F,mar=c(4,1,4,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(1,2.5,GOs[1],,adj=0,cex=1)
text(1,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=1)
text(1,7.5,"diphosphate synthase activity",adj=0,cex=1)
text(1,14.5,"oxidoreductase activity,",adj=0,cex=1)
text(1,13.5,"acting on CH or CH2 groups,",adj=0,cex=1)
text(1,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=1)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40,cex=1)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40,cex=1)
par(fig=c(0.8,1,0.2,0.8),new=TRUE,mar=c(4,1,4,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
text(rep(2,17),seq(0.5,16.5,by=1),as.character(Sort_values),adj=0,cex=0.8)
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
pdf("Enriched GO functional categories in each of pair-wise comparison.pdf",w=20,h=18)
par(fig=c(0,0.8,0,1),new=F,mar=c(4,1,4,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(1,2.5,GOs[1],,adj=0,cex=1)
text(1,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=1)
text(1,7.5,"diphosphate synthase activity",adj=0,cex=1)
text(1,14.5,"oxidoreductase activity,",adj=0,cex=1)
text(1,13.5,"acting on CH or CH2 groups,",adj=0,cex=1)
text(1,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=1)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40,cex=1)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40,cex=1)
par(fig=c(0.8,1,0.2,0.8),new=TRUE,mar=c(4,1,4,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
#text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
dev.off()
pdf("Enriched GO functional categories in each of pair-wise comparison.pdf",w=12,h=10)
par(fig=c(0,0.8,0,1),new=F,mar=c(4,1,4,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(1,2.5,GOs[1],,adj=0,cex=1)
text(1,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=1)
text(1,7.5,"diphosphate synthase activity",adj=0,cex=1)
text(1,14.5,"oxidoreductase activity,",adj=0,cex=1)
text(1,13.5,"acting on CH or CH2 groups,",adj=0,cex=1)
text(1,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=1)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40,cex=1)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40,cex=1)
par(fig=c(0.8,1,0.2,0.8),new=TRUE,mar=c(4,1,4,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
#text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
dev.off()
pdf("Enriched GO functional categories in each of pair-wise comparison.pdf",w=12,h=10)
par(fig=c(0,0.8,0,1),new=F,mar=c(2,1,2,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0,cex=1.5)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=1.5)
text(2,7.5,"diphosphate synthase activity",adj=0,cex=1.5)
text(2,14.5,"oxidoreductase activity,",adj=0,cex=1.5)
text(2,13.5,"acting on CH or CH2 groups,",adj=0,cex=1.5)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=1.5)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40,cex=1.5)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40,cex=1.5)
par(fig=c(0.8,1,0.3,0.9),new=TRUE,mar=c(2,1,2,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
#text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
dev.off()
pdf("Enriched GO functional categories in each of pair-wise comparison.pdf",w=12,h=10)
par(fig=c(0,0.8,0,1),new=F,mar=c(2,1,2,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0,cex=1.5)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=1.5)
text(2,7.5,"diphosphate synthase activity",adj=0,cex=1.5)
text(2,14.5,"oxidoreductase activity,",adj=0,cex=1.5)
text(2,13.5,"acting on CH or CH2 groups,",adj=0,cex=1.5)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=1.5)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40,cex=1.5)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40,cex=1.5)
par(fig=c(0.8,1,0.2,0.8),new=TRUE,mar=c(2,1,2,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
#text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
dev.off()
pdf("Enriched GO functional categories in each of pair-wise comparison.pdf",w=12,h=10)
par(fig=c(0,0.8,0,1),new=F,mar=c(2,1,2,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0,cex=1.5)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=1.5)
text(2,7.5,"diphosphate synthase activity",adj=0,cex=1.5)
text(2,14.5,"oxidoreductase activity,",adj=0,cex=1.5)
text(2,13.5,"acting on CH or CH2 groups,",adj=0,cex=1.5)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=1.5)
text(0.5,15.5,"GO:0052592",adj=0)
text(0.5,9.5,"GO:0046429",adj=0)
text(0.5,3.5,"GO:0019139",adj=0)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40,cex=1.5)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40,cex=1.5)
par(fig=c(0.8,1,0.2,0.78),new=TRUE,mar=c(2,1,2,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
#text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
dev.off()
group1_enrich
group1_enrich[,1:4]
pdf("Enriched GO functional categories in each of pair-wise comparison.pdf",w=12,h=10)
par(fig=c(0,0.8,0,1),new=F,mar=c(2,1,2,1))
plot(c(1,20),c(1,20),type="n",xlab="",ylab="",ylim=c(0,20),xlim=c(0,35),axes=F)
#rect(xleft_1,ybottom_1,xright_1,ytop_1,col=col_1,border="black")
rect(xleft_1+20,ybottom_1,xright_1+20,ytop_1,col=col_1,border="black")
## plot group2 gos:
xleft_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xleft_2"])
ybottom_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ybottom_2"])
xright_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"xright_2"])
ytop_2<-as.numeric(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"ytop_2"])
col_2<-as.character(Final_Plot_Data_sort[which(Final_Plot_Data_sort[["group2"]]==1),"col_2"])
#rect(xleft_2,ybottom_2,xright_2,ytop_2,col=col_2,border="black")
rect(xleft_2+20,ybottom_2,xright_2+20,ytop_2,col=col_2,border="black")
GOs<-unique(as.character(Final_Plot_Data_sort[["description"]]))
segments(0,4,34,4,lwd=1,col="black")
segments(0,10,34,10,lwd=1,col="black")
segments(0,16,34,16,lwd=1,col="black")
text(2,2.5,GOs[1],,adj=0,cex=1.5)
text(2,8.5,"4-hydroxy-3-methylbut-2-en-1-yl",adj=0,cex=1.5)
text(2,7.5,"diphosphate synthase activity",adj=0,cex=1.5)
text(2,14.5,"oxidoreductase activity,",adj=0,cex=1.5)
text(2,13.5,"acting on CH or CH2 groups,",adj=0,cex=1.5)
text(2,12.5,"with an iron-sulfur protein as acceptor",adj=0,cex=1.5)
text(0.5,15.5,"GO:0052592",adj=0,cex=1.5)
text(0.5,9.5,"GO:0046429",adj=0,cex=1.5)
text(0.5,3.5,"GO:0019139",adj=0,cex=1.5)
## 添加group名称
text(24,16.5,"XCJ2 - XCJ1",adj=0,srt=40,cex=1.5)
text(30,16.5,"XCJ3 - XCJ1",adj=0,srt=40,cex=1.5)
par(fig=c(0.8,1,0.2,0.78),new=TRUE,mar=c(2,1,2,1))
barplot(as.matrix(rep(1,length(Colors))),col=Colors,axes=F,border=NA,xlim=c(1,6))
#text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
text(rep(2,17),seq(0.5,16.5,by=1),as.character(round(Sort_values,2)),adj=0,cex=0.8)
dev.off()
?image
 require(grDevices) 
x <- y <- seq(-4*pi, 4*pi, len = 27)
x
y
 r <- sqrt(outer(x^2, y^2, "+"))
r
class(r)
dim(r)
image(z = z <- cos(r^2)*exp(-r/6), col  = gray((0:32)/32))
test<-r[1:10,1:10]
test
test<-r[1:4,1:4]
test
test[1,2:4]<-NULL
test[1,2:4]<-0
test
test[2,3:4]<-0
test[3,34]<-0
test[3,4]<-0
test
test[4,4]<-0
test[3,3]<-0
test[2,2]<-0
test
z = z <- cos(r^2)*exp(-r/6)
z
?image
Final_Plot_Data_sort
Final_Plot_Data_sort[["transcript"]]
as.character(Final_Plot_Data_sort[["transcript"]])
Final_Plot_Data_sort
write.table(Final_Plot_Data_sort,"CoEnrichedGOs_DetailInfor.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
dir()
XCJ_1SpecificExpGene<-read.delim("XCJ_1_SpecificExpGene",header=F,sep="\t")
vsAtTF_Annotation<-read.delim("vsAtTF_Annotation.txt",sep="\t",header=T)
vsAtTF_Annotation[1:4,]
vsAtTF_Annotation<-read.delim("vsAtTF_Annotation.txt",sep="\t",header=T,check.names=F)
vsAtTF_Annotation[1:4,]
names(XCJ_1SpecificExpGene)
names(XCJ_1SpecificExpGene)<-c("TranscriptName","XCJ-1_vs_XCJ-2","XCJ-1_vs_XCJ-3","XCJ-2_vs_XCJ-3","SUM","XCJ-1_count","XCJ-2_count","XCJ-3_count","XCJ-1_fpkm","XCJ-2_fpkm","XCJ-3_fpkm")
merge(XCJ_1SpecificExpGene,vsAtTF_Annotation,by="TranscriptName")->XCJ_1SpecificExpGene_TF
XCJ_1SpecificExpGene_TF[1:4,]
dim(XCJ_1SpecificExpGene_TF)
q()
dir()
XCJ_1_ExpTranscript<-read.delim("XCJ_1_ExpTranscript",header=F,sep="\t")
XCJ_1_ExpTranscript[1:4,]
dim(XCJ_1_ExpTranscript)
XCJ_1_ExpTranscript<-read.delim("XCJ_1_ExpTranscript",header=T,sep="\t",check.names=F)
XCJ_1_ExpTranscript[1:4,]
XCJ_2_ExpTranscript<-read.delim("XCJ_2_ExpTranscript",header=T,sep="\t",check.names=F)
XCJ_3_ExpTranscript<-read.delim("XCJ_3_ExpTranscript",header=T,sep="\t",check.names=F)
names(XCJ_1_ExpTranscript)
ls()
 names(XCJ_1_ExpTranscript)[1]<-names(vsAtTF_Annotation)[1]
names(XCJ_1_ExpTranscript)[1]
merge(XCJ_1_ExpTranscript,vsAtTF_Annotation,by="TranscriptName")->XCJ_1_ExpTranscript_TF
dim(XCJ_1_ExpTranscript_TF)
merge(XCJ_2_ExpTranscript,vsAtTF_Annotation,by="TranscriptName")->XCJ_2_ExpTranscript_TF
 names(XCJ_2_ExpTranscript)[1]<-names(vsAtTF_Annotation)[1]
 names(XCJ_3_ExpTranscript)[1]<-names(vsAtTF_Annotation)[1]
merge(XCJ_2_ExpTranscript,vsAtTF_Annotation,by="TranscriptName")->XCJ_2_ExpTranscript_TF
merge(XCJ_3_ExpTranscript,vsAtTF_Annotation,by="TranscriptName")->XCJ_3_ExpTranscript_TF
dim(XCJ_2_ExpTranscript_TF)
dim(XCJ_3_ExpTranscript_TF)
XCJ_1_ExpTranscript_TF[1:4,]
names(XCJ_1_ExpTranscript_TF)
XCJ_1_ExpTranscript_TF_new<-XCJ_1_ExpTranscript_TF[,c(1,6:12,17:19)]
XCJ_2_ExpTranscript_TF_new<-XCJ_2_ExpTranscript_TF[,c(1,6:12,17:19)]
XCJ_3_ExpTranscript_TF_new<-XCJ_3_ExpTranscript_TF[,c(1,6:12,17:19)]
XCJ_3_ExpTranscript_TF_new[1:4,]
XCJ_3_ExpTranscript_TF_new[1:40,]
ls()
write.table(XCJ_3_ExpTranscript_TF,"XCJ_3_ExpTranscript_TF.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(XCJ_2_ExpTranscript_TF,"XCJ_2_ExpTranscript_TF.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(XCJ_1_ExpTranscript_TF,"XCJ_1_ExpTranscript_TF.xls",sep="\t",col.names=T,row.names=F,quote=F)
ls()
XCJ_1_ExpTranscript_TF_list<-as.character(XCJ_1_ExpTranscript_TF[[1]])
XCJ_2_ExpTranscript_TF_list<-as.character(XCJ_2_ExpTranscript_TF[[1]])
XCJ_3_ExpTranscript_TF_list<-as.character(XCJ_3_ExpTranscript_TF[[1]])
length(intersect(XCJ_1_ExpTranscript_TF_list,XCJ_2_ExpTranscript_TF_list))
length(intersect(XCJ_1_ExpTranscript_TF_list,XCJ_3_ExpTranscript_TF_list))
length(intersect(XCJ_2_ExpTranscript_TF_list,XCJ_3_ExpTranscript_TF_list))
length(XCJ_1_ExpTranscript_TF_list)
length(unique(XCJ_1_ExpTranscript_TF_list))
length(intersect(intersect(XCJ_1_ExpTranscript_TF_list,XCJ_2_ExpTranscript_TF_list),XCJ_3_ExpTranscript_TF_list))
length(XCJ_1_ExpTranscript_TF_list)
ls()
write.table(XCJ_3_ExpTranscript_TF_new,"XCJ_3_ExpTranscript_TF.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(XCJ_1_ExpTranscript_TF_new,"XCJ_1_ExpTranscript_TF.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(XCJ_2_ExpTranscript_TF_new,"XCJ_2_ExpTranscript_TF.xls",sep="\t",col.names=T,row.names=F,quote=F)
XCJ_1_Specific_TFs<-setdiff(setdiff(as.character(XCJ_1_ExpTranscript_TF_new[[1]]),as.character(XCJ_2_ExpTranscript_TF_new[[1]])),as.character(XCJ_3_ExpTranscript_TF_new[[1]])))
XCJ_1_Specific_TFs<-setdiff(setdiff(as.character(XCJ_1_ExpTranscript_TF_new[[1]]),as.character(XCJ_2_ExpTranscript_TF_new[[1]])),as.character(XCJ_3_ExpTranscript_TF_new[[1]]))
length(XCJ_1_Specific_TFs)
XCJ_1_Specific_TFs
XCJ_2_Specific_TFs<-setdiff(setdiff(as.character(XCJ_2_ExpTranscript_TF_new[[1]]),as.character(XCJ_1_ExpTranscript_TF_new[[1]])),as.character(XCJ_3_ExpTranscript_TF_new[[1]]))
XCJ_3_Specific_TFs<-setdiff(setdiff(as.character(XCJ_3_ExpTranscript_TF_new[[1]]),as.character(XCJ_1_ExpTranscript_TF_new[[1]])),as.character(XCJ_2_ExpTranscript_TF_new[[1]]))
length(XCJ_2_Specific_TFs)
length(XCJ_3_Specific_TFs)
length(XCJ_1_Specific_TFs)
XCJ_1_Specific_TFs_new<-as.data.frame()
XCJ_1_Specific_TFs_new<-as.dataframe()
XCJ_1_Specific_TFs_new<-as.data.frame(XCJ_1_Specific_TFs)
names(XCJ_1_Specific_TFs_new)<-names(XCJ_1_ExpTranscript_TF_new)[1]
names(XCJ_1_Specific_TFs_new)
names(XCJ_2_Specific_TFs_new)<-names(XCJ_2_ExpTranscript_TF_new)[1]
XCJ_2_Specific_TFs_new<-as.data.frame(XCJ_2_Specific_TFs)
XCJ_3_Specific_TFs_new<-as.data.frame(XCJ_3_Specific_TFs)
names(XCJ_1_Specific_TFs_new)
names(XCJ_2_Specific_TFs_new)<-names(XCJ_2_ExpTranscript_TF_new)[1]
names(XCJ_3_Specific_TFs_new)<-names(XCJ_3_ExpTranscript_TF_new)[1]
merge(XCJ_1_Specific_TFs_new,XCJ_1_ExpTranscript_TF_new,by="TranscriptName")->XCJ_1_Specific_TFs_new_new
merge(XCJ_2_Specific_TFs_new,XCJ_2_ExpTranscript_TF_new,by="TranscriptName")->XCJ_2_Specific_TFs_new_new
merge(XCJ_3_Specific_TFs_new,XCJ_3_ExpTranscript_TF_new,by="TranscriptName")->XCJ_3_Specific_TFs_new_new
XCJ_3_Specific_TFs_new_new[1:4,]
dim(XCJ_3_Specific_TFs_new_new)
write.table(XCJ_3_Specific_TFs_new_new,"XCJ_3_Specific_TFs.xls")
write.table(XCJ_3_Specific_TFs_new_new,"XCJ_3_Specific_TFs.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(XCJ_2_Specific_TFs_new_new,"XCJ_2_Specific_TFs.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(XCJ_1_Specific_TFs_new_new,"XCJ_1_Specific_TFs.xls",sep="\t",col.names=T,row.names=F,quote=F)
Inter_XCJ12_TFs<-intersect(as.character(XCJ_1_ExpTranscript_TF_new[[1]]),as.character(XCJ_1_ExpTranscript_TF_new[[2]]))
length(Inter_XCJ12_TFs)
as.character(XCJ_1_ExpTranscript_TF_new[[1]])
as.character(XCJ_2_ExpTranscript_TF_new[[1]])
Inter_XCJ12_TFs<-intersect(as.character(XCJ_1_ExpTranscript_TF_new[[1]]),as.character(XCJ_1_ExpTranscript_TF_new[[2]]))
Inter_XCJ12_TFs
dim(XCJ_1_ExpTranscript_TF_new)
dim(XCJ_2_ExpTranscript_TF_new)
intersect(as.character(XCJ_1_ExpTranscript_TF_new[[1]]),as.character(XCJ_1_ExpTranscript_TF_new[[2]]))
intersect(as.character(XCJ_1_ExpTranscript_TF_new[[1]]),as.character(XCJ_3_ExpTranscript_TF_new[[2]]))
intersect(as.character(XCJ_2_ExpTranscript_TF_new[[1]]),as.character(XCJ_3_ExpTranscript_TF_new[[2]]))
XCJ_2_ExpTranscript_TF_new[1:4,]
class(XCJ_2_ExpTranscript_TF_new)
intersect(as.character(XCJ_2_ExpTranscript_TF_new[[1]]),as.character(XCJ_3_ExpTranscript_TF_new[[2]]))
XCJ1<-as.character(XCJ_1_ExpTranscript_TF[,1])
XCJ1
XCJ2<-as.character(XCJ_2_ExpTranscript_TF[,1])
intersect(XCJ1,XCJ2)
dim(XCJ_2_ExpTranscript_TF_new)
dim(XCJ_2_ExpTranscript_TF)
intersect(as.character(XCJ_1_ExpTranscript_TF_new[,1]),as.character(XCJ_2_ExpTranscript_TF_new[,1]))
Intersect_TFs_XCJ12<-intersect(as.character(XCJ_1_ExpTranscript_TF_new[,1]),as.character(XCJ_2_ExpTranscript_TF_new[,1]))
Intersect_TFs_XCJ13<-intersect(as.character(XCJ_1_ExpTranscript_TF_new[,1]),as.character(XCJ_3_ExpTranscript_TF_new[,1]))
Intersect_TFs_XCJ23<-intersect(as.character(XCJ_2_ExpTranscript_TF_new[,1]),as.character(XCJ_3_ExpTranscript_TF_new[,1]))
Intersect_TFs_XCJ123<-intersect(Intersect_TFs_XCJ12,as.character(XCJ_3_ExpTranscript_TF_new[,1]))
length(Intersect_TFs_XCJ123)
length(Intersect_TFs_XCJ12)
length(Intersect_TFs_XCJ13)
length(Intersect_TFs_XCJ23)
Intersect_TFs_XCJ12_new<-as.data.frame(Intersect_TFs_XCJ12)
Intersect_TFs_XCJ13_new<-as.data.frame(Intersect_TFs_XCJ13)
Intersect_TFs_XCJ23_new<-as.data.frame(Intersect_TFs_XCJ23)
Intersect_TFs_XCJ123_new<-as.data.frame(Intersect_TFs_XCJ123)
names(Intersect_TFs_XCJ12_new)<-names(XCJ_3_ExpTranscript_TF_new)[1]
names(Intersect_TFs_XCJ13_new)<-names(XCJ_3_ExpTranscript_TF_new)[1]
names(Intersect_TFs_XCJ23_new)<-names(XCJ_3_ExpTranscript_TF_new)[1]
names(Intersect_TFs_XCJ123_new)<-names(XCJ_3_ExpTranscript_TF_new)[1]
naems(vsAtTF_Annotation)
names(vsAtTF_Annotation)
merge(Intersect_TFs_XCJ12_new,vsAtTF_Annotation,by="TranscriptName")->Intersect_TFs_XCJ12_new_new
merge(Intersect_TFs_XCJ13_new,vsAtTF_Annotation,by="TranscriptName")->Intersect_TFs_XCJ13_new_new
merge(Intersect_TFs_XCJ23_new,vsAtTF_Annotation,by="TranscriptName")->Intersect_TFs_XCJ23_new_new
merge(Intersect_TFs_XCJ123_new,vsAtTF_Annotation,by="TranscriptName")->Intersect_TFs_XCJ123_new_new
names(Intersect_TFs_XCJ123_new_new)
ls()
XCJ_1_ExpTranscript[1:4,]
rbind(rbind(XCJ_1_ExpTranscript,XCJ_2_ExpTranscript),XCJ_3_ExpTranscript)->XCJ_123_ExpTranscript
unique(rbind(rbind(XCJ_1_ExpTranscript,XCJ_2_ExpTranscript),XCJ_3_ExpTranscript))->XCJ_123_ExpTranscript
dim(XCJ_123_ExpTranscript)
history()
merge(Intersect_TFs_XCJ12_new,XCJ_123_ExpTranscript,by="TranscriptName")->Intersect_TFs_XCJ12_new_new
Intersect_TFs_XCJ12_new_new[1:4,]
XCJ_1_ExpTranscript_TF_new[1:4,]
unique(rbind(rbind(XCJ_1_ExpTranscript_TF_new,XCJ_2_ExpTranscript_TF_new),XCJ_3_ExpTranscript_TF_new))->XCJ_123_ExpTranscript_TF_new
merge(Intersect_TFs_XCJ12_new,XCJ_123_ExpTranscript_TF_new,by="TranscriptName")->Intersect_TFs_XCJ12_new_new
Intersect_TFs_XCJ12_new_new[1:4,]
merge(Intersect_TFs_XCJ13_new,XCJ_123_ExpTranscript_TF_new,by="TranscriptName")->Intersect_TFs_XCJ13_new_new
merge(Intersect_TFs_XCJ23_new,XCJ_123_ExpTranscript_TF_new,by="TranscriptName")->Intersect_TFs_XCJ23_new_new
merge(Intersect_TFs_XCJ123_new,XCJ_123_ExpTranscript_TF_new,by="TranscriptName")->Intersect_TFs_XCJ123_new_new
Intersect_TFs_XCJ123_new_new[1:4,]
write.table(Intersect_TFs_XCJ12_new_new,"XCJ_12_Intersect_TFs.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Intersect_TFs_XCJ13_new_new,"XCJ_13_Intersect_TFs.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Intersect_TFs_XCJ23_new_new,"XCJ_23_Intersect_TFs.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Intersect_TFs_XCJ123_new_new,"XCJ_123_Intersect_TFs.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("genes_exp_find.xls")[,c(1,2)]   ## 转录本水平
colnames(data)=c(" PA"," PB")
data = data[,1:2] # remove the gene column since its now the rowname value
data = as.matrix(data) # convert to matrix
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted
## write.table(centered_data,"centered_data_S3S4.txt",sep="\t")
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file=" PA_PB.heatmap_gene.pdf", width=8,height=32, paper="special");
heatmap.2(centered_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2, margins=c(8,25), lhei=c(0.3,2), lwid=c(2.5,4))
dev.off()
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("genes_exp_find.xls")[,c(1,2)]
colnames(data)=c(" PA"," PB")
data = data[,1:2]
data = as.matrix(data)
data = log2(data+1)
centered_data = t(scale(t(data), scale=F))
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean")
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete")
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6)
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file=" PA_PB.heatmap_gene.pdf", width=8,height=32, paper="special")
heatmap.2(centered_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2, margins=c(8,25), lhei=c(0.3,2), lwid=c(2.5,4))
dev.off()
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("isoforms_exp_find.xls")[,c(1,2)]   ## 转录本水平
colnames(data)=c(" PA"," PB")
data = data[,1:2] # remove the gene column since its now the rowname value
data = as.matrix(data) # convert to matrix
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted
## write.table(centered_data,"centered_data_S3S4.txt",sep="\t")
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file=" PA_PB.heatmap_isoform.pdf", width=8,height=32, paper="special");
heatmap.2(centered_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2, margins=c(8,25), lhei=c(0.3,2), lwid=c(2.5,4))
dev.off()
q()
dd<-read.delim("DGGW_vs_GW_go_enrichment_go2gene.xls",sep="\t",header=T,check.names=F)
dim(dd)
dd[1:4,]
q()
dd<-read.delim("dd",sep="\t",header=T,check.names=F)
dd[1:4,]
names(dd)
for(i in 1:nrow(dd)){
goid<-dd[i,1]
doname<-dd[i,2]
diff<-unique(unlist(strsplit(as.character(dd[i,3]),";",fix=T)))
all<-unique(unlist(strsplit(as.character(dd[i,4]),";",fix=T)))
write.table(diff,paste(goid,"diffgene.txt",sep=""),sep="\n",row.names=F,col.names=F,quote=F)
write.table(all,paste(goid,"allgene.txt",sep=""),sep="\n",row.names=F,col.names=F,quote=F)
}
q()
HF_GeneLen<-read.delim("hf.mmu.gene.expr",header=T,sep="\t",check.names=F)
CB_GeneLen<-read.delim("cb.mmu.gene.expr",header=T,sep="\t",check.names=F)
names(HF_GeneLen)[1]<-"gene_id"
names(CB_GeneLen)[1]<-"gene_id"
EdgeR<-read.delim("/mnt/lustre/users/wangyan/project/liuzhihong/MJ20120929954_2ciwei/Expression_20130822/Ref_MMU/mmu/mmu.gene.mat.CB_vs_HF.edgeR.DE_results",header=F,sep="\t")
EdgeR_Result<-EdgeR[-1,]
names(EdgeR_Result)<-c("gene_id",as.character(as.matrix(EdgeR[1,1:4])))
EdgeR_Result[1:4,]
HF_EdgeR<-merge(EdgeR_Result,HF_GeneLen,by="gene_id",all.y=T)
HF_CB_EdgeR<-merge(HF_EdgeR,CB_GeneLen,by="gene_id",all.y=T)
HF_CB_EdgeR_new<-HF_CB_EdgeR[,c("gene_id","effectiveLength.x","expectCounts.x","expectCounts.y","logFC","PValue","FDR","assembleTtranscrit.x")]
names(HF_CB_EdgeR_new)<-c("gene_id","assembleTtranscrit","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR")
HF_CB_EdgeR_new[1:4,]
names(HF_CB_EdgeR_new)<-c("gene_id","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR","assembleTtranscrit")
HF_CB_EdgeR_new[1:4,]
write.table(HF_CB_EdgeR_new,"HF_CB_EdgeR_new",sep="\t",col.names=T,row.names=F,quote=F)
q()
dir()
Mmu_Gene_Pathway_GO<-read.delim("Mmu_Gene_Pathway_GO",sep="\t",header=T,check.names=F)
Mmu_Gene_Pathway_GO[1:4,]
HF_CB_EdgeR_new_new_sig<-read.delim("HF_CB_EdgeR_new_new_sig",sep="\t",header=T,check.names=F)
HF_CB_EdgeR_new_new_sig[1:4,]
names(HF_CB_EdgeR_new_new_sig)
names(HF_CB_EdgeR_new_new_sig)[1]
names(HF_CB_EdgeR_new_new_sig)[1]<-"ENSG"
merge(HF_CB_EdgeR_new_new_sig,Mmu_Gene_Pathway_GO,by="ENSG")0>HF_CB_EdgeR_new_new_sig_Annot
merge(HF_CB_EdgeR_new_new_sig,Mmu_Gene_Pathway_GO,by="ENSG")->HF_CB_EdgeR_new_new_sig_Annot
HF_CB_EdgeR_new_new_sig_Annot[1:4,]
names(HF_CB_EdgeR_new_new_sig_Annot)
HF_CB_EdgeR_new_new_sig_dir<-read.delim("HF_CB_EdgeR_new_new_sig_dir",sep="\t",header=T,check.names=F)
names(HF_CB_EdgeR_new_new_sig)[1]
names(HF_CB_EdgeR_new_new_sig_dir)[1]
names(HF_CB_EdgeR_new_new_sig_dir)[1]<-"ENSG"
merge(HF_CB_EdgeR_new_new_sig,Mmu_Gene_Pathway_GO,by="ENSG")->HF_CB_EdgeR_new_new_sig_Annot
merge(HF_CB_EdgeR_new_new_sig_dir,Mmu_Gene_Pathway_GO,by="ENSG")->HF_CB_EdgeR_new_new_sig_Annot
names(HF_CB_EdgeR_new_new_sig_dir)
names(HF_CB_EdgeR_new_new_sig_Annot)
HF_vs_CB<-HF_CB_EdgeR_new_new_sig_Annot[,c("ENSG","effectiveLength"    "HF_expectCounts"   )]
 [4] "CB_expectCounts"    "logFC"              "PValue"            
 [7] "FDR"                "Sig(FDR<=0.05)"     "Dir"               
HF_vs_CB<-HF_CB_EdgeR_new_new_sig_Annot[,c("ENSG","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR","Sig(FDR<=0.05)","Dir","assembleTtranscrit")]
HF_vs_CB[1:4,]
HF_vs_CB[1:10,]
HF_vs_CB<-unique(HF_CB_EdgeR_new_new_sig_Annot[,c("ENSG","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR","Sig(FDR<=0.05)","Dir","assembleTtranscrit")])
HF_vs_CB[1:10,]
write.table(HF_vs_CB,"HF_vs_CB",sep="\t",col.names=T,row.names=F,quote=F)
q()
HF_vs_CB[1:10,]
HF_vs_CB<-unique(HF_CB_EdgeR_new_new_sig_Annot[,c("ENSG","Name","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR","Sig(FDR<=0.05)","Dir","assembleTtranscrit")])
HF_vs_CB[1:10,]
write.table(HF_vs_CB,"HF_vs_CB",sep="\t",col.names=T,row.names=F,quote=F)
q()
dir()
Homo_1<-read.delim("HomoSpecific_ENSG2Num2Tlist",sep="\t",header=F)
Homo_1[1:4,]
names(Homo_1)
names(Homo_1)<-c("EnsembleGeneID","TranscriptNum","TranscriptList")
HomoSpecificTranscript<-read.delim("2_TranscriptsOnlyMapped2HomoGene.txt",header=F)
HomoSpecificTranscript[1:10,]
HomoMmuTranscript<-read.delim("4_TranscriptsNotOnlyMapped2HomoGeneAndMmuGene.txt",header=F)
HomoMmuTranscript[1:10,]
Homo_2<-data.frame()
for(i in 1:nrow(Homo_1)){
history()
Homo_1<-read.delim("HomoSpecific_ENSG2Num2Tlist",sep="\t",header=F)
names(Homo_1)<-c("EnsembleGeneID","TranscriptNum","TranscriptList")
HomoSpecificTranscript<-as.character(read.delim("2_TranscriptsOnlyMapped2HomoGene.txt",header=F)[[1]])
HomoMmuTranscript<-as.character(read.delim("4_TranscriptsNotOnlyMapped2HomoGeneAndMmuGene.txt",header=F)[[1]])
Homo_2<-data.frame()
Homo_1[i,3]
Homo_1[1,3]
Homo_1[2,3]
Homo_1[3,3]
Homo_1[1:10,]
getwd()
Homo_1[1:20,]
Homo_1[1:50,]
which(Homo_1[[1]]=="ENSG00000006059")
which(as.character(Homo_1[[1]])=="ENSG00000006059")
which(as.character(Homo_1[[1]])=="ENSG00000079841")
Homo_1[165,3]
as.character(Homo_1[165,3])%in%HomoSpecificTranscript
unlist(strsplit(as.character(Homo_1[165,3]),",",fix=T))
TranscriptLIST<-unlist(strsplit(as.character(Homo_1[165,3]),",",fix=T))
    TranscriptLIST%in%HomoSpecificTranscript
TranscriptLIST%in%HomoMmuTranscript
TranscriptLIST<-unlist(strsplit(as.character(Homo_1[165,3]),",",fix=T))
    TranscriptLIST[TranscriptLIST%in%HomoSpecificTranscript]
paste(TranscriptLIST[TranscriptLIST%in%HomoSpecificTranscript],collapse=";")
TranscriptLIST_1[TranscriptLIST_1%in%HomoMmuTranscript]
i
Homo_2<-data.frame()
for(i in 1:nrow(Homo_1)){
    Homo_2[i,1]<-Homo_1[i,1]
    Homo_2[i,2]<-Homo_1[i,2]
    Homo_2[i,3]<-Homo_1[i,3]
    TranscriptLIST_1<-unlist(strsplit(as.character(Homo_1[i,3]),",",fix=T))
    Homo_2[i,4]<-paste(TranscriptLIST_1[TranscriptLIST_1%in%HomoSpecificTranscript],collapse=",")
    Homo_2[i,4]<-paste(TranscriptLIST_1[TranscriptLIST_1%in%HomoMmuTranscript],collapse=",")
    print(i)
}
ls()
Homo_2[1:4,]
Homo_2<-data.frame()
for(i in 1:nrow(Homo_1)){
    Homo_2[i,1]<-Homo_1[i,1]
    Homo_2[i,2]<-Homo_1[i,2]
    Homo_2[i,3]<-Homo_1[i,3]
    TranscriptLIST_1<-unlist(strsplit(as.character(Homo_1[i,3]),",",fix=T))
    Homo_2[i,4]<-paste(TranscriptLIST_1[TranscriptLIST_1%in%HomoSpecificTranscript],collapse=",")
    Homo_2[i,5]<-paste(TranscriptLIST_1[TranscriptLIST_1%in%HomoMmuTranscript],collapse=",")
    print(i)
}
Homo_2[1:50,]
Homo_2[1:20,]
Homo_2[,5]
Homo_2[,3:4]
Homo_2[1:20,3:4]
ls()
q()
ls()
HomoMmuTranscript[1:10]
Homo_1<-read.delim("HomoSpecific_ENSG2Num2Tlist",sep="\t",header=F)
names(Homo_1)<-c("EnsembleGeneID","TranscriptNum","TranscriptList")
HomoSpecificTranscript<-as.character(read.delim("2_TranscriptsOnlyMapped2HomoGene.txt",header=F)[[1]])
HomoMmuTranscript<-as.character(read.delim("4_TranscriptsNotOnlyMapped2HomoGeneAndMmuGene.txt",header=F)[[1]])
Homo_2<-data.frame()
q()
Homo_1<-read.delim("Homo_GeneID2TranscriptNum2List",sep="\t",header=F)
names(Homo_1)<-c("EnsembleGeneID","TranscriptNum","TranscriptList")
HomoSpecificTranscript<-as.character(read.delim("2_TranscriptsOnlyMapped2HomoGene.txt",header=F)[[1]])
HomoMmuTranscript<-as.character(read.delim("4_TranscriptsNotOnlyMapped2HomoGeneAndMmuGene.txt",header=F)[[1]])
Homo_2<-data.frame()
for(i in 1:nrow(Homo_1)){
    Homo_2[i,1]<-Homo_1[i,1]
    Homo_2[i,2]<-Homo_1[i,2]
    Homo_2[i,3]<-Homo_1[i,3]
    TranscriptLIST_1<-unlist(strsplit(as.character(Homo_1[i,3]),",",fix=T))
    Homo_2[i,4]<-paste(TranscriptLIST_1[TranscriptLIST_1%in%HomoSpecificTranscript],collapse=",")
    Homo_2[i,5]<-paste(TranscriptLIST_1[TranscriptLIST_1%in%HomoMmuTranscript],collapse=",")
    print(i)
}
Homo_2[1:4,]
Homo_2[1:40,]
Homo_2[1:40,3:5]
Homo_2[1:50,3:5]
write.table(Homo_2,"Homo_2",sep="\t",col.names=F,row.names=F,quote=F)
q()
Homo_2[1:4,]
Homo_2[1:40,]
Homo_2[1:4,]
Homo_3<-Homo_2[-1,]
Homo_3[1:4,]
names(Homo_3)
names(Homo_3)<-c("EnsembleGeneID","TranscriptNum","TranscriptList","SpecificMapped2Homo","Mapped2HomoandMmu")
write.table(Homo_3,"Homo_2",sep="\t",col.names=F,row.names=F,quote=F)
q()
write.table(Homo_3,"Homo_2",sep="\t",col.names=T,row.names=F,quote=F)
q()
history()
q()
Mmu_3[1:4,]
ls()
q()
dir()
All<-as.character(read.delim("AllTranscripts",header=F)[[1]])
All[1:4]
Homo<-as.character(read.delim("Mapped2HomoTranscripts",header=F)[[1]])
Homo[1:4,]
Homo[1:4]
Homo<-paste(as.character(read.delim("Mapped2HomoTranscripts",header=F)[[1]]),sep="")
Homo<-unlist(strsplit(paste(as.character(read.delim("Mapped2HomoTranscripts",header=F)[[1]]),sep=""),",",fix=T))
Homo<-unique(unlist(strsplit(paste(as.character(read.delim("Mapped2HomoTranscripts",header=F)[[1]]),sep=""),",",fix=T)))
Homo[1:40]
Mmu<-unique(unlist(strsplit(paste(as.character(read.delim("Mapped2MmuTranscripts",header=F)[[1]]),sep=""),",",fix=T)))
Mmu[1:40]
ls()
HomoMmu<-intersect(Homo,Mmu)
length(HomoMmu)
HomoSpecific<-setdiff(Homo,HomoMmu)
intersect(HomoSpecific,HomoMmu)
HomoSpecific<-setdiff(Mmu,HomoMmu)
MmuSpecific<-setdiff(Mmu,HomoMmu)
HomoSpecific<-setdiff(Homo,HomoMmu)
ls()
length(MmuSpecific)
length(HomoSpecific)
None<-setdiff(All,unique(c(Homo,Mmu)))
length(None)
length(All)
ls()
write.table(HomoMmu,"Homo_Mmu.xls",sep="\n",col.names=F,row.names=F,quote=F)
write.table(HomoSpecific,"HomoSpecific.xls",sep="\n",col.names=F,row.names=F,quote=F)
write.table(MmuSpecific,"MmuSpecific.xls",sep="\n",col.names=F,row.names=F,quote=F)
write.table(None,"noMapped.xls",sep="\n",col.names=F,row.names=F,quote=F)
q()
dir()
Homo_Exp<-read.delim("Homo_exp",header=T,sep="\t")
Mmu_Exp<-read.delim("Mmu_exp",header=T,sep="\t")
Mmu_Exp<-read.delim("Mmu_exp",header=T,sep="\t")
Mmu_Exp<-read.delim("Mmu_exp",header=T,sep="\t",check.names=F)
ls()
Mmu_Exp<-read.delim("Mmu_exp",header=T,sep="\t")
Mmu_Exp<-read.delim("Mmu_exp",sep="\t")
q()
ls()
Mmu_Exp<-read.delim("Mmu_exp_new",sep="\t")
Mmu_Exp[1:4,]
Mmu_Exp<-read.delim("Mmu_exp_new",sep="\t",header=F)
Homo_exp[1:4,]
Homo_Exp[1:4,]
Homo_exp<-Homo_Exp[,c(1,11)]
Homo_exp[1:4,]
dim(Homo_exp)
HomoSpecificGenes<-vector()
Homo_Transcripts<-as.character(Homo_exp[[2]])
Homo_Transcripts[1:10]
Homo_Transcripts[1]
unlist(strsplit(Homo_Transcripts[1],",",fix=T)
)
unlist(strsplit(Homo_Transcripts[1],",",fix=T))
unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific
which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)=="TRUE")
which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)==TRUE)
which((unlist(strsplit(Homo_Transcripts[2],",",fix=T))%in%HomoSpecific)==TRUE)
which((unlist(strsplit(Homo_Transcripts[2],",",fix=T))%in%HomoSpecific)==FALSE)
which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)==FALSE)
length(which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)==FALSE))
length(which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)==TRUE))
length(which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)==TRUE))==length(unlist(strsplit(Homo_Transcripts[1],",",fix=T)))
length(which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)==FALSE))==length(unlist(strsplit(Homo_Transcripts[1],",",fix=T)))
length(which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)==TRUE))==length(unlist(strsplit(Homo_Transcripts[1],",",fix=T)))
sapply(Homo_Transcripts,function(x) if(length(which((unlist(strsplit(Homo_Transcripts[1],",",fix=T))%in%HomoSpecific)==TRUE))==length(unlist(strsplit(Homo_Transc
ripts[1],",",fix=T)))) y=1;else y=0)
ls()
HomoSpecificGenes
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(Homo_Transcripts[1],",",fix=T))
Flag<-1
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(Homo_Transcripts[1],",",fix=T))
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(Homo_Transcripts[i],",",fix=T))
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(Homo_exp[i,11],",",fix=T))
if(length(which(TranscriptList%in%HomoSpecific)==TRUE))==length(TranscriptList)){HomoSpecificGenes[Flag]<-Homo_exp[i,1]}
HomoSpecificGenes<-vector()
Flag
Homo_exp[1:4,]
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(Homo_exp[i,2],",",fix=T))
if(length(which((TranscriptList%in%HomoSpecific)==TRUE))==length(TranscriptList)){HomoSpecificGenes[Flag]<-Homo_exp[i,1]}
Flag<-Flag+1
print(i)}
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(as.character(Homo_exp[i,2]),",",fix=T))
if(length(which((TranscriptList%in%HomoSpecific)==TRUE))==length(TranscriptList)){HomoSpecificGenes[Flag]<-Homo_exp[i,1]}
Flag<-Flag+1
print(i)}
ls()
HomoSpecificGenes[1:4]
HomoSpecificGenes<-vector()
Flag<-1
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(as.character(Homo_exp[i,2]),",",fix=T))
if(length(which((TranscriptList%in%HomoSpecific)==TRUE))==length(TranscriptList)){HomoSpecificGenes[Flag]<-as.character(Homo_exp[i,1])}
Flag<-Flag+1
}
HomoSpecificGenes[1:10]
HomoSpecificGenes[1:100]
HomoSpecificGenes<-vector()
Flag<-1
for(i in 1:nrow(Homo_exp)){
if(length(which((TranscriptList%in%HomoSpecific)==TRUE))==length(TranscriptList)){HomoSpecificGenes[Flag]<-as.character(Homo_exp[i,1])
Flag<-Flag+1
}
print(i)
}
ls()
HomoSpecificGenes[1:10]
HomoSpecificGenes[1:100]
HomoSpecificGenes[1:1000]
HomoSpecificGenes<-vector()
Flag<-1
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(as.character(Homo_exp[i,2]),",",fix=T))
if(length(which((TranscriptList%in%HomoSpecific)==TRUE))==length(TranscriptList)){HomoSpecificGenes[Flag]<-as.character(Homo_exp[i,1]);Flag<-Flag+1}
print(i)
}
ls()
HomoSpecificGenes[1:10]
length(HomoSpecificGenes)
HomoSpecificGenes
write.table(HomoSpecificGenes,"HomoSpecificGenes.xls",sep="\n",col.names=F,row.names=F,quote=F)
history(100)
ls()
MmuSpecificGenes<-vector()
Mmu_exp[1:4,]
Mmu_Exp[1:4,]
Flag<-1
for(i in 1:nrow(Mmu_Exp)){
TranscriptList<-unlist(strsplit(as.character(Mmu_Exp[i,2]),",",fix=T))
if(length(which((TranscriptList%in%MmuSpecific)==TRUE))==length(TranscriptList)){MmuSpecificGenes[Flag]<-as.character(Mmu_Exp[i,1]);Flag<-Flag+1}
print(i)
}
ls()
MmuSpecificGenes[1:10]
ls()
length(MmuSpecificGenes)
write.table(MmuSpecificGenes,"MmuSpecificGenes.xls",sep="\n",col.names=F,row.names=F,quote=F)
q()
ls()
HomoMmu[1:4]
history()
history(1000)
Homo_exp[1:4,]
history(1000)
ls()
length(HomoMmu)
MmuandHomo_MmuGenes<-vector()
Flag<-1
for(i in 1:nrow(Mmu_Exp)){
TranscriptList<-unlist(strsplit(as.character(Mmu_Exp[i,2]),",",fix=T))
if(length(which((TranscriptList%in%HomoMmu)==TRUE))==length(TranscriptList)){MmuandHomo_MmuGenes[Flag]<-as.character(Mmu_Exp[i,1]);Flag<-Flag+1}
print(i)
}
MmuandHomo_MmuGenes[1:4]
ls()
MmuandHomo_HomoGenes<-vector()
Flag<-1
for(i in 1:nrow(Homo_Exp)){
TranscriptList<-unlist(strsplit(as.character(Homo_Exp[i,2]),",",fix=T))
if(length(which((TranscriptList%in%HomoMmu)==TRUE))==length(TranscriptList)){MmuandHomo_HomoGenes[Flag]<-as.character(Homo_Exp[i,1]);Flag<-Flag+1}
print(i)
}
hitory()
history()
history(200)
write.table(MmuandHomo_MmuGenes,"MmuandHomo_MmuGenes.xls",sep="\n",col.names=F,row.names=F,quote=F)
write.table(MmuandHomo_HomoGenes,"MmuandHomo_HomoGenes.xls",sep="\n",col.names=F,row.names=F,quote=F)
q()
ls()
length(MmuandHomo_MmuGenes)
length(MmuandHomo_HomoGenes)
MmuandHomo_HomoGenes
Homo_Exp[1:4,]
dim(Homo_Exp)
history(200)
Homo_exp[1:4,]
MmuandHomo_HomoGenes<-vector()
Flag<-1
for(i in 1:nrow(Homo_exp)){
TranscriptList<-unlist(strsplit(as.character(Homo_exp[i,2]),",",fix=T))
if(length(which((TranscriptList%in%HomoMmu)==TRUE))==length(TranscriptList)){MmuandHomo_HomoGenes[Flag]<-as.character(Homo_exp[i,1]);Flag<-Flag+1}
print(i)
}
write.table(MmuandHomo_HomoGenes,"MmuandHomo_HomoGenes.xls",sep="\n",col.names=F,row.names=F,quote=F)
length(MmuandHomo_HomoGenes)
length(MmuandHomo_MmuGenes)
MmuandHomo_HomoGenes[1:10]
MmuandHomo_MmuGenes[1:10]
q()
write.table(MmuandHomo_HomoGenes,"MmuandHomo_HomoGenes.xls",sep="\n",col.names=F,row.names=F,quote=F)
length(MmuandHomo_HomoGenes)
q()
dir()
SDEGs<-read.dleim("HF_CB_EdgeR_SDEGs",header=F)
SDEGs<-read.delim("HF_CB_EdgeR_SDEGs",header=F)
Homo_Annot<-read.delim("Homo_Gene_Pathway_GO",header=T,sep="\t",check.names=F)
Homo_Annot[1:4,]
names(SDEGs)
names(SDEGs)
names(SDEGs)<-"ENST"
merge(SDEGs,Homo_Annot,by="ENST")->SDEG_Annot
SDEG_Annot[1:4,]
SDEGs[1:4,]
names(SDEGs)<-"ENSG"
merge(SDEGs,Homo_Annot,by="ENSG")->SDEG_Annot
SDEG_Annot[1:4,]
write.table(SDEG_Annot,"SDEG_Annot",sep="\t",col.names=T,row.names=F,quote=F)
q()
AllGene<-read.delim("Homo_AllGene",header=F)
names(AllGene)
names(AllGene)<-"ENSG"
merge(AllGene,Homo_Annot,by="ENSG")->AllGene_Annot
AllGene_Annot[1:4,]
AllGeneName_<-unique(AllGene_Annot[,c("NAME","GOID")])
dim(AllGeneName_GO)
AllGeneName_GO[1:4,]
write.table(AllGeneName_GO,"AllGeneName_GO",sep="\t",col.names=T,row.names=F,quote=F)
q()
q()
write.table(AllGeneName_GO,"AllGeneName_GO",sep="\t",col.names=T,row.names=F,quote=F)
q()
HF_GeneLen<-read.delim("hf.homo.gene.expr",header=T,sep="\t",check.names=F)
CB_GeneLen<-read.delim("cb.homo.gene.expr",header=T,sep="\t",check.names=F)
HF_GeneLen[1:4,]
getwd()
EdgeR<-read.delim("/mnt/lustre/users/wangyan/project/liuzhihong/MJ20120929954_2ciwei/Expression_20130822/Ref_Homo/homo.gene.mat.CB_vs_HF.edgeR.DE_results",header=F,sep="\t")
EdgeR<-read.delim("/mnt/lustre/users/wangyan/project/liuzhihong/MJ20120929954_2ciwei/Expression_20130822/Ref_Homo/homo.gene.mat.CB_vs_HF.edgeR.DE_results",header=F,sep="\t")
EdgeR<-read.delim("/mnt/lustre/users/wangyan/project/liuzhihong/MJ20120929954_2ciwei/Expression_20130822/Ref_Homo/EdgeR/homo.gene.mat.CB_vs_HF.edgeR.DE_results",header=F,sep="\t")
EdgeR[1:4,]
EdgeR_Result<-EdgeR[-1,]
names(EdgeR_Result)<-c("gene_id",as.character(EdgeR[1,1:4]))
EdgeR_Result[1:4,]
EdgeR[1,]
EdgeR[1,1:4]
class(EdgeR[1,1:4])
as.character(EdgeR[1,1:4])
EdgeR[1,1:4]
EdgeR[1,1:4][1,]
as.character(EdgeR[1,1:4][1,])
as.character(as.matrix(EdgeR[1,1:4]))
EdgeR<-read.delim("/mnt/lustre/users/wangyan/project/liuzhihong/MJ20120929954_2ciwei/Expression_20130822/Ref_Homo/EdgeR/homo.gene.mat.CB_vs_HF.edgeR.DE_results",header=F,sep="\t")
EdgeR_Result<-EdgeR[-1,]
names(EdgeR_Result)<-c("gene_id",as.character(as.matrix(EdgeR[1,1:4])))
EdgeR_Result[1:4,]
HF_GeneLen[1:4,]
names(HF_GeneLen)<-"gene_id"
names(CB_GeneLen)<-"gene_id"
dim(HF_GeneLen)
dim(CB_GeneLen)
HF_EdgeR<-merge(EdgeR_Result,HF_GeneLen,by="gene_id",all.y=T)
HF_CB_EdgeR<-merge(HF_EdgeR,CB_GeneLen,by="gene_id",all.y=T)
HF_EdgeR[1:4,]
HF_GeneLen[1:4,]
HF_GeneLen[1:4,]
### Homo
HF_GeneLen<-read.delim("hf.homo.gene.expr",header=T,sep="\t",check.names=F)
CB_GeneLen<-read.delim("cb.homo.gene.expr",header=T,sep="\t",check.names=F)
names(HF_GeneLen)[1]<-"gene_id"
names(CB_GeneLen)[1]<-"gene_id"
EdgeR<-read.delim("/mnt/lustre/users/wangyan/project/liuzhihong/MJ20120929954_2ciwei/Expression_20130822/Ref_Homo/EdgeR/homo.gene.mat.CB_vs_HF.edgeR.DE_results",header=F,sep="\t")
EdgeR_Result<-EdgeR[-1,]
names(EdgeR_Result)<-c("gene_id",as.character(as.matrix(EdgeR[1,1:4])))
HF_EdgeR<-merge(EdgeR_Result,HF_GeneLen,by="gene_id",all.y=T)
HF_CB_EdgeR<-merge(HF_EdgeR,CB_GeneLen,by="gene_id",all.y=T)
HF_CB_EdgeR[1:10,]
HF_CB_EdgeR[1:4,]
names(### Homo)
HF_GeneLen<-read.delim("hf.homo.gene.expr",header=T,sep="\t",check.names=F)
CB_GeneLen<-read.delim("cb.homo.gene.expr",header=T,sep="\t",check.names=F)
names(HF_GeneLen)[1]<-"gene_id"
names(CB_GeneLen)[1]<-"gene_id"
EdgeR<-read.delim("/mnt/lustre/users/wangyan/project/liuzhihong/MJ20120929954_2ciwei/Expression_20130822/Ref_Homo/EdgeR/homo.gene.mat.CB_vs_HF.edgeR.DE_results",header=F,sep="\t")
EdgeR_Result<-EdgeR[-1,]
names(EdgeR_Result)<-c("gene_id",as.character(as.matrix(EdgeR[1,1:4])))
HF_EdgeR<-merge(EdgeR_Result,HF_GeneLen,by="gene_id",all.y=T)
HF_CB_EdgeR<-merge(HF_EdgeR,CB_GeneLen,by="gene_id",all.y=T)
### Homo
HF_GeneLen<-read.delim("hf.homo.gene.expr",header=T,sep="\t",check.names=F)
CB_GeneLen<-read.delim("cb.homo.gene.expr",header=T,sep="\t",check.names=F)
names(HF_GeneLen)[1]<-"gene_id"
names(CB_GeneLen)[1]<-"gene_id"
EdgeR<-read.delim("/mnt/lustre/users/wangyan/project/liuzhihong/MJ20120929954_2ciwei/Expression_20130822/Ref_Homo/EdgeR/homo.gene.mat.CB_vs_HF.edgeR.DE_results",header=F,sep="\t")
EdgeR_Result<-EdgeR[-1,]
names(EdgeR_Result)<-c("gene_id",as.character(as.matrix(EdgeR[1,1:4])))
HF_EdgeR<-merge(EdgeR_Result,HF_GeneLen,by="gene_id",all.y=T)
HF_CB_EdgeR<-merge(HF_EdgeR,CB_GeneLen,by="gene_id",all.y=T)
HF_CB_EdgeR[1:4,]
names(HF_CB_EdgeR)
HF_CB_EdgeR[1:10,c(2,8,11),]
HF_CB_EdgeR_new<-HF_CB_EdgeR[,c("gene_id","assembleTtranscrit.x","effectiveLength.x","expectCounts.x","expectCounts.y","logFC","PValue","FDR")]
HF_CB_EdgeR_new[1:4,]
names(HF_CB_EdgeR_new)<-c("gene_id","assembleTtranscrit","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR")
HF_CB_EdgeR_new[1:4,]
is.na(HF_CB_EdgeR_new[["FDR"]])
HF_CB_EdgeR_new[is.na(HF_CB_EdgeR_new[["FDR"]]),]->NA_HF_CB_EdgeR_new
NA_HF_CB_EdgeR_new[1:4,]
write.table(HF_CB_EdgeR_new,"HF_CB_EdgeR_new",sep="\t",col.names=T,row.names=F,quote=F)
q()
HF_CB_EdgeR_new<-HF_CB_EdgeR[,c("gene_id","effectiveLength.x","expectCounts.x","expectCounts.y","logFC","PValue","FDR","assembleTtranscrit.x")]
names(HF_CB_EdgeR_new)<-c("gene_id","assembleTtranscrit","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR")
write.table(HF_CB_EdgeR_new,"HF_CB_EdgeR_new",sep="\t",col.names=T,row.names=F,quote=F)
q()
dir()
ls()
HF_CB_EdgeR_new_new_sig<-read.delim("HF_CB_EdgeR_new_new_sig",sep="\t",header=T,check.names=F)
HF_CB_EdgeR_new_new_sig[1:4,]
Homo_Annot<-read.delim("Homo_Gene_Pathway_GO",sep="\t",header=T,check.names=F)
Homo_Annot[1:4,]
merge(Homo_Annot,HF_CB_EdgeR_new_new_sig,by="ENSG")->HF_CB_EdgeR_new_new_sig_Annot
names(HF_CB_EdgeR_new_new_sig)
names(HF_CB_EdgeR_new_new_sig)[1]<-"ENSG"
merge(Homo_Annot,HF_CB_EdgeR_new_new_sig,by="ENSG")->HF_CB_EdgeR_new_new_sig_Annot
HF_CB_EdgeR_new_new_sig_Annot[1:4,]
getwwd()
getwd()
HF_CB_EdgeR_new_new_sig[1:4,]
HF_CB_EdgeR_new_new_sig<-read.delim("HF_CB_EdgeR_new_new_sig_dir",sep="\t",header=T,check.names=F)
names(HF_CB_EdgeR_new_new_sig)[1]<-"ENSG"
merge(Homo_Annot,HF_CB_EdgeR_new_new_sig,by="ENSG")->HF_CB_EdgeR_new_new_sig_Annot
HF_CB_EdgeR_new_new_sig_Annot[1:4,]
write.table(HF_CB_EdgeR_new_new_sig_Annot,"HF_CB_EdgeR_new_new_sig_Annot",sep="\t",col.names=T,row.names=F,quote=F)
q()
ls()
HF_CB_EdgeR_new_new_sig_Annot[1:4,]
naems(HF_CB_EdgeR_new_new_sig_Annot)
HF_CB_EdgeR_new_new_sig_Annot_sort<-HF_CB_EdgeR_new_new_sig_Annot[order(HF_CB_EdgeR_new_new_sig_Annot[["Sig(FDR<=0.05)"]],decreasing=T),]
names(HF_CB_EdgeR_new_new_sig_Annot)
HF_CB_EdgeR_new_new_sig_Annot_sort[1:10,]
HF_CB_EdgeR_new_new_sig_Annot[1:4,]
HF_CB_EdgeR_new_new_sig_Annot_new<-HF_CB_EdgeR_new_new_sig_Annot[,-2]
HF_CB_EdgeR_new_new_sig_Annot_new[1:4,]
HF_CB_EdgeR_new_new_sig_Annot_new<-unique(HF_CB_EdgeR_new_new_sig_Annot[,-2])
HF_CB_EdgeR_new_new_sig_Annot_new[1:4,]
HF_CB_EdgeR_new_new_sig_Annot_new_sort<-HF_CB_EdgeR_new_new_sig_Annot_new[order(HF_CB_EdgeR_new_new_sig_Annot_new,decreasing=T),]
HF_CB_EdgeR_new_new_sig_Annot_new_sort[1:4,]
HF_CB_EdgeR_new_new_sig_Annot_new[1:4,]
HF_CB_EdgeR_new_new_sig_Annot_new_sort<-HF_CB_EdgeR_new_new_sig_Annot_new[order(HF_CB_EdgeR_new_new_sig_Annot_new[["Sig(FDR<=0.05)"]],decreasing=T),]
HF_CB_EdgeR_new_new_sig_Annot_new_sort[1:4,]
HF_CB_EdgeR_new_new_sig_Annot_new_sort[1:10,]

HF_CB_EdgeR_new_new_sig_Annot_new_sort[1:4,]
HF_vs_CB_diffexp_result<-HF_CB_EdgeR_new_new_sig_Annot_new_sort[,c("ENSG","NAME","effectiveLength","HF_expectCounts","CB_expectCounts",)]
HF_vs_CB_diffexp_result<-HF_CB_EdgeR_new_new_sig_Annot_new_sort[,c("ENSG","NAME","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR","Sig(FDR<=0.05)","Dir","assembleTtranscrits")]
HF_vs_CB_diffexp_result[1:4,]
HF_vs_CB_diffexp_result<-unique(HF_CB_EdgeR_new_new_sig_Annot_new_sort[,c("ENSG","NAME","effectiveLength","HF_expectCounts","CB_expectCounts","logFC","PValue","FDR","Sig(FDR<=0.05)","Dir","assembleTtranscrits")])
HF_vs_CB_diffexp_result[1:4,]
write.table(HF_vs_CB_diffexp_result,"HF_vs_CB_diffexp_result",sep="\t",col.names=T,row.names=F,quote=F)
q()
DIR()
dir()
GO<-read.delim("Homo_Mmu_Ts_Gene_GO_uniq",header=T,sep="\t",check.names=F)
GO[1:4,]
GO<-read.delim("Homo_Mmu_Ts_Gene_GO_uniq",header=F,sep="\t",check.names=F)
names(GO)<-c("")
GO[1:40,]
ls()
GO[1:4,]
Genes<-as.character(unique(GO[[1]]))
Genes
Genes<-unique(as.character(GO[[1]]))
Genes
length(Genes)
ls()
GO[1:4,]
gene_golist<-data.frame()
for(i in 1:length(Genes)){
gene_golist[i,1]<-Genes[i]
GOs<-unique(as.character(GO[which(as.character(GO[[1]])==Genes[i]),3]))
gene_golist[i,2]<-paste(GOs,collapse=";")
}
gene_golist[1:4,]
names(gene_golist)
names(gene_golist)<-c("geneID","GOID")
write.table(gene_golist,"geneID_GOlist",sep="\t",col.names=T,row.names=F,quote=F)
getwd()
q()
dir()
ENSTBEP_paths<-read.delim("ENSTBEP_paths",header=F,sep="\t")
ENSTBEP_paths[1:4,]
names(ENSTBEP_paths)
names(ENSTBEP_paths)<-c("ENSTBEP","pathID")
ENSTBEP_paths[1:40,]
dim(ENSTBEP_paths)
ls()
ENSTBEP_paths_new<-data.frame()
flag<-1
for(i in 1:nrow(ENSTBEP_paths)){
ENSTBEP_paths_new[flag,1]<-as.character(ENSTBEP_paths[i,1])
pathlist<-unlist(strsplit(as.character(ENSTBEP_paths[i,2]),";",fix=T))
N<-length(pathlist)
ENSTBEP_paths_new<-data.frame()
 flag<-1
 for(i in 1:nrow(ENSTBEP_paths)){
pathlist<-unlist(strsplit(as.character(ENSTBEP_paths[i,2]),";",fix=T))
N<-length(pathlist)
ENSTBEP_paths_new[flag:(flag+N-1),1]<-rep(as.character(ENSTBEP_paths[i,1]),N)
ENSTBEP_paths_new[flag:(flag+N-1),2]<-pathlist
flag<-flag+N
}
ENSTBEP_paths_new
dim(ENSTBEP_paths_new)
ENSTBEP_paths_new[1:4,]
ENSTBEP_paths_new[1:40,]
names(ENSTBEP_paths_new)
names(ENSTBEP_paths_new)<-c("ENSTBEP","pathID")
ENSTBEP_paths_new[1:40,]
ENSTBEP_paths_new[1:4,]
ls()
dir()
q()
ls()
dir()
pathID_pathName<-read.table("pathID_pathName",header=T,sep="\t")
pathID_pathName[1:4,]
names(pathID_pathName)
names(pathID_pathName)<-c("pathID","pathName")
names(pathID_pathName)
names(ENSTBEP_paths_new)
merge(ENSTBEP_paths_new,pathID_pathName,by="pathID")->ENSTBEP_pathID_pathName
ENSTBEP_pathID_pathName[1:4,]
ENSTBEP_paths_new[1:4,]
write.table(ENSTBEP_pathID_pathName[,c(2,1,3)],"ENSTBEP_pathID_pathName",sep="\t",col.names=T,row.names=F,quote=F)
dir()
q()
geneID_GO<-read.delim("geneID_GO",header=T,sep="\t")
dim(geneID_GO)
geneID_GO[1:4,]
q()
ls()
geneID_GO[1:4,]
geneID_GO<-read.delim("geneID_GO",header=T,sep="\t",check.names=F)
geneID_GO[1:4,]
names(geneID_GO)
names(ENSTBEP_pathID_pathName)
names(ENSTBEP_pathID_pathName)
names(ENSTBEP_pathID_pathName)[2]<-"Ensembl Protein ID"
names(ENSTBEP_pathID_pathName)
merge(geneID_GO,ENSTBEP_pathID_pathName,by="Ensembl Protein ID")->geneIDs_GO_path
geneIDs_GO_path[1:4,]
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName<-geneIDs_GO_path[,c(3,2,1,4,5,6,7)]
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName[1:4,]
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName<-geneIDs_GO_path[,c(3,2,1,4,5,6,7,8)]
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName[1:4,]
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName[1:40,]
write.table(TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName,"TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName",sep="\t",col.names=T,row.names=F,quote=F)
q()
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName[1:4,]
geneIDs_GO_path[1:4,]
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName<-geneIDs_GO_path[,c(3,2,4,5,6,7,8)]
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName[1:4,]
write.table(TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName,"TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName",sep="\t",col.names=T,row.names=F,quote=F)
q()
history(200)
quote=c("\n","\"'"),fill=T
pathID_pathName<-read.table("pathID_pathName",header=T,sep="\t",quote=c("\n","\"'"),fill=T)
history(200)
history(200)
names(ENSTBEP_paths)<-c("ENSTBEP","pathID")
ENSTBEP_paths_new<-data.frame()
 flag<-1
 for(i in 1:nrow(ENSTBEP_paths)){
pathlist<-unlist(strsplit(as.character(ENSTBEP_paths[i,2]),";",fix=T))
N<-length(pathlist)
ENSTBEP_paths_new[flag:(flag+N-1),1]<-rep(as.character(ENSTBEP_paths[i,1]),N)
ENSTBEP_paths_new[flag:(flag+N-1),2]<-pathlist
flag<-flag+N
}
history(200)
history(500)
ls()
history(500)
history(600)
history(1000)
history(1000)
history(1000)->ll
history(1000)
history(1000)
history(1000)
ENSTBEP_paths<-read.delim("ENSTBEP_paths",header=F,sep="\t")
ENSTBEP_paths_new<-data.frame()
flag<-1
 for(i in 1:nrow(ENSTBEP_paths)){
pathlist<-unlist(strsplit(as.character(ENSTBEP_paths[i,2]),";",fix=T))
N<-length(pathlist)
ENSTBEP_paths_new[flag:(flag+N-1),1]<-rep(as.character(ENSTBEP_paths[i,1]),N)
ENSTBEP_paths_new[flag:(flag+N-1),2]<-pathlist
flag<-flag+N
}
names(ENSTBEP_paths_new)<-c("ENSTBEP","pathID")
##
pathID_pathName<-read.table("pathID_pathName",header=T,sep="\t")
names(pathID_pathName)<-c("pathID","pathName")
merge(ENSTBEP_paths_new,pathID_pathName,by="pathID")->ENSTBEP_pathID_pathName
write.table(ENSTBEP_pathID_pathName[,c(2,1,3)],"ENSTBEP_pathID_pathName",sep="\t",col.names=T,row.names=F,quote=F) 
q()
pathID_pathName<-read.table("pathID_pathName",header=T,sep="\t",quote=c("\n","\"'"),fill=T)
names(pathID_pathName)<-c("pathID","pathName")
merge(ENSTBEP_paths_new,pathID_pathName,by="pathID")->ENSTBEP_pathID_pathName
write.table(ENSTBEP_pathID_pathName[,c(2,1,3)],"ENSTBEP_pathID_pathName",sep="\t",col.names=T,row.names=F,quote=F)
geneID_GO<-read.delim("geneID_GO",header=T,sep="\t",check.names=F)
names(ENSTBEP_pathID_pathName)[2]<-"Ensembl Protein ID"
merge(geneID_GO,ENSTBEP_pathID_pathName,by="Ensembl Protein ID")->geneIDs_GO_path
TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName<-geneIDs_GO_path[,c(3,2,4,5,6,7,8)]
write.table(TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName,"TranscriptID_GeneID_ProteinID_GeneName_GOID_GOName_pathID_pathName",sep="\t",col.names=
T,row.names=F,quote=F)
q()
ls()
ENSTBEP_paths_new[1:4,]
history()
geneIDs_GO_path[1:4,]
dir()
Ts_EnsembleProteinID_KO<-read.dleim("Ts_EnsembleProteinID_KO",header=T,sep="\t")
Ts_EnsembleProteinID_KO<-read.delim("Ts_EnsembleProteinID_KO",header=T,sep="\t")
names(Ts_EnsembleProteinID_KO)
names(Ts_EnsembleProteinID_KO)<-c("Ensembl Protein ID","ko")
merge(Ts_EnsembleProteinID_KO,geneIDs_GO_path,by="Ensembl Protein ID")->Ts_EnsembleIDs_KO_GO_Path
Ts_EnsembleIDs_KO_GO_Path[1:4,]
names(Ts_EnsembleIDs_KO_GO_Path)
Ts_EnsembleIDs_KO_Name_GO_Pathway<-Ts_EnsembleIDs_KO_GO_Path[,c(1,4,3,2,5,6:9)]
write.table(Ts_EnsembleIDs_KO_Name_GO_Pathway,"Ts_EnsembleIDs_KO_Name_GO_Pathway",sep="\t",col.names=T,row.names=F,quote=F)
q()
    NumberofSDEG[3]<-paste(length(UpDegLines)," / ",dim(DEG_output)[1],sep="")    ## Up
    NumberofSDEG[4]<-paste(length(DownDegLines)," / ",dim(DEG_output)[1],sep="")  ## Down
    return(NumberofSDEG)
}
NumSDEG<-c("group","all","up","down")
rbind(NumSDEG,NumSDEG)
### 添加注释：
NCBI_annotation_path<-"/mnt/lustre/users/wangyan/project/wangwei/gbk_annotation/NCBI_annotation"
stat_path<-"/mnt/lustre/users/wangyan/project/wangwei/wangwei_5_H130530049_20130626"
NumSDEG<-c("group","all","up","down")
for(i in 1:length(EdgeR_DE_result)){
    print(groupname[i])
    edgeR_DE_result<-EdgeR_DE_result[[i]]
    AA<-mergeDegResultAndNcbiPttAnnotation(DegResultPath=edgeR_DE_result,NCBI_annotation_path=NCBI_annotation_path,OutPath=out[i],logFC=1,pvalue=1e-03,group=names(EdgeR_DE_result)[i],stat_path=stat_path)
    NumSDEG<-rbind(NumSDEG,AA)
}
NumSDEG
write.table(NumSDEG,"NumofSigDiffExpGenes.xls",sep="\t",col.names=F,row.names=F,quote=F)
q()
grouplist<-as.character(read.delim("exp_paths",header=F)[[1]])
out<-as.character(read.delim("out",header=F)[[1]])
groupname<-sapply(out,function(x) unlist(strsplit(x,"/",fix=T))[10])
names(groupname)<-NULL
library(edgeR)
EdgeR_DE_result<-list()
i=1
inputpath<-grouplist[i]
outputpath<-out[i]
PValue<-1e-03
logFC<-1
FDR<-0
dir.create(outputpath)
data = read.delim(inputpath, header=T,sep="\t",check.names=F)
rnaseqMatrix = data[,2:3]
rownames(rnaseqMatrix)<-as.character(data[,1])
rnaseqMatrix = round(rnaseqMatrix)
##rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]
samples<-names(data)[2:3]
conditions = factor(c(rep(samples[1], 1), rep(samples[2], 1)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
et = exactTest(exp_study, dispersion=0.1)
tTags = topTags(et,n=NULL)
# plot
outputpath
paste(outputpath,"/MA_n_Volcano.pdf",sep="")
result_table = tTags$table
names(result_table)
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""))
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20) {
par(def.par)
dev.off()
plot_MA = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot", pch=20) {
    plot(logCounts, logFoldChange, col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);;
}
plot_Volcano = function(logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20) {
   plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);
}
plot_MA_and_Volcano = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot") {
    def.par = par(no.readonly = TRUE) # save default, for resetting...
    gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
    layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 
    plot_MA(logCounts, logFoldChange, FDR);
    plot_Volcano(logFoldChange, FDR);
    par(def.par)           
}
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""))
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
outputpath
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""))
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=6)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=5)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
getSDERole<-function(logFC,pvalue,FDR){
# FDR
if(logFC==0 && pvalue==0 && FDR!=0){SDERole<-paste("significant(FDR<=",FDR,")",sep="")}
# pvalue
if(logFC==0 && pvalue!=0 && FDR==0){SDERole<-paste("significant(pvalue<=",pvalue,")",sep="")}
# logFC
if(logFC!=0 && pvalue==0 && FDR==0){SDERole<-paste("significant(|logFC|>=",logFC,")",sep="")}
# FDR && logFC
if(logFC!=0 && pvalue==0 && FDR!=0){SDERole<-paste("significant(FDR<=",FDR,",|logFC|>=",logFC,")",sep="")}
# pvalue && logFC
if(logFC!=0 && pvalue!=0 && FDR==0){SDERole<-paste("significant(pvalue<=",pvalue,",|logFC|>=",logFC,")",sep="")}
# pvalue && FDR
if(logFC==0 && pvalue!=0 && FDR!=0){SDERole<-paste("significant(FDR<=",FDR,",pvalue>=",pvalue,")",sep="")}
# pvalue && FDR && logFC
if(logFC!=0 && pvalue!=0 && FDR!=0){SDERole<-paste("significant(FDR<=",FDR,",|logFC|>=",logFC,",pvalue<=",pvalue,")",sep="")}
return(SDERole)
}
plot_MA = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot", pch=20) {
    plot(logCounts, logFoldChange, col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);;
}
plot_Volcano = function(logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20) {
   plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);
}
plot_MA_and_Volcano = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot") {
    def.par = par(no.readonly = TRUE) # save default, for resetting...
    gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
    layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 
    plot_MA(logCounts, logFoldChange, FDR);
    plot_Volcano(logFoldChange, FDR);
    par(def.par)           
}
ww<-function(inputpath,outputpath,PValue,logFC){
dir.create(outputpath)
data = read.delim(inputpath, header=T,sep="\t",check.names=F)
rnaseqMatrix = data[,2:3]
rownames(rnaseqMatrix)<-as.character(data[,1])
rnaseqMatrix = round(rnaseqMatrix)
##rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]
samples<-names(data)[2:3]
conditions = factor(c(rep(samples[1], 1), rep(samples[2], 1)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
et = exactTest(exp_study, dispersion=0.1)
tTags = topTags(et,n=NULL)
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=5)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
# add counts:
tTags_new<-cbind(rownames(tTags$table),tTags$table)
names(tTags_new)[1]<-"Synonym"
tTags_new_new<-merge(data,tTags_new,by="Synonym")
sigLine<-intersect(which(tTags_new_new[["PValue"]]<=PValue),which(abs(tTags_new_new[["logFC"]])>=logFC))
getSDERole(logFC,PValue,FDR)
Significant<-rep(FALSE,dim(tTags_new_new)[1])
Significant[sigLine]<-"TRUE"
logFC2<--tTags_new_new[["logFC"]]
tTags_new_new_new<-cbind(tTags_new_new,logFC2,Significant)
names(tTags_new_new_new)[ncol(tTags_new_new_new)]<-getSDERole(logFC,PValue,FDR)
Name<-names(tTags_new_new_new)
names(tTags_new_new_new)[4]<-paste("logFC(",Name[3],"/",Name[2],")",sep="")
names(tTags_new_new_new)[8]<-paste("logFC(",Name[2],"/",Name[3],")",sep="")
write.table(tTags_new_new_new, file=paste(outputpath,"/edgeR.DE.result.xls",sep=""), sep='\t', quote=F, row.names=F,col.names=T)
## direction
Num_Sig<-which(tTags_new_new_new[[ncol(tTags_new_new_new)]]=="TRUE")
Sig_tTags_new_new_new<-tTags_new_new_new[Num_Sig,]
Num_up<-length(which(Sig_tTags_new_new_new[[8]]>0))
Num_down<-length(which(Sig_tTags_new_new_new[[8]]<0))
return(tTags_new_new_new)
}
mergeDegResultAndNcbiPttAnnotation<-function(DegResultPath,NCBI_annotation_path,OutPath,logFC,pvalue,group){
NumberofSDEG<-vector()
value1<-unlist(strsplit(group,"_vs_",fix=T))[1]
value2<-unlist(strsplit(group,"_vs_",fix=T))[2]
NCBI_annotation_all<-read.delim(NCBI_annotation_path,header=T,sep="\t",check.names=F)
names(NCBI_annotation_all)[5]<-"GeneName"
    Result<-unique(merge(DegResultPath,NCBI_annotation_all,by="Synonym"))
    SigDegLines<-which(Result[[ncol(DegResultPath)]]=="TRUE")
    write.table(Result,paste(OutPath,"/AllGeneNcbiAnnotationAndDegResult.xls",sep=""),sep=
    
    # all sig DEGs
    AllSigDiffExpGeneNcbiAnnotationAndDegResult_new<-as.data.frame(Result[SigDegLines,])
    row.names(AllSigDiffExpGeneNcbiAnnotationAndDegResult_new)<-NULL
    AllSigDiffExpGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult_
    write.table(AllSigDiffExpGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllSigDiffExpGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all up sig DEGs
    UpDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))>0)
    AllUpRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[UpDegLines,]
    write.table(AllUpRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllUpRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all down sig DEGs
    DownDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))<0)
    AllDownRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[DownDegLines,]
    write.table(AllDownRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllDownRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    NumberofSDEG[1]<-group
    NumberofSDEG[2]<-paste(length(SigDegLines)," / ",dim(DEG_output)[1],sep="")   ## All
    NumberofSDEG[3]<-paste(length(UpDegLines)," / ",dim(DEG_output)[1],sep="")    ## Up
    NumberofSDEG[4]<-paste(length(DownDegLines)," / ",dim(DEG_output)[1],sep="")  ## Down
    return(NumberofSDEG)
}
 ww<-function(inputpath,outputpath,PValue,logFC){
dir.create(outputpath)
data = read.delim(inputpath, header=T,sep="\t",check.names=F)
rnaseqMatrix = data[,2:3]
rownames(rnaseqMatrix)<-as.character(data[,1])
rnaseqMatrix = round(rnaseqMatrix)
##rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]
samples<-names(data)[2:3]
conditions = factor(c(rep(samples[1], 1), rep(samples[2], 1)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
et = exactTest(exp_study, dispersion=0.1)
tTags = topTags(et,n=NULL)
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=5)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
# add counts:
tTags_new<-cbind(rownames(tTags$table),tTags$table)
names(tTags_new)[1]<-"Synonym"
tTags_new_new<-merge(data,tTags_new,by="Synonym")
sigLine<-intersect(which(tTags_new_new[["PValue"]]<=PValue),which(abs(tTags_new_new[["logFC"]])>=logFC))
getSDERole(logFC,PValue,FDR)
Significant<-rep(FALSE,dim(tTags_new_new)[1])
Significant[sigLine]<-"TRUE"
logFC2<--tTags_new_new[["logFC"]]
tTags_new_new_new<-cbind(tTags_new_new,logFC2,Significant)
names(tTags_new_new_new)[ncol(tTags_new_new_new)]<-getSDERole(logFC,PValue,FDR)
Name<-names(tTags_new_new_new)
names(tTags_new_new_new)[4]<-paste("logFC(",Name[3],"/",Name[2],")",sep="")
names(tTags_new_new_new)[8]<-paste("logFC(",Name[2],"/",Name[3],")",sep="")
write.table(tTags_new_new_new, file=paste(outputpath,"/edgeR.DE.result.xls",sep=""), sep='\t', quote=F, row.names=F,col.names=T)
## direction
Num_Sig<-which(tTags_new_new_new[[ncol(tTags_new_new_new)]]=="TRUE")
Sig_tTags_new_new_new<-tTags_new_new_new[Num_Sig,]
Num_up<-length(which(Sig_tTags_new_new_new[[8]]>0))
Num_down<-length(which(Sig_tTags_new_new_new[[8]]<0))
return(tTags_new_new_new)
}
mergeDegResultAndNcbiPttAnnotation<-function(DegResultPath,NCBI_annotation_path,OutPath,logFC,pvalue,group){
NumberofSDEG<-vector()
value1<-unlist(strsplit(group,"_vs_",fix=T))[1]
value2<-unlist(strsplit(group,"_vs_",fix=T))[2]
NCBI_annotation_all<-read.delim(NCBI_annotation_path,header=T,sep="\t",check.names=F)
names(NCBI_annotation_all)[5]<-"GeneName"
    Result<-unique(merge(DegResultPath,NCBI_annotation_all,by="Synonym"))
    SigDegLines<-which(Result[[ncol(DegResultPath)]]=="TRUE")
    write.table(Result,paste(OutPath,"/AllGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    
    # all sig DEGs
    AllSigDiffExpGeneNcbiAnnotationAndDegResult_new<-as.data.frame(Result[SigDegLines,])
    row.names(AllSigDiffExpGeneNcbiAnnotationAndDegResult_new)<-NULL
    AllSigDiffExpGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult_new
    write.table(AllSigDiffExpGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllSigDiffExpGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all up sig DEGs
    UpDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))>0)
    AllUpRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[UpDegLines,]
    write.table(AllUpRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllUpRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all down sig DEGs
    DownDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))<0)
    AllDownRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[DownDegLines,]
    write.table(AllDownRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllDownRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    NumberofSDEG[1]<-group
    NumberofSDEG[2]<-paste(length(SigDegLines)," / ",dim(DEG_output)[1],sep="")   ## All
    NumberofSDEG[3]<-paste(length(UpDegLines)," / ",dim(DEG_output)[1],sep="")    ## Up
    NumberofSDEG[4]<-paste(length(DownDegLines)," / ",dim(DEG_output)[1],sep="")  ## Down
    return(NumberofSDEG)
}
grouplist<-as.character(read.delim("exp_paths",header=F)[[1]])
out<-as.character(read.delim("out",header=F)[[1]])
groupname<-sapply(out,function(x) unlist(strsplit(x,"/",fix=T))[10])
names(groupname)<-NULL
library(edgeR)
EdgeR_DE_result<-list()
for(i in 1:length(out)){
inputpath<-grouplist[i]
outputpath<-out[i]
PValue<-1e-03
logFC<-1
FDR<-0
result<-ww(inputpath,outputpath,PValue,logFC)
EdgeR_DE_result[[i]]<-result
}
names(EdgeR_DE_result)<-groupname
save(EdgeR_DE_result,file="EdgeR_DE_result.RData")
### 添加注释：
NCBI_annotation_path<-"/mnt/lustre/users/wangyan/project/wangwei/gbk_annotation/NCBI_annotation"
NumSDEG<-c("group","all","up","down")
for(i in 1:length(EdgeR_DE_result)){
    print(groupname[i])
    edgeR_DE_result<-EdgeR_DE_result[[i]]
    AA<-mergeDegResultAndNcbiPttAnnotation(DegResultPath=edgeR_DE_result,NCBI_annotation_path=NCBI_annotation_path,OutPath=out[i],logFC=1,pvalue=1e-03,group=names(EdgeR_DE_result)[i])
    NumSDEG<-rbind(NumSDEG,AA)
}
write.table(NumSDEG,"NumofSigDiffExpGenes.xls",sep="\t",col.names=F,row.names=F,quote=F)
### plot:
NumSDEG
mergeDegResultAndNcbiPttAnnotation<-function(DegResultPath,NCBI_annotation_path,OutPath,logFC,pvalue,group){
NumberofSDEG<-vector()
value1<-unlist(strsplit(group,"_vs_",fix=T))[1]
value2<-unlist(strsplit(group,"_vs_",fix=T))[2]
NCBI_annotation_all<-read.delim(NCBI_annotation_path,header=T,sep="\t",check.names=F)
names(NCBI_annotation_all)[5]<-"GeneName"
    Result<-unique(merge(DegResultPath,NCBI_annotation_all,by="Synonym"))
    SigDegLines<-which(Result[[ncol(DegResultPath)]]=="TRUE")
    write.table(Result,paste(OutPath,"/AllGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    
    # all sig DEGs
    AllSigDiffExpGeneNcbiAnnotationAndDegResult_new<-as.data.frame(Result[SigDegLines,])
    row.names(AllSigDiffExpGeneNcbiAnnotationAndDegResult_new)<-NULL
    AllSigDiffExpGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult_new
    write.table(AllSigDiffExpGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllSigDiffExpGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all up sig DEGs
    UpDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))>0)
    AllUpRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[UpDegLines,]
    write.table(AllUpRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllUpRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all down sig DEGs
    DownDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))<0)
    AllDownRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[DownDegLines,]
    write.table(AllDownRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllDownRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    #NumberofSDEG[1]<-group
    #NumberofSDEG[2]<-paste(length(SigDegLines)," / ",dim(DEG_output)[1],sep="")   ## All
    #NumberofSDEG[3]<-paste(length(UpDegLines)," / ",dim(DEG_output)[1],sep="")    ## Up
    #NumberofSDEG[4]<-paste(length(DownDegLines)," / ",dim(DEG_output)[1],sep="")  ## Down
    
    NumberofSDEG[2]<-length(SigDegLines)
    NumberofSDEG[3]<-length(UpDegLines)
    NumberofSDEG[4]<-length(DownDegLines)
    return(NumberofSDEG)
}
NCBI_annotation_path<-"/mnt/lustre/users/wangyan/project/wangwei/gbk_annotation/NCBI_annotation"
NumSDEG<-c("group","all","up","down")
NumSDEG<-vector()
for(i in 1:length(EdgeR_DE_result)){
    print(groupname[i])
    edgeR_DE_result<-EdgeR_DE_result[[i]]
    AA<-mergeDegResultAndNcbiPttAnnotation(DegResultPath=edgeR_DE_result,NCBI_annotation_path=NCBI_annotation_path,OutPath=out[i],logFC=1,pvalue=1e-03,group=names(EdgeR_DE_result)[i])
    NumSDEG<-rbind(NumSDEG,AA)
}
NumSDEG
mergeDegResultAndNcbiPttAnnotation<-function(DegResultPath,NCBI_annotation_path,OutPath,logFC,pvalue,group){
NumberofSDEG<-vector()
value1<-unlist(strsplit(group,"_vs_",fix=T))[1]
value2<-unlist(strsplit(group,"_vs_",fix=T))[2]
NCBI_annotation_all<-read.delim(NCBI_annotation_path,header=T,sep="\t",check.names=F)
names(NCBI_annotation_all)[5]<-"GeneName"
    Result<-unique(merge(DegResultPath,NCBI_annotation_all,by="Synonym"))
    SigDegLines<-which(Result[[ncol(DegResultPath)]]=="TRUE")
    write.table(Result,paste(OutPath,"/AllGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    
    # all sig DEGs
    AllSigDiffExpGeneNcbiAnnotationAndDegResult_new<-as.data.frame(Result[SigDegLines,])
    row.names(AllSigDiffExpGeneNcbiAnnotationAndDegResult_new)<-NULL
    AllSigDiffExpGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult_new
    write.table(AllSigDiffExpGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllSigDiffExpGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all up sig DEGs
    UpDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))>0)
    AllUpRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[UpDegLines,]
    write.table(AllUpRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllUpRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all down sig DEGs
    DownDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))<0)
    AllDownRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[DownDegLines,]
    write.table(AllDownRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllDownRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    NumberofSDEG[1]<-group
    #NumberofSDEG[2]<-paste(length(SigDegLines)," / ",dim(DEG_output)[1],sep="")   ## All
    #NumberofSDEG[3]<-paste(length(UpDegLines)," / ",dim(DEG_output)[1],sep="")    ## Up
    #NumberofSDEG[4]<-paste(length(DownDegLines)," / ",dim(DEG_output)[1],sep="")  ## Down
    
    NumberofSDEG[2]<-length(SigDegLines)
    NumberofSDEG[3]<-length(UpDegLines)
    NumberofSDEG[4]<-length(DownDegLines)
    return(NumberofSDEG)
}
NCBI_annotation_path<-"/mnt/lustre/users/wangyan/project/wangwei/gbk_annotation/NCBI_annotation"
NumSDEG<-c("group","all","up","down")
NumSDEG<-vector()
for(i in 1:length(EdgeR_DE_result)){
    print(groupname[i])
    edgeR_DE_result<-EdgeR_DE_result[[i]]
    AA<-mergeDegResultAndNcbiPttAnnotation(DegResultPath=edgeR_DE_result,NCBI_annotation_path=NCBI_annotation_path,OutPath=out[i],logFC=1,pvalue=1e-03,group=names(EdgeR_DE_result)[i])
    NumSDEG<-rbind(NumSDEG,AA)
}
NumSDEG
NumSDEG[,2:4]
write.table(NumSDEG,"NumofSigDiffExpGenes.xls",sep="\t",col.names=F,row.names=F,quote=F)
names(NumSDEG)<-c("group","all","up","down")
write.table(NumSDEG,"NumofSigDiffExpGenes.xls",sep="\t",col.names=T,row.names=F,quote=F)
NumSDEG
NumSDEG<-vector()
for(i in 1:length(EdgeR_DE_result)){
    print(groupname[i])
    edgeR_DE_result<-EdgeR_DE_result[[i]]
    AA<-mergeDegResultAndNcbiPttAnnotation(DegResultPath=edgeR_DE_result,NCBI_annotation_path=NCBI_annotation_path,OutPath=out[i],logFC=1,pvalue=1e-03,group=names(EdgeR_DE_result)[i])
    NumSDEG<-rbind(NumSDEG,AA)
}
NumSDEG
dim(NumSDEG)
class(NumSDEG)
?jpeg
col.names(NumSDEG)<-c("group","all","up","down")
colnames(NumSDEG)<-c("group","all","up","down")
NumSDEG
NumSDEG_new<-NumSDEG[,-1]
colnames(NumSDEG_new)<-c("all","up","down")
rownames(NumSDEG_new)<-as.character(NumSDEG[[1]])
NumSDEG_new<-NumSDEG[,-1]
NumSDEG_new
as.character(NumSDEG[[1]])
rownames(NumSDEG_new)<-as.character(NumSDEG[,1])
NumSDEG_new
write.table(NumSDEG_new,"NumofSigDiffExpGenes.xls",sep="\t",col.names=T,row.names=T,quote=F)
?pdf
?jpeg
?read.delim
?read.table
NumSDEG_new<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1, col.names=1)
NumSDEG_new
NumSDEG_new<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
NumSDEG_new
dim(NumSDEG_new)
names(NumSDEG_new)
rownames(NumSDEG_new)
q()
NSDE
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,3000),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
??brewer.pal
library(ColorBrewer)
NSDE<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
pdf("Num of DEGs Among Each Group.pdf",w=11,h=8)
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,3000),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
dev.off()
library("ColorBrewer")
library(RColorBrewer)
NSDE<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
pdf("Num of DEGs Among Each Group.pdf",w=11,h=8)
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,3000),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
dev.off()
q()
# 柱状图
library(RColorBrewer)
NSDE<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
pdf("NumofDEGsAmongEachGroup_EdgeR.pdf",w=11,h=8)
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,3000),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
dev.off()
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,3000),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
NSDE
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,1000),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
library(RColorBrewer)
NSDE<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
pdf("NumofDEGsAmongEachGroup_EdgeR.pdf",w=11,h=8)
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,800),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
dev.off()
q()
pathway<-read.delim("pathway.txt.new",sep="\t",header=F)
naems(pathway)
names(pathway)<-c("")
geneid2Synonym<-read.delim("geneid2Synonym",sep="\t",header=T)
geneid2Synonym[1:4,]
names(pathway)<-c("Gene","KOs")
merge(pathway,geneid2Synonym,by="Gene")->new
new[1:4,]
write.table(new[,c(3,2,1)],"Synonym2KOs2Gene",sep="\t",col.names=F,row.names=F,quote=F)
q()
history()
GO<-read.delim("GO.list",sep="\t",header=F)
GO[1:4,]
ls()
history()
geneid2Synonym[1:4,]
GO[1:4,]
names(GO)
names(GO)<-("","GOs")
names(GO)<-("Gene","GOs")
names(GO)<-c("Gene","GOs")
merge(geneid2Synonym,GO,by="Gene")->go_new
go_new[1:4,]
write.table(go_new,"GO.list.new",sep="\t",col.names=F,row.names=F,quote=F)
q()
zt4_matrix<-read.delim("ZT_map_gene.sam.gene.exp_multi.xls",sep="\t",header=T)
allsample_matrix<-read.delim("6sample_exp_matrix.xls",sep="\t",header=T)
zt4_matrix[1:4,]
allsample_matrix[1:4,]
names(zt4_matrix)
names(allsample_matrix)
merge(allsample_matrix,zt4_matrix,by="Gene")->All7Sample_Matrix
All7Sample_Matrix[1:4,]
names(All7Sample_Matrix)
names(All7Sample_Matrix)<-c("gene_id","16wan1","16wan2","16wan3","zt1","zt2","zt3","zt4")
names(All7Sample_Matrix)
write.table(All7Sample_Matrix,"All7Sample_Matrix.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
ls()
load("all.RData")
ls()
hc_samples
pdf("7Samples_Cluster_Result.pdf")
plot(hc_samples1,hang=-1,xlab="",ylab="",main="",sub="")
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="")
dev.off()
pdf("7Samples_Cluster_Result.pdf",w=8,h=6)
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="")
dev.off()
pdf("7Samples_Cluster_Result.pdf",w=6,h=4)
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="")
dev.off()
q()
unmapname<-scan("unmapped.fq.list",what=character())
unmapname[1:10]
?sample
x <- 1:12
     # a random permutation
     sample(x)
x
sample(x)
sample(x,4)
x <- 1:10000
sample(x,20)
randNum<-sample(1:length(unmapname),800000)
randNum[1:4]
unmapname_rand<-unmapname[randNum]
write.table(unmapname_rand,"unmapped.fq.list.rand",sep="\n",quote=F,col.names=F,row.names=F)
q()
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
# add counts:
tTags_new<-cbind(rownames(tTags$table),tTags$table)
names(tTags_new)[1]<-"Synonym"
tTags_new_new<-merge(data,tTags_new,by="Synonym")
sigLine<-intersect(which(tTags_new_new[["PValue"]]<=PValue),which(abs(tTags_new_new[["logFC"]])>=logFC))
getSDERole(logFC,PValue,FDR)
Significant<-rep(FALSE,dim(tTags_new_new)[1])
Significant[sigLine]<-"TRUE"
logFC2<--tTags_new_new[["logFC"]]
tTags_new_new_new<-cbind(tTags_new_new,logFC2,Significant)
names(tTags_new_new_new)[ncol(tTags_new_new_new)]<-getSDERole(logFC,PValue,FDR)
Name<-names(tTags_new_new_new)
names(tTags_new_new_new)[4]<-paste("logFC(",Name[3],"/",Name[2],")",sep="")
names(tTags_new_new_new)[8]<-paste("logFC(",Name[2],"/",Name[3],")",sep="")
write.table(tTags_new_new_new, file=paste(outputpath,"/edgeR.DE.result.xls",sep=""), sep='\t', quote=F, row.names=F,col.names=T)
## direction
Num_Sig<-which(tTags_new_new_new[[ncol(tTags_new_new_new)]]=="TRUE")
Sig_tTags_new_new_new<-tTags_new_new_new[Num_Sig,]
Num_up<-length(which(Sig_tTags_new_new_new[[8]]>0))
Num_down<-length(which(Sig_tTags_new_new_new[[8]]<0))
return(tTags_new_new_new)
}
mergeDegResultAndNcbiPttAnnotation<-function(DegResultPath,NCBI_annotation_path,OutPath,logFC,pvalue,group){
NumberofSDEG<-vector()
value1<-unlist(strsplit(group,"_vs_",fix=T))[1]
value2<-unlist(strsplit(group,"_vs_",fix=T))[2]
NCBI_annotation_all<-read.delim(NCBI_annotation_path,header=T,sep="\t",check.names=F)
names(NCBI_annotation_all)[5]<-"GeneName"
    Result<-unique(merge(DegResultPath,NCBI_annotation_all,by="Synonym"))
    SigDegLines<-which(
    write.table(Result,paste(OutPath,"/AllGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    
    # all sig DEGs
    AllSigDiffExpGeneNcbiAnnotationAndDegResult_new<-as.data.frame(Result[SigDegLines,])
    row.names(AllSigDiffExpGeneNcbiAnnotationAndDegResult_new)<-NULL
    AllSigDiffExpGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult
    write.table(AllSigDiffExpGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllSigDiffExpGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all up sig DEGs
    UpDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))>0)
    AllUpRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[UpDegLines,]
    write.table(AllUpRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllUpRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all down sig DEGs
    DownDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))<0)
    AllDownRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[DownDegLines,
    write.table(AllDownRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllDownRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    NumberofSDEG[1]<
    #NumberofSDEG[2]<-paste(length(SigDegLines)," / ",dim(DEG_output)[1],sep="")   ## All
    #NumberofSDEG[3]<-paste(length(UpDegLines)," / ",dim(DEG_output)[1],sep="")    ## Up
    #NumberofSDEG[4]<-paste(length(DownDegLines)," / ",dim(DEG_output)[1],sep="")  ## Down
    
    NumberofSDEG[2]<-length(SigDegLine
    NumberofSDEG[3]<-length(UpDegLines)
    NumberofSDEG[4]<-length(DownDegLines)
    return(NumberofSDEG)
}
mergeDegResultAndNcbiPttAnnotation<-function(DegResultPath,NCBI_annotation_path,OutPath,logFC,pvalue,group){
NumberofSDEG<-vector()
value1<-unlist(strsplit(group,"_vs_",fix=T))[1]
value2<-unlist(strsplit(group,"_vs_",fix=T))[2]
NCBI_annotation_all<-read.delim(NCBI_annotation_path,header=T,sep="\t",check.names=F)
names(NCBI_annotation_all)[5]<-"GeneName"
    Result<-unique(merge(DegResultPath,NCBI_annotation_all,by="Synonym"))
    SigDegLines<-which(Result[[ncol(DegResultPath)]]=="TRUE")
    write.table(Result,paste(OutPath,"/AllGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    
    # all sig DEGs
    AllSigDiffExpGeneNcbiAnnotationAndDegResult_new<-as.data.frame(Result[SigDegLines,])
    row.names(AllSigDiffExpGeneNcbiAnnotationAndDegResult_new)<-NULL
    AllSigDiffExpGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult_new
    write.table(AllSigDiffExpGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllSigDiffExpGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all up sig DEGs
    UpDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))>0)
    AllUpRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[UpDegLines,]
    write.table(AllUpRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllUpRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    # all down sig DEGs
    DownDegLines<-which(as.numeric(as.character(AllSigDiffExpGeneNcbiAnnotationAndDegResult[[8]]))<0)
    AllDownRegulatedGeneNcbiAnnotationAndDegResult<-AllSigDiffExpGeneNcbiAnnotationAndDegResult[DownDegLines,]
    write.table(AllDownRegulatedGeneNcbiAnnotationAndDegResult,paste(OutPath,"/AllDownRegulatedGeneNcbiAnnotationAndDegResult.xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
    NumberofSDEG[1]<-group
    #NumberofSDEG[2]<-paste(length(SigDegLines)," / ",dim(DEG_output)[1],sep="")   ## All
    #NumberofSDEG[3]<-paste(length(UpDegLines)," / ",dim(DEG_output)[1],sep="")    ## Up
    #NumberofSDEG[4]<-paste(length(DownDegLines)," / ",dim(DEG_output)[1],sep="")  ## Down
    
    NumberofSDEG[2]<-length(SigDegLines)
    NumberofSDEG[3]<-length(UpDegLines)
    NumberofSDEG[4]<-length(DownDegLines)
    return(NumberofSDEG)
}
#### 换edgeR方法计算差异表达
# 修改表头
grouplist<-as.character(read.delim("exp_paths_new",header=F)[[1]])
out<-as.character(read.delim("out",header=F)[[1]])
groupname<-sapply(out,function(x) unlist(strsplit(x,"/",fix=T))[10])
names(groupname)<-NULL
library(edgeR)
EdgeR_DE_result<-list()
for(i in 1:length(out)){
inputpath<-grouplist[i]
outputpath<-out[i]
PValue<-1e-03
logFC<-1
FDR<-0
result<-ww(inputpath,outputpath,PValue,logFC)
EdgeR_DE_result[[i]]<-result
}
names(EdgeR_DE_result)<-groupname
save(EdgeR_DE_result,file="EdgeR_DE_result.RData")
groupname
grouplist
unlist(strsplit(grouplist[1],"/",fix=T))
unlist(strsplit(grouplist[1],"/",fix=T))[10]
q()
groupname
grouplist<-as.character(read.delim("exp_paths_new",header=F)[[1]])
out<-as.character(read.delim("out",header=F)[[1]])
groupname<-sapply(out,function(x) unlist(strsplit(x,"/",fix=T))[9])
names(groupname)<-NULL
library(edgeR)
EdgeR_DE_result<-list()
groupname
inputpath
i
data = read.delim(inputpath, header=T,sep="\t",check.names=F,skip=1)
data[1:4,]
dim(data)
groupname[i]
strsplit(groupname[i],"-",fix=T)
unlist(strsplit(groupname[i],"-",fix=T))
dir.create(outputpath)
data = read.delim(inputpath, header=T,sep="\t",check.names=F,skip=1)
names(data)<-c("Gene",unlist(strsplit(groupname[i],"-",fix=T)))
rnaseqMatrix = data[,2:3]
rownames(rnaseqMatrix)<-as.character(data[,1])
rnaseqMatrix = round(rnaseqMatrix)
##rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]
samples<-names(data)[2:3]
conditions = factor(c(rep(samples[1], 1), rep(samples[2], 1)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
et = exactTest(exp_study, dispersion=0.1)
tTags = topTags(et,n=NULL)
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=5)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
q()
ww<-function(inputpath,outputpath,PValue,logFC,name){
dir.create(outputpath)
data = read.delim(inputpath, header=T,sep="\t",check.names=F,skip=1)
names(data)<-c("Gene",unlist(strsplit(groupname[i],"-",fix=T)))
rnaseqMatrix = data[,2:3]
rownames(rnaseqMatrix)<-as.character(data[,1])
rnaseqMatrix = round(rnaseqMatrix)
##rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]
samples<-names(data)[2:3]
conditions = factor(c(rep(samples[1], 1), rep(samples[2], 1)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
et = exactTest(exp_study, dispersion=0.1)
tTags = topTags(et,n=NULL)
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=5)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
# add counts:
tTags_new<-cbind(rownames(tTags$table),tTags$table)
names(tTags_new)[1]<-"Synonym"
tTags_new_new<-merge(data,tTags_new,by="Synonym")
sigLine<-intersect(which(tTags_new_new[["PValue"]]<=PValue),which(abs(tTags_new_new[["logFC"]])>=logFC))
getSDERole(logFC,PValue,FDR)
Significant<-rep(FALSE,dim(tTags_new_new)[1])
Significant[sigLine]<-"TRUE"
logFC2<--tTags_new_new[["logFC"]]
tTags_new_new_new<-cbind(tTags_new_new,logFC2,Significant)
names(tTags_new_new_new)[ncol(tTags_new_new_new)]<-getSDERole(logFC,PValue,FDR)
Name<-names(tTags_new_new_new)
names(tTags_new_new_new)[4]<-paste("logFC(",Name[3],"/",Name[2],")",sep="")
names(tTags_new_new_new)[8]<-paste("logFC(",Name[2],"/",Name[3],")",sep="")
write.table(tTags_new_new_new, file=paste(outputpath,"/edgeR.DE.result.xls",sep=""), sep='\t', quote=F, row.names=F,col.names=T)
## direction
Num_Sig<-which(tTags_new_new_new[[ncol(tTags_new_new_new)]]=="TRUE")
Sig_tTags_new_new_new<-tTags_new_new_new[Num_Sig,]
Num_up<-length(which(Sig_tTags_new_new_new[[8]]>0))
Num_down<-length(which(Sig_tTags_new_new_new[[8]]<0))
return(tTags_new_new_new)
}
q()
grouplist<-as.character(read.delim("exp_paths_new",header=F)[[1]])
out<-as.character(read.delim("out",header=F)[[1]])
groupname<-sapply(out,function(x) unlist(strsplit(x,"/",fix=T))[9])
names(groupname)<-NULL
library(edgeR)
EdgeR_DE_result<-list()
for(i in 1:length(out)){
inputpath<-grouplist[i]
outputpath<-out[i]
PValue<-1e-03
logFC<-1
FDR<-0
result<-ww(inputpath,outputpath,PValue,logFC,name=groupname[i])
EdgeR_DE_result[[i]]<-result
}
names(EdgeR_DE_result)<-groupname
save(EdgeR_DE_result,file="EdgeR_DE_result.RData")
i
inputpath<-grouplist[i]
outputpath<-out[i]
PValue<-1e-03
logFC<-1
FDR<-0
result<-ww(inputpath,outputpath,PValue,logFC,name=groupname[i])
EdgeR_DE_result[[i]]<-result
i
dir.create(outputpath)
data = read.delim(inputpath, header=T,sep="\t",check.names=F,skip=1)
names(data)<-c("Gene",unlist(strsplit(groupname[i],"-",fix=T)))
rnaseqMatrix = data[,2:3]
rownames(rnaseqMatrix)<-as.character(data[,1])
rnaseqMatrix = round(rnaseqMatrix)
##rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]
samples<-names(data)[2:3]
conditions = factor(c(rep(samples[1], 1), rep(samples[2], 1)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
et = exactTest(exp_study, dispersion=0.1)
tTags = topTags(et,n=NULL)
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=5)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
tTags_new<-cbind(rownames(tTags$table),tTags$table)
tTags_new[1:4,]
q()
grouplist<-as.character(read.delim("exp_paths",header=F)[[1]])
out<-as.character(read.delim("out",header=F)[[1]])
groupname<-sapply(out,function(x) unlist(strsplit(x,"/",fix=T))[9])
names(groupname)<-NULL
library(edgeR)
EdgeR_DE_result<-list()
i
inputpath<-grouplist[i]
outputpath<-out[i]
PValue<-1e-03
logFC<-1
FDR<-0
result<-ww(inputpath,outputpath,PValue,logFC,name=groupname[i])
EdgeR_DE_result[[i]]<-result
ww
dir.create(outputpath)
data = read.delim(inputpath, header=T,sep="\t",check.names=F,skip=1)
names(data)<-c("Synonym",unlist(strsplit(groupname[i],"-",fix=T)))
rnaseqMatrix = data[,2:3]
rownames(rnaseqMatrix)<-as.character(data[,1])
rnaseqMatrix = round(rnaseqMatrix)
##rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]
samples<-names(data)[2:3]
conditions = factor(c(rep(samples[1], 1), rep(samples[2], 1)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
et = exactTest(exp_study, dispersion=0.1)
tTags = topTags(et,n=NULL)
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=5)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
tTags_new<-cbind(rownames(tTags$table),tTags$table)
names(tTags_new)[1]<-"Synonym"
tTags_new_new<-merge(data,tTags_new,by="Synonym")
sigLine<-intersect(which(tTags_new_new[["PValue"]]<=PValue),which(abs(tTags_new_new[["logFC"]])>=logFC))
getSDERole(logFC,PValue,FDR)
Significant<-rep(FALSE,dim(tTags_new_new)[1])
Significant[sigLine]<-"TRUE"
logFC2<--tTags_new_new[["logFC"]]
tTags_new_new_new<-cbind(tTags_new_new,logFC2,Significant)
names(tTags_new_new_new)[ncol(tTags_new_new_new)]<-getSDERole(logFC,PValue,FDR)
Name<-names(tTags_new_new_new)
names(tTags_new_new_new)[4]<-paste("logFC(",Name[3],"/",Name[2],")",sep="")
names(tTags_new_new_new)[8]<-paste("logFC(",Name[2],"/",Name[3],")",sep="")
write.table(tTags_new_new_new, file=paste(outputpath,"/edgeR.DE.result.xls",sep=""), sep='\t', quote=F, row.names=F,col.names=T)
## direction
Num_Sig<-which(tTags_new_new_new[[ncol(tTags_new_new_new)]]=="TRUE")
Sig_tTags_new_new_new<-tTags_new_new_new[Num_Sig,]
Num_up<-length(which(Sig_tTags_new_new_new[[8]]>0))
Num_down<-length(which(Sig_tTags_new_new_new[[8]]<0))
tTags_new_new_new[1:4,]
ww<-function(inputpath,outputpath,PValue,logFC,name){
dir.create(outputpath)
data = read.delim(inputpath, header=T,sep="\t",check.names=F,skip=1)
names(data)<-c("Synonym",unlist(strsplit(name,"-",fix=T)))
rnaseqMatrix = data[,2:3]
rownames(rnaseqMatrix)<-as.character(data[,1])
rnaseqMatrix = round(rnaseqMatrix)
##rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=10,]
samples<-names(data)[2:3]
conditions = factor(c(rep(samples[1], 1), rep(samples[2], 1)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
et = exactTest(exp_study, dispersion=0.1)
tTags = topTags(et,n=NULL)
# plot
source("/share/apps/RNA-Seq/trinityrnaseq_r2013-02-25/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf(file=paste(outputpath,"/MA_n_Volcano.pdf",sep=""),w=10,h=5)
result_table = tTags$table
def.par = par(no.readonly = TRUE) # save default, for resetting...
gridlayout = matrix(c(1:4),nrow=1,ncol=2, byrow=TRUE);
layout(gridlayout, widths=c(1,1), heights=c(1,1)) 
plot_MA(logCounts=result_table$logCPM,logFoldChange=result_table$logFC,FDR=result_table$PValue,xlab="logCPM",ylab="logFC",title="MA plot",pch=20) 
plot_Volcano(logFoldChange=result_table$logFC,FDR=result_table$PValue, xlab="logFC",ylab="-1*log10(Pvalue)", title="Volcano plot", pch=20)
par(def.par)
dev.off()
# add counts:
tTags_new<-cbind(rownames(tTags$table),tTags$table)
names(tTags_new)[1]<-"Synonym"
tTags_new_new<-merge(data,tTags_new,by="Synonym")
sigLine<-intersect(which(tTags_new_new[["PValue"]]<=PValue),which(abs(tTags_new_new[["logFC"]])>=logFC))
getSDERole(logFC,PValue,FDR)
Significant<-rep(FALSE,dim(tTags_new_new)[1])
Significant[sigLine]<-"TRUE"
logFC2<--tTags_new_new[["logFC"]]
tTags_new_new_new<-cbind(tTags_new_new,logFC2,Significant)
names(tTags_new_new_new)[ncol(tTags_new_new_new)]<-getSDERole(logFC,PValue,FDR)
Name<-names(tTags_new_new_new)
names(tTags_new_new_new)[4]<-paste("logFC(",Name[3],"/",Name[2],")",sep="")
names(tTags_new_new_new)[8]<-paste("logFC(",Name[2],"/",Name[3],")",sep="")
write.table(tTags_new_new_new, file=paste(outputpath,"/edgeR.DE.result.xls",sep=""), sep='\t', quote=F, row.names=F,col.names=T)
## direction
Num_Sig<-which(tTags_new_new_new[[ncol(tTags_new_new_new)]]=="TRUE")
Sig_tTags_new_new_new<-tTags_new_new_new[Num_Sig,]
Num_up<-length(which(Sig_tTags_new_new_new[[8]]>0))
Num_down<-length(which(Sig_tTags_new_new_new[[8]]<0))
return(tTags_new_new_new)
}
result<-ww(inputpath,outputpath,PValue,logFC,name=groupname[i])
result[1:4,]
grouplist<-as.character(read.delim("exp_paths",header=F)[[1]])
out<-as.character(read.delim("out",header=F)[[1]])
groupname<-sapply(out,function(x) unlist(strsplit(x,"/",fix=T))[9])
names(groupname)<-NULL
library(edgeR)
EdgeR_DE_result<-list()
for(i in 1:length(out)){
inputpath<-grouplist[i]
outputpath<-out[i]
PValue<-1e-03
logFC<-1
FDR<-0
result<-ww(inputpath,outputpath,PValue,logFC,name=groupname[i])
EdgeR_DE_result[[i]]<-result
}
names(EdgeR_DE_result)<-groupname
save(EdgeR_DE_result,file="EdgeR_DE_result.RData")
EdgeR_DE_result[1:4,]
EdgeR_DE_result[1]
length(EdgeR_DE_result)
### 添加注释 ： 
NCBI_annotation_path<-"/mnt/lustre/users/wangyan/project/wangwei/20131014_PairCompare_EdgeR/AllGene-Annot.txt.new"
NumSDEG<-vector()
for(i in 1:length(EdgeR_DE_result)){
    print(groupname[i])
    edgeR_DE_result<-EdgeR_DE_result[[i]]
    AA<-mergeDegResultAndNcbiPttAnnotation(DegResultPath=edgeR_DE_result,NCBI_annotation_path=NCBI_annotation_path,OutPath=out[i],logFC=1,pvalue=1e-03,group=names(EdgeR_DE_result)[i])
    NumSDEG<-rbind(NumSDEG,AA)
}
NumSDEG_new<-NumSDEG[,-1]
colnames(NumSDEG_new)<-c("all","up","down")
rownames(NumSDEG_new)<-as.character(NumSDEG[,1])
write.table(NumSDEG_new,"NumofSigDiffExpGenes.xls",sep="\t",col.names=T,row.names=T,quote=F)
# 柱状图
NumSDEG_new
library(RColorBrewer)
NSDE<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
NSDE
library(RColorBrewer)
NSDE<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
pdf("NumofDEGsAmongEachGroup_EdgeR.pdf",w=11,h=8)
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
dev.off()
q()
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
library(RColorBrewer)
par(mar=c(5,5,5,3))
barplot(t(as.matrix(NSDE[,2:3])),ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
NSDE
NSDE
plot_matrix<-t(as.matrix(NSDE[,2:3]))
rownames(plot_matrix)
colnames(plot_matrix)
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
colnames(t(as.matrix(NSDE[,2:3])))
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-0.2,length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-0.5,length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-1,length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
NSDE
max(NSDE)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-(max(NSDE)/20),length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
pdf("NumofDEGsAmongEachGroup_EdgeR.pdf",w=20,h=8)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-(max(NSDE)/20),length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
dev.off()
q()
# 柱状图
library(RColorBrewer)
NSDE<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
pdf("NumofDEGsAmongEachGroup_EdgeR.pdf",w=16,h=8)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-(max(NSDE)/20),length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
dev.off()
library(RColorBrewer)
NSDE<-read.delim("NumofSigDiffExpGenes.xls",sep="\t",row.names=1)
pdf("NumofDEGsAmongEachGroup_EdgeR.pdf",w=14,h=8)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-(max(NSDE)/20),length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
dev.off()
NSDE
NSDE[[1]]
c(NSDE[[1]]+max(NSDE)/20)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-(max(NSDE)/20),length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
text(bars,c(NSDE[[1]]+max(NSDE)/20),as.character(NSDE[[1]]),cex=0.5)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-(max(NSDE)/20),length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
text(bars,c(NSDE[[1]]+max(NSDE)/25),as.character(NSDE[[1]]),cex=0.5)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-(max(NSDE)/20),length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
text(bars,c(NSDE[[1]]+max(NSDE)/30),as.character(NSDE[[1]]),cex=0.5)
pdf("NumofDEGsAmongEachGroup_EdgeR.pdf",w=14,h=8)
par(mar=c(10,5,5,3))
plot_matrix<-t(as.matrix(NSDE[,2:3]))
colnames(plot_matrix)<-NULL
bars<-barplot(plot_matrix,ylim=c(0,1500),legend.text=colnames(NSDE)[2:3],args.legend=list(bty="n",horiz=TRUE),col=brewer.pal(2,"Set1"),border="white",ylab="Number of Differently Expressed Genes (DEGs)",main="Num of DEGs Among Each Group",cex.axis=1,font.axis=2)
text(bars,rep(-(max(NSDE)/20),length(bars)),labels=colnames(t(as.matrix(NSDE[,2:3]))),srt=45,adj=1,xpd=T,cex=1)
text(bars,c(NSDE[[1]]+max(NSDE)/30),as.character(NSDE[[1]]),cex=1)
dev.off()
q()
dir()
16wan1<-read.delim("16wan1-16wan2.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
Wan<-read.delim("16wan1-16wan2.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
dim(Wan)
Wan12<-read.delim("16wan1-16wan2.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
Wan23<-read.delim("16wan2-16wan3.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
names(Wan12)
names(Wan23)
Wan23<-read.delim("16wan2-16wan3.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
names(Wan23)
merge(Wan12,Wan23,by="Synonym")->Wan123
names(Wan123)
Wan12<-read.delim("16wan1-16wan2.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
Wan23<-read.delim("16wan2-16wan3.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
names(Wan123)
names(Wan12)
names(Wan23)
merge(Wan12,Wan23,by="Synonym")->Wan123
names(Wan123)
ZT12<-read.delim("zt1-zt2.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
ZT23<-read.delim("zt2-zt3.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
ZT34<-read.delim("zt3-zt4.gene.exp_multi.xls.new",header=T,sep="\t",check.names=F)
ls()
merge(Wan123,ZT12,by="Synonym")->Wan123ZT12
merge(Wan123ZT12,ZT23,by="Synonym")->Wan123ZT123
merge(Wan123ZT123,ZT34,by="Synonym")->Wan123ZT1234
Wan123ZT1234[1:4,]
names(Wan123ZT1234)
Wan123ZT1234_new<-Wan123ZT1234[,c("Synonym","16wan1","16wan2.x","16wan3","zt1","zt2.x","zt3.x","zt4")]
names(Wan123ZT1234_new)
names(Wan123ZT1234_new)[3]<-"16wan2"
names(Wan123ZT1234_new)[6]<-"zt2"
names(Wan123ZT1234_new)[7]<-"zt3"
names(Wan123ZT1234_new)
dim(Wan123ZT1234_new)
write.table(Wan123ZT1234_new,"All3445GeneCountMatrix",sep="\t",col.names=T,row.names=F,quote=F)
q()
dir()
data1 <-read.table(file="All3445GeneCountMatrix",header=T,sep="\t",check.names=F)
data <-data1[,-c("zt1","zt2")]
exp <-data[,2:(ncol(data)-1)]
names(data1)
data <-data1[,-5:6]
data <-data1[,-(5:6)]
data[1:4,]
exp <-data[,2:(ncol(data)-1)]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp) 
# cluster
pdf("SampleCluster_UnionSDEG.pdf",w=10,h=7)
hc <- hclust(as.dist(1-cor(exp_scale, method="spearman")), method="complete") 
plot(hc,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
getwd()
data1 <-read.table(file="All3445GeneCountMatrix",header=T,sep="\t",check.names=F)
data <-data1[,-5:6]
exp <-data[,2:ncol(data)]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp) 
# cluster
pdf("SampleCluster_UnionSDEG.pdf",w=10,h=7)
hc <- hclust(as.dist(1-cor(exp_scale, method="spearman")), method="complete") 
plot(hc,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
data1[1:4,]
data <-data1[,-(5,6)]
exp <-data[,2:(ncol(data)-1)]
data <-data1[,-c(5,6)]
exp <-data[,2:(ncol(data)-1)]
exp[1:4,]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp) 
# cluster
pdf("SampleCluster_UnionSDEG.pdf",w=10,h=7)
hc <- hclust(as.dist(1-cor(exp_scale, method="spearman")), method="complete") 
plot(hc,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
exp[1:4,]
write.table(exp,"All3445GeneCountMatrix_noZT1ZT2",sep="\t",col.names=T,row.names=F,quote=F)
write.table(exp,"All3445GeneCountMatrix_noZT1ZT2",sep="\t",col.names=T,row.names=T,quote=F)
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
## 热图
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
exp <-data[,2:(ncol(data)-1)]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp) 
hc_genes = agnes(exp_scale,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(exp_scale, method="spearman")), method="complete") # cluster conditions
library(cluster)
library(gplots)
library(Biobase)
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
exp <-data[,2:(ncol(data)-1)]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp) 
hc_genes = agnes(exp_scale,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(exp_scale, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file=\"all.RData\")
myheatcol = redgreen(75)[75:1]
pdf(file="All3445GeneCountMatrix_normalized.pdf", width=12,height=20, paper=\"special\")
heatmap.2(exp_scale, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="All3445GeneCountMatrix_normalized.pdf", width=12,height=20, paper=\"special\")
heatmap.2(exp_scale, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
pdf(file="All3445GeneCountMatrix_normalized.pdf", width=12,height=20, paper="special")
heatmap.2(exp_scale, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
max(exp_scale)
min(exp_scale)
max(data1)
max(data)
max(exp)
exp[1:4,]
exp_scale[1:4,]
exp_scale_log<-log(exp_scale+1,base=10)
centered_data = t(scale(t(exp_scale_log), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="All3445GeneCountMatrix_normalized.pdf", width=12,height=20, paper="special")
heatmap.2(centered_data, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
hc_genes = agnes(exp_scale_log,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(exp_scale_log, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="All3445GeneCountMatrix_normalized.pdf", width=12,height=20, paper="special")
heatmap.2(exp_scale_log, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
exp <-data[,2:(ncol(data)-1)]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp)
exp_scale_log<-log(exp_scale+1,base=10)
centered_data = t(scale(t(exp_scale_log), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="All3445GeneCountMatrix_normalized.pdf", width=12,height=20, paper="special")
heatmap.2(centered_data, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
pdf("SampleCluster_UnionSDEG.pdf",w=10,h=7)
hc_samples <- hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
data1[1:4,]
centered_data[1:4,]
data1 <-read.table(file="Union_SDEGs_noZT1ZT2_uniq_CountMatrix_genelength",header=T,sep="\t",check.names=F)
data1[1:4,]
data1 <-read.table(file="Union_SDEGs_noZT1ZT2_uniq_CountMatrix_genelength",header=T,sep="\t",check.names=F)
data1[1:4,]
data <-data1[,-c(5,6)]
exp <-data[,2:(ncol(data)-1)]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x]))
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp)
exp_scale_log<-log(exp_scale+1,base=10)
centered_data = t(scale(t(exp_scale_log), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="Union_SDEGs_noZT1ZT2.pdf", width=12,height=20, paper="special")
heatmap.2(centered_data, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
pdf("SampleCluster_UnionSDEG_noZT1ZT2.pdf",w=10,h=7)
hc_samples <- hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
library(cluster)
library(gplots)
library(Biobase)
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
exp <-data[,2:(ncol(data)-1)]
exp[1:4,]
data1[1:4,]
library(cluster)
library(gplots)
library(Biobase)
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
exp <-data[,2:(ncol(data)-1)]
rownames(exp)<-as.character(data[[1]])
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
exp <-data1[,2:(ncol(data1)-1)]
rownames(exp)<-as.character(data1[[1]])
gene_length <-data1[,ncol(data1)] 
exp[1:4,]
library(cluster)
library(gplots)
library(Biobase)
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
exp <-data1[,2:(ncol(data1)-1)]
rownames(exp)<-as.character(data1[[1]])
gene_length <-data1[,ncol(data1)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp)
exp_scale_log<-log(exp_scale+1,base=10)
centered_data = t(scale(t(exp_scale_log), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="All3445GeneCountMatrix_Heatmap.pdf", width=12,height=20, paper="special")
heatmap.2(centered_data, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
pdf("All3445GeneCountMatrix_SampleCluster.pdf",w=10,h=7)
hc_samples <- hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
LL
ls()
###########################
## 去掉zt1和zt2
################################ 热图
library(cluster)
library(gplots)
library(Biobase)
data1 <-read.table(file="All3445GeneCountMatrix_genelength",header=T,sep="\t",check.names=F)
data <-data1[,-c(5,6)]
exp <-data[,2:(ncol(data)-1)]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp)
exp_scale_log<-log(exp_scale+1,base=10)
centered_data = t(scale(t(exp_scale_log), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="All3445GeneCountMatrix_Heatmap_noZT12.pdf", width=12,height=20, paper="special")
heatmap.2(centered_data, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
pdf("All3445GeneCountMatrix_SampleCluster_noZT12.pdf",w=10,h=7)
hc_samples <- hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
###########################
## Union基因
library(cluster)
library(gplots)
library(Biobase)
data1 <-read.table(file="Union_SDEGs_uniq_CountMatrix_new",header=T,sep="\t",check.names=F)
exp <-data1[,2:(ncol(data1)-1)]
rownames(exp)<-as.character(data1[[1]])
gene_length <-data1[,ncol(data1)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x])) 
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp)
exp_scale_log<-log(exp_scale+1,base=10)
centered_data = t(scale(t(exp_scale_log), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="UnionGeneCountMatrix_Heatmap.pdf", width=12,height=20, paper="special")
heatmap.2(centered_data, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
pdf("UnionGeneCountMatrix_SampleCluster.pdf",w=10,h=7)
hc_samples <- hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
### ## 去掉zt1和zt2
data1 <-read.table(file="Union_SDEGs_noZT1ZT2_uniq_CountMatrix_genelength",header=T,sep="\t",check.names=F)
data <-data1[,-c(5,6)]
exp <-data[,2:(ncol(data)-1)]
rownames(exp)<-as.character(data[[1]])
gene_length <-data[,ncol(data)] 
exp_scale_0<-sapply(1:ncol(exp),function(x) exp[,x]*(10^9)/sum(exp[,x])) 
exp_scale<-t(sapply(1:nrow(exp),function(x) exp_scale_0[x,]/gene_length[x]))
rownames(exp_scale) <-rownames(exp) 
colnames(exp_scale) <-colnames(exp)
exp_scale_log<-log(exp_scale+1,base=10)
centered_data = t(scale(t(exp_scale_log), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman"))
gene_partition_assignments <- cu
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[7
pdf(file="UnionGeneCountMatrix_H
heatmap.2(centered_data, dendrogram="both",Rowv=as.dendrogram(hc_genes),C
dev.off()
pdf("UnionGeneCountMatrix_SampleCluster_noZT12.pdf",w=10,h=7)
hc_sa
plot(hc_samples,hang=-1,xlab="",ylab
dev.off()
hc_genes = agnes(centered_data,diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes),k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
myheatcol = redgreen(75)[75:1]
pdf(file="UnionGeneCountMatrix_Heatmap_noZT12.pdf", width=12,height=20, paper="special")
heatmap.2(centered_data, dendrogram="both",Rowv=as.dendrogram(hc_genes),Colv=as.dendrogram(hc_samples),col=myheatcol, RowSideColors=gene_colors,scale="none",density.info="none", trace="none",cexCol=1, cexRow=0.1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(5,8))
dev.off()
pdf("UnionGeneCountMatrix_SampleCluster_noZT12.pdf",w=10,h=7)
hc_samples <- hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
plot(hc_samples,hang=-1,xlab="",ylab="",main="",sub="",cex=2)
dev.off()
q()
count2fpkm_RSEMedgeR<-function(RSEM_count_path,EdgeR_FPKM_path_gene,EdgeR_FPKM_path_isoform,outdir){
### gene:
       RSEM_count_matrix<-read.delim(paste(RSEM_count_path,"all.genes.matrix.new",sep="/"),header=T,sep="\t")
       EdgeR_fpkm_matrix<-read.delim(paste(EdgeR_FPKM_path_gene,"all.genes.matrix.new.normalized.FPKM",sep="/"),header=T,sep="\t")
       names(RSEM_count_matrix)<-c("gene",paste(names(RSEM_count_matrix)[-1],".count",sep=""))
       names(EdgeR_fpkm_matrix)<-c("gene",paste(names(EdgeR_fpkm_matrix)[-1],".fpkm",sep=""))
       genes.count2fpkm<-merge(RSEM_count_matrix,EdgeR_fpkm_matrix,by="gene")
       write.table(genes.count2fpkm,paste(outdir,"genes.count2fpkm",sep="/"),sep="\t",col.names=T,row.names=F,quote=F)
### isoform:
       RSEM_count_matrix<-read.delim(paste(RSEM_count_path,"all.isoforms.matrix.new",sep="/"),header=T,sep="\t")
       EdgeR_fpkm_matrix<-read.delim(paste(EdgeR_FPKM_path_isoform,"all.isoforms.matrix.new.normalized.FPKM",sep="/"),header=T,sep="\t")
       names(RSEM_count_matrix)<-c("isoform",paste(names(RSEM_count_matrix)[-1],".count",sep=""))
       names(EdgeR_fpkm_matrix)<-c("isoform",paste(names(EdgeR_fpkm_matrix)[-1],".fpkm",sep=""))
       isoforms.count2fpkm<-merge(RSEM_count_matrix,EdgeR_fpkm_matrix,by="isoform")
       write.table(isoforms.count2fpkm,paste(outdir,"isoforms.count2fpkm",sep="/"),sep="\t",col.names=T,row.names=F,quote=F)            
}
RSEM_count_path<-"/data/users/wangyan/project/liuzhihong_24/all_RSEM/ExpressionAnalysis"
EdgeR_FPKM_path_gene<-"/data/users/wangyan/project/liuzhihong_24/all_RSEM/ExpressionAnalysis"
EdgeR_FPKM_path_isoform<-"/data/users/wangyan/project/liuzhihong_24/all_RSEM/ExpressionAnalysis"
outdir<-"/data/users/wangyan/project/liuzhihong_24/all_RSEM/ExpressionAnalysis"
count2fpkm_RSEMedgeR(RSEM_count_path,EdgeR_FPKM_path_gene,EdgeR_FPKM_path_isoform,outdir)
q()
count2fpkm_RSEMedgeR<-function(RSEM_count_path,EdgeR_FPKM_path_gene,EdgeR_FPKM_path_isoform,outdir){
### gene:
       RSEM_count_matrix<-read.delim(paste(RSEM_count_path,"all.genes.matrix.new",sep="/"),header=T,sep="\t")
       EdgeR_fpkm_matrix<-read.delim(paste(EdgeR_FPKM_path_gene,"all.genes.matrix.new.normalized.FPKM",sep="/"),header=T,sep="\t")
       names(RSEM_count_matrix)<-c("gene",paste(names(RSEM_count_matrix)[-1],".count",sep=""))
       names(EdgeR_fpkm_matrix)<-c("gene",paste(names(EdgeR_fpkm_matrix)[-1],".fpkm",sep=""))
       genes.count2fpkm<-merge(RSEM_count_matrix,EdgeR_fpkm_matrix,by="gene")
       write.table(genes.count2fpkm,paste(outdir,"genes.count2fpkm",sep="/"),sep="\t",col.names=T,row.names=F,quote=F)
       count2fpkm<-genes.count2fpkm
       save(count2fpkm,file=paste(outdir,"gene.count2fpkm.RData",sep="/"))
### isoform:
       RSEM_count_matrix<-read.delim(paste(RSEM_count_path,"all.isoforms.matrix.new",sep="/"),header=T,sep="\t")
       EdgeR_fpkm_matrix<-read.delim(paste(EdgeR_FPKM_path_isoform,"all.isoforms.matrix.new.normalized.FPKM",sep="/"),header=T,sep="\t")
       names(RSEM_count_matrix)<-c("isoform",paste(names(RSEM_count_matrix)[-1],".count",sep=""))
       names(EdgeR_fpkm_matrix)<-c("isoform",paste(names(EdgeR_fpkm_matrix)[-1],".fpkm",sep=""))
       isoforms.count2fpkm<-merge(RSEM_count_matrix,EdgeR_fpkm_matrix,by="isoform")
       write.table(isoforms.count2fpkm,paste(outdir,"isoforms.count2fpkm",sep="/"),sep="\t",col.names=T,row.names=F,quote=F)
       count2fpkm<-isoforms.count2fpkm
       save(count2fpkm,file=paste(outdir,"isoform.count2fpkm.RData",sep="/"))       
}
RSEM_count_path<-"/data/users/wangyan/project/liuzhihong_24/all_RSEM/ExpressionAnalysis"
EdgeR_FPKM_path_gene<-"/data/users/wangyan/project/liuzhihong_24/all_RSEM/ExpressionAnalysis"
EdgeR_FPKM_path_isoform<-"/data/users/wangyan/project/liuzhihong_24/all_RSEM/ExpressionAnalysis"
outdir<-"/data/users/wangyan/project/liuzhihong_24/all_RSEM/ExpressionAnalysis"
count2fpkm_RSEMedgeR(RSEM_count_path,EdgeR_FPKM_path_gene,EdgeR_FPKM_path_isoform,outdir)
q()
load("/home/Project_Archieves/RNA-seq/liu.zhi.hong/Sheep.MJ20120929422/annotation/AllAnnotation_gene.RData")
load("/home/Project_Archieves/RNA-seq/liu.zhi.hong/Sheep.MJ20120929422/annotation/AllAnnotation_transcript.RData")
path1<-"/home/Project_Archieves/RNA-seq/liu.zhi.hong/Sheep.MJ20120929422/ExpressionAnalysis/Expression_Analysis/groups_A..M"
path2<-"/home/Project_Archieves/RNA-seq/liu.zhi.hong/Sheep.MJ20120929422/Expression_Analysis"
groups<-list.files(path1)
Big_group<-"groups_A..M"
dir.create(paste(path2,Big_group,sep="/"))
for(i in 1:length(groups)){
    print(groups[i])
    diffexp_matrix_gene<-read.table(paste(path1,groups[i],"gene.diff.exp",sep="/"),sep="\t",header=T,check.names=F)
    diffexp_matrix_isoform<-read.table(paste(path1,groups[i],"isoform.diff.exp",sep="/"),sep="\t",header=T,check.names=F)
    sample<-gsub("-",".",unlist(strsplit(groups[i],"_vs_",fix=T)))
    colname_count<-paste(sample,"count",sep=".")
    colname_fpkm<-paste(sample,"fpkm",sep=".")
    diffexp_matrix_gene_new<-diffexp_matrix_gene[,c("gene",colname_count,colname_fpkm,"logFC","PValue","FDR","Direction","Sig")]
    diffexp_matrix_isoform_new<-diffexp_matrix_isoform[,c("isoform",colname_count,colname_fpkm,"logFC","PValue","FDR","Direction","Sig")]
    # add ann
    names(diffexp_matrix_gene_new)[1]<-"gene_id"
    names(diffexp_matrix_isoform_new)[1]<-"transcript_id"
    diffexp_matrix_gene_new_new<-merge(diffexp_matrix_gene_new,AllAnnotation_gene,by="gene_id")
    diffexp_matrix_isoform_new_new<-merge(diffexp_matrix_isoform_new,AllAnnotation_transcript,by="transcript_id")
    
    # plot file copy
    dir.create(paste(path2,Big_group,groups[i],sep="/"))
    pdf_from<-paste(path1,"/",groups[i],"/",groups[i],"_gene.pdf",sep="")
    pdf_to<-paste(paste(path2,"/",Big_group,"/",groups[i],"/genes.SmearPlot.pdf",sep=""))
    file.copy(pdf_from,pdf_to)
    pdf_from<-paste(path1,"/",groups[i],"/",groups[i],"_isoform.pdf",sep="")
    pdf_to<-paste(paste(path2,"/",Big_group,"/",groups[i],"/isoforms.SmearPlot.pdf",sep=""))
    file.copy(pdf_from,pdf_to)
    write.table(diffexp_matrix_gene_new_new,paste(paste(path2,"/",Big_group,"/",groups[i],"/genes.diff_exp_ann_result",sep="")),sep="\t",col.names=T,row.names=F,quote=F)
    write.table(diffexp_matrix_isoform_new_new,paste(paste(path2,"/",Big_group,"/",groups[i],"/isoforms.diff_exp_ann_result",sep="")),sep="\t",col.names=T,row.names=F,quote=F)
     
}
load("/home/Project_Archieves/RNA-seq/liu.zhi.hong/Sheep.MJ20120929422/annotation/AllAnnotation_gene.RData")
load("/home/Project_Archieves/RNA-seq/liu.zhi.hong/Sheep.MJ20120929422/annotation/AllAnnotation_transcript.RData")
path1<-"/home/Project_Archieves/RNA-seq/liu.zhi.hong/Sheep.MJ20120929422/ExpressionAnalysis/Expression_Analysis/groups_A..M"
path2<-"/home/Project_Archieves/RNA-seq/liu.zhi.hong/Sheep.MJ20120929422/Expression_Analysis"
groups<-list.files(path1)
Big_group<-"groups_A..M"
dir.create(paste(path2,Big_group,sep="/"))
for(i in 1:length(groups)){
    print(groups[i])
    diffexp_matrix_gene<-read.table(paste(path1,groups[i],"gene.diff.exp",sep="/"),sep="\t",header=T,check.names=F)
    sample<-gsub("-",".",unlist(strsplit(groups[i],"_vs_",fix=T)))
    colname_count<-paste(sample,"count",sep=".")
    colname_fpkm<-paste(sample,"fpkm",sep=".")
    diffexp_matrix_gene_new<-diffexp_matrix_gene[,c("gene",colname_count,colname_fpkm,"logFC","PValue","FDR","Direction","Sig")]
    # add ann
    names(diffexp_matrix_gene_new)[1]<-"gene_id"
    names(diffexp_matrix_isoform_new)[1]<-"transcript_id"
    diffexp_matrix_gene_new_new<-merge(diffexp_matrix_gene_new,AllAnnotation_gene,by="gene_id")
    
    # plot file copy
    dir.create(paste(path2,Big_group,groups[i],sep="/"))
    pdf_from<-paste(path1,"/",groups[i],"/",groups[i],"_gene.pdf",sep="")
    pdf_to<-paste(paste(path2,"/",Big_group,"/",groups[i],"/genes.SmearPlot.pdf",sep=""))
    file.copy(pdf_from,pdf_to)
    write.table(diffexp_matrix_gene_new_new,paste(paste(path2,"/",Big_group,"/",groups[i],"/genes.diff_exp_ann_result",sep="")),sep="\t",col.names=T,row.names=F,quote=F)
     
}
for(i in 1:length(groups)){
    print(groups[i])
    diffexp_matrix_gene<-read.table(paste(path1,groups[i],"gene.diff.exp",sep="/"),sep="\t",header=T,check.names=F)
    sample<-gsub("-",".",unlist(strsplit(groups[i],"_vs_",fix=T)))
    colname_count<-paste(sample,"count",sep=".")
    colname_fpkm<-paste(sample,"fpkm",sep=".")
    diffexp_matrix_gene_new<-diffexp_matrix_gene[,c("gene",colname_count,colname_fpkm,"logFC","PValue","FDR","Direction","Sig")]
    # add ann
    names(diffexp_matrix_gene_new)[1]<-"gene_id"
    diffexp_matrix_gene_new_new<-merge(diffexp_matrix_gene_new,AllAnnotation_gene,by="gene_id")
    
    # plot file copy
    dir.create(paste(path2,Big_group,groups[i],sep="/"))
    pdf_from<-paste(path1,"/",groups[i],"/",groups[i],"_gene.pdf",sep="")
    pdf_to<-paste(paste(path2,"/",Big_group,"/",groups[i],"/genes.SmearPlot.pdf",sep=""))
    file.copy(pdf_from,pdf_to)
    write.table(diffexp_matrix_gene_new_new,paste(paste(path2,"/",Big_group,"/",groups[i],"/genes.diff_exp_ann_result",sep="")),sep="\t",col.names=T,row.names=F,quote=F)
     
}
q()
DNA_path<-"/home/Project_Archieves/RNA-seq/liu.zhi.hong/Erinaceus_europaeus.MJ20120929954/annotation/DNA.ann"
Protein_path<-"/home/Project_Archieves/RNA-seq/liu.zhi.hong/Erinaceus_europaeus.MJ20120929954/annotation/Protein.ann"
Out_path<-"merged_annotation"
gene2isoform<-"no mapping info"
q()
ls()
gene2isoform
gene2isoform<-"/home/Project_Archieves/RNA-seq/liu.zhi.hong/Erinaceus_europaeus.MJ20120929954/cuffdiff_out_dir/gene2isoforms"
gene2transcript<-unique(read.delim(gene2isoform,sep="\t",header=F));names(gene2transcript)<-c("gene_id","transcript_id")
DNA_path
Protein_path
ORF<-rep("no ORF",dim(DNAAnnotation)[1])
DNAAnnotation_new<-cbind(DNAAnnotation,ORF)[,names(ProteinAnnotation)]
AllAnnotation<-rbind(DNAAnnotation_new,ProteinAnnotation)[,-13]
newnames<-c("transcript_id","Qeury_length","ORF","NR_tophit_name","NR_tophit_description","NR_tophit_similarity","GOs","String_tophit_name","String_tophit_description","String_tophit_similarity","COG"                   ,"KOG","KO.Gene_ID","KEGG_GENE_NAME")
names(AllAnnotation)<-newnames
AllAnnotation_transcript<-unique(AllAnnotation[order(as.character(AllAnnotation[[1]]),decreasing=F),])
write.table(AllAnnotation_transcript,paste(Out_dir,"AllAnnotation_transcript",sep="/"),sep="\t",col.names=T,row.names=F,quote=F)
save(AllAnnotation_transcript,file=paste(Out_dir,"AllAnnotation_transcript.RData",sep="/"))
print("transcript Done !!!")
if(DNA_path=="no dna"){DNAAnnotation<-data.frame()}else{DNAAnnotation<-read.delim(DNA_path,sep="\t",header=TRUE)}
if(Protein_path=="no pep"){ProteinAnnotation<-data.frame()}else{ProteinAnnotation<-read.delim(Protein_path,sep="\t",header=TRUE)}
current_path<-getwd()
Out_dir<-paste(current_path,Out_path,sep="/")
dir.create(Out_dir)
names(DNAAnnotation)
#####\ isoform annotation
ORF<-rep("no ORF",dim(DNAAnnotation)[1])
DNAAnnotation_new<-cbind(DNAAnnotation,ORF)[,names(ProteinAnnotation)]
AllAnnotation<-rbind(DNAAnnotation_new,ProteinAnnotation)[,-13]
newnames<-c("transcript_id","Qeury_length","ORF","NR_tophit_name","NR_tophit_description","NR_tophit_similarity","GOs","String_tophit_name","String_tophit_description","String_tophit_similarity","COG"                   ,"KOG","KO.Gene_ID","KEGG_GENE_NAME")
names(AllAnnotation)<-newnames
AllAnnotation_transcript<-unique(AllAnnotation[order(as.character(AllAnnotation[[1]]),decreasing=F),])
write.table(AllAnnotation_transcript,paste(Out_dir,"AllAnnotation_transcript",sep="/"),sep="\t",col.names=T,row.names=F,quote=F)
save(AllAnnotation_transcript,file=paste(Out_dir,"AllAnnotation_transcript.RData",sep="/"))
print("transcript Done !!!")
gene_id<-vector()
transcript_num<-vector()
AllAnnotation_transcript[,1]
AllAnnotation_transcript[,1:10]
AllAnnotation_transcript[1:4,1:10]
AllAnnotation_transcript[1:4,1:3]
gene2transcript[1:4,]
gene2transcript<-unique(read.delim(gene2isoform,sep="\t",header=T));names(gene2transcript)<-c("gene_id","transcript_id")
gene2transcript[1:4,]
names(AllAnnotation_transcript)
AllAnnotation_transcript_addgene<-as.data.frame(merge(gene2transcript,AllAnnotation_transcript,by="transcript_id"))
AllAnnotation_transcript_addgene[1:4,1:3]
AllAnnotation_transcript_addgene<-as.data.frame(merge(gene2transcript,AllAnnotation_transcript,by="transcript_id"))[,c("gene_id",names(AllAnnotation_transcript))]
AllAnnotation_transcript_addgene[1:4,1:3]
row.names(AllAnnotation_transcript_addgene)<-as.character(AllAnnotation_transcript_addgene[["transcript_id"]])
gene_ids<-unique(as.character(AllAnnotation_transcript_addgene[[1]]))
nnn<-table(as.character(AllAnnotation_transcript_addgene[["gene_id"]]))
MultipleLine<-names(nnn)[which(nnn>=2)]
TobeConverted<-AllAnnotation_transcript_addgene[AllAnnotation_transcript_addgene[[1]]%in%MultipleLine,]
noNeedConverted<-AllAnnotation_transcript_addgene[(AllAnnotation_transcript_addgene[[1]]%in%MultipleLine)==FALSE,]
transcript_num<-rep(1,dim(noNeedConverted)[1])
noNeedConverted_new<-cbind(noNeedConverted[,1],transcript_num,noNeedConverted[,2:15])
names(noNeedConverted_new)[1]<-"gene_id"
noNeedConverted_new[1:4,1:3]
Genes<-as.character(TobeConverted[["gene_id"]])
TobeConverted_new<-data.frame()
for(i in 1:length(Genes)){
    TobeConverted_new[i,1]<-Genes[i]
    Line<-which(TobeConverted[["gene_id"]]==Genes[i])
    TobeConverted_new[i,2]<-length(Line)
    TobeConverted_new[i,3]<-paste(as.character(TobeConverted[Line,"transcript_id"]),collapse=";")
    TobeConverted_new[i,4]<-paste(as.character(TobeConverted[Line,"Qeury_length"]),collapse=";")
    TobeConverted_new[i,5]<-paste(as.character(TobeConverted[Line,"ORF"]),collapse=";")
    TobeConverted_new[i,6]<-paste(unique(as.character(TobeConverted[Line,"NR_tophit_name"])),collapse=";")
    TobeConverted_new[i,7]<-paste(unique(as.character(TobeConverted[Line,"NR_tophit_description"])),collapse=";")
    TobeConverted_new[i,8]<-paste(unique(as.character(TobeConverted[Line,"NR_tophit_similarity"])),collapse=";")
    TobeConverted_new[i,9]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"GOs"]),collapse=";"),";",fix=T))
),collapse=";")
    TobeConverted_new[i,10]<-paste(unique(as.character(TobeConverted[Line,"String_tophit_name"])),collapse=";")
    TobeConverted_new[i,11]<-paste(unique(as.character(TobeConverted[Line,"String_tophit_description"])),collapse=";")
    TobeConverted_new[i,12]<-paste(unique(as.character(TobeConverted[Line,"String_tophit_similarity"])),collapse=";")
    TobeConverted_new[i,13]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"COG"]),collapse=";"),";",fix=T))
),collapse=";")
    TobeConverted_new[i,14]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"KOG"]),collapse=";"),";",fix=T))
),collapse=";")
    TobeConverted_new[i,15]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"KO.Gene_ID"]),collapse=";"),";",fix=T))
),collapse=";")
    TobeConverted_new[i,16]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"KEGG_GENE_NAME"]),collapse=";"),";",fix=T))
),collapse=";")
}
length(Genes)
dim(TobeConverted)
dim(AllAnnotation_transcript_addgene)
ls（）
ls()
dim(gene2transcript)
dim(AllAnnotation_transcript)
AllAnnotation_transcript_addgene<-as.data.frame(merge(gene2transcript,AllAnnotation_transcript,by="transcript_id",all=T))[,c("gene_id",names(AllAnnotation_transcript))]
row.names(AllAnnotation_transcript_addgene)<-as.character(AllAnnotation_transcript_addgene[["transcript_id"]])
gene_ids<-unique(as.character(AllAnnotation_transcript_addgene[[1]]))
nnn<-table(as.character(AllAnnotation_transcript_addgene[["gene_id"]]))
MultipleLine<-names(nnn)[which(nnn>=2)]
TobeConverted<-AllAnnotation_transcript_addgene[AllAnnotation_transcript_addgene[[1]]%in%MultipleLine,]
noNeedConverted<-AllAnnotation_transcript_addgene[(AllAnnotation_transcript_addgene[[1]]%in%MultipleLine)==FALSE,]
transcript_num<-rep(1,dim(noNeedConverted)[1])
noNeedConverted_new<-cbind(noNeedConverted[,1],transcript_num,noNeedConverted[,2:15])
names(noNeedConverted_new)[1]<-"gene_id"
dim(TobeConverted)
dim(AllAnnotation_transcript_addgene)
dim(gene2transcript)
gene2transcript[1:4,]
length(unique(gene2transcript[[1]]))
length(unique(gene2transcript[[2]]))
length(unique(as.character(gene2transcript[[2]])))
length(unique(as.character(gene2transcript[[1]])))
noNeedConverted
TobeConverted
Genes<-as.character(TobeConverted[["gene_id"]])
TobeConverted_new<-data.frame()
for(i in 1:length(Genes)){
    TobeConverted_new[i,1]<-Genes[i]
    Line<-which(TobeConverted[["gene_id"]]==Genes[i])
    TobeConverted_new[i,2]<-length(Line)
    TobeConverted_new[i,3]<-paste(as.character(TobeConverted[Line,"transcript_id"]),collapse=";")
    TobeConverted_new[i,4]<-paste(as.character(TobeConverted[Line,"Qeury_length"]),collapse=";")
    TobeConverted_new[i,5]<-paste(as.character(TobeConverted[Line,"ORF"]),collapse=";")
    TobeConverted_new[i,6]<-paste(unique(as.character(TobeConverted[Line,"NR_tophit_name"])),collapse=";")
    TobeConverted_new[i,7]<-paste(unique(as.character(TobeConverted[Line,"NR_tophit_description"])),collapse=";")
    TobeConverted_new[i,8]<-paste(unique(as.character(TobeConverted[Line,"NR_tophit_similarity"])),collapse=";")
    TobeConverted_new[i,9]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"GOs"]),collapse=";"),";",fix=T))
),collapse=";")
    TobeConverted_new[i,10]<-paste(unique(as.character(TobeConverted[Line,"String_tophit_name"])),collapse=";")
    TobeConverted_new[i,11]<-paste(unique(as.character(TobeConverted[Line,"String_tophit_description"])),collapse=";")
    TobeConverted_new[i,12]<-paste(unique(as.character(TobeConverted[Line,"String_tophit_similarity"])),collapse=";")
    TobeConverted_new[i,13]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"COG"]),collapse=";"),";",fix=T))
),collapse=";")
    TobeConverted_new[i,14]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"KOG"]),collapse=";"),";",fix=T))
),collapse=";")
    TobeConverted_new[i,15]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"KO.Gene_ID"]),collapse=";"),";",fix=T))
),collapse=";")
    TobeConverted_new[i,16]<-paste(unique(unlist(strsplit(paste(as.character(TobeConverted[Line,"KEGG_GENE_NAME"]),collapse=";"),";",fix=T))
),collapse=";")
}
names(TobeConverted_new)<-c(names(TobeConverted)[1],"transcript_num",names(TobeConverted)[2:15])
AllAnnotation_gene<-unique(rbind(noNeedConverted_new,TobeConverted_new))
write.table(AllAnnotation_gene,paste(Out_dir,"AllAnnotation_gene",sep="/"),sep="\t",col.names=T,row.names=F,quote=F)
save(AllAnnotation_gene,file=paste(Out_dir,"AllAnnotation_gene.RData",sep="/"))
print("gene Done !!!")
q()
names(resultMatrix_isoform)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_isoform,isoform.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","isoforms.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEisoform1<-as.data.frame(DEisoform)
names(DEisoform1)<-names(isoform.fpkm.matrix)[1]
DEisoform.fpkm.matrix<-merge(DEisoform1,isoform.fpkm.matrix,by=names(isoform.fpkm.matrix)[1])
write.table(DEisoform.fpkm.matrix,paste(process.file,"/",filename,"DEisoform.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        ####### sig.detail.matrix generation
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
}
}
group_name
gene.name.list<-c("GeneName",as.character(gene.fpkm.matrix[[1]]),dim(gene.fpkm.matrix)[1])
isoform.name.list<-c("IsoformName",as.character(isoform.fpkm.matrix[[1]]),dim(isoform.fpkm.matrix)[1])
sig.detail.gene<-gene.name.list
sig.detail.isoform<-isoform.name.list
        group1<-group_name[i]
for(j in (i+1):length(group_name)){
   group2<-group_name[j]
   filename<-paste(group1,group2,sep="_vs_")
   output_dir_name<-paste(groups_out_dir,filename,sep="/")
   dir.create(output_dir_name)
######  gene
resultMatrix_gene<-find_genes(group1,group2,gene.diff)
sig_resultMatrix_gene<-find_genes(group1,group2,gene.diff.sig)
DEgene<-as.character(sig_resultMatrix_gene[,"test_id"]) 
names(resultMatrix_gene)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_gene,gene.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEgene1<-as.data.frame(DEgene)
names(DEgene1)<-names(gene.fpkm.matrix)[1]
DEgene.fpkm.matrix<-merge(DEgene1,gene.fpkm.matrix,by=names(gene.fpkm.matrix)[1])
write.table(DEgene.fpkm.matrix,paste(process.file,"/",filename,"_DEgene.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        
######  isoform
resultMatrix_isoform<-find_genes(group1,group2,isoform.diff)
sig_resultMatrix_isoform<-find_genes(group1,group2,isoform.diff.sig)
DEisoform<-as.character(sig_resultMatrix_isoform[,"test_id"])
names(resultMatrix_isoform)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_isoform,isoform.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","isoforms.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEisoform1<-as.data.frame(DEisoform)
names(DEisoform1)<-names(isoform.fpkm.matrix)[1]
DEisoform.fpkm.matrix<-merge(DEisoform1,isoform.fpkm.matrix,by=names(isoform.fpkm.matrix)[1])
write.table(DEisoform.fpkm.matrix,paste(process.file,"/",filename,"DEisoform.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        ####### sig.detail.matrix generation
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
}
for(i in 1:(length(group_name)-1)){
        group1<-group_name[i]
for(j in (i+1):length(group_name)){
   group2<-group_name[j]
   filename<-paste(group1,group2,sep="_vs_")
   output_dir_name<-paste(groups_out_dir,filename,sep="/")
   dir.create(output_dir_name)
######  gene
resultMatrix_gene<-find_genes(group1,group2,gene.diff)
sig_resultMatrix_gene<-find_genes(group1,group2,gene.diff.sig)
DEgene<-as.character(sig_resultMatrix_gene[,"test_id"]) 
names(resultMatrix_gene)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_gene,gene.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEgene1<-as.data.frame(DEgene)
names(DEgene1)<-names(gene.fpkm.matrix)[1]
DEgene.fpkm.matrix<-merge(DEgene1,gene.fpkm.matrix,by=names(gene.fpkm.matrix)[1])
write.table(DEgene.fpkm.matrix,paste(process.file,"/",filename,"_DEgene.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        
######  isoform
resultMatrix_isoform<-find_genes(group1,group2,isoform.diff)
sig_resultMatrix_isoform<-find_genes(group1,group2,isoform.diff.sig)
DEisoform<-as.character(sig_resultMatrix_isoform[,"test_id"])
names(resultMatrix_isoform)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_isoform,isoform.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","isoforms.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEisoform1<-as.data.frame(DEisoform)
names(DEisoform1)<-names(isoform.fpkm.matrix)[1]
DEisoform.fpkm.matrix<-merge(DEisoform1,isoform.fpkm.matrix,by=names(isoform.fpkm.matrix)[1])
write.table(DEisoform.fpkm.matrix,paste(process.file,"/",filename,"DEisoform.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        ####### sig.detail.matrix generation
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
print(i)
}
}
        group1<-group_name[i]
for(j in (i+1):length(group_name)){
   group2<-group_name[j]
   filename<-paste(group1,group2,sep="_vs_")
   output_dir_name<-paste(groups_out_dir,filename,sep="/")
   dir.create(output_dir_name)
######  gene
resultMatrix_gene<-find_genes(group1,group2,gene.diff)
sig_resultMatrix_gene<-find_genes(group1,group2,gene.diff.sig)
DEgene<-as.character(sig_resultMatrix_gene[,"test_id"]) 
names(resultMatrix_gene)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_gene,gene.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEgene1<-as.data.frame(DEgene)
names(DEgene1)<-names(gene.fpkm.matrix)[1]
DEgene.fpkm.matrix<-merge(DEgene1,gene.fpkm.matrix,by=names(gene.fpkm.matrix)[1])
write.table(DEgene.fpkm.matrix,paste(process.file,"/",filename,"_DEgene.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        
######  isoform
resultMatrix_isoform<-find_genes(group1,group2,isoform.diff)
sig_resultMatrix_isoform<-find_genes(group1,group2,isoform.diff.sig)
DEisoform<-as.character(sig_resultMatrix_isoform[,"test_id"])
names(resultMatrix_isoform)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_isoform,isoform.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","isoforms.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEisoform1<-as.data.frame(DEisoform)
names(DEisoform1)<-names(isoform.fpkm.matrix)[1]
DEisoform.fpkm.matrix<-merge(DEisoform1,isoform.fpkm.matrix,by=names(isoform.fpkm.matrix)[1])
write.table(DEisoform.fpkm.matrix,paste(process.file,"/",filename,"DEisoform.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        ####### sig.detail.matrix generation
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
print(i)
}
gene.name.list<-c("GeneName",as.character(gene.fpkm.matrix[[1]]),dim(gene.fpkm.matrix)[1])
isoform.name.list<-c("IsoformName",as.character(isoform.fpkm.matrix[[1]]),dim(isoform.fpkm.matrix)[1])
sig.detail.gene<-gene.name.list
sig.detail.isoform<-isoform.name.list
length(sig.detail.gene)
length(tmp_gene)
length(sig.detail.isoform)
length(tmp_isoform)
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
gene.name.list<-c("GeneName",as.character(gene.fpkm.matrix[[1]]),dim(gene.fpkm.matrix)[1])
isoform.name.list<-c("IsoformName",as.character(isoform.fpkm.matrix[[1]]),dim(isoform.fpkm.matrix)[1])
sig.detail.gene<-gene.name.list
sig.detail.isoform<-isoform.name.list
group1<-group_name[i]
for(j in (i+1):length(group_name)){
   group2<-group_name[j]
   filename<-paste(group1,group2,sep="_vs_")
   output_dir_name<-paste(groups_out_dir,filename,sep="/")
   dir.create(output_dir_name)
######  gene
resultMatrix_gene<-find_genes(group1,group2,gene.diff)
sig_resultMatrix_gene<-find_genes(group1,group2,gene.diff.sig)
DEgene<-as.character(sig_resultMatrix_gene[,"test_id"]) 
names(resultMatrix_gene)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_gene,gene.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEgene1<-as.data.frame(DEgene)
names(DEgene1)<-names(gene.fpkm.matrix)[1]
DEgene.fpkm.matrix<-merge(DEgene1,gene.fpkm.matrix,by=names(gene.fpkm.matrix)[1])
write.table(DEgene.fpkm.matrix,paste(process.file,"/",filename,"_DEgene.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        
######  isoform
resultMatrix_isoform<-find_genes(group1,group2,isoform.diff)
sig_resultMatrix_isoform<-find_genes(group1,group2,isoform.diff.sig)
DEisoform<-as.character(sig_resultMatrix_isoform[,"test_id"])
names(resultMatrix_isoform)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_isoform,isoform.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","isoforms.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEisoform1<-as.data.frame(DEisoform)
names(DEisoform1)<-names(isoform.fpkm.matrix)[1]
DEisoform.fpkm.matrix<-merge(DEisoform1,isoform.fpkm.matrix,by=names(isoform.fpkm.matrix)[1])
write.table(DEisoform.fpkm.matrix,paste(process.file,"/",filename,"DEisoform.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        ####### sig.detail.matrix generation
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
print(i)
}
dim(sig.detail.isoform)
   group2<-group_name[j]
   filename<-paste(group1,group2,sep="_vs_")
   output_dir_name<-paste(groups_out_dir,filename,sep="/")
   dir.create(output_dir_name)
######  gene
resultMatrix_gene<-find_genes(group1,group2,gene.diff)
sig_resultMatrix_gene<-find_genes(group1,group2,gene.diff.sig)
DEgene<-as.character(sig_resultMatrix_gene[,"test_id"]) 
names(resultMatrix_gene)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_gene,gene.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEgene1<-as.data.frame(DEgene)
names(DEgene1)<-names(gene.fpkm.matrix)[1]
DEgene.fpkm.matrix<-merge(DEgene1,gene.fpkm.matrix,by=names(gene.fpkm.matrix)[1])
write.table(DEgene.fpkm.matrix,paste(process.file,"/",filename,"_DEgene.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        
######  isoform
resultMatrix_isoform<-find_genes(group1,group2,isoform.diff)
sig_resultMatrix_isoform<-find_genes(group1,group2,isoform.diff.sig)
DEisoform<-as.character(sig_resultMatrix_isoform[,"test_id"])
names(resultMatrix_isoform)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_isoform,isoform.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction,DetailMatrix[,10:12])
write.table(DetailMatrix.new,paste(output_dir_name,"/","isoforms.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEisoform1<-as.data.frame(DEisoform)
names(DEisoform1)<-names(isoform.fpkm.matrix)[1]
DEisoform.fpkm.matrix<-merge(DEisoform1,isoform.fpkm.matrix,by=names(isoform.fpkm.matrix)[1])
write.table(DEisoform.fpkm.matrix,paste(process.file,"/",filename,"DEisoform.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        ####### sig.detail.matrix generation
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
print(i)
dim(DetailMatrix)
dim(Direction)
length(Direction)
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction)
DetailMatrix.new.new<-cbind(DetailMatrix.new,DetailMatrix[,10:12])
DetailMatrix.new.new[1:4,]
gene.name.list<-c("GeneName",as.character(gene.fpkm.matrix[[1]]),dim(gene.fpkm.matrix)[1])
isoform.name.list<-c("IsoformName",as.character(isoform.fpkm.matrix[[1]]),dim(isoform.fpkm.matrix)[1])
sig.detail.gene<-gene.name.list
sig.detail.isoform<-isoform.name.list
for(i in 1:(length(group_name)-1)){
        group1<-group_name[i]
for(j in (i+1):length(group_name)){
   group2<-group_name[j]
   filename<-paste(group1,group2,sep="_vs_")
   output_dir_name<-paste(groups_out_dir,filename,sep="/")
   dir.create(output_dir_name)
######  gene
resultMatrix_gene<-find_genes(group1,group2,gene.diff)
sig_resultMatrix_gene<-find_genes(group1,group2,gene.diff.sig)
DEgene<-as.character(sig_resultMatrix_gene[,"test_id"]) 
names(resultMatrix_gene)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_gene,gene.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction)
DetailMatrix.new.new<-cbind(DetailMatrix.new,DetailMatrix[,10:12])
write.table(DetailMatrix.new.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEgene1<-as.data.frame(DEgene)
names(DEgene1)<-names(gene.fpkm.matrix)[1]
DEgene.fpkm.matrix<-merge(DEgene1,gene.fpkm.matrix,by=names(gene.fpkm.matrix)[1])
write.table(DEgene.fpkm.matrix,paste(process.file,"/",filename,"_DEgene.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        
######  isoform
resultMatrix_isoform<-find_genes(group1,group2,isoform.diff)
sig_resultMatrix_isoform<-find_genes(group1,group2,isoform.diff.sig)
DEisoform<-as.character(sig_resultMatrix_isoform[,"test_id"])
names(resultMatrix_isoform)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_isoform,isoform.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction)
DetailMatrix.new.new<-cbind(DetailMatrix.new,DetailMatrix[,10:12])
write.table(DetailMatrix.new.new,paste(output_dir_name,"/","isoforms.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEisoform1<-as.data.frame(DEisoform)
names(DEisoform1)<-names(isoform.fpkm.matrix)[1]
DEisoform.fpkm.matrix<-merge(DEisoform1,isoform.fpkm.matrix,by=names(isoform.fpkm.matrix)[1])
write.table(DEisoform.fpkm.matrix,paste(process.file,"/",filename,"DEisoform.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        ####### sig.detail.matrix generation
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
print(i)
}
}
gene.name.list<-c("GeneName",as.character(gene.fpkm.matrix[[1]]),dim(gene.fpkm.matrix)[1])
isoform.name.list<-c("IsoformName",as.character(isoform.fpkm.matrix[[1]]),dim(isoform.fpkm.matrix)[1])
sig.detail.gene<-gene.name.list
sig.detail.isoform<-isoform.name.list
group1<-group_name[i]
 group2<-group_name[j]
   filename<-paste(group1,group2,sep="_vs_")
   output_dir_name<-paste(groups_out_dir,filename,sep="/")
   dir.create(output_dir_name)
######  gene
resultMatrix_gene<-find_genes(group1,group2,gene.diff)
sig_resultMatrix_gene<-find_genes(group1,group2,gene.diff.sig)
DEgene<-as.character(sig_resultMatrix_gene[,"test_id"]) 
names(resultMatrix_gene)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_gene,gene.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction)
DetailMatrix.new.new<-cbind(DetailMatrix.new,DetailMatrix[,10:12])
write.table(DetailMatrix.new.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
dim(DetailMatrix.new)
dim(DetailMatrix[,10:12])
dim(DetailMatrix)
dim(resultMatrix_gene)
length(Direction)
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
length(log2FC)
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
length(log2FC)
length(Direction)
p
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
Direction<-vector()
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
length(Direction)
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction)
DetailMatrix.new.new<-cbind(DetailMatrix.new,DetailMatrix[,10:12])
write.table(DetailMatrix.new.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
for(i in 1:(length(group_name)-1)){
        group1<-group_name[i]
for(j in (i+1):length(group_name)){
   group2<-group_name[j]
   filename<-paste(group1,group2,sep="_vs_")
   output_dir_name<-paste(groups_out_dir,filename,sep="/")
   dir.create(output_dir_name)
######  gene
resultMatrix_gene<-find_genes(group1,group2,gene.diff)
sig_resultMatrix_gene<-find_genes(group1,group2,gene.diff.sig)
DEgene<-as.character(sig_resultMatrix_gene[,"test_id"]) 
names(resultMatrix_gene)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_gene,gene.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
Direction<-vector()
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction)
DetailMatrix.new.new<-cbind(DetailMatrix.new,DetailMatrix[,10:12])
write.table(DetailMatrix.new.new,paste(output_dir_name,"/","genes.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEgene1<-as.data.frame(DEgene)
names(DEgene1)<-names(gene.fpkm.matrix)[1]
DEgene.fpkm.matrix<-merge(DEgene1,gene.fpkm.matrix,by=names(gene.fpkm.matrix)[1])
write.table(DEgene.fpkm.matrix,paste(process.file,"/",filename,"_DEgene.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        
######  isoform
resultMatrix_isoform<-find_genes(group1,group2,isoform.diff)
sig_resultMatrix_isoform<-find_genes(group1,group2,isoform.diff.sig)
DEisoform<-as.character(sig_resultMatrix_isoform[,"test_id"])
names(resultMatrix_isoform)[1]<-"tracking_id"      
DetailMatrix<-merge(resultMatrix_isoform,isoform.count.matrix,by="tracking_id")[,c("tracking_id","gene_id","gene","locus",Count_col,"value_1","value_2","log2.fold_change.","p_value","q_value","significant")]
names(DetailMatrix)[7:8]<-FPKM_col
# direction
log2FC<-as.numeric(DetailMatrix[,"log2.fold_change."])
Direction<-vector()
for(p in 1:length(log2FC)){
    if(log2FC[p]<0){Direction[p]<-"down"}
    if(log2FC[p]>0){Direction[p]<-"up"}
    if(log2FC[p]==0){Direction[p]<-"-"}
}
DetailMatrix.new<-cbind(DetailMatrix[,1:8],Direction)
DetailMatrix.new.new<-cbind(DetailMatrix.new,DetailMatrix[,10:12])
write.table(DetailMatrix.new.new,paste(output_dir_name,"/","isoforms.diff_exp_result",sep=""),sep="\t",col.names=T,row.names=F,quote=F) 
### FPKM_matrix:
DEisoform1<-as.data.frame(DEisoform)
names(DEisoform1)<-names(isoform.fpkm.matrix)[1]
DEisoform.fpkm.matrix<-merge(DEisoform1,isoform.fpkm.matrix,by=names(isoform.fpkm.matrix)[1])
write.table(DEisoform.fpkm.matrix,paste(process.file,"/",filename,"DEisoform.fpkm_matrix",sep=""),sep="t",col.names=T,row.names=F,quote=F)
        ####### sig.detail.matrix generation
tmp_gene<-c(filename,as.character(resultMatrix_gene[["significant"]]),length(DEgene))
sig.detail.gene<-as.data.frame(cbind(sig.detail.gene,tmp_gene))
        tmp_isoform<-c(filename,as.character(resultMatrix_isoform[["significant"]]),length(DEisoform))
sig.detail.isoform<-as.data.frame(cbind(sig.detail.isoform,tmp_isoform))
print(i)
}
}
q()
vs_matrix<-read.delim("CB_vs_HF_DEgene.fpkm_matrix",sep="\t",header=T,check.names=F)
union_matrix<-read.delim("union.DEGene.fpkm.matrix",sep="\t",header=T,check.names=F)
inter_matrix<-read.delim("1_inter.DEisoform.fpkm.matrix",sep="\t",header=T,check.names=F)
vs_isoforms<-as.character(vs_matrix[[1])
vs_isoforms<-as.character(vs_matrix[[1]]
)
union_isoforms<-as.character(union_matrix[[1]])
inter_isoforms<-as.character(inter_matrix[[1]])
length(vs_isoforms)
length(union_isoforms)
length(inter_isoforms)
intersect(vs_isoforms,union_isoforms)
intersect(vs_isoforms,inter_isoforms)
ls()
inter_isoforms
intersect(vs_isoforms,inter_isoforms)
vs_isoforms
inter_matrix<-read.delim("1_inter.DEgene.fpkm.matrix",sep="\t",header=T,check.names=F)
inter_isoforms<-as.character(inter_matrix[[1]])
intersect(vs_isoforms,inter_isoforms)
intersect(vs_isoforms,union_isoforms)
intersect(inter_isoforms,union_isoforms)
vs_matrix<-read.delim("CB_vs_HF_DEisoform.fpkm_matrix",sep="\t",header=T,check.names=F)
inter_matrix<-read.delim("1_inter.DEisoform.fpkm.matrix",sep="\t",header=T,check.names=F)
union_matrix<-read.delim("union.DEIsoform.fpkm.matrix",sep="\t",header=T,check.names=F)
union_isoforms<-as.character(union_matrix[[1]])
inter_isoforms<-as.character(inter_matrix[[1]])
vs_isoforms<-as.character(vs_matrix[[1]]
)
intersect(union_isoforms,inter_isoforms)
q()
SSR_result<-read.delim("SSR_Barplot.txt",header=F,sep="\t")
SSR_result
SSR_result<-read.delim("SSR_Barplot.txt",header=F,sep="\t")
SSR_result
SSR_result<-read.delim("SSR_Barplot.txt",header=F,sep="\t")
SSR_result
numbers<-as.numeric(SSR_result[[2]])
jump<-max(numbers)/10
ylims=c(0,(max(numbers)+jump))
jump
ylims
numbers
bar<-barplot(numbers,col="darkred",horiz=TRUE,ylim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
xlims=c(0,(max(numbers)+jump))
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
bar
labels<-as.character(SSR_result[[1]])
labl~es
labels
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
par(mar=c(4,20,4,4))
numbers<-as.numeric(SSR_result[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
par(mar=c(4,15,4,4))
numbers<-as.numeric(SSR_result[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=0.8,font=1)
jump
numbers
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.8,font=1)
SSR_result
?barplot
SSR_result
SSR_result_new<-SSR_result[nrow(SSR_result):1]
SSR_result_new<-SSR_result[nrow(SSR_result):1,]
SSR_result_new
SSR_result
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.8,font=1)
box()
?axis
xlims
axis_s<-seq(xlims[1],xlims[2],length=1000)
axis_s
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis_s
axis_s<-seq(xlims[1],xlims[2],by=500)
axis(side=3,pos=c(0,0),at=axis_s,labels=as.character(axis_s),cex.axis=2,las=1) 
axis(side=3,pos=c(0,0),at=axis_s,labels=as.character(axis_s),cex.axis=0.8) 
?axis
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.8) 
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.8,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.8) 
minor.tick(nx=2,ny=2,tick.ratio=0.5)
library(Hmisc)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.8,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.8) 
minor.tick(nx=2,ny=2,tick.ratio=0.5)
?minor.tick
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.8,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.8) 
minor.tick(nx=2,ny=2,tick.ratio=0.2)
?axis
?axis
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="number of ...",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.8,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.8) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
par(mar=c(2,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.8,font=1)
axis_s<-seq(xlims[1],xlims[2],by=500)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.8) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
pdf(file="ssr_bar.pdf",width=12,height=8)
par(mar=c(2,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.8,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.8) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
pdf(file="ssr_bar.pdf",width=12,height=8)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.9,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.9) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=0.8,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=0.9,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=0.9) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
pdf(file="ssr_bar.pdf",width=12,height=8)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col="darkred",horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,558,43,652,31,610,477,588,99,81,503,562,76,96,495)
length(mycol)
SSR_result
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,558,43,652,31,610,477,588,99,81,503,562,76,96,495)
mycol <-colors()[rep(mycol,20)]
mycol
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,558,43,652,31,610,477,588,99,81,503,562,76,96,495)
mycol <-colors()[mycol]
mycol
pdf(file="ssr_bar.pdf",width=12,height=8)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/2),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
pdf(file="ssr_bar.pdf",width=12,height=8)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
rainbow()
?rainbow
plot(1:20,col=rainbow(20))
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=rainbow(nrow(SSR_result_new)),horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
pdf(file="ssr_bar_rainbow.pdf",width=12,height=8)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=rainbow(nrow(SSR_result_new)),horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
pdf(file="ssr_bar_random.pdf",width=12,height=8)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
q()
history
history()
ls()
SSR_result_new
dir()
SSR_result<-read.delim("SSR_Barplot.txt",header=F,sep="\t")
SSR_result
SSR_result_new<-SSR_result[nrow(SSR_result):1,]
SSR_result_new
pdf(file="ssr_bar_random_0416.pdf",width=12,height=8)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
q()
?tiff
tiff(file="ssr_bar_random_0416.tiff",width=12,height=8)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
q()
?tiff
tiff(filename="ssr_bar_random_0416.tiff",width=12,height=10)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
?tiff
tiff(filename="ssr_bar_random_0416.tiff")
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
tiff(filename="ssr_bar_random_0416.tiff",width=600,height=480)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
?barplot
tiff(filename="ssr_bar_random_0416.tiff",width=600,height=480)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,borders=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
pdf(filename="ssr_bar_random_0416.pdf",width=12,height=10)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,borders=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
pdf(file="ssr_bar_random_0416.pdf",width=12,height=10)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,borders=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
?barplot
pdf(file="ssr_bar_random_0416.pdf",width=12,height=10)
par(mar=c(4,10,4,4))
numbers<-as.numeric(SSR_result_new[[2]])
jump<-max(numbers)/10
xlims=c(0,(max(numbers)+jump))
labels<-as.character(SSR_result_new[[1]])
bar<-barplot(numbers,col=mycol,border=mycol,horiz=TRUE,xlim=ylims,axisnames=FALSE,ylab="",xlab="numbers.",cex.lab=1,font.lab=1,axes=FALSE)
text(rep(-(jump/3),length(bar)),bar,labels=labels,adj=1,xpd=T,cex=1,font=1)
text((numbers+jump/3),bar,labels=as.character(numbers),cex=1,font=1)
axis_s<-seq(xlims[1],xlims[2],by=1000)
axis(side=3,at=axis_s,labels=as.character(axis_s),cex.axis=1) 
#minor.tick(nx=2,ny=2,tick.ratio=0.2)
box()
dev.off()
q()
#----
# ----
getwd()
setwd("/mnt/lustre/users/wangyan/project/sunshengming/RiBenZhaoXia_MJ20120912594/all_RSEM/20130809_buchong")
getwd()
dir()
setwd("/mnt/lustre/users/wangyan/project/sunshengming/RiBenZhaoXia_MJ20120912594/all_RSEM/20130809_buchong/Venn/")
dir()
dir()
GO<-read.delim("GO.list",sep="\t",header=F)
GO[1:4,]
naems(GO)
names(GO)
names(GO)<-c("transcript_id","gene_id")
transcript2gene<-read.delim("contig_singlet_gene2transcript",sep="\t",header=F)
transcript2gene[1:4,]
names(transcript2gene)<-c("gene_id","transcript_id")
names(GO)<-c("transcript_id","GO")
merge(GO,transcript2gene,by="transcript_id")->ll
ll[1:4,]
gene_GO<-ll[,3:2]
gene_GO[1:4,]
write.table(gene_GO,"genes.GO",sep="\t",col.names=F,row.names=F,quote=F)
q()
level2_GO<-read.delim("3h_24h_DEGs_GOs_new_level2.txt",sep="\t",header=T,skip=2,check.names=F)
level2_GO[1:4,]
level2_GO<-read.delim("3h_24h_DEGs_GOs_new_level2.txt",sep="\t",header=T,skip=1,check.names=F)
level2_GO[1:4,]
class_term_num<-unique(level2_GO[,1:3])
class_term_num_sort<-class_term_num[order(as.character(class_term_num[[1]]),as.numeric(class_term_num[[3]]),decreasing=T),]
plot_data<-as.matrix(class_term_num_sort[,3])
rownames(plot_data)<-as.character(class_term_num_sort[[2]])
plot_data
class_term_num
level2_GO<-read.delim("3h_24h_DEGs_GOs_new_level2.txt",sep="\t",header=T,skip=1,check.names=F)
class_term_num<-unique(level2_GO[,1:3])
class_term_num_sort<-class_term_num[order(as.character(class_term_num[[1]]),as.numeric(class_term_num[[3]]),decreasing=T),]
plot_data<-as.matrix(class_term_num_sort)
plot_data
class_term_num_sort
plot_data<-as.matrix(class_term_num_sort[,3])
rownames(plot_data)<-as.character(class_term_num_sort[[2]])
Classes<-as.character(class_term_num_sort[[1]])
plot_data
Classes
N_MF<-length(which(Classes=="molecular_function"))
N_CC<-length(which(Classes=="cellular_component"))
N_BP<-length(which(Classes=="biological_process"))
N_MF
N_CC
N_BP
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
N_MF<-which(Classes=="molecular_function")
N_MF
N_CC<-which(Classes=="cellular_component")
N_BP<-which(Classes=="biological_process")
N_CC
N_BP
N_MF<-length(which(Classes=="molecular_function"))
N_CC<-length(which(Classes=="cellular_component"))
N_BP<-length(which(Classes=="biological_process"))
abline(h=round(bars[N_MF])+1.2)
abline(h=round(bars[N_MF+N_CC])+1.2)
plot_data
MAX(plot_data)
max(plot_data)
class_term_num_sort
level2_GO[1:4,]
class_term_num<-unique(level2_GO[,1:4])
class_term_num_sort<-class_term_num[order(as.character(class_term_num[[1]]),as.numeric(class_term_num[[3]]),decreasing=T),]
#plot_data<-as.matrix(class_term_num_sort[,3])
plot_data<-as.matrix(class_term_num_sort[,4])
plot_data
plot_data<-as.matrix(class_term_num_sort[,4])*100
plot_data
class_term_num_sort[,4]*100
round(class_term_num_sort[,4]*100,1)
round(class_term_num_sort[,4]*100)
as.matrix(round(class_term_num_sort[,4]*100))
round(class_term_num_sort[,4])
plot_data<-as.matrix(round(class_term_num_sort[,4]*100))
plot_data
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,100),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,100,by=20),as.character(seq(0,100,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[N_MF])+1.2)
abline(h=round(bars[N_MF+N_CC])+1.2)
paste(as.character(seq(0,100,by=20)),"%",sep="")
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,100),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,100,by=20),paste(as.character(seq(0,100,by=20)),"%",sep=""),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[N_MF])+1.2)
abline(h=round(bars[N_MF+N_CC])+1.2)
round(class_term_num_sort[,4]
)
round(class_term_num_sort[,4],1)
round(class_term_num_sort[,4],2)
round(class_term_num_sort[,4],3)
seq(0,100,by=20)
seq(0,1,by=0.2)
seq(0,1,by=0.2)*123
round(seq(0,1,by=0.2)*123)
?scan
AllNum<-scan("3h_24h_DEGs_GOs_new_level2.txt",what = character(),nlines=1)
AllNum
AllNum
FirstLine<-scan("3h_24h_DEGs_GOs_new_level2.txt",what = character(),nlines=1)
AllNum<-as.numeric(FirstLine[length(FirstLine)])
AllNum
seq(0,1,by=.2)
round(seq(0,1,by=0.2)*AllNum)
axis(1,seq(0,1,by=.2),as.character(round(seq(0,1,by=0.2)*AllNum)),tcl=1,tck=0.02,font=2)
axis(1,seq(0,100,by=20),as.character(round(seq(0,1,by=0.2)*AllNum)),tcl=1,tck=0.02,font=2)
text(x=98,y=round(bars[N_MF-1]),"Molecular Function",adj=c(1,0),font=2,cex=1.6)
text(x=98,y=round(bars[N_MF+N_CC-1]),"Cellular Component",adj=c(1,0),font=2,cex=1.6)
text(x=98,y=round(bars[N_MF+N_CC+N_BP-1]),"Biological Process",adj=c(1,0),font=2,cex=1.6)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)
bars
rownames(plot_data)
rownames(plot_data)<-as.character(class_term_num_sort[[2]])
rownames(plot_data)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)
###plot
pdf("GO_DiffExpGene.pdf",w=9,h=10)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,100),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,100,by=20),paste(as.character(seq(0,100,by=20)),"%",sep=""),tcl=1,tck=0.02,font=2)
#axis(1,seq(0,100,by=20),as.character(round(seq(0,1,by=0.2)*AllNum)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[N_MF])+1.2)
abline(h=round(bars[N_MF+N_CC])+1.2)
text(x=98,y=round(bars[N_MF-1]),"Molecular Function",adj=c(1,0),font=2,cex=1.6)
text(x=98,y=round(bars[N_MF+N_CC-1]),"Cellular Component",adj=c(1,0),font=2,cex=1.6)
text(x=98,y=round(bars[N_MF+N_CC+N_BP-1]),"Biological Process",adj=c(1,0),font=2,cex=1.6)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)
dev.off()
q()
## BP
par(mar=c(0.5,0.5,4,4))
bars<-barplot(BPs_new_new,horiz=T,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(0.5,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F)
text(rep(max(bars),length(bars)),bars,rownames(BPs_new_new),col="black",adj=1)
## MF
par(mar=c(0.5,0.5,0.5,4))
bars<-barplot(MFs_new_new,horiz=T,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
## 差异基因 GO-annotation 圆饼图.jpg
layout(matrix(c(rep(2,2),rep(1,2),rep(4,2),rep(3,2),rep(6,2),rep(5,2)),3,4,byrow=T))
## BP
par(mar=c(0.5,0.5,4,4))
bars<-barplot(BPs_new_new,horiz=T,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(0.5,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F)
text(rep(max(bars),length(bars)),bars,rownames(BPs_new_new),col="black",adj=1)
## MF
par(mar=c(0.5,0.5,0.5,4))
bars1<-barplot(MFs_new_new,horiz=T,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
##axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(0.5,4,0.5,0.5))
barplot(bars1,beside=T,horiz=T,col="white",border=F,axes=F)
text(rep(max(bars1),length(bars1)),bars1,rownames(MFs_new_new),col="black",adj=1)
## CC
par(mar=c(0.5,0.5,0.5,4))
bars2<-barplot(CCs_new_new,horiz=T,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
##axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(0.5,4,0.5,0.5))
barplot(bars2,beside=T,horiz=T,col="white",border=F,axes=F)
text(rep(max(bars2),length(bars2)),bars2,rownames(CCs_new_new),col="black",adj=1)
?BARPLOT
?barplot
?barplot
 
## 差异基因 GO-annotation 圆饼图.jpg
layout(matrix(c(rep(2,2),rep(1,2),rep(4,2),rep(3,2),rep(6,2),rep(5,2)),3,4,byrow=T))
## BP
par(mar=c(0.5,0.5,4,4))
bars<-barplot(BPs_new_new,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(0.5,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F)
text(rep(max(bars),length(bars)),bars,rownames(BPs_new_new),col="black",adj=1)
## MF
par(mar=c(0.5,0.5,0.5,4))
bars1<-barplot(MFs_new_new,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
##axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(0.5,4,0.5,0.5))
barplot(bars1,beside=T,horiz=T,col="white",border=F,axes=F)
text(rep(max(bars1),length(bars1)),bars1,rownames(MFs_new_new),col="black",adj=1)
## CC
par(mar=c(0.5,0.5,0.5,4))
bars2<-barplot(CCs_new_new,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
##axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(0.5,4,0.5,0.5))
barplot(bars2,beside=T,horiz=T,col="white",border=F,axes=F)
text(rep(max(bars2),length(bars2)),bars2,rownames(CCs_new_new),col="black",adj=1)
?barplot
dim(BPs)
dim(MFs)
dim(CCs)
ls()
SDEG_GO_Class_Des
ls()
SDEG_GOs
ls()
MFs_new_new
MFs_new_new_sort<-MFs_new_new[order(MFs_new_new[,1],decreasing=F),]
MFs_new_new_sort
class(MFs_new_new_sort)
ls()
SDEG_GO_Class_Des[1:4,]
SDEG_GOs_new[1:4,]
ls()
level2_GO[1:4,]
class_term_num<-level2_GO[,1:3]
class_term_num[1:4,]
class_term_num_sort<-class_term_num[order(as.character(class_term_num[[1]]),as.numeric(class_term_num[[3]]),decreasing=F),]
class_term_num_sort
class_term_num<-unique(level2_GO[,1:3])
class_term_num_sort<-class_term_num[order(as.character(class_term_num[[1]]),as.numeric(class_term_num[[3]]),decreasing=F),])
class_term_num_sort<-class_term_num[order(as.character(class_term_num[[1]]),as.numeric(class_term_num[[3]]),decreasing=F),]
class_term_num_sort
plot_data<-as.matrix(class_term_num_sort[,3])
rownames(plot_data)
rownames(plot_data)<-as.character(class_term_num_sort[[2]])
rownames(plot_data)
plot_data
history()
level2_GO[1:2,]
par(mar=c(0.5,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
class_term_num<-unique(level2_GO[,1:3])
class_term_num_sort<-class_term_num[order(as.character(class_term_num[[1]]),as.numeric(class_term_num[[3]]),decreasing=T),]
plot_data<-as.matrix(class_term_num_sort[,3])
rownames(plot_data)<-as.character(class_term_num_sort[[2]])
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(4,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F)
plot_data
rownames(plot_data)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1)
plot_data
rep(max(bars),length(bars))
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(4,4,4,0.5))
bars
axis(2,bars,bars)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
par(mar=c(4,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1)
round(bars[10])
abline(h=round(bars[9]))
abline(h=round(bars[10]))
abline(h=round(bars[11]))
abline(h=round(bars[10]))
abline(h=round(bars[10])+0.5)
abline(h=round(bars[10])+0.8)
abline(h=round(bars[21])+0.8)
abline(h=round(bars[22])+0.8)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
abline(h=round(bars[10])+0.8)
abline(h=round(bars[22])+0.8)
par(mar=c(4,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
abline(h=round(bars[10])+0.9)
abline(h=round(bars[22])+0.9)
par(mar=c(4,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
abline(h=round(bars[10])+1)
abline(h=round(bars[22])+1)
par(mar=c(4,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
abline(h=round(bars[10])+1.1)
abline(h=round(bars[22])+1.1)
par(mar=c(4,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1)
###plot
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
par(mar=c(4,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1)
###plot
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
par(mar=c(4,4,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
###plot
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,8,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
par(mar=c(4,4,8,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[10]),"Molecular Function")
text(x=110,y=round(bars[10]),"Molecular Function",adj=c(1,0))
text(x=110,y=round(bars[10]),"Molecular Function",adj=c(1,0))
###plot
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0))
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0))
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0))
par(mar=c(4,4,6,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
?pdf
?pdf
pdf("GO_DiffExpGene.pdf",w=12,h=16)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,4,6,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
dev.off()
getwd()
pdf("GO_DiffExpGene.pdf",w=12,h=16)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,0.5,6,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
dev.off()
bars
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
barplot(bars,beside=T,horiz=T,border=F,axes=F,space=0.5)
strsplit(rownames(plot_data))
strwidth(rownames(plot_data))
?strwidth
strwidth(rownames(plot_data),cex=1,font=2)
max(strwidth(rownames(plot_data),cex=1,font=2))
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,0.5,6,0.5))
#barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
barplot(1:max(strwidth(rownames(plot_data),cex=1,font=2)),beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,0.5,6,0.5))
barplot(1:max(strwidth(rownames(plot_data),cex=1,font=2)),beside=T,horiz=T,col="black",border=F,axes=F,space=0.5)
1:max(strwidth(rownames(plot_data)
)
max(strwidth(rownames(plot_data)))
strwidth(rownames(plot_data),cex=1,font=2)
max(strwidth(rownames(plot_data),cex=1,font=2))
rownames(plot_data)
bars
max(strwidth(rownames(plot_data),cex=1))
max(strwidth(rownames(plot_data)))
rownames(plot_data))
max(strwidth(as.character(rownames(plot_data))))
?strwidth
?strwidth
 sumex <- expression(sum(x[i], i=1,n), e^{i * pi} == -1)
sumex
strwidth(sumex)
?strwidth
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,0.5,6,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
strwidth(rownames(plot_data),cex=1,font=2)
max(strwidth(rownames(plot_data),cex=1,font=2))
box()
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
pdf("GO_DiffExpGene.pdf",w=12,h=16)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,0.5,6,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
dev.off()
pdf("GO_DiffExpGene.pdf",w=10,h=12)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,6,4))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,0.5,6,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
dev.off()
###plot
pdf("GO_DiffExpGene.pdf",w=10,h=12)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
dev.off()
###plot
pdf("GO_DiffExpGene.pdf",w=9,h=10)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1,font=2)
dev.off()
pdf("GO_DiffExpGene.pdf",w=9,h=10)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2,cex=2)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2.cex=2)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2,cex=2)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=2)
dev.off()
pdf("GO_DiffExpGene.pdf",w=9,h=10)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2,cex=1.6)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2.cex=1.6)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2,cex=1.6)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)
dev.off()
###plot
pdf("GO_DiffExpGene.pdf",w=9,h=10)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2,cex=1.6)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2.cex=1)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2,cex=1.6)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)
dev.off()
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2,cex=1.6)
text(x=110,y=round(bars[20]),"Cellular Component",adj=c(1,0),font=2.cex=1.6)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2,cex=1.6)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)pdf("GO_DiffExpGene.pdf",w=9,h=10)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2,cex=1.6)
text(x=110,y=round(bars[20]),"Cellular Component",adj=c(1,0),font=2,cex=1.6)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2,cex=1.6)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)
dev.off()
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=110,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2,cex=1.6)
text(x=110,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2,cex=1.6)
text(x=110,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2,cex=1.6)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)
pdf("GO_DiffExpGene.pdf",w=9,h=10)
layout(matrix(c(rep(2,2),rep(1,2)),2,2,byrow=F))
par(mar=c(4,0.5,4,6))
bars<-barplot(plot_data,horiz=T,width=1,xlim=c(0,120),beside=T,col="steelblue",space=0.5,border=T,axes=F,bty="o")
axis(3,seq(0,120,by=20),as.character(seq(0,120,by=20)),tcl=1,tck=0.02,font=2)
box()
abline(h=round(bars[10])+1.2)
abline(h=round(bars[22])+1.2)
text(x=115,y=round(bars[9]),"Molecular Function",adj=c(1,0),font=2,cex=1.6)
text(x=115,y=round(bars[21]),"Cellular Component",adj=c(1,0),font=2,cex=1.6)
text(x=115,y=round(bars[43]),"Biological Process",adj=c(1,0),font=2,cex=1.6)
par(mar=c(4,0,4,0.5))
barplot(bars,beside=T,horiz=T,col="white",border=F,axes=F,space=0.5)
text(rep(max(bars),length(bars)),bars,rownames(plot_data),col="black",adj=1,cex=1.2,font=1)
dev.off()
ls()
level2_GO[1:4,]
q()
######### 差异基因 KEGG PATHWAY 富集分析 饼图
N<-20
path_enrich<-read.delim("pathway_enrich.xls",header=T,sep="\t",check.names=F)
path2num<-path_enrich[,c(1,4)]
path2num_sort<-path2num[order(path2num[[2]],decreasing=T),]
num_sum<-sum(path2num_sort[[2]])
paste(path_enrich[["Genes"]],collapse="|")->enrich_gene
enrich_genes<-unlist(strsplit(enrich_gene,";",fix=T))
enrich_genes<-unique(unlist(strsplit(enrich_gene,"|",fix=T)))
Percents<-sapply(path2num_sort[[2]],function(x) x*100/length(enrich_genes))
cbind(path2num_sort,Percents)->path2num2percent
path_enrich_sort<-path_enrich[order(path_enrich[[4]],decreasing=T),]
path_enrich_sort_top10<-path_enrich_sort[1:N,]
top10path_genes<-unique(unlist(strsplit(paste(path_enrich_sort_top10[["Genes"]],collapse="|"),"|",fix=T)))
top10path_num<-path_enrich_sort_top10[,c(1,4)]
percent<-sapply(top10path_num[[2]],function(x) 100*x/length(top10path_genes))
cbind(top10path_num,percent)->top10path_num_percent
paste(top10path_num_percent[[1]],"(",top10path_num_percent[[2]],",",round(top10path_num_percent[[3]],2),"%)",sep="")->lbls_new
## plot
library(RColorBrewer)
library(plotrix)
lbls<-top10path_num_percent[[1]]  
slices_NR<-top10path_num_percent[[3]]
col=rainbow(length(lbls))
#col=brewer.pal(length(lbls_new),"Set1")
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=10)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.75,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
## plot
library(RColorBrewer)
library(plotrix)
lbls<-top10path_num_percent[[1]]  
slices_NR<-top10path_num_percent[[3]]
col=rainbow(length(lbls))
#col=brewer.pal(length(lbls_new),"Set1")
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=10)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.6,cex=0.75,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
## plot
library(RColorBrewer)
library(plotrix)
lbls<-top10path_num_percent[[1]]  
slices_NR<-top10path_num_percent[[3]]
col=rainbow(length(lbls))
#col=brewer.pal(length(lbls_new),"Set1")
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=10)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.75,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
## plot
library(RColorBrewer)
library(plotrix)
lbls<-top10path_num_percent[[1]]  
slices_NR<-top10path_num_percent[[3]]
col=rainbow(length(lbls))
#col=brewer.pal(length(lbls_new),"Set1")
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=10)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.75,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
Q()
q()
transcript_gene<-read.delim("contig_singlet_gene2transcript",sep="\t",header=F)
transcript_gene[1:4,]
names(transcript_gene)
names(transcript_gene)<-c("gene_id","transcript_id")
q()
ls()
transcript_gene[1:4,]
pathway<-read.delim("pathway.txt.new",header=F,sep="\t")
names(pathway)
names(pathway)<-c("transcript_id","KOP")
names(pathway)<-c("transcript_id","KO")
merge(transcript_gene,pathway,by="transcript_id")->ll
head(ll)
gene_pathway<-ll[,2:3]
gene_pathway[1:10,]
write.table(gene_pathway,"gene_pathway.txt",sep="\t",col.names=F,row.names=F,quote=F)
q()
ls()
dir()
path_enrich<-read.delim("pathway_enrichment.xls",header=T,sep="\t")
path_enrich[1:4,]
path_enrich<-read.delim("pathway_enrichment.xls",header=T,sep="\t",check.names=F)
path_enrich[1:4,]
path2num<-path_enrich[,c(1,5)]
path2num[1:4,]
path2num_sort<-path2num[order(path2num[[2]],decreasing=T),]
path2num_sort[1:4,]
path2num<-path_enrich[,c(1,4)]
path2num[1:4,]
path2num_sort<-path2num[order(path2num[[2]],decreasing=T),]
path2num_sort[1:4,]
path2num_sort[1:40,]
num_sum<-sum(path2num_sort[[2]])
num_sum
path_enrich[1:4,]
paste(path_enrich[["Genes"]],collapse="|")->enrich_gene
enrich_gene
enrich_genes<-unlist(strsplit(enrich_gene,";",fix=T))
enrich_genes<-unique(unlist(strsplit(enrich_gene,";",fix=T)))
enrich_genes
enrich_genes<-unique(unlist(strsplit(enrich_gene,"|",fix=T)))
enrich_genes
ls()
path2num_sort[1:4,]
Percents<-sapply(path2num[[2]],function(x) x/length(enrich_genes))
Percents
Percents<-sapply(path2num[[2]],function(x) x*100/length(enrich_genes))
Percents
path2num_sort
Percents<-sapply(path2num_sort[[2]],function(x) x*100/length(enrich_genes))
Percents
path2num_sort
29/length(enrich_genes)
path2num_sort
cbind(path2num_sort,Percents)
cbind(path2num_sort,Percents)[1:4,]
cbind(path2num_sort,Percents)[1:40,]
cbind(path2num_sort,Percents)->path2num2percent
path2num2percent
path2num2percent[1:4,]
pie(path2num2percent[[3]])
history()
history(100)
ls()
path_enrich[1:4,]
path_enrich[1:10,]
names(path_enrich)
path_enrich_sort<-path_enrich[order(path_enrich[[4]],decreasing=T),]
path_enrich_sort[1:4,]
path_enrich_sort[1:10,]
path_enrich_sort_top10<-path_enrich_sort[1:10,]
path_enrich_sort_top10<-path_enrich_sort[1:10,c(1,4)]
path_enrich_sort_top10<-path_enrich_sort[1:10,]
path_enrich_sort_top10
top10path_genes<-unique(unlist(strsplit(paste(path_enrich_sort_top10[["Genes"]],collapse="|"),"|",fix=T)))
top10path_genes
?read.table
ls()
top10path_num
top10path_genes
top10path
ls()
top10path_num<-path_enrich_sort_top10[,c(1,4)]
top10path_num
percent<-sapply(top10path_num[[2]],function(x) 100*x/length(top10path_genes))
percent
cbind(top10path_num,percent)
cbind(top10path_num,percent)->top10path_num_percent
top10path_num_percent
pie(top10path_num_percent[[3]])
top10path_num_percent
library(plotrix)
lbls<-top10path_num_percent[[1]]  
slices_NR<-top10path_num_percent[[3]]
col=rainbow(length(lbls))
pie(slices_NR,labels=Lab2Num_NR[[3]],col=col)
pie(slices_NR,labels=lbls,col=col)
?pie
require(grDevices)
     pie(rep(1, 24), col = rainbow(24), radius = 0.9)
pie.sales <- c(0.12, 0.3, 0.26, 0.16, 0.04, 0.12)
names(pie.sales) <- c("Blueberry", "Cherry",
         "Apple", "Boston Cream", "Other", "Vanilla Cream")
     pie(pie.sales) # default colours
pie(pie.sales, col = c("purple", "violetred1", "green3",
                            "cornsilk", "cyan", "white"))
 pie(pie.sales, col = gray(seq(0.4, 1.0, length = 6)))
pie(pie.sales, density = 10, angle = 15 + 10 * 1:6)
pie(pie.sales, clockwise = TRUE, main = "pie(*, clockwise = TRUE)")
segments(0, 0, 0, 1, col = "red", lwd = 2)
text(0, 1, "init.angle = 90", col = "red")
     
     n <- 200
pie(rep(1, n), labels = "", col = rainbow(n), border = NA,
         main = "pie(*, labels=\"\", col=rainbow(n), border=NA,..")
?pie
pie(slices_NR,labels=lbls,col=col,radius=0.9)
pie(slices_NR,labels=lbls,col=col,radius=0.5)
pie(slices_NR,labels=lbls,col=col,radius=0.1)
pie(slices_NR,labels=lbls,col=col,radius=0.2)
pie(slices_NR,labels=lbls,col=col,radius=0.3)
pie(slices_NR,labels=lbls,col=col,radius=0.5)
pie(slices_NR,labels=lbls,col=col,radius=0.5,cex=1)
pie(slices_NR,labels=lbls,col=col,radius=0.5,cex=0.8)
par(mar=c(4,4,4,4))
pie(slices_NR,labels=lbls,col=col,radius=0.9)
par(mar=c(8,8,8,8))
pie(slices_NR,labels=lbls,col=col,radius=0.9)
pie(slices_NR,labels=lbls,col=col,radius=0.8)
pie(slices_NR,labels=lbls,col=col,radius=0.8)
par(mar=c(8,8,8,8))
pie(slices_NR,labels=lbls,col=col,radius=0.7,cex=0.8)
pie(slices_NR,labels=lbls,col=col,radius=0.7,cex=0.8,font=2)
pie(slices_NR,labels=lbls,col=col,radius=0.7,cex=0.8,font=2)->l
l
summary(l)
??brewer.pal
top10path_num_percent
lbls
paste(top10path_num_percent[[1]],top10path_num_percent[[2]],sep="=")
paste(top10path_num_percent[[1]],"(",top10path_num_percent[[2]],")",sep="")
paste(top10path_num_percent[[1]],"(",top10path_num_percent[[2]],")",sep="")->lbls_new
library(RColorBrewer)
library(plotrix)
lbls_new<-top10path_num_percent[[1]]  
slices_NR<-top10path_num_percent[[3]]
# col=rainbow(length(lbls))
col=brewer.pal(length(lbls_new),"Set1")
#pdf(file="Diff Exp Gene Num in Pathways.pdf",height=5,width=10)
par(mar=c(8,8,8,8))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.8,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways")
lbls_new
paste(top10path_num_percent[[1]],"(",top10path_num_percent[[2]],")",sep="")->lbls_new
library(RColorBrewer)
library(plotrix)
lbls<-top10path_num_percent[[1]]  
slices_NR<-top10path_num_percent[[3]]
# col=rainbow(length(lbls))
col=brewer.pal(length(lbls_new),"Set1")
#pdf(file="Diff Exp Gene Num in Pathways.pdf",height=5,width=10)
par(mar=c(8,8,8,8))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.8,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways")
?pie
par(mar=c(8,8,8,8))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.8,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise)
#pdf(file="Diff Exp Gene Num in Pathways.pdf",height=5,width=10)
par(mar=c(8,8,8,8))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.8,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
par(mar=c(8,8,8,8))
pie(slices_NR,labels=lbls_new,col=col,radius=0.6,cex=0.8,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
par(mar=c(8,8,8,6))
pie(slices_NR,labels=lbls_new,col=col,radius=0.6,cex=0.7,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
par(mar=c(8,8,8,6))
pie(slices_NR,labels=lbls_new,col=col,radius=0.6,cex=0.6,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=10,width=16)
par(mar=c(8,8,8,6))
pie(slices_NR,labels=lbls_new,col=col,radius=0.6,cex=0.8,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
getwd()
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=8,width=16)
par(mar=c(4,4,4,4))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=1,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
clockwise
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=8,width=12)
par(mar=c(4,4,4,2))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=1,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=8,width=14)
par(mar=c(4,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=1,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=8,width=14)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=1,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=12)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=1,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
col=brewer.pal(length(lbls_new),"Set1")
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=10)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=1,font=1,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=10)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.8,font=1,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
percent
paste(top10path_num_percent[[1]],"(",top10path_num_percent[[2]],",",top10path_num_percent[[3]],")",sep="")->lbls_new
lbls_new
paste(top10path_num_percent[[1]],"(",top10path_num_percent[[2]],",",round(top10path_num_percent[[3]],2),"%)",sep="")->lbls_new
lbls_new
library(RColorBrewer)
library(plotrix)
lbls<-top10path_num_percent[[1]]  
slices_NR<-top10path_num_percent[[3]]
# col=rainbow(length(lbls))
col=brewer.pal(length(lbls_new),"Set1")
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=10)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.8,font=1,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
pdf(file="Diff Exp Gene Num in Pathways.pdf",height=6,width=10)
par(mar=c(1,4,4,0))
pie(slices_NR,labels=lbls_new,col=col,radius=0.7,cex=0.75,font=2,border="white",main="Percentage of DEGs in TOP10 Enriched Pathways",clockwise=T)
dev.off()
top10path_num_percent
q()
dir()
geneList<-read.delim("AllSDEGene",header=F)
geneList
names(geneList)
names(geneList)<-"gene_id"
ls()
getwd()
Matrix_3h<-read.delim("/mnt/lustre/users/wangyan/project/sunshengming/RiBenZhaoXia_MJ20120912594/all_RSEM/20130809_buchong/Expression_Analysis/Normoxia_vs_Hypoxia_3h/genes.diff_exp_result",header=T,sep="\t",check.names=F)
Matrix_3h[1:4,]
Matrix_24h<-read.delim("/mnt/lustre/users/wangyan/project/sunshengming/RiBenZhaoXia_MJ20120912594/all_RSEM/20130809_buchong/Expression_Analysis/Normoxia_vs_Hypoxia_24h/genes.diff_exp_result",header=T,sep="\t",check.names=F)
Matrix_O3h<-read.delim("/mnt/lustre/users/wangyan/project/sunshengming/RiBenZhaoXia_MJ20120912594/all_RSEM/20130809_buchong/Expression_Analysis/Normoxia_vs_Reoxygenation_3h/genes.diff_exp_result",header=T,sep="\t",check.names=F)
ls()
merge(geneList,Matrix_3h,by="gene_id")[,c("gene_id","log2FC_2(Hypoxia_3h/Normoxia)")]->logFC1
logFC1[1:4,]
merge(geneList,Matrix_3h,by="gene_id")[,c("gene_id","Normoxia_fpkm","Hypoxia_3h_fpkm","log2FC_2(Hypoxia_3h/Normoxia)")]->logFC1
merge(geneList,Matrix_3h,by="gene_id",all.x=T)[,c("gene_id","Normoxia_fpkm","Hypoxia_3h_fpkm","log2FC_2(Hypoxia_3h/Normoxia)")]->logFC1
merge(logFC1,Matrix_24h,by="gene_id",all.x=T)[,c("gene_id","Normoxia_fpkm","Hypoxia_3h_fpkm","Hypoxia_24h_fpkm","log2FC_2(Hypoxia_3h/Normoxia)","log2FC_2(Hypoxia_24h/Normoxia)")]->logFC2
names(Matrix_24h)
merge(logFC1,Matrix_24h,by="gene_id",all.x=T)[,c("gene_id","Normoxia_fpkm","Hypoxia_3h_fpkm","Hypoxia_24h_fpkm","log2FC_2(Hypoxia_3h/Normoxia)","log2FC_2(Hypoxia_24h/Normoxia)")]->logFC2
merge(logFC1,Matrix_24h,by="gene_id",all.x=T)->logFC
names(logFC)
merge(logFC1,Matrix_24h,by="gene_id",all.x=T)[,c("gene_id","Normoxia_fpkm.x","Hypoxia_3h_fpkm","Hypoxia_24h_fpkm","log2FC_2(Hypoxia_3h/Normoxia)","log2FC_2(Hypoxia_24h/Normoxia)")]->logFC2
names(logFC2)[2]
names(Matrix_O3h)
merge(logFC2,Matrix_O3h,by="gene_id",all.x=T)[,c("gene_id","Normoxia_fpkm","Hypoxia_3h_fpkm","Hypoxia_24h_fpkm","Reoxygenation_3h_fpkm","log2FC_2(Hypoxia_3h/Normoxia)","log2FC_2(Hypoxia_24h/Normoxia)","log2FC_2(Reoxygenation_3h/Normoxia)")]->logFC3
logFC3[1:4,]
dim(logFC3)
dim(geneList)
write.table(logFC3,"SDEG_fpkm_logFC.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
LiangZu<-read.delim("logFC_Matrix.xls",header=T,sep="\t",check.names=F)
logFCMatrix<-read.delim("logFC_Matrix.xls",header=T,sep="\t",check.names=F)
LiangZu<-read.delim("3h_24h_DEGs",header=F,check.names=F)
LiangZu[1:10]
LiangZu[1:10,]
names(LiangZu)
names(LiangZu)<-names(logFCMatrix)[1]
names(LiangZu)
names(logFCMatrix)
merge(LiangZu,logFCMatrix,by="gene_id",all.x=T)->LiangZu_logFCMatrix
dim(LiangZu_logFCMatrix)
dim(LiangZu)
LiangZu_logFCMatrix[1:4,]
write.table(LiangZu_logFCMatrix,"LiangZu_logFCMatrix.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
ls()
names(LiangZu_logFCMatrix)
LiangZu_logFCMatrix_new<-LiangZu_logFCMatrix[,c(1,2,3)]
LiangZu_logFCMatrix_new[1:4,]
names(LiangZu_logFCMatrix_new)
names(LiangZu_logFCMatrix_new)
names(LiangZu_logFCMatrix_new)[2:3]
names(LiangZu_logFCMatrix_new)[2:3]<-c("Hypoxia 3h_vs_Normoxia 3h","Hypoxia 24h_vs_Normoxia 24h")
write.table(LiangZu_logFCMatrix_new,"LiangZu_logFCMatrix.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
dir()
transcript2gene<-read.delim("contig_singlet_gene2transcript",sep="\t",header=F)
transcript2gene[1:4,]
names(transcript2gene)
names(transcript2gene)<-c("gene_id","transcript_id")
dir()
annot<-read.delim("annotation.table.xls",header=T,sep="\t",check.names=F)
annot[1:4,]
names(annot)
names(annot)[1]
names(annot)[1]<-"transcript_id"
ll
merge(annot,transcript2gene,by="transcript_id")->ll
names(ll)
gene_annot<-ll[,c(15,3,4,5,6,7,8,9,10,11,13,14)]
gene_annot<-unique(ll[,c(15,3,4,5,6,7,8,9,10,11,13,14)])
write.table(gene_annot,"gene_annot.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("cluster_diff_knownMir.matrix", header=T, com='', sep="\t")
rownames(data) = data[,1] # set rownames to gene identifiers
data = data[,2:length(data[1,])] # remove the gene column since its now the rowname value
data = as.matrix(data) # convert to matrix
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
myheatcol = redgreen(75)[75:1]
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="cluster_diff_knownMir.matrix.heatmap.pdf", width=8,height=43.05, paper="special");
heatmap.2(centered_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE,keysize=1.2,cexCol=2.5, margins=c(8,8), lhei=c(0.3,2), lwid=c(2.5,4))
dev.off()
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("cluster_diff_novoMir.matrix", header=T, com='', sep="\t")
rownames(data) = data[,1] # set rownames to gene identifiers
data = data[,2:length(data[1,])] # remove the gene column since its now the rowname value
data = as.matrix(data) # convert to matrix
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") # cluster genes
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") # cluster conditions
myheatcol = redgreen(75)[75:1]
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="cluster_diff_novoMir.matrix.heatmap.pdf", width=8,height=38.85, paper="special");
heatmap.2(centered_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(8,8), lhei=c(0.3,2), lwid=c(2.5,4))
dev.off()
q
q()
miRanda_rno<-read.delim("rn4_predictions_S_C_aug2010.txt",header=T,sep="\t",checknames=F)
miRanda_rno<-read.delim("rn4_predictions_S_C_aug2010.txt",header=Tsep="\t",check.names=F)
miRanda_rno<-read.delim("rn4_predictions_S_C_aug2010.txt",header=Tsep="\t",check.names=FALSE)
miRanda_rno<-read.delim("rn4_predictions_S_C_aug2010.txt",header=T,sep="\t",check.names=FALSE)
q()
ls()
history()
miRanda_rno<-read.delim("rn4_predictions_S_C_aug2010.txt",header=T,sep="\t",check.names=FALSE)
rno_geneid_ensgid<-read.delim("rno_rno_gene_id_rno_ensg_id.txt",sep="\t",header=T,check.names=T)
miRanda_rno[1:4,]
rno_geneid_ensgid[1:4,]
names(rno_geneid_ensgid)[1]<-"gene_id"
miRanda_rno_new<-merge(miRanda_rno,rno_geneid_ensgid,by="gene_id")
miRanda_rno_new[1:4,]
names(miRanda_rno_new)
miRanda_rno<-miRanda_rno_new[,names(miRanda_rno_new)[c(2:4,1,5:20)]]
miRanda_rno[1,]
miRanda_rno<-miRanda_rno_new[,names(miRanda_rno_new)[c(2:4,1,20,5:19)]]
miRanda_rno[1:4,]
miRanda_rno_detail_info<-miRanda_rno
miRanda_rno_brief_info<-miRanda_rno[,c(2,3,5)]
write.table(miRanda_rno_detail_info,"miRanda_rno_detail_info",sep="\t",col.names=T,row.names=F,quote=F)
write.table(miRanda_rno_brief_info,"miRanda_rno_brief_info",sep="\t",col.names=T,row.names=F,quote=F)
q()
dir()
S1_GO<-read.delim("S1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
W1_GO<-read.delim("W1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
S1_GO[1:4,]
co_GO<-merge(S1_GO,W1_GO,by="id")
co_GO[1:4,]
names(co_GO)
q()
co_GO_sort
## 两组差异表达基因 GO 富集分析结果统计及展示（筛选后的结果：fdr<=0.05）
S1_GO<-read.delim("S1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
W1_GO<-read.delim("W1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
GOterm_class_id<-read.delim("/mnt/lustre/users/wangyan/Data/GO_Related/GOterm_3class_GOID",sep="\t",header=T)
S1_GO_class<-merge(S1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
W1_GO_class<-merge(W1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
names(S1_GO_class)[2]<-"description"
names(W1_GO_class)[2]<-"description"
##交集
co_GO<-merge(S1_GO_class,W1_GO_class,by="id")[,c("id","description.x","class.x","p_uncorrected.x","p_uncorrected.y")]
names(co_GO)<-c("id","description","class","p_uncorrected.S1","p_uncorrected.W1")
log10pvalue_S1<--log10(co_GO[["p_uncorrected.S1"]])
log10pvalue_W1<--log10(co_GO[["p_uncorrected.W1"]])
co_GO_NEW<-cbind(co_GO[,2:3],log10pvalue_S1,log10pvalue_W1)
co_GO_sort<-co_GO_NEW[order(co_GO_NEW[["class"]],as.numeric(co_GO_NEW[["log10pvalue_S1"]]),decreasing=T),]
co_GO_sort
matrixp<-t(as.matrix(co_GO_sort[,3:4]))
matrixp
GOtype<-as.character(co_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
mycol
c(0,max(matrixp)+2)
round(13.28233)
barplot(matrixp,col=mycol,ylim=c(0,round(max(matrixp)+2)),density=rep(c(20,0),dim(matrixp)[1]))
barplot(matrixp,col=mycol,ylim=c(0,round(max(matrixp)+2)),density=rep(c(20,0),dim(matrixp)[1]),beside=T)
barplot(matrixp,col=mycol,ylim=c(0,round(max(matrixp)+2)),density=rep(c(100,20),dim(matrixp)[1]),beside=T)
pdf("test.pdf",w=12,h=10)
barplot(matrixp,col=mycol,ylim=c(0,round(max(matrixp)+2)),density=rep(c(100,20),dim(matrixp)[1]),beside=T)
dev.off()
gwtwd()
getwd()
rep(mycol,2)
mycol
?rep
rep(mycol,each=2)
mycol
pdf("test.pdf",w=12,h=10)
barplot(matrixp,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp)+2)),density=rep(c(100,20),dim(matrixp)[1]),beside=T)
dev.off()
?barplot
barplot(matrixp,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp)+2)),density=rep(c(100,20),dim(matrixp)[1]),beside=T,border=NA)
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp_co)+2)),density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
matrixp_co<-t(as.matrix(co_GO_sort[,3:4]))
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp_co)+2)),density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
dev.off()
grid
?grid
plot(1:3)
grid(NA, 5, lwd = 2) 
pdf("ddd.pdf")
plot(1:3)
grid(NA, 5, lwd = 2) 
dev.off()
?grid
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp_co)+2)),density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
grid()
dev.off()
ylims<-c(0,round(max(matrixp_co)+2))
ylims
seq(0,ylims,by=1)
seq(0,ylims[2],by=1)
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
ylims<-c(0,round(max(matrixp_co)+2))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
abline(h=seq(1,ylims[2],by=1),col="darkgrey",lwd=1,lty=2)
dev.off()
bar
mean(bar)
bar[1,]+0.5
co_GO_sort
ylims<-c(0,round(max(matrixp_co)+2))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
abline(h=seq(1,ylims[2],by=1),col="darkgrey",lwd=1,lty=2)
axis(2,ylims,as.character(ylims),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=90,adj=1,xpd=T,cex=0.5,font=1,col="black")
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
ylims<-c(0,round(max(matrixp_co)+2))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
abline(h=seq(1,ylims[2],by=1),col="darkgrey",lwd=1,lty=2)
axis(2,ylims,as.character(ylims),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=90,adj=1,xpd=T,cex=0.5,font=1,col="black")
dev.off()
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
ylims<-c(0,round(max(matrixp_co)+2))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
abline(h=seq(1,ylims[2],by=1),col="darkgrey",lwd=1,lty=2)
axis(2,seq(1,ylims[2],by=1),as.character(seq(1,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=90,adj=1,xpd=T,cex=0.5,font=1,col="black")
dev.off()
?layout
S1_GOids<-as.data.frame(setdiff(as.character(S1_GO_class[["id"]]),as.character(co_GO[["id"]])))
names(S1_GOids)<-"id"
S1_specific_GO<-merge(S1_GOids,S1_GO_class,by="id")
log10pvalue_S1<--log10(S1_specific_GO[["p_uncorrected"]])
S1_specific_GO_new<-cbind(S1_specific_GO[,2:3],log10pvalue_S1)
S1_specific_GO_sort<-S1_specific_GO_new[order(S1_specific_GO_new[["class"]],S1_specific_GO_new[["log10pvalue_S1"]],decreasing=T),]
S1_specific_GO_sort
q()
dir()
S1_GO<-read.delim("S1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
W1_GO<-read.delim("W1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
S1_GO[1:4,]
co_GO<-merge(S1_GO,W1_GO,by="id")
co_GO[1:4,]
names(co_GO)
q()
co_GO_sort
## 两组差异表达基因 GO 富集分析结果统计及展示（筛选后的结果：fdr<=0.05）
S1_GO<-read.delim("S1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
W1_GO<-read.delim("W1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
GOterm_class_id<-read.delim("/mnt/lustre/users/wangyan/Data/GO_Related/GOterm_3class_GOID",sep="\t",header=T)
S1_GO_class<-merge(S1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
W1_GO_class<-merge(W1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
names(S1_GO_class)[2]<-"description"
names(W1_GO_class)[2]<-"description"
##交集
co_GO<-merge(S1_GO_class,W1_GO_class,by="id")[,c("id","description.x","class.x","p_uncorrected.x","p_uncorrected.y")]
names(co_GO)<-c("id","description","class","p_uncorrected.S1","p_uncorrected.W1")
log10pvalue_S1<--log10(co_GO[["p_uncorrected.S1"]])
log10pvalue_W1<--log10(co_GO[["p_uncorrected.W1"]])
co_GO_NEW<-cbind(co_GO[,2:3],log10pvalue_S1,log10pvalue_W1)
co_GO_sort<-co_GO_NEW[order(co_GO_NEW[["class"]],as.numeric(co_GO_NEW[["log10pvalue_S1"]]),decreasing=T),]
co_GO_sort
matrixp<-t(as.matrix(co_GO_sort[,3:4]))
matrixp
GOtype<-as.character(co_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
mycol
c(0,max(matrixp)+2)
round(13.28233)
barplot(matrixp,col=mycol,ylim=c(0,round(max(matrixp)+2)),density=rep(c(20,0),dim(matrixp)[1]))
barplot(matrixp,col=mycol,ylim=c(0,round(max(matrixp)+2)),density=rep(c(20,0),dim(matrixp)[1]),beside=T)
barplot(matrixp,col=mycol,ylim=c(0,round(max(matrixp)+2)),density=rep(c(100,20),dim(matrixp)[1]),beside=T)
pdf("test.pdf",w=12,h=10)
barplot(matrixp,col=mycol,ylim=c(0,round(max(matrixp)+2)),density=rep(c(100,20),dim(matrixp)[1]),beside=T)
dev.off()
gwtwd()
getwd()
rep(mycol,2)
mycol
?rep
rep(mycol,each=2)
mycol
pdf("test.pdf",w=12,h=10)
barplot(matrixp,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp)+2)),density=rep(c(100,20),dim(matrixp)[1]),beside=T)
dev.off()
?barplot
barplot(matrixp,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp)+2)),density=rep(c(100,20),dim(matrixp)[1]),beside=T,border=NA)
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp_co)+2)),density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
matrixp_co<-t(as.matrix(co_GO_sort[,3:4]))
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp_co)+2)),density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
dev.off()
grid
?grid
plot(1:3)
grid(NA, 5, lwd = 2) 
pdf("ddd.pdf")
plot(1:3)
grid(NA, 5, lwd = 2) 
dev.off()
?grid
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=c(0,round(max(matrixp_co)+2)),density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
grid()
dev.off()
ylims<-c(0,round(max(matrixp_co)+2))
ylims
seq(0,ylims,by=1)
seq(0,ylims[2],by=1)
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
ylims<-c(0,round(max(matrixp_co)+2))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
abline(h=seq(1,ylims[2],by=1),col="darkgrey",lwd=1,lty=2)
dev.off()
bar
mean(bar)
bar[1,]+0.5
co_GO_sort
ylims<-c(0,round(max(matrixp_co)+2))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
abline(h=seq(1,ylims[2],by=1),col="darkgrey",lwd=1,lty=2)
axis(2,ylims,as.character(ylims),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=90,adj=1,xpd=T,cex=0.5,font=1,col="black")
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
ylims<-c(0,round(max(matrixp_co)+2))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
abline(h=seq(1,ylims[2],by=1),col="darkgrey",lwd=1,lty=2)
axis(2,ylims,as.character(ylims),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=90,adj=1,xpd=T,cex=0.5,font=1,col="black")
dev.off()
pdf("GO Enrichment Distribution.pdf",w=15,h=10)
ylims<-c(0,round(max(matrixp_co)+2))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(200,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black")
abline(h=seq(1,ylims[2],by=1),col="darkgrey",lwd=1,lty=2)
axis(2,seq(1,ylims[2],by=1),as.character(seq(1,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=90,adj=1,xpd=T,cex=0.5,font=1,col="black")
dev.off()
?layout
S1_GOids<-as.data.frame(setdiff(as.character(S1_GO_class[["id"]]),as.character(co_GO[["id"]])))
names(S1_GOids)<-"id"
S1_specific_GO<-merge(S1_GOids,S1_GO_class,by="id")
log10pvalue_S1<--log10(S1_specific_GO[["p_uncorrected"]])
S1_specific_GO_new<-cbind(S1_specific_GO[,2:3],log10pvalue_S1)
S1_specific_GO_sort<-S1_specific_GO_new[order(S1_specific_GO_new[["class"]],S1_specific_GO_new[["log10pvalue_S1"]],decreasing=T),]
S1_specific_GO_sort
q()
ls()
co_GO_sort[1:4,]
co_GO[1:4,]
co_GO_S1_W1<-merge(S1_GO_class,W1_GO_class,by="id")
co_GO_S1_W1[1:4,]
S1_GO_class<-merge(S1_GO,GOterm_class_id,by="id")
W1_GO_class<-merge(W1_GO,GOterm_class_id,by="id")
names(S1_GO_class)[2]<-"description"
names(W1_GO_class)[2]<-"description"
co_GO<-merge(S1_GO_class,W1_GO_class,by="id")
co_GO[1:4,]
write.table(co_GO,"co_GOs_S1_W1.xls",sep="\t",col.names=T,row.names=F,quote=F)
gwtwd()
getwd()
S1_specific_GO<-merge(S1_GOids,S1_GO_class,by="id")
S1_specific_GO[1:4,]
write.table(S1_specific_GO,"S1_specific_GO<-merge(S1_GOids,S1_GO_class,by="id")")
write.table(S1_specific_GO,"S1_specific_GO.xls",sep="\t",col.names=T,row.names=F,quote=F)
W1_specific_GO<-merge(W1_GOids,W1_GO_class,by="id")
W1_GOids<-as.data.frame(setdiff(as.character(W1_GO_class[["id"]]),as.character(co_GO[["id"]])))
names(W1_GOids)<-"id"
W1_GOids
W1_specific_GO<-merge(W1_GOids,W1_GO_class,by="id")
W1_specific_GO
write.table(W1_specific_GO,"W1_specific_GO.xls",sep="\t",col.names=T,row.names=F,quote=F)
getwd()
setwd("/mnt/lustre/users/wangyan/project/liuwei/DGE/buchongfenxi/KEGG")
getwd()
dir()
logFC2color<-read.delim("log2FC_颜色对应关系",sep=" ",header=F)
logFC2color)
plot(1:10,rep(2,10),col=as.character(logFC2color[[3]]),pch=15,cex=5)
lables<-paste(logFC2color[[1]],logFC2color[[2]],sep=" ~ ")
lables
legend("top",lables,pch=15)
plot(1,0,xlim=c(0,10),ylim=c(0,10),cex=0)
legend("right",lables,pch=15,col=as.character(logFC2color[[3]]),cex=5)
plot(1,0,xlim=c(0,10),ylim=c(0,10),cex=0)
lables<-paste(logFC2color[[1]],logFC2color[[2]],sep=" ~ ")
legend("right",lables,pch=15,col=as.character(logFC2color[[3]]),cex=2)
plot(1,0,xlim=c(0,10),ylim=c(0,10),cex=0)
lables<-paste(logFC2color[[1]],logFC2color[[2]],sep=" ~ ")
legend("top",lables,pch=15,col=as.character(logFC2color[[3]]),cex=2)
dir()
pdf("logFC_color.pdf")
plot(1,0,xlim=c(0,10),ylim=c(0,10),cex=0)
lables<-paste(logFC2color[[1]],logFC2color[[2]],sep=" ~ ")
legend("top",lables,pch=15,col=as.character(logFC2color[[3]]),cex=2)
dev.off()
ls()
q()
dir()
transcript2gene2name<-read.delim("transcript2gene2name",header=F,sep="\t")
transcript2gene2name[1:4,]
names(transcript2gene2name)<-c("transcript","gene","name")
ncbi1<-read.delim("ncbi_ENSRT_rno_KO_logFCs_20131024.txt",header=T,sep="\t",check.names=F)
ncbi1[1:4,]
ncbi2<-read.delim("ncbi_logFCs_colors_20131024.txt",header=T,sep="\t",check.names=F)
ncbi2[1:4,]
names(transcript2gene2name)
names(transcript2gene2name)[1]<-ENSRT
names(transcript2gene2name)[1]<-"ENSRT"
merge(ncbi1,transcript2gene2name,by="ENSRT")->merge1
dim(merge1)
dim(ncbi1)
merge1[1:4,]
merge(merge1,ncbi2,by="ENSRT")->merge2
merge(merge1,ncbi2,by="ncbi")->merge2
dim(merge2)
merge2[1:4,]
dim(ncbi2)
write.table(merge2,"merge2.xls",sep="\t",col.names=T,row.names=F,quote=F)
names(merge2)
merge2_sort<-merge2[order(as.character(merge2[["name"]]),as.character(merge2[["rno"]]),as.character(merge2[["KO"]],as.character(merge2[["ENSRT"]]))),]
write.table(merge2_sort,"merge2_sort.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
W1_C1<-read.delim("W1_vs_C1_PathwayEnrich",header=T,sep="\t",check.names=F)
S1_C1<-read.delim("S1_vs_C1_PathwayEnrich",header=T,sep="\t",check.names=F)
names(W1_C1)
W1_C1_sig<-W1_C1[which(W1_C1[["P-Value"]]<=0.1),] 
S1_C1_sig<-S1_C1[which(S1_C1[["P-Value"]]<=0.1),] 
W1_C1_sig[[1]]
S1_C1_sig[[1]]
as.charactrer(W1_C1_sig[[1]])
as.character(W1_C1_sig[[1]])
as.character(S1_C1_sig[[1]])
as.character(S1_C1_sig[[1]])unique(c(as.character(S1_C1_sig[[1]]),as.character(W1_C1_sig[[1]])))
unique(c(as.character(S1_C1_sig[[1]]),as.character(W1_C1_sig[[1]])))
W1_C1_sig<-W1_C1[which(W1_C1[["P-Value"]]<=0.05),] 
W1_C1_sig<-W1_C1[which(W1_C1[["P-Value"]]<=0.05),] 
S1_C1_sig<-S1_C1[which(S1_C1[["P-Value"]]<=0.05),] 
as.character(S1_C1_sig[[1]])
intersect(as.character(W1_C1[[1]]),as.character(S1_C1[[1]]))
?merge
coPath<-merge(S1_C1,W1_C1,by="Id")
dim(coPath)
write.table(coPath,"coPathway.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
method<-"cuffdiff"
cuffdiffdir<-"/mnt/lustre/users/wangyan/project/liuwei/DGE/ExpressionAnalysis"
DEGseqdirs<-"no"
rsemdir<-"no"
ergdir<-"no"
ertdir<-"no"
FDR<-0.05
pvalue<-0
logFC<-1
width<-12
height<-10
geneAnnPath<-"no"
isoformAnnPath<-"no"
hmp_groups<-"all"
# output dir name prepare
if(logFC!=0 && pvalue==0 && FDR!=0){SDERole<-paste("FDR:",FDR,"-logFC:",logFC,sep="")}
if(logFC==0 && pvalue==0 && FDR!=0){SDERole<-paste("FDR:",FDR,sep="")}
if(logFC==0 && pvalue!=0 && FDR==0){SDERole<-paste("pvalue:",pvalue,sep="")}
if(logFC!=0 && pvalue==0 && FDR==0){SDERole<-paste("logFC:",logFC,sep="")}
if(logFC!=0 && pvalue!=0 && FDR==0){SDERole<-paste("pvalue:",pvalue,"-logFC:",logFC,sep="")}
if(logFC==0 && pvalue!=0 && FDR!=0){SDERole<-paste("FDR:",FDR,"-pvalue:",pvalue,sep="")}
if(logFC!=0 && pvalue!=0 && FDR!=0){SDERole<-paste("FDR:",FDR,"logFC:",logFC,"pvalue:",pvalue,sep="")}
outfilename<-paste(getwd(),"/","diffexpstat","-",method,"-",SDERole,sep="")
q()
## 两组差异表达基因 GO 富集分析结果统计及展示（筛选后的结果：fdr<=0.05）
S1_GO<-read.delim("S1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
W1_GO<-read.delim("W1_vs_C1_go_enrichment.xls",sep="\t",header=T,check.names=F)
GOterm_class_id<-read.delim("/mnt/lustre/users/wangyan/Data/GO_Related/GOterm_3class_GOID",sep="\t",header=T)
S1_GO_class<-merge(S1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
W1_GO_class<-merge(W1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
names(S1_GO_class)[2]<-"description"
names(W1_GO_class)[2]<-"description"
## 两组差异表达基因 GO 富集分析结果统计及展示（筛选后的结果：fdr<=0.05）
S1_GO<-read.delim("W1_vs_C1_DEgene.list.enrichment.detail.xls",sep="\t",header=T,check.names=F)
W1_GO<-read.delim("S1_vs_C1_DEgene.list.enrichment.detail.xls",sep="\t",header=T,check.names=F)
GOterm_class_id<-read.delim("/mnt/lustre/users/wangyan/Data/GO_Related/GOterm_3class_GOID",sep="\t",header=T)
S1_GO_class<-merge(S1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
W1_GO_class<-merge(W1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
names(S1_GO_class)[2]<-"description"
names(W1_GO_class)[2]<-"description"
co_GO<-merge(S1_GO_class,W1_GO_class,by="id")[,c("id","description.x","class.x","p_uncorrected.x","p_uncorrected.y")]
names(co_GO)<-c("id","description","class","p_uncorrected.S1","p_uncorrected.W1")
log10pvalue_S1<--log10(co_GO[["p_uncorrected.S1"]])
log10pvalue_W1<--log10(co_GO[["p_uncorrected.W1"]])
co_GO_NEW<-cbind(co_GO[,2:3],log10pvalue_S1,log10pvalue_W1)
co_GO_sort<-co_GO_NEW[order(co_GO_NEW[["class"]],as.numeric(co_GO_NEW[["log10pvalue_S1"]]),decreasing=T),]
GOtype<-as.character(co_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_co<-t(as.matrix(co_GO_sort[,3:4]))
ylims<-c(0,round(max(matrixp_co)+2))
pdf("Co Enriched GO Terms for S1 and W1 groups.pdf",w=16,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(1000,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Co Enrichment GO Terms for S1 Group and W1 Group")
matrixp_co
dim(matrixp_co)
S1_GO[1:4,]
dim(S1_GO)
S1_GO[1:10,1:5]
S1_GO[1:10,1:8]
S1_GO[1:10,1:9]
S1_GO[1:10,1:10]
S1_GO[10:20,1:10]
S1_GO[10:60,1:10]
S1_GO[60:100,10]
names(S1_GO)
## 两组差异表达基因 GO 富集分析结果统计及展示（筛选后的结果：fdr<=0.05）
S1_GO_raw<-read.delim("W1_vs_C1_DEgene.list.enrichment.detail.xls",sep="\t",header=T,check.names=F)
W1_GO_raw<-read.delim("S1_vs_C1_DEgene.list.enrichment.detail.xls",sep="\t",header=T,check.names=F)
S1_GO<-S1_GO_raw[which(as.numeric(S1_GO_raw[["p_fdr"]])<=0.05),]
dim(S1_GO)
dim(S1_GO_raw)
## 两组差异表达基因 GO 富集分析结果统计及展示（筛选后的结果：fdr<=0.05）
S1_GO_raw<-read.delim("W1_vs_C1_DEgene.list.enrichment.detail.xls",sep="\t",header=T,check.names=F)
W1_GO_raw<-read.delim("S1_vs_C1_DEgene.list.enrichment.detail.xls",sep="\t",header=T,check.names=F)
S1_GO<-S1_GO_raw[which(as.numeric(S1_GO_raw[["p_fdr"]])<=0.05),]
W1_GO<-W1_GO_raw[which(as.numeric(W1_GO_raw[["p_fdr"]])<=0.05),]
GOterm_class_id<-read.delim("/mnt/lustre/users/wangyan/Data/GO_Related/GOterm_3class_GOID",sep="\t",header=T)
S1_GO_class<-merge(S1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
W1_GO_class<-merge(W1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
names(S1_GO_class)[2]<-"description"
names(W1_GO_class)[2]<-"description"
##layout(matrix(c(rep(2,12),rep(c(1,1,3,3,4,4),2)),4,6,byrow=T))
##交集
co_GO<-merge(S1_GO_class,W1_GO_class,by="id")[,c("id","description.x","class.x","p_uncorrected.x","p_uncorrected.y")]
names(co_GO)<-c("id","description","class","p_uncorrected.S1","p_uncorrected.W1")
log10pvalue_S1<--log10(co_GO[["p_uncorrected.S1"]])
log10pvalue_W1<--log10(co_GO[["p_uncorrected.W1"]])
co_GO_NEW<-cbind(co_GO[,2:3],log10pvalue_S1,log10pvalue_W1)
co_GO_sort<-co_GO_NEW[order(co_GO_NEW[["class"]],as.numeric(co_GO_NEW[["log10pvalue_S1"]]),decreasing=T),]
GOtype<-as.character(co_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_co<-t(as.matrix(co_GO_sort[,3:4]))
ylims<-c(0,round(max(matrixp_co)+2))
pdf("Co Enriched GO Terms for S1 and W1 groups.pdf",w=16,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(1000,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Co Enrichment GO Terms for S1 Group and W1 Group")
matrixp_co
dim(matrixp_co)
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(1000,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Co Enrichment GO Terms for S1 Group and W1 Group")
pdf("Co Enriched GO Terms for S1 and W1 groups.pdf",w=16,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(1000,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Co Enrichment GO Terms for S1 Group and W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
getwd()
pdf("Co Enriched GO Terms for S1 and W1 groups.pdf",w=18,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(1000,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Co Enrichment GO Terms for S1 Group and W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
## S1特异的
S1_GOids<-as.data.frame(setdiff(as.character(S1_GO_class[["id"]]),as.character(co_GO[["id"]])))
names(S1_GOids)<-"id"
S1_specific_GO<-merge(S1_GOids,S1_GO_class,by="id")
log10pvalue_S1<--log10(S1_specific_GO[["p_uncorrected"]])
S1_specific_GO_new<-cbind(S1_specific_GO[,2:3],log10pvalue_S1)
S1_specific_GO_new_psort<-S1_specific_GO_new[order(S1_specific_GO_new[["log10pvalue_S1"]],decreasing=T),][1:60,]
S1_specific_GO_sort<-S1_specific_GO_new_psort[order(S1_specific_GO_new_psort[["class"]],S1_specific_GO_new_psort[["log10pvalue_S1"]],decreasing=T),]
GOtype<-as.character(S1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_S1<-S1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_S1)+2))
pdf("S1 Specific Enriched GO Terms.pdf",w=16,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_S1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only S1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=S1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
## W1特异的
W1_GOids<-as.data.frame(setdiff(as.character(W1_GO_class[["id"]]),as.character(co_GO[["id"]])))
names(W1_GOids)<-"id"
W1_specific_GO<-merge(W1_GOids,W1_GO_class,by="id")
log10pvalue_W1<--log10(W1_specific_GO[["p_uncorrected"]])
W1_specific_GO_new<-cbind(W1_specific_GO[,2:3],log10pvalue_W1)
W1_specific_GO_sort<-W1_specific_GO_new[order(W1_specific_GO_new[["class"]],W1_specific_GO_new[["log10pvalue_W1"]],decreasing=T),]
GOtype<-as.character(W1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_W1<-W1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_W1)+2))
pdf("W1 Specific Enriched GO Terms.pdf",w=10,h=6)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_W1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=W1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
matrixp_S1
names(S1_GO_class)
S1_specific_GO_sort
S1_GOids
S1_specific_GO
S1_specific_GO_sort
S1_specific_GO_new_psort
S1_specific_GO_new
S1_specific_GO_new_psort<-S1_specific_GO_new[order(S1_specific_GO_new[["log10pvalue_S1"]],decreasing=T),]
S1_specific_GO_sort<-S1_specific_GO_new_psort[order(S1_specific_GO_new_psort[["class"]],S1_specific_GO_new_psort[["log10pvalue_S1"]],decreasing=T),]
S1_specific_GO_sort
GOtype<-as.character(S1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_S1<-S1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_S1)+2))
pdf("S1 Specific Enriched GO Terms.pdf",w=16,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_S1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only S1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=S1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
W1_GOids
length(W1_GOids)
dim(W1_GOids)
dim(W1_GO_class)
W1_GO_class[1:4,1:10]
W1_GO_class[1:4,]
W1_GO[1:4,1:10]
W1_GO[,10]
sort(W1_GO[,10])
## 两组差异表达基因 GO 富集分析结果统计及展示（筛选后的结果：fdr<=0.05）
W1_GO_raw<-read.delim("W1_vs_C1_DEgene.list.enrichment.detail.xls",sep="\t",header=T,check.names=F)
S1_GO_raw<-read.delim("S1_vs_C1_DEgene.list.enrichment.detail.xls",sep="\t",header=T,check.names=F)
S1_GO<-S1_GO_raw[which(as.numeric(S1_GO_raw[["p_fdr"]])<=0.05),]
W1_GO<-W1_GO_raw[which(as.numeric(W1_GO_raw[["p_fdr"]])<=0.05),]
GOterm_class_id<-read.delim("/mnt/lustre/users/wangyan/Data/GO_Related/GOterm_3class_GOID",sep="\t",header=T)
S1_GO_class<-merge(S1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
W1_GO_class<-merge(W1_GO,GOterm_class_id,by="id")[,c("id","description.x","class","p_uncorrected")]
names(S1_GO_class)[2]<-"description"
names(W1_GO_class)[2]<-"description"
##layout(matrix(c(rep(2,12),rep(c(1,1,3,3,4,4),2)),4,6,byrow=T))
##交集
co_GO<-merge(S1_GO_class,W1_GO_class,by="id")[,c("id","description.x","class.x","p_uncorrected.x","p_uncorrected.y")]
names(co_GO)<-c("id","description","class","p_uncorrected.S1","p_uncorrected.W1")
log10pvalue_S1<--log10(co_GO[["p_uncorrected.S1"]])
log10pvalue_W1<--log10(co_GO[["p_uncorrected.W1"]])
co_GO_NEW<-cbind(co_GO[,2:3],log10pvalue_S1,log10pvalue_W1)
co_GO_sort<-co_GO_NEW[order(co_GO_NEW[["class"]],as.numeric(co_GO_NEW[["log10pvalue_S1"]]),decreasing=T),]
GOtype<-as.character(co_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_co<-t(as.matrix(co_GO_sort[,3:4]))
ylims<-c(0,round(max(matrixp_co)+2))
pdf("Co Enriched GO Terms for S1 and W1 groups.pdf",w=18,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_co,beside=T,col=rep(mycol,each=2),ylim=ylims,density=rep(c(1000,20),dim(matrixp_co)[2]),axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Co Enrichment GO Terms for S1 Group and W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar[1,]+0.5,rep(-0.2,length(bar)),labels=co_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
dim(matrixp_co)
## S1特异的
S1_GOids<-as.data.frame(setdiff(as.character(S1_GO_class[["id"]]),as.character(co_GO[["id"]])))
names(S1_GOids)<-"id"
S1_specific_GO<-merge(S1_GOids,S1_GO_class,by="id")
log10pvalue_S1<--log10(S1_specific_GO[["p_uncorrected"]])
S1_specific_GO_new<-cbind(S1_specific_GO[,2:3],log10pvalue_S1)
S1_specific_GO_new_psort<-S1_specific_GO_new[order(S1_specific_GO_new[["log10pvalue_S1"]],decreasing=T),][1:100,]
S1_specific_GO_sort<-S1_specific_GO_new_psort[order(S1_specific_GO_new_psort[["class"]],S1_specific_GO_new_psort[["log10pvalue_S1"]],decreasing=T),]
GOtype<-as.character(S1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_S1<-S1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_S1)+2))
pdf("S1 Specific Enriched GO Terms.pdf",w=16,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_S1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only S1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=S1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
## W1特异的
W1_GOids<-as.data.frame(setdiff(as.character(W1_GO_class[["id"]]),as.character(co_GO[["id"]])))
names(W1_GOids)<-"id"
W1_specific_GO<-merge(W1_GOids,W1_GO_class,by="id")
log10pvalue_W1<--log10(W1_specific_GO[["p_uncorrected"]])
W1_specific_GO_new<-cbind(W1_specific_GO[,2:3],log10pvalue_W1)
W1_specific_GO_sort<-W1_specific_GO_new[order(W1_specific_GO_new[["class"]],W1_specific_GO_new[["log10pvalue_W1"]],decreasing=T),]
GOtype<-as.character(W1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_W1<-W1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_W1)+2))
pdf("W1 Specific Enriched GO Terms.pdf",w=10,h=6)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_W1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=W1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
## S1特异的
S1_GOids<-as.data.frame(setdiff(as.character(S1_GO_class[["id"]]),as.character(co_GO[["id"]])))
names(S1_GOids)<-"id"
S1_specific_GO<-merge(S1_GOids,S1_GO_class,by="id")
log10pvalue_S1<--log10(S1_specific_GO[["p_uncorrected"]])
S1_specific_GO_new<-cbind(S1_specific_GO[,2:3],log10pvalue_S1)
S1_specific_GO_new_psort<-S1_specific_GO_new[order(S1_specific_GO_new[["log10pvalue_S1"]],decreasing=T),][1:80,]
S1_specific_GO_sort<-S1_specific_GO_new_psort[order(S1_specific_GO_new_psort[["class"]],S1_specific_GO_new_psort[["log10pvalue_S1"]],decreasing=T),]
GOtype<-as.character(S1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_S1<-S1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_S1)+2))
pdf("S1 Specific Enriched GO Terms.pdf",w=16,h=10)
par(mar=c(15,5,5,5))
bar<-barplot(matrixp_S1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only S1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=S1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
S1_specific_GO_new_psort<-S1_specific_GO_new[order(S1_specific_GO_new[["log10pvalue_S1"]],decreasing=T),][1:100,]
S1_specific_GO_sort<-S1_specific_GO_new_psort[order(S1_specific_GO_new_psort[["class"]],S1_specific_GO_new_psort[["log10pvalue_S1"]],decreasing=T),]
GOtype<-as.character(S1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_S1<-S1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_S1)+2))
pdf("S1 Specific Enriched GO Terms.pdf",w=16,h=10)
par(mar=c(20,10,5,1))
bar<-barplot(matrixp_S1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only S1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=S1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
S1_specific_GO_new[1:20,1:10]
names(S1_specific_GO_new)
S1_specific_GO_new[1:10,]
S1_specific_GO_new[1:50,]
S1_specific_GO_new_bp<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="biological_process"),]
dim(S1_specific_GO_new_bp)
S1_specific_GO_new_bp<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="molecular_function"),]
# CC
S1_specific_GO_new_bp<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="cellular_component"),]
# BP
S1_specific_GO_new_bp<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="biological_process"),]
# MF
S1_specific_GO_new_mf<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="molecular_function"),]
# CC
S1_specific_GO_new_cc<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="cellular_component"),]
dim(S1_specific_GO_new_bp)
dim(S1_specific_GO_new_cc)
dim(S1_specific_GO_new_mf)
S1_specific_GO_new_bp[1:10,]
S1_specific_GO_new_bp_sort<-S1_specific_GO_new_bp[order(S1_specific_GO_new_bp[["log10pvalue_S1"]],decreasing=T),][1:30]
S1_specific_GO_new_bp[1:4,]
S1_specific_GO_new_bp[order(as.numeric(S1_specific_GO_new_bp[["log10pvalue_S1"]]),decreasing=T),][1:30]
as.numeric(S1_specific_GO_new_bp[["log10pvalue_S1"]]
)
S1_specific_GO_new_bp[order(as.numeric(S1_specific_GO_new_bp[["log10pvalue_S1"]]),decreasing=T),][1:30,]
S1_specific_GO_new_bp[order(as.numeric(S1_specific_GO_new_bp[["log10pvalue_S1"]]),decreasing=T),][,3]
# BP
S1_specific_GO_new_bp<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="biological_process"),]
S1_specific_GO_new_bp_sort<-S1_specific_GO_new_bp[order(as.numeric(S1_specific_GO_new_bp[["log10pvalue_S1"]]),decreasing=T),][1:30,]
# MF
S1_specific_GO_new_mf<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="molecular_function"),]
S1_specific_GO_new_mf_sort<-S1_specific_GO_new_mf[order(as.numeric(S1_specific_GO_new_mf[["log10pvalue_S1"]]),decreasing=T),][1:30,]
# CC
S1_specific_GO_new_cc<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="cellular_component"),]
S1_specific_GO_new_cc_sort<-S1_specific_GO_new_cc[order(as.numeric(S1_specific_GO_new_cc[["log10pvalue_S1"]]),decreasing=T),][1:30,]
S1_specific_GO_new_mf_sort
S1_specific_GO_new_cc_sort
# BP
S1_specific_GO_new_bp<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="biological_process"),]
S1_specific_GO_new_bp_sort<-S1_specific_GO_new_bp[order(as.numeric(S1_specific_GO_new_bp[["log10pvalue_S1"]]),decreasing=T),][1:20,]
# MF
S1_specific_GO_new_mf<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="molecular_function"),]
S1_specific_GO_new_mf_sort<-S1_specific_GO_new_mf[order(as.numeric(S1_specific_GO_new_mf[["log10pvalue_S1"]]),decreasing=T),][1:20,]
# CC
S1_specific_GO_new_cc<-S1_specific_GO_new[which(S1_specific_GO_new[["class"]]=="cellular_component"),]
S1_specific_GO_new_cc_sort<-S1_specific_GO_new_cc[order(as.numeric(S1_specific_GO_new_cc[["log10pvalue_S1"]]),decreasing=T),][1:20,]
S1_specific_GO_new_cc_sort
S1_specific_GO_new_mf_sort
S1_specific_GO_new_psort<-rbind(S1_specific_GO_new_bp_sort,S1_specific_GO_new_mf_sort,S1_specific_GO_new_cc_sort)
S1_specific_GO_new_psort
S1_specific_GO_new_psort)adim(
dim(S1_specific_GO_new_psort)
S1_specific_GO_new_psort<-S1_specific_GO_new1[order(S1_specific_GO_new1[["log10pvalue_S1"]],decreasing=T),]
S1_specific_GO_sort<-S1_specific_GO_new_psort[order(S1_specific_GO_new_psort[["class"]],S1_specific_GO_new_psort[["log10pvalue_S1"]],decreasing=T),]
S1_specific_GO_new1<-rbind(S1_specific_GO_new_bp_sort,S1_specific_GO_new_mf_sort,S1_specific_GO_new_cc_sort)
S1_specific_GO_new_psort<-S1_specific_GO_new1[order(S1_specific_GO_new1[["log10pvalue_S1"]],decreasing=T),]
S1_specific_GO_sort<-S1_specific_GO_new_psort[order(S1_specific_GO_new_psort[["class"]],S1_specific_GO_new_psort[["log10pvalue_S1"]],decreasing=T),]
S1_specific_GO_sort
GOtype<-as.character(S1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_S1<-S1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_S1)+2))
pdf("S1 Specific Enriched GO Terms.pdf",w=16,h=10)
par(mar=c(20,10,5,1))
bar<-barplot(matrixp_S1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only S1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=S1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
pdf("W1 Specific Enriched GO Terms.pdf",w=6,h=6)
par(mar=c(15,8,5,5))
bar<-barplot(matrixp_W1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=W1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
W1_specific_GO_sort
GOtype<-as.character(W1_specific_GO_sort[["class"]])
bpl<-which(GOtype=="biological_process")
ccl<-which(GOtype=="cellular_component")
mfl<-which(GOtype=="molecular_function")
mycol<-GOtype
mycol[bpl]<-"#6B8E23"
mycol[ccl]<-"#7EC0EE"
mycol[mfl]<-"#FF69B4"
matrixp_W1<-W1_specific_GO_sort[[3]]
ylims<-c(0,round(max(matrixp_W1)+2))
matrixp_W1
mycol
matrixp_W1
pdf("W1 Specific Enriched GO Terms.pdf",w=6,h=6)
par(mar=c(15,8,5,5))
bar<-barplot(matrixp_W1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=W1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
GOtype
bar<-barplot(matrixp_W1,beside=T,col=mycol,ylim=ylims)
mycol
plot(1:5,col=mycol)
pdf("W1 Specific Enriched GO Terms.pdf",w=8,h=6)
par(mar=c(15,8,5,5))
bar<-barplot(matrixp_W1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=W1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
dev.off()
dev.off()
dev.off()
pdf("W1 Specific Enriched GO Terms.pdf",w=8,h=6)
par(mar=c(15,8,5,5))
bar<-barplot(matrixp_W1,beside=T,col=mycol,ylim=ylims,axes=FALSE,axisnames=FALSE,ylab=paste("Enrichment Score : -log10(Pvalue)",sep=""),cex.lab=1,font.lab=1,col.lab="black",main="Enrichment GO Terms for only W1 Group")
axis(2,seq(0,ylims[2],by=1),as.character(seq(0,ylims[2],by=1)),cex.axis=1,font.axis=1,col.axis="black")
text(bar,rep(-0.2,length(bar)),labels=W1_specific_GO_sort[[1]],srt=45,adj=1,xpd=T,cex=1,font=1,col="black")
dev.off()
q()
co_GO[1:4,]
co_GO_sort
S1_GO[1:4,]
S1_GO_class<-merge(S1_GO,GOterm_class_id,by="id")
S1_GO_class[1:10,1:10]
W1_GO_class<-merge(W1_GO,GOterm_class_id,by="id")
co_GO<-merge(S1_GO_class,W1_GO_class,by="id")
co_GO[1:4,]
names(co_GO)
co_GO_short<-co_GO[,c(1,3,4,5,6,7,8,9,10,11,18:24)]
co_GO_short[1:10,]
write.table(co_GO_short,"co GOs for S1 and W1.xls",sep="\t",col.names=T,row.names=F,quote=F)
S1_specific_GO<-merge(S1_GOids,S1_GO_class,by="id")
S1_specific_GO[1:4,]
names(S1_specific_GO)
S1_specific_GO_short<-S1_specific_GO[,c(1,3,11,4:10)]
S1_specific_GO_short[1:4,]
write.table(S1_specific_GO_short,"W1_specific_GO.xls",sep="\t",col.names=T,row.names=F,quote=F)
W1_specific_GO<-merge(W1_GOids,W1_GO_class,by="id")
W1_specific_GO_short<-W1_specific_GO[,c(1,3,11,4:10)]
W1_specific_GO_short[1:4,]
write.table(W1_specific_GO_short,"W1_specific_GO.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(W1_specific_GO_short,"W1_specific_GO.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(S1_specific_GO_short,"S1_specific_GO.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
cc<-unique(as.character(read.delim("Colors")[[1]]))
cc
plot(,pch=15,col=cc,cex=2)
q()
cc<-read.delim("Colors",sep="\t",header=F)
cc
names(cc)
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]))
q()
cc<-read.delim("Colors",sep="\t",header=F)
cc
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]))
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=5)
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=5,xlim=c(0,10))
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=5,xlim=c(0,10),ylim=c(-2,12))
library(gplot)
library(gplots)
plot(1:10,pch=15,cex=2,col=redblue(20)[1:10])
redblue(20)[1:10]
ls()
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=5,xlim=c(0,10),ylim=c(-2,12))
cc
plot(1:10,pch=15,cex=2,col=redblue(20)[1:10])
cc[4,2]
cc[4,2]<-"#FFE3E3"
q()
cc<-read.delim("Colors",sep="\t",header=F)
cc
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=5,xlim=c(0,10),ylim=c(-2,12))
plot(1:10,pch=15,cex=5,col=c("#FF0000","#FF3939"))
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=5,xlim=c(0,10),ylim=c(-2,12))
plot(1:10,pch=15,cex=2,col=redblue(20)[1:10])
library(gplots)
plot(1:10,pch=15,cex=2,col=redblue(20)[1:10])
redblue(20)[1:10]
q()
cc<-read.delim("Colors",sep="\t",header=F)
cc
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=5,xlim=c(0,10),ylim=c(-2,12))
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=3,xlim=c(0,10),ylim=c(-2,12))
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=4,xlim=c(0,10),ylim=c(-2,12))
text(x=rep(4,9),y=1:9,lables=as.character(cc[[1]]))
as.character(cc[[1]])
text(x=rep(4,9),y=1:9,as.character(cc[[1]]))
text(x=rep(4,9),y=1:9,as.character(cc[[1]]),adj=1)
text(x=rep(4,9),y=1:9,as.character(cc[[1]]),adj=0)
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=4,xlim=c(0,10),ylim=c(-2,12))
text(x=rep(4,9),y=1:9,as.character(cc[[1]]),adj=0)
text(x=rep(3,9),y=1:9,as.character(cc[[1]]),adj=0)
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=4,xlim=c(0,10),ylim=c(-2,12))
text(x=rep(3,9),y=1:9,as.character(cc[[1]]),adj=0)
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=4,xlim=c(0,6),ylim=c(-2,12))
text(x=rep(3,9),y=1:9,as.character(cc[[1]]),adj=0)
pdf("ColorCard.pdf")
plot(x=rep(2,9),y=1:9,col=as.character(cc[[2]]),pch=15,cex=4,xlim=c(0,6),ylim=c(-2,12))
text(x=rep(3,9),y=1:9,as.character(cc[[1]]),adj=0)
dev.off()
q()
dir()
ENSG2Ncbi<-read.delim("ENSG2Ncbi.txt.1",header=F,sep="\t")
ncbiID2rnoID<-read.delim("ncbiID2rnoID",header=F,sep="\t")
ENSG2Ncbi[1:4,]
ncbiID2rnoID[1:4,]
names(ENSG2Ncbi)<-c("ENSG","Ncbi")
names(ncbiID2rnoID)<-c("Ncbi","rnoID")
ENSG2Ncbi[1:4,]
ncbiID2rnoID[1:4,]
merge(ENSG2Ncbi,ncbiID2rnoID,by="Ncbi")->ENSG2Ncbi2rnoID
ENSG2Ncbi2rnoID[1:4,]
ENSG2Ncbi2rnoID_uniq<-unique(ENSG2Ncbi2rnoID)
write.table(ENSG2Ncbi2rnoID_uniq,"ENSG2Ncbi2rnoID_uniq",sep="\t",col.names=T,row.names=F,quote=F)
q()
diff_gene_list<-read.delim("**")
diff_gene_list<-read.delim("TRM_vs_TRMS_DEisoform.list")
names(diff_gene_list)<-"gene_id"
all_gene_COG<-read.delim
diff_gene_COG<-merge
("TRM_vs_TRMS_DEisoform.list",header=T,sep="\t",check.names=F)
q90
q()
diff_gene_list<-read.delim("TRM_vs_TRMS_DEisoform.list")
names(diff_gene_list)<-"gene_id"
all_gene_COG<-read.delim("gene_id,COG_ID,COG_Description",header=T,sep="\t",check.names=F)
all_gene_COG<-read.delim("T341_COG.txt",header=T,sep="\t",check.names=F)
names(all_gene_COG)<-c("gene_id","COG_ID","COG_Description")
diff_gene_COG<-merge(diff_gene_list,all_gene_COG,by="gene_id")
T1<-table(diff_gene_COG[["COG_Description"]])
COG2Num<-cbind(names(T1),T1)
pdf(file="pie_COG.pdf",height=20,width=28)
slices<-as.numeric(COG2Num[[2]])
slices<-as.numeric(COG2Num[[2]])
lbls<-as.character(COG2Num[[1]])
COG2Num[1:2,]
GETWD()
getwd()
all_gene_COG<-read.delim("T341_COG.txt",header=T,sep="\t",check.names=F)
all_gene_COG[1:4,]
names(all_gene_COG)<-c("gene_id","COG_ID","COG_Description")
diff_gene_COG<-merge(diff_gene_list,all_gene_COG,by="gene_id")
T1<-table(diff_gene_COG[["COG_Description"]])
COG2Num<-cbind(names(T1),T1)
COG2Num
dim(COG2Num)
COG2Num[1:4,]
names(all_gene_COG)
names(diff_gene_list)
diff_gene_list[1:4,]
names(diff_gene_COG)
dim(diff_gene_COG)
diff_gene_COG[1:4,]
T1<-table(diff_gene_COG[["COG_Description"]])
T1
COG2Num<-cbind(names(T1),T1)
COG2Num
dim(COG2Num)
slices<-as.numeric(COG2Num[[2]])
lbls<-as.character(COG2Num[[1]])
as.numeric(COG2Num[[2]])
COG2Num
dim(COG2Num)
COG2Num[,2]
as.numeric(COG2Num[,2])
slices<-as.numeric(COG2Num[,2])
lbls<-as.character(COG2Num[,1])
pdf(file="pie_COG.pdf",height=20,width=28)
pct<-round(slices/sum(slices)*100,2)
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(slices,cex=1.5,labels=lbls,col=rainbow(length(lbls)),main="Pie Chart of KEGG pathways")
dev.off()
history()
history
ls
q()
source("Rfam.R")
q()
source("Rfam.R")
q()
source("Rfam.R")
q()
dir()
pathway_table<-read.delim("pathway_table.xls",header=T,sep=:)
pathway_table<-read.delim("pathway_table.xls",header=T,sep="\t",check.names=T)
SDEG<-read.delim("S_0-5M_vs_S_1-5M_DEisoform",header=F)
pathway_table[1:4,]
q()
pathway_table<-read.delim("pathway_table.xls.new",header=T,sep="\t",check.names=T)
pathway_table[1:4,]
ls()
d

SDEG<-scan("S_0-5M_vs_S_1-5M_DEisoform",what=character())
SDEG
ls()
pathway_table[1:4,]
SDEG
pathway_table[1:4,]
pathway_table_new<-data.frame()
for(i in 1:nrow(pathway_table)){
pathway_table_new[i,1]<-pathway_table[i,1]
pathway_table_new[i,2]<-pathway_table[i,2]
pathway_table_new[i,3]<-pathway_table[i,3]
pathway_table_new[i,4]<-pathway_table[i,4]
all_genes<-unlist(strsplit(as.character(pathway_table[i,4]),";",fix=T))
diff_genes<-intersect(all_genes,SDEG)
pathway_table_new[i,5]<-paste(diff_genes,collapse=";")
pathway_table_new[i,6]<-length(diff_genes)
}
pathway_table_new[1:4,]
pathway_table[1:4,]
names(pathway_table_new)
names(pathway_table)
names(pathway_table_new)
names(pathway_table_new)<-c("PathWay","Pathway_definition","number_of_seqs","seqs_kos.genes_list","seqs_kos.diffgenes.list","number_of_diffseqs")
pathway_table_new_new<-pathway_table_new[,c(1:2,6,5)]
names(pathway_table_new_new)
names(pathway_table_new_new)[3:4]<-c("number_of_seqs","seqs_kos.genes.list")
write.table(pathway_table_new_new,"pathway_table_sdeg.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
isoformID_geneName<-read.delim("isoformID_geneName",sep="\t",header=F)
top20path<-read.delim("170_diff_genes_top20_pathway.xls",sep="\t",header=T)
top20path[1:4,]
genename<-vector()
isoformID_geneName[1:4,]
for(i in 1:nrow(top20path)){
genelist<-unlist(strsplit(as.character(top20path[i,"Genes"]),"|",fix=T))
}
i=1
genelist<-unlist(strsplit(as.character(top20path[i,"Genes"]),"|",fix=T))
genelist
genelist<-as.data.frame(unlist(strsplit(as.character(top20path[i,"Genes"]),"|",fix=T)))
genelist
names(genelist)<-"V1"
genelist
merge(genelist,isoformID_geneName,by="V1")
?merge
merge(genelist,isoformID_geneName,by="V1",sort=F)
genelist
as.character(genelist)
as.character(genelist[[1]])
merge(genelist,isoformID_geneName,by="V1",sort=F)->ll
as.character(ll[[2]])
history
history()
history()
genename<-vector()
for(i in 1:nrow(top20path)){
genelist<-as.data.frame(unlist(strsplit(as.character(top20path[i,"Genes"]),"|",fix=T)))
names(genelist)<-"V1"
aa<-as.character(merge(genelist,isoformID_geneName,by="V1",sort=F)[[2]])
genename[i]<-paste(aa,collapse="|")
}
genename
top20path
genename
cbind(top20path,genename)->ll
ll[1:4,]
write.table(ll,"top20path.xls",sep="\t",col.names=T,row.names=T,quote=F)
write.table(ll,"top20path.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
kegg<-read.delim("170_diff_genes_kegg_enrichment.xls",header=T,sep="\t")
kegg[1:4,]
history()
kegg<-read.delim("170_diff_genes_kegg_enrichment.xls",header=T,sep="\t")
isoformID_geneName<-read.delim("isoformID_geneName",header=T,sep="\t")
GeneName<-vector()
for(i in 1:nrow(kegg)){
genelist<-as.data.frame(unlist(strsplit(as.character(kegg[i,"Genes"]),"|",fix=T)))
names(genelist)<-"V1"
aa<-as.character(merge(genelist,isoformID_geneName,by="V1",sort=F)[[2]])
GeneName[i]<-paste(aa,collapse="|")
}
kegg_new<-cbind(kegg,GeneName)
i
kegg<-read.delim("170_diff_genes_kegg_enrichment.xls",header=T,sep="\t")
isoformID_geneName<-read.delim("isoformID_geneName",header=T,sep="\t")
GeneName<-vector()
genelist<-as.data.frame(unlist(strsplit(as.character(kegg[i,"Genes"]),"|",fix=T)))
names(genelist)<-"V1"
aa<-as.character(merge(genelist,isoformID_geneName,by="V1",sort=F)[[2]])
GeneName[i]<-paste(aa,collapse="|")
genelist
isoformID_geneName[1:4,]
kegg<-read.delim("170_diff_genes_kegg_enrichment.xls",header=T,sep="\t")
isoformID_geneName<-read.delim("isoformID_geneName",header=F,sep="\t")
GeneName<-vector()
for(i in 1:nrow(kegg)){
genelist<-as.data.frame(unlist(strsplit(as.character(kegg[i,"Genes"]),"|",fix=T)))
names(genelist)<-"V1"
aa<-as.character(merge(genelist,isoformID_geneName,by="V1",sort=F)[[2]])
GeneName[i]<-paste(aa,collapse="|")
}
kegg_new<-cbind(kegg,GeneName)
kegg_new[1:4,]
dim(kegg_new)
kegg_new<-cbind(kegg,GeneName)[,c(1:8,10,9)]
kegg_new[1:4,]
write.table(kegg_new,"170_diff_genes_kegg_enrichment_addname.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("13069_vs_49.sigDEgene.fpkm.new", header=T, com='', sep="\t")
rownames(data) = data[,1] 
data = data[,2:length(data[1,])] 
data = as.matrix(data) 
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) 
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") 
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=25,height=26, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=2, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(15,20))
dev.off()
q()
library(cluster)
library(gplots)
library(Biobase)
load("all.RData")
outdir = "clusters.K_6"
dir.create(outdir)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6)
gene_colors = partition_colors[gene_partition_assignments]
pdf(file="union.DEgene.clusters.K_6.heatmap.pdf", width=25, height=26, paper="special");
heatmap.2(centered_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE,keysize=1.2,cexRow=0.5, cexCol=3,margins=c(15,15),lhei=c(0.3,1.5), lwid=c(2.5,5))
dev.off()
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("13069_vs_49.sigDEgene.fpkm.new", header=T, com='', sep="\t")
rownames(data) = data[,1] 
data = data[,2:length(data[1,])] 
data = as.matrix(data) 
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) 
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") 
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=30,height=30, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexRow=0.5, cexCol=2, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(15,20))
dev.off()
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("13069_vs_49.sigDEgene.fpkm.new", header=T, com='', sep="\t")
rownames(data) = data[,1] 
data = data[,2:length(data[1,])] 
data = as.matrix(data) 
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) 
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean") 
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=30,height=30, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexRow=0.8, cexCol=3, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(15,20))
dev.off()
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=3, cexRow=1, lhei=c(0.3,2), lwid=c(2.5,4),margins=c(15,20))
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("13069_vs_49.sigDEgene.fpkm.new", header=T, com='', sep="\t")
rownames(data) = data[,1] 
data = data[,2:length(data[1,])] 
data = as.matrix(data) 
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) 
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean")
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=30,height=30, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=3, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(15,20))
dev.off()
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=30,height=35, paper="special");
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=40,height=45, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(15,20))
dev.off()
?heatmap.2
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(15,20),las=1)
dev.off()
?las
??las
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("13069_vs_49.sigDEgene.fpkm.new", header=T, com='', sep="\t")
rownames(data) = data[,1] 
data = data[,2:length(data[1,])] 
data = as.matrix(data) 
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) 
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean")
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=30,height=30, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=3, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(20,20))
dev.off()
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=35,height=40, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=4, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(20,20))
dev.off()
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=45,height=50, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=4, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,4),margins=c(30,15))
dev.off()
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=5, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(32,30))
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=5, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(30,30))
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=5, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(30,25))
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=5, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(30,20))
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=5, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(30,15))
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=5, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(25,20))
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", cexCol=5, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(20,20))
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("13069_vs_49.sigDEgene.fpkm.new", header=T, com='', sep="\t")
rownames(data) = data[,1] 
data = data[,2:length(data[1,])] 
data = as.matrix(data) 
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) 
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean")
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=30,height=40, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(25,20))
dev.off()
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,2), lwid=c(2.5,3),margins=c(25,20),lhei=c(0.3,1.5), lwid=c(2.5,5))
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,1), lwid=c(2.5,3),margins=c(25,20))
dev.off()
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=25,height=45, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,1), lwid=c(2.5,3),margins=c(25,20))
dev.off()
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,0.5), lwid=c(2.5,4),margins=c(25,20))
dev.off()
q()
library(cluster)
library(gplots)
library(Biobase)
data = read.table("13069_vs_49.sigDEgene.fpkm.new", header=T, com='', sep="\t")
rownames(data) = data[,1] 
data = data[,2:length(data[1,])] 
data = as.matrix(data) 
data = log2(data+1)
centered_data = t(scale(t(data), scale=F)) 
hc_genes = agnes(centered_data, diss=FALSE, metric="euclidean")
hc_samples = hclust(as.dist(1-cor(centered_data, method="spearman")), method="complete") 
myheatcol = redgreen(75)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=10);
partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]
save(list=ls(all=TRUE), file="all.RData")
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=30,height=40, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,1), lwid=c(2.5,3),margins=c(25,20))
dev.off()
pdf(file="13069_vs_49.sigDEgene.heatmap.pdf", width=30,height=35, paper="special");
heatmap.2(centered_data, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=NA, col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none",cexCol=4, cexRow=1,lhei=c(0.3,1), lwid=c(2.5,3),margins=c(25,20))
dev.off()
q()
library("DEGseq")
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(2))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(3))
outputDir <- file.path("13069_vs_49")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir)
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(4))
outputDir <- file.path("13069_vs_DK")
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(4))
outputDir <- file.path("13069_vs_DK")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="DK",method="MARS",outputDir=outputDir)
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(3))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(4))
outputDir <- file.path("49_vs_DK")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="49",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="DK",method="MARS",outputDir=outputDir)
q()
library("DEGseq")
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(2))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(3))
outputDir <- file.path("13069_vs_49_0508") 
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,pValue=1e-3,zScore=4, qValue=1e-3, foldChange=4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,pValue=1e-3,zScore=4,qValue=5e-4,foldChange=4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,pValue=1e-3,zScore=4,qValue=1e-3,foldChange=4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,pValue=1e-3,zScore=4,qValue=1e-3,foldChange=4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,pValue=1e-3,zScore=4,qValue=1e-3,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=1e-3,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-3,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-3,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(4))
outputDir <- file.path("13069_vs_DK_0508")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="DK",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(3))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(4))
outputDir <- file.path("49_vs_DK_0508")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="49",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="DK",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
q()
?DEGseq
?DEGexp
??DEGexp
?DEGseq::DEGexp
q()
library("DEGseq")
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(2,3,4))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(2,3,4))
outputDir <- file.path("test") 
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2,3,4),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2,3,4),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
q()
library("DEGseq")
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(2))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(3))
outputDir <- file.path("13069_vs_49.0516")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=5e-2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-3)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=1)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=1)
?DEGseq::DEGexp
q()
library("DEGseq")
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(2))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(3))
outputDir <- file.path("13069_vs_49.0516")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=1)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1 =1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=3)
?DEGseq::DEGexp
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method=c("MARS","FC"),outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=1)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=1)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-5)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=1e-3,,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=105)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=1.5)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=1)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-5)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=1e-3)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=5,qValue=5e-2,foldChange=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=3)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=2)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=1)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-4)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-5)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-5)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-6)
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1, valCol=c(4))
outputDir <- file.path("13069_vs_DK")
outputDir <- file.path("13069_vs_DK.0516")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="DK",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-6)
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(3))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(4))
outputDir <- file.path("49_vs_DK.0516")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="49",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="DK",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-6)
q()
library("DEGseq")
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(2))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(3))
outputDir <- file.path("13069_vs_49.MATR") 
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MATR",outputDir=outputDir,thresholdKind=4,qValue=1e-3)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="LRT",outputDir=outputDir,thresholdKind=4,qValue=1e-3)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="FET",outputDir=outputDir,thresholdKind=4,qValue=1e-3)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-3)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MATR",outputDir=outputDir,thresholdKind=4,qValue=1e-3)
?DEGseq::DEGexp
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MATR",replicateExpMatrix1=NULL,replicateExpMatrix2=NULL,outputDir=outputDir,thresholdKind=4,qValue=1e-3)
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="LRT",outputDir=outputDir,thresholdKind=4,qValue=1e-6)
outputDir <- file.path("13069_vs_49.0517") 
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-8)
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(2))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(3))
outputDir <- file.path("13069_vs_49.0516") 
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="49",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-8)
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(4))
outputDir <- file.path("13069_vs_DK.0516") 
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="13069",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="DK",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-8)
geneExpMatrix1 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(3))
geneExpMatrix2 <- readGeneExp(file="13069-vs-49-vs-DK.gene.exp_uniq.xls",geneCol=1,valCol=c(4))
outputDir <- file.path("49_vs_DK.0516") 
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2),groupLabel1="49",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2),groupLabel2="DK",method="MARS",outputDir=outputDir,thresholdKind=4,qValue=1e-8)
q()
??DEGseq
?DEGseq::DEGseq
q()
spl5<-read.delim("spl5.junctions.bed.novel.classify.newfilter.xls",header=T,sep="\t",check.names=T)
spl5[1:4,]
splice_type<-unique(as.character(spl5[["splice_type"]]))
splice_type
q()
spl5<-read.delim("spl5.junctions.bed.novel.classify.filtered.new",header=T,sep="\t",check.names=T)
splice_type<-unique(as.character(spl5[["splice_type"]]))
splice_type
spl5[1:4,]
delLines<-unique(c(which(spl5[["splice_type"]]=="Other"),which(spl5[["splice_type"]]=="Intergenic")))
delLines
## 一个基因发生两种或两种以上的个数统计
spl5_raw<-read.delim("spl5.junctions.bed.novel.classify.filtered",header=T,sep="\t",check.names=T)
## del Other and Intergenic
OtherLines<-which(spl5[["splice_type"]]=="Other")
IntergenicLines<-which(spl5[["splice_type"]]=="Intergenic")
delLines<-unique(c(OtherLines,IntergenicLines))
spl5<-spl5_raw[-delLines,]
splice_type<-unique(as.character(spl5[["splice_type"]]))
dim(spl5)
spl5[1:4,]
spl5[1:40,]
OtherLines
spl5_raw<-read.delim("spl5.junctions.bed.novel.classify.filtered",header=T,sep="\t",check.names=T)
## del Other and Intergenic
OtherLines<-which(spl5_raw[["splice_type"]]=="Other")
IntergenicLines<-which(spl5_raw[["splice_type"]]=="Intergenic")
delLines<-unique(c(OtherLines,IntergenicLines))
spl5<-spl5_raw[-delLines,]
splice_type<-unique(as.character(spl5[["splice_type"]]))
spl5[1:40,]
SpliceSummary<-data.frame()
Genes<-unique(as.character(spl5[["gene_id"]]))
   SpliceSummary[i,1]<-Genes[i]
i=1
   SpliceSummary[i,1]<-Genes[i]
SpliceSummary
   nline<-which(spl5[["gene_id"]]==Genes[i])
nline
 spl5[nline,]
  SpliceSummary[i,1]<-Genes[i]
   nline<-which(spl5[["gene_id"]]==Genes[i])
   Stype<-as.character(spl5[nline,"splice_type"])
   
<Stype
Stype
"5S"%inStype
"5S"%in%Stype
grep("5S",Stype)
grep("5sS",Stype)
grep("5sS",Stype)
length(grep("5sS",Stype))
length(grep("5S",Stype))
SpliceSummary[i,1]<-Genes[i]
   nline<-which(spl5[["gene_id"]]==Genes[i])
   Stype<-as.character(spl5[nline,"splice_type"])
   SpliceSummary[i,2]<-length(grep("A3SS",Stype))
   SpliceSummary[i,3]<-length(grep("3UTR",Stype))
   SpliceSummary[i,4]<-length(grep("A5SS",Stype))
   SpliceSummary[i,5]<-length(grep("5UTR",Stype))
   SpliceSummary[i,6]<-length(grep("SE",Stype))
   SpliceSummary[i,7]<-length(grep("IR",Stype))
   
SpliceSummary
Stype
   SpliceSummary[i,1]<-Genes[i]
   nline<-which(spl5[["gene_id"]]==Genes[i])
   Stype<-as.character(spl5[nline,"splice_type"])
   SpliceSummary[i,2]<-length(grep("3S",Stype))
   SpliceSummary[i,3]<-length(grep("3UTR",Stype))
   SpliceSummary[i,4]<-length(grep("5S",Stype))
   SpliceSummary[i,5]<-length(grep("5UTR",Stype))
   SpliceSummary[i,6]<-length(grep("ES",Stype))
   SpliceSummary[i,7]<-length(grep("IR",Stype))
   
   
SpliceSummary
 which(SpliceSummary[i,]==1)
   
names(which(SpliceSummary[i,]==1))
colname<-c("3S","3UTR","5S","5UTR","ES","IR")
  colname[which(SpliceSummary[i,]==1)]
   
      colname[which(SpliceSummary[i,]==1)+1,]
   
   colname[which(SpliceSummary[i,]==1)+1]
colname<-c("gene_id","3S","3UTR","5S","5UTR","ES","IR","Num","Type")
   colname[which(SpliceSummary[i,]==1),]
   
   colname[which(SpliceSummary[i,]==1)]
Stype
 paste(colname[which(SpliceSummary[i,]==1)],collapse="-")
   
which(SpliceSummary[i,]==1)
SpliceSummary<-data.frame()
colname<-c("gene_id","3S","3UTR","5S","5UTR","ES","IR","Num","Type")
Genes<-unique(as.character(spl5[["gene_id"]]))
for(i in 1:length(Genes)){
   SpliceSummary[i,1]<-Genes[i]
   nline<-which(spl5[["gene_id"]]==Genes[i])
   Stype<-as.character(spl5[nline,"splice_type"])
   SpliceSummary[i,2]<-length(grep("3S",Stype))
   SpliceSummary[i,3]<-length(grep("3UTR",Stype))
   SpliceSummary[i,4]<-length(grep("5S",Stype))
   SpliceSummary[i,5]<-length(grep("5UTR",Stype))
   SpliceSummary[i,6]<-length(grep("ES",Stype))
   SpliceSummary[i,7]<-length(grep("IR",Stype))
   SpliceSummary[i,8]<-paste(colname[which(SpliceSummary[i,]==1)],collapse="-")
   SpliceSummary[i,9]<-length(which(SpliceSummary[i,]==1))
}
names(SpliceSummary)<-colname
q()
SpliceSummary[1:10,]
i
i=4
Genes[4]
nline<-which(spl5[["gene_id"]]==Genes[i])
nline
   Stype<-as.character(spl5[nline,"splice_type"])
Stype
paste(colname[which(SpliceSummary[i,]==1)],collapse="-")
paste(colname[which(SpliceSummary[i,]>=1)],collapse="-")
SpliceSummary[i,]
which(SpliceSummary[i,]>=1)
which(SpliceSummary[i,2:7]>=1)
which(SpliceSummary[i,2:7]>=1)
i=1
which(SpliceSummary[i,2:7]>=1)
SpliceSummary[1,]
which(SpliceSummary[1,]>=1)
which(SpliceSummary[1,2:7]>=1)
which(SpliceSummary[i,2:7]>=1)
which(SpliceSummary[i,2:7]>=1)+1
i=4
paste(colname[(which(SpliceSummary[i,2:7]>=1)+1)],collapse="-")
length(which(SpliceSummary[i,]==1))
length(which(SpliceSummary[i,]>=1))
SpliceSummary[i,]
SpliceSummary[1:10,]
i
   SpliceSummary[i,1]<-Genes[i]
   nline<-which(spl5[["gene_id"]]==Genes[i])
   Stype<-as.character(spl5[nline,"splice_type"])
   SpliceSummary[i,2]<-length(grep("3S",Stype))
   SpliceSummary[i,3]<-length(grep("3UTR",Stype))
   SpliceSummary[i,4]<-length(grep("5S",Stype))
   SpliceSummary[i,5]<-length(grep("5UTR",Stype))
   SpliceSummary[i,6]<-length(grep("ES",Stype))
   SpliceSummary[i,7]<-length(grep("IR",Stype))
   SpliceSummary[i,8]<-length(which(SpliceSummary[i,]==1))
   SpliceSummary[i,9]<-paste(colname[(which(SpliceSummary[i,2:7]>=1)+1)],collapse="-")
SpliceSummary[i,]
length(which(SpliceSummary[i,]>=1))
   SpliceSummary[i,2]<-length(grep("3S",Stype))
   SpliceSummary[i,3]<-length(grep("3UTR",Stype))
   SpliceSummary[i,4]<-length(grep("5S",Stype))
   SpliceSummary[i,5]<-length(grep("5UTR",Stype))
   SpliceSummary[i,6]<-length(grep("ES",Stype))
   SpliceSummary[i,7]<-length(grep("IR",Stype))
   SpliceSummary[i,8]<-length(which(SpliceSummary[i,2:7]>=1))
   SpliceSummary[i,9]<-paste(colname[(which(SpliceSummary[i,2:7]>=1)+1)],collapse="-")
 SpliceSummary[i,]
sum(SpliceSummary[i,2:7])
Genes<-unique(as.character(spl5[["gene_id"]]))
for(i in 1:length(Genes)){
   SpliceSummary[i,1]<-Genes[i]
   nline<-which(spl5[["gene_id"]]==Genes[i])
   Stype<-as.character(spl5[nline,"splice_type"])
   SpliceSummary[i,2]<-length(grep("3S",Stype))
   SpliceSummary[i,3]<-length(grep("3UTR",Stype))
   SpliceSummary[i,4]<-length(grep("5S",Stype))
   SpliceSummary[i,5]<-length(grep("5UTR",Stype))
   SpliceSummary[i,6]<-length(grep("ES",Stype))
   SpliceSummary[i,7]<-length(grep("IR",Stype))
   SpliceSummary[i,8]<-sum(SpliceSummary[i,2:7])
   SpliceSummary[i,9]<-paste(colname[(which(SpliceSummary[i,2:7]>=1)+1)],collapse="-")
}
names(SpliceSummary)<-colname
SpliceSummary[1:4,]
SpliceSummary[1:40,]
unique(as.character(SpliceSummary[[9]]))
sort(unique(as.character(SpliceSummary[[9]])))
SpliceSummary[1:40,]
j=1
     length(strsplit(as.character(SpliceSummary[j,"Type"]),"-",fix=T))
SpliceSummary[1:10,]
as.character(SpliceSummary[j,"Type"])
strsplit(as.character(SpliceSummary[j,"Type"]),"-",fix=T)
 length(unlist(strsplit(as.character(SpliceSummary[j,"Type"]),"-",fix=T)))
    
SpliceSummary1<-SpliceSummary
for(j in 1:nrow(SpliceSummary)){
    SpliceSummary1[j,10]<-length(unlist(strsplit(as.character(SpliceSummary[j,"Type"]),"-",fix=T)))
}
SpliceSummary1[1:10,]
SpliceSummary1[1:40,]
names(SpliceSummary1)[c(8,10)]<-c("ASEventNum","ASTypeNum")
SpliceSummary1[1:40,]
SpliceSummary2<-SpliceSummary1[,c(1:8,10,9)]
SpliceSummary2[1:40,]
dim(SpliceSummary2)
SpliceSummary2[1:40,]
SpliceSummary2[1:20,]
ASTypeNum<-unique(as.numeric(SpliceSummary2[["ASTypeNum"]]))
ASTypeNum
ASTypeNum<-sort(unique(as.numeric(SpliceSummary2[["ASTypeNum"]])))
ASTypeNum
ASTypeNum<-sort(unique(as.numeric(SpliceSummary2[["ASTypeNum"]])))
sapply(ASTypeNum,function(x) length(which(SpliceSummary2[["ASTypeNum"]]==x)))
ASType<-sort(unique(as.character(SpliceSummary2[["Type"]])))
ASType
ASType<-sort(unique(as.character(SpliceSummary2[["Type"]])))
AS54Num<-sapply(ASType,function(x) length(which(SpliceSummary2[["Type"]]==x)))
AS54Num
ASType<-sort(unique(as.character(SpliceSummary2[["Type"]])))
AS54Num<-sapply(ASType,function(x) length(which(SpliceSummary2[["Type"]]==x)))
ASType_GeneNum<-cbind(names(AS54Num),AS54Num)
ASType_GeneNum
cbind(names(AS54Num),as.numeric(AS54Num))
dim(SpliceSummary2)
ASType_GeneNum[[2]]
ASType_GeneNum[,2]
as.numeric(ASType_GeneNum[,2])
sapply(as.numeric(ASType_GeneNum[,2]),function(x) x/dim(SpliceSummary2)[1])
sapply(as.numeric(ASType_GeneNum[,2]),function(x) paste(round(100*x/dim(SpliceSummary2)[1],2),"%",sep=""))
sapply(as.numeric(ASType_GeneNum[,2]),function(x) round(100*x/dim(SpliceSummary2)[1],2))
ratio<-sapply(as.numeric(ASType_GeneNum[,2]),function(x) round(100*x/dim(SpliceSummary2)[1],2))
ASType_GeneNum_Ratio<-cbind(ASType_GeneNum,ratio)
ASType_GeneNum_Ratio
ASType_GeneNum_Ratio_sort<-ASType_GeneNum_Ratio[order(as.numeric(ASType_GeneNum_Ratio[,3]),decreasing=T),]
ASType_GeneNum_Ratio_sort
ASType_GeneNum_Ratio_sort<-ASType_GeneNum_Ratio[order(as.numeric(ASType_GeneNum_Ratio[,2]),decreasing=T),]
ASType_GeneNum_Ratio_sort
write.table(ASType_GeneNum_Ratio_sort,"ASType_GeneNum_Ratio_sort.xls",sep="\t",col.names=T,row.names=F,quote=F)
write.table(SpliceSummary2,"SpliceSummaryForEachGene",sep="\t",col.names=T,row.names=F,quote=F)
SixNum
ASTypeNum<-sort(unique(as.numeric(SpliceSummary2[["ASTypeNum"]])))
SixNum<-sapply(ASTypeNum,function(x) length(which(SpliceSummary2[["ASTypeNum"]]==x)))
SixNum
barplot(SixNum,col=heat.colors(6))
round(max(SixNum)/10)
x<-barplot(SixNum,col=heat.colors(6),ylim=c(0,max(SixNum)+round(max(SixNum)/10)))
?text
text(x,c(SixNum+round(max(SixNum)/10)),labels=1:6)
x<-barplot(SixNum,col=heat.colors(6),ylim=c(0,max(SixNum)+round(max(SixNum)/10)))
text(x,c(SixNum+round(max(SixNum)/10)),labels=SixNum)
text(x,rep(-5,6),labels=SixNum)
?text
text(x,rep(-5,6),labels=SixNum,adj=1)
text(x,rep(-5,6),labels=SixNum,adj=0)
text(x,rep(100,6),labels=SixNum,adj=0)
text(x,rep(-100,6),labels=SixNum,adj=0)
skip<-round(max(SixNum)/10)
x<-barplot(SixNum,col=heat.colors(6),ylim=c(0,max(SixNum)+skip))
text(x,c(SixNum+skip),labels=SixNum)
text(x,rep(-skip,6),labels=SixNum,adj=1)
x<-barplot(SixNum,col=heat.colors(6),ylim=c(0,max(SixNum)+skip))
text(x,c(SixNum+skip/2),labels=SixNum)
text(x,rep(-skip/2,6),labels=SixNum,adj=1)
SixNum
text(x,rep(-skip/2,6),labels=1:6,adj=1)
rep(-skip/2,6)
x
?text
?text
?text
?text
x<-barplot(SixNum,col=heat.colors(6),ylim=c(0,max(SixNum)+skip))
text(x,c(SixNum+skip/2),labels=SixNum)
text(x,rep(-skip/2,6),labels=1:6,adj=1,xpd=T)
x<-barplot(SixNum,col=heat.colors(6),ylim=c(0,max(SixNum)+skip),xlab="Number of AS Type",ylab="Number of Gene")
text(x,c(SixNum+skip/2),labels=SixNum)
text(x,rep(-skip/2,6),labels=1:6,adj=1,xpd=T)
pdf(file="Distribution of alternative splicing type number.pdf",w=10,h=9)
x<-barplot(SixNum,col=heat.colors(6),ylim=c(0,max(SixNum)+skip),xlab="Number of AS Type",ylab="Number of Gene")
text(x,c(SixNum+skip/2),labels=SixNum)
text(x,rep(-skip/2,6),labels=1:6,adj=1,xpd=T)
dev.off()
q()
write.table(SpliceSummary2,"SpliceSummaryForEachGene",sep="\t",col.names=T,row.names=F,quote=F)
q()
spl5_raw<-read.delim("spl5_ox.junctions.bed.novel.classify.filtered",header=T,sep="\t",check.names=T)
## del "Other" and "Intergenic"
OtherLines<-which(spl5_raw[["splice_type"]]=="Other")
IntergenicLines<-which(spl5_raw[["splice_type"]]=="Intergenic")
delLines<-unique(c(OtherLines,IntergenicLines))
spl5<-spl5_raw[-delLines,]
splice_type<-unique(as.character(spl5[["splice_type"]]))
## 哪几种类型的可变剪接在这个基因上发生
##  genename3S3UTR5S5UTRESIRNumType
##  gene110110003A3SS-A5SS-5UTR
##  gene210000001A3SS
##  ...
SpliceSummary<-data.frame()
colname<-c("gene_id","3S","3UTR","5S","5UTR","ES","IR","Num","Type")
Genes<-unique(as.character(spl5[["gene_id"]]))
for(i in 1:length(Genes)){
   SpliceSummary[i,1]<-Genes[i]
   nline<-which(spl5[["gene_id"]]==Genes[i])
   Stype<-as.character(spl5[nline,"splice_type"])
   SpliceSummary[i,2]<-length(grep("3S",Stype))
   SpliceSummary[i,3]<-length(grep("3UTR",Stype))
   SpliceSummary[i,4]<-length(grep("5S",Stype))
   SpliceSummary[i,5]<-length(grep("5UTR",Stype))
   SpliceSummary[i,6]<-length(grep("ES",Stype))
   SpliceSummary[i,7]<-length(grep("IR",Stype))
   SpliceSummary[i,8]<-sum(SpliceSummary[i,2:7])
   SpliceSummary[i,9]<-paste(colname[(which(SpliceSummary[i,2:7]>=1)+1)],collapse="-")
}
names(SpliceSummary)<-colname
SpliceSummary1<-SpliceSummary
for(j in 1:nrow(SpliceSummary)){
    SpliceSummary1[j,10]<-length(unlist(strsplit(as.character(SpliceSummary[j,"Type"]),"-",fix=T)))
}
names(SpliceSummary1)[c(8,10)]<-c("ASEventNum","ASTypeNum")
SpliceSummary2<-SpliceSummary1[,c(1:8,10,9)]
write.table(SpliceSummary2,"Spl5_ox_SpliceSummaryForEachGene",sep="\t",col.names=T,row.names=F,quote=F)
## 发生1种可变剪切数目、发生2种可变剪切数目、...
ASTypeNum<-sort(unique(as.numeric(SpliceSummary2[["ASTypeNum"]])))
SixNum<-sapply(ASTypeNum,function(x) length(which(SpliceSummary2[["ASTypeNum"]]==x)))
skip<-round(max(SixNum)/10)
pdf(file="Distribution of alternative splicing type number.pdf",w=10,h=9)
x<-barplot(SixNum,col=heat.colors(6),ylim=c(0,max(SixNum)+skip),xlab="Number of AS Type",ylab="Number of Gene")
text(x,c(SixNum+skip/2),labels=SixNum)
text(x,rep(-skip/2,6),labels=1:6,adj=1,xpd=T)
dev.off()
ASType<-sort(unique(as.character(SpliceSummary2[["Type"]])))
AS54Num<-sapply(ASType,function(x) length(which(SpliceSummary2[["Type"]]==x)))
ASType_GeneNum<-cbind(names(AS54Num),as.numeric(AS54Num))
ratio<-sapply(as.numeric(ASType_GeneNum[,2]),function(x) round(100*x/dim(SpliceSummary2)[1],2))
ASType_GeneNum_Ratio<-cbind(ASType_GeneNum,ratio)
ASType_GeneNum_Ratio_sort<-ASType_GeneNum_Ratio[order(as.numeric(ASType_GeneNum_Ratio[,2]),decreasing=T),]
write.table(ASType_GeneNum_Ratio_sort,"spl5_ox_ASType_GeneNum_Ratio_sort.xls",sep="\t",col.names=T,row.names=F,quote=F)
q()
for(i in 1:6){
file[[i]]<-read.table(dir()[i],header=T,sep="\t")
newfile[[i]]<-file[[i]][1:5]
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,7]))
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,6]))
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,8:9]))
newfile[[i]]<-as.data.frame(newfile[[i]])
colnames(newfile[[i]][6])<-c("p_bonferroni")
colnames(newfile[[i]][7])<-c("p_uncorrected")
}
for(i in 1:6){
file[[i]]<-read.table(dir()[i],header=T,sep="\t")
newfile[[i]]<-file[[i]][1:5]
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,7]))
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,6]))
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,8:9]))
}
file<-list()
newfile<-list()
history()
for(i in 1:6){
file[[i]]<-read.table(dir()[i],header=T,sep="\t")
newfile[[i]]<-file[[i]][1:5]
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,7]))
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,6]))
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,8:9]))
newfile[[i]]<-as.data.frame(newfile[[i]])
colnames(newfile[[i]][6])<-c("p_bonferroni")
colnames(newfile[[i]][7])<-c("p_uncorrected")
}
for(i in 1:6){
file[[i]]<-read.table(dir()[i],header=T,sep="\t")
newfile[[i]]<-file[[i]][1:5]
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,7]))
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,6]))
newfile[[i]]<-cbind(as.data.frame(newfile[[i]]),as.data.frame(file[[i]][,8:9]))
}
dir()
quit()
source("cmd.r")
    count2fpkm_gene<-mergeCountFpkm(count_matrix_path=rsemdir,type=type)
count_matrix_path=rsemdir;type=type
count_matrix_path
 countMatrix<-read.delim(paste(count_matrix_path,"/genes.mat",sep=""),header=T,sep="\t",check.names=F)
 Files1<-colnames(countMatrix)[-1]
 #Sample_name<-gsub(".RSEM.genes.results","",Files1)
 Sample_name<-gsub(".gene","",Files1,perl=T)
Sample_name
 RSEM_File<-paste(count_matrix_path,Files1,sep="/")
 matrix1<-read.delim(RSEM_File[1],check.names=F,header=T,sep="\t")[,c("gene_id","length","expected_count","FPKM")]
         matrix1.sort<-matrix1[order(as.character(matrix1[[1]]),decreasing=F),]
 names(matrix1.sort)[3]<-paste(Sample_name[1],"_count",sep="")
 names(matrix1.sort)[4]<-paste(Sample_name[1],"_fpkm",sep="")
 count2fpkm_matrix<-matrix1.sort
 for(i in 2:length(Sample_name)){
         matrix1<-read.delim(RSEM_File[i],check.names=F,header=T,sep="\t")[,c("gene_id","length","expected_count","FPKM")]
                 matrix1.sort<-matrix1[order(as.character(matrix1[[1]]),decreasing=F),]
 names(matrix1.sort)[3]<-paste(Sample_name[i],"_count",sep="")
 names(matrix1.sort)[4]<-paste(Sample_name[i],"_fpkm",sep="")
 count2fpkm_matrix<-cbind(count2fpkm_matrix,matrix1.sort[,3:4])
 }
i
countMatrix[1:4,]
q()
library(VennDiagram)
install.packages('VennDiagram')
library(VennDiagram)
draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c('red','green'))
venn.diagram(area1=18586,area2=18352,cross.area=17145,fill=c('red','green'))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"),filename="venn_1.pdf")
q()
library(VennDiagram)
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"),cat.col=c("dodgerblue", "goldenrod1"))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"),cat.col=c("dodgerblue", "goldenrod1"))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE)
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=1.5,cat.cex=1.5)
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.05, 3),ext.pos = rep(0, 2), ext.dist = rep(0, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.1, 3),ext.pos = rep(0.5, 2), ext.dist = rep(0.5, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.1, 3),ext.pos = rep(0.05, 2), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.1, 3),ext.pos = rep(1, 2), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.1, 3),ext.pos = c(90, 180), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.01, 3),ext.pos = c(270, 90), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.1, 3),ext.pos = c(270, 90), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.06, 3),ext.pos = c(270, 90), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = rep(0.09, 3),ext.pos = c(270, 90), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0,0.08),ext.pos = c(270, 90), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.05,0.05,0),ext.pos = c(270, 90), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.05,0.1,0),ext.pos = c(270, 90), ext.dist = rep(0.05, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.09,0.09,0),ext.pos = c(270, 90), ext.dist = rep(0.01, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.08,0),ext.pos = c(270, 90), ext.dist = rep(0.01, 2),ext.length = rep(0.5, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.07,0.07,0),ext.pos = c(270, 90), ext.dist = rep(0.01, 2),ext.length = rep(0.9, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.07,0.07,0),ext.pos = c(270, 90), ext.dist = rep(0.02, 2),ext.length = rep(0.9, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.07,0.07,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.9, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.06,0.06,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.9, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.07,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.9, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.9, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.8, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.8, 2),cat.pos=c(0,0))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.8, 2),cat.pos=c(50,-50))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.8, 2),cat.pos=c(-45,45))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.8, 2),cat.pos=c(-45,45),cat.dist = rep(0.025, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.8, 2),cat.pos=c(-45,45),cat.dist = rep(0.03, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0, 2),ext.length = rep(0.8, 2),cat.pos=c(-45,45),cat.dist = rep(0.1, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(0.1, 2),ext.length = rep(0.8, 2),cat.pos=c(-45,45),cat.dist = rep(0.1, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.05, 2),ext.length = rep(0.8, 2),cat.pos=c(-45,45),cat.dist = rep(0.1, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.07, 2),ext.length = rep(0.7, 2),cat.pos=c(-45,45),cat.dist = rep(0.1, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.06, 2),ext.length = rep(0.65, 2),cat.pos=c(-45,45),cat.dist = rep(0.1, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.08, 2),ext.length = rep(0.65, 2),cat.pos=c(-45,45),cat.dist = rep(0.1, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.09, 2),ext.length = rep(0.65, 2),cat.pos=c(-45,45),cat.dist = rep(0.08, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.1, 2),ext.length = rep(0.65, 2),cat.pos=c(-45,45),cat.dist = rep(0.08, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.2, 2),ext.length = rep(0.65, 2),cat.pos=c(-45,45),cat.dist = rep(0.08, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.15, 2),ext.length = rep(0.65, 2),cat.pos=c(-45,45),cat.dist = rep(0.08, 2))
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.13, 2),ext.length = rep(0.5, 2),cat.pos=c(-45,45),cat.dist = rep(0.08, 2))
pdf(file="venn.pdf")
q()
library(VennDiagram)
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.13, 2),ext.length = rep(0.5, 2),cat.pos=c(-45,45),cat.dist = rep(0.08, 2))
pdf(venn.plot,file="venn.pdf")
pdf(venn.plot,file="venn.pdf",w=7,h=8)
pdf(file="venn.pdf",w=7,h=8)
q()
library(VennDiagram)
pdf(file="venn.pdf",w=7,h=8)
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.13, 2),ext.length = rep(0.5, 2),cat.pos=c(-45,45),cat.dist = rep(0.08, 2))
dev.off()
q()
library(VennDiagram)
pdf(file="venn.pdf",w=10,h=10)
venn.plot<-draw.pairwise.venn(area1=18586,area2=18352,cross.area=17145,fill=c("dodgerblue", "goldenrod1"),margin = 0.05,category = c("Sample_W","Sample_K"), ext.line.lty = "solid",ext.text=TRUE,cex=2,cat.cex=2,ext.percent = c(0.08,0.065,0),ext.pos = c(270, 90), ext.dist = rep(-0.13, 2),ext.length = rep(0.5, 2),cat.pos=c(-45,45),cat.dist = rep(0.08, 2))
dev.off()
q()
options(warn=-100)
 library('VennDiagram')
?venn.plot
venn.plot?
draw.pairwise.venn(area1=18586, area2=18352,cross.area=17145,category=c('Sample_W','Sample_K'),lwd=rep(1,1),lty=rep(2,2),,col=c('black','black'),fill=c('dodgerblue','goldenrod1'),scale=TRUE)
pdf(file = paste("venn_1.pdf"),width=100,height=100)
draw.pairwise.venn(area1=18586, area2=18352,cross.area=17145,category=c('Sample_W','Sample_K'),lwd=rep(1,1),lty=rep(2,2),,col=c('black','black'),fill=c('dodgerblue','goldenrod1'),scale=TRUE)
q()
data <- read.table("2.txt", header=TRUE)
data <- read.table("2.txt", header=TRUE,sep="\t",check.names=TRUE)
pdf("tmp.pdf")
boxplot(data, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"))
dev.off()
bye
q
()
q()
pdf("tmp.1.pdf")
data1 <- read.table("2.txt", header=TRUE,sep="\t",check.names=TRUE)
boxplot(Sample~w,data=data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"))
boxplot(data=data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"))
boxplot(data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"))
dev.off()
q()
data1 <- read.table("2.txt", header=TRUE,sep="\t",check.names=TRUE)
pdf("tmp.2.pdf")
boxplot(data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"))
q()
data1 <- read.table("3.txt", header=TRUE,sep="\t",check.names=TRUE)
pdf("tmp.pdf")
boxplot(data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"))
q()
data1 <- read.table("3.txt", header=TRUE,sep="\t",check.names=TRUE)
pdf("tmp.pdf")
boxplot(data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"))
q()
data1 <- read.table("3.txt", header=TRUE,sep="\t",check.names=TRUE)
pdf("tmp.pdf")
boxplot(data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"))
q()
data1 <- read.table("3.txt", header=TRUE,sep="\t",check.names=TRUE)
pdf("tmp.pdf")
boxplot(data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"),xlab="Samples",ylab="dN/dS")
q()
data1 <- read.table("3.txt", header=TRUE,sep="\t",check.names=TRUE)
pdf("boxplot.pdf")
boxplot(data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"),xlab="Samples",ylab="dN/dS",cex=2)
q()
data1 <- read.table("3.txt", header=TRUE,sep="\t",check.names=TRUE)
pdf("boxplot.pdf")
boxplot(data1, col=c("steelblue","mediumturquoise","sandybrown" ,"hotpink","mediumslateblue","wheat"),xlab="Samples",ylab="dN/dS")
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
plot(data1[,2],data[,1],col="red")
pdf("tmp.pdf")
plot(data1[,2],data[,1],col="red")
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
bb<-vector()
for (n in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
}
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
}
bb
plot(aa[,1],aa[,2],col=bb)
plot(data1[,1],data1[,2],col=bb)
pdf("tmp.pdf")
plot(data1[,1],data1[,2],col=bb)
dev.off()
dev.off()
dev.off()
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
pdf("tmp.pdf")
bb<-vector()
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
}
plot(data1[,1],data1[,2],col=bb,xlab="dS",ylab="dN",pch=2)
dev.off()
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
pdf("points.pdf")
bb<-vector()
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
}
plot(data1[,1],data1[,2],col=bb,xlab="dS",ylab="dN",pch=20)
dev.off()
 q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
pdf("points.pdf")
bb<-vector()
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
}
plot(data1[,1],data1[,2],col=bb,xlab="dS",ylab="dN",main="Sample Tk",pch=20)
legend(0.3,0.3,pch=20,legend=c("dN/dS>1","dN/dS<0.1","0.1<=dN/dS<=1"), col=c("red","blue","black"))
dev.off()
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
pdf("points.pdf")
bb<-vector()
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
}
plot(data1[,1],data1[,2],col=bb,xlab="dS",ylab="dN",main="Sample Tk",pch=20)
legend(0.11,0.021,pch=20,legend=c("dN/dS>1","dN/dS<0.1","0.1<=dN/dS<=1"), col=c("red","blue","black"))
dev.off()
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
bb<-vector()
pdf("points.pdf")
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
 }
plot(data1[,1],data1[,2],col=bb,xlab="dS",ylab="dN",main="Sample Tk",pch=20)
legend(0.11,0.026,pch=20,legend=c("dN/dS>1","dN/dS<0.1","0.1<=dN/dS<=1"), col=c("red","blue","black"))
dev.off()
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
bb<-vector()
pdf("points.pdf")
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
}
plot(data1[,1],data1[,2],col=bb,xlab="dS",ylab="dN",main="Sample Tk",pch=20)
legend(0.11,0.029,pch=20,legend=c("dN/dS>1","dN/dS<0.1","0.1<=dN/dS<=1"), col=c("red","blue","black"))
dev.off()
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
bb<-vector()
pdf("points.pdf")
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
}
plot(data1[,1],data1[,2],col=bb,xlab="dS",ylab="dN",main="Sample Tk",pch=20)
legend(0.11,0.0275,pch=20,legend=c("dN/dS>1","dN/dS<0.1","0.1<=dN/dS<=1"), col=c("red","blue","black"))
dev.off()
q()
data1 <- read.table("Tk.dS_dN.txt",sep="\t")
bb<-vector()
pdf("points.pdf")
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
 }
pdf("Tk.points.pdf")
for (i in 1:nrow(data1)){
if(data1[i,2]/data1[i,1]<0.1){bb[i]<-"blue"}
if(data1[i,2]/data1[i,1]>1){bb[i]<-"red"}
if(data1[i,2]/data1[i,1]>=0.1 && data1[i,2]/data1[i,1]<1){bb[i]<-"black"}
}
plot(data1[,1],data1[,2],col=bb,xlab="dS",ylab="dN",main="Sample Tk",pch=20)
legend(0.11,0.028,pch=20,legend=c("dN/dS>1","dN/dS<0.1","0.1<=dN/dS<=1"), col=c("red","blue","black"))
dev.off()
dev.off()
q()
datafile = system.file( "gene_exp_count", package="pasilla" )
pasillaCountTable = read.table( datafile, header=TRUE, row.names=1 )
datafile
pasillaCountTable = read.table( "gene_exp_count", header=TRUE, row.names=1 )
ll
head( pasillaCountTable )
pasillaDesign = data.frame(row.names=colnames(pasillaCountTable),condition=c("group_M","group_M","group_M","group_M","group_T","group_T","group_T","group_T"),libType=c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end"))
pasillaDesign
pairedSamples = pasillaDesign$libType == "paired-end"
countTable = pasillaCountTable[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]
head(countTable)
condition = factor( c( "group_M", "group_M", "group_M", "group_M" ,"group_T", "group_T", "group_T", "group_T") )
library(DESeq)
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds )
res = nbinomTest( cds, "group_M", "group_T" )
head(res)
plotMA(res)
resSig = res[ res$padj < 0.05, ]
resSig1 = res[ resSig$foldChange > 1, ]
ll
head( resSig[ order(resSig1$pval), ] )
head( resSig1[ order(resSig1$pval), ] )
head( resSig1[ order( resSig1$foldChange, -resSig1$baseMean ), ] )
write.csv( resSig1, file="DESseq_result.csv" )
q()
write.table( res, file="DESseq_result.txt", sep="\t" )
q()
source("cmd.r")
q()
library(DESeq)
x <- read.delim("gene_exp_count",row.names="gene_id")
group <- factor(c(1,1,1,1,2,2,2,2))
y <- DGEList(counts=x,group=group)
library(edgeR)
y <- DGEList(counts=x,group=group)
head(y)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)
write.table("et",file="edgeR_result.txt",sep="\t")
q()
write.table(et,file="edgeR_result.txt",sep="\t")
q()
write.table(et,file="edgeR_result.txt",sep="\t")
head(et)
et <- exactTest(y, pair=c("1","2"))
head (et)
topTags(et)
tTags = topTags(et,n=NULL)
write.table(tTags$table, file='edgeR_result.txt', sep='\t', quote=F, row.names=T)
ll
head(edgeR_result.txt)
q()
library("DEGseq")
geneExpMatrix1<=readGeneExp(file=geneExpfile,geneCol=1,valCal=c(2,3,4))
geneExpMatrix1 <- readGeneExp(file="exp_new.txt",geneCol=1,valCol=c(2,3,4))
 geneExpMatrix2 <- readGeneExp(file="exp_new.txt",geneCol=1,valCol=c(5,6,7))
outputDir<-file.path("AvsB")
DEGexp(geneExpMatrix1=geneExpMatrix1,geneCol1=1,expCol1=c(2,3,4),groupLabel1="A",geneExpMatrix2=geneExpMatrix2,geneCol2=1,expCol2=c(2,3,4),groupLabel2="B",method="MARS",outputDir=outputDir)
q()
