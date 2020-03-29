#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %opts;
my $VERSION="1.0";
GetOptions( \%opts,"e=s","o=s","a=s","m=s","h!");

my $usage = <<"USAGE";
       Descript: Used for WGCNA result arrangement
       Version : $VERSION
       Contact : meng.luo\@majorbio.com
       Modify  :2020-03-30
       Usage   : $0 [options]
       Notice  : Run script in WGCNA result dir
       Options :
	        -e*            file            exp.matrix
                -a*            file            annotation file used for run WGCNA
                -o             string          output dir name,default : WGCNA_result
		-m	       num	       num of edges to draw, default : 200
                -h                             Display this usage information
     
       "*" --> must be given Argument      
USAGE
die $usage if ( !( $opts{a} && $opts{e} ) || $opts{h} );
$opts{o}=$opts{o}? $opts{o}:"WGCNA_result";
$opts{m}=$opts{m}? $opts{m}:200;
mkdir($opts{o},"493") or die "Can't create dir at $opts{o}\n" unless -e $opts{o};

#step1 store exp matrix
my %exp;
open EXP,"<$opts{e}" || die "can't open $opts{e}";
my $first_exp=<EXP>;
chomp $first_exp;
my @samples=split /\t/,$first_exp;
while (<EXP>){
	chomp;
	my @line=split /\t/;
	for (my $i=1;$i<=$#line;$i++){
		$exp{$line[0]}{$samples[$i]}=$line[$i];
	}
}
close EXP;
shift @samples;

#step2 store node module infor
my %node_module;
my %module_node;
my %tmp_col;
my $total_gene=0;
open NM,"<all_nodes.txt" || die "cann't open all_nodes.txt, please check your pwd!\n";
<NM>;
map {
	chomp;
	my @line=split /\t/;
	$node_module{$line[0]}=$line[2];
	$module_node{$line[2]}{$line[0]}=1;
	$tmp_col{$line[2]}=1;
	$total_gene+=1;
} <NM>;
close NM;
delete $tmp_col{'grey'};
my @colors=sort (keys %tmp_col);

#step3 out whole stat infor
my %num_node;
open STAT,">$opts{o}/WGCNA_stat.xls";
print STAT "Module\tNode_num\n";
foreach (sort keys %module_node){
	my $nn=scalar (keys %{$module_node{$_}});
	$num_node{$_}=$nn;
	print STAT "$_\t$nn\n";
}
print STAT "Total\t$total_gene\n";
close STAT;
open STAT_R,">$opts{o}/WGCNA_stat.r";
print STAT_R '
library("ggplot2")
a<-read.table("WGCNA_stat.xls",sep="\t",header=T)
a<-a[-which (a$Module == "Total" | a$Module == "grey"),]
add<-max(a$Node_num)/100
a$Module<-factor(a$Module,levels=a$Module)
theme_set(theme_bw())
pdf("WGCNA_stat.pdf",width=16,height=9)
par(mar=c(1,2,2,1))
ggplot(data=a,aes(x=Module,y=Node_num,fill=Module))+scale_fill_manual(values=as.vector(a$Module))+theme(axis.text.x=element_text(size=8,angle=45,hjust=1),plot.title = element_text(hjust=0.5))+ ggtitle("Nodes numeber in Modules")+annotate("text",x=a$Module,y= a$Node_num+add,parse=T,label=a$Node_num,size=3) +geom_bar(stat="identity")
dev.off()
';

system ("cd $opts{o} && Rscript WGCNA_stat.r && cd ..");

#step3 store kMEs
my %kmes;
open KME,"<kME.xls" || die "can't open kME.xls, please check your dir\n";
my $first_kme=<KME>;
chomp $first_kme;
my @raw_modules=split /\t/,$first_kme;
while (<KME>){
	chomp;
	my @line=split /\t/;
	for (my $i=1;$i<=$#line;$i++){
		my $s_module=$raw_modules[$i];
		$s_module=~s/kME//g;
		$kmes{$s_module}{$line[0]}=$line[$i];
	}
}
close KME;

#step4 store kMEs_p
my %kmes_p;
open KMEP,"<kME_pvalue.xls" || die "can't open kME_pvalue.xls, please check your dir";
<KMEP>;
while (<KMEP>){
        chomp;
        my @line=split /\t/;
        for (my $i=1;$i<=$#line;$i++){
		my $s_module=$raw_modules[$i];
		$s_module=~s/kME//g;
                $kmes_p{$s_module}{$line[0]}=$line[$i];
        }
}       
close KMEP;

#step5 store adjacency
my %adjacency;
open ADJ,"<adjacency.xls" || die "can't open adjacency.xls, please check your dir\n";
my $first_adj=<ADJ>;
chomp $first_adj;
my @genes=split /\t/,$first_adj;
while (<ADJ>){
	chomp;
	my @line=split /\t/;
	for (my $i=1;$i<=$#line;$i++){
		$adjacency{$line[0]}{$genes[$i]}=$line[$i];
	}
	delete $adjacency{$line[0]}{$line[0]}; #remove self
}
close ADJ;

#step6 read ano.txt
my %ano;
open ANO,"<$opts{a}" || die "cannot open $opts{a}\n";
while (<ANO>){
	chomp;
	my @line=split /\t/;
	$ano{$line[0]}=$line[1];
}
close ANO;

#step7 output node_stat and hub_stat
my %hub;
open WNST,">$opts{o}/WGCNA_node_stat.xls";
open HUB,">$opts{o}/WGCNA_hubgene_stat.xls";
print WNST "GeneID\tConn-total\tConn-intramodular\tkME\tkME-pvalue\tModule_color\tAnnotation\t".(join "\t",@samples)."\n";
print HUB "GeneID\tConn-total\tConn-intramodular\tkME\tkME-pvalue\tModule_color\tAnnotation\t".(join "\t",@samples)."\n";
#print @colors;exit;
for my $c (@colors){
	my $total_g=0;
	for my $g (sort {$kmes{$c}{$b} <=> $kmes{$c}{$a}} keys %{$kmes{$c}}){
		next if(! exists $module_node{$c}{$g});
		$total_g+=1;
		my $ct=0;
		my $ci=0;
		my $exp_g;
		foreach (@samples){
			$exp_g.="\t".$exp{$g}{$_};
		}
		for my $g2(keys %{$adjacency{$g}}){
			if (exists $module_node{$c}{$g2}){
				$ci+=$adjacency{$g}{$g2};
			}
			$ct+=$adjacency{$g}{$g2};
		}
		print WNST "$g\t$ct\t$ci\t$kmes{$c}{$g}\t$kmes_p{$c}{$g}\t$c\t$ano{$g}"."$exp_g\n";
		# defined head 10% kMEs genes as hub genes
		if ($total_g <= 0.1*$num_node{$c}){
			print HUB "$g\t$ct\t$ci\t$kmes{$c}{$g}\t$kmes_p{$c}{$g}\t$c\t$ano{$g}"."$exp_g\n"
		}
		#select first 3 hub gene for display
		if ($total_g <=3){
			$hub{$c}{$g}=1;
		}
	}
}
close WNST;
close HUB;

#step8 prepare enrichment files
mkdir ("enrichment");
#print scalar (keys %module_node);exit;
for my $c (@colors){
	open DE,">enrichment/module_$c.DE.list";
	foreach (keys %{$module_node{$c}}){
		print DE $_."\n";
	}
	close DE;
}

#step9 auto network
mkdir ("network");
for my $c (@colors){
	my %edge_old;
	open EDGE_OLD, "<CytoscapeInput-edges-$c.txt" || die "cann't find $c edge file";
	<EDGE_OLD>;
	while (<EDGE_OLD>){
		my @line=split /\t/;
		$edge_old{$line[0]."\t".$line[1]}=$line[2];
	}
	close EDGE_OLD;
	my %num_hub;
	my $count=0;
	open EDGE,">network/$c.edge.txt";
	open NODE,">network/$c.node.txt";
	print EDGE "From\tTo\tWeight\n";
	print NODE "Node\tDegree\n";
	for my $pairs (sort {$edge_old{$b} <=> $edge_old{$a}} keys %edge_old){
		my @line=split /\t/,$pairs;
		if (exists $hub{$c}{$line[0]}  or exists $hub{$c}{$line[1]}){
			$count++;
			if ($count <= $opts{m}){
				print EDGE "$pairs\t$edge_old{$pairs}\n";
				if (exists $num_hub{$line[0]}){
					$num_hub{$line[0]}+=1;
				}else{
					$num_hub{$line[0]}=1;
				}
				if (exists $num_hub{$line[1]}){
					$num_hub{$line[1]}+=1;
				}else{
					$num_hub{$line[1]}=1;
				}
			}
		}
	}
	close EDGE;

	foreach (keys %num_hub){
		print NODE "$_\t$num_hub{$_}\n";
	}
	close NODE;
}

chdir ("network");
open NET,">net.r";
print NET "library(\"igraph\")\n";
for my $c (@colors){
	print NET "
### module $c ###
links <- read.table(\"$c.edge.txt\",header=T,sep=\"\\t\")
net <- graph_from_data_frame(d=links, directed=F)
net <- simplify(net, remove.multiple = T, remove.loops = T)

#degree
k <- degree(net, mode=\"all\")
count <- as.data.frame(k)
colnames(count) <- \"degree\"

#gene size
gene_sep <- quantile(count\$degree, probs = seq(0, 1, by = 0.2))
gene_size <- count\$degree
base_size <- 3
gene_size[which(count\$degree < gene_sep[2])] <- base_size
gene_size[which(count\$degree < gene_sep[3] & count\$degree >= gene_sep[2])] <- base_size + 1
gene_size[which(count\$degree < gene_sep[4] & count\$degree >= gene_sep[3])] <- base_size + 2
gene_size[which(count\$degree < gene_sep[5] & count\$degree >= gene_sep[4])] <- base_size + 3
gene_size[which(count\$degree < gene_sep[6] & count\$degree >= gene_sep[5])] <- base_size + 4
gene_size[which(count\$degree >= gene_sep[6])] <- base_size + 10
V(net)\$size <- gene_size

#label size
label_size <- count\$degree
base_size <- 0.2
label_size[which(count\$degree < gene_sep[6])] <- base_size
label_size[which(count\$degree >= gene_sep[6])] <- base_size * 3
V(net)\$label.cex <- label_size

#frame color
fram_col <- count\$degree
fram_col[which(count\$degree < gene_sep[6])] <- \"$c\"
fram_col[which(count\$degree >= gene_sep[6])] <- \"darkgreen\"
V(net)\$frame.color=fram_col

#edge width
edge_wide <- links\$Weight
edge_sep <- quantile(edge_wide, probs = seq(min(edge_wide), max(edge_wide), by = ((max(edge_wide)-min(edge_wide))/5)))
edge_wide[which(links\$weight < edge_sep[2])] <- 0.05
edge_wide[which(links\$weight < edge_sep[3] & links\$weight >= edge_sep[2])] <- 0.15
edge_wide[which(links\$weight < edge_sep[4] & links\$weight >= edge_sep[3])] <- 0.25
edge_wide[which(links\$weight < edge_sep[5] & links\$weight >= edge_sep[4])] <- 0.35
edge_wide[which(links\$weight < edge_sep[6] & links\$weight >= edge_sep[5])] <- 0.45
edge_wide[which(links\$weight >= edge_sep[6])] <- 0.55
E(net)\$width <- edge_wide

pdf(\"net_$c.pdf\", h=12, w=12)
plot(net,
		edge.color=\"grey80\", 
		vertex.color=\"$c\", 
		vertex.label.color=\"black\",
		edge.curved=0,
		main=\"Module $c network\",
		vertex.label.font=1,
		layout=layout_in_circle
)
###################################################################
";
}
close NET;
system ("Rscript net.r && cd ..");
