

qs|grep step02.uniform.s |cut -d "." -f1|tr '\n' ' ' kill all the procees READ THE RELATE FILE IN READ

nohup perl ~dna/Environment/qsub-sge.pl --Resource mem=500G --CPU 32 

ls *.sh|perl -ne 'chomp;print "$_ "if(!-f"$_.check")' ###read all the check.sh file 
ex..cat step03.20180807.ustack.sh step03.20180807.ustack.sh.reqsub.sh step03.ustacks.sh_00024.sh
 step03.ustacks.sh_00025.sh step03.ustacks.sh_00028.sh step03.ustacks.sh_00033.sh step03.ustacks.sh_00049.sh
 > step03.20180808.ustack.sh

le total.gt|cut -f 1|uniq -c   read uniq chromosome
le total.gt|cut -f 1|uniq |wc -l
tar -zxvf filename.tar.gz

tar -zcvf filename.tar.gz newfolder here


BAIDU Cloud 
DNA@majorbaio.com
majorDNA1120
 #tree


lfs quota /mnt/ilustre/ -h


Rscript /mnt/ilustre/users/long.huang/Pipeline/GIT/bioinformatics/Pipeline/07.POP/bin/bin/tree.R --infile pop.nj.tree --group g.list --outfile ./
 
 
##before BSA
## make the fq.list
 
ls *1.fastq.gz|perl -e 'print "Sample Initial Name\tSample Analysis Name\tBatch\tLibrary\tFile Name\n";while(<>){chomp;@name=split(/\./,$_);$ne=$name[0];$ny=$_;$ny=~s/1.fastq.gz/2.fastq.gz/g; print "$ne\t$ne\tbatch1\tWGS\t$_,$ny\n"}' >fr.list

RF *1.fastq.gz|perl -ne 'chomp;@name=split(/\//,$_);@ne=split(/\./,$name[-1]);$nea=$ne[0];$ny=$_;$ny=~s/1.fastq.gz/2.fastq.gz/g; print "$nea\t$_\t$ny\tWGS\n"' >fqa.list 
le  fqa.list |cut -f1,2|perl -ne 'chomp;($a,$b)=split/\s+/,$_;$c=`less $b|head -n1`;print "$a\t$c\n"' > stat.xls

###sangerupload 
 
SangerUpload -i /mnt/ilustre/users/meng.luo/project/maoguangzhi_BSA/maoguangzhi_bsa/ -c 'VXDKZZ|0bd8a0aa4e81018a015f7cd1aa515211' -m sanger -l maoguangzhi
 
##download data

 python /mnt/ilustre/users/dna/Environment/cloud/dataexchange-V2/download.py -c '' -t ./ -m sanger
 ##print gene in mutmap on linux use the perl 
le region.threshold.gene.total|perl -ne 'chomp;($n1,$n2)=split(/\s+/,$_,2);if ($n1 =~ /@/){print "$_\n"}' 



####write the ID sample name > the group.list
head -n1000 step01.vcf-filter/pop.recode.vcf|grep "#C"|cut -f 10-100|sed 's/\t/\n/g'|perl -ne 'chomp;print "$_\t1\n"' >group.list

##liangzhiqiang
##remove the de.list

grep -v -f de.list total.list >11.list

##take out the de.list 

grep  -f de.list total.list >22.list

cat total.list 22.list >all.list 

cut -f1 all.list|sort -u|wc -l
le all.list |sort|wc -l
cut -f1 all.list|sort -u


###liangzhiqiang no reference genome

cut -f1 all.list|sort|uniq|perl -ne 'chomp;$m=rand()*1;print "$_\t$m\n"'|sort -k2,2|head -n50|cut -f1 >sample.list

cut -f1 all.list|sort|uniq|perl -ne 'chomp;$m++;print "$_\t$m\n"' >group.list 

less step02.uniform.sh|perl -ne 'chomp;$na=(split(/\s+/))[-1];$nb=(split(/\./,$na))[2]; print "$nb.check\n"'

less step02.uniform.sh|perl -ne 'chomp;$na=(split(/\s+/))[-1];$nb=(split(/\./,$na))[2]; print "touch /mnt/ilustre/users/meng.luo/project/liangzhiqiang_NOREF_20180831/data/02.$nb.check\n"'

less step02.uniform.sh|perl -ne 'chomp;$na=(split(/\s+/))[-1];$nb=(split(/\./,$na))[2]; print "touch /mnt/ilustre/users/meng.luo/project/liangzhiqiang_NOREF_20180831/data/02.$nb.check\n"' >02.sh



###collinearrity

map
##map not right with empty

Fragaria_vesca_v4.0.a1.fasta|perl -ne 'chomp;@a=split; next if(/^$/ ||$a[0] eq ""); print "$_\n"' >Fragaria_vesca.fa

samtools faidx Fragaria_vesca.fa

head -n7 Fragaria_vesca.fa.fai|cut -f1,2 >chrlength.map

le mappp.mp|cut -f1,4,5|grep Fvb|perl -ne '@a=split;@b=split(/\:/,$a[0]);print "$b[0]\t$a[1]\t$a[2]\n"' >mappp.map

perl /mnt/ilustre/users/long.huang/bin/Bin/old/genetic_map/v1.2.1/graph/drawAligmentRalationMap.pl -m ../mapp1.map -a mapp1.map -k map1 -o ./

cairosvg -f png map1.svg >map1.png




###comsumetimes

RAxML_info.pop
le RAxML_info.pop|grep ": Time"|cut -d " " -f3|awk '{s+=$1}END{print s}'|le

###delete duplicatble

le mapp2.map|cut -f1|uniq -d >de2.list

grep -v -f dep.list mm2.map >de_mm2.map

perl /mnt/ilustre/users/long.huang/bin/Bin/old/genetic_map/v1.2.1/graph/drawAligmentRalationMap.pl -m ../mapp1.map -a mapp1.map -k map1 -o ./

cairosvg -f png mm1.svg >mm1.png



##check all users/space

lfs quota -u meng.luo /mnt/ilustre/

perl snp.indel.number.pl -infile index-calc.result.index -outfile 3-10.csv
grep '@chr' region.threshold.vcf.total|sed 's/\t/,/g' >3-11.csv
grep '@chr' region.threshold.gene.total|sed 's/\t/,/g'|sed 's/#@chr/Chromosome ID/g'|sed 's/@//g'|sed 's/pos1/Start/g'|sed 's/pos2/end/g' >3-12.csv

grep '@chr' region.threshold.gene.total|sed 's/\t/,/g' >



###all work use the tools for make more work effectively


perl -e '@a=glob("*/05.annovar/*/work_flow/");print join("\ncd ",@a),"\n"' > work.sh

less -S work.sh |sed 's/snp\.indel\.number\.pl/"/mnt/ilustre/users/meng.luo/project/liujianxiang_20180903/01MJ20171213041_bsa_mutmap/05.annovar/C10XC/work_flow/snp.indel.number.pl"/g'|less -S 

readlink -f ../*/*/*/work_flow|perl -ne '@a=split(/\//,$_);print $a[-2],"\n"'

*/05.annovar/*/work_flow/


##variation step01 to make the GCTA figure

for i in `ls` ;do mv -f $i `echo $i|sed 's/:/-/'`;done

sed -i 's/:/-/g' *.stat
le 01.fastq-qc.sh|perl -ne 'chomp;@st=split/\&&/;print "$st[-1]\n"' >01fig.sh 

lfs quota -u meng.luo /mnt/ilustre/
Disk quotas for user meng.luo (uid 804):
     Filesystem  kbytes   quota   limit   grace   files    quota   limit   grace
  /mnt/ilustre/ 14989243096       0       0       -  134076       0       0       -

亲本平均测序深度为16.11X(总的数据量（多少G/基因组的大小）)，子代平均测序深度为10.33X。样品与参考基因组平均比对效率为93.37 %，
平均覆盖深度为16.10X(最后一列)，基因组覆盖度为%（至少一个碱基覆盖（第三列））（3-5.table）。

运行过的命令要记得保存一份！

##change the Total.lg file name #Gmap

le pop.filtered.marker|cut -f1|grep -v "#" |perl -ne 'chomp;my @a=split(/\_/,$_);if(!exists $hash{$a[0]}){print "\n>$a[0]\n$_";$hash{$a[0]}=1}else{print " $_ ";}'|less -S >LG12

##check error

ls *.sh|perl -ne 'chomp;print "$_ "if(!-f"$_.check")'


新的上传数据

python /mnt/ilustre/users/meng.luo/Pipeline/V3/get_file_list.py -i /mnt/ilustre/users/meng.luo/project/meiyongxian_MJ20180807035/data -l outfile

python /mnt/ilustre/users/meng.luo/Pipeline/V3/upload.py -i /mnt/ilustre/users/meng.luo/project/meiyongxian_MJ20180807035/data -l outfile -c "NTMPAW|618d3db615ca3d000cbbdc3bbb6f0ede" -m "nsanger"

##new

python /mnt/ilustre/users/meng.luo/Pipeline/V3/get_file_list.py -i ./ -l outfile

python /mnt/ilustre/users/meng.luo/Pipeline/V3/upload.py -i ./ -l outfile -c "DGQTDR|7532229f97ee86116f47331df107c43a" -m "nsanger"

##download data
python ~chunxiang.xue/scripts/i-sanger/V3/download.py -c 'OGHFOJ|bc4638361a7614f0fc5c01fcb2c1f1d2' -t ./ -m nsanger


python /mnt/ilustre/users/dna/Environment/cloud/dataexchange-V2/download.py -c '' -t ./ -m sanger

python /mnt/ilustre/users/qingmei.cui/dataexchange/V3/download.py -c '' -t ./ -m nsanger


python ~chunxiang.xue/scripts/i-sanger/V3/download.py -c '' -t ./ -m nsanger 



qstat|grep rna. |cut -d "." -f1|tr '\n' ' '

qdel kill the processs

# list the mechine
lfs df

qsource 
qcheck
qstatf
qhost
##PCA change the name for the GCTA
le pop.recode.vcf |perl -ne 'chomp;if(/^#/ || /^##/){print "$_\n";}else{my ($chr,$all)=split/\s+/,$_,2;print "Un\t$all\n"}'>pop.un.vcf

##提取表中的结果!两个表中是一样的结果!
less -S  ~long.huang/newmdt/WTH/2018-11-27/result/gatk/fst-thetapi-tajima/pop2_pop3.fst.thetapi.tajima.pop1.select.snp ~long.huang/newmdt/WTH/2018-11-27/vcf/gatk.snp.table|perl -ne '@a=split;$b=$a[0]."-".$a[1];$stat{$b}++;print $_ if($stat{$b} > 1);' >pop2_pop3.fst.thetapi.tajima.pop1.select.snp.table && less -S  ~long.huang/newmdt/WTH/2018-11-27/result/gatk/fst-thetapi-tajima/pop2_pop3.fst.thetapi.tajima.pop2.select.snp ~long.huang/newmdt/WTH/2018-11-27/vcf/gatk.snp.table|perl -ne '@a=split;$b=$a[0]."-".$a[1];$stat{$b}++;print $_ if($stat{$b} > 1);' >pop2_pop3.fst.thetapi.tajima.pop2.select.snp.table

##修改提取的区域
le sgRNA.region |grep "chr7"|cut -f 1,2,3|perl -ne 'chomp;my ($a,$b,$c)=split/\s+/,$_;my $start=$b-5000;my $end=$c+5000;print "$a\t$start\t$end\n";'>5K.S2.chr7.pos

#get fa.list for qc

RF */*1_001.fastq.gz|perl -ne 'chomp;@name=split(/\//,$_);@ne=split(/\./,$name[-1]);$nea=$ne[0];$ny=$_;$ny=~s/1_001.fastq.gz/2_001.fastq.gz/g; print "$nea\t$_\t$ny\n"' >fqa.list

# bwa
bwa index -a bwtsw M.clean.1.fa
bwa mem -t 4 /mnt/ilustre/users/minghao.zhang/newmdt/Project/MJ20180612007_wanhuihua_Gmap/var_20180619/02.ref-config/ref query.fq > query.sam
bwa mem  -M -a -t 8 -R "@RG\tID:2\tLG:M\tLB:1\tPL:illumina\tSM:M\tPU:run_barcode\tCN:MajorBioDS:reseq" \
/mnt/ilustre/centos7users/meng.luo/project/yuanjie_MJ20190619038/bwa/1/M.clean.1.fa \
/mnt/ilustre/centos7users/meng.luo/project/yuanjie_MJ20190619038/04.bam-sort/Z.sort.mapped.bam.R1.fastq.gz\
/mnt/ilustre/centos7users/meng.luo/project/yuanjie_MJ20190619038/04.bam-sort/Z.sort.mapped.bam.R2.fastq.gz \
| samtools view -bS - > /mnt/ilustre/centos7users/meng.luo/project/yuanjie_MJ20190619038/bwa/1/Z.b1.bam


#blast
blat /mnt/ilustre/users/minghao.zhang/newmdt/Genome/newGenome/Rosa_chinensis/01.newref/ref.fa all.fa -t=dna -q=dna -out=blast ssr3
makeblastdb -in /mnt/ilustre/users/minghao.zhang/newmdt/Genome/newGenome/Rosa_chinensis/01.newref/ref.fa \
-dbtype nucl -title ssr1 -parse_seqids -out ssr1 -logfile ssr1.log && blastn -query /mnt/ilustre/users/meng.luo/project/wanhuihhua/SSR/query1.fa \
-db /mnt/ilustre/users/meng.luo/project/wanhuihhua/SSR/ssr1 -evalue 1e-5 -num_threads 8  -outfmt 7 -out ssr1 


#bowtie
bowtie-build M.clean.1.fa M.clean.1.fa
bowtie -f -n 1 -e 80 -l 18 -a -m 5 --best --strata --al M.clean.1.mapped --un M.clean.1.unmapped M.clean.1.fa \
/mnt/ilustre/centos7users/meng.luo/project/yuanjie_MJ20190619038/04.bam-sort/Z.sort.mapped.bam.fa M.clean.1.bwt 2> M.clean.1.log



#install 

python setup.py install --prefix=/mnt/ilustre/centos7users/meng.luo/.env/bin


## remove empty with start

le populations.new.tag |grep -v '^\s*$' > population.tag

