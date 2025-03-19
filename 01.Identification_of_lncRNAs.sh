########### using NGS data to identify lncRNAs
########### Preparation: Reference genome sequence, reference genome annotation, NGS raw data
genome=Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
gtf=Oryza_sativa.IRGSP-1.0.56.gtf
fq1=fastq1.gz 
fq2=fastq2.gz 
sample=example

##################################################################
########### Quality control of raw data
########### fastqc and trim-galore
trim_galore -q 25 --phred33 --fastqc --length 35 --stringency 3 --paired -o ./  $fq1 $fq2

##################################################################
########### Mapping to the genome
########### hisat2 
extract_exons.py $gtf > genome.exon
extract_splice_sites.py $gtf > genome.ss
hisat2-build $genome --ss genome.ss --exon genome.exon index
hisat2 -p 10 -x index -1 $fq1 -2 $fq2|samtools sort -@ 10 -o ./${sample}.bam
hisat2 -p 10 -x index -U $fq1|samtools sort -@ 10 -o ./${sample}.bam # If the raw data is Single-Read

##################################################################
########### stringtie the transcripts
stringtie -B -p 10 -G $gtf -c 2 -o ./${sample}_stringtie.gtf -A ./${sample}.tab -l $sample ${sample}.bam
########### When there is more than one sample
stringtie --merge -p 8 -G $gtf -o ./merged.gtf list_of_sample

##################################################################
########### Classification of transcripts
########### gffcompare  
gffcompare candidate.gtf -R -r $gtf -o gffcompare
awk '($3 ~ /^[xjoiup]/) && ($10 >= 200) {print $5}' gffcompare.merge.gtf.tmap > candidate.id 
### x means NAT lncRNAs, j and o mean sense lncRNAs, i means intronic lncRNAs, u and p mean lincRNAs.
cat merged.gtf|grep -w -Ff candidate.id > candidate.gtf
gffread  candidate.gtf -w candidate.fa -g $genome

########### The 3 files (candidate.id, candidate.gtf and candidate.fa) are informations of candidate lncRNAs.

##################################################################
########### Identify open reading frames (ORF)
########### TransDecoder
TransDecoder.LongOrfs -t candidate.fa --output_dir transdecoder_out
cat transdecoder_out/longest_orfs.gff3|grep -v '^\s*$'|cut -f1|sort -u > orf_101.ID
cat candidate.id |grep -v -w -Ff orf_101.ID > filtered_by_ORF.ID

##################################################################
########### Filter out the housekeeping RNA
########### Rfam database
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.4/Rfam.cm.gz
cmpress Rfam.cm # indexing
cmscan -Z 216349 --cut_ga --rfam --nohmmonly  --tblout my-genome.tblout --fmt 2 \
        --clanin Rfam.clanin Rfam.cm candidate.fa > my-genome.cmscan
cat my-genome.tblout |awk '{if($20=="=")print $4}'|sort -u > ncRNA.ID
cat candidate.id |grep -v -w -Ff ncRNA.ID > filtered_by_Rfam.ID

##################################################################
########### Predict coding ability
########### CPC2
python /public/home/xujw/app/CPC2-master/bin/CPC2.py -i $fa -o cpc2_result
cat cpc2_result |awk '{if($8=="noncoding") print $1}' > filtered_by_CPC2.ID

##################################################################
########### Filter protein functional threshold
########### Pfam database
transeq candidate.fa filter_protein.fa -frame=6
pfam_scan.pl -fasta ./filter_protein.fa -dir /public3/labmember/xujw/biosoft/database/Pfam -out ./Pfam_scan.out
grep -v '^#' Pfam_scan.out | grep -v '^\s*$' | awk -F "_" '{print$1}' | sort | uniq > protein.ID
cat candidate.id |grep -v -w -Ff protein.ID > filtered_by_Pfam.ID

##################################################################
########### Final identification result
cat filtered_by_ORF.ID filtered_by_Rfam.ID filtered_noncoding.ID filtered_by_Pfam.ID|sort -n|uniq -c|awk '{if($1=="4")print $2}' > final.id
cat candidate.gtf |grep -w -Ff final.id > final.gtf
gffread -w final.fasta -g $genome ./final.gtf
seqtk seq -l 0 final.fasta > final.fa

########### Merge lncRNA and mRNA into the same GTF file
cat mRNA_gtf lncRNA_gtf > cat.gtf;
cat cat.gtf | awk 'BEGIN{FS="[\t;]";OFS="\t"} {for(i=1; i<=NF; i++){ if( $i ~ "gene_id"){print $1,$4,$i}}}' |awk -F'["\t]' 'BEGIN{OFS="\t"}{print $1,$2,$4}' > chr_start_gene_id
# In R language, order all gene ids according to chromosomes and loci
Rscript -e '
library(dplyr)
df <- read.table("./chr_start_gene_id",sep= "\t")
df <- df %>% group_by(V3) %>% filter(V2==min(V2))
df <- df[!duplicated(df$V3),]
df <- df[order(df$V1,df$V2),]
write.table(df,file = "chr_start_gene_id.sort_uniq",col.names=F,row.names= F,sep="\t")
'
cat $mRNA_gtf |grep "#" >  mRNA_lncRNA.gtf
cat chr_start_gene_id.sort_uniq|cut -f3|while read id;
do
cat cat.gtf|grep -w -F $id >> mRNA_lncRNA.gtf
done