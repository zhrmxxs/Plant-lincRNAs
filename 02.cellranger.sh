##############################################################
# make reference for cellranger
cd /public3/labmember/xujw/scRNA/genome/rice
        cellranger mkref --genome=refdata_all_transcript \
        --fasta=Oryza_sativa.IRGSP-1.0.dna.toplevel.fa \
	--genes=Oryza_sativa_mRNR_lincRNA.gtf;
cd /public3/labmember/xujw/scRNA/genome/arabidopsis
        cellranger mkref --genome=refdata_all_transcript \
        --fasta=Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
        --genes=Arabidopsis_thaliana_mRNA_lincRNA.gtf;
cd /public3/labmember/xujw/scRNA/genome/tomato
        cellranger mkref --genome=refdata_all_transcript \
        --fasta=S_lycopersicum_chromosomes.4.00.fa \
        --genes=Solanum_lycopersicum_mRNA_lincRNA.gtf ;
cd /public3/labmember/xujw/scRNA/genome/maize
       cellranger mkref --genome=refdata_all_transcript \
       --fasta=Zea_mays.B73_RefGen_v4.dna.toplevel.fa \
       --genes=Zea_mays_mRNA_lincRNA.gtf;

##############################################################
# generate the expression matrix containing both lncRNAs and mRNAs
plant=rice
cellranger count --id ${b}_all_transcript \
	--fastqs /public3/labmember/xujw/lncRNA/RICE/PRJNA706435/lncRNA_project/01.raw_data/ \
	--sample ${b} \
	--transcriptome /public3/labmember/xujw/scRNA/genome/${plant}/refdata_all_transcript/
done