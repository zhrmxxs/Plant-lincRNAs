for plant in rice arabidopsis maize tomato
do
	cd /public3/labmember/xujw/lincRNA/outs/08.cis-acting/${plant}
	cat /public3/labmember/xujw/scRNA/genome/${plant}/information/mRNA.gtf|awk '$3=="gene"{print}' |awk 'BEGIN{FS="[\t;]";OFS="\t"} {for(i=1; i<=NF; i++){ if( $i ~ "gene_id"){print $1,$4,$5,$7,$i}}}' |awk -F'["\t]' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6}' > mRNA 
	cat /public3/labmember/xujw/scRNA/genome/${plant}/information/lincRNA.gtf|awk '$3=="gene"{print}' |awk 'BEGIN{FS="[\t;]";OFS="\t"} {for(i=1; i<=NF; i++){ if( $i ~ "gene_id"){print $1,$4,$5,$7,$i}}}' |awk -F'["\t]' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6}' > lincRNA 

	cat mRNA|while read mRNA;do (arr=($mRNA); cat lincRNA|awk -v CHR=${arr[0]} -v START=${arr[1]} -v End=${arr[2]} -v STRAND=${arr[3]} -v ID=${arr[4]} 'BEGIN{FS="\t";OFS="\t"}{if($1==CHR && (($2>End && $2-End<=5000)||($3<START && START-$3<=5000 ))){print $5,ID}}' );done > cis-acting_5000_no_strand

done

conda activate hdWGCNA
Rscript -e '
library(Seurat)
library(tidyverse)
library(cowplot)
library(igraph)
library(WGCNA)
library(hdWGCNA)
library(pheatmap)
library(openxlsx)
###########################################
list <- read.table("/public3/labmember/xujw/lincRNA/config")
list <- as.character(list$V1)
for(id in list)
# id=list[10]
{
plant <- strsplit(id,"_")[[1]][1]
setwd("/public3/labmember/xujw/lincRNA/outs/08.cis-acting/")
if(!dir.exists(id)){dir.create(id)}
setwd(id);  output_path <- "./"

cis <- read.table(paste0("/public3/labmember/xujw/lincRNA/outs/08.cis-acting/",
            plant,"/cis-acting_5000_no_strand"))
colnames(cis) <- c("lincRNA","mRNA")

seurat_obj <- readRDS(paste0("/public3/labmember/xujw/lincRNA/outs/03.hdWGCNA/",id,"/hdWGCNA.rds")) # hdWGCNA rds
dataExpr <- GetDatExpr(seurat_obj) # Metacell
###### Calculate correlations using Metacell
cis <- cis[cis$mRNA %in% colnames(dataExpr) & cis$lincRNA %in% colnames(dataExpr),]
COR <- c()
for(i in 1:nrow(cis)){
    cor <- cor(dataExpr[,cis[i,1]],dataExpr[,cis[i,2]],method="pearson")
    COR <- c(COR, cor)
}
cis$cor_Metacell <- COR

findallmarkers <- read.table(paste0("/public3/labmember/xujw/lincRNA/outs/02.markers/",id,"_findallmarkers.xlsx"),header=T,sep="\t")
findallmarkers <- findallmarkers[order(findallmarkers$gene,findallmarkers$avg_log2FC,decreasing=T),]
findallmarkers <- findallmarkers %>% group_by(gene) %>% 
                summarize(celltype = paste0(cluster," (",avg_log2FC,")",collapse=", "))
colnames(findallmarkers) <- c("lincRNA","celltype_lincRNA")
cis <- merge(cis,findallmarkers,by="lincRNA",all.x=T)
colnames(findallmarkers) <- c("mRNA","celltype_mRNA")
cis <- merge(cis,findallmarkers,by="mRNA",all.x=T)
write.xlsx(cis,"./cis_5k.xlsx")
}
'