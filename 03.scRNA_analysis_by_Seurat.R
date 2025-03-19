library(devtools);
library(Seurat);
library(cowplot);
library(harmony);
library(patchwork);
library(ggplot2);
library(ppEffect);
library(rmarkdown);
library(markdown);
library(dplyr);
library(plotly);
library(ggvenn);
library(RColorBrewer)
library(clustree)
library(reshape2)
###############################################################
plant = "rice"
project = "PRJNA609092"
sample = "SRR11194114"
setwd(dir = paste0("/public3/labmember/xujw/scRNA/annotation/",plant,"/",project,"/",sample))
###############################################################
counts <- Read10X(data.dir= paste0("/public3/labmember/xujw/scRNA/matrix/",plant,"/",project,"/", sample, "_all_transcript/outs/filtered_feature_bc_matrix", gene.column= 1))
seurat_obj <- CreateSeuratObject(counts= counts, project= project, min.cells=3, min.features=200)
###############################################################
# seurat_obj <- NormalizeData(seurat_obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(seurat_obj , verbose = FALSE)
seurat_obj <- SCTransform(seurat_obj) %>% RunPCA(seurat_obj ,npcs = 50, verbose = FALSE)
p <- ElbowPlot(seurat_obj[[time]], ndims = 50, reduction = "pca") 
        ggsave(p, file = "ElboewPlot.png")
seurat_obj <- FindNeighbors(dims = 1:30, verbose = FALSE) %>% 
        FindClusters(resolution = 0.5, verbose = FALSE) %>% 
        RunUMAP(dims = 1:30) %>% 
        RunTSNE(dims = 1:30)
######## Multiple sample integration
seurat_obj <- merge(obj1, y = c(obj2,obj3,obj4))
seurat_obj <-  SCTransform(seurat_obj) %>%
        RunPCA(npcs = 50, verbose = F) %>%
        RunHarmony(group.by.vars = "orig.ident") %>%
        FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
        FindClusters(resolution = 0.5) %>%
        RunUMAP(reduction = "harmony", dims = 1:30) %>%
        RunTSNE(reduction = "harmony", dims = 1:30)

###############################################################  
######## clustree
library(clustree)
pdf(file = "./outs/04_clustree.pdf", width = 9, height = 9)
obj <- seurat_obj
seq <-seq(0.1,1, by=0.1)
for(res in seq){obj <-FindClusters(obj, resolution=res )}
p <- clustree(obj@meta.data, prefix= 'SCT_snn_res.')
print(p)
dev.off()
##### Change resolution
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE) 
###############################################################
######## annotation, according to the marker gene or the author's barcodes
seurat_obj$celltype <- "x"
######### save
saveRDS(seurat_obj, file = "annotation.RDS")
