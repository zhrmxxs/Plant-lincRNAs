## conda activate hdWGCNA
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(devtools)
library(igraph)
library(corrplot)
library(WGCNA)
library(hdWGCNA)
library(pheatmap)
library(openxlsx)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)
id = "example"  # id = "rice_root"
setwd(paste0("/public3/labmember/xujw/lincRNA/outs/03.hdWGCNA/",id))
output_path <- "./"

seurat_obj <- readRDS(paste0("/public3/labmember/xujw/lincRNA/outs/01.rds/",id ,".rds"))

assay="SCT"
Idents(seurat_obj) <- seurat_obj$celltype

######################################################################
### step 1
#创建对象，选择基因， 根据细胞分组选择基因
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  group.by = "celltype",
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
length(seurat_obj@misc$tutorial$wgcna_genes)

######################################################################
### step 2
#创建metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
 # assay = assay,
  group.by = "celltype", # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  reduction = "harmony", # 
  slot='counts',
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype' # set the Idents of the metacell seurat object
)

seurat_obj <- NormalizeMetacells(seurat_obj)
metacell_obj <- GetMetacellObject(seurat_obj)
clusters<- names(table(seurat_obj@misc$tutorial$wgcna_metacell_obj$celltype))

######################################################################
#### step 3
##### select celltypes, here choose all
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = clusters ,
  group.by='celltype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
 # assay = assay, # using SCT assay
  slot = 'data' # using normalized data
)

# t<-GetDatExpr(seurat_obj)
#  t<-cbind(rownames(t),t)
#  write.table(t,file=paste0(output_path,"01.DatExpr.xlsx"),sep="\t",quote=F,row.names=F,col.names=T)

# metacell clusters
runMetaUMAP <- function(obj,title= "Metacell"){
my20color <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
               "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
               "#008b8b", "#e6ab02", "#a6761d", "#b15928", "#8dd3c7",
               "#d73027", "#4575b4", "#313695", "#fec44f", "#800080")     
metacell_obj <- GetMetacellObject(obj)
metacell_obj <- SCTransform(metacell_obj,assay="SCT")
metacell_obj <- RunPCA(metacell_obj, npcs = 30, verbose = FALSE)
metacell_obj <- RunUMAP(metacell_obj, reduction = "pca", dims = 1:30)
p <- DimPlot(metacell_obj, reduction = "umap", group.by = "celltype", 
    label = TRUE, repel = TRUE, cols = my20color) + 
    theme_void() +
    theme(
      legend.key.size = unit(2, "lines"),
      text = element_text(size = 18),  # 调整主要文本的字体大小
#      axis.text = element_text(size = 12),  # 调整坐标轴文本的字体大小
#      axis.title = element_text(size = 14),  # 调整坐标轴标题的字体大小
      legend.text = element_text(size = 18),  # 调整图例文本的字体大小
      plot.title = element_text(hjust = 0.5, size = 20)  # 调整标题的水平位置和字体大小
    ) +
    ggtitle(title)
#     theme(legend.position = "none") # 不要图例
return(p)  
}

p <- runMetaUMAP(seurat_obj)
ggsave(p,filename= "./01_metacell_UMAP.png",width = 12, height =9 )

######################################################################
###### STEP4
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed', 
#  networkType = 'unsigned',
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2))) # 
  # networkType = 'signed' # "unsigned" or "signed hybrid"
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj,
                            point_size = 5,
                            text_size = 3)
                            # selected_power = NULL
# assemble with patchwork
plot3 <-wrap_plots(plot_list, ncol=2) 
ggsave(plot3,filename=paste0(output_path, "02_Wrap_signed.pdf"),width=12,height=9)
power_table <- GetPowerTable(seurat_obj)

##############
softPower<-power_table[power_table$SFT.R.sq>=0.8,1][1]

######################################################################
# step5
     #TOM
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power=softPower,
  setDatExpr=FALSE,
  corType = "pearson",
  networkType = "signed", #区分正相关、负相关， "unsigned"不区分
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 50,  #模块的最少基因数
  mergeCutHeight = 0.2, #合并模块的阈值
  tom_outdir = "TOM", # 输出文件夹
  tom_name = "hdWGCNA" # name of the topoligical overlap matrix written to disk
)

 pdf(file=paste0(output_path, "03_Dendrogram.pdf"),width=1200,height=900)
 PlotDendrogram(seurat_obj, main=paste0('hdWGCNA Dendrogram'))
 dev.off()

#t <- as.data.frame(GetTOM(seurat_obj))
#t<-cbind(rownames(t),t)
#write.table(t,file=paste0(output_path,"07.tom.xlsx"),quote=F,sep="\t",row.names=F)

#dataExpr<-GetDatExpr(seurat_obj)
#TOM<-TOMsimilarityFromExpr(dataExpr, power = softPower,TOMType = 'signed',networkType = 'signed')
#write.table(TOM,file=paste0(output_path,"07.TOM.xlsx"),sep="\t",quote=F,row.names=F,col.names=F)

######################################################################
###### STEP6
################### 模块特征值(ME)
     #计算ME值
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj)

# compute all MEs in the full single-cell dataset
if( length(unique(seurat_obj$orig.ident)) > 1 )
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars= "orig.ident" # 根据样本去批次化 harmonize
) else {seurat_obj <- ModuleEigengenes(seurat_obj)}

### 统计并输出各模块的mRNA，lncRNA数量
#df <- as.data.frame(seurat_obj@misc$tutorial$wgcna_net$colors)
#mRNA<-table(df[!grepl("LNC",rownames(df)),])
#lncRNA<-table(df[grepl("LNC",rownames(df)),])
#mRNA <- data.frame(color=names(mRNA),mRNA=as.numeric(mRNA))
#lncRNA <- data.frame(color=names(lncRNA),lncRNA=as.numeric(lncRNA))
#t <- merge(mRNA,lncRNA,all=T)
#  t$allRNA<-(t$mRNA+t$lncRNA)
#   t<-cbind(rownames(t),t)
#    write.xlsx(t,file=paste0(output_path,"02.color.xlsx"),quote=F,sep="\t",row.names=F)

# t<-GetMEs(seurat_obj)
# t<-cbind(rownames(t),t)
# write.table(t,file=paste0(output_path,"03.MEs.xlsx"),sep="\t",quote=F,row.names=F,col.names=T)

######################################################################
###### STEP7
      #计算连接度
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'celltype', 
  group_name = clusters
)

plot5 <- PlotKMEs(seurat_obj, 
              ncol=5,
              n_hubs = 10, # number of hub genes to display
              text_size = 2,
              plot_widths = c(3, 2) # the relative width between the kME rank plot and the hub gene text
              )
ggsave(plot5,filename=paste0(output_path, "04_KMEs.pdf"),width=12,height=9)

t <- GetHubGenes(seurat_obj, n_hubs = 10)
t<-cbind(rownames(t),t)
write.table(t,file=paste0(output_path,"04.HubGene.xlsx"),sep="\t",quote=F,row.names=F,col.names=T)

t<-GetModules(seurat_obj)
t<-cbind(rownames(t),t)
write.table(t,file=paste0(output_path,"05.Modules.xlsx"),sep="\t",quote=F,row.names=F,col.names=T)

######################################################################
###### STEP8
################ Visualization
# FeaturePlot
pdf(file=paste0(output_path, "05_FeaturesPlot.pdf"));
ModuleFeaturePlot(
  seurat_obj,
  reduction = "umap",
  features = "hMEs",
  order_points = TRUE, # order so the points with highest hMEs are on top
  point_size = 0.5,
  raster_scale = 1)
dev.off();  
#####################

#############################
# Intermodule correlation, heatmap
pdf(file=paste0(output_path, "06_heatmap.pdf"))
pheatmap((ModuleCorrelogram(seurat_obj))$corr,display_numbers=T)
dev.off()    

MEs <- GetMEs(seurat_obj, harmonized=T)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# DotPlot
p <- DotPlot(seurat_obj, features = mods, group.by = 'celltype') +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high = 'red', mid = 'grey95', low = 'white') +
  labs(x = "Module", y = "Celltype") +
  theme(axis.title = element_text(size = 18),  # 调整标题字体大小
        axis.text = element_text(size = 12),   # 调整轴标签字体大小
        plot.title = element_text(size = 18))  # 调整整体标题字体大小     
ggsave(p, filename = paste0(output_path, "07_bubblePlot.png"), width = 12, height = 9)

# Vlnplot
plot_list <- lapply(mods, function(x) {
  print(x)
  p <- VlnPlot(
    seurat_obj,
    features = x,
    group.by = 'celltype',
    pt.size = 0 # don't show actual data points
  )
  # add box-and-whisker plots on top:
  p <- p + geom_boxplot(width=.25, fill='white')
  
  # change axis labels and remove legend:
  p <- p + xlab('') + ylab('hMEs') + NoLegend()
  p
})
pdf(file=paste0(output_path, "08_VlnPlot.pdf"))
plot_list
dev.off()

# single moduel plot
#remotes::install_github("igraph/rigraph@master")

# Visualizes the top hub genes for selected modules as a circular network plot
ModuleNetworkPlot(
  seurat_obj,
  mods = "all", # all modules are plotted.
  outdir = paste0(output_path,"ModuleNetworks"), # The directory where the plots will be stored.
  plot_size = c(6, 6),
  label_center = FALSE,
  edge.alpha = 0.25,
  vertex.label.cex = 0.5, # text size
  vertex.size = 8 # point size
)

###############################
# Composite network diagram
# hubgene network
pdf(file=paste0(output_path, "11_HubGeneNetworkPlot.pdf"))
HubGeneNetworkPlot(
  seurat_obj,
  mods = "all", # all modules are plotted.
  n_hubs = 3, 
  n_other= 5,
  edge_prop = 0.75,
  vertex.label.cex = 0.4
)
dev.off()

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color, # color each point by WGCNA module
   size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
ggsave(p,filename=paste0(output_path, "09_NET_UMAP.pdf"),width=12,height=9)

pdf(file=paste0(output_path, "10_Network.pdf"))
  ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=1 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)
dev.off()

saveRDS(seurat_obj,file=paste0(output_path,"hdWGCNA.rds"))