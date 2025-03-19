library(monocle3)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(SeuratData)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(patchwork)
cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999","#FF34B3", "#00F5FF", "#FF6A6A", "#7FFFD4", "#AB82FF", "#90EE90", "#00CD00", "#008B8B", "#6495ED", "#FFC1C1", "#CD5C5C", "#8B008B", "#FF3030", "#7CFC00", "#708090", "#CC66FF", "#CC6600", "#990033", "#9900CC","#660000", "#000000")

setwd("D:/Desktop/lincRNA修改补充/文稿/重投稿/拟时序/拟南芥/根皮层/")
sample_name <- "Arabidopsis_root_cortex"

############################################ 读取文件
obj1 <- readRDS("D:/Desktop/lincRNA修改补充/文稿/重投稿/拟时序/RDS/Arabidopsis_root.rds")

DimPlot(obj1, group.by="celltype",cols = cols, pt.size = 1 , label =T)

####################### Create cds
counts <- GetAssayData(object = obj1, assay = "SCT", slot = "counts")
cell.metadata <- obj1[[]]
feature.metadata <- obj1[["SCT"]][[]]
feature.metadata$gene_short_name <- rownames(x = feature.metadata)
cds <- new_cell_data_set(counts,
                         cell_metadata = cell.metadata,
                         gene_metadata = feature.metadata)

#################### target cell types
valid_cells <- row.names(subset(pData(cds), is.element(celltype,  c("Quiescent center","Cortex","Endodermis"))))
cds_sub <- cds[,valid_cells]

## preprocess
cds_sub <- preprocess_cds(cds_sub, num_dim = 50)
cds_sub <- align_cds(cds_sub, alignment_group = "orig.ident")
cds_sub <- reduce_dimension(cds_sub)

######### 
plot_cells(cds_sub, group_label_size = 6, cell_size = 0.8, label_cell_groups = F,color_cells_by = "celltype") +
    theme_void() + theme(legend.position = c(0.80, 0.2)) + 
    theme(legend.text = element_text(size = 18),legend.title = element_text(size = 18)) + 
    scale_color_manual(values = cols)
ggsave("1_celltype.png",width = 7,height = 7);

#######
cds_sub <- cluster_cells(cds_sub,  cluster_method = "louvain")
cds_sub <- learn_graph(cds_sub, use_partition = F, close_loop = F, learn_graph_control = list(minimal_branch_len = 15))
# genes = c("Os12g0207000") ### Some of the genes specified in the article are used to evaluate the root node
# plot_cells(cds_sub, genes= genes, show_trajectory_graph=FALSE, #
#            label_cell_groups=FALSE,  label_leaves=FALSE , cell_size = 1) +  scale_color_gradientn(colors = c('grey','yellow','orange','red'))
# ggsave("2_genes_DimPlot.png",width = 7,height = 7)

cds_sub <- order_cells(cds_sub) # choose the root node

plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4,
           cell_size = 1,
           label_principal_points = F, #TRUE
           trajectory_graph_segment_size = 1) + 
  theme_void() + theme(legend.position = c(0.8, 0.2)) + 
  theme(legend.text = element_text(size = 18),legend.title = element_text(size = 18))
ggsave("2_pseudotime.png",width = 7,height = 7);

saveRDS(cds_sub,file= paste0(sample_name,"_pseudotime.rds"))
######## Track_genes
Track_genes <- graph_test(cds_sub, neighbor_graph="principal_graph", cores=32)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
Track_genes <- Track_genes[order(Track_genes$morans_I,decreasing = T),]
table(grepl("LNC",Track_genes$gene_short_name))
write.csv(Track_genes, file = paste0(sample_name ,"_Trajectory_genes.csv"), row.names = F) ## 

####################################################################################################
######################### Visualization of trajetory-dependent lincRNAs
plot_gene <- function(cds_sub , gene, group.by = "celltype", dir = "plot_gene" ){
  if(!dir.exists(dir)){dir.create(dir)}
  p <- plot_genes_in_pseudotime(cds_sub[gene,], color_cells_by = group.by, min_expr = 0.5) 
  p1 <- p + xlab("Pseudotime") +
    theme(text = element_text(size = 20, hjust = 0.5),
          plot.title = element_text(size = 20, hjust = 0.5, face = "italic"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
           legend.position = c("none"),  # 图例字体大小及位置
          axis.title = element_text(size = 20),               # 设置坐标轴标题字体大小
          axis.text = element_text(size = 20)) +              # 设置坐标轴刻度标签字体大小) 
    guides(color = guide_legend(override.aes = list(size = 5))) + # 图例点大小
    scale_color_manual(values = cols)   
  p2 <- plot_cells(cds_sub, genes = gene, show_trajectory_graph = FALSE,
                   label_cell_groups = FALSE, label_leaves = FALSE, cell_size = 1) + 
    theme_void() + # 隐藏坐标轴及其标签
    theme(
      text = element_text(size = 20, hjust = 0.5),
      axis.title = element_blank(),  # 隐藏坐标轴标题
      axis.text = element_blank(),   # 隐藏坐标轴刻度标签
      legend.text = element_text(size = 20), # 设置图例文本大小
      legend.title = element_text(size = 20), # 设置图例标题大小
      legend.position = "right" +  # 可选：调整图例的位置
        theme(plot.title = element_text(size = 20, hjust = 0.5, face = "italic"))
        ) +
    scale_color_gradientn(colors = c('grey', 'green', 'yellow', 'orange', 'red'))
  ggsave(p1,filename = paste0(dir,"/",gene,"_p1.png"), width = 7, height = 7)
  ggsave(p2,filename = paste0(dir,"/",gene,"_p2.png"), width = 7 ,height = 7)
}

#############
setwd("D:/Desktop/lincRNA修改补充/文稿/重投稿/拟时序/Plot_genes_in_pseudoptime/")

df <- data.frame(data1 = c("../水稻/根毛表皮/Rice_root_hair_pseudotime.rds","../水稻/根毛表皮/Rice_root_hair_Trajectory_genes.csv","rice_root_hair"),
                 data2 = c("../水稻/皮层厚壁/Rice_root_cortex_pseudotime.rds","../水稻/皮层厚壁/Rice_root_cortex_Trajectory_genes.csv","rice_root_cortex"),
                 data3 = c("../水稻/雌蕊/Rice_pistils_ovule_pseudotime.rds","../水稻/雌蕊/Rice_pistils_ovule_Trajectory_genes.csv","rice_pistils"),
                 data4 = c("../水稻/叶肉/Rice_leaf_pseudotime.rds","../水稻/叶肉/Rice_leaf_Trajectory_genes.csv","rice_leaf"),
                 data5 = c("../拟南芥/根毛/Arabidopsis_root_hair_pseudotime.rds","../拟南芥/根毛/Arabidopsis_root_hair_Trajectory_genes.csv","Arabidopsis_root_hair"),
                 data6 = c("../拟南芥/韧皮部/Arabidopsis_root_phloem_pseudotime.rds","../拟南芥/韧皮部/Arabidopsis_root_phloem_Trajectory_genes.csv","Arabidopsis_phloem"),
                 data7 = c("../番茄/leaf_initiation_and_vasculature2/tomato_pseudotime.rds","../番茄/leaf_initiation_and_vasculature2/tomato_Trajectory_genes.csv","tomatos")) 
rownames(df) <- c("rds_path","gene_path","data_name")
df <- as.data.frame(t(df))
                 
rds <- track_genes <- list();
for(i in 1:7){
  name = df$data_name[i];
  rds[[name]] <- readRDS(file = df$rds_path[i]);
  track_genes[[name]] <- read.csv(file = df$gene_path[i])
  }

####################################################################  plots
for(name in df$data_name){
  data <- track_genes[[name]] 
  data <- data[data$morans_I >= 0.2, ] # morans_I >= 0.2
  data <- data[grepl("LNC",data$gene_short_name),] # 
  lncRNAs <- data$gene_short_name
    #############
  if(length(lncRNAs) >= 0 ){
    for(gene in lncRNAs){
      plot_gene(cds_sub = rds[[name]], gene =  gene, group.by = "celltype", dir = name) 
    }
  }
}