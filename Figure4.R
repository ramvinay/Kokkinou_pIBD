library(Seurat)
packageVersion("Seurat")
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(gridExtra)
# remotes::install_github("czarnewski/nicerplots")
#install.packages("remotes")
library(niceRplots)
library(slingshot)
setwd("/Users/rampan/OneDrive-KI.SE/OneDrive-KI.SE/Mac/Documents/proj/KI/Jenny/10Xseq_data/github")

outdir="Results"
dir.create(file.path(outdir), showWarnings = FALSE)
#setwd(outdir)

output_prefix <- paste(outdir, "/Kokkinou_pIBD_etal_2023",sep = "")

figure4_seurat_object <- "Data/Figure4_ILC_subcluster/New_Cluster_setting_with_JMJ09_10n_2022_03_28_cluster11_res_0.6.rds"
subset_data_sub <- readRDS(figure4_seurat_object,refhook = NULL)

plot_file1 <- paste(output_prefix,"Figure_4A.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 7, height = 7, units = "in",res=400)

print(DimPlot(subset_data_sub, reduction = "umap",label = T,label.size = 8, repel=TRUE,pt.size=0.5)+
        guides(colour=guide_legend(override.aes=list(size=8)))+
        theme(panel.background = element_rect(fill = "white", colour = "black"),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text=element_text(size=18, face = "bold"),
              axis.title=element_text(size=20, face = "bold"),
              legend.text=element_text(size=16,face="bold"),
              legend.title=element_text(size=20,face="bold")))
dev.off()




gene_list <- c("IL22", "CSF2", "IL4I1", "TSC22D3", "DUSP1", "ZFP36L2",
               "HLA-DPA1","HLA-DRA", "HLA-DQB1", "KLF2", "TCF7", "SELL")

gene_list1 <- gene_list[length(gene_list):1]
gene_list1 <- unique(gene_list1)
#plot_file1 <- paste(out_prefix2,"_DotPlot_ILC_cluster_p3.pdf",sep = "")
#pdf(file = plot_file1, width = 10, height = 7, family = "Helvetica") # defaults to 8 x 8 inches

plot_file1 <- paste(output_prefix,"Figure_4B.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 10, height = 7, units = "in",res=400)


print(DotPlot(subset_data_sub, dot.scale = 15, features = gene_list1,
)+ RotatedAxis()+ coord_flip()+
  #guides(colour=guide_legend(override.aes=list(size=8)))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text=element_text(size=18, face = "plain",color = "#000000"),
        legend.text=element_text(size=18,face="plain",color = "#000000"),
        legend.title=element_text(size=20,face="plain",color = "#000000"),
        axis.text.y = element_text(face = "italic")))
dev.off()


############### Figure 3C to 3G
gene_list <- c("IL22", "CSF2", "IL4I1", "NCR2", "HLA-DPA1","HLA-DRA", "HLA-DQB1",
               "TSC22D3", "DUSP1", "CD69", "SELL")

for (g1 in gene_list){
  plot_file1 <- paste(output_prefix,"Figure_4C_to_4G_",g1,".jpeg",sep = "")
  jpeg(filename = plot_file1, width = 4, height = 4, units = "in",res=400)
  print(plot_feat(subset_data_sub,feat = g1, add_legend=T,cex = 0.4))
  dev.off()
}


########## Figure 4H
seu <- subset_data_sub

head(seu@meta.data)
cl <- seu$seurat_clusters
rd <- Embeddings(seu, "umap")
pto <- slingshot(rd, cl, start.clus = '4')

sds <- SlingshotDataSet(pto)
reducedDim(sds)

dim(rd)
length(cl)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

#cell_colors <- cell_pal(seu$cell.types, brewer_pal("qual", "Set3"))
cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())



plot_file1 <- paste(output_prefix,"_Figure_4H.jpeg",sep = "")
jpeg(filename = plot_file1, width = 4, height = 4, units = "in",res=400)
#plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.2)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
dev.off()


#### Suppl. Figure S4-G

pto <- slingshot(rd, cl, start.clus = '1')

sds <- SlingshotDataSet(pto)
reducedDim(sds)

dim(rd)
length(cl)
#cell_colors <- cell_pal(seu$cell.types, brewer_pal("qual", "Set3"))
cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())



plot_file1 <- paste(output_prefix,"_Suppl_Figure_4S_G.jpeg",sep = "")
jpeg(filename = plot_file1, width = 4, height = 4, units = "in",res=400)
#plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.2)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
dev.off()





#### Suppl. Figure S4-B

gene_list <- c("RORC", "TBX21", "EOMES","NCR2", "KIT", "TNFSF13B", "SOX4")
gene_list1 <- gene_list[length(gene_list):1]
gene_list1 <- unique(gene_list1)
plot_file1 <- paste(output_prefix,"_Suppl_Figure_4S_B.jpeg",sep = "")
jpeg(filename = plot_file1, width = 9, height = 5, units = "in",res=400)

print(DotPlot(subset_data_sub, dot.scale = 15, features = gene_list1,
)+ RotatedAxis()+ coord_flip()+
  #guides(colour=guide_legend(override.aes=list(size=8)))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text=element_text(size=18, face = "plain",color = "#000000"),
        legend.text=element_text(size=18,face="plain",color = "#000000"),
        legend.title=element_text(size=20,face="plain",color = "#000000"),
        axis.text.y = element_text(face = "italic")))
dev.off()


#### Suppl. Figure S4-F
g1 <- "RORC"
plot_file1 <- paste(output_prefix,"_Suppl_Figure_4S_F.jpeg",sep = "")
jpeg(filename = plot_file1, width = 6, height = 6, units = "in",res=400)

print(plot_feat(subset_data_sub,feat = g1, add_legend=T,cex = 0.5))
dev.off()

#### Suppl. Figure S4-C & D
