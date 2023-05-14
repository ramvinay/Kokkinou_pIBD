library(Seurat)
packageVersion("Seurat")
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(gridExtra)

outdir="Results"
dir.create(file.path(outdir), showWarnings = FALSE)
#setwd(outdir)

output_prefix <- paste(outdir, "/Kokkinou_pIBD_etal_2023",sep = "")

main_seurat_object <- "Data/New_Cluster_setting_one_noninflamed_sample_add_selected_clusterscluster7_removed1_1_40_k10_5_subset_data1.rds"

subset_data <- readRDS(main_seurat_object,refhook = NULL)

DimPlot(subset_data,label = T)

subset_data1 <- subset_data
subset_data1 <- ScaleData(subset_data1,assay='RNA', features = rownames(subset_data1))
dim(subset_data1@assays$RNA@scale.data)

#1-8 - CD4+ Tn/cm
#2-2 - NK pop 1
#3-12 - Activ. CD4+T/Tregs
#4-4 - CD8+ Tn/cm
#5-5 - CD8+ Teff
#6-6 - CD8+ Tem
#7-10 - CD4+ Trm
#8-7 - CD8+ Trm
#9-9 - CD4+ Teff
#10-11 - Tfh
#11-1 - ILC
#12-3 - NK pop 2

cluster_order <- c("11","2",	"12",	"4",	"5",	"6",	"8",	"1",	"9",	"7",	"10",	"3")


#### Annotating clusters
cell.types <- vector("logical", length = ncol(subset_data1))
names(cell.types) <- colnames(subset_data1)

cell.types[subset_data1@meta.data$clusters_louvain_subset1=="11"] <- "ILC"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="2"] <- "NK pop 1"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="12"] <- "NK pop 2"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="4"] <- "CD8+ Tn/cm"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="5"] <- "CD8+ Teff"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="6"] <- "CD8+ Tem"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="8"] <- "CD8+ Trm"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="1"] <- "CD4+ Tn/cm"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="9"] <- "CD4+ Teff"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="7"] <- "CD4+ Trm"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="10"] <- "Tfh"
cell.types[subset_data1@meta.data$clusters_louvain_subset1=="3"] <- "Activ. CD4+T/Tregs"

subset_data1[["cell.types"]] <- cell.types

Idents(subset_data1) <- subset_data1@meta.data$cell.types

##### Figure 2B
plot_file1 <- paste(output_prefix,"Figure_2B.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 8, height = 6, units = "in",res=400)
DimPlot(subset_data1,label = T, label.size = 3)
dev.off()



##### Figure 2C

gene_list <- c("ID2","IL7R","IL22","ZBTB16","PRF1","CD8B","CD27",
               "CCR7","SELL","TCF7","IFNG","GZMA", "KLRC2","ZNF683",
               "ITGA1","ITGAE","CD4","ICOS","CD69","CXCR5","PDCD1",
               "FAS","FOXP3", "CTLA4","IL10")


gene_list1 <- gene_list[length(gene_list):1]
gene_list1 <- unique(gene_list1)

unique(subset_data1@meta.data$cell.types)

cell.types.orders <- c("ILC","NK pop 1","NK pop 2","CD8+ Tn/cm",
                       "CD8+ Teff","CD8+ Tem","CD8+ Trm","CD4+ Tn/cm",
                       "CD4+ Teff","CD4+ Trm","Tfh","Activ. CD4+T/Tregs")
subset_data1@meta.data$cell.types <- factor(subset_data1@meta.data$cell.types, levels = cell.types.orders)

Idents(subset_data1) <- subset_data1@meta.data$cell.types

plot_file1 <- paste(output_prefix,"Figure_2C.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 9, height = 11, units = "in",res=400)


#pdf(file = plot_file1, width = 9, height = 11, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches
print(DotPlot(subset_data1, dot.scale = 12, features = gene_list1, col.min= -0.5, col.max=1.5,
              scale.max = 50, scale.min = 0)+ RotatedAxis()+ coord_flip()+
        #guides(colour=guide_legend(override.aes=list(size=8)))+
        theme(panel.background = element_rect(fill = "white", colour = "black"),
              axis.text=element_text(size=18, face = "plain",color = "#000000"),
              legend.text=element_text(size=18,face="plain",color = "#000000"),
              legend.title=element_text(size=20,face="plain",color = "#000000"),
              axis.text.y = element_text(face = "italic")))
dev.off()




######## Figure 3
module_genes_file <- "Data/Modules chategorized genes for fig.xlsx"
gene_df <- openxlsx::read.xlsx(module_genes_file,sheet = 5)
head(gene_df)


######## Figure 3C
colnames(gene_df)
column_list <- c(1,4,11)

j<-1
plot_list <- list()
for (i in column_list){
  #i<-5
  #i<-1
  print(i)
  module_name <- colnames(gene_df)[i]
  module_name <- gsub("\\."," ",module_name)
  #out_prefix2 <- paste(out_prefix1,"_",module_name,sep = "")
  gene_list <- subset(gene_df,  !is.na(gene_df[,i]))
  gene_list <- as.character(gene_list[,i])
  
  #gene_list <- append(gene_list, "CD83")
  #gene_list <- append(gene_list, "CD300LF")
  heatmap_genes1 <- gene_list[length(gene_list):1]
  heatmap_genes1 <- unique(heatmap_genes1)
  #plot_file1 <- paste(output_prefix,"_DotPlot_fixed_expre_scaled2.pdf",sep = "")
  #pdf(file = plot_file1, width = 12, height = 5, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches
  
  p1 <- DotPlot(subset_data1, dot.scale = 15, features = heatmap_genes1, col.min= -0.5, col.max=1.5,
                scale.max = 50, scale.min = 0)+ RotatedAxis()+ coord_flip()+
    xlab(module_name)+ ylab('')+

          #guides(colour=guide_legend(override.aes=list(size=8)))+
          theme(panel.background = element_rect(fill = "white", colour = "black"),
                axis.text=element_text(size=18, face = "plain",color = "#000000"),
                legend.text=element_text(size=18,face="plain",color = "#000000"),
                legend.title=element_text(size=20,face="plain",color = "#000000"),
                axis.text.y = element_text(face = "italic"))
  #dev.off()
  #axis.title.x = element_blank()
  plot_list[[j]] <- p1
  j<-j+1
}



plot_file1 <- paste(output_prefix,"Figure_3C.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 11, height = 14, units = "in",res=400)
ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol=1)
dev.off()



######## Figure 3D
colnames(gene_df)
column_list <- c(3,6)

j<-1
plot_list <- list()
for (i in column_list){
  #i<-5
  #i<-1
  print(i)
  module_name <- colnames(gene_df)[i]
  module_name <- gsub("\\."," ",module_name)
  #out_prefix2 <- paste(out_prefix1,"_",module_name,sep = "")
  gene_list <- subset(gene_df,  !is.na(gene_df[,i]))
  gene_list <- as.character(gene_list[,i])
  
  #gene_list <- append(gene_list, "CD83")
  #gene_list <- append(gene_list, "CD300LF")
  heatmap_genes1 <- gene_list[length(gene_list):1]
  heatmap_genes1 <- unique(heatmap_genes1)
  #plot_file1 <- paste(output_prefix,"_DotPlot_fixed_expre_scaled2.pdf",sep = "")
  #pdf(file = plot_file1, width = 12, height = 5, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches
  
  p1 <- DotPlot(subset_data1, dot.scale = 15, features = heatmap_genes1, col.min= -0.5, col.max=1.5,
                scale.max = 50, scale.min = 0)+ RotatedAxis()+ coord_flip()+
    xlab(module_name)+ ylab('')+
    
    #guides(colour=guide_legend(override.aes=list(size=8)))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text=element_text(size=18, face = "plain",color = "#000000"),
          legend.text=element_text(size=18,face="plain",color = "#000000"),
          legend.title=element_text(size=20,face="plain",color = "#000000"),
          axis.text.y = element_text(face = "italic"))
  #dev.off()
  #axis.title.x = element_blank()
  plot_list[[j]] <- p1
  j<-j+1
}



plot_file1 <- paste(output_prefix,"Figure_3D.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 11, height = 8, units = "in",res=400)
ggarrange(plot_list[[1]], plot_list[[2]], ncol=1)
dev.off()



######## Figure 3E
colnames(gene_df)
column_list <- c(9)

j<-1
plot_list <- list()
for (i in column_list){
  #i<-5
  #i<-1
  print(i)
  module_name <- colnames(gene_df)[i]
  module_name <- gsub("\\."," ",module_name)
  #out_prefix2 <- paste(out_prefix1,"_",module_name,sep = "")
  gene_list <- subset(gene_df,  !is.na(gene_df[,i]))
  gene_list <- as.character(gene_list[,i])
  
  #gene_list <- append(gene_list, "CD83")
  #gene_list <- append(gene_list, "CD300LF")
  heatmap_genes1 <- gene_list[length(gene_list):1]
  heatmap_genes1 <- unique(heatmap_genes1)
  #plot_file1 <- paste(output_prefix,"_DotPlot_fixed_expre_scaled2.pdf",sep = "")
  #pdf(file = plot_file1, width = 12, height = 5, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches
  
  p1 <- DotPlot(subset_data1, dot.scale = 15, features = heatmap_genes1, col.min= -0.5, col.max=1.5,
                scale.max = 50, scale.min = 0)+ RotatedAxis()+ coord_flip()+
    xlab(module_name)+ ylab('')+
    
    #guides(colour=guide_legend(override.aes=list(size=8)))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text=element_text(size=18, face = "plain",color = "#000000"),
          legend.text=element_text(size=18,face="plain",color = "#000000"),
          legend.title=element_text(size=20,face="plain",color = "#000000"),
          axis.text.y = element_text(face = "italic"))
  #dev.off()
  #axis.title.x = element_blank()
  plot_list[[j]] <- p1
  j<-j+1
}



plot_file1 <- paste(output_prefix,"Figure_3E.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 11, height = 8, units = "in",res=400)
ggarrange(plot_list[[1]], ncol=1)
dev.off()


######## Figure 3F
colnames(gene_df)
column_list <- c(2)

j<-1
plot_list <- list()
for (i in column_list){
  #i<-5
  #i<-1
  print(i)
  module_name <- colnames(gene_df)[i]
  module_name <- gsub("\\."," ",module_name)
  #out_prefix2 <- paste(out_prefix1,"_",module_name,sep = "")
  gene_list <- subset(gene_df,  !is.na(gene_df[,i]))
  gene_list <- as.character(gene_list[,i])
  
  #gene_list <- append(gene_list, "CD83")
  #gene_list <- append(gene_list, "CD300LF")
  heatmap_genes1 <- gene_list[length(gene_list):1]
  heatmap_genes1 <- unique(heatmap_genes1)
  #plot_file1 <- paste(output_prefix,"_DotPlot_fixed_expre_scaled2.pdf",sep = "")
  #pdf(file = plot_file1, width = 12, height = 5, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches
  
  p1 <- DotPlot(subset_data1, dot.scale = 15, features = heatmap_genes1, col.min= -0.5, col.max=1.5,
                scale.max = 50, scale.min = 0)+ RotatedAxis()+ coord_flip()+
    xlab(module_name)+ ylab('')+
    
    #guides(colour=guide_legend(override.aes=list(size=8)))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text=element_text(size=18, face = "plain",color = "#000000"),
          legend.text=element_text(size=18,face="plain",color = "#000000"),
          legend.title=element_text(size=20,face="plain",color = "#000000"),
          axis.text.y = element_text(face = "italic"))
  #dev.off()
  #axis.title.x = element_blank()
  plot_list[[j]] <- p1
  j<-j+1
}



plot_file1 <- paste(output_prefix,"Figure_3F.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 11, height = 8, units = "in",res=400)
ggarrange(plot_list[[1]], ncol=1)
dev.off()


######## Figure 3G
colnames(gene_df)
column_list <- c(7)

j<-1
plot_list <- list()
for (i in column_list){
  #i<-5
  #i<-1
  print(i)
  module_name <- colnames(gene_df)[i]
  module_name <- gsub("\\."," ",module_name)
  #out_prefix2 <- paste(out_prefix1,"_",module_name,sep = "")
  gene_list <- subset(gene_df,  !is.na(gene_df[,i]))
  gene_list <- as.character(gene_list[,i])
  
  #gene_list <- append(gene_list, "CD83")
  #gene_list <- append(gene_list, "CD300LF")
  heatmap_genes1 <- gene_list[length(gene_list):1]
  heatmap_genes1 <- unique(heatmap_genes1)
  #plot_file1 <- paste(output_prefix,"_DotPlot_fixed_expre_scaled2.pdf",sep = "")
  #pdf(file = plot_file1, width = 12, height = 5, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches
  
  p1 <- DotPlot(subset_data1, dot.scale = 15, features = heatmap_genes1, col.min= -0.5, col.max=1.5,
                scale.max = 50, scale.min = 0)+ RotatedAxis()+ coord_flip()+
    xlab(module_name)+ ylab('')+
    
    #guides(colour=guide_legend(override.aes=list(size=8)))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text=element_text(size=18, face = "plain",color = "#000000"),
          legend.text=element_text(size=18,face="plain",color = "#000000"),
          legend.title=element_text(size=20,face="plain",color = "#000000"),
          axis.text.y = element_text(face = "italic"))
  #dev.off()
  #axis.title.x = element_blank()
  plot_list[[j]] <- p1
  j<-j+1
}



plot_file1 <- paste(output_prefix,"Figure_3G.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 11, height = 8, units = "in",res=400)
ggarrange(plot_list[[1]], ncol=1)
dev.off()



######## Figure 3H
colnames(gene_df)
column_list <- c(5,8)

j<-1
plot_list <- list()
for (i in column_list){
  #i<-5
  #i<-1
  print(i)
  module_name <- colnames(gene_df)[i]
  module_name <- gsub("\\."," ",module_name)
  #out_prefix2 <- paste(out_prefix1,"_",module_name,sep = "")
  gene_list <- subset(gene_df,  !is.na(gene_df[,i]))
  gene_list <- as.character(gene_list[,i])
  
  #gene_list <- append(gene_list, "CD83")
  #gene_list <- append(gene_list, "CD300LF")
  heatmap_genes1 <- gene_list[length(gene_list):1]
  heatmap_genes1 <- unique(heatmap_genes1)
  #plot_file1 <- paste(output_prefix,"_DotPlot_fixed_expre_scaled2.pdf",sep = "")
  #pdf(file = plot_file1, width = 12, height = 5, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches
  
  p1 <- DotPlot(subset_data1, dot.scale = 15, features = heatmap_genes1, col.min= -0.5, col.max=1.5,
                scale.max = 50, scale.min = 0)+ RotatedAxis()+ coord_flip()+
    xlab(module_name)+ ylab('')+
    
    #guides(colour=guide_legend(override.aes=list(size=8)))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text=element_text(size=18, face = "plain",color = "#000000"),
          legend.text=element_text(size=18,face="plain",color = "#000000"),
          legend.title=element_text(size=20,face="plain",color = "#000000"),
          axis.text.y = element_text(face = "italic"))
  #dev.off()
  #axis.title.x = element_blank()
  plot_list[[j]] <- p1
  j<-j+1
}



plot_file1 <- paste(output_prefix,"Figure_3H.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 11, height = 10, units = "in",res=400)
ggarrange(plot_list[[1]],plot_list[[2]], ncol=1)
dev.off()




######## Figure 3I
colnames(gene_df)
column_list <- c(10)

j<-1
plot_list <- list()
for (i in column_list){
  #i<-5
  #i<-1
  print(i)
  module_name <- colnames(gene_df)[i]
  module_name <- gsub("\\."," ",module_name)
  #out_prefix2 <- paste(out_prefix1,"_",module_name,sep = "")
  gene_list <- subset(gene_df,  !is.na(gene_df[,i]))
  gene_list <- as.character(gene_list[,i])
  
  #gene_list <- append(gene_list, "CD83")
  #gene_list <- append(gene_list, "CD300LF")
  heatmap_genes1 <- gene_list[length(gene_list):1]
  heatmap_genes1 <- unique(heatmap_genes1)
  #plot_file1 <- paste(output_prefix,"_DotPlot_fixed_expre_scaled2.pdf",sep = "")
  #pdf(file = plot_file1, width = 12, height = 5, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches
  
  p1 <- DotPlot(subset_data1, dot.scale = 15, features = heatmap_genes1, col.min= -0.5, col.max=1.5,
                scale.max = 50, scale.min = 0)+ RotatedAxis()+ coord_flip()+
    xlab(module_name)+ ylab('')+
    
    #guides(colour=guide_legend(override.aes=list(size=8)))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text=element_text(size=18, face = "plain",color = "#000000"),
          legend.text=element_text(size=18,face="plain",color = "#000000"),
          legend.title=element_text(size=20,face="plain",color = "#000000"),
          axis.text.y = element_text(face = "italic"))
  #dev.off()
  #axis.title.x = element_blank()
  plot_list[[j]] <- p1
  j<-j+1
}



plot_file1 <- paste(output_prefix,"Figure_3I.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 11, height = 8, units = "in",res=400)
ggarrange(plot_list[[1]], ncol=1)
dev.off()