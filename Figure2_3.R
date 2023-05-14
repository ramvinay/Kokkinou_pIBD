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



#### Suppl. Figure S2-B

subset_data1 <- subset_data
clusters_louvain_subset2 <- vector("logical", length = ncol(subset_data1))
names(clusters_louvain_subset2) <- colnames(subset_data1)
head(subset_data1@meta.data)

clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="1"] <- "8"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="2"] <- "2"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="3"] <- "12"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="4"] <- "4"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="5"] <- "5"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="6"] <- "6"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="7"] <- "10"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="8"] <- "7"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="9"] <- "9"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="10"] <- "11"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="11"] <- "1"
clusters_louvain_subset2[subset_data1@meta.data$clusters_louvain_subset1=="12"] <- "3"

subset_data1[["clusters_louvain_subset2"]] <- clusters_louvain_subset2
head(subset_data1@meta.data)


cluster_cell <- as.data.frame(subset_data1@meta.data)
head(cluster_cell)

cluster_cell$Cell_ID <- row.names(cluster_cell)
head(cluster_cell)

unique(cluster_cell$clusters_louvain_subset2)
cluster1_cell_summary1 <- ddply(cluster_cell, .(clusters_louvain_subset2), summarize, Total_Cell_Count_In_Cluster = length(Cell_ID))
head(cluster1_cell_summary1)


head(cluster_cell)

######################################
cluster1_cell_summary3 <- ddply(cluster_cell, .(clusters_louvain_subset2,donor1), summarize, Each_Sample_Cell_Count_In_Cluster = length(Cell_ID))

cluster1_cell_summary3

all_cluster_cell_df1 <- merge(cluster1_cell_summary1,cluster1_cell_summary3,by="clusters_louvain_subset2")
head(all_cluster_cell_df1)
colnames(all_cluster_cell_df1)[3]<-"Donor"
all_cluster_cell_df1$Percent <- (all_cluster_cell_df1$Each_Sample_Cell_Count_In_Cluster/all_cluster_cell_df1$Total_Cell_Count_In_Cluster)*100



all_cluster_cell_df1$clusters_louvain_subset2 <- factor(all_cluster_cell_df1$clusters_louvain_subset2, levels = c("1","2","3","4","5","6","7","8","9","10","11","12"))
all_cluster_cell_df1$Donor <- factor(x = all_cluster_cell_df1$Donor, levels = c("D1", "D6","D10","D11","D16","D19")) # change the order of the factor levels

p4 <- ggplot() + geom_bar(aes(y = Percent, x = as.factor(clusters_louvain_subset2), fill = Donor), data = all_cluster_cell_df1,
                          stat="identity")+
  #scale_fill_manual(name="Sample",values=sample_col)+
  xlab("Cluster")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=18, face = "plain"),
        axis.title=element_text(size=20, face = "plain"),
        legend.text=element_text(size=16,face="plain"),
        legend.title=element_text(size=20,face="plain"))

p4


plot_file1 <- paste(output_prefix,"_Suppl_Figure_2S_B1.jpeg",sep = "")

jpeg(filename = plot_file1, width = 8, height = 6, units = "in",res=400)
print(p4)
dev.off()

all_cluster_cell_df <- cluster1_cell_summary3
colnames(all_cluster_cell_df) <- c("Cluster","Donor", "Cell_Count")
write.table(all_cluster_cell_df, file = paste(output_prefix,"_Suppl_Figure_2S_B1.txt",sep = ""), append = F, quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, qmethod = c("escape", "double"))

#################### Suppl_Figure_2S_B2
cluster1_cell_summary1 <- ddply(cluster_cell, .(donor1), summarize, Total_Cell_Count_In_Cluster = length(Cell_ID))
head(cluster1_cell_summary1)

cluster1_cell_summary3 <- ddply(cluster_cell, .(donor1,clusters_louvain_subset2), summarize, Each_Sample_Cell_Count_In_Cluster = length(Cell_ID))
cluster1_cell_summary3

all_cluster_cell_df1 <- merge(cluster1_cell_summary1,cluster1_cell_summary3,by="donor1")
head(all_cluster_cell_df1)
colnames(all_cluster_cell_df1)[1]<-"Donor"
colnames(all_cluster_cell_df1)[3]<-"Cluster"
all_cluster_cell_df1$Percent <- (all_cluster_cell_df1$Each_Sample_Cell_Count_In_Cluster/all_cluster_cell_df1$Total_Cell_Count_In_Cluster)*100



all_cluster_cell_df1$Cluster <- factor(all_cluster_cell_df1$Cluster, levels = c("1","2","3","4","5","6","7","8","9","10","11","12"))
all_cluster_cell_df1$Donor <- factor(x = all_cluster_cell_df1$Donor, levels = c("D1", "D6","D10","D11","D16","D19")) # change the order of the factor levels

p4 <- ggplot() + geom_bar(aes(y = Percent, x = as.factor(Donor), fill = Cluster), data = all_cluster_cell_df1,
                          stat="identity")+
  #scale_fill_manual(name="Sample",values=sample_col)+
  xlab("Cluster")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=18, face = "plain"),
        axis.title=element_text(size=20, face = "plain"),
        legend.text=element_text(size=16,face="plain"),
        legend.title=element_text(size=20,face="plain"))

p4


plot_file1 <- paste(output_prefix,"_Suppl_Figure_2S_B2.jpeg",sep = "")

jpeg(filename = plot_file1, width = 8, height = 6, units = "in",res=400)
print(p4)
dev.off()

all_cluster_cell_df <- cluster1_cell_summary3
colnames(all_cluster_cell_df) <- c("Donor","Cluster", "Cell_Count")
write.table(all_cluster_cell_df, file = paste(output_prefix,"_Suppl_Figure_2S_B2.txt",sep = ""), append = F, quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, qmethod = c("escape", "double"))


###################################### Suppl_Figure_2S_B3
cluster1_cell_summary3 <- ddply(cluster_cell, .(clusters_louvain_subset2,orig.ident), summarize, Each_Sample_Cell_Count_In_Cluster = length(Cell_ID))

cluster1_cell_summary3

all_cluster_cell_df1 <- merge(cluster1_cell_summary1,cluster1_cell_summary3,by="clusters_louvain_subset2")
head(all_cluster_cell_df1)
colnames(all_cluster_cell_df1)[3]<-"Sample"
all_cluster_cell_df1$Percent <- (all_cluster_cell_df1$Each_Sample_Cell_Count_In_Cluster/all_cluster_cell_df1$Total_Cell_Count_In_Cluster)*100




sample_col <- c("JMJ01_1n"="#feb743",
                "JMJ02_1i"="#cae63b",
                "JMJ07_19n"="#9c79b4",
                "JMJ08_19i"="#ff58c0",
                "JMJ03_16n"="#4ad2d6",
                "JMJ04_16i"="#823acb",
                "JMJ05_6n"="#e46c71",
                "JMJ06_6i"="#121163",
                "JMJ09_10n"="#f342d8",
                "JMJ10_10i"="#669966",
                "JMJ11_11n"="#ff954e",
                "JMJ12_11i"="#326ada")



all_cluster_cell_df1$clusters_louvain_subset2 <- factor(all_cluster_cell_df1$clusters_louvain_subset2, levels = c("1","2","3","4","5","6","7","8","9","10","11","12"))


p4 <- ggplot() + geom_bar(aes(y = Percent, x = as.factor(clusters_louvain_subset2), fill = Sample), data = all_cluster_cell_df1,
                          stat="identity")+
  scale_fill_manual(name="Sample",values=sample_col)+
  xlab("Cluster")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=18, face = "plain"),
        axis.title=element_text(size=20, face = "plain"),
        legend.text=element_text(size=16,face="plain"),
        legend.title=element_text(size=20,face="plain"))

p4



plot_file1 <- paste(output_prefix,"_Suppl_Figure_2S_B3.jpeg",sep = "")

jpeg(filename = plot_file1, width = 8, height = 6, units = "in",res=400)
print(p4)
dev.off()

all_cluster_cell_df <- cluster1_cell_summary3
colnames(all_cluster_cell_df) <- c("Cluster","Sample", "Cell_Count")
write.table(all_cluster_cell_df, file = paste(output_prefix,"_Suppl_Figure_2S_B3.txt",sep = ""), append = F, quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, qmethod = c("escape", "double"))



(B) Stacked bar plot showing the contribution of each donor (n = 6) to
the clusters shown in Figure 2B (left and middle). Stacked bar plot showing the contribution of each colon biopsy
sample (n = 12) to the clusters shown in Figure 2B (right). Cells from each donor were analyzed in independent
experiments.
