library(Seurat)
packageVersion("Seurat")
library(ggplot2)

outdir="Results"
dir.create(file.path(outdir), showWarnings = FALSE)
#setwd(outdir)

output_prefix <- paste(outdir, "/Kokkinou_pIBD_etal_2023",sep = "")

main_seurat_object <- "Data/New_Cluster_setting_one_noninflamed_sample_add_selected_clusterscluster7_removed1_1_40_k10_5_subset_data1.rds"

subset_data <- readRDS(main_seurat_object,refhook = NULL)
subset_data <- ScaleData(subset_data,assay='RNA', features = rownames(subset_data))

DimPlot(subset_data,label = T)



jakob_table_df <- as.data.frame(read.csv("Data/Complete_metadata.csv",sep = 
                                           ",",header = T))

row.names(jakob_table_df) <- jakob_table_df$X
head(jakob_table_df)
colnames(jakob_table_df)

jakob_table_df1 <- jakob_table_df[,c(16,17,18,19)]

colnames(jakob_table_df)
head(jakob_table_df1)
colnames(jakob_table_df1)



subset_data <- AddMetaData(object=subset_data, metadata=jakob_table_df1)


head(subset_data@meta.data)
subset_data1 <- subset_data
Idents(subset_data1) <- subset_data1@meta.data$neighbor_smooth_CHASTea


head(jakob_table_df1)
jakob_table_df2 <- jakob_table_df1
jakob_table_df2$CHASTea<-NULL
jakob_table_df2$inflam<-NULL
jakob_table_df2$NonInflam<-NULL
jakob_table_df2$above5samplesInNeighbors<-NULL
colnames(jakob_table_df2) <- "hist.score"
head(jakob_table_df2)
row.names(t(jakob_table_df2))


subset_data1 <- subset_data
subset_data1[["HIST"]] <- CreateAssayObject(counts = t(jakob_table_df2))

dim(jakob_table_df2)

DefaultAssay(subset_data1) <- "HIST"




11-1
2-2
12-3
8-7
5-5
6-6
4-4
1-8
10-11
3-12
7-10
9-9

#### Annotating clusters
cell.types <- vector("logical", length = ncol(subset_data1))
names(cell.types) <- colnames(subset_data1)
head(subset_data1@meta.data)
subset_data1

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

Idents(subset_data1) <- subset_data1@meta.data$clusters_louvain_subset2
#row.names(subset_data1)
library(RColorBrewer)
#FeaturePlot(subset_data1, features = "hist.score", order = T,label=T) & 
#  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#plot_file1 <- paste(output_prefix,"_hist.score1.pdf",sep = "")
#pdf(file = plot_file1, width = 9, height = 8, family = "Helvetica",useDingbats = F) # defaults to 8 x 8 inches

plot_file1 <- paste(output_prefix,"_Figure_5C.jpeg",sep = "")
jpeg(filename = plot_file1, width = 9, height = 8, units = "in",res=400)


FeaturePlot(subset_data1, features = "hist.score", order = T,label=T, label.size = 6)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
  
  #guides(colour=guide_legend(override.aes=list(size=8)))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=18, face = "plain"),
        axis.title=element_text(size=20, face = "plain"),
        legend.text=element_text(size=16,face="plain"),
        legend.title=element_text(size=20,face="plain"))
dev.off()


#### Figure 5G-H


plot_dotplot <- function(subset_data1,gene_list,out_prefix1,op,height1, type1, cluster_number){
  
  jakob_table_df <- as.data.frame(read.csv("Data/Complete_metadata.csv",sep = 
                                             ",",header = T))
  
  
  
  row.names(jakob_table_df) <- jakob_table_df$X
  head(jakob_table_df)
  colnames(jakob_table_df)
  #jakob_table_df1 <- jakob_table_df[,c(16,17,18,19,20)]
  jakob_table_df1 <- jakob_table_df[,c(16,17,18,19,20)]
  
  colnames(jakob_table_df)
  head(jakob_table_df1)
  colnames(jakob_table_df1)
  jakob_table_df1$Hist.score <- jakob_table_df1$neighbor_smooth_CHASTea
  dim(jakob_table_df1)
  
  
  subset_data1 <- AddMetaData(object=subset_data1, metadata=jakob_table_df1)
  head(subset_data1@meta.data)
  heatmap_genes1 <- gene_list[length(gene_list):1]
  heatmap_genes1 <- unique(heatmap_genes1)
  
  subset_data2 <- subset(subset_data1, idents=c(cluster_number), invert=FALSE)
  subset_data2@meta.data$barcode <- row.names(subset_data2@meta.data)
  unique(subset_data2@meta.data$clusters_louvain_subset1)
  unique(subset_data2@meta.data$NonInflam)
  
  table(subset_data2@meta.data$NonInflam)
  if (type1=="Noninflamed"){
    
    subset_data2@meta.data$NonInflam <- factor(subset_data2@meta.data$NonInflam, levels=c("TRUE","FALSE"))
    subset_data2@meta.data$NonInflam <- ifelse(subset_data2@meta.data$NonInflam=="TRUE","LEAST","REST")
    
    Idents(subset_data2) <- subset_data2@meta.data$NonInflam

    #plot_file1 <- paste(out_prefix1,"_VlnPlot_",op,"_TRUE_FALSE.pdf",sep = "")
    plot_file1 <- paste(output_prefix,"_Figure_5_G_H_",cluster_number,"_Noninflam.jpeg",sep = "")
    
    #pdf(file = plot_file1, width = 6, height = height1, family = "Helvetica") # defaults to 8 x 8 inches
    jpeg(filename = plot_file1, width = 4, height = height1, units = "in",res=400)
    
    
    
    print(VlnPlot(
      subset_data2,
      features=gene_list,
      stack = TRUE,
      combine = TRUE,
      flip = TRUE
    )+theme(legend.position = 'none'))
    dev.off()
    
  } else {
    subset_data2@meta.data$Inflam <- factor(subset_data2@meta.data$Inflam, levels=c("TRUE","FALSE"))
    subset_data2@meta.data$Inflam <- ifelse(subset_data2@meta.data$Inflam=="TRUE","MOST","REST")
    subset_data2@meta.data$Inflam <- factor(subset_data2@meta.data$Inflam, levels = c("MOST","REST") )
    Idents(subset_data2) <- subset_data2@meta.data$Inflam
    
    #paste(out_prefix1,"_Figure_5_G_H_",op,"_TRUE_FALSE.pdf",sep = "")
    
    plot_file1 <- paste(output_prefix,"_Figure_5_G_H_",cluster_number,"_Inflam.jpeg",sep = "")
    jpeg(filename = plot_file1, width = 4, height = height1, units = "in",res=400)
    
    print(VlnPlot(
      subset_data2,
      features=gene_list,
      stack = TRUE,
      combine = TRUE,
      flip = TRUE
    )+theme(legend.position = 'none'))
    dev.off()
    
  }
}

#plot_dotplot <- function(subset_data1,gene_list,out_prefix1,op,height1, type1, cluster_number){

subset_data1 <- subset_data

out_prefix1 <- output_prefix
op <- "CL1_noninflamed_new_gene_order"
gene_list <- c("PTGER4","IL16","CD164","LEF1","SELL","CCR7","TCF7")
plot_dotplot(subset_data1,gene_list,out_prefix1,op,7, "Noninflamed",1)

op <- "CL4_noninflamed_new_gene_order"
gene_list <- c("NR4A2","PTGER4","IL16","CD164","LEF1","SELL","CCR7","ANXA6","ITGB2","SOCS3")
plot_dotplot(subset_data1,gene_list,out_prefix1,op,9, "Noninflamed",4)




op <- "CL8_noninflamed_new_gene_order"
gene_list <- c("NR4A2",
               "JUND",
               "PTGER4",
               "ICOS",
               "CCL4",
               "KLRB1",
               "NFKBID",
               "NFKBIZ",
               "RORA",
               "MYC",
               "TNF",
               "IFNG",
               "TNFSF9",
               "REL")

plot_dotplot(subset_data1,gene_list,out_prefix1,op,12, "Noninflamed",8)


op <- "CLl11_noninflamed_1_new_gene_order"
gene_list <- c("NR4A2",
               "JUND","ICOS","CTNNB1","CD83","MAFF")
plot_dotplot(subset_data1,gene_list,out_prefix1,op,6, "Noninflamed",11)


op <- "CL7_noninflamed_new_gene_order"
gene_list <- c("NR4A2","JUND","PTGER4","ICOS")
plot_dotplot(subset_data1,gene_list,out_prefix1,op,4, "Noninflamed",7)





op <- "Cl2_inflamed"
gene_list <- c("GNLY", "CD81", "CLEC2B")
plot_dotplot(subset_data1,gene_list,out_prefix1,op,4, "Inflamed",2)

op <- "Cl3_inflamed"
gene_list <- c("LAIR2","CARD16","SOD1","PRDX2","PRDX1","SELK","TNFRSF18")
plot_dotplot(subset_data1,gene_list,out_prefix1,op,6, "Inflamed",3)

op <- "Cl5_inflamed"
gene_list <- c("SELK", "CREM")
plot_dotplot(subset_data1,gene_list,out_prefix1,op,3, "Inflamed",5)


op <- "CL9_inflamed"
gene_list <- c("GZMK","GZMH")
plot_dotplot(subset_data1,gene_list,out_prefix1,op,3, "Inflamed",9)


