library(Seurat)
library(rafalib)
library(igraph)
library(scales)
library(RColorBrewer)
library(MASS)
library(spatstat)
library(magick)
library(ggplot2)
library(raster)
library(pheatmap)
library(circlize)
library(rafalib)
library(dendextend)
library(WGCNA)
library(RANN)
library(plyr)
library(dplyr)


outdir="Results"
dir.create(file.path(outdir), showWarnings = FALSE)
#setwd(outdir)

output_prefix <- paste(outdir, "/Kokkinou_pIBD_etal_2023",sep = "")


m <- as.matrix(read.csv2("Data/50_modules_genes.csv", row.names = 1))[,1]
a <- readRDS("Data/50_pbmc.sub2.module.rds")

a@meta.data$Cluster <- a@meta.data$clusters_louvain_subset1
head(a@meta.data)
jakob_table_df <- as.data.frame(read.csv("Data/Complete_metadata.csv",sep = 
                                           ",",header = T))



row.names(jakob_table_df) <- jakob_table_df$X
head(jakob_table_df)
colnames(jakob_table_df)
#jakob_table_df1 <- jakob_table_df[,c(16,17,18,19,20)]
jakob_table_df1 <- jakob_table_df[,c(16,17,18,19)]

colnames(jakob_table_df)
head(jakob_table_df1)
colnames(jakob_table_df1)
jakob_table_df1$Hist.score <- jakob_table_df1$neighbor_smooth_CHASTea
dim(jakob_table_df1)


a <- AddMetaData(object=a, metadata=jakob_table_df1)
head(a@meta.data)


#rowmax <- apply(a@assays[[ opt$assay ]]@data[names(m),],1,function(x) sort(x[x>0],T)[1] )
#module_means <- rowsum(as.matrix(a@assays[[ opt$assay ]]@data[names(m),]/(rowmax+1) ), m)

rowmax <- apply(a@assays[[ "RNA" ]]@data[names(m),],1,function(x) sort(x[x>0],T)[1] )
module_means <- rowsum(as.matrix(a@assays[[ "RNA" ]]@data[names(m),]/(rowmax+1) ), m)

module_means <- t(module_means)/c(table(m))
a@assays[["modules"]] <- Seurat::CreateAssayObject(data = t(as.matrix(module_means)), min.cells = 0,min.features = 0)

n <-  nrow(a@assays$modules@data)
# 
# # Compute the kNN graph for the top 5 gene modules
cors <- (1 - WGCNA::cor( t(a@assays$modules@data) , nThreads = 0)) / 2
NN <- RANN::nn2(cors, k = 5, eps = 0)
NN <- data.frame( rep(NN$nn.idx[,1], ncol(NN$nn.idx)-1), c(NN$nn.idx[,-1]), c(NN$nn.dists[,-1]) )
colnames(NN) <- c("from","to","weight")
NN <- NN[ NN$weight < 0.65 , ]
NN$scaled_weight <- ((1-NN$weight) - min(1-NN$weight) ) / (max(1-NN$weight) - min(1-NN$weight) )
dim(NN)
# 
# #dend <- hclust( dist( a@assays$modules@data ) ,method = "complete" )
# # dend2 <- hclust( as.dist( 1-cor( t(a@assays$modules@data) )) ,method = "complete" )
# 
# # Defining the relevance of each module by calculating how much variability is explained by each metadata factor
head(a@meta.data)
colnames(a@meta.data)[12] <- "Donor"
unique(a@meta.data$Donor)
unique(a@meta.data$donor1)

meta_data1 <- as.data.frame(a@meta.data)
meta_data1$barcode <- row.names(meta_data1)
Donor_sum <- ddply(meta_data1, .(Donor), summarise, Cell_Count=length(barcode))
#write.table(Donor_sum,"Donor_to_cell_count.txt", row.names=F, sep="\t", quote=F)

Donor_sum <- ddply(meta_data1, .(Donor,Cluster), summarise, Cell_Count=length(barcode))
#write.table(Donor_sum,"Donor+Cluster_to_cell_count.txt", row.names=F, sep="\t", quote=F)



meta_data <- c("RNA_snn_res.2.4","Tissue","Celltype","Donor","Plate")
meta_data <- c("Cluster","Donor")

colnames(a@meta.data)
a@meta.data$Tissue
# #Defining the relative importance of each module for explaining the clustering
for(x in meta_data){
  message(x)
  dev <- sapply(1:n, x=x, function(i,x){mod <- glm( a@assays$modules@data[i,] ~ a@meta.data[,x] ); c(mod$deviance, mod$null.deviance) } )
  assign( paste0("deltadeviance","_",x) , (dev[2,]-dev[1,])/dev[2,])
}
# 
#modules_interest <- sort(unique(c(32,23,28,2,3,13,60,4,26,70,7,6,50,17,65,41,11,23,28,34,54,25)))
modules_interest <- sort(unique(c(1:50)))
# 
# pdf("~/Box/Projects/J_Mjosberg/20200110/gene_relevance_plot.pdf",width = 12,height = 6,useDingbats = F)
rafalib::mypar(length(meta_data),1,mar=c(2,15,1,1),mgp = c(2,0.5,0) )
for(x in meta_data){
  i <- get(paste0("deltadeviance","_",x))
  barplot(i,las=2, pch= 16,col="grey80",
          xaxs="i",names.arg = 1:n,border = NA,ylim=c(0,1.1),cex.main=.7)
  # text( (modules_interest*1.2-0.5), i[ modules_interest ], labels = modules_interest,pos = 3 )
  text( (modules_interest*1.2-0.5), i[ modules_interest ], labels = "*",pos = 3 )
  mtext(x,side = 2,las=1,cex=.8,adj = 1,line = 2)
  lines( c(0,n*1.2), c(0,0) )
  lines( c(0,n*1.2), 0.1*c(1,1),lty=2,lwd=.5 )
  points( (1:n)*1.2-0.5 ,i, las=1, cex=.8, pch= 21,bg="grey70")
}
# dev.off()
# 
# 
# # Plot the circos plot

#pdf("Gene_Module_circular_plot_2022_04_22.pdf",width = 7,height = 7,useDingbats = F)
# pdf("Gene_Module_circular_plot_2022_08_28.pdf",width = 7,height = 7,useDingbats = F)
# 
# rafalib::mypar(mar=c(0,0,0,0))
# circos.par(start.degree = 90)
# 
# circos.initialize("foo",xlim = c(0,n+10 ))
# pos = circlize(x = 0:( n+7 ), 0:( n+7 ), sector.index = "foo")
# 
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   circos.text((1:n)-0.25, rep(0, n), ifelse(1:n %in% modules_interest,"*",""), col = "black",
#               facing = "outside", adj = c(0, 0.5),cex = .9)
# }, bg.border = NA, track.height = 0.05,track.margin=c(0.00,0.00),cell.padding=c(0,0,0,0))
# 
# 
# 
# #Add annotation to the modules
# deltadeviance_Cluster <- (deltadeviance_Cluster - 0) / (max(deltadeviance_Cluster) - 0)
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   circos.rect(1:n-.1, rep(0, n), 1:n-.8, deltadeviance_Cluster , col = colorRampPalette(c("grey95","black"))(100)[ round(deltadeviance_Cluster*98)+1 ]
#               , border = NA)
# }, bg.border = NA, track.height=0.05,track.margin=c(0,.03),cell.padding=c(0,0,0,0))
# circos.yaxis(side="left",at = c(0.5),labels = expression( paste("Cluster ", Delta, " deviance" ) ), labels.cex = .7,col = "white")
# 
# deltadeviance_Tissue <- (deltadeviance_Tissue - 0) / (max(deltadeviance_Tissue) - 0)
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   circos.rect(1:n-.1, rep(0, n), 1:n-.8, deltadeviance_Tissue , col = colorRampPalette(c("grey95","black"))(100)[ round(deltadeviance_Tissue*98)+1 ]
#               , border = NA)
# }, bg.border = NA, track.height=0.05,track.margin=c(0,.01),cell.padding=c(0,0,0,0))
# circos.yaxis(side="left",at = c(0.5),labels = expression( paste("Tissue ", Delta, " deviance" ) ), labels.cex = .7,col = "white")
# 
# deltadeviance_Donor <- (deltadeviance_Donor - 0) / (max(deltadeviance_Donor) - 0)
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   circos.rect(1:n-.1, rep(0, n), 1:n-.8, deltadeviance_Donor , col = colorRampPalette(c("grey95","black"))(100)[ round(deltadeviance_Donor*98)+1 ]
#               , border = NA)
# }, bg.border = NA, track.height=0.05,track.margin=c(0,.01),cell.padding=c(0,0,0,0))
# circos.yaxis(side="left",at = c(0.5),labels = expression( paste("Donor ", Delta, " deviance") ), labels.cex = .7,col = "white")
# 
# deltadeviance_Hist.score <- (deltadeviance_Hist.score - 0) / (max(deltadeviance_Hist.score) - 0)
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   circos.rect(1:n-.1, rep(0, n), 1:n-.8, deltadeviance_Hist.score , col = colorRampPalette(c("grey95","black"))(100)[ round(deltadeviance_Hist.score*98)+1 ]
#               , border = NA)
# }, bg.border = NA, track.height=0.05,track.margin=c(0,.01),cell.padding=c(0,0,0,0))
# circos.yaxis(side="left",at = c(0.5),labels = expression( paste("Hist.score ", Delta, " deviance") ), labels.cex = .7,col = "white")
# 
# 
# #Adds the gene module ID
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   circos.rect(1:n-.1, rep(0, n), 1:n-.8,
#               (table(m)/max(table(m))) , col = scales::hue_pal()(n), border = NA)
# }, bg.border = NA, track.height=0.05,track.margin=c(0,.01),cell.padding=c(0,0,0,0))
# circos.yaxis(side="left",at = c(0,1),labels = c(0,max(table(m))), labels.cex = .7)
# circos.yaxis(side="left",at = c(0.5),labels = "no. of genes        ", labels.cex = .7)
# 
# circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#   circos.text(1:n-0.5, rep(0, n), 1:n, col = scales::hue_pal()(n),
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = .7)
# }, bg.border = NA, track.height = 0.05,track.margin=c(0,0.02),cell.padding=c(0,0,0,0))
# circos.yaxis(side="left",at = c(0.5),labels = "module ID", labels.cex = .6,col = "white")
# 
# 
# #Adds the lines linking modules that are related to others
# for(i in 1:nrow(NN)){
#   circos.link("foo", NN[i,1]-0.5, 'foo', NN[i,2]-.5,h.ratio = 1,
#               lwd = NN[i,4]*2+.2,col = paste0( scales::hue_pal()(n)[min(NN[i,1:2])],"60")  )}
# 
# #Adds the number of genes in each module
# circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,track.margin=c(0.00,0.04),cell.padding=c(0,0,0,0))
# circos.yaxis(side="left",at = c(0.5),labels = "module kNN", labels.cex = .7,col = "white")
# 
# circos.clear()
# dev.off()



###############################

#pdf("Gene_Module_circular_plot_no_Tissue_no_Hist.Score_2022_04_23.pdf",width = 7,height = 7,useDingbats = F)
plot_file1 <- paste(output_prefix,"Figure_3B.jpeg",sep = "_")
jpeg(filename = plot_file1, width = 7, height = 7, units = "in",res=400)


rafalib::mypar(mar=c(0,0,0,0))
circos.par(start.degree = 90)

circos.initialize("foo",xlim = c(0,n+10 ))
pos = circlize(x = 0:( n+7 ), 0:( n+7 ), sector.index = "foo")

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text((1:n)-0.25, rep(0, n), ifelse(1:n %in% modules_interest,"*",""), col = "black",
              facing = "outside", adj = c(0, 0.5),cex = .9)
}, bg.border = NA, track.height = 0.05,track.margin=c(0.00,0.00),cell.padding=c(0,0,0,0))



#Add annotation to the modules
deltadeviance_Cluster <- (deltadeviance_Cluster - 0) / (max(deltadeviance_Cluster) - 0)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(1:n-.1, rep(0, n), 1:n-.8, deltadeviance_Cluster , col = colorRampPalette(c("grey95","black"))(100)[ round(deltadeviance_Cluster*98)+1 ]
              , border = NA)
}, bg.border = NA, track.height=0.05,track.margin=c(0,.03),cell.padding=c(0,0,0,0))
circos.yaxis(side="left",at = c(0.5),labels = expression( paste("Cluster ", Delta, " deviance" ) ), labels.cex = .7,col = "white")

deltadeviance_Donor <- (deltadeviance_Donor - 0) / (max(deltadeviance_Donor) - 0)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(1:n-.1, rep(0, n), 1:n-.8, deltadeviance_Donor , col = colorRampPalette(c("grey95","black"))(100)[ round(deltadeviance_Donor*98)+1 ]
              , border = NA)
}, bg.border = NA, track.height=0.05,track.margin=c(0,.01),cell.padding=c(0,0,0,0))
circos.yaxis(side="left",at = c(0.5),labels = expression( paste("Donor ", Delta, " deviance") ), labels.cex = .7,col = "white")


#Adds the gene module ID
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(1:n-.1, rep(0, n), 1:n-.8,
              (table(m)/max(table(m))) , col = scales::hue_pal()(n), border = NA)
}, bg.border = NA, track.height=0.05,track.margin=c(0,.01),cell.padding=c(0,0,0,0))
circos.yaxis(side="left",at = c(0,1),labels = c(0,max(table(m))), labels.cex = .7)
circos.yaxis(side="left",at = c(0.5),labels = "no. of genes        ", labels.cex = .7)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:n-0.5, rep(0, n), 1:n, col = scales::hue_pal()(n),
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = .7)
}, bg.border = NA, track.height = 0.05,track.margin=c(0,0.02),cell.padding=c(0,0,0,0))
circos.yaxis(side="left",at = c(0.5),labels = "module ID", labels.cex = .6,col = "white")


#Adds the lines linking modules that are related to others
for(i in 1:nrow(NN)){
  circos.link("foo", NN[i,1]-0.5, 'foo', NN[i,2]-.5,h.ratio = 1,
              lwd = NN[i,4]*2+.2,col = paste0( scales::hue_pal()(n)[min(NN[i,1:2])],"60")  )}

#Adds the number of genes in each module
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,track.margin=c(0.00,0.04),cell.padding=c(0,0,0,0))
circos.yaxis(side="left",at = c(0.5),labels = "module kNN", labels.cex = .7,col = "white")

circos.clear()
dev.off()


######### Suppl. figure S3 A-G
pbmc.sub2.module <- readRDS("Data/50_pbmc.sub2.module.rds")

feats <- rownames(pbmc.sub2.module@assays$modules@data)
feats <- c("1","11","43","8","19","26","7","20","18","22","30")
for(i in feats){ 
  #pdf( paste(outdir2,"/",k,"_module",i,"_genes_umap_blue1.pdf",sep = ""),width = 6,height = 6, useDingbats=FALSE)
  plot_file1 <- paste(output_prefix,"_Figure_S3_A_to_G_module_",i,".jpeg",sep = "")
  jpeg(filename = plot_file1, width = 4, height = 4, units = "in",res=400)

  plot_feat(pbmc.sub2.module, feat = i , assay="modules", 
            col = c("grey95","grey70","#99B1FF","#003CFF"),
            add_legend=T,cex = 0.2)
  #  print(plot_feat(subset_data_sub,feat = g1, add_legend=T,cex = 0.4))

  dev.off()
}

