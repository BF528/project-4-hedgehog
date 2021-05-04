suppressPackageStartupMessages({
  install.packages("patchwork")
  install.packages("Seurat")
  install.packages("R.filesets")
  install.packages("ggsci")
  install.packages("eulerr")
  install.packages('UpSetR')
  install.packages("shiny")
})


suppressPackageStartupMessages({
  library(ggplot2)
  library(Rtsne)
  library(dplyr)
  library(ggrepel)
  library(pheatmap)
  library(R.filesets)  
  library(dplyr)
  library(magrittr)
  library(glue)
  library(patchwork)
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(cowplot)
  library(RColorBrewer)
  library(ggsci)
  library(eulerr)
  library(UpSetR)
  library(Rtsne)
  library(shiny)
})

#1 Load the  dataset
cellsnew <- loadRDS("/projectnb2/bf528/users/hedgehog/project_4/Analyst/cellsnew.rds")
saveRDS(cellsnew, "/projectnb2/bf528/users/hedgehog/project_4/Analyst/cellsnew.rds")

#clustering, tried several min.pct and logfc.thershold options 15% to 90% & 25% to 90% respectively. Below we only list 33% 
markersDE <- FindAllMarkers(object = cellsnew, min.pct = .33,
                            only.pos = TRUE,
                            logfc.threshold = .33, return.thresh = 0.05) 

#Save initial file as CSV
write.csv(markersDE,'/projectnb/bf528/users/hedgehog/project_4/Analyst/allgenes2_pan.csv')


##########################
#2 Label the clusters with cell types
#Create new dataframe for cell type assignment
markDE<-markersDE

#Add labels to clusters
DE_labels<-markDE%>%mutate(celltype = 
                             case_when(markDE$cluster == "0" ~ "Mixed",
                                       markDE$cluster == "1" ~ "Alpha1", 
                                       markDE$cluster == "2" ~ "Acinar", 
                                       markDE$cluster == "3" ~ "Alpha2", 
                                       markDE$cluster == "4" ~ "Gamma", 
                                       markDE$cluster == "5" ~ "Delta-Beta", 
                                       markDE$cluster == "6" ~ "Ductal-Epsilon", 
                                       markDE$cluster == "7" ~ "Delta-Gamma", 
                                       markDE$cluster== "8" ~ "Activated Stellate", 
                                       markDE$cluster == "9" ~ "Beta", 
                                       markDE$cluster == "10" ~ "Ductal-Acinar", 
                                       markDE$cluster == "11" ~ "Schwann", 
                                       markDE$cluster == "12" ~ "Endothelial", 
                                       markDE$cluster == "13" ~ "Macrophage",
                                       markDE$cluster == "14" ~ "Quiescent Stellate-Epsilon"))


####################
#6 Identified marker genes with celltype
write.csv(DE_labels,'/projectnb/bf528/users/hedgehog/project_4/Analyst/DE_genes_all.csv')


####################
#3 Visualization of clustered data using UMAP
cellp <- readRDS("/projectnb/bf528/users/hedgehog/project_4/Programmer/panc.rda")

##renaming levels to label celltype
names(newid)<-levels(cellp)
cellp<-RenameIdents(cellp,newid)

#run umap
pan_umap <- RunUMAP(cellp, dims = 1:10, LabelClusters = newid)
#plot Seurat object
pan_plot<-DimPlot(pan_umap, reduction = 'umap', group.by = "pan_umap@active.ident")
#custom labels
pan_plot + labs(title = "Pancreas Gene Clustering", repel = TRUE, label_size = 2)+DarkTheme()



#########################
#4 Visualize top 5, and top 3marker genes per cluster
sig_genes<-DE_labels
top5 <- sig_genes %>% group_by(celltype)%>% top_n(5, wt = avg_log2FC)
top3 <- sig_genes %>% group_by(celltype)%>% top_n(3, wt = avg_log2FC)

# Save top 5 
write.csv(top5,'/projectnb/bf528/users/hedgehog/project_4/Analyst/top5_genes.csv')

#Read-in Seurat Object
cellx <- readRDS("/projectnb/bf528/users/hedgehog/project_4/Analyst/cellsnew.rds")

#rename clusters levels
newid <- c("Mixed", "Alpha1", "Acinar", "Alpha2", "Gamma", "Delta-Beta", "Ductal-Epsilon", "Delta-Gamma", 
           "Activated stellate", "Beta", "Ductal-Acinar", "Schwann", "Endothelial", "Macrophage", "Quiescent Stellate-Epsilon")
names(newid)<-levels(cellx)
cellx<-RenameIdents(cellx,newid)

#Create Heatmap, using TOP 5 
DoHeatmap(cellx, features = top5$gene, size = 3.5)+NoLegend()

#create violin plots with the most DE genes, as expressed on heatmap and paper - 
VlnPlot(cellsnew, features = c("SST", "INS", "IAPP", "REG3A", "GCG", "KRT19", "KRT8", "PPY","COL1A2", "VWF", "CD68", "MT-ND1"), slot = "counts", log = TRUE)

#create Featureplots based on top DE genes, as expressed on heatmap and paper
FeaturePlot(cellx, features = c("SST", "INS", "IAPP", "REG3A", "GCG", "KRT19", "KRT8", "PPY","COL1A2", "VWF", "CD68","MT-ND1"))


