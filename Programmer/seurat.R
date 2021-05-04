# Taylor Falk
# tfalk@bu.edu
# BF528 - Project 

# bioconductor installs:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("fishpond")

library(biomaRt)
library(dplyr)
library(fishpond)
library(patchwork)
library(Seurat)
library(tximport)

data_1 <- paste0("/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant/",
                 "alevin/quants_mat.gz")
data_1 <- paste0("/projectnb2/bf528/users/hedgehog/project_4/DataCurator/", 
                 "alevin_output_all_1k_v2.2/alevin/quants_mat.gz")
file.exists(data_1)

txi <- tximport(data_1, type="alevin")

# ensembl names
splitnames <- sapply(strsplit(as.character(txi$counts@Dimnames[[1]]), "\\."), "[[", 1)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
newNames <- getBM(filters= "ensembl_gene_id", 
                  attributes= c("ensembl_gene_id", "hgnc_symbol"),
                  values=splitnames, 
                  mart= mart)
newNames <- newNames[-which(duplicated(newNames$ensembl_gene_id)),]
names = c()
for (i in c(1:length(splitnames))) {
  row = which(newNames$ensembl_gene_id == splitnames[i])
  if (newNames$hgnc_symbol[row] != "") {
    names <- append(names, newNames$hgnc_symbol[row])
  }
  else {
    names <- append(names, splitnames[i])
  }
}

sprintf("Pre-filtering: %i gene %i cell", 
        length(txi$counts@Dimnames[[1]]), length(txi$counts@Dimnames[[2]]))
txi$counts@Dimnames[[1]] <- names
panc <- CreateSeuratObject(counts = txi$counts , min.cells = 3, 
                           min.features = 200, project = "panc")

sprintf("There are %i cells and %i genes after loading.", 
        panc@assays$RNA@counts@Dim[2], panc@assays$RNA@counts@Dim[1])


VlnPlot(panc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

plot1 <- FeatureScatter(panc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

panc <- subset(panc, subset = nFeature_RNA > 300 & nFeature_RNA < 2500)
sprintf("Post-post-filtering: %i gene %i cell", 
        panc@assays$RNA@counts@Dim[1], panc@assays$RNA@counts@Dim[2])
panc <- NormalizeData(panc)

# set threshold of variable features
panc <- FindVariableFeatures(panc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(panc), 10)

# plot variable features
plot1 <- VariableFeaturePlot(panc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, labels = top10)
plot2

# return only 2,000 features
panc <- subset(panc, features = VariableFeatures(panc))
sprintf("Dims after variable filter: %i gene %i cell", 
        panc@assays$RNA@counts@Dim[1], panc@assays$RNA@counts@Dim[2])

# scaling the data
all.genes <- rownames(panc)
panc <- ScaleData(panc, features = all.genes)

# linear dim reduction
panc <- RunPCA(panc, features = VariableFeatures(object = panc))
VizDimLoadings(panc, dims = 1:2, reduction = "pca")
DimPlot(panc, reduction = "pca")
DimHeatmap(panc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(panc, dims = 1:3, cells = 500, balanced = TRUE)

# determine dimensionality
panc <- JackStraw(panc, num.replicate = 100)
panc <- ScoreJackStraw(panc, dims = 1:20)
JackStrawPlot(panc, dims = 1:20) +
  ElbowPlot(panc)
# looks to be about 8-9 dims

# clustering
panc <- FindNeighbors(panc, dims = 1:9)
panc <- FindClusters(panc, resolution = 1)

panc <- RunUMAP(panc, dims = 1:9)
panc <- RunTSNE(panc, dims = 1:9)
DimPlot(panc, reduction = "tsne") +
  DimPlot(panc, reduction = "umap")

# diff exp clusters
cluster1.markers <- FindMarkers(panc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(panc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

panc.markers <- FindAllMarkers(panc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topN <- panc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
VlnPlot(panc, features = topN[1:3,7]$gene)
VlnPlot(panc, features = topN[1:2,7]$gene, slot = "counts", log = TRUE)
FeaturePlot(panc, features = topN[1:9,7]$gene)
top10 <- panc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(panc, features = top10$gene) + NoLegend()

sprintf("Dims saved: %i gene %i cell", 
        panc@assays$RNA@counts@Dim[1], panc@assays$RNA@counts@Dim[2])
saveRDS(panc, file = "/projectnb2/bf528/users/hedgehog/project_4/Programmer/panc.rda")
