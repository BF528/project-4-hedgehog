library(dplyr)
library(Seurat)
install.packages("patchwork")
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "//projectnb/bf528/project_4_scrnaseq/GSM2230760_seurat.rda")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc