#setwd("/projectnb/bf528/users/hedgehog/project_4/Biologist")
library(dplyr)

# Read in the data
marker_genes <- read.csv("/projectnb/bf528/users/hedgehog/project_4/Analyst/DE_genes_all.csv")

# gL function
# takes the marker genes df, returns list of genes in each cluster
gL <- function(df){
  geneList <- list()
  clust <- levels(factor(df$cluster))
  for(i in clust){
    genes <- df %>% 
      filter(cluster==i) %>% 
      pull(gene)
    genes <- list(genes)
    names(genes) <- paste0("C_",i) 
    geneList <- c(geneList, genes)
  }
  return(geneList)
}


# Finding threshold values ------------------------------------------------

# adjusted p-value
l <- 250
vals <- seq(max(marker_genes$p_val_adj), min(marker_genes$p_val_adj), length.out = l)
lengths <- rep(0,l) # vector of lengths
for(i in 1:length(vals)){
  test <- which(marker_genes$p_val_adj<=vals[i])
  lengths[i] <- length(test)
}
thresh_pval <- vals[which.min(lengths)-1]
#jpeg(file = "suppFig1.jpeg")
plot(vals, lengths,
     xlab = "Adjusted P-Values",
     ylab = "Number of observations",
     main = "Distribution of observation number vs p-value")
text(y=lengths[which.min(lengths)-1]-200, x=vals[which.min(lengths)]-0.005,
     label="p-value\n=0.004",
     cex=0.7,
     pos=4)
abline(v=thresh_pval, col="red", lty=3)
#dev.off()

# log2 fold change
l <- 1000
vals <- seq(min(marker_genes$avg_log2FC), max(marker_genes$avg_log2FC), length.out = l)
lengths <- rep(0,l) # vector of lengths
for(i in 1:length(vals)){
  test <- which(marker_genes$avg_log2FC>vals[i])
  lengths[i] <- length(test)
}

obv <- dim(marker_genes)[1]*0.5 # want to retain at least 50% of data
thresh_logfc <- vals[min(which(lengths<=obv))]

#jpeg(file = "suppFig2.jpeg")
plot(vals, lengths,
     xlab = "Log2 Fold Change",
     ylab = "Number of observations",
     main = "Distribution of observation number vs fold change")
abline(h=lengths[min(which(lengths<=obv))],
       col="red", lty=3)
abline(v=vals[min(which(lengths<=obv))],
       col="red", lty=3)
text(y=lengths[min(which(lengths<=obv))]+75, x=vals[min(which(lengths<=obv))],
     label="log2(FC)=0.733",
     cex=0.7,
     pos=4)
#dev.off()



# Getting the gene lists --------------------------------------------------

# No filters applied
geneList <- gL(marker_genes)

# With filters (threshold values)
marker_genes_filt <- marker_genes %>%
  filter(p_val_adj<thresh_pval) %>% 
  filter(avg_log2FC>thresh_logfc)
geneList_filt <- gL(marker_genes_filt)

# text files --------------------------------------------------------------
for(i in 1:length(geneList)){
  filename <- names(geneList[i])
  filename <- paste0('no_filt/', filename, ".txt")
  fileConn <- file(filename)
  writeLines(geneList[i][[1]], fileConn)
  close(fileConn)
}

for(i in 1:length(geneList_filt)){
  filename <- names(geneList_filt[i])
  filename <- paste0('filt/', filename, ".txt")
  fileConn <- file(filename)
  writeLines(geneList_filt[i][[1]], fileConn)
  close(fileConn)
}

