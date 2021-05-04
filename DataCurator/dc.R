gc() #for memory management
rm(list=ls()) #Clear environment

install.packages("BiocManager")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("gplots")
library("tidyverse")
library("samtools")

#Metadata
metapath<-"/projectnb/bf528/users/hedgehog/project_4/DataCurator/p4_sample_metadata.txt"
metadata<-read_csv(metapath)
spec(metadata)#Metadata fields and type

#Sample codes are under column "Run"
SRR51years <- subset(metadata, metadata$AGE==51 & metadata$sex=="female", select=c("Run","AGE","sex"))



#Reading barcodes and their frequecy counts
SRR3879604_counts<- read.csv('/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/results/count_SRR3879604.csv', col.names=c('barcode', 'count'), sep='')
SRR3879605_counts<- read.csv('/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/results/count_SRR3879605.csv', col.names=c('barcode', 'count'), sep='')
SRR3879606_counts<- read.csv('/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/results/count_SRR3879606.csv', col.names=c('barcode', 'count'), sep='')

#Cleaning data => Removing null values (nas) from our data 
SRR3879604_counts <- SRR3879604_counts %>% filter(!is.na(count))
SRR3879605_counts <- SRR3879605_counts %>% filter(!is.na(count))
SRR3879606_counts <- SRR3879606_counts %>% filter(!is.na(count))

#Computing the mean for each sample group
SRR3879604_mean <- mean(SRR3879604_counts$count)
SRR3879605_mean <- mean(SRR3879605_counts$count)
SRR3879606_mean <- mean(SRR3879606_counts$count)

#sorting
SRR3879604_counts <-SRR3879604_counts %>% arrange(desc(count))
SRR3879605_counts<- SRR3879605_counts %>% arrange(desc(count))
SRR3879606_counts <- SRR3879606_counts %>% arrange(desc(count))

#filtering out those values greater than the mean
SRR3879604_whitelisted <- SRR3879604_counts %>% filter(SRR3879604_counts$count>SRR3879604_mean)
SRR3879605_whitelisted <- SRR3879605_counts %>% filter(SRR3879605_counts$count>SRR3879605_mean)
SRR3879606_whitelisted <- SRR3879606_counts %>% filter(SRR3879606_counts$count>SRR3879606_mean)

#Saving data version for saving- it avoids general lost of data when somthing goes wrong with this dataframe
SRR3879604_whitelisted<-data.frame(SRR3879604_whitelisted$barcode)
SRR3879605_whitelisted<-data.frame(SRR3879605_whitelisted$barcode)
SRR3879606_whitelisted<-data.frame(SRR3879606_whitelisted$barcode)

write_csv(SRR3879604_whitelisted, '/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/salmon/SRR3879604_whitelist.txt', col_names = FALSE)
write_csv(SRR3879605_whitelisted, '/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/salmon/SRR3879605_whitelist.txt', col_names = FALSE)
write_csv(SRR3879606_whitelisted, '/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/salmon/SRR3879606_whitelist.txt', col_names = FALSE)





png('/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/plots/commulative_plot_SRR3879604.png',  width = 900, height = 600 )
plot(ecdf(SRR3879604_counts$count), cex=0)
dev.off()

png('/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/plots/commulative_plot_SRR3879605.png',  width = 900, height = 600 )
plot(ecdf(SRR3879605_counts$count), cex=0)
dev.off()


png('/projectnb/bf528/users/hedgehog/project_4/DataCurator/a2/plots/commulative_plot_SRR3879606.png',  width = 900, height = 600 )
plot(ecdf(SRR3879606_counts$count), cex=0)
dev.off()



