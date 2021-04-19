setwd("/projectnb2/bf528/users/hedgehog/project_4/DataCurator")

m <- read.table("counts.txt", header=FALSE)

plot(ecdf(m$V1), main = "Cumulative Distribution Plot for scRNA-seq Barcodes")
par(mfrow = c(1,1))
hist(m$V1, xlim=c(0,250000), breaks=30)
summary(m$V1)

filter <- subset(m, V1 > 396)

summary(filter$V1)

write.csv(filter, "whitlist.csv", row.names = F)
