rm(list = ls())

library(data.table)
library(phateR)
library(ggplot2)

# Prepare input data
count.df = as.data.frame(fread("/mnt/e/JH/Research/MSO_matrix/rna.txt", header = T))
count.df[,1] = NULL
gene.df = read.table("/mnt/e/JH/Research/MSO_matrix/gene.txt", header = T, as.is = T, sep = ",")
meta.df = read.table("/mnt/e/JH/Research/MSO_matrix/meta.txt", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df$X
colnames(count.df) = gene.df[,2]
rownames(count.df) = rownames(meta.df)

output.dir = "/mnt/e/JH/Research/MSO_matrix/PHATE"
dir.create(output.dir, showWarnings = F)
input.df = count.df # Row as Cell, Column as Gene

# Keep genes expressed in at least 10 cells
keep_cols <- colSums(input.df > 0) > 10
input.df <- input.df[,keep_cols]

# Keep cells with at least 1000 UMIs
keep_rows <- rowSums(input.df) > 1000
input.df <- input.df[keep_rows,]

input.df <- library.size.normalize(input.df)
input.df <- sqrt(input.df)

input_PCA <- as.data.frame(prcomp(input.df)$x)
write.table(input_PCA, paste0(output.dir, "/PCA.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

# Run PHATE
input_PHATE <- phate(input_PCA)
output.df = as.data.frame(input_PHATE[[1]])
write.table(output.df, paste0(output.dir, "/PHATE.output.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

# Cell type plot (Merge)
pdf(paste0(output.dir, "/PHATE.Cluster.pdf"), width=7.5, height=5)
ggplot(output.df) +
  geom_point(aes(PHATE1, PHATE2, colour = leiden), size = 0.1) +
  labs(colour="Type") +
  xlim(c(min(output.df$PHATE1), max(output.df$PHATE1))) +
  ylim(c(min(output.df$PHATE2), max(output.df$PHATE2))) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
dev.off()