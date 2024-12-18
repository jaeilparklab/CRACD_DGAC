rm(list = ls())
library(Seurat)
library(data.table)
library(CytoTRACE)

rna.df = as.data.frame(fread("/mnt/e/JH/Research/MSO_matrix/rna.txt"))
rna.df[,1] = NULL
rna.df = rna.df[2:nrow(rna.df),]
meta.df = read.csv("/mnt/e/JH/Research/MSO_matrix/meta.txt", sep = ",", header = T, as.is = T)
rownames(meta.df) = meta.df[,1]
gene.df = read.csv("/mnt/e/JH/Research/MSO_matrix/gene.txt", as.is = T)
colnames(rna.df) = gene.df[,2]
rownames(rna.df) = meta.df[,1]
t.rna.df = as.data.frame(t(rna.df))

output.dir = "/mnt/e/JH/Research/MSO_matrix/CytoTRACE"
dir.create(output.dir, showWarnings = F)

rownames(meta.df) = meta.df[,1]
colnames(t.rna.df) = rownames(meta.df)
input.mat = as.matrix(t.rna.df)  

results <- CytoTRACE(input.mat)
saveRDS(results, file = paste0(output.dir, "/CytoTRACE_output.RDS"))
results <- readRDS(paste0(output.dir, "/CytoTRACE_output.RDS"))

cyto.array = results$CytoTRACE
rank.array = results$CytoTRACErank
result.df = data.frame(Cell = names(cyto.array), Cyto = cyto.array, Rank = rank.array)
write.table(result.df, paste0(output.dir, "/CytoTRACE_output.txt"), sep = "\t", quote = F, col.names = T, row.names = F)