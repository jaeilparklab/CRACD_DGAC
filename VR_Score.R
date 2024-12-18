rm(list = ls())

phate.file = "/mnt/e/JH/Research/MSO_matrix/PHATE/PHATE.output.txt"
dist.file = "/mnt/e/JH/Research/MSO_matrix/graphtools.txt"
ccat.file = "/mnt/e/JH/Research/MSO_matrix/SCENT/CCAT.txt"

scvelo.files = Sys.glob("/mnt/e/JH/Research/MSO_matrix/scVelo/*/*.scvelo.txt")
scvelo.df = data.frame()
for(i in 1:length(scvelo.files)) {
  scvelo.file = scvelo.files[i]
  tmp.type = strsplit(basename(scvelo.file), "\\.")[[1]][1]
  tmp.scvelo.df = read.csv(scvelo.file, header = T, as.is = T, sep = ",", row.names = 1)
  tmp.scvelo.df = data.frame(Cell = paste0(rownames(tmp.scvelo.df), tmp.scvelo.df$sample_batch), Type = tmp.type, velocity_length = tmp.scvelo.df$velocity_length)
  scvelo.df =  rbind(scvelo.df, tmp.scvelo.df)
}
rownames(scvelo.df) = scvelo.df$Cell
meta.df = read.table("/mnt/e/JH/Research/MSO_matrix/meta.txt", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df[,1]

phate.df = read.delim(phate.file, header = T, as.is = T)
dist.df = read.csv(dist.file, header = T, sep = ",", row.names = 1)
ccat.df = read.delim(ccat.file, header = T, as.is = T)
rownames(ccat.df) = as.character(sapply(ccat.df[,1], function(x){strsplit(x,"_")[[1]][1]}))

inter.cells = intersect(rownames(phate.df), rownames(ccat.df))
inter.cells = intersect(inter.cells, rownames(scvelo.df))
inter.cells = intersect(inter.cells, rownames(dist.df))
output.df = data.frame(PHATE1 = phate.df[inter.cells,"PHATE1"], PHATE2 = phate.df[inter.cells,"PHATE2"], CCAT = ccat.df[inter.cells,"CCAT"], VelocityLength = scvelo.df[inter.cells,"velocity_length"], Dist = dist.df[inter.cells,"dist_scaled"], Cluster = meta.df[inter.cells, "ClusterType"])
rownames(output.df) = inter.cells
output.df$Scaled_VelocityLength = (output.df$VelocityLength)^(-1)/max((output.df$VelocityLength)^(-1))
cluster.array = unique(output.df$Cluster)
for(i in 1:length(cluster.array)) {
  tmp.cluster = cluster.array[i]
  tmp.cells = rownames(output.df)[which(output.df$Cluster == tmp.cluster)]
  output.df[tmp.cells,"Med_CCAT"] = median(output.df[tmp.cells,"CCAT"])
  output.df[tmp.cells,"Med_Scaled_VelocityLength"] = median(output.df[tmp.cells,"Scaled_VelocityLength"])
}
output.df$Valley = output.df$Med_CCAT
output.df$Ridge = output.df$Med_Scaled_VelocityLength * output.df$Dist
output.df$VR_score = 0.9*output.df$Valley + 0.1*output.df$Ridge

output.df$dataset = as.character(sapply(rownames(output.df), function(x){strsplit(x,"-")[[1]][3]}))

write.table(output.df, "/mnt/e/JH/Research/MSO_matrix/PHATE/Waddington_VR_score_PHATE_gamma.txt", row.names = T, col.names = T, sep = "\t", quote = F)