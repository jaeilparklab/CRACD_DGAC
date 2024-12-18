rm(list = ls())
library(SCENT)
library(data.table)

count.df = as.data.frame(fread("/mnt/e/JH/Research/MSO_matrix/rna.txt", header = T))
count.df[,1] = NULL
gene.df = read.table("/mnt/e/JH/Research/MSO_matrix/gene.txt", header = T, as.is = T, sep = ",")
meta.df = read.table("/mnt/e/JH/Research/MSO_matrix/meta.txt", header = T, as.is = T, sep = ",")
meta.df$Cell_Ident = paste0(meta.df[,1], "_", meta.df$Sample, "_", meta.df$n_genes)
rownames(meta.df) = meta.df$Cell_Ident
colnames(count.df) = gene.df[,2]
rownames(count.df) = rownames(meta.df)

output.dir = "/mnt/e/JH/Research/MSO_matrix/SCENT"
dir.create(output.dir, showWarnings = F)

# Normlize to Count Per Million
cpm.df = sapply(as.data.frame(t(count.df)), function(x) {1000000 * (x/sum(x))})
rownames(cpm.df) = colnames(count.df)

# Gene symbol to entrez id
library(org.Mm.eg.db)
Mm <- org.Mm.eg.db
my.symbols <- rownames(cpm.df)
my.df = select(Mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
b.index = which(is.na(my.df[,2]))
if(length(b.index) > 0) {
  id.cpm.df = cpm.df[-b.index,]
  rownames(id.cpm.df) = my.df[-b.index,2]
}

or.file = "/home/jjang2/gene_orthologs"
or.df = as.data.frame(fread(or.file))
hm.or.df = subset(or.df, Other_tax_id == "10090" & tax_id == "9606")

id.array = rownames(id.cpm.df)
hm.array = sapply(id.array, function(x){return(subset(hm.or.df, Other_GeneID == x)[1,"GeneID"])})

b.index = which(is.na(hm.array))
if(length(b.index) > 0) {
  f.hm.array = hm.array[-b.index]
}
f.id.cpm.df = id.cpm.df[-b.index,]
rownames(f.id.cpm.df) = f.hm.array

# Run SCENT
data(net13Jun12)
log.f.id.cpm.df <- log2(f.id.cpm.df+1);
ccat.v <- CompCCAT(exp = log.f.id.cpm.df, ppiA = net13Jun12.m)
output.df = data.frame(Cell = colnames(log.f.id.cpm.df), CCAT = ccat.v)
write.table(output.df, paste0(output.dir, "/CCAT.txt"), row.names = F, col.names = T, sep = "\t", quote = F)