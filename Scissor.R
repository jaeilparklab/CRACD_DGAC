rm(list = ls())

library(dplyr)
library(patchwork)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(biomaRt)
library(Scissor)

# memory.limit(512000)
# options(future.globals.maxSize = 256000*1024^2)

output.dir = "/mnt/e/JH/Research/MSO_matrix/scissor"
dir.create(output.dir, showWarnings = F)

# 1. prepare mouse gene symbol annotated TCGA-STAD RNA-seq gene matrix
# download TCGA-STAD gene matrix including 'Primary Tumor' and 'Solid Tissue Normal'
query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  access= 'open', data.format='TXT', experimental.strategy='RNA-Seq')
GDCdownload(query)
RNA_seq <- GDCprepare(query=query, save=TRUE, save.filename='TCGA_STAD_Exp.rda')
GeneMatrix <- assay(RNA_seq)
write.csv(GeneMatrix, file=paste0(output.dir, "/TCGA_STAD_GeneMatrix.csv"), sep = ",", row.names=TRUE, col.names = TRUE)

# convert human ensembl id to mouse gene symbol
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host='https://dec2021.archive.ensembl.org/')
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host='https://dec2021.archive.ensembl.org/')
# Use follows if makes error
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
####

human_gene <- rownames(GeneMatrix)
mouse_gene <- getLDS(attributes = c("ensembl_gene_id_version"), filters = "ensembl_gene_id_version", values = human_gene, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=TRUE)
head(mouse_gene)

mouse_gene <- mouse_gene[!(mouse_gene$MGI.symbol==""),]
mouse_gene <- mouse_gene[!duplicated(mouse_gene$Gene.stable.ID.version),]
mouse_gene <- mouse_gene[!duplicated(mouse_gene$MGI.symbol),]

# make the mouse gene symbol annotated gene matrix
class(GeneMatrix)
GeneMatrix_df<- as.data.frame(GeneMatrix)
class(GeneMatrix_df)

GeneMatrix_MGI <- GeneMatrix_df
GeneMatrix_MGI$genes <- mouse_gene$MGI.symbol[match(rownames(GeneMatrix_df), mouse_gene$Gene.stable.ID.version)]
GeneMatrix_MGI <- na.omit(GeneMatrix_MGI)
rownames(GeneMatrix_MGI) <- GeneMatrix_MGI$genes
GeneMatrix_MGI = GeneMatrix_MGI[,!(colnames(GeneMatrix_MGI) %in% 'genes')]
write.csv(GeneMatrix_MGI, file=paste0(output.dir, "/GeneMatrix_MGI.csv"), sep = ",", row.names=TRUE, col.names = TRUE)

# 2. prepare phenotype list
# download TCGA-STAD clinical index data
clinical_index <- GDCquery_clinic(project = "TCGA-STAD", type = "Clinical", save.csv = TRUE)

# 2. get clinical info

query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Clinical",
  data.format = "bcr xml"
)
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
datatable(clinical, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)
write.csv(clinical, file = paste0(output.dir, "/clinical_STAD.csv"), sep = ",", row.names=TRUE, col.names = TRUE)
write.table(clinical, file = paste0(output.dir, "/clinical.STAD.txt"), sep="\t")


# Subsetting Diffuse Type

DGAC_clinical = subset(clinical, histological_type %in% "Stomach, Adenocarcinoma, Diffuse Type")
GeneMatrix_DGAC <- subset(GeneMatrix_MGI, select = substr(colnames(GeneMatrix_MGI),1,12) %in% DGAC_clinical$bcr_patient_barcode)
dim(GeneMatrix_STAD)

# Normal
GeneMatrix_Normal <- subset(GeneMatrix_MGI, select = substr(colnames(GeneMatrix_MGI),14,15) %in% '11')
write.csv(GeneMatrix_Normal, file="/mnt/e/JH/Research/MSO_matrix/scissor/GeneMatrix_Normal.csv", sep = ",", row.names=TRUE, col.names = TRUE)
Normal_phenotype <- as.data.frame(colnames(GeneMatrix_Normal))
Normal_phenotype$Disease <- 'Normal'
Normal_phenotype$Status <- 0
colnames(Normal_phenotype)<-  c("Patient", 'Disease', "Status")

# DGAC
DGAC_phenotype <- as.data.frame(colnames(GeneMatrix_DGAC))
DGAC_phenotype$Disease <- 'DGAC'
DGAC_phenotype$Status <- 1
colnames(DGAC_phenotype)<-  c("Patient", 'Disease', "Status")
Phenotype_normal_DGAC <- rbind(Normal_phenotype, DGAC_phenotype)
write.table(Phenotype_normal_DGAC, file="/mnt/e/JH/Research/MSO_matrix/scissor/Phenotype_normal_DGAC.csv", sep = ",", row.names=FALSE)

# 3. match gene matrix and phenotype list
# only include patients that also exist in the phenotype list
GeneMatrix_normal_DGAC <- subset(GeneMatrix_MGI, select = colnames(GeneMatrix_MGI) %in% Phenotype_normal_DGAC$Patient)
patient = as.data.frame(colnames(GeneMatrix_normal_DGAC))

duplicated(Phenotype_normal_DGAC$Patient)
Phenotype_normal_DGAC_1 <- Phenotype_normal_DGAC[!duplicated(Phenotype_normal_DGAC$Patient),]
Phenotype_normal_DGAC <- Phenotype_normal_DGAC_1

# make gene matrix and phenotype list the same order
GeneMatrix_normal_DGAC <- GeneMatrix_normal_DGAC[ , order(match(colnames(GeneMatrix_normal_DGAC), Phenotype_normal_DGAC$Patient))]
all(colnames(GeneMatrix_normal_DGAC) == Phenotype_normal_DGAC$Patient)

# make the sample id annotated phenotype list
rownames(Phenotype_normal_DGAC) <- Phenotype_normal_DGAC$Patient
Phenotype_normal_DGAC_input = Phenotype_normal_DGAC[,!(colnames(Phenotype_normal_DGAC) %in% c('Patient','Disease'))]
Phenotype_normal_DGAC_input
names(Phenotype_normal_DGAC_input) = rownames(Phenotype_normal_DGAC)
Phenotype_normal_DGAC_input
table(Phenotype_normal_DGAC_input)

# 4. scRNA-seq data preprocess before
library(data.table)
rna.df = as.data.frame(fread("/mnt/e/JH/Research/MSO_matrix/rna.txt", header = T))
rna.df[,1] = NULL
gene.df = read.table("/mnt/e/JH/Research/MSO_matrix/gene.txt", header = T, as.is = T, sep = ",")
meta.df = read.table("/mnt/e/JH/Research/MSO_matrix/meta.txt", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df$X
colnames(rna.df) = gene.df[,2]
rownames(rna.df) = rownames(meta.df)
t.rna.df = as.data.frame(t(rna.df))

# 5. run Scissor
phenotype <- Phenotype_normal_DGAC_input
tag <- c('Normal', 'DGAC')
GeneMatrix_normal_DGAC_input  <- as.matrix(GeneMatrix_normal_DGAC)

infos <- Scissor(bulk_dataset=GeneMatrix_normal_DGAC_input, sc_dataset=t.rna.df, phenotype=phenotype, tag = tag, alpha = 0.32, family = "binomial", Save_file = "/mnt/e/JH/Research/MSO_matrix/scissor/normal_DGAC.RData")
Scissor_select <- rep(0, ncol(t.rna.df))
names(Scissor_select) <- colnames(t.rna.df)
Scissor_select[infos$Scissor_pos] <- 1
Scissor_select[infos$Scissor_neg] <- 2
write.table(as.data.frame(Scissor_select), "/mnt/e/JH/Research/MSO_matrix/scissor/DGAC/scissor_a0.32.csv", row.names = T, col.names = T, quote = F)    # export meta.data of sc_dataset, you can check Scissor information