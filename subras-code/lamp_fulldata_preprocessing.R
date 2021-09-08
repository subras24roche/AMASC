#!/usr/bin/env Rscript
library(Matrix)
#library(scran)
#library(scuttle)
library(dplyr)

ar_path <- "../mountArvados/complete_dataset_6149/hg38/"
samples <- c("Sample_1a","Sample_1b","Sample_1c","Sample_2a","Sample_2b","Sample_2c","Sample_2d",
             "Sample_3a","Sample_3b","Sample_3c","Sample_3d","Sample_4a","Sample_4b","Sample_4c",
             "Sample_4d","Sample_5a","Sample_5b","Sample_5c","Sample_5d")

full_data <- readMM(paste0(ar_path,paste0(samples[1],"/outs/filtered_gene_bc_matrices/hg38/matrix.mtx")))
genes.names = read.delim(paste0(ar_path,paste0(samples[1],"/outs/filtered_gene_bc_matrices/hg38/genes.tsv")),
                         header = FALSE,
                         stringsAsFactors = FALSE)
barcode.names = read.delim(paste0(ar_path,paste0(samples[1],"/outs/filtered_gene_bc_matrices/hg38/barcodes.tsv")),
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(full_data) = barcode.names$V1
rownames(full_data) = genes.names$V2

for (i in 2:19){
  new_data <- readMM(paste0(ar_path,paste0(samples[i],"/outs/filtered_gene_bc_matrices/hg38/matrix.mtx")))
  genes.names = read.delim(paste0(ar_path,paste0(samples[i],"/outs/filtered_gene_bc_matrices/hg38/genes.tsv")),
                           header = FALSE,
                           stringsAsFactors = FALSE)
  barcode.names = read.delim(paste0(ar_path,paste0(samples[i],"/outs/filtered_gene_bc_matrices/hg38/barcodes.tsv")),
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(new_data) = barcode.names$V1
  rownames(new_data) = genes.names$V2
  full_data <- cbind(full_data, new_data)
}

mat <- as.matrix(full_data)

#mat <- logNormCounts(mat)
#sce <- calculateSumFactors(mat)
#mat <- normalizeCounts(mat,sce)

metadata <- read.csv("LC_metadata.csv")
metadata$Cell <- gsub(".*_","",metadata$Cell)
colnames(mat) <- gsub("-1","",colnames(mat))

metadata <- metadata[!duplicated(metadata$Cell),]
rownames(metadata) <- metadata$Cell
metadata <- metadata[,-1]
mat <- mat[,!duplicated(colnames(mat))]
metadata_t <- metadata[rownames(metadata) %in% colnames(mat),]

mat2 <- mat[,colnames(mat) %in% rownames(metadata_t)]
mat2 <- mat2[,order(colnames(mat2))]
metadata_t <- metadata_t[order(rownames(metadata_t)),]
annot <- select(metadata_t, starts_with("CellType"))
colnames(annot) <- "truth"

#write.csv(mat3, "Lamprecht_full_normalized_filt.csv")
#write.csv(annot, "Lamprecht_full_celltype.csv")

annot2 <- subset(annot, annot$truth != 'Cancer')
annot2 <- subset(annot2, annot2$truth != 'Fibroblast')
annot2 <- subset(annot2, annot2$truth != 'Alveolar')
annot2 <- subset(annot2, annot2$truth != 'Mast_cell')
annot2 <- subset(annot2, annot2$truth != 'EC')
annot2 <- subset(annot2, annot2$truth != 'Erythroblast')
annot2 <- subset(annot2, annot2$truth != 'Epithelial')

bad_cells <- rownames(annot)[!(rownames(annot) %in% rownames(annot2))]
mat2 <- mat2[ ,which((colnames(mat2) %in% bad_cells) == FALSE)]

annot2$truth <- replace(annot2$truth, annot2$truth == "Myeloid", 2)
annot2$truth <- replace(annot2$truth, annot2$truth == "B_cell", 4)
annot2$truth <- replace(annot2$truth, annot2$truth == "T_cell", 6)

#write.csv(annot, "Lamprecht_full_celltype_clean.csv")
write.csv(annot2, "Lamprecht_full_celltype_clean.csv")
write.csv(mat2, "Lamprecht_full_normalized_filt.csv")
