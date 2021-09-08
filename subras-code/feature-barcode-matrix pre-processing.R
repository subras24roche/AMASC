library(Matrix)
library(cellranger)
library(scran)
library(dplyr)

matrix_dir = "Lamprecht sample1a/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
genes.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")

mat <- readMM(file = matrix.path)
genes.names = read.delim(genes.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

colnames(mat) = barcode.names$V1
rownames(mat) = genes.names$V2

mat <- as.matrix(mat)
write.csv(mat, "Lamprecht_1a_counts.csv")

sce <- calculateSumFactors(mat)
mat2 <- normalizeCounts(mat,sce)
write.csv(mat2, "Lamprecht_1a_normalized.csv")

metadata <- read.csv("LC_metadata.csv")
metadata$Cell <- gsub(".*_","",metadata$Cell)
colnames(mat2) <- gsub("-1","",colnames(mat2))

metadata <- metadata[!duplicated(metadata$Cell),]
rownames(metadata) <- metadata$Cell
metadata <- metadata[,-1]
metadata_t <- metadata[rownames(metadata) %in% colnames(mat2),]

mat3 <- mat2[,colnames(mat2) %in% rownames(metadata_t)]
mat3 <- mat3[,order(colnames(mat3))]
metadata_t <- metadata_t[order(rownames(metadata_t)),]
annot <- select(metadata_t, starts_with("CellType"))
colnames(annot) <- "truth"

#write.csv(mat3, "Lamprecht_1a_normalized_filt.csv")
write.csv(annot, "Lamprecht_1a_celltype.csv")

annot2 <- subset(annot, annot$truth != 'Cancer')
annot2 <- subset(annot2, annot2$truth != 'Fibroblast')
annot2 <- subset(annot2, annot2$truth != 'Alveolar')
annot2 <- subset(annot2, annot2$truth != 'Mast_cell')

bad_cells <- rownames(annot)[!(rownames(annot) %in% rownames(annot2))]
mat3 <- mat3[ ,which((colnames(mat3) %in% bad_cells) == FALSE)]

annot2$truth <- replace(annot2$truth, annot2$truth == "Myeloid", 2)
annot2$truth <- replace(annot2$truth, annot2$truth == "B_cell", 4)
annot2$truth <- replace(annot2$truth, annot2$truth == "T_cell", 6)

write.csv(annot, "Lamprecht_1a_celltype_clean.csv")
write.csv(annot2, "Lamprecht_1a_celltype_clean_new.csv")
write.csv(mat3, "Lamprecht_1a_normalized_filt_new.csv")
