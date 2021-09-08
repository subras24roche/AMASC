library('scRNAseq')
library('SingleCellExperiment')
library('SingleR')
library('scater')

baron <- scRNAseq::BaronPancreasData()
baron_se <- logNormCounts(baron)

mur <- scRNAseq::MuraroPancreasData()
mur_se <- logNormCounts(mur)
rownames(mur_se@assays@data$logcounts) <- gsub("_.*","",rownames(mur_se@assays@data$logcounts))
rownames(mur_se@assays@data$counts) <- gsub("_.*","",rownames(mur_se@assays@data$counts))

'
matrix_dir <- "byrnes_mouse_data/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2
#mat <- log(mat+1) #scater::logNormCounts(as.matrix(summary(mat)))
sce <- SingleCellExperiment(assays=list(counts=mat))
sce_lognorm <- logNormCounts(sce)
'

pred.mur <- SingleR(test = mur_se@assays@data$counts, ref = baron_se, labels = baron_se$label)

write.table(pred.mur, file="muraro_singler_baron_log.csv",sep=",")
