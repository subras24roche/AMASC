library("scPred")
library("Seurat")
library("scater")
#library("beachmat")
library("scuttle")
library("magrittr")
library("SingleR")
library("Matrix")
library("celldex")
library("SingleCellExperiment")
library("scRNAseq")

#scPred tutorial workflow
reference <- scPred::pbmc_1
query <- scPred::pbmc_2

reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)

reference <- getFeatureSpace(reference, "cell_type")
reference <- trainModel(reference)

get_scpred(reference)
plot_probabilities(reference)

reference <- trainModel(reference, model = "mda", reclassify = c("cMono", "ncMono"))

query <- NormalizeData(query)
query <- scPredict(query, reference)
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")

query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(query, group.by = "cell_type", label = TRUE, repel = TRUE)

FeaturePlot(query, c("scpred_B.cell", "scpred_CD4.T.cell", "scpred_CD8.T.cell", 
                     "scpred_cMono", "scpred_ncMono", "scpred_Plasma.cell", 
                     "scpred_cDC", "scpred_pDC"))

crossTab(query, "cell_type", "scpred_prediction")



# scPred PBMC workflow
matrix_dir <- "pbmc10k_filtered/pbmc10k_filtered/"
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
mat <- mat[!duplicated(rownames(mat)),]
sce <- SingleCellExperiment(assays=list(counts=mat))
sce <- logNormCounts(sce)
sce@assays@data$counts <- as(sce@assays@data$counts, 'dgCMatrix')
sce@assays@data$logcounts <- as(sce@assays@data$logcounts, 'dgCMatrix')

pbmc_seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")
dice.se <- celldex::DatabaseImmuneCellExpressionData()


dice.sce <- as(dice.se, "SingleCellExperiment")
dice_seurat <- as.Seurat(dice.sce, counts = "logcounts", data = "logcounts")
dice_seurat <- dice_seurat %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)


# scPred pancreas (Mur) workflow
baron <- scRNAseq::BaronPancreasData()
baron <- logNormCounts(baron)
baron_seurat <- as.Seurat(baron, counts = "counts", data = "logcounts")
mur <- scRNAseq::MuraroPancreasData()
mur <- logNormCounts(mur)

mur_t <- which(!is.na(mur$label))
#mur$label <- mur$label.replace("duct","ductal")
mur <- mur[,mur_t]

#rownames(mur) <- gsub("_.*","",rownames(mur))
mur <- mur[!duplicated(rownames(mur)),]
mur_seurat <- as.Seurat(mur, counts = "counts", data = "logcounts")

baron_seurat <- baron_seurat %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

DimPlot(baron_seurat, group.by = "label", label = TRUE, repel = TRUE)

baron_seurat <- getFeatureSpace(baron_seurat, "label")
baron_seurat <- trainModel(baron_seurat)

get_scpred(baron_seurat)
plot_probabilities(baron_seurat)

mur_seurat <- NormalizeData(mur_seurat)
mur_seurat <- scPredict(mur_seurat, baron_seurat)
DimPlot(mur_seurat, group.by = "scpred_prediction", reduction = "scpred")

mur_seurat <- RunUMAP(mur_seurat, reduction = "scpred", dims = 1:30)
DimPlot(mur_seurat, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(mur_seurat, group.by = "label", label = TRUE, repel = TRUE)

FeaturePlot(mur_seurat, c("scpred_acinar", "scpred_alpha", "scpred_beta", 
                     "scpred_delta", "scpred_gamma", "scpred_endothelial", 
                     "scpred_epsilon", "scpred_schwann", "scpred_mast", 
                     "scpred_activated_stellate", "scpred_macrophage",
                     "scpred_t_cell", "scpred_ductal", "scpred_quiescent_stellate"))

crossTab(mur_seurat, "label", "scpred_prediction")

mur_label <- mur_seurat@meta.data$label
mur_pred <- mur_seurat@meta.data$scpred_prediction

acc <- mur_label == mur_pred
sum(acc)/length(acc)


# scPred pancreas (Seg) workflow
baron <- scRNAseq::BaronPancreasData()
baron <- logNormCounts(baron)
baron_seurat <- as.Seurat(baron, counts = "counts", data = "logcounts")
seg <- scRNAseq::SegerstolpePancreasData()
seg <- logNormCounts(seg)

seg_t <- which(!is.na(seg$`cell type`))
seg <- seg[,seg_t]

seg <- seg[!duplicated(rownames(seg)),]
seg_seurat <- as.Seurat(seg, counts = "counts", data = "logcounts")

baron_seurat <- baron_seurat %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

DimPlot(baron_seurat, group.by = "label", label = TRUE, repel = TRUE)

baron_seurat <- getFeatureSpace(baron_seurat, "label")
baron_seurat <- trainModel(baron_seurat)

get_scpred(baron_seurat)
plot_probabilities(baron_seurat)

seg_seurat <- NormalizeData(seg_seurat)
seg_seurat <- scPredict(seg_seurat, baron_seurat)
DimPlot(seg_seurat, group.by = "scpred_prediction", reduction = "scpred")

seg_seurat <- RunUMAP(seg_seurat, reduction = "scpred", dims = 1:30)
DimPlot(seg_seurat, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(seg_seurat, group.by = "cell.type", label = TRUE, repel = TRUE)

FeaturePlot(seg_seurat, c("scpred_acinar", "scpred_alpha", "scpred_beta", 
                          "scpred_delta", "scpred_gamma", "scpred_endothelial", 
                          "scpred_epsilon", "scpred_schwann", "scpred_mast", 
                          "scpred_activated_stellate", "scpred_macrophage",
                          "scpred_t_cell", "scpred_ductal", "scpred_quiescent_stellate"))

crossTab(seg_seurat, "cell.type", "scpred_prediction")

seg_label <- seg_seurat@meta.data$cell.type
seg_label <- gsub(" .*","",seg_label)
seg_label <- gsub("unclassified","unassigned",seg_label)
seg_pred <- seg_seurat@meta.data$scpred_prediction

acc <- seg_label == seg_pred
sum(acc)/length(acc)
