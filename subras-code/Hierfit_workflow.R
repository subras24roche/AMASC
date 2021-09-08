install.packages("devtools")
devtools::install_github("yasinkaymaz/HieRFIT")
library(HieRFIT)
library(celldex)
library(Seurat)
library(DiagrammeR)
pbmc <- readRDS("pbmc3k_final.rds")
pbmc <- UpdateSeuratObject(pbmc)

table(pbmc@meta.data$ClusterNames_0.6)
treeTable <- read.delim("pbmc3k_taxa.txt", header = F)
PlotTopoTree(obj = treeTable)

refmod <- CreateHieR(RefData = pbmc[["RNA"]]@data,
                     ClassLabels = pbmc@meta.data$ClusterNames_0.6,
                     Tree = treeTable,
                     species = "hsapiens")

SaveHieRMod(refMod = refmod, filePrefix = "PBMC3K_HierMod")

PlotTopoNodeAcc(refMod = refmod)

new.pbmc.data <- Read10X("pbmc_10k_v3/filtered_feature_bc_matrix/")

newPBMC <- CreateSeuratObject(counts = new.pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)

newPBMC <- NormalizeData(newPBMC)

newPBMC <- FindVariableFeatures(newPBMC, selection.method = "vst", nfeatures = 2000)
newPBMC <- ScaleData(newPBMC)
newPBMC <- RunPCA(newPBMC)
newPBMC <- FindNeighbors(newPBMC, dims = 1:10)
newPBMC <- FindClusters(newPBMC, resolution = 1)

refmod <- readRDS("PBMC3K_HierMod.RDS")

hierObj <- HieRFIT(Query = newPBMC[["RNA"]]@data, refMod = refmod)

PlotBarStats(HieRobj = hierObj)

PlotTopoStats(HieRobj = hierObj)

CrossCheck(HieRobj = hierObj, Prior = newPBMC@meta.data$res.1)



## PBMC workflow ##
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

table(pbmc_seurat)
treeTable <- read.delim("pbmc3k_taxa.txt", header = F)
PlotTopoTree(obj = treeTable)

refmod <- CreateHieR(RefData = pbmc[["RNA"]]@data,
                     ClassLabels = pbmc@meta.data$ClusterNames_0.6,
                     Tree = treeTable,
                     species = "hsapiens")

SaveHieRMod(refMod = refmod, filePrefix = "PBMC3K_HierMod")

PlotTopoNodeAcc(refMod = refmod)

new.pbmc.data <- Read10X("pbmc_10k_v3/filtered_feature_bc_matrix/")

newPBMC <- CreateSeuratObject(counts = new.pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)

newPBMC <- NormalizeData(newPBMC)

newPBMC <- FindVariableFeatures(newPBMC, selection.method = "vst", nfeatures = 2000)
newPBMC <- ScaleData(newPBMC)
newPBMC <- RunPCA(newPBMC)
newPBMC <- FindNeighbors(newPBMC, dims = 1:10)
newPBMC <- FindClusters(newPBMC, resolution = 1)

refmod <- readRDS("PBMC3K_HierMod.RDS")

hierObj <- HieRFIT(Query = newPBMC[["RNA"]]@data, refMod = refmod)

PlotBarStats(HieRobj = hierObj)

PlotTopoStats(HieRobj = hierObj)

CrossCheck(HieRobj = hierObj, Prior = newPBMC@meta.data$res.1)
