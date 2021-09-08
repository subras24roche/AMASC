library("Matrix")
#library("DelayedArray")
#print('loaded matrix')
#library("scater")
#print('loaded scater')
#library("SingleCellExperiment")
#print('loaded sce')
library("SingleR")

blueprint.se <- BlueprintEncodeData(rm.NA = "rows")
dice.se <- DatabaseImmuneCellExpressionData()

print('loaded data')

start_time_load <- Sys.time()
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
#mat <- log(mat+1) #scater::logNormCounts(as.matrix(summary(mat)))
sce <- SingleCellExperiment(assays=list(counts=mat))
sce_lognorm <- logNormCounts(sce)
end_time_load <- Sys.time()

load_time <- end_time_load-start_time_load
print(load_time)

start_time_blue <- Sys.time()
pred.blueprint <- SingleR(test = sce_lognorm@assays@data$logcounts, ref = blueprint.se, labels = blueprint.se$label.main)
end_time_blue <- Sys.time()

blue_time <- end_time_blue-start_time_blue
print(blue_time)

start_time_dice <- Sys.time()
pred.dice <- SingleR(test = sce_lognorm@assays@data$logcounts, ref = dice.se, labels = dice.se$label.main)
end_time_dice <- Sys.time()

dice_time <- end_time_dice-start_time_dice
print(dice_time)

save(pred.blueprint,pred.dice,load_time,blue_time,dice_time,file="pbmc10k_singler.rda")

write.table(pred.blueprint, file="pbmc10k_singler_blue_log.csv",sep=",")
write.table(pred.dice, file="pbmc10k_singler_dice_log.csv",sep=",")
