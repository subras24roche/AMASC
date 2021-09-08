library("Matrix")
library("SingleR")

dice.se <- DatabaseImmuneCellExpressionData()
dice_normalized_filt <- assay(dice.se)
dice_celltype_clean <- colData(dice.se)

#remove columns 2 and 3 from labels
dice_celltype_clean <- subset(dice_celltype_clean, select = label.main)
dice_celltype_clean <- as.data.frame(dice_celltype_clean)
colnames(dice_celltype_clean) <- "truth"

write.csv(dice_celltype_clean, "dice_celltype_clean.csv")
write.csv(dice_normalized_filt, "dice_normalized_filt.csv")
