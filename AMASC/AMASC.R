#AMASC MAIN

DATA_PATH <- "./data"
OUTPUT_PATH <- "./data/amasc/"
threshold_protein <- 7 
threshold_cluster <- 0.40
pg_neighbors <- 50
fs_iteration <- 10
fs_selected <- 9 

files_rna <- dir(path=DATA_PATH, pattern="_rna.csv")
files_protein <- gsub( "_rna.csv" , "_pe.csv", files_rna) 

library("Matrix")
library("Rphenograph")
source("util.R")

start_time <- format(Sys.time(), "%Y-%m-%d_%H:%M")
amasc <- list() 
files_groundtruth <- gsub( "_rna.csv" , "_celltype.csv", files_rna) 
files_fs <- gsub( "_rna.csv" , "_features.csv", files_rna) 
files_accuracy <- gsub( "_rna.csv" , "_accuracy.csv", files_rna) 

print(files_fs)

for(i in 1:length(files_rna)){
    amasc[[i]] <- list()

    amasc[[i]]$protein  <- data.matrix(read.table(paste( DATA_PATH, files_protein[i], sep="/"),sep=",",row.name=1,header=TRUE, stringsAsFactors=FALSE))
    amasc[[i]]$groundtruth <- label_celltype(amasc[[i]]$protein, threshold_protein=threshold_protein, threshold_cluster=threshold_cluster, neighbors=pg_neighbors)
   
    df <- data.frame(truth=amasc[[i]]$groundtruth)
    write.table(df, file = paste(OUTPUT_PATH, files_groundtruth[i], sep="/"), sep=",",quote=FALSE )
    system(paste("python fs.py", paste( DATA_PATH, files_rna[i], sep="/"), 
                               paste(OUTPUT_PATH, files_groundtruth[i], sep="/"), 
                               paste(OUTPUT_PATH, files_fs[i], sep="/"), 
                               paste(OUTPUT_PATH, files_accuracy[i], sep="/"), 
                               fs_iteration), wait=TRUE )

    cat("Done\n")
}

geneset_list <- list()
for(fs_itr in 1:length(files_fs)){
    amasc[[fs_itr]]$features<-read.table(paste(OUTPUT_PATH, files_fs[fs_itr], sep="/"), sep=",", stringsAsFactors=FALSE)
    geneset_list[[fs_itr]] <- names( which(table(amasc[[fs_itr]]$features$V3) > fs_selected) )
}

feature_all <- names( which( table( unlist(geneset_list) ) >=2) )
write(feature_all, file=paste0("features_",start_time,".txt"), sep="\t")

cat("Done\n")
