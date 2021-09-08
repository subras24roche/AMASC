#ml R/4.0.1-intel-2018b
#GCC compiler builtin_FILE error in 3.6

#require(devtools) 
#install_version("xgboost", version = "0.90.0.2", repos = "http://cran.us.r-project.org")
library('readr')

g <- read_csv('seg_genelist.csv', col_names = FALSE)
genes <- as.character(g$X2)

ge <- read.csv('../home_jupyter/segerstolpe_pancreas_counts_clean.csv', row.names = 1)
label_gex <- read.csv('../home_jupyter/segerstolpe_pancreas_annotations_clean.csv', row.names = 1)
ge_gex <- data.matrix( ge[genes,] )


minmaxscaler <- function(mat, axis= 1){
  apply(mat, axis, function(x) {
    (x-min(x))/(max(x)-min(x)) } )
}


binarizer <- function(mat, bin_threshold ){
  ifelse(mat > bin_threshold, 1, 0)
}

normalize_dataset <- function(ge, log1ptransform = FALSE, minmaxscale = TRUE, bin_threshold = 0.01){

  if(log1ptransform == TRUE)
    ge <- log2(ge+1)

  if(minmaxscale == TRUE)
    ge <- minmaxscaler(ge)

  ge <- binarizer(ge, bin_threshold)

  return(ge)
}


#precision <- posPredValue(predictions, y, positive="1")
#recall <- sensitivity(predictions, y, positive="1")

#F1 <- (2 * precision * recall) / (precision + recall)

f1_score <- function(predicted, expected, positive.class="1") {
    predicted <- factor(as.character(predicted), levels=unique(as.character(expected)))
    expected  <- as.factor(expected)
    cm = as.matrix(table(expected, predicted))

    precision <- diag(cm) / colSums(cm)
    recall <- diag(cm) / rowSums(cm)
    f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))

    #Assuming that F1 is zero when it's not possible compute it
    f1[is.na(f1)] <- 0

    #Binary F1 or Multi-class macro-averaged F1
    #ifelse(nlevels(expected) == 2, f1[positive.class], mean(f1))
	mean(f1)
}

calculate.accuracy <- function(predictions, ref.labels) {
    return(length(which(predictions == ref.labels)) / length(ref.labels))
}

ge_gex_norm <- normalize_dataset(ge_gex)

ge_train <- ge_gex_norm
label_train <- label_gex
label_train$truth <- as.character(label_train$truth)

#ge_train <- ge_gex_norm
#label_train <- label_gex

set.seed(20200819)
split_train = sample(nrow(label_train), nrow(label_train)*.8)

library("caret")
library("xgboost")

mat_maxperf <- matrix(0,length(genes), length(genes))
rownames(mat_maxperf) <- genes
colnames(mat_maxperf) <- 1:length(genes)


mat_maxperf_valid <- matrix(0,length(genes), length(genes))
rownames(mat_maxperf_valid) <- genes
colnames(mat_maxperf_valid) <- 1:length(genes)

print("progress")

gene_max_prev <- NULL

for(j in 1:length(genes)){
	for(i in 1:length(genes)){
		if(genes[i] %in% gene_max_prev == FALSE ){
			feature_list_this <- unique(c(gene_max_prev, genes[i]))
			if(j>1){
				df_this <- ge_train[split_train,feature_list_this ]
			}else{
				df_this <-t(t(ge_train[split_train,feature_list_this ]))
			}
			bstSparse <- xgboost(data = df_this, label = label_train$truth[split_train], 
			max.depth = 3, nthread = 4, nrounds = 50, objective = "multi:softmax", num_class = 15, verbose = 0 )#eta = 1 

			pred <- predict(bstSparse,t(t(ge_train[split_train, feature_list_this])))
			mat_maxperf[genes[i],j]<- calculate.accuracy(pred, label_train$truth[split_train])

			pred_valid <- predict(bstSparse,t(t(ge_train[-split_train, feature_list_this])))
			mat_maxperf_valid[genes[i],j]<- calculate.accuracy(pred_valid, label_train$truth[-split_train])

			cat("Round ",j, " gene ", i, " ",feature_list_this, " ", round(mat_maxperf[genes[i],j],2), "\n")
		}else{
			cat(genes[i], "has been selected\n")
		}

	}
    
	if(j>1){
		#if( max(mat_maxperf[,j]) > max(mat_maxperf[,j-1]) ){
			gene_max_prev<-unique(c( gene_max_prev, names(which.max(mat_maxperf[,j])) ))
			cat("Round ",j, " max ", round(max(mat_maxperf[,j]),2)," ", gene_max_prev,"\n")
		#}else{
			"Saturation\n"
		#}
	}else{
		gene_max_prev<-unique(c( gene_max_prev, names(which.max(mat_maxperf[,j])) ))
		cat("Round ",j, " max ", round(max(mat_maxperf[,j]),2)," ", gene_max_prev,"\n")
	}
}


save(split_train, genes, gene_max_prev, mat_maxperf, mat_maxperf_valid, file="geneselection_xgboost_accuracy_10kgex_valid_80_3_01_6.rda")



#plot(1:ncol(mat_maxperf), apply(mat_maxperf, 2,  max) )
#abline(h=18)
df <- data.frame(feature_num=1:ncol(mat_maxperf),accuracy_train=apply(mat_maxperf, 2,  max),accuracy_test=apply(mat_maxperf_valid, 2,  max),
				feature=factor(gene_max_prev, levels=gene_max_prev))

#library("ggrepel")

p <-ggplot(df, aes(x= feature, group = 1) ) + #colour="green",
   #geom_point(size = 2,alpha = 0.6) +
   geom_line(aes(y = accuracy_train), color = "darkred") + 
   geom_line(aes(y = accuracy_test), color = "steelblue") + 
    #geom_text_repel(aes(label = feature ),  #color = 'white', fill = factor(cyl), geom_label_repel
    #               angle=90, min.segment.length = 0, size = 3, force=1) +
  geom_hline(yintercept=0.9074879, linetype="dashed", color = "red") +
  geom_vline(xintercept=18, linetype="dashed", color = "red") +
  #geom_text(x=20, y=0.95, label="0.91", size=5) + 
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Features") + ylab("Accuracy") + ylim(0, 1) # +
  #geom_text(aes(label=Name),hjust=0, vjust=0)

png("xgboost_progressive_selection_seg.png",width=1000,height=600,unit="px")

print(p)

dev.off()