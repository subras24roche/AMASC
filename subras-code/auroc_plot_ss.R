library("ggplot2")
library("reshape2")

print('loaded libraries')

df_raw <- read.table("mountArvados/by_id/2de2948b2d44edcd781487cdfcbaea8c+1220/auclist_bygene_20201021.csv", sep=",", header=TRUE)

df_raw$feature_name<-paste(df_raw$Celltype, df_raw$Feature)


df <- melt(df_raw,id.vars = c("feature_name"), 
		measure.vars=c("PBMC10KNG","PBMC10KV3","CBMC","PBMC1K","PBMC5K","Sort","MALT" )  )

feature_color<-rep("red",nrow(df_raw))
feature_order<-sort(df_raw$feature_name)
feature_color[grep("CD19", feature_order )] <- "black"
feature_color[grep("CD14", feature_order )] <- "black"
feature_color[grep("NCAM1", feature_order )] <- "black"



p<-ggplot(df, aes(x=feature_name, y=value, fill = variable)) + 
  #facet_wrap(~Celltype) + 
  geom_bar(position = "dodge", stat="identity") +
   theme_bw(base_size=25) +
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color=feature_color),
   		 legend.position="bottom") +
   xlab("") + ylab("AUROC") +labs(fill="Data Set")


png("auroc_gene_18_20201021.png",width=1000,height=600,unit="px")
print(p)
dev.off()
print('saved file')
