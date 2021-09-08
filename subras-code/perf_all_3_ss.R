df_raw <- read.table("mountArvados/by_id/2de2948b2d44edcd781487cdfcbaea8c+1220/accuracy_all_methods_3.csv", sep=",", header=TRUE)
df_raw$dataset<-toupper(df_raw$dataset)


df <- reshape2::melt(df_raw[,-2],id.vars = c("dataset") ) #, measure.vars=c("PBMC10KV3","PBMC10KNG","CBMC","PBMC1K","PBMC5K","PBMCSORT","MALT" ) 

df$split <- rep("Test",nrow(df))

#df$split[which(df$dataset=="PBMC10KV3")] <- "Training"
df$split[which(df$dataset=="PBMC10KNG")] <- "Training"
#df$split[which(df$dataset=="CBMC")] <- "Training"

df$dataset <- factor(df$dataset, level=c("PBMC10KNG","PBMC10KV3","CBMC","PBMC1K","PBMC5K","PBMCSORT","MALT" ) )

df$split <- factor(df$split, level=c("Training", "Test" ) )

df$value_round <- round(df$value,2)
library("ggplot2")
p<-ggplot(df, aes(x=dataset, y=value, fill = split) ) + 
  facet_wrap(~variable, ncol=3) + 
  geom_bar(position = "dodge", stat="identity") +
  geom_text(aes(label=value_round), vjust=1.6, color="white", size=6)+
   theme_bw(base_size=25) +
   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), 
   legend.position="bottom") + # , color=split
   			xlab("") + ylab("Accuracy") +labs(fill="Data Set")


png("perf_all_18g3_accu.png",width=1000,height=1000,unit="px")
print(p)
dev.off()