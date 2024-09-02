#sizes of myelonets and functional marker expressions for myelonets for further analysis
#new set of IDS samples

data_ss <- read.csv("E:/cellcycle/cellcycle_phenotypes_TIM3.csv")
#clinical <- read.csv("P:/h345/afarkkilab/Projects/Sciset/clinical_data.csv")

obs_all <- data_ss

#first RCN rich in Macrophages
#load in results from delaunay triangulation computed in cpu1

macrophage_blobs_all <- list()
for (i in unique(data_ss$imageid)){
  
  macrophage_blob <- read.csv(paste0("E:/cellcycle/Giotto/",i,"_spatial_kmeans_13.csv"))
  
  
  macrophage_blob <- unique(macrophage_blob[,-c(1)])
  
  macrophage_blobs_all[[i]] <- macrophage_blob
  
}


macrophage_blobs_all <- do.call(rbind, macrophage_blobs_all)
library(dplyr)
blobs <- list()
for (i in unique(macrophage_blobs_all$imageid)){
  
  blob <- macrophage_blobs_all[which(macrophage_blobs_all$imageid == i),] %>% 
    group_by(components.membership) %>% 
    summarise(number = n())
  blob$imageid <- i
  
  blobs[[i]] <- blob
}


blobs <- do.call(rbind, blobs)

blobs <- blobs[-which(blobs$number<9),]

#then boxplot mean size per sample

#mean per sample

blobs_mean <- blobs[,-c(1)] %>%
  group_by(imageid)%>%
  summarize(mean_n = mean(number))

blobs_median <- blobs[,-c(1)] %>%
  group_by(imageid)%>%
  summarize(median_n = median(number))


#add clinical

blobs_mean <- merge(blobs_mean, clinical, by.x="imageid", by.y="block name")

blobs_median <- merge(blobs_median, clinical, by.x="imageid", by.y="block name")

library(ggplot2)
library(ggpubr)
pdf("E:/cellcycle/plots/boxplots/boxplot_CD163_blob.pdf", width = 4, height = 4)
ggplot(blobs_mean, aes(x=PFIcat_12, y=mean_n)) + geom_boxplot(aes(fill=PFIcat_12)) + geom_point() + scale_fill_manual(values=c("lightblue", "blue")) + stat_compare_means(method="wilcox") + theme_bw()
ggplot(blobs_median, aes(x=PFIcat_12, y=median_n)) + geom_boxplot(aes(fill=PFIcat_12)) + geom_point() + scale_fill_manual(values=c("lightblue", "blue")) + stat_compare_means(method="wilcox") + theme_bw()
dev.off()

#repeat the previous for other RCNs