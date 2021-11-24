###############################################################################################################################
######################### BOX PLOTS CELLS PERCENTAGES #########################################################################
###############################################################################################################################

#For Figures 2 and 3 and Supplementary Figures 2 and 3 

library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)


#ANNOTATION DATA
all_data_TMA <- read.csv("TMA_clinicaldata.csv")
ann <- all_data_TMA[, c("Identifier","Category", "PFI_time")]
ann$PFI_time <- as.numeric(ann$PFI_time)
ann$PFI_time[which(ann$`PFI_time`> 365)] <- "long"
ann$PFI_time[which(ann$`PFI_time`> 182 & ann$`PFI_time`< 365)] <- "short"
ann$PFI_time[which(ann$`PFI_time`< 182)] <- "short"
ann$Category <- as.character(ann$Category)
ann$Category[which(ann$`Category`== "HR")] <- "HRwt"
ann$Category[which(ann$`Category`== "BRCA1")] <- "BRCA1/2 mutated"
ann$Category[which(ann$`Category`== "BRCA2")] <- "BRCA1/2 mutated"


#SINGLE CELL DATA
all_celltypes <- read.csv("TMA_annotated_single_cell_data.csv")

ann <- ann[is.element(ann$Identifier, all_celltypes$Sample),]

all_clusters_sorted <- all_celltypes %>%
  group_by(Sample, GlobalCellType) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))

all_clusters_sorted$freq <- all_clusters_sorted$freq*100
all_clusters_sorted <- all_clusters_sorted[, -3]
all_clusters_sorted_wide <- pivot_wider(all_clusters_sorted, names_from="Sample", values_from="freq")
all_clusters_sorted_wide <- as.data.frame(all_clusters_sorted_wide)
rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide$GlobalCellType
all_clusters_sorted_wide <- all_clusters_sorted_wide[,-1]

all_clusters_sorted_wide <- as.data.table(t(all_clusters_sorted_wide), keep.colnames = T, keep.rownames = T)

rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide$rn
colnames(all_clusters_sorted_wide)[1] <- "Patient"
all_clusters_sorted_wide[is.na(all_clusters_sorted_wide)] <- 0
all_clusters_sorted_wide <- as.data.frame(all_clusters_sorted_wide)
rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide[, 1]
all_clusters_sorted_wide <- all_clusters_sorted_wide[, -1]

ann <- ann[is.element(ann$Identifier, all_celltypes_24092020$Sample),]
all_clusters_sorted_wide <- cbind(all_clusters_sorted_wide, ann)


theme <- theme(panel.border = element_rect(colour = "black", size=1, fill=NA), 
               panel.background = element_blank(), plot.title=element_text(hjust=0.5)) 


colnames(all_clusters_sorted_wide) <- make.names(colnames(all_clusters_sorted_wide))
for (i in colnames(all_clusters_sorted_wide)[1:26]){
  p <- ggplot(all_clusters_sorted_wide, aes_string(x="Type", y=i)) + geom_boxplot(lwd=1.5, aes(fill=`Type`, alpha=0.9)) + xlab("HR status")+theme+
    geom_jitter(position=position_jitter(0.2), size=1, alpha=0.9) + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), aspect.ratio=1, legend.position="none")+ 
    scale_fill_manual(values=c("#3366CC", "#CC0033")) + ylab(i) + stat_compare_means(method = "wilcox.test",label.x=1.65,label.y.npc = "top", size=5)
  
  print(p)
}



#NOW THE SAME BUT IMNMUNE CELL PROPORTIONS OUT OF IMMUNE CELLS

all_celltypes_immune <- all_celltypes[which(all_celltypes$Subtype == "Immune"),]

all_clusters_sorted <- all_celltypes_immune %>%
  group_by(Sample, GlobalCellType) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))

all_clusters_sorted$freq <- all_clusters_sorted$freq*100
all_clusters_sorted <- all_clusters_sorted[, -3]
all_clusters_sorted_wide <- pivot_wider(all_clusters_sorted, names_from="Sample", values_from="freq")
all_clusters_sorted_wide <- as.data.frame(all_clusters_sorted_wide)
rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide$GlobalCellType
all_clusters_sorted_wide <- all_clusters_sorted_wide[,-1]

all_clusters_sorted_wide <- as.data.table(t(all_clusters_sorted_wide), keep.colnames = T, keep.rownames = T)

rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide$rn
colnames(all_clusters_sorted_wide)[1] <- "Patient"
all_clusters_sorted_wide[is.na(all_clusters_sorted_wide)] <- 0
all_clusters_sorted_wide <- as.data.frame(all_clusters_sorted_wide)
rownames(all_clusters_sorted_wide) <- all_clusters_sorted_wide[, 1]
all_clusters_sorted_wide <- all_clusters_sorted_wide[, -1]

ann <- ann[is.element(ann$Identifier, all_celltypes_24092020$Sample),]
all_clusters_sorted_wide <- cbind(all_clusters_sorted_wide, ann)


theme <- theme(panel.border = element_rect(colour = "black", size=1, fill=NA), 
               panel.background = element_blank(), plot.title=element_text(hjust=0.5)) 


colnames(all_clusters_sorted_wide) <- make.names(colnames(all_clusters_sorted_wide))
for (i in colnames(all_clusters_sorted_wide)[1:10]){
  p <- ggplot(all_clusters_sorted_wide, aes_string(x="Type", y=i)) + geom_boxplot(lwd=1.5, aes(fill=`Type`, alpha=0.9)) + xlab("HR status")+theme+
    geom_jitter(position=position_jitter(0.2), size=1, alpha=0.9) + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), aspect.ratio=1, legend.position="none")+ 
    scale_fill_manual(values=c("#3366CC", "#CC0033")) + ylab(i) + stat_compare_means(method = "wilcox.test",label.x=1.65,label.y.npc = "top", size=5)
  
  print(p)
}





