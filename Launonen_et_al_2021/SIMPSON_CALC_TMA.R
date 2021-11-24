###################################################################################################################
####################### SIMPSON DIVERSITY INDEX CALCULATION #######################################################
###################################################################################################################

#For Figures 2, 3 and Supplementary Figures 1, 2 and 3

library(Hmisc)
library(vegan)
library(sm)
library(ISLR)
library(ggpubr)
library(dplyr)
library(tidyr)

all_celltypes <- read.csv("TMA_annotated_single_cell_data.csv")

metaclusters <- all_celltypes[which(all_celltypes$Subtype == "Tumor"),]

metacluster_numbers <- metaclusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n())
metacluster_numbers <- as.data.frame(metacluster_numbers)


metacluster_numbers <- metacluster_numbers %>% pivot_wider(names_from = Metacluster, values_from = n)

metacluster_numbers[is.na(metacluster_numbers)] <- 0

simpson <- diversity(metacluster_numbers[,2:8], "simpson",MARGIN = 1)
simpson_tumor <- data.frame(sample=metacluster_numbers$Sample, simpson)


## immune
immuneclusters <- all_celltypes[which(all_celltypes$Subtype == "Immune"),]

immune_numbers <- immuneclusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n())
immune_numbers <- as.data.frame(immune_numbers)

immune_numbers <- immune_numbers %>% pivot_wider(names_from = GlobalCellType, values_from = n)

immune_numbers[is.na(immune_numbers)] <- 0


simpson.i <- diversity(immune_numbers[2:9],"simpson" ,MARGIN = 1)
simpson.i <- data.frame(sample=immune_numbers$Sample, simpson.i)

##stroma


stromalclusters <- all_celltypes[which(all_celltypes$Subtype == "Stroma_Endothelia"),]

stroma_numbers <- stromalclusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n())
stroma_numbers <- as.data.frame(stroma_numbers)

stroma_numbers <- stroma_numbers %>% pivot_wider(names_from = GlobalCellType, values_from = n)

stroma_numbers[is.na(stroma_numbers)] <- 0


simpson.s <- diversity(stroma_numbers[2:9],"simpson" ,MARGIN = 1)
simpson.s <- data.frame(sample=stroma_numbers$Sample, simpson.s)

