#heatmap figure 3f

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(viridis)
#E:/sciset/all_myeloid_RCN_myelonets/All_RCN_Macrophage_blobs_size_mean_functional_marker_expression_and_proportions.csv

Macrophages <- read.csv("D:/Sciset/scimap/plots/dalaunay/Macrophage_blobs_size_mean_functional_marker_expression.csv")
Macrophages2 <- read.csv("D:/Sciset/scimap/plots/dalaunay/Immune_blobs_size_mean_functional_marker_expression.csv")
Macrophages3 <- read.csv("D:/Sciset/scimap/plots/dalaunay/IBA1CD163_Macrophage_blobs_size_mean_functional_marker_expression.csv")
Macrophages4 <- read.csv("D:/Sciset/scimap/plots/dalaunay/CD11c_blobs_size_mean_functional_marker_expression.csv")

#names

Macrophages$blob <- "Macrophages"
Macrophages2$blob <- "Immune"
Macrophages3$blob <- "IBA1CD163_Macrophages"
Macrophages4$blob <- "CD11c"


Macrophages_all <- rbind(Macrophages, Macrophages2, Macrophages3, Macrophages4)

# Rowann <- HeatmapAnnotation(df=Macrophages_all[, c(1, 3, 13, 21, 11, 27)], which="row",  
#                             annotation_name_gp = gpar(fontsize=10,fontface = "bold"))
# Heatmap(scale(Macrophages_all[,c(4:9)]), height = unit(10.5, "cm") , width = unit(7.3, "cm"), row_dend_width = unit(3, "cm"),
#         left_annotation = Rowann, row_split = 6)
# 




Macrophages_all$number_log <- log2(Macrophages_all$number)
Var = circlize::colorRamp2(seq(0, 300, length = 100), rev(hcl.colors(100,"Heat")))
colours<- list("Stage"=c("primary"="royalblue","interval"="yellow"), "blob"=c("CD11c" = "#0a2463", "IBA1CD163_Macrophages"="#3e92cc", "Immune"="#fffaff", "Macrophages" = "#d8315b"), "HRD_status" = c("HRD" = "blue", "HRP" = "red"), "number_log" = Var)
#colours<- list("Stage"=c("primary"="royalblue","interval"="yellow"), "blob"=c("CD11c" = "#ff9898ff", "IBA1CD163_Macrophages"="#d9636cff", "Immune"="#a91e45ff", "Macrophages" = "#691238ff"), "HRD_status" = c("HRD" = "blue", "HRP" = "red"), "number" = Var)
Rowann <- HeatmapAnnotation(df=Macrophages_all[which(Macrophages_all$number>9), c(1, 3, 13, 21, 27)], which="row",  
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), col=colours)
col_c <- colorRamp2(c(-1.5, -0.5, 0.5, 1.5), c( "#0000FF","#AAAAFF", "#F9AAAA", "#EE0000"))
#col_c <- colorRamp2(c(0, 0.3, 0.6,1), c( "#0000FF","#AAAAFF", "#F9AAAA", "#EE0000"))

hmap1 <- Heatmap(scale(Macrophages_all[which(Macrophages_all$number>9),c(4, 6:9)]), height = unit(10.5, "cm") , width = unit(3.3, "cm"), row_dend_width = unit(3, "cm"),
                 left_annotation = Rowann, show_row_names = F, col=col_c)
hmap1


