
#############################################################################################################################################################
#################################### ALL METACLUSTERS AND IMMUNE CELL SUBTYPES ###############################################################################
############################################################################################################################################################

#For Figures 2, 3, 4 and Supplementary Figure 1

library(tidyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(dplyr)

ann <- read.csv("TMA_clinicaldata.csv")
all_celltypes <- read.csv("TMA_annotated_single_cell_data.csv")
all_celltypes$GlobalCellType <- as.character(all_celltypes$GlobalCellType)
all_celltypes$Subtype <- as.character(all_celltypes$Subtype)

all_celltypes[which(all_celltypes$GlobalCellType == "CD11c+CD163+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD163+CD11c+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "CD11c+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD11c+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "CD163+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD163+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "CD4"), "GlobalCellType"] <- "CD4+ Effector T-cells"
all_celltypes[which(all_celltypes$GlobalCellType == "CD8"), "GlobalCellType"] <- "CD8+ T-cells"
all_celltypes[which(all_celltypes$GlobalCellType == "CD163+Macrophages"), "GlobalCellType"] <- "CD163+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+ Macrophages"
all_celltypes[which(all_celltypes$GlobalCellType == "CD11c+APC"), "GlobalCellType"] <- "CD11c+APC"
all_celltypes[which(all_celltypes$GlobalCellType == "FOXP3+CD4+Tregs"), "GlobalCellType"] <- "FOXP3+CD4+ T-regs"
all_celltypes[which(all_celltypes$GlobalCellType == "High-PDL1"), "GlobalCellType"] <- "Functional stroma"
all_celltypes[which(all_celltypes$GlobalCellType == "Non-proliferative_Stroma"), "GlobalCellType"] <- "Non-proliferative Stroma"
all_celltypes[which(all_celltypes$GlobalCellType == "Low_eccentricity_medium_vimentin"), "GlobalCellType"] <- "Low eccentricity"
all_celltypes[which(all_celltypes$GlobalCellType == "Proliferative_Stroma"), "GlobalCellType"] <- "Proliferative Stroma"
all_celltypes[which(all_celltypes$GlobalCellType == "High-proliferative_Stroma"), "GlobalCellType"] <- "High-proliferative Stroma"
all_celltypes[which(all_celltypes$GlobalCellType == "Hyperfunctional epithelial"), "GlobalCellType"] <- "Functional epithelial"


all_clusters <- all_celltypes[, c(1, 46)]


metacluster_percentages <- all_clusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
metacluster_percentages <- as.data.frame(metacluster_percentages)

metacluster_percentages <- metacluster_percentages[, -c(3)]
metacluster_percentages$proportion <- metacluster_percentages$proportion*100


metacluster_percentages <- metacluster_percentages %>% pivot_wider(names_from = GlobalCellType, values_from = proportion)

metacluster_percentages <- as.data.frame(metacluster_percentages)
rownames(metacluster_percentages) <- metacluster_percentages[, 1]
metacluster_percentages <- metacluster_percentages[, -1]

metacluster_percentages[is.na(metacluster_percentages)] <- 0


ann <- ann[rownames(metacluster_percentages),]

colours<- list("simpson_immune"=col_het,"simpson_tumor"=col_het_t,"Type"=c("BRCA1"="royalblue","BRCA2"="lightblue", "HRwt"="red3"), "PFI_time"=c("long" = "#3399CC", "short"="#003366"))

rownames(ann) <- ann[, 1]
ann <- ann[,-1]

ann <- as.data.frame(ann)

ann[which(ann$PFI_time == "medium"), "PFI_time"] <- "short"

Rowann <- HeatmapAnnotation(df=ann[, c(1:2)], which="row", col=colours, 
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(1, "mm"))
col_c <- colorRamp2(c(-1.5, -0.5, 0.5, 1.5), c( "#0000FF","#AAAAFF", "#F9AAAA", "#EE0000"))
hmap <- Heatmap(scale(as.matrix(metacluster_percentages)), clustering_method_rows ="ward.D",name="subtype %\nper patient (Z-score)", 
                row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                left_annotation = Rowann, row_dend_width = unit(1.3, "cm"),show_row_names = F, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                show_row_dend=T, show_column_dend=T, column_title = "subtypes", row_title = "patients",
                height = unit(10.5, "cm") , width = unit(7.3, "cm"),border="white",
                rect_gp = gpar(col = "white", lwd = 2), col=col_c, row_split=3, column_split = 2)


draw(hmap)


#################################################################################################################################
################################ ALL CELLS HEATMAP ################################################################################
#################################################################################################################################

all_celltypes_24092020 <- sample_n(all_celltypes_24092020, 20000)

all_celltypes_24092020$Subtype <- as.character(all_celltypes_24092020$Subtype)
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Endothelia"), "Subtype"] <- "Endothelia"
all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Stroma_Endothelia"), "Subtype"] <- "Stroma"


colours<- list("df"=c("Tumor"="#666666", "Immune"="#7570B3","Stroma"="#E6AB02", "Endothelia"="#1B9E77"), "HR_defect"=c("1" = "blue", "0" = "red"))
all_celltypes_24092020$Subtype <- as.factor(all_celltypes_24092020$Subtype)
#, col=colours
Rowann <- HeatmapAnnotation(df=all_celltypes_24092020[, c(47)], which="row",col=colours,
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(1, "mm"))
png("all_cells_heatmap_wardD.png", res=300, width=5, height=7, units = "in")
Heatmap(scale(as.matrix(all_celltypes_24092020[, c(8,10, 17, 18, 21, 29, 11:14, 16, 22:28, 30:31)])),name="marker expression per cell", 
        column_names_gp = gpar(fontsize = 6),column_title_gp = gpar(fontsize = 10, fontface = "bold"),clustering_method_rows ="ward.D",
        row_dend_width = unit(1.3, "cm"),show_row_names = F, row_title_gp = gpar(fontsize = 10, fontface = "bold"),left_annotation = Rowann,
        show_row_dend=T, show_column_dend=T, column_title = "markers", row_title = "cells",
        height = unit(10.5, "cm") , width = unit(3.3, "cm"),border="white")
dev.off()



###################################################################################################################################################
########################################### IMMUNE CELL SUBTYPE ABUNDANCE ########################################################################
###################################################################################################################################################

all_cell_types_immune <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Immune"),]
all_clusters <- all_cell_types_immune[, c(1, 46)]



immune_from_immune <- all_clusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
immune_from_immune <- as.data.frame(immune_from_immune)

immune_from_immune <- immune_from_immune[, -c(3)]
immune_from_immune$proportion <- immune_from_immune$proportion*100


immune_from_immune <- immune_from_immune %>% pivot_wider(names_from = GlobalCellType, values_from = proportion)

immune_from_immune <- as.data.frame(immune_from_immune)
rownames(immune_from_immune) <- immune_from_immune[, 1]
immune_from_immune <- immune_from_immune[, -1]

immune_from_immune[is.na(immune_from_immune)] <- 0



ann <- ann[rownames(immune_from_immune),]

pal <- colorRampPalette(c("thistle1", "plum4"))
pal(3)

pal <- colorRampPalette(c("lightblue", "lightblue4"))
pal(3)

col_het = colorRamp2(c(0.4, 0.6, 0.8), c("thistle1", "#C4A3C4", "plum4"))
col_het_t = colorRamp2(c(0.3, 0.55, 0.80), c("lightblue", "#8AADB8", "lightblue4"))


colours<- list("simpson_immune"=col_het,"simpson_tumor"=col_het_t,"Type"=c("BRCA1"="royalblue","BRCA2"="lightblue", "HRwt"="red3"), "PFI_time"=c("long" = "#3399CC", "short"="#003366"))

rownames(ann) <- ann[, 1]
ann <- ann[,-1]

ann <- as.data.frame(ann)

ann[which(ann$PFI_time == "medium"), "PFI_time"] <- "short"

Rowann <- HeatmapAnnotation(df=ann[, c(1:2)], which="row", col=colours, 
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(1, "mm"))
col_c <- colorRamp2(c(-1.5, -0.5, 0.5, 1.5), c( "#0000FF","#AAAAFF", "#F9AAAA", "#EE0000"))
hmap <- Heatmap(scale(as.matrix(immune_from_immune)), name="subtype %\nper patient (Z-score)", 
                row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                left_annotation = Rowann, row_dend_width = unit(1.3, "cm"),show_row_names = F, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                show_row_dend=T, show_column_dend=T, column_title = "subtypes", row_title = "patients",
                height = unit(10.5, "cm") , width = unit(7.3, "cm"),border="white",
                rect_gp = gpar(col = "white", lwd = 2), col=col_c)


draw(hmap)


###################################################################################################################################################
########################################### STROMAL CELL SUBTYPE ABUNDANCE ########################################################################
###################################################################################################################################################

all_cell_types_stroma <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Stroma_Endothelia"),]
all_clusters <- all_cell_types_stroma[, c(1, 46)]



stroma_from_stroma <- all_clusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
stroma_from_stroma <- as.data.frame(stroma_from_stroma)

stroma_from_stroma <- stroma_from_stroma[, -c(3)]
stroma_from_stroma$proportion <- stroma_from_stroma$proportion*100


stroma_from_stroma <- stroma_from_stroma %>% pivot_wider(names_from = GlobalCellType, values_from = proportion)

stroma_from_stroma <- as.data.frame(stroma_from_stroma)
rownames(stroma_from_stroma) <- stroma_from_stroma[, 1]
stroma_from_stroma <- stroma_from_stroma[, -1]

stroma_from_stroma[is.na(stroma_from_stroma)] <- 0



hmap <- Heatmap(scale(as.matrix(stroma_from_stroma)), name="subtype %\nper patient (Z-score)", 
                row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                left_annotation = Rowann, row_dend_width = unit(1.3, "cm"),show_row_names = F, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                show_row_dend=T, show_column_dend=T, column_title = "subtypes", row_title = "patients",
                height = unit(10.5, "cm") , width = unit(7.3, "cm"),border="white",
                rect_gp = gpar(col = "white", lwd = 2), col=col_c)


draw(hmap)



###################################################################################################################################################
########################################### TUMOR CELL SUBTYPE ABUNDANCE ########################################################################
###################################################################################################################################################

all_cell_types_tumor <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]
all_clusters <- all_cell_types_tumor[, c(1, 46)]



tumor_from_tumor <- all_clusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
tumor_from_tumor <- as.data.frame(tumor_from_tumor)

tumor_from_tumor <- tumor_from_tumor[, -c(3)]
tumor_from_tumor$proportion <- tumor_from_tumor$proportion*100


tumor_from_tumor <- tumor_from_tumor %>% pivot_wider(names_from = GlobalCellType, values_from = proportion)

tumor_from_tumor <- as.data.frame(tumor_from_tumor)
rownames(tumor_from_tumor) <- tumor_from_tumor[, 1]
tumor_from_tumor <- tumor_from_tumor[, -1]

tumor_from_tumor[is.na(tumor_from_tumor)] <- 0


hmap <- Heatmap(scale(as.matrix(tumor_from_tumor)), name="subtype %\nper patient (Z-score)", 
                row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                left_annotation = Rowann, row_dend_width = unit(1.3, "cm"),show_row_names = F, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                show_row_dend=T, show_column_dend=T, column_title = "subtypes", row_title = "patients",
                height = unit(10.5, "cm") , width = unit(7.3, "cm"),border="white",
                rect_gp = gpar(col = "white", lwd = 2), col=col_c)


draw(hmap)



###################################################################################################################################################
########################################### IMMUNE CELL MARKER EXPRESSION ########################################################################
###################################################################################################################################################

immune <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Immune"),]

data <-  immune %>% 
  group_by(GlobalCellType) %>%
  summarise(
    cCasp3 = mean(cCasp3),
    pSTAT1 = mean(pSTAT1),
    Ki67 = mean(Ki67),
    PDL1 = mean(PDL1),
    PD1 = mean(PD1),
    P21 = mean(P21)
  )


data_scaled <- scale(data[, c(2:7)])
data <- cbind(data[, c(1)], data_scaled)



data <- as.data.frame(data)
rownames(data) <- data[, c(1)]
data <- data[, -c(1)]


data <- data[c(8, 7, 9, 5, 3, 10, 4, 6, 2, 1),]

hmap <- Heatmap(as.matrix(data), name="Z-score", 
                show_row_names=T, show_column_names=T, cluster_columns=T,cluster_rows=F,show_row_dend=T,
                show_column_dend=T,  row_dend_reorder=F, column_dend_reorder=F,
                border="white", rect_gp = gpar(col = "white", lwd = 3), width = unit(4, "cm"), 
                height = unit(7, "cm"), column_names_rot = 45, column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10))
draw(hmap)


###################################################################################################################################################
########################################### STROMAL CELL SUBTYPE MARKER EXPRESSION ########################################################################
###################################################################################################################################################


stroma <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Stroma" | all_celltypes_24092020$Subtype == "Endothelia"),]

data <-  stroma %>% 
  group_by(GlobalCellType) %>%
  summarise(
    cCasp3 = mean(cCasp3),
    pSTAT1 = mean(pSTAT1),
    Ki67 = mean(Ki67),
    PDL1 = mean(PDL1),
    Vimentin = mean(vimentin),
    CD31 = mean(CD31),
    Eccentricity = mean(Eccentricity),
    P21 = mean(P21)
  )


data_scaled <- scale(data[, c(2:9)])
data <- cbind(data[, c(1)], data_scaled)



data <- as.data.frame(data)
rownames(data) <- data[, c(1)]
data <- data[, -c(1)]

subtype_numbers <- stroma %>% group_by(GlobalCellType) %>% summarise(n=n())
subtype_numbers <- as.data.frame(subtype_numbers)
rownames(subtype_numbers) <- subtype_numbers[, 1]
rows <- rowAnnotation(cells=anno_barplot(subtype_numbers$n, height = unit(1, "cm")))



data <- data[c(4,9, 8, 5, 7, 6, 2, 3, 1),]

hmap <- Heatmap(as.matrix(data), name="Z-score", 
                show_row_names=T, show_column_names=T, cluster_columns=T,cluster_rows=T,show_row_dend=F,
                show_column_dend=T,  row_dend_reorder=F, column_dend_reorder=F,
                border="white", rect_gp = gpar(col = "white", lwd = 3), width = unit(6, "cm"), 
                height = unit(7, "cm"), column_names_rot = 45, column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10))
draw(hmap)

###################################################################################################################################################
########################################### TUMOR CELL MARKER EXPRESSION ########################################################################
###################################################################################################################################################

tumor <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]
subtype_numbers <- tumor %>% group_by(GlobalCellType) %>% summarise(n=n())
subtype_numbers <- as.data.frame(subtype_numbers)
rownames(subtype_numbers) <- subtype_numbers[, 1]
subtype_numbers <- subtype_numbers[c(3, 7, 2, 6, 5, 4, 1),]
rows <- rowAnnotation(cells=anno_barplot(subtype_numbers$n, height = unit(1, "cm")))


tumor <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]

data <-  tumor %>% 
  group_by(GlobalCellType) %>%
  summarise(
    cCasp3 = mean(cCasp3),
    pSTAT1 = mean(pSTAT1),
    Ki67 = mean(Ki67),
    PDL1 = mean(PDL1),
    Ecadherin = mean(Ecadherin),
    Vimentin = mean(vimentin),
    P21 = mean(P21),
    CK7 = mean(CK7)
  )


data_scaled <- scale(data[, c(2:9)])
data <- cbind(data[, c(1)], data_scaled)



data <- as.data.frame(data)
rownames(data) <- data[, c(1)]
data <- data[, -c(1)]

#right_annotation = rows,
#own code

data <- data[c(3, 7, 2, 6, 5, 4, 1),]

hmap <- Heatmap(as.matrix(data), name="Z-score", 
                show_row_names=T, show_column_names=T, cluster_columns=T,cluster_rows=F,show_row_dend=T,
                show_column_dend=T,  row_dend_reorder=F, column_dend_reorder=F,
                border="white", rect_gp = gpar(col = "white", lwd = 3), width = unit(6, "cm"), 
                height = unit(6, "cm"), column_names_rot = 45, column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10))
draw(hmap)












