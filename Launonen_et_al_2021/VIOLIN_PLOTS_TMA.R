###########################################################################################################################
############################################# VIOLIN PLOTS ################################################################
###########################################################################################################################

#For Figure 5 and Supplementary Figure 5

library(ggpubr)
library(ggplot2)


ann <- read.csv("TMA_clinicaldata.csv")
ann <- ann[, c("Identifier","Category", "PFI_time")]

ann$PFI_time[which(ann$`PFI_time`> 365)] <- "long"
ann$PFI_time[which(ann$`PFI_time`> 182 & ann$`PFI_time`< 365)] <- "medium"
ann$PFI_time[which(ann$`PFI_time`< 182)] <- "short"
ann$Category <- as.character(ann$Category)
ann$Category[which(ann$`Category`== "HR")] <- "HRwt"
ann$Category[which(ann$`Category`== "BRCA1")] <- "BRCA1/2 mutated"
ann$Category[which(ann$`Category`== "BRCA2")] <- "BRCA1/2 mutated"
rownames(ann) <- ann$Identifier

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


#PD1

cells_of_interest_2 <- all_celltypes[which(all_celltypes$GlobalCellType == "CD4+ Effector T-cells" | all_celltypes$GlobalCellType == "CD8+ T-cells" |all_celltypes$GlobalCellType == "FOXP3+CD4+ T-regs"),]

cells_of_interest_2$HR_defect <- as.factor(cells_of_interest_2$HR_defect)
cells_of_interest_2$Sample <- as.factor(cells_of_interest_2$Sample)

cells_of_interest_2 <- merge(cells_of_interest_2, ann, by.x="Sample", by.y="Identifier")


cells_of_interest_2$combined <- paste0(cells_of_interest_2$GlobalCellType, cells_of_interest_2$HR_defect)

my_comparisons <- list(c("CD8+ T-cells0", "CD8+ T-cells1"), c("CD4+ Effector T-cells0", "CD4+ Effector T-cells1"), c("FOXP3+CD4+ T-regs0", "FOXP3+CD4+ T-regs1"))

pvals <- compare_means(PD1~combined, cells_of_interest_2, method = "wilcox", comparisons=my_comparisons, p.adjust.method = "BH")

dodge <- position_dodge(width = 0.5)

p <- ggplot(cells_of_interest_2, aes(x=combined, y=PD1)) + geom_violin(aes(fill=Type), position=dodge) + geom_boxplot(width=0.1, aes(fill=Type), position=dodge, outlier.shape = NA) + theme(axis.text.x = element_text(angle = 25, hjust=1)) + xlab("subtype") + scale_fill_manual(values=c("blue","red"))+ stat_compare_means(pvals, label = "p.signif", method = "wilcox", comparisons=my_comparisons, label.y = 8.2) + scale_y_continuous(limits=c(7.4, 8.4))
p


##########################################################################################################################################



#PDL1

cells_of_interest <- all_celltypes[which(all_celltypes$GlobalCellType == "CD163+ Macrophages" | all_celltypes$GlobalCellType == "IBA1+CD163+ Macrophages" |all_celltypes$GlobalCellType == "IBA1+ Macrophages" |all_celltypes$GlobalCellType == "IBA1+CD163+CD11c+ Macrophages"|all_celltypes$GlobalCellType == "IBA1+CD11c+ Macrophages"|all_celltypes$GlobalCellType == "CD11c+APC" | all_celltypes$Subtype == "Tumor" | all_celltypes$Subtype == "Stroma_Endothelia"),]

cells_of_interest$HR_defect <- as.factor(cells_of_interest$HR_defect)
cells_of_interest$Sample <- as.factor(cells_of_interest$Sample)

ann$Sample <- rownames(ann)

cells_of_interest <- merge(cells_of_interest, ann, by.x="Sample", by.y="Identifier")


cells_of_interest[which(cells_of_interest$Subtype == "Tumor"), "GlobalCellType"] <- "Tumor"
cells_of_interest[which(cells_of_interest$Subtype == "Stroma_Endothelia"), "GlobalCellType"] <- "Stromal"



cells_of_interest$combined <- paste0(cells_of_interest$GlobalCellType, cells_of_interest$HR_defect)

my_comparisons <- list(c("CD11c+APC0", "CD11c+APC1"), c("CD163+ Macrophages0", "CD163+ Macrophages1"), c("IBA1+ Macrophages0", "IBA1+ Macrophages1"), c("IBA1+CD11c+ Macrophages0", "IBA1+CD11c+ Macrophages1"), c("IBA1+CD163+ Macrophages0", "IBA1+CD163+ Macrophages1"), c("IBA1+CD163+CD11c+ Macrophages0", "IBA1+CD163+CD11c+ Macrophages1"), c("Stromal0", "Stromal1"), c("Tumor0", "Tumor1"))


dodge <- position_dodge(width = 0.5)

pvals <- compare_means(PDL1~combined, cells_of_interest3, method = "wilcox", comparisons=my_comparisons, p.adjust.method = "BH")

all_macrophages <- cells_of_interest[-which(cells_of_interest$Subtype == "Stromal" | cells_of_interest$Subtype == "Tumor"),]
all_macrophages$GlobalCellType <- "all_macrophages"
all_macrophages$combined <- paste0(all_macrophages$GlobalCellType, all_macrophages$HR_defect)

cells_of_interest3 <- rbind(cells_of_interest, all_macrophages)

#add all macrophages vs stroma and tumor


order_x <- c( "IBA1+CD11c+ Macrophages0","IBA1+CD11c+ Macrophages1","CD11c+APC0","CD11c+APC1", "IBA1+CD163+CD11c+ Macrophages0", "IBA1+CD163+CD11c+ Macrophages1", "IBA1+CD163+ Macrophages0","IBA1+CD163+ Macrophages1", "CD163+ Macrophages0","CD163+ Macrophages1", "IBA1+ Macrophages0","IBA1+ Macrophages1","Stromal0","Stromal1", "Tumor0", "Tumor1")


p <- ggplot(cells_of_interest, aes(x=combined, y=PDL1)) + geom_violin(aes(fill=Type), position=dodge) + geom_boxplot(width=0.1, aes(fill=Type), position=dodge, outlier.shape = NA) + theme(axis.text.x = element_text(angle = 25, hjust=1)) + xlab("Macrophage/tumor/stroma subtype") + scale_fill_manual(values=c("blue","red"))+ stat_compare_means(pvals,label = "p.signif", method = "wilcox", comparisons=my_comparisons, label.y = 8.2) + scale_y_continuous(limits=c(7.4, 8.4))
p




########################################################################################################################################################



