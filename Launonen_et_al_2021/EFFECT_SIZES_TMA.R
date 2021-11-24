########################################################################################################################
####################################  EFFECT SIZES #####################################################################
########################################################################################################################

#For figures 2, 3 as well as Supplementary Figures 1, 2 and 3

library(tidyr)
library(rstatix)
library(dplyr)

ann <- read.csv("TMA_clinicaldata.csv")
ann <- ann[, c("Identifier","Category", "PFI_time")]

ann$PFI_time[which(ann$`PFI_time`> 365)] <- "long"
ann$PFI_time[which(ann$`PFI_time`> 182 & ann$`PFI_time`< 365)] <- "medium"
ann$PFI_time[which(ann$`PFI_time`< 182)] <- "short"
ann$PFI_time[which(ann$`PFI_time`== "56")] <- "short"
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

all_clusters <- all_celltypes[, c(1, 46)]

metacluster_percentages <- all_clusters %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
metacluster_percentages <- as.data.frame(metacluster_percentages)

metacluster_percentages <- metacluster_percentages[, -c(3)]
metacluster_percentages$proportion <- metacluster_percentages$proportion*100
#to wide format

metacluster_percentages <- metacluster_percentages %>% pivot_wider(names_from = GlobalCellType, values_from = proportion)

metacluster_percentages <- as.data.frame(metacluster_percentages)
rownames(metacluster_percentages) <- metacluster_percentages[, 1]
metacluster_percentages <- metacluster_percentages[, -1]

metacluster_percentages[is.na(metacluster_percentages)] <- 0


ann <- ann[rownames(metacluster_percentages),]
BRCA_patients <- rownames(ann[which(ann$`Category` == "BRCA1/2 mutated"),])
HRwt_patients <- rownames(ann[which(ann$`Category` == "HRwt"),])

metacluster_percentages <- cbind(metacluster_percentages, ann)



colnames(metacluster_percentages) <- make.names(colnames(metacluster_percentages))

profepi <- wilcox_effsize(
  metacluster_percentages,
  Proliferating.epithelial ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

profEMT <- wilcox_effsize(
  metacluster_percentages,
  Proliferating.EMT ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

EMT <- wilcox_effsize(
  metacluster_percentages,
  EMT ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

metacluster_percentages$total_EMT <- metacluster_percentages$EMT + metacluster_percentages$Proliferating.EMT

metacluster_percentages$total_macrophages <- metacluster_percentages$CD163..Macrophages + metacluster_percentages$IBA1.CD163..Macrophages + metacluster_percentages$IBA1..Macrophages

total_mac <- wilcox_effsize(
  metacluster_percentages,
  total_macrophages ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

total_EMT <- wilcox_effsize(
  metacluster_percentages,
  total_EMT ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

CD163 <- wilcox_effsize(
  metacluster_percentages,
  CD163..Macrophages ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

IBA1CD163 <- wilcox_effsize(
  metacluster_percentages,
  IBA1.CD163..Macrophages ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)


IBA1CD163CD11c <- wilcox_effsize(
  metacluster_percentages,
  IBA1.CD163.CD11c..Macrophages ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

metacluster_percentages$total_tumor <- metacluster_percentages[, 1] + metacluster_percentages[, 7] +metacluster_percentages[, 9] +metacluster_percentages[, 11] +metacluster_percentages[, 22] +metacluster_percentages[, 24] +metacluster_percentages[, 25]

metacluster_percentages$total_stroma <- metacluster_percentages[, 8] +metacluster_percentages[, 12] +metacluster_percentages[, 13] +metacluster_percentages[, 14] +metacluster_percentages[, 15] +metacluster_percentages[, 20] +metacluster_percentages[, 21] +metacluster_percentages[,23 ] +metacluster_percentages[, 26]


metacluster_percentages$tumor_stroma_ratio <- metacluster_percentages$total_tumor/metacluster_percentages$total_stroma


tumor_stroma <- wilcox_effsize(
  metacluster_percentages,
  tumor_stroma_ratio ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

IBA1 <- wilcox_effsize(
  metacluster_percentages,
  IBA1..Macrophages ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)


eff_sizes <- rbind(profepi, profEMT, EMT, total_EMT, total_mac, CD163, IBA1CD163, IBA1CD163CD11c, tumor_stroma, IBA1)



#effect sizes immune percentages and ratios

immune_from_immune <- read.csv("immune_percentages_from_immune.csv")


immune_from_immune <- cbind(immune_from_immune, ann)

CD8 <- wilcox_effsize(
  metacluster_percentages,
  CD8..T.cells ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

CD11APC <- wilcox_effsize(
  immune_from_immune,
  CD11c..APC ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

IBA1CD163CD11c <- wilcox_effsize(
  immune_from_immune,
  IBA1.CD163.CD11c..Macrophages ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

immune_from_immune$CD8_IBA1CD163CD11 <- immune_from_immune$CD8..T.cells/immune_from_immune$IBA1.CD163.CD11c..Macrophages

immune_from_immune$CD8_IBA1CD11 <- immune_from_immune$CD8..T.cells/immune_from_immune$IBA1.CD11c..Macrophages

CD8_IBA1CD11c <- wilcox_effsize(
  immune_from_immune,
  CD8_IBA1CD11 ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

CD8_IBA1CD163CD11c <- wilcox_effsize(
  immune_from_immune,
  CD8_IBA1CD163CD11 ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)

simpson_immune <- wilcox_effsize(
  immune_from_immune,
  simpson_immune ~ Category,
  paired = FALSE,
  alternative = "two.sided",
  ci = TRUE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000)


eff_sizes_immune <- rbind(simpson_immune, CD8_IBA1CD11c, CD8_IBA1CD163CD11c, CD8, CD11APC, IBA1CD163CD11c)

