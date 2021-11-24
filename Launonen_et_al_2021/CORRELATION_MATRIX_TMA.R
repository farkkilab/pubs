####################################################################################################################
######################## CORRELATION DOTPLOT #######################################################################
####################################################################################################################

# For Figure 4

library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(corrplot)
library(RColorBrewer)

#without fdr adjustment

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
BRCA_patients <- rownames(ann[which(ann$`Category` == "BRCA1/2 mutated"),])
HRwt_patients <- rownames(ann[which(ann$`Category` == "HRwt"),])




all_data <- metacluster_percentages


all_data <- all_data[, c(25, 9, 7, 24, 22, 1, 11,10, 5, 6, 18, 19, 16, 17, 3, 4, 2, 23, 20, 15, 8, 26, 14, 12, 21, 13)]

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#matrix of the p-value of the correlation
#p.mat <- cor.mtest(res[BRCA_patients,])
#p.mat <- cor.mtest(res[HRwt_patients,])

BRCA_patients <- as.character(BRCA_patients)
res <- all_data
res[is.na(res)] <- 0
p.mat <- cor.mtest(res[BRCA_patients,])
res <- cor(res[BRCA_patients,], method = "spearman")
diag(res) = NA


col_c<- colorRampPalette(c("blue", "white", "red"))(20)
corrplot(res, method="circle", col=col_c, tl.col="black", tl.srt=45,
         p.mat = p.mat, sig.level = 0.05, insig = "blank", type="upper",
         addgrid.col="grey",  tl.cex=0.7, tl.pos = "tl", bg="white", na.label = " ")

res <- all_data
HRwt_patients <- as.character(HRwt_patients)
p.mat <- cor.mtest(res[HRwt_patients,])
res <- cor(res[HRwt_patients,], method = "spearman")
diag(res) = NA


corrplot(res, method="circle", col=col_c, tl.col="black", tl.srt=45,
         p.mat = p.mat, sig.level = 0.05, insig = "blank", type="lower",
         addgrid.col="grey",  tl.cex=0.7, tl.pos = "n",add=T, cl.pos="n",
         bg="lavenderblush", rect.lwd = 2, na.label = "*")


#### NOW FDR CORRECTION

BRCA_patients <- as.character(BRCA_patients)
res <- all_data
res[is.na(res)] <- 0
p.mat <- cor.mtest(res[BRCA_patients,])
res <- cor(res[BRCA_patients,], method = "spearman")
diag(res) = NA

res <- all_data
HRwt_patients <- as.character(HRwt_patients)
p.mat_hr <- cor.mtest(res[HRwt_patients,])

p.mat[1,] <- p.mat_hr[1,]
p.mat[2,c(2:26)] <- p.mat_hr[2,c(2:26)]
p.mat[3,c(3:26)] <- p.mat_hr[3,c(3:26)]
p.mat[4,c(4:26)] <- p.mat_hr[4,c(4:26)]
p.mat[5,c(5:26)] <- p.mat_hr[5,c(5:26)]
p.mat[6,c(6:26)] <- p.mat_hr[6,c(6:26)]
p.mat[7,c(7:26)] <- p.mat_hr[7,c(7:26)]
p.mat[8,c(8:26)] <- p.mat_hr[8,c(8:26)]
p.mat[9,c(9:26)] <- p.mat_hr[9,c(9:26)]
p.mat[10,c(10:26)] <- p.mat_hr[10,c(10:26)]
p.mat[11,c(11:26)] <- p.mat_hr[11,c(11:26)]
p.mat[12,c(12:26)] <- p.mat_hr[12,c(12:26)]
p.mat[13,c(13:26)] <- p.mat_hr[13,c(13:26)]
p.mat[14,c(14:26)] <- p.mat_hr[14,c(14:26)]
p.mat[15,c(15:26)] <- p.mat_hr[15,c(15:26)]
p.mat[16,c(16:26)] <- p.mat_hr[16,c(16:26)]
p.mat[17,c(17:26)] <- p.mat_hr[17,c(17:26)]
p.mat[18,c(18:26)] <- p.mat_hr[18,c(18:26)]
p.mat[19,c(19:26)] <- p.mat_hr[19,c(19:26)]
p.mat[20,c(20:26)] <- p.mat_hr[20,c(20:26)]
p.mat[21,c(21:26)] <- p.mat_hr[21,c(21:26)]
p.mat[22,c(22:26)] <- p.mat_hr[22,c(22:26)]
p.mat[23,c(23:26)] <- p.mat_hr[23,c(23:26)]
p.mat[24,c(24:26)] <- p.mat_hr[24,c(24:26)]
p.mat[25,c(25:26)] <- p.mat_hr[25,c(25:26)]
p.mat[26,26] <- p.mat_hr[26,26]

p.mat <- p.mat %>% 
  as.vector %>% 
  p.adjust(method='fdr')  %>%
  matrix(ncol=26)


BRCA_patients <- as.character(BRCA_patients)
res <- all_data_try
res[is.na(res)] <- 0
res <- cor(res[BRCA_patients,], method = "spearman")
diag(res) = NA



col_c<- colorRampPalette(c("blue", "white", "red"))(20)
corrplot(res, method="circle", col=col_c, tl.col="black", tl.srt=45,
         p.mat = p.mat, sig.level = 0.1, insig = "blank", type="lower",
         addgrid.col="grey",  tl.cex=0.7, tl.pos = "tl", bg="white", na.label = " ")



res <- all_data_try
HRwt_patients <- as.character(HRwt_patients)
res <- cor(res[HRwt_patients,], method = "spearman")
diag(res) = NA


corrplot(res, method="circle", col=col_c, tl.col="black", tl.srt=45,
         p.mat = p.mat, sig.level = 0.1, insig = "blank", type="upper",
         addgrid.col="grey",  tl.cex=0.7, tl.pos = "n",add=T, cl.pos="n",
         bg="lavenderblush", rect.lwd = 2, na.label = "*")


