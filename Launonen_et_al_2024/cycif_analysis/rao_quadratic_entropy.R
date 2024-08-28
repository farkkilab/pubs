#Rao's quadratic entropy for Figure 3

# SYNCSA package rao.diversity inputs:
# comm
# Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data. Alternatively comm can be an object of class metacommunity.data, an alternative way to set all data.frames/matrices. When you use the class metacommunity.data the arguments traits, phylodist and put.together must be null. See details.
# 
# traits
# Data frame or matrix data of species described by traits, with traits as columns and species as rows 


#install.packages('SYNCSA')
library(SYNCSA)
library(dplyr)
library(readxl)
library(tidyr)
#as comm, calculate proportion of cell types per sample
#as traits, calculate median functional marker expression per cell type

#per sample
#per neighborhood per sample

data <- read.csv("D:/sciset/data_20231125.csv")
#first per RCN

comm <- data %>% group_by(neighbordood_cluster2, GlobalCellType2) %>%  summarize(n=n())%>% mutate(proportion = n / sum(n))

#cell types as columns

comm <- pivot_wider(comm[,-3], values_from = proportion, names_from = GlobalCellType2)
comm[is.na(comm)] <- 0
comm <- as.data.frame(comm)
rownames(comm) <- comm[,1]
comm <- comm[,-1]

#expression data also

data_functional <- data[ ,c("GlobalCellType2",'TIM3','CD45RO' ,'yH2AX','Annexin','MHCI','MHCII','pSTAT1','FOXOA3','TAZ','HE4','Ki67')]

traits <- data_functional %>% group_by(GlobalCellType2) %>% summarize(mean_TIM3 = mean(TIM3),
                                                           mean_CD45RO = mean(CD45RO),
                                                           mean_yH2AX = mean(yH2AX),
                                                           mean_Annexin = mean(Annexin),
                                                           mean_MHCI = mean(MHCI),
                                                           mean_MHCII = mean(MHCII),
                                                           mean_pSTAT1 = mean(pSTAT1),
                                                           mean_FOXOA3 = mean(FOXOA3),
                                                           mean_TAZ = mean(TAZ),
                                                           mean_HE4 = mean(HE4),
                                                           mean_Ki67 = mean(Ki67))

traits <- as.data.frame(traits)
rownames(traits) <- traits[,1]
traits <- traits[,-1]


res <- rao.diversity(
  comm,
  traits,
  phylodist = NULL,
  checkdata = TRUE,
  ord = "metric",
  put.together = NULL,
  standardize = TRUE
)

#Fun Rao in the output has the diversity values

results <- as.data.frame(res$FunRao)
colnames(results) <- "Rao"
results$neighbordood_cluster2 <- rownames(results)
ggplot(results, aes(y=Rao, x=factor(neighbordood_cluster2, levels=c("Proliferating.epithelial","epithelial_and_proliferating_epithelial", "Epithelial", "Epithelial_EMT", "EMT", "Proliferating.EMT", "tumor-stroma-interface", "Desmin.positive","SMA.Desmin.myofibroblast",  "stroma",  "SMA.CD31.positive",  "Myofibroblast",  "Fibroblast", "Macrophage","IBA1.CD163.Macrophages","CD11c.myeloid", "CD8_CD4_T.cells","Immune")))) + geom_bar(stat='identity') + theme_bw()


#then separately for each patient to compare pre and post stromal cell diversity as suggested by reviewer 1

comm <- data %>% group_by(Sample_code, neighbordood_cluster2, GlobalCellType2) %>%  summarize(n=n())%>% mutate(proportion = n / sum(n))

comm$neighbordood_cluster2 <- as.character(comm$neighbordood_cluster2)
comm[which(comm$neighbordood_cluster2 == "CD8_CD4_T.cells"), "neighbordood_cluster2"] <- "CD8.CD4.T.cells"
comm[which(comm$neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"), "neighbordood_cluster2"] <- "epithelial.proliferating.epithelial"
comm[which(comm$neighbordood_cluster2 == "Epithelial_EMT"), "neighbordood_cluster2"] <- "Epithelial.EMT"


#cell types as columns

comm <- pivot_wider(comm[,-4], values_from = proportion, names_from = c(Sample_code, neighbordood_cluster2))
comm[is.na(comm)] <- 0
comm <- as.data.frame(comm)
rownames(comm) <- comm[,1]
comm <- comm[,-1]
comm <- t(comm)

#rows should be patient + neighborhood
#traits 

res <- rao.diversity(
  comm,
  traits,
  phylodist = NULL,
  checkdata = TRUE,
  ord = "metric",
  put.together = NULL,
  standardize = TRUE
)

#Fun Rao has the entropy values

results <- as.data.frame(res$FunRao)
colnames(results) <- "Rao"
results$sample_neighbordood_cluster <- rownames(results)

#have to separate sample code from neighborhood

df <- results %>% separate(sample_neighbordood_cluster, c("Patient", "Sample_code", "neighbordood_cluster2"), sep="_")
df$Sample_code <- paste0(df$Patient, "_", df$Sample_code)

clinical <- read.csv("P:/h345/afarkkilab/Projects/Sciset/clinical_data.csv")
df <- merge(df, clinical, by.x="Sample_code", by.y="Sample_code")

#p-values

df_paired <- df[which(df$Paired_sample == TRUE),]

wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "tumor-stroma-interface" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "tumor-stroma-interface" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "CD11c.myeloid" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "CD11c.myeloid" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Fibroblast" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Fibroblast" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Proliferating.epithelial" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Proliferating.epithelial" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "SMA.CD31.positive" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "SMA.CD31.positive" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "SMA.Desmin.myofibroblast" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "SMA.Desmin.myofibroblast" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Myofibroblast" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Myofibroblast" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "stroma" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "stroma" & df_paired$Stage == "primary"), "Rao"], paired = T)

wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "EMT" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "EMT" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Epithelial" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Epithelial" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Epithelial.EMT" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Epithelial.EMT" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "epithelial.proliferating.epithelial" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "epithelial.proliferating.epithelial" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Proliferating.EMT" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Proliferating.EMT" & df_paired$Stage == "primary"), "Rao"], paired = T)

wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "CD8.CD4.T.cells" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "CD8.CD4.T.cells" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Immune" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Immune" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Macrophages" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Macrophages" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "IBA1.CD163.Macrophages" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "IBA1.CD163.Macrophages" & df_paired$Stage == "primary"), "Rao"], paired = T)
wilcox.test(df_paired[which(df_paired$neighbordood_cluster2 == "Desmin.positive" & df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$neighbordood_cluster2 == "Desmin.positive" & df_paired$Stage == "primary"), "Rao"], paired = T)



#then only stromal cells

stroma <- data[which(data$GlobalCellType2 == "Desmin.positive.cell" | data$GlobalCellType2 == "SMA.Desmin.positive.cell" | data$GlobalCellType2 == "SMA.CD31.positive.cell" | data$GlobalCellType2 == "Fibroblast" | data$GlobalCellType2 == "Myofibroblast" | data$GlobalCellType2 == "Endothelial.cell"),]



comm <- stroma %>% group_by(Sample_code, GlobalCellType2) %>%  summarize(n=n())%>% mutate(proportion = n / sum(n))

#cell types as columns

comm <- pivot_wider(comm[,-3], values_from = proportion, names_from = c(Sample_code))
comm[is.na(comm)] <- 0
comm <- as.data.frame(comm)
rownames(comm) <- comm[,1]
comm <- comm[,-1]
comm <- t(comm)

#rows should be patient + neighborhood

traits <- stroma %>% group_by(GlobalCellType2) %>% summarize(mean_TIM3 = mean(TIM3),
                                                                      mean_CD45RO = mean(CD45RO),
                                                                      mean_yH2AX = mean(yH2AX),
                                                                      mean_Annexin = mean(Annexin),
                                                                      mean_MHCI = mean(MHCI),
                                                                      mean_MHCII = mean(MHCII),
                                                                      mean_pSTAT1 = mean(pSTAT1),
                                                                      mean_FOXOA3 = mean(FOXOA3),
                                                                      mean_TAZ = mean(TAZ),
                                                                      mean_HE4 = mean(HE4),
                                                                      mean_Ki67 = mean(Ki67))

traits <- as.data.frame(traits)
rownames(traits) <- traits[,1]
traits <- traits[,-1]


res <- rao.diversity(
  comm,
  traits,
  phylodist = NULL,
  checkdata = TRUE,
  ord = "metric",
  put.together = NULL,
  standardize = TRUE
)

results <- as.data.frame(res$FunRao)
colnames(results) <- "Rao"
results$Sample_code <- rownames(results)

clinical <- read.csv("P:/h345/afarkkilab/Projects/Sciset/clinical_data.csv")
df <- merge(results, clinical, by.x="Sample_code", by.y="Sample_code")

#p-values

df_paired <- df[which(df$Paired_sample == TRUE),]

wilcox.test(df_paired[which(df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$Stage == "primary"), "Rao"], paired = T)

#

#same for tumor

tumor <- data[which(data$GlobalCellType2 == "Proliferating.epithelial" | data$GlobalCellType2 == "Proliferating.EMT" | data$GlobalCellType2 == "EMT" | data$GlobalCellType2 == "Epithelial"),]

comm <- tumor %>% group_by(Sample_code, GlobalCellType2) %>%  summarize(n=n())%>% mutate(proportion = n / sum(n))

#cell types as columns

comm <- pivot_wider(comm[,-3], values_from = proportion, names_from = c(Sample_code))
comm[is.na(comm)] <- 0
comm <- as.data.frame(comm)
rownames(comm) <- comm[,1]
comm <- comm[,-1]
comm <- t(comm)

#rows should be patient + neighborhood

traits <- tumor %>% group_by(GlobalCellType2) %>% summarize(mean_yH2AX = mean(yH2AX),
                                                             mean_Annexin = mean(Annexin),
                                                             mean_MHCI = mean(MHCI),
                                                             mean_MHCII = mean(MHCII),
                                                             mean_pSTAT1 = mean(pSTAT1),
                                                             mean_FOXOA3 = mean(FOXOA3),
                                                             mean_TAZ = mean(TAZ),
                                                             mean_HE4 = mean(HE4),
                                                             mean_Ki67 = mean(Ki67))

traits <- as.data.frame(traits)
rownames(traits) <- traits[,1]
traits <- traits[,-1]


res <- rao.diversity(
  comm,
  traits,
  phylodist = NULL,
  checkdata = TRUE,
  ord = "metric",
  put.together = NULL,
  standardize = TRUE
)

results <- as.data.frame(res$FunRao)
colnames(results) <- "Rao"
results$Sample_code <- rownames(results)

clinical <- read.csv("P:/h345/afarkkilab/Projects/Sciset/clinical_data.csv")
df <- merge(results, clinical, by.x="Sample_code", by.y="Sample_code")

#p-values

df_paired <- df[which(df$Paired_sample == TRUE),]

wilcox.test(df_paired[which(df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$Stage == "primary"), "Rao"], paired = T)



#same for immune


immune <- data[which(data$GlobalCellType2 == "IBA1.CD11c.Macrophages" | data$GlobalCellType2 == "IBA1.CD163.Macrophages" | data$GlobalCellType2 == "CD163.Macrophages" | data$GlobalCellType2 == "CD11c.myeloid"| data$GlobalCellType2 == "FOXP3.CD4.Tregs"| data$GlobalCellType2 == "CD8.T.cells"| data$GlobalCellType2 == "CD4.T.cells"),]

comm <- immune %>% group_by(Sample_code, GlobalCellType2) %>%  summarize(n=n())%>% mutate(proportion = n / sum(n))

#cell types as columns

comm <- pivot_wider(comm[,-3], values_from = proportion, names_from = c(Sample_code))
comm[is.na(comm)] <- 0
comm <- as.data.frame(comm)
rownames(comm) <- comm[,1]
comm <- comm[,-1]
comm <- t(comm)

#rows should be patient + neighborhood

traits <- immune %>% group_by(GlobalCellType2) %>% summarize(mean_yH2AX = mean(yH2AX),
                                                            mean_Annexin = mean(Annexin),
                                                            mean_MHCI = mean(MHCI),
                                                            mean_MHCII = mean(MHCII),
                                                            mean_pSTAT1 = mean(pSTAT1),
                                                            mean_TIM3 = mean(TIM3),
                                                            mean_TAZ = mean(TAZ),
                                                            mean_CD45RO = mean(CD45RO),
                                                            mean_Ki67 = mean(Ki67))

traits <- as.data.frame(traits)
rownames(traits) <- traits[,1]
traits <- traits[,-1]


res <- rao.diversity(
  comm,
  traits,
  phylodist = NULL,
  checkdata = TRUE,
  ord = "metric",
  put.together = NULL,
  standardize = TRUE
)

results <- as.data.frame(res$FunRao)
colnames(results) <- "Rao"
results$Sample_code <- rownames(results)

clinical <- read.csv("P:/h345/afarkkilab/Projects/Sciset/clinical_data.csv")
df <- merge(results, clinical, by.x="Sample_code", by.y="Sample_code")

#p-values

df_paired <- df[which(df$Paired_sample == TRUE),]

wilcox.test(df_paired[which(df_paired$Stage == "interval"), "Rao"], df_paired[which(df_paired$Stage == "primary"), "Rao"], paired = T)


























