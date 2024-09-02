library(readxl)
library(ggplot2)
library(ggpubr)

clinical <- read_excel("D:/validation_cycif/clinical_data_validation_cycif.xlsx")
clinical <- read_excel("E:/cellcycle/clinical_data_validation_cycif.xlsx")

clinical <- read_excel("D:/validation_cycif/clinical_data_validation_cycif.xlsx")
clinical <- read_excel("/media/oncosys/Expansion/cellcycle/clinical_data_validation_cycif.xlsx")


clinical <- clinical[, c("block name", "PFIcat_12")]
clinical <- clinical[c(1:30),]

data_ss <- read.csv('/media/oncosys/Expansion/cellcycle/cellcycle_phenotypes_RCNs.csv')
data_ss <- read.csv('/media/oncosys/Expansion/cellcycle/cellcycle_phenotypes_TIM3.csv')

data_ss[which(data_ss$phenotype == "CD8_new"), "phenotype"] <- "CD8Tcells"
data_ss[which(data_ss$phenotype == "Tumor_new2"), "phenotype"] <- "Tumor"
data_ss[which(data_ss$phenotype == "Tumor_new"), "phenotype"] <- "Tumor"
data_ss[which(data_ss$phenotype == "CD11c_new"), "phenotype"] <- "CD11c.myeloid"

#sit data_ss

data_ss <- read.csv("E:/cellcycle/cellcycle_phenotypes_TIM3.csv")

#CD163 expression in all macrophages

CD163_median <- data_ss[which(data_ss$phenotype == "IBA1.CD11.mac" |data_ss$phenotype == "IBA1.CD163.mac" |data_ss$phenotype == "CD163.mac" ),] %>%
  group_by(imageid) %>% summarize(CD163 = median(CD163))


pfi <- merge(CD163_median, clinical, by.x="imageid", by.y="block name")

#library(MoMacolors)

p <- ggplot(pfi, aes(x=PFIcat_12, y=CD163, fill=PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue"))
p
pdf("E:/cellcycle/plots/boxplots/CD163_PFI.pdf", width=4, height = 4)
p
dev.off()

#TIM3 CD8Tcells
TIM3_median <- data_ss[which(data_ss$phenotype == "CD8Tcells"),] %>%
  group_by(imageid) %>% summarize(TIM3 = median(TIM.3))

pfi <- merge(TIM3_median, clinical, by.x="imageid", by.y="block name")

p <- ggplot(pfi, aes(x=PFIcat_12, y=TIM3, fill=PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue"))
p
pdf("/media/oncosys/Expansion/cellcycle/plots/boxplots/TIM3_PFI.pdf", width=4, height = 4)
p
dev.off()


CD163_median <- data_ss[which(data_ss$phenotype == "IBA1.CD11.mac" |data_ss$phenotype == "IBA1.CD163.mac" |data_ss$phenotype == "CD163.mac"| data_ss$phenotype == "CD11c.myeloid"),] %>%
  group_by(imageid) %>% summarize(CD163 = median(CD163))


pfi <- merge(CD163_median, clinical, by.x="imageid", by.y="block name")

library(ggplot2)
library(ggpubr)
#library(MoMacolors)

p <- ggplot(pfi, aes(x=PFIcat_12, y=CD163, fill=PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue"))
p
pdf("E:/cellcycle/plots/boxplots/median_CD163_from_all_myeloid_PFI.pdf", width=4, height = 4)
p
dev.off()

#proportion of CD163 from all myeloid

CD163_prop <- data_ss[which(data_ss$phenotype == "IBA1.CD11.mac" |data_ss$phenotype == "IBA1.CD163.mac" |data_ss$phenotype == "CD163.mac"| data_ss$phenotype == "CD11c.myeloid"),] %>%
  group_by(imageid, phenotype) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))

CD163_prop <- CD163_prop[which(CD163_prop$phenotype == "CD163.mac"),]

pfi <- merge(CD163_prop, clinical, by.x="imageid", by.y="block name")

library(ggplot2)
library(ggpubr)
#library(MoMacolors)

p <- ggplot(pfi, aes(x=PFIcat_12, y=freq, fill=PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue"))
p
pdf("E:/cellcycle/plots/boxplots/CD163_prop_from_all_myeloid_PFI.pdf", width=4, height = 4)
p
dev.off()

#entä proportion of CD163 from all macs

CD163_prop <- data_ss[which(data_ss$phenotype == "IBA1.CD11c.mac" |data_ss$phenotype == "IBA1.CD163.mac" |data_ss$phenotype == "CD163.mac"),] %>%
  group_by(imageid, phenotype) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))

CD163_prop <- CD163_prop[which(CD163_prop$phenotype == "CD163.mac"),]

pfi <- merge(CD163_prop, clinical, by.x="imageid", by.y="block name")



p <- ggplot(pfi, aes(x=PFIcat_12, y=freq, fill=PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue"))
p
pdf("E:/cellcycle/plots/boxplots/CD163_prop_from_all_macs_PFI.pdf", width=4, height = 4)
p
dev.off()

#add IBA1+CD163
data_ss$phenotype_mac <- data_ss$phenotype
data_ss[which(data_ss$phenotype == "IBA1.CD11.mac" |data_ss$phenotype == "IBA1.CD163.mac" |data_ss$phenotype == "CD163.mac"| data_ss$phenotype == "CD11c.myeloid"),"phenotype_mac"] <- "Macrophages"

data_ss$phenotype_CD163 <- data_ss$phenotype
data_ss[which(data_ss$phenotype == "IBA1.CD163.mac" |data_ss$phenotype == "CD163.mac"),'phenotype_CD163'] <- "CD163_pos_mac"

CD163_prop <- data_ss[which(data_ss$phenotype_mac == "Macrophages"),] %>%
  group_by(imageid, phenotype_CD163) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))

CD163_prop <- CD163_prop[which(CD163_prop$phenotype_CD163 == "CD163_pos_mac"),]

pfi <- merge(CD163_prop, clinical, by.x="imageid", by.y="block name")


p <- ggplot(pfi, aes(x=PFIcat_12, y=freq, fill=PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue")) + ylab("CD163+Macrophage proportion from myeloids")
p
pdf("E:/cellcycle/plots/boxplots/CD163_all_prop_from_all_myeloid_PFI.pdf", width=4, height = 4)
p
dev.off()


#download Giotto results

for (j in c("RCN7", "RCN9", "RCN6", "RCN5", "RCN1", "RCN2")){
giotto <- list()
for (i in unique(data_ss$imageid)){

  test <- read.csv(paste0("E:/cellcycle/Giotto/Giotto_phenotypes_",i,"_",j,".csv"))
#take genes: TIM.3
#cell_type: CD8Tcells
#int_cell_type are macrophages
#columns log2fc (target cell type interacting vs not interacting with interacting cell type)

  test <- test[which(test$genes == "TIM.3" & test$cell_type == "CD8Tcells"), c("int_cell_type", "log2fc")]

  test <- as.data.frame(t(test))
  colnames(test) <- test[1,]
  test <- test[-1,]
  test$imageid <- i
  giotto[[i]] <- test

}

#giotto_all <- do.call(rbind, giotto)
giotto_all <- bind_rows(giotto)

giotto_all$Tumor <- as.numeric(giotto_all$Tumor)
giotto_all$Stromal <- as.numeric(giotto_all$Stromal)
giotto_all$IBA1.CD163.mac <- as.numeric(giotto_all$IBA1.CD163.mac)
giotto_all$IBA1.CD11.mac <- as.numeric(giotto_all$IBA1.CD11.mac)
giotto_all$CD8Tcells <- as.numeric(giotto_all$CD8Tcells)
giotto_all$CD4Tcells <- as.numeric(giotto_all$CD4Tcells)
giotto_all$CD11c.myeloid <- as.numeric(giotto_all$CD11c.myeloid)
giotto_all$CD163.mac <- as.numeric(giotto_all$CD163.mac)
giotto_all$Other <- as.numeric(giotto_all$Other)
giotto_all$Unknown <- as.numeric(giotto_all$Unknown)

giotto_all <- merge(giotto_all, clinical, by.x="imageid", by.y="block name")
library(tidyr)
giotto_all <- pivot_longer(giotto_all, cols = 2:11, names_to='celltype', values_to='log2fc_TIM3')
pdf(paste0("E:/cellcycle/plots/boxplots/giotto/boxplots_", j, "_PFI.pdf"), width=20, height=4)
print(ggplot(giotto_all, aes(x=PFIcat_12, y=log2fc_TIM3, fill = PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue")) + ylab(paste0("log2fc TIM3 in CD8Tcells")) + facet_grid(. ~ celltype))
dev.off()

}



#download interactions
#per RCN


test <- read.csv(paste0("E:/cellcycle/interactions/zscore/","all","_interactions.csv"))
test <- test[,-c(1, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43)]
test$interactions <- paste0(test$phenotype, "_", test$neighbour_phenotype)
test <- test[,-c(1:2)]
test <- test[, c(21, 1:20)]
test <- t(test)
test <- as.data.frame(test)
colnames(test) <- test[1,]
test <- test[-1,]
test$imageid <- rownames(test)

test <- merge(test, clinical, by.x="imageid", by.y="block name")
#test$IBA1.CD163.mac_CD8Tcells <- as.numeric(test$IBA1.CD163.mac_CD8Tcells)
columns <- colnames(test)[2:73]
test[, columns] <- lapply(columns, function(x) as.numeric(test[[x]]))

keep <- c("imageid","CD11c.myeloid_CD8Tcells", "CD163.mac_CD8Tcells", "CD8Tcells_CD11c.myeloid", "CD8Tcells_CD163.mac", "CD8Tcells_IBA1.CD11.mac", "CD8Tcells_IBA1.CD163.mac","CD8Tcells_Tumor", "IBA1.CD11.mac_CD8Tcells","IBA1.CD163.mac_CD8Tcells", "Tumor_CD8Tcells", "PFIcat_12")

test <- test[, keep]

test <- pivot_longer(test, cols=2:11, names_to="interaction", values_to = "Interaction_zscore")

pdf(paste0("E:/cellcycle/plots/boxplots/interaction/zscore/PFI_interaction_","all",".pdf"), width=20, height=4)
print(ggplot(test, aes(x=PFIcat_12, y=Interaction_zscore, fill = PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue")) + ylab(paste0("Interaction"))) + facet_grid(.~ interaction)
dev.off()

#
#pscore

# for (i in c("CD11c", "CD163", "IBA1CD11c", "IBA1CD163")){
# test <- read.csv(paste0("E:/cellcycle/interactions/pscore/all_spatial_pscore_",i,".csv"))
# test <- merge(test, clinical, by.x="imageid", by.y="block name")
# 
# library(ggplot2)
# library(ggpubr)
# pdf(paste0("E:/cellcycle/plots/boxplots/interaction/pscore/PFI_interaction_pscore_all_",i,".pdf"), width=4, height=4)
# print(ggplot(test, aes(x=PFIcat_12, y=Proximity.Volume, fill = PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue")) )
# print(ggplot(test, aes(x=PFIcat_12, y=Proximity.Density, fill = PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue")) )
# dev.off()
# 
# }

#myelonet sizes
#see other R script

#myeloids from all

data_ss$phenotype_mac <- data_ss$phenotype
data_ss[which(data_ss$phenotype == "IBA1.CD11.mac" |data_ss$phenotype == "IBA1.CD163.mac" |data_ss$phenotype == "CD163.mac"),"phenotype_mac"] <- "Macrophages"

mac_prop <- data_ss %>%
  group_by(imageid, phenotype_mac) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))

mac_prop <- mac_prop[which(mac_prop$phenotype_mac == "Macrophages"),]

pfi <- merge(mac_prop, clinical, by.x="imageid", by.y="block name")


p <- ggplot(pfi, aes(x=PFIcat_12, y=freq, fill=PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue"))
p
pdf("E:/cellcycle/plots/boxplots/mac_prop_from_all_PFI.pdf", width=4, height = 4)
p
dev.off()

#CD163mac prop from all

data_ss$phenotype_mac <- data_ss$phenotype
data_ss[which(data_ss$phenotype == "IBA1.CD11.mac" |data_ss$phenotype == "IBA1.CD163.mac" |data_ss$phenotype == "CD163.mac"),"phenotype_mac"] <- "Macrophages"

mac_prop <- data_ss %>%
  group_by(imageid, phenotype) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))

mac_prop <- mac_prop[which(mac_prop$phenotype == "CD163.mac"),]

pfi <- merge(mac_prop, clinical, by.x="imageid", by.y="block name")



p <- ggplot(pfi, aes(x=PFIcat_12, y=freq, fill=PFIcat_12)) + geom_boxplot() + geom_point() + stat_compare_means(method="wilcox") + theme_bw() + scale_fill_manual(values=c("lightblue", "blue"))
p
pdf("E:/cellcycle/plots/boxplots/CD163mac_prop_from_all_PFI.pdf", width=4, height = 4)
p
dev.off()


#






