#CellID's connected to myelonets
#add them to the original data data


myelonet_macs <- read.csv("D:/sciset/delaunay/macrophage_myelonet_CellIDs.csv")
myelonet_IBA1CD163 <- read.csv("D:/sciset/delaunay/IBA1.CD163.macrophage_myelonet_CellIDs.csv")
myelonet_CD11c <- read.csv("D:/sciset/delaunay/CD11c_myelonet_CellIDs.csv")
myelonet_mixed_Immune <- read.csv("D:/sciset/delaunay/Immune_myelonet_CellIDs.csv")

#add to data
test <- read.csv("D:/sciset/data_20231125.csv")
test$ID_Sample_code <- paste0(test$ID, ".", test$Sample_code)

myelonet_macs$ID_Sample_code <- paste0(myelonet_macs$ID, ".", myelonet_macs$Sample_code)
myelonet_IBA1CD163$ID_Sample_code <- paste0(myelonet_IBA1CD163$ID, ".", myelonet_IBA1CD163$Sample_code)
myelonet_CD11c$ID_Sample_code <- paste0(myelonet_CD11c$ID, ".", myelonet_CD11c$Sample_code)
myelonet_mixed_Immune$ID_Sample_code <- paste0(myelonet_mixed_Immune$ID, ".", myelonet_mixed_Immune$Sample_code)

test$myelonet_macs <- NA
test[which(test$ID_Sample_code %in% myelonet_macs$ID_Sample_code), "myelonet_macs"] <- "myelonet_mac"

test$myelonet_IBA1CD163 <- NA
test[which(test$ID_Sample_code %in% myelonet_IBA1CD163$ID_Sample_code), "myelonet_IBA1CD163"] <- "myelonet_IBA1CD163"

test$myelonet_CD11c <- NA
test[which(test$ID_Sample_code %in% myelonet_CD11c$ID_Sample_code), "myelonet_CD11c"] <- "myelonet_CD11c"

test$myelonet_mixed_immune <- NA
test[which(test$ID_Sample_code %in% myelonet_mixed_Immune$ID_Sample_code), "myelonet_mixed_immune"] <- "myelonet_mixed_immune"

test$myelonet_any <- NA
test[which(test$myelonet_macs == "myelonet_mac"), "myelonet_any"] <- "myelonet_any"
test[which(test$myelonet_IBA1CD163 == "myelonet_IBA1CD163"), "myelonet_any"] <- "myelonet_any"
test[which(test$myelonet_CD11c == "myelonet_CD11c"), "myelonet_any"] <- "myelonet_any"
test[which(test$myelonet_mixed_immune == "myelonet_mixed_immune"), "myelonet_any"] <- "myelonet_any"


#export data to be able to load into scimap

write.csv(test, "D:/sciset/data_with_myelonets.csv", row.names = F)

test$myelonet_any_celltype <- paste0(test$GlobalCellType2, '.',test$myelonet_any)

#TIM3 expression in myelonets vs elsewhere

clinical <- read.csv("L:/ltdk_farkkila/Data/Sciset/clinical_data/clinical_data.csv")

library(ggplot2)

test <- merge(test, clinical[, c("Patient.x", "Sample_code", "Stage", "HRD_status")], by="Sample_code")

#entä ota vaan paired

ggplot(test[which(test$GlobalCellType2 == "CD8.T.cells" & test$Paired == TRUE),], aes(x=myelonet_any, y=TIM3, fill=Stage)) + geom_violin(position = position_dodge(0.9)) + geom_boxplot(width=0.2, position = position_dodge(0.9))

ggplot(test[which(test$GlobalCellType2 == "CD8.T.cells"& test$Paired == TRUE),], aes(x=myelonet_macs, y=TIM3, fill=Stage)) + geom_violin(position = position_dodge(0.9)) + geom_boxplot(width=0.2, position = position_dodge(0.9))

ggplot(test[which(test$GlobalCellType2 == "CD8.T.cells"& test$Paired == TRUE),], aes(x=myelonet_IBA1CD163, y=TIM3, fill=Stage)) + geom_violin(position = position_dodge(0.9)) + geom_boxplot(width=0.2, position = position_dodge(0.9))

ggplot(test[which(test$GlobalCellType2 == "CD8.T.cells"& test$Paired == TRUE),], aes(x=myelonet_CD11c, y=TIM3, fill=Stage)) + geom_violin(position = position_dodge(0.9)) + geom_boxplot(width=0.2, position = position_dodge(0.9))

ggplot(test[which(test$GlobalCellType2 == "CD8.T.cells"& test$Paired == TRUE),], aes(x=myelonet_mixed_immune, y=TIM3, fill=Stage)) + geom_violin(position = position_dodge(0.9)) + geom_boxplot(width=0.2, position = position_dodge(0.9))



#next
#how big part the myeloid cells are in the myelonets?

#in myelonets / elsewhere

test$myelonet_all <- NA
test[which(test$myelonet_macs == "myelonet_mac"), "myelonet_all"] <- "myelonet_RCN14"
test[which(test$myelonet_IBA1CD163 == "myelonet_IBA1CD163"), "myelonet_all"] <- "myelonet_RCN15"
test[which(test$myelonet_CD11c == "myelonet_CD11c"), "myelonet_all"] <- "myelonet_RCN16"
test[which(test$myelonet_mixed_immune == "myelonet_mixed_immune"), "myelonet_all"] <- "myelonet_RCN18"

#x myeloid cell types
#y location

macs <- test[which(test$GlobalCellType2 == "IBA1.CD163.Macrophages" | test$GlobalCellType2 == "IBA1.CD11c.Macrophages" | test$GlobalCellType2 == "CD163.Macrophages" | test$GlobalCellType2 == "CD11c.myeloid"),]
p <- ggplot(macs, 
       aes(x=factor(GlobalCellType2, levels=c("CD163.Macrophages", "CD11c.myeloid", "IBA1.CD163.Macrophages", "IBA1.CD11c.Macrophages")), fill=myelonet_all)) + geom_bar(position="fill", stat="count") + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("myeloid cell type")
  
pdf("E:/Sciset/scimap/plots/revision/myeloids_inside_myelonets.pdf", width=3, height=5)
p
dev.off()
#3e92ccff , #d8315bff,#0a2463ff, #fffaffff, seashell3

#then stacked

p <- ggplot(macs, aes(x=factor(GlobalCellType2, levels=c("CD163.Macrophages","IBA1.CD11c.Macrophages", "CD11c.myeloid", "IBA1.CD163.Macrophages")), fill=myelonet_all)) +
  geom_bar() + scale_fill_manual(values=c(  '#3e92ccff' ,  '#d8315bff', '#0a2463ff', '#fffaffff', 'seashell3')) +theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("myeloid cell type")
p

pdf("E:/Sciset/scimap/plots/revision/myeloids_inside_myelonets_counts.pdf", width=3, height=5)
p
dev.off()



#then violin plot

library(ggplot2)
library(ggpubr)

test_all$myelonet_merged <- test_all$myelonet_all
test_all[which(test_all$myelonet_merged == "myelonet_RCN14" | test_all$myelonet_merged == "myelonet_RCN15"| test_all$myelonet_merged== "myelonet_RCN16" | test_all$myelonet_merged=="myelonet_RCN18"), "myelonet_merged"] <- "myelonet"

my_comparisons <- list(c("myelonet", "no_interactions"), c("no_interactions", "non_myelonet"), c("myelonet", "non_myelonet"))

p <- ggplot(test_all, aes(x=factor(myelonet_merged, levels=c("myelonet", "non_myelonet", "no_interactions")), y=TIM3)) + 
  geom_violin(aes(fill=myelonet_merged)) + geom_boxplot(width=0.2, outlier.shape = NA, aes(fill=myelonet_merged)) + 
  theme_bw() + scale_fill_manual(values=c("orange", "pink", "lightblue")) + 
  stat_compare_means(comparisons=my_comparisons)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + xlab("location of interaction")
p

pdf("D:/sciset/revision/TIM3_in_CD8T_in_myelonets_and_elsewhere_all.pdf", width=3, height=4)
p
dev.off()

# primary ja interval
test_primary[which(test_primary$spatial_pscore_IBA1.CD163.Macrophages == "other" & test_primary$spatial_pscore_CD163.Macrophages == "other" & test_primary$spatial_pscore_CD11c.myeloid == "other" & test_primary$spatial_pscore_IBA1.CD11c.Macrophages == "other"), "myelonet_all"] <- "no_interactions"
test_primary$myelonet_merged <- test_primary$myelonet_all
test_primary[which(test_primary$myelonet_merged == "myelonet_RCN14" | test_primary$myelonet_merged == "myelonet_RCN15"| test_primary$myelonet_merged== "myelonet_RCN16" | test_primary$myelonet_merged=="myelonet_RCN18"), "myelonet_merged"] <- "myelonet"
p <- ggplot(test_primary, aes(x=factor(myelonet_merged, levels=c("myelonet", "non_myelonet", "no_interactions")), y=TIM3)) + 
  geom_violin(aes(fill=myelonet_merged)) + geom_boxplot(width=0.2, outlier.shape = NA, aes(fill=myelonet_merged)) + 
  theme_bw() + scale_fill_manual(values=c("orange", "pink", "lightblue")) + 
  stat_compare_means(comparisons=my_comparisons)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + xlab("location of interaction")
p

pdf("D:/sciset/revision/TIM3_in_CD8T_in_myelonets_and_elsewhere_chemonaive.pdf", width=3, height=4)
p
dev.off()

#interval

test_interval[which(test_interval$spatial_pscore_IBA1.CD163.Macrophages == "other" & test_interval$spatial_pscore_CD163.Macrophages == "other" & test_interval$spatial_pscore_CD11c.myeloid == "other" & test_interval$spatial_pscore_IBA1.CD11c.Macrophages == "other"), "myelonet_all"] <- "no_interactions"
test_interval$myelonet_merged <- test_interval$myelonet_all
test_interval[which(test_interval$myelonet_merged == "myelonet_RCN14" | test_interval$myelonet_merged == "myelonet_RCN15"| test_interval$myelonet_merged== "myelonet_RCN16" | test_interval$myelonet_merged=="myelonet_RCN18"), "myelonet_merged"] <- "myelonet"
p <- ggplot(test_interval, aes(x=factor(myelonet_merged, levels=c("myelonet", "non_myelonet", "no_interactions")), y=TIM3)) + 
  geom_violin(aes(fill=myelonet_merged)) + geom_boxplot(width=0.2, outlier.shape = NA, aes(fill=myelonet_merged)) + 
  theme_bw() + scale_fill_manual(values=c("orange", "pink", "lightblue")) + 
  stat_compare_means(comparisons=my_comparisons)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + xlab("location of interaction")
p

pdf("D:/sciset/revision/TIM3_in_CD8T_in_myelonets_and_elsewhere_interval.pdf", width=3, height=4)
p
dev.off()

#sitten dotplot post ve pre

test_interval$stage <- "IDS"
test_interval <- test_interval[, c(1, 2, 3, 4, 5, 6, 7, 9, 8)]

test_interval_dotplot <- test_interval[-which(test_interval$myelonet_merged == "no_interactions"),]
test_primary_dotplot <- test_primary[-which(test_primary$myelonet_merged == "no_interactions"),]

test_all_dotplot <- rbind(test_interval_dotplot, test_primary_dotplot)

test_all_dotplot <- test_all_dotplot[, c(1, 3, 4, 5, 6, 7, 8)]

library(tidyr)
data_bulk_long <- pivot_longer(test_all_dotplot, cols=spatial_pscore_IBA1.CD163.Macrophages:spatial_pscore_CD11c.myeloid, values_to='interaction', names_to='pair')
metacluster_percentages <- data_bulk_long

metacluster_percentages <- metacluster_percentages[, -4]

library(dplyr)
test <- metacluster_percentages %>%
  group_by(myelonet_all, interaction) %>%
  summarise(
    Stage_pval = wilcox.test(TIM3[stage == 'IDS'], 
                             TIM3[stage == 'chemo-naive'], paired=FALSE)$p.value)


df2 <- metacluster_percentages %>%
  group_by(myelonet_all, interaction, stage) %>%
  summarize(TIM3 = mean(TIM3))

df2

df3 <- df2 %>%
  pivot_wider(names_from = 'stage', values_from = 'TIM3')
df3

df3 <- df3 %>%
  mutate(foldChange = log2(`IDS`/`chemo-naive`))

df2_bulk <- merge(test, df3)


df2_bulk <- df2_bulk[,-c(4, 5)]

colnames(df2_bulk)[3] <- "pvalue"
colnames(df2_bulk)[4] <- "log_foldchange"

truncate.df <- function(df, na.cutoff=0.1, na.var="pvalue", na.var.boundary=1e-500, range.lims=c(-2, 2), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(df2_bulk) 
#to.plot$comparison <- "bulk"
#row_order <- rev(c("tcycif", "scRNAseq", "bulkRNAseq"))
to.plot <- to.plot[-which(to.plot$interaction == "other"),]
to.plot[which(to.plot$pvalue == 0), "pvalue"] <- 1e-100

p <- ggplot(to.plot,  aes(x=factor(myelonet_all), y=factor(interaction))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log2) IDS vs chemo-naive") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90",  breaks=waiver()) +
  scale_size_area("p-values", trans="log10",max_size = 1.5, breaks=c(1e-50, 1e-20, 0.05), limits=c(0,0.05)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=90, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)

pdf("D:/sciset/revision/dotplot_TIM3_in_RCN_myelonets.pdf", width=8, height=5)
p
dev.off()





#where do myeloid cells reside in percentages

proportion <- macs[which(macs$GlobalCellType2 == "CD163.Macrophages"),] %>% group_by(myelonet_any) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))
#17,1
proportion <- macs[which(macs$GlobalCellType2 == "IBA1.CD163.Macrophages"),] %>% group_by(myelonet_any) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))
#42,2
proportion <- macs[which(macs$GlobalCellType2 == "IBA1.CD11c.Macrophages"),] %>% group_by(myelonet_any) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))
#44,6%
proportion <- macs[which(macs$GlobalCellType2 == "CD11c.myeloid"),] %>% group_by(myelonet_any) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))
#30,6

#then interactions
cd8t <- test[which(test$GlobalCellType2 == "CD8.T.cells"),]

proportion <- cd8t[which(cd8t$spatial_pscore_IBA1.CD163.Macrophages == "IBA1.CD163.Macrophages_CD8.T.cells"),] %>% group_by(myelonet_any) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))
proportion
#23,6

proportion <- cd8t[which(cd8t$spatial_pscore_IBA1.CD11c.Macrophages == "IBA1.CD11c.Macrophages_CD8.T.cells"),] %>% group_by(myelonet_any) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))
proportion
#27,5

proportion <- cd8t[which(cd8t$spatial_pscore_CD163.Macrophages == "CD163.Macrophages_CD8.T.cells"),] %>% group_by(myelonet_any) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))
proportion
#18,8

proportion <- cd8t[which(cd8t$spatial_pscore_CD11c.myeloid == "CD11c.myeloid_CD8.T.cells"),] %>% group_by(myelonet_any) %>% summarize(n=n()) %>% mutate(freq=n/sum(n))
proportion
#23,7

#



















