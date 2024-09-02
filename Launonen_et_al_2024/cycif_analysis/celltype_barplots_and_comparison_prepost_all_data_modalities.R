#Figures 1f-g and supplementary figure 1 f, i and j

#plot all data types with matching cell types


library(dplyr)
library(tidyr)
library(RColorBrewer)
library(circlize)
library(viridis)
library(ComplexHeatmap)
library(ggplot2)
#cycif

data_sciset <- read.csv("E:/sciset/data_20231125.csv")

data_sciset$Global_celltype_tumor_stroma_merged <- data_sciset$GlobalCellType2
data_sciset[which(data_sciset$GlobalCellType2 == "EMT" | data_sciset$GlobalCellType2 == "Epithelial" | data_sciset$GlobalCellType2 == "Proliferating.EMT" | data_sciset$GlobalCellType2 == "Proliferating.epithelial"), "Global_celltype_tumor_stroma_merged"] <- "tumor"
data_sciset[which(data_sciset$GlobalCellType2 == "Fibroblast" | data_sciset$GlobalCellType2 == "SMA.CD31.positive.cell" | data_sciset$GlobalCellType2 == "Myofibroblast" | data_sciset$GlobalCellType2 == "Endothelial.cell" | data_sciset$GlobalCellType2 == "Desmin.positive.cell" | data_sciset$GlobalCellType2 == "SMA.Desmin.positive.cell" ), "Global_celltype_tumor_stroma_merged"] <- "stroma"

data_sciset[which(data_sciset$GlobalCellType2 == "CD11c.myeloid" | data_sciset$GlobalCellType2 == "IBA1.CD11c.Macrophages" | data_sciset$GlobalCellType2 == "IBA1.CD163.Macrophages" | data_sciset$GlobalCellType2 == "CD163.Macrophages"), "Global_celltype_tumor_stroma_merged"] <- "Myeloid cells"

clinical <- read.csv("D:/Sciset/clinical_data.csv")
colnames(clinical)[2] <- "Patient"
clinical <- unique(clinical[, c("Sample_code", "Stage", "Patient")])

TSEI_percentages <- data_sciset %>% group_by(Sample_code, Global_celltype_tumor_stroma_merged) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))%>%
  ungroup() %>%
  complete(Sample_code, Global_celltype_tumor_stroma_merged,
           fill = list(n = 0, proportion = 0))

TSEI_percentages <- merge(TSEI_percentages, clinical[, c("Sample_code", "Stage", "Patient", "Paired")], by="Sample_code")

metacluster_percentages <- as.data.frame(TSEI_percentages)


metacluster_percentages <- metacluster_percentages[which(metacluster_percentages$Paired == TRUE),]

test <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    Stage_pval = wilcox.test(proportion[Stage == 'interval'], 
                             proportion[Stage == 'primary'], paired=TRUE)$p.value)


df2 <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged, Stage) %>%
  summarize(proportion = mean(proportion))

df2

df3 <- df2 %>%
  pivot_wider(names_from = 'Stage', values_from = 'proportion')
df3

df3 <- df3 %>%
  mutate(foldChange = log2(`interval`/`primary`))

#
df2_sciset <- merge(test, df3)

#scRNAseq

data_sc <- read.table("D:/Sciset/scimap/Figures/for_figure_1/scrnaseq_cell_fractions_celltype_high_v2.txt", header=TRUE, sep="\t")

data_sc <- data_sc[, c("Other_immune","Myeloid","B_cell","CD8_T_cell","CD4_T_cell" ,"FOXP3_CD4_Treg","Tumor", "Stroma")]
colnames(data_sc) <- c("Other_immune","Myeloid","B_cell","CD8.T.cells","CD4.T.cells" ,"FOXP3.CD4.Tregs","tumor", "stroma")
colnames(data_sc)[2] <- "Myeloid cells"
data_sc$Sample_code <- rownames(data_sc)
data_sc$Paired <- TRUE

data_sc <- data_sc[which(data_sc$Paired == TRUE),]

data_sc$Stage <- "primary"
data_sc[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43), "Stage"] <- "interval"

data_sc_long <- pivot_longer(data_sc, cols=Other_immune:stroma, values_to='proportion', names_to='Global_celltype_tumor_stroma_merged')
metacluster_percentages <- data_sc_long

test <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    Stage_pval = wilcox.test(proportion[Stage == 'interval'], 
                             proportion[Stage == 'primary'], paired=TRUE)$p.value)


df2 <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged, Stage) %>%
  summarize(proportion = mean(proportion, na.rm=T))

df2

df3 <- df2 %>%
  pivot_wider(names_from = 'Stage', values_from = 'proportion')
df3

df3 <- df3 %>%
  mutate(foldChange = log2(`interval`/`primary`))

df2_sc <- merge(test, df3)


#bulkRNAseq

data_bulk <- read.table("E:/sciset/bulk/preds_50_samples.tsv")

data_bulk <- t(data_bulk)
data_bulk <- as.data.frame(data_bulk)

data_bulk$Bcell <- data_bulk$B_cells
data_bulk$CD8.T.cells <- data_bulk$CD8_T_cells
data_bulk$CD4.T.cells <- data_bulk$T_helpers
data_bulk$stroma <- data_bulk$Fibroblasts + data_bulk$Endothelium
data_bulk$myeloids <- data_bulk$Macrophages + data_bulk$Monocytes
data_bulk$FOXP3.CD4.Tregs <- data_bulk$Tregs
data_bulk$tumor <- data_bulk$Other

data_bulk <- data_bulk[, c("Bcell","myeloids", "CD8.T.cells","CD4.T.cells" ,"FOXP3.CD4.Tregs","tumor", "stroma")]
colnames(data_bulk)[2] <- "Myeloid cells"
data_bulk$Sample_code <- rownames(data_bulk)
data_bulk$Paired <- TRUE

data_bulk <- data_bulk[which(data_bulk$Paired == TRUE),]

data_bulk$Stage <- "primary"
data_bulk[c(1:25), "Stage"] <- "interval"

data_bulk_long <- pivot_longer(data_bulk, cols=Bcell:stroma, values_to='proportion', names_to='Global_celltype_tumor_stroma_merged')
metacluster_percentages <- data_bulk_long

test <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    Stage_pval = wilcox.test(proportion[Stage == 'interval'], 
                             proportion[Stage == 'primary'], paired=TRUE)$p.value)


df2 <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged, Stage) %>%
  summarize(proportion = mean(proportion))

df2

df3 <- df2 %>%
  pivot_wider(names_from = 'Stage', values_from = 'proportion')
df3

df3 <- df3 %>%
  mutate(foldChange = log2(`interval`/`primary`))

df2_bulk <- merge(test, df3)



#then merge data 

df2_sc <- df2_sc[-which(df2_sc$Global_celltype_tumor_stroma_merged == "Other_immune"),]
df2_sc[which(df2_sc$Global_celltype_tumor_stroma_merged == "B_cell"),"Global_celltype_tumor_stroma_merged"] <- "Bcell"


df2_sciset$comparison <- "tcycif"
df2_sc$comparison <- "scRNAseq"
df2_bulk$comparison <- "bulkRNAseq"


df2 <- rbind(df2_sciset[, c(1, 2, 5, 6)], df2_sc[, c(1, 2, 5, 6)], df2_bulk[, c(1, 2, 5, 6)])
colnames(df2)[2] <- "pvalue"
colnames(df2)[3] <- "log_foldchange"

truncate.df <- function(df, na.cutoff=0.05, na.var="pvalue", na.var.boundary=1e-50, range.lims=c(-2, 2), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(df2) 
row_order <- rev(c("tcycif", "scRNAseq", "bulkRNAseq"))
p <- ggplot(to.plot,  aes(x=factor(Global_celltype_tumor_stroma_merged), y=factor(comparison, level=row_order))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90",  breaks=waiver()) +
  scale_size_area("p-values", trans="log10",max_size = 5.5, breaks=c(1e-6, 1e-2, 0.05), limits=c(1e-10,0.05)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)


pdf("E:/sciset/bulk/comparison_pre_post_all_data_modalities_20231030.pdf", width = 10, height = 2.5)
print(p)
dev.off()


#same but as a heatmap


df2_wide <- pivot_wider(df2[, c(1, 3, 4)], names_from = "Global_celltype_tumor_stroma_merged", values_from="log_foldchange")
df2_wide <- as.data.frame(df2_wide)
rownames(df2_wide) <- df2_wide$comparison


col_c <- colorRamp2(c(-1.3, 0, 1.3), hcl_palette = "Blue-Red 3")
df2_wide <- df2_wide[, -c(7)]
hmap <- Heatmap(as.matrix(df2_wide[, -1]), name="log2fc interval vs primary", cluster_rows = F,
                row_names_gp = gpar(fontsize = 8),clustering_method_rows ="ward.D2",column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),cluster_columns = F,
                left_annotation = NULL, row_dend_width = unit(1.3, "cm"),show_row_names = T, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                show_row_dend=T, show_column_dend=T, column_title = "subtypes", row_title = "Data modality",
                border="black", width = unit(6.3, "cm"),height = unit(3, "cm"),
                rect_gp = gpar(col = "white", lwd = 2), col=col_c)


draw(hmap)

pdf("E:/sciset/bulk/comparison_pre_post_all_data_modalities_heatmap_20240123.pdf", width = 10, height = 2.5)
draw(hmap)
dev.off()


#then correlations

#cycif proportions
TSEI_percentages <- data_sciset %>% group_by(Sample_code, Global_celltype_tumor_stroma_merged) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))%>%
  ungroup() %>%
  complete(Sample_code, Global_celltype_tumor_stroma_merged,
           fill = list(n = 0, proportion = 0))


TSEI_percentages <- merge(TSEI_percentages, clinical_sciset[, c("Sample_code", "Stage", "Patient", "Paired")], by="Sample_code")

metacluster_percentages <- as.data.frame(TSEI_percentages)


#scRNAseq proportions


data_sc <- read.table("D:/Sciset/scimap/Figures/for_figure_1/scrnaseq_cell_fractions_celltype_high_v3.txt", header=TRUE, sep="\t")
colnames(data_sc) <- c("Bcell", "CD4.T.cells","CD8.T.cells","FOXP3.CD4.Tregs", "Myeloid cells", "Other","stroma","tumor" )
data_sc$Sample_code <- rownames(data_sc)
data_sc_long <- pivot_longer(data_sc, cols=Bcell:tumor, values_to='proportion', names_to='Global_celltype_tumor_stroma_merged')
metacluster_percentages2 <- data_sc_long

#then bulkRNAseq

data_bulk <- read.csv("E:/sciset/bulk/predictions-27-samples.tsv", sep="\t")
colnames(data_bulk)
data_bulk <- as.data.frame(data_bulk)
rownames(data_bulk) <- data_bulk$X
data_bulk <- data_bulk[,-1]
data_bulk <- t(data_bulk)
data_bulk <- as.data.frame(data_bulk)



data_bulk$Bcell <- data_bulk$B_cells
data_bulk$CD8.T.cells <- data_bulk$CD8_T_cells
data_bulk$CD4.T.cells <- data_bulk$T_helpers
data_bulk$Endothelial.cell <- data_bulk$Endothelium
data_bulk$myeloids <- data_bulk$Monocytes + data_bulk$Macrophages
data_bulk$FOXP3.CD4.Tregs <- data_bulk$Tregs
data_bulk$other_immune <- data_bulk$NK_cells + data_bulk$Neutrophils
data_bulk$tumor <- data_bulk$Other
data_bulk$stroma <- data_bulk$Fibroblasts

data_bulk <- data_bulk[, c("Bcell","myeloids", "CD8.T.cells","CD4.T.cells" ,"Endothelial.cell","FOXP3.CD4.Tregs","tumor", "stroma")]
colnames(data_bulk)[2] <- "Myeloid cells"
data_bulk$Sample_code <- rownames(data_bulk)
data_bulk$stroma <- data_bulk$Endothelial.cell + data_bulk$stroma
data_bulk <- data_bulk[,-5]

data_bulk_long <- pivot_longer(data_bulk, cols=Bcell:stroma, values_to='proportion', names_to='Global_celltype_tumor_stroma_merged')
metacluster_percentages3 <- data_bulk_long


#keep first 8 characters of Sample_code

metacluster_percentages$Sample_code <- substr(metacluster_percentages$Sample_code, start = 1, stop = 8)
metacluster_percentages2$Sample_code <- substr(metacluster_percentages2$Sample_code, start = 1, stop = 8)
metacluster_percentages3$Sample_code <- substr(metacluster_percentages3$Sample_code, start = 1, stop = 8)

metacluster_percentages2[which(metacluster_percentages2$Global_celltype_tumor_stroma_merged == "B_cell"), "Global_celltype_tumor_stroma_merged"] <- "Bcell"

cycif_scRNAseq <- merge(metacluster_percentages[, c("Sample_code", "Global_celltype_tumor_stroma_merged", "proportion")], metacluster_percentages2, by=c("Sample_code", "Global_celltype_tumor_stroma_merged"))

#then correlations and coefficients
#estimate and p.value
test <- cycif_scRNAseq %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    cor_pval = cor.test(proportion.x, 
                             proportion.y, method="spearman")$p.value)



#coeffifients

cor <- cycif_scRNAseq %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    cor = cor.test(proportion.x, 
                        proportion.y, method="spearman")$estimate)



#same cycif ja bulk

cycif_bulk <- merge(metacluster_percentages[, c("Sample_code", "Global_celltype_tumor_stroma_merged", "proportion")], metacluster_percentages3, by=c("Sample_code", "Global_celltype_tumor_stroma_merged"))

#sit correlations and coefficients
#estimate and p.value
test2 <- cycif_bulk %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    cor_pval = cor.test(proportion.x, 
                        proportion.y, method="spearman")$p.value)



#coeffifients

cor2 <- cycif_bulk %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    cor = cor.test(proportion.x, 
                   proportion.y, method="spearman")$estimate)

test <- as.data.frame(test)
cor <- as.data.frame(cor)

test2 <- as.data.frame(test2)
cor2 <- as.data.frame(cor2)


#lis?? comparison

test <- cbind(test, cor)
test2 <- cbind(test2, cor2)

test <- test[, -c(3)]
test2 <- test2[,-c(3)]

test$comparison <- "tcycif_scRNAseq"
test2$comparison <- "tcycif_bulkRNAseq"

test3 <- rbind(test, test2)

colnames(test3)[2] <- "pvalue"

row_order <- c("tcycif_scRNAseq", "tcycif_bulkRNAseq")

to.plot <- test3

to.plot <- to.plot[-which(to.plot$Global_celltype_tumor_stroma_merged == "Other"),]

p <- ggplot(to.plot,  aes(x=factor(Global_celltype_tumor_stroma_merged), y=factor(comparison, level=row_order))) +
  geom_point(aes(color=cor, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="correlation coef") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90",  breaks=waiver()) +
  scale_size_area("p-values", trans="log10",max_size = 2) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))
p

pdf("D:/Sciset/scimap/Figures/for_figure_1/dotplot_correlations_all_20240129.pdf", height = 2.5, width=5)
p
dev.off()



#then barplots from cell types

#metacluster percentages1 2 3 have proportions

#metacluster_percentages - cycif
#metacluster_percentages2 - scRNAseq
#metacluster_percentages3 - bulk

# data_sc <- read.table("D:/Sciset/scimap/Figures/for_figure_1/scrnaseq_cell_fractions_celltype_high.txt", header=TRUE, sep="\t")
# data_sc <- data_sc[, c("Other_immune","Myeloid","B_cell", "CD8_T_cell","CD4_T_cell" ,"FOXP3_CD4_Treg","Tumor", "Stroma")]
# colnames(data_sc) <- c("Other_immune","Myeloid","B_cell","CD8.T.cells","CD4.T.cells" ,"FOXP3.CD4.Tregs","tumor", "stroma")
# 
# colnames(data_sc)[2] <- "Myeloid cells"
# data_sc$Sample_code <- rownames(data_sc)
# 
# 
# data_sc_long <- pivot_longer(data_sc, cols=Other_immune:stroma, values_to='proportion', names_to='Global_celltype_tumor_stroma_merged')
# metacluster_percentages2 <- data_sc_long
metacluster_percentages$Sample_code <- substr(metacluster_percentages$Sample_code, start = 1, stop = 8)
data_cycif_wide <- metacluster_percentages

data_cycif_wide <- pivot_wider(data_cycif_wide, names_from="Global_celltype_tumor_stroma_merged", values_from="proportion")


#data_sc_wide <- metacluster_percentages2[-which(metacluster_percentages2$Global_celltype_tumor_stroma_merged == "Other"),]
data_sc_wide <- metacluster_percentages2
data_sc_wide[which(data_sc_wide$Global_celltype_tumor_stroma_merged == "Other"), "Global_celltype_tumor_stroma_merged"] <- "other_immune"
data_sc_wide <- pivot_wider(data_sc_wide, names_from="Global_celltype_tumor_stroma_merged", values_from="proportion")



#tumor cell percentage ordering
# plot_order <- c("H158_pPer", "H098_pOme", "H102_pMes", "H142_pPer", "H110_pOvaR", 
#                 "H142_iOme", "H144_pPer", "H166_pPer", "H114_pOme", "H173_pOme", 
#                 "H087_iOvaR", "H095_pOvaR", "H092_pPer", "H144_iOme", "H114_iOvaR", 
#                 "H103_pPer", "H173_iOme", "H116_iOme")
# plot_order <- as.data.frame(plot_order)
# colnames(plot_order)[1] <- "Sample_code"
#plot_order$Sample_code <- substr(plot_order$Sample_code, start = 1, stop = 8)

data_sc_wide$Sample_code <- substr(data_sc_wide$Sample_code, start = 1, stop = 6)
plot_order$Sample_code <- substr(plot_order$Sample_code, start = 1, stop = 6)


#data_sc_wide is long
data_sc_wide <- merge(data_sc_wide, plot_order, by="Sample_code")



df <- data_sc_wide

df <- pivot_longer(df, cols=Bcell:tumor, names_to = "GlobalCellType3")

df[which(df$GlobalCellType3 == "Myeloid cells"), "GlobalCellType3"] <- "myeloid"

pdf("D:/Sciset/scimap/plots/barplots/barplot_scRNAseq_celltypes.pdf", width=7, height=4)
ggplot(df, aes(x = factor(Sample_code, levels=plot_order$Sample_code), y = value, fill = GlobalCellType3))+
  geom_bar(stat = "identity", colour="black", position = "stack") + scale_fill_manual(values=group.colors)+ xlab("Sample") + ylab("Proportion")+ 
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_discrete(expand = c(0,0)) + theme(legend.position="right")

dev.off()

group.colors <- c('tumor' = "#627D87",'CD8.T.cells' = "#ff5a5f", 'CD4.T.cells'="#57cc99", 'B_cell'="#0b3954", 'FOXP3.CD4.Tregs'="#bfd7ea", 'Myeloid cells'="#996fd6",'stroma'="#d6ccc2", "Other_immune" = "gray69")

#RNAseq barplot matched samples
data_bulk <- read.csv("E:/sciset/bulk/predictions-27-samples.tsv", sep="\t")
colnames(data_bulk)
data_bulk <- as.data.frame(data_bulk)
rownames(data_bulk) <- data_bulk$X
data_bulk <- data_bulk[,-1]
data_bulk <- t(data_bulk)
data_bulk <- as.data.frame(data_bulk)


data_bulk$Bcell <- data_bulk$B_cells
data_bulk$CD8.T.cells <- data_bulk$CD8_T_cells
data_bulk$CD4.T.cells <- data_bulk$T_helpers
data_bulk$stroma <- data_bulk$Fibroblasts + data_bulk$Endothelium
data_bulk$myeloids <- data_bulk$Monocytes + data_bulk$Macrophages
data_bulk$FOXP3.CD4.Tregs <- data_bulk$Tregs
data_bulk$other_immune <- data_bulk$NK_cells + data_bulk$Neutrophils
data_bulk$tumor <- data_bulk$Other

data_bulk <- data_bulk[, c("Bcell","myeloids", "CD8.T.cells","CD4.T.cells" ,"other_immune","FOXP3.CD4.Tregs","tumor", "stroma")]
colnames(data_bulk)[2] <- "myeloid"
data_bulk$Sample_code <- rownames(data_bulk)

data_bulk_long <- pivot_longer(data_bulk, cols=Bcell:stroma, values_to='proportion', names_to='Global_celltype_tumor_stroma_merged')
metacluster_percentages3 <- data_bulk_long

metacluster_percentages3 <- metacluster_percentages3[!metacluster_percentages3$Sample_code %in% c("H092_pPer1_RNA1_Data2", "H095_pOvaR1_RNA2", "H116_iOme1_RNA1", "H142_iOme2_RNA1_Data2.1", "H142_iOme2_RNA1_Data2","H144_iOme1_RNA1_Data2"),]


metacluster_percentages3$Sample_code <- substr(metacluster_percentages3$Sample_code, start = 1, stop = 6)



metacluster_percentages3 <- as.data.frame(metacluster_percentages3)

plot_order$Sample_code <- substr(plot_order$Sample_code, start = 1, stop = 6)
metacluster_percentages3$Sample_code <- substr(metacluster_percentages3$Sample_code, start = 1, stop = 6)

metacluster_percentages3 <- metacluster_percentages3[metacluster_percentages3$Sample_code %in% plot_order$Sample_code,]

group.colors <- c('tumor' = "#627D87",'CD8.T.cells' = "#ff5a5f", 'CD4.T.cells'="#57cc99", 'Bcell'="#0b3954", 'FOXP3.CD4.Tregs'="#bfd7ea", 'myeloid'="#996fd6",'stroma'="#d6ccc2", "other_immune" = "gray69")
pdf("D:/Sciset/scimap/plots/barplots/barplot_RNAseq_celltypes.pdf", width=7, height=4)
ggplot(metacluster_percentages3, aes(x = factor(Sample_code, levels=plot_order$Sample_code), y = proportion, fill = Global_celltype_tumor_stroma_merged))+
  geom_bar(stat = "identity", colour="black", position = "stack") + scale_fill_manual(values=group.colors)+ xlab("Sample") + ylab("Proportion")+ 
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_discrete(expand = c(0,0)) + theme(legend.position="right")

dev.off()

group.colors <- c('Epithelial' = "#627D87", 'EMT' = "#467181", 'Proliferating.epithelial' ="#747483", 'Proliferating.EMT' = "#235A6F", 
                  'CD8.T.cells' = "#ff5a5f", 'CD4.T.cells'="#57cc99", 'Bcell'="#0b3954", 'FOXP3.CD4.Tregs'="#bfd7ea", 'CD11c.myeloid'="#996fd6", 
                  'IBA1.Macrophages'="#a786db", 'IBA1.CD11c.Macrophages'="#b59ce0", 'IBA1.CD163.Macrophages'="#c2b3e5", 'CD163.Macrophages'="#d0c9ea",
                  'SMA.Desmin.positive.cell'="#edede9", 'Desmin.positive.cell'="#d6ccc2", 'SMA.CD31.positive.cell'="#f5ebe0",
                  'Fibroblast'="#e3d5ca", 'Myofibroblast'="#d5bdaf", 'Endothelial.cell'="#9E8E6E")


#then annotaatiot
colours<- list("PFI at outcome update when no prog_Days"=col_fun2,"PFI after Primary therapy if Prog_Days"=col_fun1,"CUD.Progression"=col_prog,"CUD.Treatment.strategy"=col_tr,"CRS Omental" = col_CRS, "Primary.therapy.outcome" = col_response,"tcycif"= col_tissue,"tcycif_GeoMx" = col_tissue,"bulkRNAseq" = col_tissue,"scRNAseq" = col_tissue,"Stage"=c("primary"="royalblue","interval"="#e8c547"), "HRD_status"=c("HRD" = "#3399CC", "HRP"="red3"),"SBS3"=c("HRD" = "#3399CC", "HRP"="red3") )

Rowann <- HeatmapAnnotation(df=Annotations[, c(1:14)], which="column", col=colours, 
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(1, "mm"), na_col = "grey80")
col_c <- colorRamp2(c(-1.5, -0.5, 0.5, 1.5), c( "#0000FF","#AAAAFF", "#F9AAAA", "#EE0000"))
hmap <- Heatmap(t(as.matrix(Annotations[, c("PFI after Primary therapy if Prog_Days")])), cluster_columns = F, top_annotation = Rowann, width = 10, height = 10)


#sama cycif

data_sciset <- read.csv("E:/sciset/data_20230811.csv")
data_sciset$Global_celltype_tumor_stroma_merged <- data_sciset$GlobalCellType2
data_sciset[which(data_sciset$GlobalCellType2 == "EMT" | data_sciset$GlobalCellType2 == "Epithelial" | data_sciset$GlobalCellType2 == "Proliferating.EMT" | data_sciset$GlobalCellType2 == "Proliferating.epithelial"), "Global_celltype_tumor_stroma_merged"] <- "tumor"
data_sciset[which(data_sciset$GlobalCellType2 == "Fibroblast" | data_sciset$GlobalCellType2 == "SMA.CD31.positive.cell" | data_sciset$GlobalCellType2 == "Myofibroblast" | data_sciset$GlobalCellType2 == "Endothelial_cell" | data_sciset$GlobalCellType2 == "Desmin.positive.cell" | data_sciset$GlobalCellType2 == "SMA.Desmin.positive.cell" | data_sciset$GlobalCellType2 == "Other"), "Global_celltype_tumor_stroma_merged"] <- "stroma"
data_sciset[which(data_sciset$GlobalCellType2 == "CD11c.myeloid" | data_sciset$GlobalCellType2 == "IBA1.CD11c.Macrophages" | data_sciset$GlobalCellType2 == "IBA1.CD163.Macrophages" | data_sciset$GlobalCellType2 == "CD163.Macrophages"), "Global_celltype_tumor_stroma_merged"] <- "Myeloid cells"

clinical <- unique(data_sciset[, c("Sample_code", "Stage", "Patient")])

data_sciset[which(data_sciset$Global_celltype_tumor_stroma_merged == "Endothelial.cell"), "Global_celltype_tumor_stroma_merged"] <- "stroma"
TSEI_percentages <- data_sciset %>% group_by(Sample_code, Global_celltype_tumor_stroma_merged) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))%>%
  ungroup() %>%
  complete(Sample_code, Global_celltype_tumor_stroma_merged,
           fill = list(n = 0, proportion = 0))

#clinical[which(clinical$Sample_code == "H110_pOva"), "Sample_code"] <- "H110_pOvaR"
TSEI_percentages <- merge(TSEI_percentages, clinical[, c("Sample_code", "Stage", "Patient")], by="Sample_code")

metacluster_percentages <- as.data.frame(TSEI_percentages)

order <- metacluster_percentages[which(metacluster_percentages$Global_celltype_tumor_stroma_merged == "tumor"),]

plot_order <- order[order(order$proportion, decreasing = TRUE),]

#only matching samples with scRNAseq and bulk
plot_order <- plot_order[-which(plot_order$Sample_code == "H166_iOth" | plot_order$Sample_code == "H091_pOva" | plot_order$Sample_code == "H122_pPer"),]

data_cycif_wide <- metacluster_percentages[metacluster_percentages$Sample_code %in% plot_order$Sample_code,]

df <- data_cycif_wide
df[which(df$Global_celltype_tumor_stroma_merged == "Myeloid cells"), "Global_celltype_tumor_stroma_merged"] <- "myeloid"


pdf("D:/Sciset/scimap/plots/barplots/barplot_tcycif_celltypes.pdf", width=7, height=4)
p <- ggplot(df, aes(x = factor(Sample_code, levels=plot_order$Sample_code), y = proportion, fill = Global_celltype_tumor_stroma_merged))+
  geom_bar(stat = "identity", colour="black", position = "stack") + xlab("Sample") + ylab("Proportion")+  scale_fill_manual(values=group.colors)+
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_discrete(expand = c(0,0)) + theme(legend.position="right")


print(p)
dev.off()



#barplots of the number of cells
TSEI_percentages <- data_sciset %>% group_by(Sample_code) %>% summarise(n = n())
library(ggplot2)

TSEI_percentages$Sample_code <- substr(TSEI_percentages$Sample_code, start = 1, stop = 6)

TSEI_percentages <- TSEI_percentages[which(TSEI_percentages$Sample_code %in% plot_order$Sample_code),]

pdf("D:/Sciset/scimap/plots/barplots/barplot_cycif_cells_celltypes.pdf", width=7, height=2)
p <- ggplot(TSEI_percentages, aes(x = factor(Sample_code, levels=plot_order$Sample_code), y = n))+
  geom_bar(stat = "identity", position = "stack", fill="grey") + xlab("Sample") + ylab("n")+ 
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(legend.position="none")
p
dev.off()

#barplot scRNAseq
counts_scRNAseq <- read.csv("D:/Sciset/scimap/Figures/for_figure_1/scrnaseq_cell_counts_celltype_high_v3.txt", sep="\t")
counts_scRNAseq$Sample_code <- rownames(counts_scRNAseq)
counts_scRNAseq$Sample_code <- substr(counts_scRNAseq$Sample_code, start = 1, stop = 6)
plot_order$Sample_code <- substr(plot_order$Sample_code, start = 1, stop = 9)
counts_scRNAseq <- counts_scRNAseq[counts_scRNAseq$Sample_code %in% plot_order$Sample_code, ]
#

counts_scRNAseq$total <- rowSums(counts_scRNAseq[, c(1:8)])

pdf("D:/Sciset/scimap/plots/barplots/barplot_scRNAseq_counts_celltypes.pdf", width=7, height=2)
p <- ggplot(counts_scRNAseq, aes(x = factor(Sample_code, levels=plot_order$Sample_code), y = total))+
  geom_bar(stat = "identity", position = "stack", fill="grey") + xlab("Sample") + ylab("n")+ 
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(legend.position="none")
p
dev.off()



#scRNAseq changes in omentum, also all paired

data_sc <- read.table("D:/Sciset/scimap/Figures/for_figure_1/scrnaseq_cell_fractions_celltype_low.txt", header=TRUE, sep="\t")
data_sc <- read.table("D:/Sciset/scimap/Figures/for_figure_1/scrnaseq_cell_fractions_celltype_high.txt", header=TRUE, sep="\t")


data_sc$Sample_code <- rownames(data_sc)
data_sc$Paired <- FALSE
data_sc[c(1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 21, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51), "Paired"] <- TRUE

data_sc <- data_sc[which(data_sc$Paired == TRUE),]
data_sc <- data_sc[-c(23),]

data_sc$Stage <- "primary"
data_sc[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41), "Stage"] <- "interval"

library(tidyr)
library(dplyr)
data_sc_long <- pivot_longer(data_sc, cols=CD16..NK.cells:Type.17.helper.T.cells, values_to='proportion', names_to='Global_celltype_tumor_stroma_merged')
metacluster_percentages <- data_sc_long

test <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    Stage_pval = wilcox.test(proportion[Stage == 'interval'], 
                             proportion[Stage == 'primary'], paired=TRUE)$p.value)


df2 <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged, Stage) %>%
  summarize(proportion = mean(proportion))

df2

df3 <- df2 %>%
  pivot_wider(names_from = 'Stage', values_from = 'proportion')
df3

df3 <- df3 %>%
  mutate(foldChange = log2(`interval`/`primary`))

test <- read.table("D:/Sciset/scimap/Figures/for_figure_1/propeller_cell_abundance_allpairs_cell_type_low.txt", header=TRUE, sep="\t")
test$Global_celltype_tumor_stroma_merged <- rownames(test)
test$Global_celltype_tumor_stroma_merged <- sub(" ", ".", test$Global_celltype_tumor_stroma_merged)
test[which(test$Global_celltype_tumor_stroma_merged == "Memory.B cells"), "Global_celltype_tumor_stroma_merged"] <- 'Memory.B.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Tem/Trm.cytotoxic T cells"), "Global_celltype_tumor_stroma_merged"] <- 'Tem.Trm.cytotoxic.T.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Tcm/Naive.helper T cells"), "Global_celltype_tumor_stroma_merged"] <- 'Tcm.Naive.helper.T.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Type.17 helper T cells"), "Global_celltype_tumor_stroma_merged"] <- 'Type.17.helper.T.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Regulatory.T cells"), "Global_celltype_tumor_stroma_merged"] <- 'Regulatory.T.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Naive.B cells"), "Global_celltype_tumor_stroma_merged"] <- 'Naive.B.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "CD16-.NK cells"), "Global_celltype_tumor_stroma_merged"] <- 'CD16..NK.cells.1'
test[which(test$Global_celltype_tumor_stroma_merged == "CD16+.NK cells"), "Global_celltype_tumor_stroma_merged"] <- 'CD16..NK.cells'


df2_sc <- merge(test, df3, by="Global_celltype_tumor_stroma_merged")

#make a dot plot

rownames(df2_sc) <- df2_sc$Global_celltype_tumor_stroma_merged

# col_c <- colorRamp2(c(-1.3, 0, 1.3), hcl_palette = "Blue-Red 3")
# hmap <- Heatmap(as.matrix(df2_sc[, -c(1, 2, 3, 4)]), name="log2fc IDS \nvs chemonaive", cluster_rows = T,
#                 row_names_gp = gpar(fontsize = 8),clustering_method_rows ="ward.D2",column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),cluster_columns = F,
#                 left_annotation = NULL, row_dend_width = unit(1.3, "cm"),show_row_names = T, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
#                 show_row_dend=T, show_column_dend=T, column_title = "subtypes", row_title = "cell type",
#                 border="black", width = unit(6.3, "cm"),height = unit(3, "cm"),
#                 rect_gp = gpar(col = "white", lwd = 2), col=col_c)
# 
# 
# draw(hmap)
# 
# pdf("E:/sciset/bulk/comparison_pre_post_high_sc_20231219.pdf", width = 10, height = 2.5)
# draw(hmap)
# dev.off()

colnames(df2_sc)[c(5,9)] <- c("pvalue", "log_foldchange")

to.plot <- df2_sc
to.plot$comparison <- "high_level_all_pairs"

#row_order <- rev(c("tcycif", "scRNAseq", "bulkRNAseq"))
p <- ggplot(to.plot,  aes(x=factor(Global_celltype_tumor_stroma_merged), y=factor(comparison))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90",  breaks=waiver()) +
  scale_size_area("p-values", trans="log10",max_size = 5, breaks=c(1e-4, 1e-2, 0.05), limits=c(1e-6, 1e1)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)


pdf("D:/Sciset/scimap/Figures/for_figure_1/comparison_pre_post_high_sc_20231220_v3.pdf", width = 10, height = 4.5)
print(p)
dev.off()

df2_sc$pvalue


df2_sc[which(df2_sc$pvalue<0.05), "Global_celltype_tumor_stroma_merged"]

#then only omentum pairs

metacluster_percentages$Sample_code

metacluster_percentages <- metacluster_percentages[metacluster_percentages$Sample_code %like% "Ome", ]

unique(metacluster_percentages$Sample_code)

metacluster_percentages <- metacluster_percentages[-which(metacluster_percentages$Sample_code =="H086_pOme" |metacluster_percentages$Sample_code == "H114_pOme2"|metacluster_percentages$Sample_code == "H122_iOme1"|metacluster_percentages$Sample_code == "H142_iOme2"|metacluster_percentages$Sample_code == "H144_iOme1"|metacluster_percentages$Sample_code == "H147_iOme"|metacluster_percentages$Sample_code =="H185_iOme1" |metacluster_percentages$Sample_code == "H190_iOme1"|metacluster_percentages$Sample_code =="H256_iOme1"  ),]

test <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged) %>%
  summarise(
    Stage_pval = wilcox.test(proportion[Stage == 'interval'], 
                             proportion[Stage == 'primary'], paired=TRUE)$p.value)


df2 <- metacluster_percentages %>%
  group_by(Global_celltype_tumor_stroma_merged, Stage) %>%
  summarize(proportion = mean(proportion))

df2

df3 <- df2 %>%
  pivot_wider(names_from = 'Stage', values_from = 'proportion')
df3

df3 <- df3 %>%
  mutate(foldChange = log2(`interval`/`primary`))


test <- read.table("D:/Sciset/scimap/Figures/for_figure_1/propeller_cell_abundance_omepairs_cell_type_low.txt", header=TRUE, sep="\t")
test$Global_celltype_tumor_stroma_merged <- rownames(test)
test$Global_celltype_tumor_stroma_merged <- sub(" ", ".", test$Global_celltype_tumor_stroma_merged)
test[which(test$Global_celltype_tumor_stroma_merged == "Memory.B cells"), "Global_celltype_tumor_stroma_merged"] <- 'Memory.B.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Tem/Trm.cytotoxic T cells"), "Global_celltype_tumor_stroma_merged"] <- 'Tem.Trm.cytotoxic.T.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Tcm/Naive.helper T cells"), "Global_celltype_tumor_stroma_merged"] <- 'Tcm.Naive.helper.T.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Type.17 helper T cells"), "Global_celltype_tumor_stroma_merged"] <- 'Type.17.helper.T.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Regulatory.T cells"), "Global_celltype_tumor_stroma_merged"] <- 'Regulatory.T.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "Naive.B cells"), "Global_celltype_tumor_stroma_merged"] <- 'Naive.B.cells'
test[which(test$Global_celltype_tumor_stroma_merged == "CD16-.NK cells"), "Global_celltype_tumor_stroma_merged"] <- 'CD16..NK.cells.1'
test[which(test$Global_celltype_tumor_stroma_merged == "CD16+.NK cells"), "Global_celltype_tumor_stroma_merged"] <- 'CD16..NK.cells'


df2_sc <- merge(test, df3, by="Global_celltype_tumor_stroma_merged")



to.plot <- df2_sc
to.plot$comparison <- "high_level_omentum_pairs"

colnames(to.plot)[c(5, 9)] <- c("pvalue", "log_foldchange")
p <- ggplot(to.plot,  aes(x=factor(Global_celltype_tumor_stroma_merged), y=factor(comparison))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90",  breaks=waiver()) +
  scale_size_area("p-values", trans="log10",max_size = 5, breaks=c(1e-4, 1e-2, 0.05), limits=c(1e-6, 1e1)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)


pdf("D:/Sciset/scimap/Figures/for_figure_1/comparison_pre_post_high_sc_omentum_20231220_v3.pdf", width = 10, height = 4.5)
print(p)
dev.off()

df2_sc[which(df2_sc$P.Value<0.05), "Global_celltype_tumor_stroma_merged"]


