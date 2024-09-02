#tcell functional markers in RCNs
#Figure 3d
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(data.table)

data <- read.csv("D:/sciset/data_20231125.csv")

clinical <- read.csv("P:/h345/afarkkilab/Projects/Sciset/clinical_data.csv")

data <- merge[data, clinical, by="imageid"]
data <- data[which(data$Paired==TRUE),]

#dotplot of the change in CD8+T-cell states chemo-naive to interval with FRD correction

tcells <- data[, c("pSTAT1", "MHCI", "SNAT1", "Ki67", "CD45RO","TIM3", "Stage", "Sample_code", "GlobalCellType2", "neighbordood_cluster2")]

p_values <- tcells[which(tcells$GlobalCellType2 == "CD8.T.cells"),] %>% group_by(neighbordood_cluster2)%>%
  summarise_each(funs(wilcox.test(.[Stage == "primary"], .[Stage == "interval"])$p.value), vars = pSTAT1:TIM3)
p_values <- as.data.frame(p_values)

rownames(p_values) <- p_values[, 1]
p_values <- p_values[, -c(1)]
colnames(p_values) <- colnames(tcells)[1:6]


median_expression <- tcells[which(tcells$GlobalCellType2 == "CD8.T.cells"),] %>% 
  group_by(Stage, neighbordood_cluster2) %>%
  summarize(
    median_pSTAT1 = median(pSTAT1),
    median_MHCI = median(MHCI),
    median_SNAT1 = median(SNAT1),
    median_Ki67 = median(Ki67),
    median_CD45RO = median(CD45RO),
    median_TIM3 = median(TIM3))


median_expression_long <- pivot_longer(median_expression, cols = 3:8, names_to = "Marker", values_to = "Expression")

median_expression_wide <- pivot_wider(median_expression_long, names_from = "Stage", values_from = "Expression")

#fold change

median_expression_wide$fold_change <- median_expression_wide$`interval`/median_expression_wide$primary
median_expression_wide$log_fold_change <- log2(median_expression_wide$fold_change)

fold_change_data <- median_expression_wide[, c("neighbordood_cluster2", "Marker", "log_fold_change")]


log_fold_change_data_wide <- pivot_wider(fold_change_data, names_from = "neighbordood_cluster2", values_from = "log_fold_change")

log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]


log_fold_change_data_wide <- as.data.table(t(log_fold_change_data_wide), keep.colnames = T, keep.rownames = T)
log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]

dotplot_data <- log_fold_change_data_wide

colnames(dotplot_data) <- c("foldchange_pSTAT1", "foldchange_MHCI", "foldchange_SNAT1", "foldchange_Ki67", "foldchange_CD45RO", "foldchange_TIM3")

dotplot_data <- as.data.frame(dotplot_data)

setDT(dotplot_data, keep.rownames = "Neighbordood_cluster2")[]
dotplot_data_long <- pivot_longer(dotplot_data, cols = 2:7, names_to = "Foldchange")
colnames(dotplot_data_long) <- c("Neighbordood_cluster2", "Marker", "log_foldchange")

setDT(p_values, keep.rownames = "Neighbordood_cluster2")[]
p_values_long <- pivot_longer(p_values, cols = 2:7, names_to = "pvalue")
colnames(p_values_long) <- c("Neighbordood_cluster2", "Marker", "pvalue")


poly.results <- cbind(dotplot_data_long, p_values_long, by="Neighbordood_cluster2")
poly.results <- poly.results[, c(1, 2, 3, 6)]

#row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface"))
#row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT"))
row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"))

#poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == ""),]
#poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == "Epithelial" | poly.results$Neighbordood_cluster2 == "EMT" | poly.results$Neighbordood_cluster2 == "Epithelial_EMT"| poly.results$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"| poly.results$Neighbordood_cluster2 == "Proliferating.epithelial" | poly.results$Neighbordood_cluster2 =="Proliferating.EMT"),]
poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == "Epithelial" | poly.results$Neighbordood_cluster2 == "EMT" | poly.results$Neighbordood_cluster2 == "Epithelial_EMT"| poly.results$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"| poly.results$Neighbordood_cluster2 == "Proliferating.epithelial" | poly.results$Neighbordood_cluster2 =="Proliferating.EMT" | poly.results$Neighbordood_cluster2 == "Fibroblast"| poly.results$Neighbordood_cluster2 ==  "stroma" | poly.results$Neighbordood_cluster2 == "SMA.Desmin.myofibroblast" | poly.results$Neighbordood_cluster2 == "SMA.CD31.positive" | poly.results$Neighbordood_cluster2 =="Desmin.positive" | poly.results$Neighbordood_cluster2 == "Myofibroblast"),]

poly.results$pvalue <- p.adjust(poly.results$pvalue)

truncate.df <- function(df, na.cutoff=0.05, na.var="pvalue", na.var.boundary=1e-50, range.lims=c(-0.1, 0.1), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(poly.results) 


to.plot[which(to.plot$Neighbordood_cluster2 == "EMT"), "Neighbordood_cluster2"] <- "RCN5"
to.plot[which(to.plot$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"), "Neighbordood_cluster2"] <- "RCN2"
to.plot[which(to.plot$Neighbordood_cluster2 == "Fibroblast"), "Neighbordood_cluster2"] <- "RCN13"
to.plot[which(to.plot$Neighbordood_cluster2 == "stroma"), "Neighbordood_cluster2"] <- "RCN10"
to.plot[which(to.plot$Neighbordood_cluster2 == "CD11c.myeloid"), "Neighbordood_cluster2"] <- "RCN16"
to.plot[which(to.plot$Neighbordood_cluster2 == "SMA.Desmin.myofibroblast"), "Neighbordood_cluster2"] <- "RCN9"
to.plot[which(to.plot$Neighbordood_cluster2 == "Epithelial_EMT"), "Neighbordood_cluster2"] <- "RCN4"
to.plot[which(to.plot$Neighbordood_cluster2 == "SMA.CD31.positive"), "Neighbordood_cluster2"] <- "RCN11"
to.plot[which(to.plot$Neighbordood_cluster2 == "tumor-stroma-interface"), "Neighbordood_cluster2"] <- "RCN7"
to.plot[which(to.plot$Neighbordood_cluster2 == "Immune"), "Neighbordood_cluster2"] <- "RCN18"
to.plot[which(to.plot$Neighbordood_cluster2 == "CD8_CD4_T.cells"), "Neighbordood_cluster2"] <- "RCN17"
to.plot[which(to.plot$Neighbordood_cluster2 == "Proliferating.EMT"), "Neighbordood_cluster2"] <- "RCN6"
to.plot[which(to.plot$Neighbordood_cluster2 == "Desmin.positive"), "Neighbordood_cluster2"] <- "RCN8"
to.plot[which(to.plot$Neighbordood_cluster2 == "Proliferating.epithelial"), "Neighbordood_cluster2"] <- "RCN1"
to.plot[which(to.plot$Neighbordood_cluster2 == "Epithelial"), "Neighbordood_cluster2"] <- "RCN3"
to.plot[which(to.plot$Neighbordood_cluster2 == "IBA1.CD163.Macrophages"), "Neighbordood_cluster2"] <- "RCN15"
to.plot[which(to.plot$Neighbordood_cluster2 == "Myofibroblast"), "Neighbordood_cluster2"] <- "RCN12"
to.plot[which(to.plot$Neighbordood_cluster2 == "Macrophages"), "Neighbordood_cluster2"] <- "RCN14"


row_order = rev(c("RCN16","RCN14","RCN15","RCN18","RCN17","RCN7","RCN1","RCN2","RCN3","RCN4","RCN5","RCN6","RCN8","RCN9","RCN10","RCN11","RCN12","RCN13"))

p <- ggplot(to.plot,  aes(x=factor(Marker), y=factor(Neighbordood_cluster2, level=row_order))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "gray39", high="#cc2127", mid = "gray90", breaks=seq(-0.1,0.1,0.05), limits=c(-0.1, 0.1)) +
  scale_size_area("p-values", trans="log10",max_size = 1.5, breaks=c(1e-50, 1e-20, 1e-10, 1e-5, 1e-1, 0.05), limits=c(1e-50,0.05)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)

pdf("D:/Sciset/scimap/plots/dotplot/dotplot_CD8.T.cells_in_all_neighborhoods_gray_20231126.pdf", width=4)
print(p)
dev.off()



#heatmaps of functional marker expression in CD8+T-cells in different neighborhoods in chemo-naive and IDS samples



# use the previous data: median expression wide

median_expression_wide_primary <- pivot_wider(median_expression_wide[c(1, 2, 4)], names_from = Marker, values_from = primary)

median_expression_wide_interval <- pivot_wider(median_expression_wide[c(1, 2, 3)], names_from = Marker, values_from = interval)

median_expression_wide_primary <- as.data.frame(median_expression_wide_primary)
median_expression_wide_interval <- as.data.frame(median_expression_wide_interval)


rownames(median_expression_wide_primary) <- median_expression_wide_primary$neighbordood_cluster2

rownames(median_expression_wide_interval) <- median_expression_wide_interval$neighbordood_cluster2

median_expression_wide_primary <- median_expression_wide_primary[, c(1, 6, 5, 3, 2, 4, 7)]
median_expression_wide_interval <- median_expression_wide_interval[, c(1, 6, 5, 3, 2, 4, 7)]


median_expression_wide_primary <- median_expression_wide_primary[c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"),]


median_expression_wide_interval <- median_expression_wide_interval[c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"),]

library(circlize)
col_fun<-colorRamp2(c(-4,-2,0,2,4), c("darkblue","skyblue","white", "brown1", "darkred"))

pdf("D:/Sciset/scimap/plots/heatmaps/heatmap_cd8_states_neighborhoods_primary_20231126.pdf", width=6)
Heatmap(scale(median_expression_wide_primary[,-1]), cluster_rows = F, col=col_fun, show_row_names = T, cluster_columns = F, width = unit(3, "cm"),height = unit(9, "cm"), heatmap_legend_param = list(col_fun = col_fun, at = c(-4.5, -2, 0, 2, 4.5)))
dev.off()


pdf("D:/Sciset/scimap/plots/heatmaps/heatmap_cd8_states_neighborhoods_interval_20231126.pdf", width=6)
Heatmap(scale(median_expression_wide_interval[,-1]), cluster_rows = F, show_row_names = T,col=col_fun, cluster_columns = F, width = unit(3, "cm"),height = unit(9, "cm"), heatmap_legend_param = list(col_fun = col_fun, at = c(-4.5, -2, 0, 2, 4.5)))
dev.off()



#same for cd4

p_values <- tcells[which(tcells$GlobalCellType2 == "CD4.T.cells"),] %>% group_by(neighbordood_cluster2)%>%
  summarise_each(funs(wilcox.test(.[Stage == "primary"], .[Stage == "interval"])$p.value), vars = pSTAT1:TIM3)
p_values <- as.data.frame(p_values)

rownames(p_values) <- p_values[, 1]
p_values <- p_values[, -c(1)]
colnames(p_values) <- colnames(tcells)[1:6]


median_expression <- tcells[which(tcells$GlobalCellType2 == "CD4.T.cells"),] %>% 
  group_by(Stage, neighbordood_cluster2) %>%
  summarize(
    median_pSTAT1 = median(pSTAT1),
    median_MHCI = median(MHCI),
    median_SNAT1 = median(SNAT1),
    median_Ki67 = median(Ki67),
    median_CD45RO = median(CD45RO),
    median_TIM3 = median(TIM3))

library(tidyr)
median_expression_long <- pivot_longer(median_expression, cols = 3:8, names_to = "Marker", values_to = "Expression")

median_expression_wide <- pivot_wider(median_expression_long, names_from = "Stage", values_from = "Expression")

#fold change

median_expression_wide$fold_change <- median_expression_wide$`interval`/median_expression_wide$primary
median_expression_wide$log_fold_change <- log2(median_expression_wide$fold_change)

fold_change_data <- median_expression_wide[, c("neighbordood_cluster2", "Marker", "log_fold_change")]


log_fold_change_data_wide <- pivot_wider(fold_change_data, names_from = "neighbordood_cluster2", values_from = "log_fold_change")

log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]


log_fold_change_data_wide <- as.data.table(t(log_fold_change_data_wide), keep.colnames = T, keep.rownames = T)
log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]

dotplot_data <- log_fold_change_data_wide

colnames(dotplot_data) <- c("foldchange_pSTAT1", "foldchange_MHCI", "foldchange_SNAT1", "foldchange_Ki67", "foldchange_CD45RO", "foldchange_TIM3")

dotplot_data <- as.data.frame(dotplot_data)

setDT(dotplot_data, keep.rownames = "Neighbordood_cluster2")[]
dotplot_data_long <- pivot_longer(dotplot_data, cols = 2:7, names_to = "Foldchange")
colnames(dotplot_data_long) <- c("Neighbordood_cluster2", "Marker", "log_foldchange")

setDT(p_values, keep.rownames = "Neighbordood_cluster2")[]
p_values_long <- pivot_longer(p_values, cols = 2:7, names_to = "pvalue")
colnames(p_values_long) <- c("Neighbordood_cluster2", "Marker", "pvalue")


poly.results <- cbind(dotplot_data_long, p_values_long, by="Neighbordood_cluster2")
poly.results <- poly.results[, c(1, 2, 3, 6)]

#row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface"))
#row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT"))
row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"))

#poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == ""),]
#poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == "Epithelial" | poly.results$Neighbordood_cluster2 == "EMT" | poly.results$Neighbordood_cluster2 == "Epithelial_EMT"| poly.results$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"| poly.results$Neighbordood_cluster2 == "Proliferating.epithelial" | poly.results$Neighbordood_cluster2 =="Proliferating.EMT"),]
poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == "Epithelial" | poly.results$Neighbordood_cluster2 == "EMT" | poly.results$Neighbordood_cluster2 == "Epithelial_EMT"| poly.results$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"| poly.results$Neighbordood_cluster2 == "Proliferating.epithelial" | poly.results$Neighbordood_cluster2 =="Proliferating.EMT" | poly.results$Neighbordood_cluster2 == "Fibroblast"| poly.results$Neighbordood_cluster2 ==  "stroma" | poly.results$Neighbordood_cluster2 == "SMA.Desmin.myofibroblast" | poly.results$Neighbordood_cluster2 == "SMA.CD31.positive" | poly.results$Neighbordood_cluster2 =="Desmin.positive" | poly.results$Neighbordood_cluster2 == "Myofibroblast"),]

poly.results$pvalue <- p.adjust(poly.results$pvalue)

truncate.df <- function(df, na.cutoff=0.05, na.var="pvalue", na.var.boundary=1e-50, range.lims=c(-0.1, 0.1), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(poly.results) 


to.plot[which(to.plot$Neighbordood_cluster2 == "EMT"), "Neighbordood_cluster2"] <- "RCN5"
to.plot[which(to.plot$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"), "Neighbordood_cluster2"] <- "RCN2"
to.plot[which(to.plot$Neighbordood_cluster2 == "Fibroblast"), "Neighbordood_cluster2"] <- "RCN13"
to.plot[which(to.plot$Neighbordood_cluster2 == "stroma"), "Neighbordood_cluster2"] <- "RCN10"
to.plot[which(to.plot$Neighbordood_cluster2 == "CD11c.myeloid"), "Neighbordood_cluster2"] <- "RCN16"
to.plot[which(to.plot$Neighbordood_cluster2 == "SMA.Desmin.myofibroblast"), "Neighbordood_cluster2"] <- "RCN9"
to.plot[which(to.plot$Neighbordood_cluster2 == "Epithelial_EMT"), "Neighbordood_cluster2"] <- "RCN4"
to.plot[which(to.plot$Neighbordood_cluster2 == "SMA.CD31.positive"), "Neighbordood_cluster2"] <- "RCN11"
to.plot[which(to.plot$Neighbordood_cluster2 == "tumor-stroma-interface"), "Neighbordood_cluster2"] <- "RCN7"
to.plot[which(to.plot$Neighbordood_cluster2 == "Immune"), "Neighbordood_cluster2"] <- "RCN18"
to.plot[which(to.plot$Neighbordood_cluster2 == "CD8_CD4_T.cells"), "Neighbordood_cluster2"] <- "RCN17"
to.plot[which(to.plot$Neighbordood_cluster2 == "Proliferating.EMT"), "Neighbordood_cluster2"] <- "RCN6"
to.plot[which(to.plot$Neighbordood_cluster2 == "Desmin.positive"), "Neighbordood_cluster2"] <- "RCN8"
to.plot[which(to.plot$Neighbordood_cluster2 == "Proliferating.epithelial"), "Neighbordood_cluster2"] <- "RCN1"
to.plot[which(to.plot$Neighbordood_cluster2 == "Epithelial"), "Neighbordood_cluster2"] <- "RCN3"
to.plot[which(to.plot$Neighbordood_cluster2 == "IBA1.CD163.Macrophages"), "Neighbordood_cluster2"] <- "RCN15"
to.plot[which(to.plot$Neighbordood_cluster2 == "Myofibroblast"), "Neighbordood_cluster2"] <- "RCN12"
to.plot[which(to.plot$Neighbordood_cluster2 == "Macrophages"), "Neighbordood_cluster2"] <- "RCN14"


row_order = rev(c("RCN16","RCN14","RCN15","RCN18","RCN17","RCN7","RCN1","RCN2","RCN3","RCN4","RCN5","RCN6","RCN8","RCN9","RCN10","RCN11","RCN12","RCN13"))


p <- ggplot(to.plot,  aes(x=factor(Marker), y=factor(Neighbordood_cluster2, level=row_order))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "gray39", high="#cc2127", mid = "gray90", breaks=seq(-0.1,0.1,0.05), limits=c(-0.1, 0.1)) +
  scale_size_area("p-values", trans="log10",max_size = 1.5, breaks=c(1e-50, 1e-20, 1e-10, 1e-5, 1e-1, 0.05), limits=c(1e-50,0.05)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)

pdf("D:/Sciset/scimap/plots/dotplot/dotplot_CD4.T.cells_in_all_neighborhoods_gray_20231126.pdf", width=4)
print(p)
dev.off()


#heatmaps of functional marker expression in different neighborhoods in chemo-naive and IDS samples

#median expression wide

median_expression_wide_primary <- pivot_wider(median_expression_wide[c(1, 2, 4)], names_from = Marker, values_from = primary)

median_expression_wide_interval <- pivot_wider(median_expression_wide[c(1, 2, 3)], names_from = Marker, values_from = interval)

median_expression_wide_primary <- as.data.frame(median_expression_wide_primary)
median_expression_wide_interval <- as.data.frame(median_expression_wide_interval)

rownames(median_expression_wide_primary) <- median_expression_wide_primary$neighbordood_cluster2

rownames(median_expression_wide_interval) <- median_expression_wide_interval$neighbordood_cluster2

median_expression_wide_primary <- median_expression_wide_primary[, c(1, 6, 5, 3, 2, 4, 7)]
median_expression_wide_interval <- median_expression_wide_interval[, c(1, 6, 5, 3, 2, 4, 7)]


median_expression_wide_primary <- median_expression_wide_primary[c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"),]


median_expression_wide_interval <- median_expression_wide_interval[c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"),]

col_fun<-colorRamp2(c(-4,-2,0,2,4), c("darkblue","skyblue","white", "brown1", "darkred"))

pdf("D:/Sciset/scimap/plots/heatmaps/heatmap_cd4_states_neighborhoods_primary_20231126.pdf", width=6)
Heatmap(scale(median_expression_wide_primary[,-1]), cluster_rows = F, col=col_fun, show_row_names = T, width = unit(3, "cm"),height = unit(9, "cm"), cluster_columns = F,  heatmap_legend_param = list(col_fun = col_fun, at = c(-4.5, -2, 0, 2, 4.5)))
dev.off()


pdf("D:/Sciset/scimap/plots/heatmaps/heatmap_cd4_states_neighborhoods_interval_20231126.pdf", width=6)
Heatmap(scale(median_expression_wide_interval[,-1]), cluster_rows = F, show_row_names = T,col=col_fun, width = unit(3, "cm"),height = unit(9, "cm"), cluster_columns = F, heatmap_legend_param = list(col_fun = col_fun, at = c(-4.5, -2, 0, 2, 4.5)))
dev.off()


#sama for T-regs


p_values <- tcells[which(tcells$GlobalCellType2 == "FOXP3.CD4.Tregs"),] %>% group_by(neighbordood_cluster2)%>%
  summarise_each(funs(wilcox.test(.[Stage == "primary"], .[Stage == "interval"])$p.value), vars = pSTAT1:TIM3)
p_values <- as.data.frame(p_values)

rownames(p_values) <- p_values[, 1]
p_values <- p_values[, -c(1)]
colnames(p_values) <- colnames(tcells)[1:6]


median_expression <- tcells[which(tcells$GlobalCellType2 == "FOXP3.CD4.Tregs"),] %>% 
  group_by(Stage, neighbordood_cluster2) %>%
  summarize(
    median_pSTAT1 = median(pSTAT1),
    median_MHCI = median(MHCI),
    median_SNAT1 = median(SNAT1),
    median_Ki67 = median(Ki67),
    median_CD45RO = median(CD45RO),
    median_TIM3 = median(TIM3))

library(tidyr)
median_expression_long <- pivot_longer(median_expression, cols = 3:8, names_to = "Marker", values_to = "Expression")

median_expression_wide <- pivot_wider(median_expression_long, names_from = "Stage", values_from = "Expression")

#fold change

median_expression_wide$fold_change <- median_expression_wide$`interval`/median_expression_wide$primary
median_expression_wide$log_fold_change <- log2(median_expression_wide$fold_change)

fold_change_data <- median_expression_wide[, c("neighbordood_cluster2", "Marker", "log_fold_change")]


log_fold_change_data_wide <- pivot_wider(fold_change_data, names_from = "neighbordood_cluster2", values_from = "log_fold_change")

log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]


log_fold_change_data_wide <- as.data.table(t(log_fold_change_data_wide), keep.colnames = T, keep.rownames = T)
log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]

dotplot_data <- log_fold_change_data_wide

colnames(dotplot_data) <- c("foldchange_pSTAT1", "foldchange_MHCI", "foldchange_SNAT1", "foldchange_Ki67", "foldchange_CD45RO", "foldchange_TIM3")

dotplot_data <- as.data.frame(dotplot_data)

setDT(dotplot_data, keep.rownames = "Neighbordood_cluster2")[]
dotplot_data_long <- pivot_longer(dotplot_data, cols = 2:7, names_to = "Foldchange")
colnames(dotplot_data_long) <- c("Neighbordood_cluster2", "Marker", "log_foldchange")

setDT(p_values, keep.rownames = "Neighbordood_cluster2")[]
p_values_long <- pivot_longer(p_values, cols = 2:7, names_to = "pvalue")
colnames(p_values_long) <- c("Neighbordood_cluster2", "Marker", "pvalue")


poly.results <- cbind(dotplot_data_long, p_values_long, by="Neighbordood_cluster2")
poly.results <- poly.results[, c(1, 2, 3, 6)]

#row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface"))
#row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT"))
row_order <- rev(c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"))

#poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == ""),]
#poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == "Epithelial" | poly.results$Neighbordood_cluster2 == "EMT" | poly.results$Neighbordood_cluster2 == "Epithelial_EMT"| poly.results$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"| poly.results$Neighbordood_cluster2 == "Proliferating.epithelial" | poly.results$Neighbordood_cluster2 =="Proliferating.EMT"),]
poly.results <- poly.results[which(poly.results$Neighbordood_cluster2 == "CD11c.myeloid" | poly.results$Neighbordood_cluster2 == "Macrophages" | poly.results$Neighbordood_cluster2 == "IBA1.CD163.Macrophages" | poly.results$Neighbordood_cluster2 == "Immune" | poly.results$Neighbordood_cluster2 == "CD8_CD4_T.cells" | poly.results$Neighbordood_cluster2 == "tumor-stroma-interface" | poly.results$Neighbordood_cluster2 == "Epithelial" | poly.results$Neighbordood_cluster2 == "EMT" | poly.results$Neighbordood_cluster2 == "Epithelial_EMT"| poly.results$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"| poly.results$Neighbordood_cluster2 == "Proliferating.epithelial" | poly.results$Neighbordood_cluster2 =="Proliferating.EMT" | poly.results$Neighbordood_cluster2 == "Fibroblast"| poly.results$Neighbordood_cluster2 ==  "stroma" | poly.results$Neighbordood_cluster2 == "SMA.Desmin.myofibroblast" | poly.results$Neighbordood_cluster2 == "SMA.CD31.positive" | poly.results$Neighbordood_cluster2 =="Desmin.positive" | poly.results$Neighbordood_cluster2 == "Myofibroblast"),]

poly.results$pvalue <- p.adjust(poly.results$pvalue)

truncate.df <- function(df, na.cutoff=0.05, na.var="pvalue", na.var.boundary=1e-50, range.lims=c(-0.1, 0.1), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(poly.results) 


to.plot[which(to.plot$Neighbordood_cluster2 == "EMT"), "Neighbordood_cluster2"] <- "RCN5"
to.plot[which(to.plot$Neighbordood_cluster2 == "epithelial_and_proliferating_epithelial"), "Neighbordood_cluster2"] <- "RCN2"
to.plot[which(to.plot$Neighbordood_cluster2 == "Fibroblast"), "Neighbordood_cluster2"] <- "RCN13"
to.plot[which(to.plot$Neighbordood_cluster2 == "stroma"), "Neighbordood_cluster2"] <- "RCN10"
to.plot[which(to.plot$Neighbordood_cluster2 == "CD11c.myeloid"), "Neighbordood_cluster2"] <- "RCN16"
to.plot[which(to.plot$Neighbordood_cluster2 == "SMA.Desmin.myofibroblast"), "Neighbordood_cluster2"] <- "RCN9"
to.plot[which(to.plot$Neighbordood_cluster2 == "Epithelial_EMT"), "Neighbordood_cluster2"] <- "RCN4"
to.plot[which(to.plot$Neighbordood_cluster2 == "SMA.CD31.positive"), "Neighbordood_cluster2"] <- "RCN11"
to.plot[which(to.plot$Neighbordood_cluster2 == "tumor-stroma-interface"), "Neighbordood_cluster2"] <- "RCN7"
to.plot[which(to.plot$Neighbordood_cluster2 == "Immune"), "Neighbordood_cluster2"] <- "RCN18"
to.plot[which(to.plot$Neighbordood_cluster2 == "CD8_CD4_T.cells"), "Neighbordood_cluster2"] <- "RCN17"
to.plot[which(to.plot$Neighbordood_cluster2 == "Proliferating.EMT"), "Neighbordood_cluster2"] <- "RCN6"
to.plot[which(to.plot$Neighbordood_cluster2 == "Desmin.positive"), "Neighbordood_cluster2"] <- "RCN8"
to.plot[which(to.plot$Neighbordood_cluster2 == "Proliferating.epithelial"), "Neighbordood_cluster2"] <- "RCN1"
to.plot[which(to.plot$Neighbordood_cluster2 == "Epithelial"), "Neighbordood_cluster2"] <- "RCN3"
to.plot[which(to.plot$Neighbordood_cluster2 == "IBA1.CD163.Macrophages"), "Neighbordood_cluster2"] <- "RCN15"
to.plot[which(to.plot$Neighbordood_cluster2 == "Myofibroblast"), "Neighbordood_cluster2"] <- "RCN12"
to.plot[which(to.plot$Neighbordood_cluster2 == "Macrophages"), "Neighbordood_cluster2"] <- "RCN14"


row_order = rev(c("RCN16","RCN14","RCN15","RCN18","RCN17","RCN7","RCN1","RCN2","RCN3","RCN4","RCN5","RCN6","RCN8","RCN9","RCN10","RCN11","RCN12","RCN13"))


p <- ggplot(to.plot,  aes(x=factor(Marker), y=factor(Neighbordood_cluster2, level=row_order))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "gray39", high="#cc2127", mid = "gray90", breaks=seq(-0.1,0.1,0.05), limits=c(-0.1, 0.1)) +
  scale_size_area("p-values", trans="log10",max_size = 1.5, breaks=c(1e-50, 1e-20, 1e-10, 1e-5, 1e-1, 0.05), limits=c(1e-50,0.05)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)

pdf("D:/Sciset/scimap/plots/dotplot/dotplot_Tregs_in_all_neighborhoods_gray_20231126.pdf", width=4)
print(p)
dev.off()


#heatmaps of functional marker expression in different neighborhoods in IDS and chemo-naive samples

#use median expression wide

median_expression_wide_primary <- pivot_wider(median_expression_wide[c(1, 2, 4)], names_from = Marker, values_from = primary)

median_expression_wide_interval <- pivot_wider(median_expression_wide[c(1, 2, 3)], names_from = Marker, values_from = interval)

median_expression_wide_primary <- as.data.frame(median_expression_wide_primary)
median_expression_wide_interval <- as.data.frame(median_expression_wide_interval)

rownames(median_expression_wide_primary) <- median_expression_wide_primary$neighbordood_cluster2

rownames(median_expression_wide_interval) <- median_expression_wide_interval$neighbordood_cluster2

median_expression_wide_primary <- median_expression_wide_primary[, c(1, 6, 5, 3, 2, 4, 7)]
median_expression_wide_interval <- median_expression_wide_interval[, c(1, 6, 5, 3, 2, 4, 7)]


median_expression_wide_primary <- median_expression_wide_primary[c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"),]


median_expression_wide_interval <- median_expression_wide_interval[c("CD11c.myeloid", "Macrophages", "IBA1.CD163.Macrophages", "Immune", "CD8_CD4_T.cells", "tumor-stroma-interface", "Epithelial", "EMT", "Epithelial_EMT", "Proliferating.epithelial", "epithelial_and_proliferating_epithelial", "Proliferating.EMT", "Fibroblast", "stroma", "SMA.Desmin.myofibroblast", "SMA.CD31.positive", "Desmin.positive", "Myofibroblast"),]

col_fun<-colorRamp2(c(-4,-2,0,2,4), c("darkblue","skyblue","white", "brown1", "darkred"))

pdf("D:/Sciset/scimap/plots/heatmaps/heatmap_treg_states_neighborhoods_primary_20231126.pdf", width=6)
Heatmap(scale(median_expression_wide_primary[,-1]), cluster_rows = F, col=col_fun, show_row_names = T, width = unit(3, "cm"),height = unit(9, "cm"), cluster_columns = F,  heatmap_legend_param = list(col_fun = col_fun, at = c(-4.5, -2, 0, 2, 4.5)))
dev.off()


pdf("D:/Sciset/scimap/plots/heatmaps/heatmap_treg_states_neighborhoods_interval_20231126.pdf", width=6)
Heatmap(scale(median_expression_wide_interval[,-1]), cluster_rows = F, show_row_names = T,col=col_fun, width = unit(3, "cm"),height = unit(9, "cm"), cluster_columns = F, heatmap_legend_param = list(col_fun = col_fun, at = c(-4.5, -2, 0, 2, 4.5)))
dev.off()


