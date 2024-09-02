#RCN interactions

# spatial interaction dotplot for RCN7
# Figure 2g

library(dplyr)
library(rstatix)
library(data.table)
library(coin)
library(ggplot2)
library(tidyr)

#unique(data_expr$neighbordood_cluster2)[-1]

#download interactions computed in scimap

for(i in unique(data_expr$neighbordood_cluster2)[12:18]){

RCN_interactions <- read.csv(paste0("D:/Sciset/",i,"_interactions.csv"))

RCN_interactions <- RCN_interactions[, -1]

#here are the interacting pairs

#then a wilcox test per pair across samples
#remove extra columns

RCN_interactions <- RCN_interactions[, -c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48)]

#
#need data to two columns - other has the groups, other the values (many columns of values for each pair)

rownames(RCN_interactions) <- paste0(RCN_interactions$phenotype, "_",RCN_interactions$neighbour_phenotype)

RCN_interactions <- t(RCN_interactions)
RCN_interactions <- as.data.frame(RCN_interactions)

#remove first two rows with redundant information

RCN_interactions <- RCN_interactions[-c(1, 2),]
RCN_interactions$imageid <- rownames(RCN_interactions)

#obs has clinical info
obs <- info

RCN_interactions$imageid <- rownames(RCN_interactions)
RCN_interactions <- merge(RCN_interactions, obs[,c("Sample_code", 'imageid', 'Stage', 'Patient')], by="imageid")

#interaction data to numeric
n <- ncol(RCN_interactions)
n <- n-3
RCN_interactions[, c(2:n)] <- lapply(RCN_interactions[, c(2:n)], function(x) as.numeric(as.character(x)))
#wilcoxon tests - too few samples with B-cells so I omitted those

modelList<-list()
for(m in c(2:290)){
  fmla <- formula(paste(names(RCN_interactions)[m], " ~ Stage"))
  try(modelList[[m]]<-wilcox.test(fmla, data = RCN_interactions, paired = FALSE, conf.int=T), silent=TRUE)
}


#then effect size

modelList2<-list()
for(m in c(2:290)){
  fmla <- formula(paste(names(RCN_interactions)[m], " ~ Stage"))
  modelList2[[m]]<-wilcox_effsize(fmla, data = RCN_interactions, paired = FALSE)
}


#then bind to data frames

#modelList - p.value and data.name
#modelList2 - effsize and .y.


df <- data.frame(matrix(unlist(modelList), nrow=289, byrow=TRUE),stringsAsFactors=FALSE)
colnames(df) <- c("statistic", "p.val", "null", "type", "test","celltypes", "conf.int1", "conf.int2", "estimate")


df2 <- data.frame(matrix(unlist(modelList2), nrow=289, byrow=TRUE),stringsAsFactors=FALSE)
colnames(df2) <- c("celltype.pair", "stage.i", "stage.p", "effect.size", "n1","n2", "magnitude")


#then cbind

df <- cbind(df, df2)

#separate celltypes

#estimate value gives the direction of change

colnames(df)

#df_heatmap <- df

df <- df[, c("p.val", "estimate", "effect.size", "celltype.pair")]

df <- df %>% separate("celltype.pair", c("celltype", "neighbor"), sep="_")

df[, c(1:3)] <- lapply(df[, c(1:3)], function(x) as.numeric(as.character(x)))

df$effect.size <- df$effect.size * sign(df$estimate)

dim(df[which(df$p.val < 0.05),])


#plot the dotplot

# ggplot(df, aes(x=celltype, y=neighbor, color=effect.size)) + geom_point(aes(size=p.val)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="effect size with direction") +
#   scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90", breaks=seq(-1,0,1), limits=c(-1, 1)) +
#   scale_size_area("p.val", trans="log10",max_size = 2.5, breaks=c( 1e-5, 1e-1, 0.05), limits=c(1e-5,0.05)) +
#   theme(axis.text=element_text(size=rel(1)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
#   theme(plot.title = element_text(size = 12, face = "bold"), aspect.ratio=1, panel.border = element_rect(colour = "black", size=1.5, fill=NA))
# 
# 
# 
# #pval correction
# df$p.val <- p.adjust(df$p.val, method = "BH")
# 
# ggplot(df, aes(x=celltype, y=neighbor, color=effect.size)) + geom_point(aes(size=p.val)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="effect size with direction") +
#   scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90", breaks=seq(-1,0,1), limits=c(-1, 1)) +
#   scale_size_area("p.val", trans="log10",max_size = 2.5, breaks=c( 1e-5, 1e-1, 0.05), limits=c(1e-5,0.05)) +
#   theme(axis.text=element_text(size=rel(1)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
#   theme(plot.title = element_text(size = 12, face = "bold"), aspect.ratio=1, panel.border = element_rect(colour = "black", size=1.5, fill=NA))
# 

#take only immuuni neighbors
df_immune <- df[which(df$neighbor== "CD11c.myeloid" | df$neighbor== "CD8.T.cells" |df$neighbor== "CD4.T.cells" |df$neighbor== "FOXP3.CD4.Tregs" |df$neighbor== "CD163.Macrophages" |df$neighbor== "IBA1.CD163.Macrophages" | df$neighbor== "IBA1.CD11c.Macrophages" ),]
df_immune$p.val <- p.adjust(df_immune$p.val, method = "BH")


x_order <- c("Epithelial", "Proliferating.epithelial", "EMT", "Proliferating.EMT", "CD163.Macrophages", "IBA1.CD163.Macrophages", "IBA1.CD11c.Macrophages", "CD11c.myeloid", "CD8.T.cells", "CD4.T.cells", "FOXP3.CD4.Tregs", "Fibroblast", "Myofibroblast", "Endothelial.cell", "SMA.CD31.positive.cell", "SMA.Desmin.positive.cell", "Desmin.positive.cell")
y_order <- c("CD11c.myeloid", "FOXP3.CD4.Tregs", "CD4.T.cells", "CD8.T.cells", "IBA1.CD11c.Macrophages", "IBA1.CD163.Macrophages", "CD163.Macrophages")


df_immune$celltype <- factor(df_immune$celltype, levels=x_order)
df_immune$neighbor <- factor(df_immune$neighbor, levels=rev(y_order))

p <- ggplot(df_immune, aes(x=celltype, y=neighbor, color=effect.size)) + geom_point(aes(size=p.val)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="effect size with direction") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90", breaks=seq(-1,0,1), limits=c(-1, 1)) +
  scale_size("p.val", trans="log10", breaks=c(0.05, 1), range=c(5, 2)) +
  theme(axis.text=element_text(size=rel(1)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

p1 <- ggplot(df_immune[which(df_immune$p.val <= 0.05),], aes(x=celltype, y=neighbor, color=effect.size)) + geom_point(aes(size=p.val)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="effect size with direction") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90", breaks=seq(-1,0,1), limits=c(-1, 1)) +
  scale_size("p.val", trans="log10", breaks=c(0.01, 0.049), range=c(5, 2)) +
  theme(axis.text=element_text(size=rel(1)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"),  panel.border = element_rect(colour = "black", size=1.5, fill=NA))

pdf(paste0("D:/Sciset/scimap/plots/revision/only_immune_neighbors_stage_",i,".pdf"), width = 8, height = 4.5)
print(p)
print(p1)
dev.off()


#show all interactions but highlight only the significant ones with boxes


#now heatmaps separately for primary and interval to show the actual values


mean_int <- RCN_interactions[which(RCN_interactions$Stage == "interval"),]
mean_prim <- RCN_interactions[which(RCN_interactions$Stage == "primary"),]

mean_int <- colMeans(mean_int[, c(2:290)], na.rm=T) 
mean_int <- as.data.frame(mean_int)


mean_prim <- colMeans(mean_prim[, c(2:290)], na.rm=T) 
mean_prim <- as.data.frame(mean_prim)

mean_int$celltype.pair <- rownames(mean_int)
mean_prim$celltype.pair <- rownames(mean_prim)

mean_int <- mean_int %>% separate("celltype.pair", c("celltype", "neighbor"), sep="_")
mean_prim <- mean_prim %>% separate("celltype.pair", c("celltype", "neighbor"), sep="_")


mean_int <- pivot_wider(mean_int, names_from = 'celltype', values_from='mean_int')

mean_prim <- pivot_wider(mean_prim, names_from = 'celltype', values_from='mean_prim')

mean_int <- as.data.frame(mean_int)
rownames(mean_int) <- mean_int$neighbor
mean_int <- as.data.frame(mean_int)

mean_prim <- as.data.frame(mean_prim)
rownames(mean_prim) <- mean_prim$neighbor
mean_prim <- as.data.frame(mean_prim)


library(ComplexHeatmap)


neighbors <- c("CD11c.myeloid", "FOXP3.CD4.Tregs", "CD4.T.cells", "CD8.T.cells", "IBA1.CD11c.Macrophages", "IBA1.CD163.Macrophages", "CD163.Macrophages")

mean_int <- mean_int[which(mean_int$neighbor %in% neighbors),]
mean_prim <- mean_prim[which(mean_prim$neighbor %in% neighbors),]

mean_int$neighbor <- factor(mean_int$neighbor, levels=neighbors)
mean_prim$neighbor <- factor(mean_prim$neighbor, levels=neighbors)

library(RColorBrewer)
library(viridis)
library(circlize)
col_c = colorRamp2(c(-1, 0, 1), hcl_palette = "Blue-Red 3")

x_order <- c("neighbor","Epithelial", "Proliferating.epithelial", "EMT", "Proliferating.EMT", "CD163.Macrophages", "IBA1.CD163.Macrophages", "IBA1.CD11c.Macrophages", "CD11c.myeloid", "CD8.T.cells", "CD4.T.cells", "FOXP3.CD4.Tregs", "Fibroblast", "Myofibroblast", "Endothelial.cell", "SMA.CD31.positive.cell", "SMA.Desmin.positive.cell", "Desmin.positive.cell")
mean_int <- mean_int[,x_order]
mean_prim <- mean_prim[,x_order]

pdf(paste0("D:/Sciset/scimap/plots/revision/heatmap_interaction_values_",i,".pdf"))

hmap1 <- Heatmap(mean_int[,-1], name="Spatial interaction", 
        row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_row_names = T, row_title_gp = gpar(fontsize = 10, fontface = "bold"), cluster_rows = F, cluster_columns = F,
        show_row_dend=F, show_column_dend=F, column_title = "IDS samples - celltypes", row_title = "neighbors",
        height = unit(3.5, "cm") , width = unit(7.3, "cm"),border="white",
        rect_gp = gpar(col = "white", lwd = 2), col=col_c)



hmap2 <- Heatmap(mean_prim[,-1], name="Spatial interaction", 
        row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_row_names = T, row_title_gp = gpar(fontsize = 10, fontface = "bold"), cluster_rows = F, cluster_columns = F,
        show_row_dend=F, show_column_dend=F, column_title = "chemo-naive samples - celltypes", row_title = "neighbors",
        height = unit(3.5, "cm") , width = unit(7.3, "cm"),border="white",
        rect_gp = gpar(col = "white", lwd = 2), col=col_c)

pdf(paste0("D:/Sciset/scimap/plots/revision/heatmap_interaction_values_",i,".pdf"))

draw(hmap1)
draw(hmap2)

dev.off()


}
