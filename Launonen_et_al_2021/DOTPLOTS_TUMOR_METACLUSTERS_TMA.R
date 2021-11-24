#######################################################################################################################
############################### FUNCTIONAL MARKER DOTPLOTS ############################################################
#######################################################################################################################

#For Figure 3 and Supplementary Figure 3

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(data.table)

####################################### HR status ####################################################################################
tumor_metaclusters <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]
colnames(tumor_metaclusters)[46] <- "Metacluster"

tumor_metaclusters[which(tumor_metaclusters$HR_defect == "1"), "HR_defect"] <- "BRCA1/2 mutated"
tumor_metaclusters[which(tumor_metaclusters$HR_defect == "0"), "HR_defect"] <- "HRwt"

p_values <- tumor_metaclusters %>% group_by(Metacluster)%>%
  summarise_each(funs(wilcox.test(.[HR_defect == "BRCA1/2 mutated"], .[HR_defect == "HRwt"])$p.value), vars = CD11c:CD45)
p_values <- as.data.frame(p_values)

rownames(p_values) <- p_values[, 1]
p_values <- p_values[, -c(1)]
colnames(p_values) <- colnames(tumor_metaclusters)[8:31]
p_values <- p_values[, c("cCasp3", "pSTAT1","Ki67", "PDL1", "P21")]

ann <- read.csv("TMA_clinicaldata.csv")
ann <- ann[, c("Identifier","Category", "PFI_time")]
ann$PFI_time[which(ann$`PFI_time`> 365)] <- "long"
ann$PFI_time[which(ann$`PFI_time`> 182 & ann$`PFI_time`< 365)] <- "short"
ann$PFI_time[which(ann$`PFI_time`< 182)] <- "short"


#tumor metaclusters
tumor_metaclusters <- merge(tumor_metaclusters, ann, by.x = "Sample", by.y="Identifier")


tumor_metaclusters_new <- exp(tumor_metaclusters[, c(8:31)])
tumor_metaclusters_new <- cbind(tumor_metaclusters_new, tumor_metaclusters[, c(1, 45,46, 49, 50)])

median_expression <- tumor_metaclusters_new %>% 
  group_by(Metacluster, HR_defect) %>%
  summarize(
            median_cCasp3 = median(cCasp3),
            median_pSTAT1 = median(pSTAT1),
            median_Ki67 = median(Ki67),
            median_PDL1 = median(PDL1),
            median_P21 = median(P21))

library(tidyr)
median_expression_long <- pivot_longer(median_expression, cols = 3:7, names_to = "Marker", values_to = "Expression")

median_expression_wide <- pivot_wider(median_expression_long, names_from = "HR_defect", values_from = "Expression")

#fold change

median_expression_wide$fold_change <- median_expression_wide_new$`BRCA1/2 mutated`/median_expression_wide_new$HRwt
median_expression_wide$log_fold_change <- log2(median_expression_wide$fold_change)

fold_change_data <- median_expression_wide[, c("Metacluster", "Marker", "log_fold_change")]


log_fold_change_data_wide <- pivot_wider(fold_change_data, names_from = "Metacluster", values_from = "log_fold_change")

log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]


log_fold_change_data_wide <- as.data.table(t(log_fold_change_data_wide), keep.colnames = T, keep.rownames = T)
log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]

dotplot_data <- log_fold_change_data_wide

colnames(dotplot_data) <- c("foldchange_cCasp3", "foldchange_pSTAT1", "foldchange_Ki67", "foldchange_PDL1", "foldchange_P21")

dotplot_data <- as.data.frame(dotplot_data)

setDT(dotplot_data, keep.rownames = "Subtype")[]
colnames(dotplot_data) <- c("Metacluster", "cCasp3", "pSTAT1", "Ki67", "PDL1", "P21")
dotplot_data_long <- pivot_longer(dotplot_data, cols = 2:6, names_to = "Foldchange")
colnames(dotplot_data_long) <- c("Metacluster", "Marker", "log_foldchange")

setDT(p_values, keep.rownames = "Metacluster")[]
colnames(p_values) <- c("Metacluster", "cCasp3", "pSTAT1", "Ki67", "PDL1","P21")
p_values_long <- pivot_longer(p_values, cols = 2:6, names_to = "pvalue")
colnames(p_values_long) <- c("Metacluster", "Marker", "pvalue")


poly.results <- cbind(dotplot_data_long, p_values_long, by="Metacluster")
poly.results <- poly.results[, c(1, 2, 3, 6)]

row_order <- rev(c("Epithelial", "Proliferating epithelial", "EMT", "Proliferating EMT", "Mesenchymal",
                              "Functional epithelial", "Apoptotic"))



truncate.df <- function(df, na.cutoff=0.05, na.var="pvalue", na.var.boundary=1e-50, range.lims=c(-0.1, 0.1), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(poly.results) 
p <- ggplot(to.plot,  aes(x=factor(Marker), y=factor(Metacluster, level=row_order))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90", breaks=seq(-0.1,0.1,0.05), limits=c(-0.1, 0.1)) +
  scale_size_area("p-values", trans="log10",max_size = 2.5, breaks=c(1e-50, 1e-20, 1e-10, 1e-5, 1e-1, 0.05), limits=c(1e-50,0.05)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"), aspect.ratio=1, panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)


################################ WITH PFI ######################################################################################

### FIRST BRCA1/2MUT TUMORS

tumor_metaclusters <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]
colnames(tumor_metaclusters)[46] <- "Metacluster"

tumor_metaclusters[which(tumor_metaclusters$HR_defect == "1"), "HR_defect"] <- "BRCA1/2 mutated"
tumor_metaclusters[which(tumor_metaclusters$HR_defect == "0"), "HR_defect"] <- "HRwt"

tumor_metaclusters <- tumor_metaclusters[which(tumor_metaclusters$HR_defect == "BRCA1/2 mutated"),]


#tumor metaclusters
tumor_metaclusters <- merge(tumor_metaclusters, ann, by.x = "Sample", by.y="Identifier")


p_values <- tumor_metaclusters %>% group_by(Metacluster)%>%
  summarise_each(funs(wilcox.test(.[PFI_time == "long"], .[PFI_time == "short"])$p.value), vars = CD11c:CD45)
p_values <- as.data.frame(p_values)

rownames(p_values) <- p_values[, 1]
p_values <- p_values[, -c(1)]
colnames(p_values) <- colnames(tumor_metaclusters)[8:31]
p_values <- p_values[, c("cCasp3", "pSTAT1","Ki67", "PDL1", "P21")]


#tumor metaclusters
tumor_metaclusters <- merge(tumor_metaclusters, ann, by.x = "Sample", by.y="Identifier")


tumor_metaclusters_new <- exp(tumor_metaclusters[, c(8:31)])
tumor_metaclusters_new <- cbind(tumor_metaclusters_new, tumor_metaclusters[, c(1, 45,46, 49, 50)])

median_expression <- tumor_metaclusters_new %>% 
  group_by(Metacluster, PFI_time) %>%
  summarize(
    median_cCasp3 = median(cCasp3),
    median_pSTAT1 = median(pSTAT1),
    median_Ki67 = median(Ki67),
    median_PDL1 = median(PDL1),
    median_P21 = median(P21))

library(tidyr)
median_expression_long <- pivot_longer(median_expression, cols = 3:7, names_to = "Marker", values_to = "Expression")

median_expression_wide <- pivot_wider(median_expression_long, names_from = "PFI_time", values_from = "Expression")

#fold change

median_expression_wide$fold_change <- median_expression_wide_new$`long`/median_expression_wide_new$short
median_expression_wide$log_fold_change <- log2(median_expression_wide$fold_change)

fold_change_data <- median_expression_wide[, c("Metacluster", "Marker", "log_fold_change")]


log_fold_change_data_wide <- pivot_wider(fold_change_data, names_from = "Metacluster", values_from = "log_fold_change")

log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]


log_fold_change_data_wide <- as.data.table(t(log_fold_change_data_wide), keep.colnames = T, keep.rownames = T)
log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]

dotplot_data <- log_fold_change_data_wide

colnames(dotplot_data) <- c("foldchange_cCasp3", "foldchange_pSTAT1", "foldchange_Ki67", "foldchange_PDL1", "foldchange_P21")

dotplot_data <- as.data.frame(dotplot_data)

setDT(dotplot_data, keep.rownames = "Subtype")[]
colnames(dotplot_data) <- c("Metacluster", "cCasp3", "pSTAT1", "Ki67", "PDL1", "P21")
dotplot_data_long <- pivot_longer(dotplot_data, cols = 2:6, names_to = "Foldchange")
colnames(dotplot_data_long) <- c("Metacluster", "Marker", "log_foldchange")

setDT(p_values, keep.rownames = "Metacluster")[]
colnames(p_values) <- c("Metacluster", "cCasp3", "pSTAT1", "Ki67", "PDL1","P21")
p_values_long <- pivot_longer(p_values, cols = 2:6, names_to = "pvalue")
colnames(p_values_long) <- c("Metacluster", "Marker", "pvalue")


poly.results <- cbind(dotplot_data_long, p_values_long, by="Metacluster")
poly.results <- poly.results[, c(1, 2, 3, 6)]

row_order <- rev(c("Epithelial", "Proliferating epithelial", "EMT", "Proliferating EMT", "Mesenchymal",
                   "Functional epithelial", "Apoptotic"))



truncate.df <- function(df, na.cutoff=0.05, na.var="pvalue", na.var.boundary=1e-50, range.lims=c(-0.1, 0.1), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(poly.results) 
p <- ggplot(to.plot,  aes(x=factor(Marker), y=factor(Metacluster, level=row_order))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90", breaks=seq(-0.1,0.1,0.05), limits=c(-0.1, 0.1)) +
  scale_size_area("p-values", trans="log10",max_size = 2.5, breaks=c(1e-50, 1e-20, 1e-10, 1e-5, 1e-1, 0.05), limits=c(1e-50,0.05)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"), aspect.ratio=1, panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)




#HRWT TUMORS



tumor_metaclusters <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]
colnames(tumor_metaclusters)[46] <- "Metacluster"

tumor_metaclusters[which(tumor_metaclusters$HR_defect == "1"), "HR_defect"] <- "BRCA1/2 mutated"
tumor_metaclusters[which(tumor_metaclusters$HR_defect == "0"), "HR_defect"] <- "HRwt"

tumor_metaclusters <- tumor_metaclusters[which(tumor_metaclusters$HR_defect == "HRwt"),]


#tumor metaclusters
tumor_metaclusters <- merge(tumor_metaclusters, ann, by.x = "Sample", by.y="Identifier")


p_values <- tumor_metaclusters %>% group_by(Metacluster)%>%
  summarise_each(funs(wilcox.test(.[PFI_time == "long"], .[PFI_time == "short"])$p.value), vars = CD11c:CD45)
p_values <- as.data.frame(p_values)

rownames(p_values) <- p_values[, 1]
p_values <- p_values[, -c(1)]
colnames(p_values) <- colnames(tumor_metaclusters)[8:31]
p_values <- p_values[, c("cCasp3", "pSTAT1","Ki67", "PDL1", "P21")]


#tumor metaclusters
tumor_metaclusters <- merge(tumor_metaclusters, ann, by.x = "Sample", by.y="Identifier")


tumor_metaclusters_new <- exp(tumor_metaclusters[, c(8:31)])
tumor_metaclusters_new <- cbind(tumor_metaclusters_new, tumor_metaclusters[, c(1, 45,46, 49, 50)])

median_expression <- tumor_metaclusters_new %>% 
  group_by(Metacluster, PFI_time) %>%
  summarize(
    median_cCasp3 = median(cCasp3),
    median_pSTAT1 = median(pSTAT1),
    median_Ki67 = median(Ki67),
    median_PDL1 = median(PDL1),
    median_P21 = median(P21))

library(tidyr)
median_expression_long <- pivot_longer(median_expression, cols = 3:7, names_to = "Marker", values_to = "Expression")

median_expression_wide <- pivot_wider(median_expression_long, names_from = "PFI_time", values_from = "Expression")

#fold change

median_expression_wide$fold_change <- median_expression_wide_new$`long`/median_expression_wide_new$short
median_expression_wide$log_fold_change <- log2(median_expression_wide$fold_change)

fold_change_data <- median_expression_wide[, c("Metacluster", "Marker", "log_fold_change")]


log_fold_change_data_wide <- pivot_wider(fold_change_data, names_from = "Metacluster", values_from = "log_fold_change")

log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]


log_fold_change_data_wide <- as.data.table(t(log_fold_change_data_wide), keep.colnames = T, keep.rownames = T)
log_fold_change_data_wide <- as.data.frame(log_fold_change_data_wide)
rownames(log_fold_change_data_wide) <- log_fold_change_data_wide[, 1]
log_fold_change_data_wide <- log_fold_change_data_wide[, -c(1)]

dotplot_data <- log_fold_change_data_wide

colnames(dotplot_data) <- c("foldchange_cCasp3", "foldchange_pSTAT1", "foldchange_Ki67", "foldchange_PDL1", "foldchange_P21")

dotplot_data <- as.data.frame(dotplot_data)

setDT(dotplot_data, keep.rownames = "Subtype")[]
colnames(dotplot_data) <- c("Metacluster", "cCasp3", "pSTAT1", "Ki67", "PDL1", "P21")
dotplot_data_long <- pivot_longer(dotplot_data, cols = 2:6, names_to = "Foldchange")
colnames(dotplot_data_long) <- c("Metacluster", "Marker", "log_foldchange")

setDT(p_values, keep.rownames = "Metacluster")[]
colnames(p_values) <- c("Metacluster", "cCasp3", "pSTAT1", "Ki67", "PDL1","P21")
p_values_long <- pivot_longer(p_values, cols = 2:6, names_to = "pvalue")
colnames(p_values_long) <- c("Metacluster", "Marker", "pvalue")


poly.results <- cbind(dotplot_data_long, p_values_long, by="Metacluster")
poly.results <- poly.results[, c(1, 2, 3, 6)]

row_order <- rev(c("Epithelial", "Proliferating epithelial", "EMT", "Proliferating EMT", "Mesenchymal",
                   "Functional epithelial", "Apoptotic"))



truncate.df <- function(df, na.cutoff=0.05, na.var="pvalue", na.var.boundary=1e-50, range.lims=c(-0.1, 0.1), range.var="log_foldchange"){
  df[which(df[,na.var] > na.cutoff), na.var] <- NA
  df[which(df[,na.var] < na.var.boundary), na.var] <- na.var.boundary
  df[,range.var] <- pmax( range.lims[1], pmin( df[,range.var], range.lims[2]))
  return(df)
}
to.plot <- truncate.df(poly.results) 
p <- ggplot(to.plot,  aes(x=factor(Marker), y=factor(Metacluster, level=row_order))) +
  geom_point(aes(color=log_foldchange, size=pvalue)) + theme_classic() + xlab(NULL) + ylab(NULL) + labs(color="fold change (log)") +
  scale_color_gradient2(low = "#0b71b0", high="#cc2127", mid = "gray90", breaks=seq(-0.1,0.1,0.05), limits=c(-0.1, 0.1)) +
  scale_size_area("p-values", trans="log10",max_size = 2.5, breaks=c(1e-50, 1e-20, 1e-10, 1e-5, 1e-1, 0.05), limits=c(1e-50,0.05)) +
  theme(axis.text=element_text(size=rel(1.3)), axis.text.x = element_text(angle=45, hjust=1), strip.placement = "outside",strip.background = element_blank())+
  theme(plot.title = element_text(size = 12, face = "bold"), aspect.ratio=1, panel.border = element_rect(colour = "black", size=1.5, fill=NA))

print(p)



