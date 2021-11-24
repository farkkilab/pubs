####################################################################################################################
################################ BARPLOTS ##########################################################################
####################################################################################################################

# For Figure 1

library(dplyr)
library(ComplexHeatmap)
library(ggplot2)

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


all_celltypes[which(all_celltypes$GlobalCellType == "Endothelia"), "Subtype"] <- "Endothelia"
all_celltypes[which(all_celltypes$Subtype == "Stroma_Endothelia"), "Subtype"] <- "Stroma"


TSEI_percentages <- all_celltypes %>% group_by(Sample, Subtype) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
plotorder <- TSEI_percentages[which(TSEI_percentages$Subtype == "Tumor"),c("Sample","proportion")]
TSEI_percentages <- merge(TSEI_percentages, plotorder, by="Sample")

df <- TSEI_percentages
#axis.text.x=element_blank(), axis.ticks.x = element_blank(),

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(4)
ggplot(df, aes(x = as.factor(reorder(Sample,-proportion.y)), y = proportion.x, fill = Subtype))+
  geom_bar(stat = "identity", colour="black", position = "stack") + scale_fill_manual(values=mycolors)+ xlab("Sample") + ylab("Proportion")+ 
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_discrete(expand = c(0,0)) + theme(legend.position="none")



stroma <- all_celltypes[which(all_celltypes$Subtype == "Stroma" | all_celltypes$Subtype == "Endothelia"),]

SE_percentages <- stroma %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
plotorder <- SE_percentages[which(SE_percentages$GlobalCellType == "Non-proliferative Stroma"),c("Sample","proportion")]
SE_percentages <- merge(SE_percentages, plotorder, by="Sample")

df <- SE_percentages


mycolors <- colorRampPalette(brewer.pal(9, "Set3"))(9)
mycolors <- c("#D9D9D9","#FB8072","#B3DE69","#FDB462","#80B1D3","#BEBADA","#FFFFB3","#8DD3C7","#FCCDE5")

ggplot(df, aes(x = as.factor(reorder(Sample,-proportion.y)), y = proportion.x, fill = factor(GlobalCellType, levels=c("Proliferative Stroma", "High-proliferative Stroma", "Low eccentricity", "Low-vimentin", "High-Vimentin","High-P21","Functional stroma", "Endothelia","Non-proliferative Stroma"))))+
  geom_bar(stat = "identity", position = "stack", colour="black", size=0.25) + scale_fill_manual(values=mycolors)+ xlab("Sample") + ylab("Proportion")+ 
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5), axis.ticks.x = element_blank(),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+ scale_y_discrete(expand = c(0,0)) + theme(legend.position="none")



immune <- all_celltypes[which(all_celltypes$Subtype == "Immune"),]


I_percentages <- immune %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
plotorder <- I_percentages[which(I_percentages$GlobalCellType == "IBA1+CD163+ Macrophages"),c("Sample","proportion")]
I_percentages <- merge(I_percentages, plotorder, by="Sample")

df <- I_percentages


mycolors <- colorRampPalette(brewer.pal(10, "Paired"))(10)
mycolors <- c("#A6CEE3", "#FB9A99", "#E31A1C","#33A02C", "#1F78B4", "#B2DF8A", "#FDBF6F", "#FF7F00", "#6A3D9A", "#CAB2D6")

ggplot(df, aes(x = as.factor(reorder(Sample,-proportion.y)), y = proportion.x, fill = factor(GlobalCellType, levels=rev(c("IBA1+CD163+ Macrophages","IBA1+CD163+CD11c+ Macrophages", "IBA1+CD11c+ Macrophages", "IBA1+ Macrophages", "CD163+ Macrophages", "CD11c+ APC", "CD4+ Effector T-cells", "FOXP3+CD4+ T-regs", "CD8+ T-cells","B-cells")))))+
  geom_bar(stat = "identity", position = "stack", colour="black", size=0.25) + scale_fill_manual(values=mycolors)+ xlab("Sample") + ylab("Proportion")+ 
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5), axis.ticks.x = element_blank(),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "white"))+ scale_y_discrete(expand = c(0,0))





#tumor
tumor <- all_celltypes[which(all_celltypes$Subtype == "Tumor"),]


T_percentages <- tumor %>% group_by(Sample, GlobalCellType) %>% summarise(n = n()) %>% mutate(proportion = n / sum(n))
plotorder <- T_percentages[which(T_percentages$GlobalCellType == "Proliferating epithelial"),c("Sample","proportion")]
plotorder <- as.data.frame(plotorder)

#60417 has no proliferating epithelial, add the the smeelst proportion to get it to plotorder
extra <- c("60417", 0.001)
extra <- as.data.frame(t(extra))
colnames(extra) <- c("Sample", "proportion")

plotorder <- rbind(plotorder, extra)

T_percentages <- merge(T_percentages, plotorder, by="Sample")

df <- T_percentages
df$proportion.y <- as.numeric(df$proportion.y)

mycolors <- colorRampPalette(brewer.pal(7, "Accent"))(7)
mycolors <- c("#FFFF99","#BF5B17", "#BEAED4","#FDC086","#F0027F", "#7FC97F", "#386CB0")

ggplot(df, aes(x = as.factor(reorder(Sample,-proportion.y)), y = proportion.x, fill = factor(GlobalCellType, levels=rev(c("Proliferating epithelial","Epithelial", "EMT", "Mesenchymal", "Proliferating EMT", "Functional epithelial" ,"Apoptotic")))))+
  geom_bar(stat = "identity", position = "stack", colour="black", size=0.25) + scale_fill_manual(values=mycolors)+ xlab("Sample") + ylab("Proportion")+ 
  guides(fill = guide_legend(title = "Cell type")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5), axis.ticks.x = element_blank(),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "white"))+ scale_y_discrete(expand = c(0,0))+ theme(legend.position="none")

#+ theme(legend.position="none")



#complex heatmap annotation bar for the barplots


all_data_TMA <- read.csv("TMA_clinicaldata.csv")
data_heatmap_TMA <- all_data_TMA[, c("Identifier","Category", "PFI_time")]
ann <- data_heatmap_TMA
ann$PFI_time[which(ann$`PFI_time`> 365)] <- "long"
ann$PFI_time[which(ann$`PFI_time`> 182 & ann$`PFI_time`< 365)] <- "short"
ann$PFI_time[which(ann$`PFI_time`< 182)] <- "short"
ann$PFI_time[which(ann$`PFI_time`== "56")] <- "short"
ann$PFI_time[which(ann$`PFI_time`== "55")] <- "short"
ann$Category <- as.character(ann$Category)
ann$Category[which(ann$`Category`== "HR")] <- "HRwt"

clinical <- ann
colnames(clinical) <- c("Sample", "HR status", "PFI")
clinical$PFI <- factor(clinical$PFI)

clinical <- clinical[is.element(clinical$Sample, tumor$Sample),]

colours<- list("HR status"=c("BRCA1"="royalblue", "BRCA2" = "blue","HRwt"="red3"), "PFI"=c("long" = "#3399CC", "medium"="#006699", "short"="#003366"))



T_percentages$proportion.y <- as.numeric(T_percentages$proportion.y)
T_percentages$Sample <- as.integer(T_percentages$Sample)

order_T = as.factor(reorder(T_percentages$Sample,-T_percentages$proportion.y))

order_I = as.factor(reorder(I_percentages$Sample,-I_percentages$proportion.y))

order_SE = as.factor(reorder(SE_percentages$Sample,-SE_percentages$proportion.y))

order_TSEI = as.factor(reorder(TSEI_percentages$Sample,-TSEI_percentages$proportion.y))


#tumor
#clinical$Sample <- clinical$Sample[names(levels(order_T))]

#clinical <- clinical[order(match(clinical$Sample, levels(order_T))), ]


rownames(clinical) <- clinical$Sample

clinical <- clinical[levels(order_T),]

Rowann <- HeatmapAnnotation(df=clinical[, c(2:3)], which="column", col=colours, 
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(0, "mm"))




#immune
#clinical$Sample <- clinical$Sample[as.factor(levels(order_I))]
clinical <- clinical[levels(order_I),]


Rowann <- HeatmapAnnotation(df=clinical[, c(2:3)], which="column", col=colours, 
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(0, "mm"))



#stroma
#clinical$Sample <- clinical$Sample[as.factor(levels(order_SE))]

clinical <- clinical[levels(order_SE),]

Rowann <- HeatmapAnnotation(df=clinical[, c(2:3)], which="column", col=colours, 
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(0, "mm"))


clinical <- clinical[levels(order_TSEI),]

Rowann <- HeatmapAnnotation(df=clinical[, c(2:3)], which="column", col=colours, 
                            annotation_name_gp = gpar(fontsize=10,fontface = "bold"), gap = unit(0, "mm"))


